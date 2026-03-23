import os
import glob
import numpy as np
import wrf
import xarray as xr
import metpy.calc as mpcalc
from metpy.units import units
from netCDF4 import Dataset
from wrf import getvar, latlon_coords, to_np, interplevel
import multiprocessing as mp
import time
import pandas as pd
from datetime import datetime
import warnings

warnings.filterwarnings('ignore')

# --- CONFIGURAÇÃO ---
dia = datetime.now().strftime('%d')

DIR_ENTRADA = '/home/gilberto/data/WRF'
DIR_SAIDA   = '/home/gilberto/data/WRF/indices'
N_CORES     = 64

os.makedirs(DIR_SAIDA, exist_ok=True)


# --- FÍSICA ---

def scalar_moist_lapse(t_ref_val: float, p_ref_val: float) -> float:
    """Eleva uma parcela da pressão p_ref até 500 hPa pela adiabática saturada.
    Retorna a temperatura da parcela a 500 hPa (K), ou NaN em caso de falha."""
    try:
        if np.isnan(t_ref_val) or np.isnan(p_ref_val) or p_ref_val < 500.0:
            return np.nan
        res = mpcalc.moist_lapse(
            500.0 * units.hPa,
            t_ref_val * units.kelvin,
            reference_pressure=p_ref_val * units.hPa,
        )
        return float(res.m)
    except Exception:
        return np.nan


vec_moist_lapse = np.vectorize(scalar_moist_lapse)


def lift_parcels(lcl_t: np.ndarray, lcl_p: np.ndarray) -> np.ndarray:
    """Aplica vec_moist_lapse apenas nos pontos válidos, evitando overhead desnecessário."""
    result = np.full(lcl_t.shape, np.nan, dtype=np.float64)
    valid = ~(np.isnan(lcl_t) | np.isnan(lcl_p) | (lcl_p < 500.0))
    if valid.any():
        result[valid] = vec_moist_lapse(lcl_t[valid], lcl_p[valid])
    return result


def processar_arquivo(caminho: str):
    t_start = time.time()
    nome = os.path.basename(caminho)

    try:
        ncfile = Dataset(caminho)

        # 1. Variáveis 3D do WRF
        p   = getvar(ncfile, 'pres', units='hPa')
        tk  = getvar(ncfile, 'tk')                  # temperatura (K)
        td  = getvar(ncfile, 'td', units='K')        # ponto de orvalho (K)
        qv  = getvar(ncfile, 'QVAPOR')
        z   = getvar(ncfile, 'z', units='m')
        psfc = getvar(ncfile, 'PSFC')               # pressão superficial (Pa)
        ter = getvar(ncfile, 'ter', units='m')
        ua  = getvar(ncfile, 'ua', units='m s-1')
        va  = getvar(ncfile, 'va', units='m s-1')

        # Timestamp para o nome de saída
        v_time   = getvar(ncfile, 'times', timeidx=0)
        ts       = pd.to_datetime(str(v_time.values))
        nome_saida  = ts.strftime(f'wrf_indices_i{dia}_%Y_%m_%d_%H.nc')
        caminho_saida = os.path.join(DIR_SAIDA, nome_saida)

        # 2. Interpolação para níveis de pressão
        t500  = interplevel(tk, p, 500.0)
        t700  = interplevel(tk, p, 700.0)
        t850  = interplevel(tk, p, 850.0)
        td850 = interplevel(td, p, 850.0)
        td700 = interplevel(td, p, 700.0)

        # Superfície (nível mais baixo do modelo)
        t_sfc  = tk[0, :, :]
        td_sfc = td[0, :, :]
        p_sfc  = to_np(psfc) / 100.0  # Pa → hPa (2D)

        # 3. K-Index e Total Totals (°C)
        t500c, t700c, t850c = t500 - 273.15, t700 - 273.15, t850 - 273.15
        td850c = td850 - 273.15
        td700c = td700 - 273.15

        k_index     = (t850c - t500c) + td850c - (t700c - td700c)
        total_totals = (t850c - t500c) + (td850c - t500c)

        # 4. Showalter e Lifted via MetPy (consistente com ICON/MPAS)
        t850_flat  = to_np(t850).flatten()
        td850_flat = to_np(td850).flatten()
        t_sfc_flat = to_np(t_sfc).flatten()
        td_sfc_flat = to_np(td_sfc).flatten()
        p_sfc_flat = p_sfc.flatten()
        t500_flat  = to_np(t500).flatten()

        lcl_p_si, lcl_t_si = mpcalc.lcl(
            850 * units.hPa, t850_flat * units.kelvin, td850_flat * units.kelvin
        )
        lcl_p_li, lcl_t_li = mpcalc.lcl(
            p_sfc_flat * units.hPa, t_sfc_flat * units.kelvin, td_sfc_flat * units.kelvin
        )

        t_parcel_si = lift_parcels(lcl_t_si.m, lcl_p_si.m)
        shw_index   = (t500_flat - t_parcel_si).reshape(to_np(t500).shape)

        t_parcel_li = lift_parcels(lcl_t_li.m, lcl_p_li.m)
        lifted_index = (t500_flat - t_parcel_li).reshape(to_np(t500).shape)

        # 5. CAPE e CIN via wrf-python
        cape_3d = wrf.cape_2d(p, tk, qv, z, ter, psfc, ter_follow=True)
        cape = cape_3d[0, :]
        cin  = cape_3d[1, :]

        # 6. BRN (cisalhamento 0–6 km)
        u6km = interplevel(ua, z, 6000.0)
        v6km = interplevel(va, z, 6000.0)

        u_shear  = u6km - ua[0, :, :]
        v_shear  = v6km - va[0, :, :]
        shear_sq = u_shear**2 + v_shear**2

        brn = np.where(shear_sq > 0.1, to_np(cape) / (0.5 * to_np(shear_sq)), 0.0)
        brn = np.clip(brn, 0, 500)

        # 7. Salvar NetCDF
        lats, lons = latlon_coords(p)
        ds = xr.Dataset(
            data_vars={
                'CAPE':         (('lat', 'lon'), to_np(cape).astype(np.float32)),
                'CIN':          (('lat', 'lon'), to_np(cin).astype(np.float32)),
                'K_INDEX':      (('lat', 'lon'), to_np(k_index).astype(np.float32)),
                'TOTAL_TOTALS': (('lat', 'lon'), to_np(total_totals).astype(np.float32)),
                'SHOWALTER':    (('lat', 'lon'), shw_index.astype(np.float32)),
                'LIFTED_INDEX': (('lat', 'lon'), lifted_index.astype(np.float32)),
                'BRN':          (('lat', 'lon'), brn.astype(np.float32)),
            },
            coords={
                'time':      ts,
                'latitude':  (('lat', 'lon'), to_np(lats)),
                'longitude': (('lat', 'lon'), to_np(lons)),
            },
        )

        encoding = {v: {'zlib': True, 'complevel': 4} for v in ds.data_vars}
        ds.to_netcdf(caminho_saida, encoding=encoding)
        ncfile.close()

        return True, nome, nome_saida, time.time() - t_start

    except Exception as e:
        return False, nome, str(e), 0.0


if __name__ == '__main__':
    t0 = time.time()

    arquivos = sorted(glob.glob(os.path.join(DIR_ENTRADA, f'wrfout_i{dia}_d02_*.nc')))
    total = len(arquivos)

    if total == 0:
        print(f'Nenhum arquivo wrfout encontrado em {DIR_ENTRADA}')
        raise SystemExit(0)

    print(f'WRF: {total} arquivos | {N_CORES} núcleos | saída: {DIR_SAIDA}\n')

    with mp.Pool(N_CORES) as pool:
        for i, (ok, orig, novo, tempo) in enumerate(
            pool.imap_unordered(processar_arquivo, arquivos), 1
        ):
            status = '✓' if ok else '✗ ERRO'
            print(f'[{i}/{total}] {status} {orig} -> {novo} | {tempo:.2f}s')

    print(f'\nConcluído em {(time.time() - t0) / 60:.2f} minutos.')
