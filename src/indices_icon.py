import xarray as xr
import metpy.calc as mpcalc
from metpy.units import units
import numpy as np
import multiprocessing as mp
import time
import pandas as pd
import warnings
import os
import glob
import sys
import argparse
from datetime import datetime

warnings.filterwarnings('ignore')

# --- CONFIGURAÇÃO ---
dia = datetime.now().strftime('%d')
dir_date = datetime.now().strftime('%Y%m%d')

dir_dados = f'/home/gilberto/data/ICON/{dir_date}00'
dir_saida = '/home/gilberto/data/ICON/indices'
N_CORES_MAX = 64

os.makedirs(dir_saida, exist_ok=True)


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


def calc_indices(ds, tag: str) -> tuple:
    """Calcula todos os índices a partir de um dataset xarray com níveis de pressão."""

    t500 = ds['temp'].sel(plev=50000).values.flatten()
    t700 = ds['temp'].sel(plev=70000).values.flatten()
    t850 = ds['temp'].sel(plev=85000).values.flatten()

    rh850 = ds['rh'].sel(plev=85000).values.flatten() / 100.0
    rh700 = ds['rh'].sel(plev=70000).values.flatten() / 100.0

    t2m  = ds['t_2m'].values.flatten()
    td2m = ds['td_2m'].values.flatten()
    psfc = ds['pres_sfc'].values.flatten() / 100.0  # Pa → hPa

    # Unidades MetPy
    t850_u  = t850  * units.kelvin
    t700_u  = t700  * units.kelvin
    t2m_u   = t2m   * units.kelvin
    td2m_u  = td2m  * units.kelvin
    psfc_u  = psfc  * units.hPa

    td850_u = mpcalc.dewpoint_from_relative_humidity(t850_u, rh850)
    td700_u = mpcalc.dewpoint_from_relative_humidity(t700_u, rh700)

    # LCL
    lcl_p_si, lcl_t_si = mpcalc.lcl(850 * units.hPa, t850_u, td850_u)
    lcl_p_li, lcl_t_li = mpcalc.lcl(psfc_u, t2m_u, td2m_u)

    # Showalter e Lifted (gargalo — só processa pontos válidos)
    print(f"  -> [{tag}] Calculando Showalter...", flush=True)
    t_parcel_si = lift_parcels(lcl_t_si.m, lcl_p_si.m)
    sw = t500 - t_parcel_si

    print(f"  -> [{tag}] Calculando Lifted...", flush=True)
    t_parcel_li = lift_parcels(lcl_t_li.m, lcl_p_li.m)
    li = t500 - t_parcel_li

    # Índices empíricos (°C)
    t500c, t700c, t850c = t500 - 273.15, t700 - 273.15, t850 - 273.15
    td850c = td850_u.to('degC').m
    td700c = td700_u.to('degC').m

    k_idx = (t850c - t500c) + td850c - (t700c - td700c)
    tt_idx = (t850c - t500c) + (td850c - t500c)

    # BRN
    u_shear = ds['u'].sel(plev=50000).values.flatten() - ds['u_10m'].values.flatten()
    v_shear = ds['v'].sel(plev=50000).values.flatten() - ds['v_10m'].values.flatten()
    shear_sq = u_shear**2 + v_shear**2

    with np.errstate(divide='ignore', invalid='ignore'):
        brn = np.where(shear_sq > 0.1,
                       ds['cape_mu'].values.flatten() / (0.5 * shear_sq),
                       0.0)

    return k_idx, tt_idx, sw, li, brn


def process_file(file_pl: str):
    t_start = time.time()
    file_sfc = file_pl.replace('_PL', '')
    tag = os.path.basename(file_pl).split('_')[-1].replace('.nc', '')

    if not os.path.exists(file_sfc):
        print(f"[{tag}] AVISO: arquivo de superfície não encontrado, pulando.", flush=True)
        return

    try:
        print(f"[{tag}] Iniciando...", flush=True)
        with xr.open_dataset(file_pl) as ds_pl, xr.open_dataset(file_sfc) as ds_sfc:
            ds = xr.merge([ds_pl, ds_sfc]).load()

        ts = pd.to_datetime(np.atleast_1d(ds.time.values)[0])
        fname_out = ts.strftime(f"icon_indices_i{dia}_%Y_%m_%d_%H.nc")
        path_out  = os.path.join(dir_saida, fname_out)

        k, tt, sw, li, brn = calc_indices(ds, tag)

        n_lat, n_lon = ds.dims['lat'], ds.dims['lon']
        shape = (n_lat, n_lon)

        out_ds = xr.Dataset(
            data_vars={
                'CAPE':      (('lat', 'lon'), ds['cape_mu'].values.reshape(shape).astype(np.float32)),
                'CIN':       (('lat', 'lon'), ds['cin_mu'].values.reshape(shape).astype(np.float32)),
                'K_INDEX':   (('lat', 'lon'), k.reshape(shape).astype(np.float32)),
                'TOTALS':    (('lat', 'lon'), tt.reshape(shape).astype(np.float32)),
                'SHOWALTER': (('lat', 'lon'), sw.reshape(shape).astype(np.float32)),
                'LIFTED':    (('lat', 'lon'), li.reshape(shape).astype(np.float32)),
                'BRN':       (('lat', 'lon'), brn.reshape(shape).astype(np.float32)),
            },
            coords={'time': ds.time, 'lat': ds.lat, 'lon': ds.lon},
        )

        encoding = {v: {'zlib': True, 'complevel': 4} for v in out_ds.data_vars}
        out_ds.to_netcdf(path_out, encoding=encoding)

        print(f"✓ {fname_out} ({time.time() - t_start:.1f}s)", flush=True)
        ds.close()

    except Exception as e:
        print(f"!!! ERRO [{tag}]: {e}", flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calcula índices de instabilidade do ICON (distribuído).')
    parser.add_argument('--rank', type=int, default=0, help='ID do nó atual (0-based)')
    parser.add_argument('--size', type=int, default=1, help='Número total de nós')
    args = parser.parse_args()

    all_files = sorted(glob.glob(os.path.join(dir_dados, 'ICON_LAM_DOM02_PL_*.nc')))
    total = len(all_files)

    if total == 0:
        print(f'Nenhum arquivo encontrado em {dir_dados}')
        sys.exit(0)

    my_files = np.array_split(all_files, args.size)[args.rank].tolist()
    n_workers = min(N_CORES_MAX, len(my_files)) if my_files else 1

    print('=' * 60)
    print(f'NÓ {args.rank + 1}/{args.size} | Arquivos: {len(my_files)} | Workers: {n_workers}')
    print('=' * 60)

    t0 = time.time()
    with mp.Pool(n_workers) as pool:
        pool.map(process_file, my_files)

    print('=' * 60)
    print(f'NÓ {args.rank + 1} CONCLUÍDO EM {(time.time() - t0) / 60:.2f} MIN')
    print('=' * 60)
