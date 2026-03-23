import xarray as xr
import metpy.calc as mpcalc
from metpy.units import units
import numpy as np
import multiprocessing as mp
import time
import pandas as pd
import warnings
import os
import argparse
import sys
from datetime import datetime

warnings.filterwarnings('ignore')

# --- CONFIGURAÇÃO ---
dia   = datetime.now().strftime('%d')
agora = datetime.now()

INPUT_FILE  = agora.strftime('/home/gilberto/data/MPAS/mpas_regrid3km_%Y-%m-%d-00.nc')
OUTPUT_DIR  = '/home/gilberto/data/MPAS/indices'
N_CORES_MAX = 64

os.makedirs(OUTPUT_DIR, exist_ok=True)


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


def calc_indices(ds_t, tag: str) -> tuple:
    """Calcula todos os índices a partir de um dataset xarray de um único timestep."""

    t500  = ds_t['temperature_500hPa'].values.flatten()
    t700  = ds_t['temperature_700hPa'].values.flatten()
    t850  = ds_t['temperature_850hPa'].values.flatten()
    tsfc  = ds_t['temperature_surface'].values.flatten()

    td850 = ds_t['dewpoint_850hPa'].values.flatten()
    td700 = ds_t['dewpoint_700hPa'].values.flatten()
    tdsfc = ds_t['dewpoint_surface'].values.flatten()

    psfc  = ds_t['mslp'].values.flatten() / 100.0  # Pa → hPa

    # LCL
    lcl_p_si, lcl_t_si = mpcalc.lcl(
        850 * units.hPa, t850 * units.kelvin, td850 * units.kelvin
    )
    lcl_p_li, lcl_t_li = mpcalc.lcl(
        psfc * units.hPa, tsfc * units.kelvin, tdsfc * units.kelvin
    )

    # Showalter e Lifted (só processa pontos válidos)
    print(f"  -> [{tag}] Calculando parcelas...", flush=True)
    t_parcel_si = lift_parcels(lcl_t_si.m, lcl_p_si.m)
    sw = t500 - t_parcel_si

    t_parcel_li = lift_parcels(lcl_t_li.m, lcl_p_li.m)
    li = t500 - t_parcel_li

    # Índices empíricos (°C)
    t500c, t700c, t850c = t500 - 273.15, t700 - 273.15, t850 - 273.15
    td850c, td700c = td850 - 273.15, td700 - 273.15

    k_idx  = (t850c - t500c) + td850c - (t700c - td700c)
    tt_idx = (t850c - t500c) + (td850c - t500c)

    # BRN
    u_shear  = ds_t['uzonal_6km'].values.flatten()    - ds_t['uzonal_surface'].values.flatten()
    v_shear  = ds_t['umeridional_6km'].values.flatten() - ds_t['umeridional_surface'].values.flatten()
    shear_sq = u_shear**2 + v_shear**2

    with np.errstate(divide='ignore', invalid='ignore'):
        brn = np.where(shear_sq > 0.1,
                       ds_t['cape'].values.flatten() / (0.5 * shear_sq),
                       0.0)

    return k_idx, tt_idx, sw, li, brn


def process_timestep(time_idx: int):
    t_start = time.time()

    try:
        with xr.open_dataset(INPUT_FILE, chunks={}) as ds:
            ds_t = ds.isel(time=time_idx).load()

        ts    = pd.to_datetime(str(ds_t.time.values))
        tag   = f"H{ts.strftime('%H')}"

        print(f"[{tag}] Carregado. Iniciando cálculos...", flush=True)

        k, tt, sw, li, brn = calc_indices(ds_t, tag)

        fname_out = ts.strftime(f"mpas_indices_i{dia}_%Y_%m_%d_%H.nc")
        path_out  = os.path.join(OUTPUT_DIR, fname_out)

        n_lat = ds_t.latitude.size
        n_lon = ds_t.longitude.size
        shape = (n_lat, n_lon)

        out_ds = xr.Dataset(
            data_vars={
                'CAPE':      (('latitude', 'longitude'), ds_t['cape'].values.astype(np.float32)),
                'CIN':       (('latitude', 'longitude'), ds_t['cin'].values.astype(np.float32)),
                'K_INDEX':   (('latitude', 'longitude'), k.reshape(shape).astype(np.float32)),
                'TOTALS':    (('latitude', 'longitude'), tt.reshape(shape).astype(np.float32)),
                'SHOWALTER': (('latitude', 'longitude'), sw.reshape(shape).astype(np.float32)),
                'LIFTED':    (('latitude', 'longitude'), li.reshape(shape).astype(np.float32)),
                'BRN':       (('latitude', 'longitude'), brn.reshape(shape).astype(np.float32)),
            },
            coords={'time': ds_t.time, 'latitude': ds_t.latitude, 'longitude': ds_t.longitude},
        )

        encoding = {v: {'zlib': True, 'complevel': 4} for v in out_ds.data_vars}
        out_ds.to_netcdf(path_out, encoding=encoding)

        print(f"✓ {fname_out} ({time.time() - t_start:.1f}s)", flush=True)

    except Exception as e:
        import traceback
        print(f"!!! ERRO no timestep {time_idx}: {e}", flush=True)
        traceback.print_exc()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calcula índices de instabilidade do MPAS (distribuído).')
    parser.add_argument('--rank', type=int, default=0, help='ID do nó atual (0-based)')
    parser.add_argument('--size', type=int, default=1, help='Número total de nós')
    args = parser.parse_args()

    if not os.path.exists(INPUT_FILE):
        print(f'ERRO: arquivo de entrada não encontrado: {INPUT_FILE}')
        sys.exit(1)

    with xr.open_dataset(INPUT_FILE) as ds_meta:
        total_steps = ds_meta.sizes['time']

    if total_steps == 0:
        print('Nenhum timestep encontrado no arquivo de entrada.')
        sys.exit(0)

    my_indices = np.array_split(np.arange(total_steps), args.size)[args.rank].tolist()
    n_workers  = min(N_CORES_MAX, len(my_indices)) if my_indices else 1

    print('=' * 60)
    print(f'MPAS NÓ {args.rank + 1}/{args.size} | Timesteps: {my_indices} | Workers: {n_workers}')
    print(f'Arquivo: {os.path.basename(INPUT_FILE)}')
    print('=' * 60)

    t0 = time.time()
    with mp.Pool(n_workers) as pool:
        pool.map(process_timestep, my_indices)

    print('=' * 60)
    print(f'NÓ {args.rank + 1} CONCLUÍDO EM {(time.time() - t0) / 60:.2f} MIN')
    print('=' * 60)
