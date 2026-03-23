# Instability Indices Calculator

**Author:** Gilberto Ricardo Bonatti — Especialista em Modelagem Numérica

Batch processing system for computing atmospheric convective instability indices from three operational NWP models: **ICON**, **MPAS**, and **WRF**. Designed for distributed execution on HPC clusters (tested on Thor cluster, 64 cores/node).

---

## Indices Computed

| Index | Acronym | Levels Used | Method |
|-------|---------|-------------|--------|
| Convective Available Potential Energy | CAPE | Surface → EL | From model output |
| Convective Inhibition | CIN | Surface → LFC | From model output |
| K-Index | KI | 850, 700, 500 hPa | Empirical formula |
| Total Totals | TT | 850, 500 hPa | Empirical formula |
| Showalter Index | SI | 850 → 500 hPa | Parcel theory (MetPy) |
| Lifted Index | LI | Surface → 500 hPa | Parcel theory (MetPy) |
| Bulk Richardson Number | BRN | Surface → 6 km | CAPE / wind shear |

For complete physical definitions and equations, see [Documentation.md](Documentation.md).

---

## Models Supported

### ICON (Icosahedral Nonhydrostatic)
- **Input pattern:** `ICON_LAM_DOM02_PL_*.nc` (pressure levels) + surface file
- **Distribution:** by file, across multiple cluster nodes
- **Launcher:** `src/roda_indices_icon.sh`

### MPAS (Model for Prediction Across Scales)
- **Input pattern:** `mpas_regrid3km_YYYY-MM-DD-00.nc` (single daily file, multiple timesteps)
- **Distribution:** by timestep, across multiple cluster nodes
- **Launcher:** `src/roda_indices_mpas.sh`

### WRF (Weather Research and Forecasting)
- **Input pattern:** `wrfout_i{day}_d02_*.nc` (one file per timestep)
- **Distribution:** by file, across multiple cluster nodes
- **Launcher:** `src/roda_indices_wrf.sh`

---

## Project Structure

```
instability-indices/
├── README.md
├── Documentation.md
├── .gitignore
└── src/
    ├── indices_icon.py        # ICON index calculator
    ├── indices_mpas.py        # MPAS index calculator
    ├── indices_wrf.py         # WRF index calculator
    ├── roda_indices_icon.sh   # ICON cluster launcher (4 nodes)
    ├── roda_indices_mpas.sh   # MPAS cluster launcher (4 nodes)
    └── roda_indices_wrf.sh    # WRF cluster launcher (4 nodes)
```

---

## Requirements

```
python >= 3.8
xarray
metpy
numpy
pandas
netCDF4
wrf-python   # WRF only
```

Install with conda (recommended):
```bash
conda install -c conda-forge xarray metpy numpy pandas netCDF4
conda install -c conda-forge wrf-python  # for WRF
```

---

## Usage

### ICON

**Single node:**
```bash
python src/indices_icon.py --rank 0 --size 1
```

**Distributed across 4 nodes (via launcher):**
```bash
bash src/roda_indices_icon.sh
```

### MPAS

**Single node:**
```bash
python src/indices_mpas.py --rank 0 --size 1
```

**Distributed across 4 nodes (via launcher):**
```bash
bash src/roda_indices_mpas.sh
```

### WRF

**Single node:**
```bash
python src/indices_wrf.py --rank 0 --size 1
```

**Distributed across 4 nodes (via launcher):**
```bash
bash src/roda_indices_wrf.sh
```

---

## Configuration

Edit the configuration block at the top of each script:

| Parameter | Description |
|-----------|-------------|
| `dir_dados` / `INPUT_FILE` | Path to model input data |
| `dir_saida` / `OUTPUT_DIR` | Path to output directory |
| `N_CORES_MAX` / `n_cores` | Max parallel workers per node |

Cluster nodes are configured in the launcher scripts (`NOS` array).

---

## Output Format

All scripts produce **NetCDF4** files with zlib compression (level 4).

**Naming convention:**
```
{model}_indices_i{day}_{YYYY}_{MM}_{DD}_{HH}.nc
```
Example: `icon_indices_i23_2024_03_15_12.nc`

**Variables in each output file:**

| Variable | Type | Description |
|----------|------|-------------|
| `CAPE` | float32 | Convective Available Potential Energy (J/kg) |
| `CIN` | float32 | Convective Inhibition (J/kg) |
| `K_INDEX` | float32 | K-Index (°C) |
| `TOTALS` / `TOTAL_TOTALS` | float32 | Total Totals Index (°C) |
| `SHOWALTER` | float32 | Showalter Index (°C) |
| `LIFTED` / `LIFTED_INDEX` | float32 | Lifted Index (°C) |
| `BRN` | float32 | Bulk Richardson Number (dimensionless) |

**Coordinates:** `time`, `latitude`, `longitude` (2D arrays for WRF).

---

## Architecture

### Parallel Processing Strategy

```
Cluster Node 0 ─┐
Cluster Node 1 ─┼─► SSH launcher ──► Python script (--rank N --size M)
Cluster Node 2 ─┤                         │
Cluster Node 3 ─┘                    mp.Pool (64 workers)
                                          │
                                     Per file/timestep:
                                     Load → Compute → Save NetCDF
```

- **ICON:** distributes files across nodes; each node runs a local pool
- **MPAS:** distributes timesteps across nodes; each node runs a local pool
- **WRF:** distributes files across nodes; each node runs a local pool

### Moist Adiabatic Lapse Rate

The Showalter and Lifted indices require lifting a parcel along the moist adiabat from the LCL to 500 hPa. This is implemented via `metpy.calc.moist_lapse()` wrapped in a scalar function and vectorized with `numpy.vectorize` to handle NaN propagation and physical constraints robustly.

---

## Paths (Cluster)

| Resource | Path |
|----------|------|
| ICON input/output | `/home/gilberto/data/ICON/` |
| MPAS input/output | `/home/gilberto/data/MPAS/` |
| WRF input/output | `/home/gilberto/data/WRF/` |
| Scripts | `/scratch/gilberto/scripts/` |
| Logs | `/scratch/gilberto/logs/` |

---

## References

See [Documentation.md](Documentation.md) for scientific references and detailed index equations.
