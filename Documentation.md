# Scientific Documentation — Atmospheric Instability Indices

**Author:** Gilberto Ricardo Bonatti — Especialista em Modelagem Numérica

This document describes the physical basis, mathematical formulation, implementation details, and interpretation thresholds for each index computed by this system.

---

## Table of Contents

1. [CAPE — Convective Available Potential Energy](#1-cape)
2. [CIN — Convective Inhibition](#2-cin)
3. [K-Index](#3-k-index)
4. [Total Totals Index](#4-total-totals-index)
5. [Showalter Index](#5-showalter-index)
6. [Lifted Index](#6-lifted-index)
7. [BRN — Bulk Richardson Number](#7-bulk-richardson-number)
8. [Implementation Notes](#8-implementation-notes)
9. [References](#9-references)

---

## 1. CAPE

### Definition

CAPE is the positive buoyancy energy available to a lifted air parcel. It represents the maximum kinetic energy a rising parcel can acquire through buoyancy alone, and is directly related to the potential updraft speed of deep moist convection.

### Equation

$$\text{CAPE} = g \int_{z_{LFC}}^{z_{EL}} \frac{T_{v,\text{parcel}} - T_{v,\text{env}}}{T_{v,\text{env}}} \, dz \quad \text{[J kg}^{-1}\text{]}$$

Or equivalently in pressure coordinates:

$$\text{CAPE} = -R_d \int_{p_{LFC}}^{p_{EL}} \left(T_{v,\text{parcel}} - T_{v,\text{env}}\right) d\ln p \quad \text{[J kg}^{-1}\text{]}$$

Where:
- $g$ = gravitational acceleration (9.81 m s⁻²)
- $T_{v,\text{parcel}}$ = virtual temperature of the lifted parcel (K)
- $T_{v,\text{env}}$ = virtual temperature of the environment (K)
- $z_{LFC}$ = height of the Level of Free Convection (m)
- $z_{EL}$ = height of the Equilibrium Level (m)
- $R_d$ = gas constant for dry air (287.05 J kg⁻¹ K⁻¹)

### Implementation

CAPE is extracted directly from model output (`cape_mu` for ICON/MPAS, `cape_2d` for WRF). ICON and MPAS provide most-unstable CAPE; WRF uses `wrf.cape_2d()` with terrain-following option enabled.

### Interpretation

| CAPE (J/kg) | Convective Potential |
|-------------|----------------------|
| 0–500 | Weak |
| 500–1500 | Moderate |
| 1500–3000 | Large |
| > 3000 | Extreme |

---

## 2. CIN

### Definition

CIN is the negative buoyancy energy that must be overcome for a parcel to reach the LFC. It acts as a "cap" that inhibits or delays convection initiation.

### Equation

$$\text{CIN} = g \int_{z_0}^{z_{LFC}} \frac{T_{v,\text{parcel}} - T_{v,\text{env}}}{T_{v,\text{env}}} \, dz \quad \text{[J kg}^{-1}\text{]} \quad (\leq 0)$$

Or in pressure coordinates:

$$\text{CIN} = -R_d \int_{p_0}^{p_{LFC}} \left(T_{v,\text{parcel}} - T_{v,\text{env}}\right) d\ln p \quad \text{[J kg}^{-1}\text{]}$$

Where $z_0$ and $p_0$ denote the parcel's initial level.

### Implementation

CIN is extracted directly from model output (`cin_mu` for ICON/MPAS, `cape_2d[1]` for WRF). Values are negative (or zero).

### Interpretation

| CIN (J/kg) | Description |
|------------|-------------|
| 0 to -25 | Weak cap, convection easily triggered |
| -25 to -100 | Moderate inhibition |
| -100 to -200 | Strong cap, needs strong forcing |
| < -200 | Very strong cap, convection unlikely |

---

## 3. K-Index

### Definition

The K-Index (George 1960) is an empirical thermal-moisture index for assessing thunderstorm potential in tropical and sub-tropical environments. It accounts for atmospheric instability (temperature lapse rate) and moisture content at mid-levels.

### Equation

$$\text{KI} = (T_{850} - T_{500}) + T_{d850} - (T_{700} - T_{d700}) \quad \text{[°C]}$$

Where:
- $T_{850}$, $T_{700}$, $T_{500}$ = temperature at 850, 700, and 500 hPa (°C)
- $T_{d850}$ = dewpoint temperature at 850 hPa (°C)
- $T_{d700}$ = dewpoint temperature at 700 hPa (°C)
- $(T_{700} - T_{d700})$ = dewpoint depression at 700 hPa (°C), a measure of dryness

### Physical Interpretation

- **$(T_{850} - T_{500})$**: thermal lapse rate — larger values indicate greater instability
- **$T_{d850}$**: low-level moisture — higher dewpoints favor deep convection
- **$(T_{700} - T_{d700})$**: mid-level moisture deficit — small values indicate more moisture in the mid-troposphere

### Implementation

Computed from pressure-level temperatures and dewpoints. For ICON and MPAS, dewpoints at 850/700 hPa are derived from relative humidity via `metpy.calc.dewpoint_from_relative_humidity()`. For WRF, dewpoints are extracted via `wrf.getvar("td")` and interpolated to pressure levels.

```python
k_idx = (T850c - T500c) + Td850c - (T700c - Td700c)
```

### Interpretation

| KI (°C) | Thunderstorm Probability |
|---------|--------------------------|
| < 20 | None |
| 20–25 | Isolated (< 20%) |
| 26–30 | Widely scattered (20–40%) |
| 31–35 | Scattered (40–60%) |
| 36–40 | Numerous (60–80%) |
| > 40 | Very numerous (> 80%) |

---

## 4. Total Totals Index

### Definition

The Total Totals Index (Miller 1972) combines the Cross Totals (moisture contribution) and Vertical Totals (lapse rate contribution) into a single measure of convective potential.

### Equation

$$\text{TT} = \underbrace{(T_{850} - T_{500})}_{\text{Vertical Totals}} + \underbrace{(T_{d850} - T_{500})}_{\text{Cross Totals}} \quad \text{[°C]}$$

Which simplifies to:

$$\text{TT} = T_{850} + T_{d850} - 2 \cdot T_{500}$$

Where all temperatures are in °C.

### Component Breakdown

| Component | Equation | Physical Meaning |
|-----------|----------|-----------------|
| Vertical Totals (VT) | $T_{850} - T_{500}$ | Dry thermal instability |
| Cross Totals (CT) | $T_{d850} - T_{500}$ | Moisture contribution |
| Total Totals (TT) | VT + CT | Combined instability |

### Implementation

```python
tt_idx = (T850c - T500c) + (Td850c - T500c)
```

### Interpretation

| TT (°C) | Thunderstorm Potential |
|---------|------------------------|
| < 44 | None |
| 44–50 | Isolated, weak |
| 50–55 | Scattered, moderate |
| 55–60 | Numerous, some severe |
| > 60 | Severe thunderstorms possible |

---

## 5. Showalter Index

### Definition

The Showalter Index (Showalter 1953) assesses the instability of the mid-troposphere by lifting a parcel from 850 hPa to 500 hPa along the moist adiabat and comparing its temperature with the environmental temperature at 500 hPa.

### Equation

$$\text{SI} = T_{500,\text{env}} - T_{500,\text{parcel}} \quad \text{[°C or K]}$$

Where:
- $T_{500,\text{env}}$ = environmental temperature at 500 hPa
- $T_{500,\text{parcel}}$ = temperature of a parcel lifted from 850 hPa to 500 hPa along the moist adiabat, starting at the Lifted Condensation Level (LCL)

### Parcel Lifting Procedure

1. Compute the LCL from the 850 hPa parcel using:
   $$p_{LCL}, T_{LCL} = \text{LCL}(p_{850}, T_{850}, T_{d850})$$

2. Lift the parcel from $(p_{LCL}, T_{LCL})$ to 500 hPa along the saturated adiabat:
   $$\frac{dT}{d\ln p} = \frac{R_d T + L_v r_s}{c_{pd} + \frac{L_v^2 r_s}{R_v T^2}}$$

   where $L_v$ is the latent heat of vaporization and $r_s$ is the saturation mixing ratio.

3. The parcel temperature at 500 hPa gives $T_{500,\text{parcel}}$.

### Implementation

Uses `metpy.calc.lcl()` for the LCL and `metpy.calc.moist_lapse()` for the saturated adiabatic ascent, wrapped in a scalar function and vectorized:

```python
lcl_p, lcl_t = mpcalc.lcl(850 * units.hPa, T850, Td850)
T_parcel_500 = vec_moist_lapse(lcl_t.m, lcl_p.m)  # lifts to 500 hPa
SI = T500 - T_parcel_500
```

### Interpretation

| SI (°C) | Interpretation |
|---------|----------------|
| > 3 | Stable, convection unlikely |
| 1 to 3 | Slightly stable |
| -1 to 1 | Near-neutral, convection possible |
| -3 to -1 | Unstable, thunderstorms likely |
| -6 to -3 | Very unstable, severe storms possible |
| < -6 | Extremely unstable |

---

## 6. Lifted Index

### Definition

The Lifted Index (Galway 1956) is similar to the Showalter Index but uses a surface-based parcel. It measures the temperature excess of the environment over a parcel lifted from the surface to 500 hPa.

### Equation

$$\text{LI} = T_{500,\text{env}} - T_{500,\text{parcel}} \quad \text{[°C or K]}$$

Where $T_{500,\text{parcel}}$ is the temperature of a parcel lifted from the surface to 500 hPa along the moist adiabat.

### Parcel Lifting Procedure

1. Compute the LCL from the surface parcel:
   $$p_{LCL}, T_{LCL} = \text{LCL}(p_{\text{sfc}}, T_{\text{sfc}}, T_{d,\text{sfc}})$$

2. Lift the parcel from the LCL to 500 hPa along the saturated adiabat (same moist lapse rate equation as in SI).

### Difference from Showalter Index

| Feature | Showalter | Lifted |
|---------|-----------|--------|
| Parcel origin | 850 hPa | Surface |
| Sensitivity | Mid-level instability | Near-surface instability |
| Use case | General convection | Surface-based convection |

### Implementation

For ICON and MPAS:
```python
lcl_p, lcl_t = mpcalc.lcl(p_sfc, T_sfc, Td_sfc)
T_parcel_500 = vec_moist_lapse(lcl_t.m, lcl_p.m)
LI = T500 - T_parcel_500
```

For WRF (approximation using Bolton's formula):
```python
LI = (T500 - 273.15) - ((T_sfc - 273.15) - ((T_sfc - Td_sfc) * 0.4 + 22))
```

### Interpretation

| LI (°C) | Stability |
|---------|-----------|
| > 2 | Stable |
| 0 to 2 | Marginally stable |
| -2 to 0 | Slightly unstable |
| -4 to -2 | Unstable, thunderstorms likely |
| -6 to -4 | Very unstable, severe potential |
| < -6 | Extremely unstable |

---

## 7. Bulk Richardson Number

### Definition

The Bulk Richardson Number (Weisman & Klemp 1982) is a dimensionless parameter that relates the buoyancy energy available for convection (CAPE) to the kinetic energy of the vertical wind shear. It is used to distinguish between different convective storm modes (supercells vs. multicell vs. ordinary).

### Equation

$$\text{BRN} = \frac{\text{CAPE}}{\frac{1}{2} \|\Delta \mathbf{V}\|^2}$$

Where:
$$\|\Delta \mathbf{V}\|^2 = (u_{6\text{km}} - u_{\text{sfc}})^2 + (v_{6\text{km}} - v_{\text{sfc}})^2 \quad \text{[m}^2 \text{s}^{-2}\text{]}$$

- CAPE is in J/kg = m² s⁻²
- $u_{6\text{km}}$, $v_{6\text{km}}$ = wind components at ~6 km AGL (m/s)
- $u_{\text{sfc}}$, $v_{\text{sfc}}$ = surface wind components (m/s)
- BRN is dimensionless

### Implementation

```python
u_shear = u_6km - u_sfc
v_shear = v_6km - v_sfc
shear_sq = u_shear**2 + v_shear**2

BRN = CAPE / (0.5 * shear_sq)  # where shear_sq > 0.1 m²/s²
```

A threshold of `shear_sq > 0.1` avoids division by near-zero wind shear. WRF applies an additional clip `BRN = clip(BRN, 0, 500)` to remove physically unreasonable values.

**ICON:** Shear computed from 500 hPa wind vs. 10 m wind (as a proxy for the 0–6 km layer).
**MPAS/WRF:** Shear computed from wind at 6 km AGL vs. surface wind.

### Interpretation

| BRN | Storm Mode |
|-----|------------|
| < 10 | Supercell favored (strong shear, moderate CAPE) |
| 10–45 | Supercell likely |
| 45–50 | Transition zone (supercell to multicell) |
| 50–100 | Multicell storms likely |
| > 100 | Weak shear, disorganized or pulse convection |

---

## 8. Implementation Notes

### LCL Computation

The Lifted Condensation Level is computed analytically using Bolton's (1980) approximation, as implemented in MetPy:

$$T_{LCL} = \frac{1}{\frac{1}{T_d - 56} + \frac{\ln(T/T_d)}{800}} + 56 \quad \text{[K]}$$

$$p_{LCL} = p \left(\frac{T_{LCL}}{T}\right)^{c_{pd}/R_d}$$

### Moist Adiabatic Lapse Rate

MetPy's `moist_lapse()` integrates the pseudoadiabatic lapse rate numerically:

$$\Gamma_s = g \frac{1 + \frac{L_v r_s}{R_d T}}{c_{pd} + \frac{L_v^2 r_s}{R_v T^2}}$$

Where:
- $L_v$ = 2.501 × 10⁶ J kg⁻¹ (latent heat of vaporization)
- $r_s$ = saturation mixing ratio (kg/kg)
- $R_v$ = 461.5 J kg⁻¹ K⁻¹ (gas constant for water vapor)
- $c_{pd}$ = 1005.7 J kg⁻¹ K⁻¹ (specific heat of dry air)

### NaN Handling

The scalar moist lapse function includes explicit guards:
- Returns `NaN` if input temperature or pressure is `NaN`
- Returns `NaN` if starting pressure < 500 hPa (parcel already above target level)
- Any exception during MetPy computation returns `NaN`

### Output Compression

All output NetCDF files use `zlib` compression at level 4, which provides a good balance between compression ratio and I/O speed for float32 spatial fields.

---

## 9. References

- **Bolton, D. (1980).** The Computation of Equivalent Potential Temperature. *Monthly Weather Review*, 108, 1046–1053.

- **Galway, J. G. (1956).** The Lifted Index as a Predictor of Latent Instability. *Bulletin of the American Meteorological Society*, 37(10), 528–529.

- **George, J. J. (1960).** *Weather Forecasting for Aeronautics*. Academic Press, New York.

- **Miller, R. C. (1972).** *Notes on Analysis and Severe-Storm Forecasting Procedures of the Air Force Global Weather Central*. Tech. Report 200 (Rev.), Air Weather Service, Scott AFB, IL.

- **Showalter, A. K. (1953).** A Stability Index for Thunderstorm Forecasting. *Bulletin of the American Meteorological Society*, 34(6), 250–252.

- **Weisman, M. L., & Klemp, J. B. (1982).** The Dependence of Numerically Simulated Convective Storms on Vertical Wind Shear and Buoyancy. *Monthly Weather Review*, 110(6), 504–520.

- **Skamarock, W. C., et al. (2019).** *A Description of the Advanced Research WRF Model Version 4*. NCAR Technical Note NCAR/TN-556+STR.

- **Zängl, G., et al. (2015).** The ICON (ICOsahedral Non-hydrostatic) modelling framework of DWD and MPI-M: Description of the non-hydrostatic dynamical core. *Quarterly Journal of the Royal Meteorological Society*, 141, 563–579.

- **Skamarock, W. C., et al. (2012).** A Multiscale Nonhydrostatic Atmospheric Model Using Centroidal Voronoi Tesselations and C-Grid Staggering. *Monthly Weather Review*, 140, 3090–3105.

- **May, R. M., et al. (2022).** MetPy: A Meteorological Python Library for Data Analysis and Visualization. *Bulletin of the American Meteorological Society*, 103(10), E2273–E2284.
