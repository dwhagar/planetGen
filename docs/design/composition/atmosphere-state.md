# The air's density, scale height and greenhouse warming from its mixture (GEN.423)

Module: `src/planetgen/physics/atmosphere_state.py`. Tests: `src/tests/test_atmosphere_state.py`. Lane L, 2026-10-24.

## What it does

Given the mixture (the dict of `physics/mixture.py`: gas-table formula to partial pressure in Pa), a temperature and a gravity:

| Quantity | Formula |
|---|---|
| mean molar mass | `M = sum x_i M_i` |
| surface density | `rho = P M / (R T)` |
| scale height | `H = R T / (M g)` (so `rho g H = P`) |
| column mass | `P / g` |
| mass of gas i above 1 m2 | `m_i = x_i M_i P / (M g)` (they add up to `P / g`) |
| infrared optical depth | `tau0 = sum kappa_i m_i`, `kappa_i` the grey mass absorption coefficient, m2/kg |
| warming | `T_s` from `physics/greenhouse.py` (Robinson and Catling) with `tau0`, `gamma` of the mixture (NASA fits of `species_thermo`), `alpha` and `n`; `G = (T_s/T_e)^4 - 1` |

A change of the mix therefore moves all of them: swapping N2 for CO2 raises `M` and `rho` and lowers `H`; adding water vapour adds
`kappa_H2O m_H2O` to `tau0` and raises `G`.

## Stand-ins and limits (replaced by GEN.703, lane E, when it lands)

- `kappa_i` are read from the research's statements, not tuned to any planet: Koll 2019 (arXiv:1907.13145, Fig. 4 text, part 03
  section 6.2) finds a pure H2O air thick above about 1e-2 bar (`kappa` 1e-2 m2/kg), a pure CO2 air of order unity at 1 bar (taken as
  3 optical depths, `kappa` 3e-4) and CO thin even at 1 bar (1e-5). Any other gas the table flags as an infrared absorber gets 1e-3
  (between CO2 and H2O), a "weak" one 1e-5, collision-induced-only gases and noble gases zero. All are stand-ins.
- No collision-induced absorption (dense N2, H2 or CO2 air at tens of bar is under-warmed), no condensation (the stored H2O is
  taken as it is), one grey band, the pressure exponent `n` = 1, the lapse scaling `alpha` = 0.8 for every body.
- Untuned results as a check, not as inputs: Earth air with 1 kPa of water vapour gives `tau0` 0.64 and `G` 0.41 (measured 0.625),
  a Venus-like 92 bar CO2 air gives `tau0` 304 and `G` 113 (99.6), a Mars-like 600 Pa air gives `G` 0.03.

## Wiring (not done here)

`physics/planets.py` still draws `atm_molar_density`, the atmosphere density and the greenhouse multiplier apart. Wiring this module in
is part of lane I's GEN.411 rewrite of that path (the mixture is drawn first, then the density, the scale height, the pressure and
the greenhouse come from it), which takes a generator epoch; `scripts/test.py` and GEN.703's tau0 check read the same functions.
