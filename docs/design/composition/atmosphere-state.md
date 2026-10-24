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
| infrared optical depth | `tau0` of the spectrum on a wavenumber grid (`infrared_opacity.py`), at the surface temperature |
| warming | `T_s` from `physics/greenhouse.py` (Robinson and Catling) with `tau0`, `gamma` of the mixture (NASA fits of `species_thermo`), `alpha` and `n`; `G = (T_s/T_e)^4 - 1` |

A change of the mix therefore moves all of them: swapping N2 for CO2 raises `M` and `rho` and lowers `H`; adding water vapour adds
`kappa_H2O m_H2O` to `tau0` and raises `G`.

## Limits

The infrared optical depth comes from `physics/infrared_opacity.py` (GEN.703, lane E, design note `infrared-opacity.md`): a
spectrum on a 10 cm-1 grid from the gases' bands, their far wings, the water vapour continuum and the collision-induced pairs, taken
as the equivalent grey thickness at the surface temperature (the two are solved together). The band and wing opacities are
STAND-INS until real line data are in; the stored water vapour is taken as it is; the closure's pressure exponent is 2.
Untuned checks: Earth air with 1 kPa of well-mixed vapour gives `tau0` 1.81 and `G` 0.65 (1.82 and 0.625), a Venus-like 92 bar CO2 air
gives `G` 76 (99.6), a Mars-like 600 Pa air gives `G` 0.075 (0.10 to 0.21).

## Wiring (not done here)

`physics/planets.py` still draws `atm_molar_density`, the atmosphere density and the greenhouse multiplier apart. Wiring this module in
is part of lane I's GEN.411 rewrite of that path (the mixture is drawn first, then the density, the scale height, the pressure and
the greenhouse come from it), which takes a generator epoch; `scripts/test.py` and GEN.703's tau0 check read the same functions.
