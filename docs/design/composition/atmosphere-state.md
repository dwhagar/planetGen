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
| warming | `T_s` from the greenhouse table (`physics/greenhouse_table.py`, GEN.850, design note `greenhouse-table.md`): the surface temperature the radiative-convective model gives for the mixture, the gravity and the absorbed flux (and a haze's share of it taken aloft); `G = (T_s/T_e)^4 - 1` |

A change of the mix therefore moves all of them: swapping N2 for CO2 raises `M` and `rho` and lowers `H`; adding water vapour raises `G`.

## Limits

The warming is the fitted table of `physics/radiative_convective.py` (GEN.850), which reads the downloaded ExoMolOP and HITRAN line data (GEN.715):
its accuracy, range and stand-ins are those of `greenhouse-table.md` and `line-opacity.md`. The stored water vapour is the surface's, thinned with
height by Manabe and Wetherald's humidity profile. The grey closure of GEN.703 (`infrared_opacity.py`, Robinson and Catling's closure in
`greenhouse.py`) is gone.

## Wiring

`surface_air.solve` calls `greenhouse_warming` at each step of the fixed point (GEN.850), with the haze's share of the sunlight and its level when the
air hazes (GEN.851); `physics/planets.py` stores the answer.
