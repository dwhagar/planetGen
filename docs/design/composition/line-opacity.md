# The infrared optical depth of an air, band by band (GEN.715)

Modules: `src/planetgen/physics/line_opacity.py` (the optical depth of a layer), `src/planetgen/physics/radiative_convective.py` (the column that
uses it), `src/planetgen/datasets/ir_ktables.py` and `src/planetgen/datasets/hitran_cia.py` (the downloaded data). Tests:
`test_ir_datasets.py`, `test_line_opacity.py`, `test_radiative_convective.py` (all offline, on fixtures built by
`ir_datasets_support.py`) and `test_radiative_convective_real.py` (against measured numbers, on the real files when the machine has them).
Math lane, 2026-10-24. Research: part 03 of `atmosphere-composition-parts/` (sections 2, 3, 5.4, 6.2).

## What it replaces

GEN.703 gave an air one grey optical depth from band centres whose strengths were stand-ins. It could not say what a 15 micron band does when
the water vapour thins aloft, put Titan twice too warm, and warmed Earth by a third to a half of what a doubling of carbon dioxide is measured to.
Boss's decision of 4:01 PM EDT on 2026-10-24 was to use the public line-list data, downloaded at setup and reduced to small tables. This is that.

## The data (downloaded by `update.sh`, about 1.4 GB once, then 6.5 MB kept)

| Dataset | Fill | What | Source |
|---|---|---|---|
| `ir_ktables` | AD10 | Correlated-k tables of H2O, CO2 and CH4 at resolving power 1000, 0.3 to 50 micron, 22 pressures (1e-5 to 100 bar), 27 temperatures (100 to 3400 K), 20 g points | ExoMolOP [R0795] from POKAZATEL, UCL-4000 and YT34to10 |
| `ir_cia` | AD11 | Collision-induced absorption of 11 pairs: N2-N2, O2-O2, N2-O2, CO2-CO2, N2-CH4, N2-H2, CO2-CH4, CO2-H2, H2-H2, H2-CH4, H2-He | HITRAN CIA [R2267] |

Setup downloads each file, checks its sha256, reduces it and deletes it (`transient` pieces). The reduction keeps the band mean of the k-distribution:
all the (wavenumber, g) values of a band are sorted by k, cut at the new g edges, and each interval's mean is the new k. The result is 129 bands (20 cm-1
wide near 10 to 5 micron, 50, 100 and 250 cm-1 wide further out) by 10 g points (six under g = 0.9, four over it so the strongest lines are resolved),
on the 3400 K tables cut at 2000 K. The tables are in the cache folder, never in the repository. **Generation requires them** (GEN.590): a missing file
makes `fetch.DataMissing` name the dataset and the setup command, and a run told `--disregard-missing-data` takes the labelled stand-ins instead.

## The optical depth of a layer

`line_opacity.layer_depth(prepared, temperature, pressure, water_fraction, molecules_cm2)` returns the depth of every layer on 130 bands (the first
is 0 to 200 cm-1, below the tables) by 10 g points:

1. **Line absorption** of H2O, CO2 and CH4: the table read at the layer's temperature (linear in 1/T) and pressure (linear in ln P), times the molecules of that gas.
2. **Random overlap** of gases in a band: all `10 x 10` sums with their joint weights, sorted again and cut at the g edges (`overlap`). Where one gas is thin (largest depth under 0.05) it adds only its mean.
3. **The water vapour continuum** (MT_CKD 4.3): self and foreign coefficients at eight wavenumbers between 200 and 3600 cm-1, temperature exponent on the self part.
4. **Collision-induced absorption**: `k(T) x1 x2 n N` for each pair of the air, k the band mean at the layer's temperature.
5. **Gases without a table** (stand-ins IR1, IR2): a plateau at each band of the gas table with exponential wings, sorted into a distribution.
6. **The far wings** (rows FW1 to FW7): the tables cut every line at 25 cm-1, and a Venus-like air needs what lies beyond.
   `physics/co2_far_wing.py` (GEN.1019, rows FW1 to FW5; it replaced the flat stand-in IR3): each carbon dioxide band's Lorentz tail beyond 25 cm-1 times the [R3490] sub-Lorentzian factor, out to 500 cm-1, with the band strengths taken from the carbon dioxide table at its lowest pressure.
   GEN.1020, wired into the column through `solve(..., cloud=(tau_550, radius_um))`: the acid cloud's infrared opacity (`cloud_infrared`, Mie droplets from the HITRAN sulfuric acid indices, rows CF1 to CF5) and the far wings of water, sulfur dioxide and carbon monoxide in a carbon dioxide air (rows FW6, FW7). With both, the clear-sky Venus-like test air reaches 757 K (measured 735 K); the wing alone gave 539 K.

## The column

`radiative_convective.solve(partial_pressures, gravity, absorbed_flux, aloft_share=0, aloft_pressure=None)` cuts the air into 40 layers equally spaced in ln P
from a ten-thousandth of the surface pressure, and finds the surface temperature at which

* the profile follows the adiabat of the mixture (heat capacity from `species_thermo` at each temperature; the water pseudo-adiabat once the vapour
  saturates, or in an air with no water the methane pseudo-adiabat, IR14; cold hydrogen at its ortho and para equilibrium, IR15) up to the tropopause,
  and is in radiative equilibrium above it, and
* the infrared flux to space equals the absorbed flux.

The flux is the two-stream form with a diffusivity factor of 1.66 and a source linear in optical depth between level temperatures, solved band by band and
g point by g point. The water vapour has Manabe and Wetherald's (1967) relative humidity, `RH = RH_s (p/p_s - 0.02) / 0.98`, with the surface humidity the
mixture's own vapour over the saturation pressure, and none above the coldest level below (the cold trap). `aloft_share` puts that share of the absorbed flux
at `aloft_pressure_pa` as a Gaussian in ln P instead of at the ground: the sunlight a haze takes (the anti-greenhouse).

GEN.1030: with the dry adiabat alone the ice giants came out 20 percent too warm at 1 bar. Methane condenses near 80 K and 1.3 bar there and the
condensation shallows the lapse rate, and the heat capacity of hydrogen below 200 K is the equilibrium one (4.4 R at 60 K against the 3.5 R the fit of
`species_thermo` holds from 200 K up), which shallows it again. At the measured effective temperatures the column gives Uranus 79 K (76 measured), Neptune 78 K
(72) and Jupiter 167 K (166); Saturn comes out 149 K against 134 K measured, a miss this change does not close. The check is
`test_the_giants_sit_near_their_measured_temperature_at_one_bar` (Jupiter, Uranus, Neptune).

A step of the adiabat that began unsaturated and crossed the saturation corner used to take the saturated slope for its whole length; the branch is now
fixed for a step and the step is cut at the crossing (`test_the_adiabat_is_continuous_and_rises_with_the_surface_temperature_across_the_saturation_corner`).

## Checks (not inputs)

| Case | Model | Measured or published |
|---|---|---|
| Doubling CO2 from 280 ppm, vapour held (GEN.717) | +1.04 K | about 1.2 K (3.7 W/m2 forcing over 3.2 W/m2/K) |
| Same, vapour free | +1.72 K | larger than held, as water vapour feedback is |
| Earth, F = 240 W/m2, 410 ppm, 1 kPa of vapour, with the cloud share (IR12) | about 284 K, tropopause 13 kPa (275.6 K and 20 kPa clear sky); 288 K at its own fixed point | 288 K with clouds; clear-sky models 275 to 280 K |
| Mars, 600 Pa, F = 110 W/m2 | 214.6 K | 210 to 215 K |
| Venus, 92 bar, F = 160 W/m2, sourced far wing and acid cloud (depth 25, 1 micron) | 757 K | 735 K |
| the same with the far wing alone | 539 K | 735 K |
| Titan, 1.47 bar, F = 2.61 W/m2, no haze (GEN.718) | 106.8 K | 94 K: a clear sky is too warm |
| Titan with 60 percent of the sunlight taken at 30 mbar | 93.3 K | 94 K |

GEN.717's miss (the old model gave 0.3 to 0.6 K for a doubling) was the missing line structure; GEN.718's (Titan twice too warm) was the missing absorption of sunlight
by the haze, which the old closure had no way to say. Earth's clear-sky target is rebased: the old `0.625` for G assumed a closure tuned to a cloudy 288 K.

## Limits (registry rows IR1 to IR13, `gap-filling-plan.md` 6.6b)

- The cloud is a share of the outgoing flux (IR12), not a layer, and there is no condensation but water. The humidity profile is Earth's. Sunlight is deposited by the caller, never absorbed by the gases here.
- The stratosphere is clear-sky radiative equilibrium without ozone heating. A gas without a table has stand-in bands.
- The tables hold 100 K as the lowest temperature and 100 bar as the highest pressure (edge value beyond); H2 and He broadening is under-represented for H2O and CO2.
- Venus passes only with the acid cloud's infrared opacity (the sourced far wing alone leaves it near 540 K; the droplet radius within the published modes spans about 620 to 820 K). Very dense or pure-steam airs take 20 to 30 seconds.
- Speed: one to several seconds a call; the generator reads the fitted table (GEN.850).
- A drier Earth (1000 Pa of vapour against 400 Pa) gives a noisier doubling (1.5 to 2.0 K held) because the tropopause moves one layer at a time (a follow-up).

## Wiring

`radiative_convective.warming` takes one to several seconds a call (up to half a minute for a dense steam or carbon dioxide air), too slow for the several
calls a planet's fixed point makes. GEN.850 puts a fitted table in front of it (`greenhouse_table.py`, `greenhouse-table.md`): the generator reads the table, the
model is what the table is fitted to, and `infrared_opacity.py` and the grey closure of `greenhouse.py` are deleted. This module's band stand-ins for gases
without a table (IR1, IR2) stay in `line_opacity.py` alone. The cloud of the Earth's longwave effect is IR12 and the haze's sunlight is IR8 and IR13.
