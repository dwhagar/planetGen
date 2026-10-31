# The greenhouse table: the radiative-convective model in front of the fixed point (GEN.850, GEN.851)

Modules: `src/planetgen/physics/greenhouse_table.py` (the table), `physics/radiative_convective.py` and `physics/line_opacity.py` (the model it is fitted to,
`line-opacity.md`), `physics/aerosols.py` (the haze it reads), `physics/surface_air.py` (the caller). Fit: `scripts/fit_greenhouse_table.py`. File:
`physics/data/greenhouse_table.json`. Tests: `src/tests/test_greenhouse_table.py`, `test_haze_in_the_air.py`, `test_aerosols.py`,
`test_radiative_convective.py` and `test_radiative_convective_real.py`. Lane Air 1, 2026-10-30.

## What it does

`surface_air.solve` finds an air and its surface temperature together: it asks for the greenhouse of the air it has five to ten times a body. The model that
knows the answer (`radiative_convective.warming`: 40 layers, 130 bands by 10 g points of the ExoMolOP correlated-k tables of water, carbon dioxide and methane,
the HITRAN collision-induced pairs, the MT_CKD water continuum) takes one to twelve seconds a call and up to a few minutes for a dense air, so the
generator reads a table fitted to it:

    ln(T_s / T_e) = net(features)          T_e = (F / sigma)^(1/4), F the absorbed flux, G = (T_s / T_e)^4 - 1

`net` is a network of two tanh layers; its weights are the file. The grey closure of GEN.703 (`infrared_opacity.py`, Robinson and Catling's closure in
`greenhouse.py`) is deleted. The downloaded line data (1.4 GB) are still downloaded at setup, because the table is fitted to them and checked against them:
the file carries a digest of the dataset pieces it was fitted to (`data_digest`), and `greenhouse_table.table()` refuses it when they differ ("fit it again").

## The features

Sixteen numbers describe an air (`greenhouse_table.FEATURE_NAMES`): the logarithms of the effective temperature, the gravity and the surface pressure; of the mole
fractions of H2O, CO2, CH4, N2, O2, H2 and He (an absent gas sits at `FRACTION_FLOOR`); of the share of the other infrared-active gases (those the gas table flags
`ir`) and of the weak ones (`weak`), which the model gives stand-in bands (IR1, IR2); of the mean molar mass; the heat-capacity ratio at `T_e`; and the share of the
absorbed sunlight a haze takes aloft with the logarithm of the pressure it takes it at over the surface pressure (GEN.851 below). A feature outside the fitted range is
held at its edge, as the k-tables hold their edge beyond 100 K and 100 bar; the table is fitted up to a surface pressure of 1e8 Pa (1000 bar), above which the model
does not converge, and the output is capped at 4000 K.

## The fit

`scripts/fit_greenhouse_table.py` (steps `sample`, `designed`, `solve`, `fit`; the work folder is outside the repository):

- **The airs.** The greenhouse calls of the generator over seeds 1000 to 1800 (27,000 calls, 14,900 distinct airs, haze included since GEN.851), one perturbed copy of each (the
  pressure scaled by up to a factor of e, every gas by up to about 1.5, the gravity by 1.35, the flux by 2.2, and in one in five a haze share log-uniform from 0.003 to 0.9
  at a pressure between 3e-6 and a tenth of the surface's), and 8,000 designed airs of four families the generator seldom makes but the table must be right for
  (terrestrial, carbon-dioxide, cold methane-nitrogen, hydrogen-helium).
- **The solves.** The model on each air, in parallel and resumable; airs the model does not converge on are dropped.
- **The loss.** The error of `ln(T_s / T_e)` weighted by `1 / (0.004 + 0.03 y)^2` (a few percent of the warming, with a floor), so a thin air is fitted as closely in
  absolute terms and a thick one in relative; L-BFGS-B on numpy; the validation is on fresh seeds the fit never saw.
- **The file** holds the weights, the feature ranges, the digest and the measured error (`summary`). `tests/test_greenhouse_table.py::test_the_file_states_an_accuracy_the_fit_reached` holds the file to the
  accuracy it states.

Fitting again is needed when the model or its data change (the digest says so): `python scripts/fit_greenhouse_table.py sample|designed|solve|fit /path/to/work`.

## Accuracy

The network has three hidden layers (128, 128, 64 tanh units) and was fitted in two rounds: the first on the generator's airs and the designed families, the second adding airs the first table's own fixed point visits
(DAgger-style, so the table is right where it will be asked). Airs where the model "converged" to a surface colder than its effective temperature (a numerical artefact of dense hot airs) are dropped. The validation is
1,306 airs from seeds the fit never saw (`validation2`, taken under the fitted table); the error is the relative error of `T_s` and the error in kelvin:

| Surface temperature | Airs | Median error | 90th percentile | 99th percentile |
|---|---|---|---|---|
| 0 to 300 K | 1,000 | 0.2 K (0.20 %) | 5.6 K (3.3 %) | 20.5 K (9 %) |
| 300 to 500 K | 213 | 7.1 K (2.1 %) | 35.4 K (9.1 %) | 78.2 K (17 %) |
| 500 to 1000 K | 50 | 62.6 K (11 %) | 256 K (28 %) | 319 K (49 %) |
| over 1000 K | 43 | 39.3 K (3.5 %) | 156 K (13 %) | 270 K (23 %) |
| all | 1,306 | 0.54 K (0.41 %) | 16.8 K (5.3 %) | 170 K (18 %) |

The table was fitted again for GEN.1030 and GEN.1032, on the model with the methane pseudo-adiabat, the ortho-para heat capacity of cold hydrogen (IR14, IR15) and the sourced far wings of carbon dioxide (GEN.1019): about 22,000 training airs. Against the model it stands, in a clear-sky carbon dioxide and nitrogen air at Venus's gravity, 2 to 8 percent warm from 5 to 92 bar (336 K against 329 K at 5 bar and 160 W/m2; 20 bar and 300 W/m2 is the worst, 509 K against 475 K), on the Earth-like airs between 225 and 255 W/m2 within about 3 K, and 1 to 2 percent warm on the three giants at their measured effective temperatures. It is blind to the methane share of a hydrogen air: the model cools by 7 percent as methane goes from 0.5 to 4 percent of a Uranus-like air (condensation shallows the lapse rate), and the table by under 1, so a drawn ice giant with a few percent of methane is read about 5 percent warm. A further round with 3,000 designed one-bar hydrogen and methane airs did not teach it, and it moved the Earth-like airs 5 to 8 K cold (a fit of this size is sensitive to its seed there), so it was not shipped. The 500 to 1000 K row is dominated by a few dense steam and carbon dioxide airs at 10 to 200 bar; the table is weakest there.

The loss weights (an error floor of 0.003 on `ln(T_s/T_e)` and a relative term of 0.025) make the fit follow the cool airs the habitability and the life read; the hot dense airs (over 500 K: runaway greenhouses and
Venus-like airs) are the model's own weakest and the table's: a hot dense air is right to the first tens of kelvin, not to the degree. The weights are stated in `scripts/fit_greenhouse_table.py` (`NOISE_FLOOR`, `RELATIVE_ERROR`).

## The guards (stand-in IR11)

The network was fitted on airs that always hold an absorber and are at least 10 Pa thick, and left alone it gave a noble-gas air a greenhouse of 0.4. Two guards multiply `ln(T_s / T_e)`: a surface pressure under
`THIN_AIR_PA` (10 Pa) scales it by `P / 10 Pa` (an optically thin air warms with its column), and a share of gases that absorb or collide in the infrared (everything the gas table does not flag `none`) under
`ABSORBER_SHARE_MIN` (1 percent) scales it by the share over that. Both are registry row IR11 (tier T5, tweak L0).

## The cloud (stand-in IR12)

The model is a clear-sky one and a body's albedo is its all-sky one, so a world that keeps the Earth's clouds would lose their longwave warming (26.5 W/m2 of the 266 W/m2
the clear-sky Earth sends out, CERES EBAF, 0.10). The sky sends out `(1 - s) * clear`, `s = CLOUD_LONGWAVE_SHARE * h / 0.77` for a surface relative humidity `h` (at most 1.3
times), nil without water vapour. A first version put the cloud into the column as a grey layer; it hid the whole of the vapour's feedback (the doubling of carbon dioxide
went from +1.72 to +0.99 K), so the cloud is a share of the flux, which masks forcing and feedback alike by a tenth. The CERES share leaves the Earth twin of the fixed point at 284.6 K
(the clear-sky model's own bias of 3 K); 0.13 puts it at the measured 288 K (tier T9, tweak L3 in the registry). The cloud falls with the surface humidity, so a world whose vapour
is held while it warms loses cloud: the fixed point's vapour follows the temperature and does not.

## The haze (GEN.851, stand-in IR13)

The aerosols were made after the air, from its final temperature. They are part of it now. The three haze draws (K_zz, yield, grain radius) come from a stream of the body's
own (`aerosols.haze_stream_of`, seeded by the mass, radius and distance, which are fixed before the air), carried in the `Setting`; at each step `surface_air.solve` reads
`aerosols.haze_of` for the air it has: the depth at 550 nm from the CH4/CO2 gate, the diffusion-limited supply and the fall integral, the share of the absorbed sunlight the haze takes
(the tholin's absorption at each wavelength over the star's blackbody, on a slant path of 2) and the pressure where half of it is taken from the top, in the haze's mass profile. Those two
go to the table. `generate_aerosols` makes the same haze of the final air (the life stage may have changed the air since, and the stored depth then follows the air, not the temperature).
A hazy air is cooler at the ground than the same air clear (the anti-greenhouse of Titan): `test_haze_in_the_air.py`.

## What moves

Every planet and moon with an air: the surface temperature, and with it the class, the habitability and the life draw. 28 percent of the airs the generator makes shift by more than 5 K against the
grey closure, the dense ones by 60 to 190 K (the closure was too cold); the worlds with a haze are cooler. The Earth twin of the fixed point stays at 288 K. Airs above 1000 bar are held at the table's edge.

## Limits

- The table is as good as its training airs: an air of a kind the generator and the designed families never made is the network's extrapolation, held at the feature ranges' edge.
- The derivative of the temperature with the flux (the climate sensitivity of GEN.972 to GEN.974) is smooth but not exact: the table's slope error is a few hundredths of the slope.
- The line data are the ones of GEN.715; the stand-ins of the model (IR1 to IR13) are in the table with it.
