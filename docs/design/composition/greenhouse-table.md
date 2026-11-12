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

Twenty-one numbers describe an air (`greenhouse_table.FEATURE_NAMES`): the logarithms of the effective temperature, the gravity and the surface pressure; of the mole
fractions of H2O, CO2, CH4, NH3, SO2, H2S, CO, HCN (GEN.852), N2, O2, H2 and He (an absent gas sits at `FRACTION_FLOOR`); of the share of the other infrared-active gases (those the gas table flags
`ir`) and of the weak ones (`weak`), which the model gives stand-in bands (IR1, IR2); of the mean molar mass; the heat-capacity ratio at `T_e`; and the share of the
absorbed sunlight a haze takes aloft with the logarithm of the pressure it takes it at over the surface pressure (GEN.851 below). A feature outside the fitted range is
held at its edge, as the k-tables hold their edge beyond 100 K and 100 bar; the table is fitted up to a surface pressure of 1e8 Pa (1000 bar), above which the model
does not converge. The model's surface stops at 4000 K, so the fit ends there; since GEN.1580 the table no longer clamps its output at 4000 K: above it the warming is the same
log ratio of the air and the flux, an extrapolation (registry IR18), and the surface of an air stops where its main gas is half ionized (`ionization.ceiling_k`, CC7).

## The fit

`scripts/fit_greenhouse_table.py` (steps `sample`, `designed`, `solve`, `fit`; the work folder is outside the repository):

- **The airs.** The greenhouse calls of the generator over seeds 1000 to 1800 (27,000 calls, 14,900 distinct airs, haze included since GEN.851), one perturbed copy of each (the
  pressure scaled by up to a factor of e, every gas by up to about 1.5, the gravity by 1.35, the flux by 2.2, and in one in five a haze share log-uniform from 0.003 to 0.9
  at a pressure between 3e-6 and a tenth of the surface's), and 8,000 designed airs of four families the generator seldom makes but the table must be right for
  (terrestrial, carbon-dioxide, cold methane-nitrogen, hydrogen-helium; since GEN.852 also cold ammonia, volcanic sulphur and reducing carbon-monoxide airs, and the hydrogen family with hydrogen sulphide and carbon monoxide).
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
| 0 to 300 K | 635 | 0.6 K (0.37 %) | 7.0 K (3.0 %) | 17.4 K (7.0 %) |
| 300 to 500 K | 274 | 6.5 K (1.9 %) | 32.3 K (8.2 %) | 117 K (28 %) |
| 500 to 1000 K | 73 | 20.5 K (3.5 %) | 98 K (15 %) | 188 K (32 %) |
| over 1000 K | 34 | 48.9 K (4.7 %) | 180 K (15 %) | 257 K (18 %) |
| all | 1,016 | 1.7 K (0.74 %) | 22.1 K (5.6 %) | 162 K (18 %) |

(GEN.852: the figures above are for the table fitted again to eight tabled gases and the line lists of GEN.852; the validation is 1,016 airs from seeds 5000 to 5100 the fit never saw, a larger share of them with the sulphur, carbon-monoxide and ammonia airs the older table did not carry. The table fitted before it read 0.54 K, 16.8 K and 170 K on 1,306 airs of the model before the line lists changed; the two sets are not the same airs, and the old table on the 570 first airs of this set read 5.7 K median and 52 K at the 90th percentile, because the model moved under it.) The third fit round added 1,805 airs around the airs the checks use (the Earth grid of 225 to 255 W/m2 and 600 to 4,000 Pa of vapour, a carbon-dioxide sweep from 100 to 2,400 ppm, Mars, Venus, a 5 bar carbon-dioxide air, Titan clear and under haze, a cold nitrogen-methane air, hydrogen-helium, and Jupiter, Uranus and Neptune from 0.03 to 10 bar, each three times over), because the first six fits, which saw only the generator's airs and the designed families, each missed one of those checks by 5 to 40 percent. Against the model on the checks the file stands within 2 K on the Earth-like grid, 0.2 K on Titan, 14 K (2 percent) warm on Venus and within 5 K on the three giants. A carbon-dioxide doubling gives +0.70 to +0.99 K in the table and +0.86 K in the model.

The table was fitted for GEN.1030 and GEN.1032, on the model with the methane pseudo-adiabat, the ortho-para heat capacity of cold hydrogen (IR14, IR15) and the sourced far wings of carbon dioxide (GEN.1019): about 22,000 training airs. Against the model it stands, in a clear-sky carbon dioxide and nitrogen air at Venus's gravity, 2 to 8 percent warm from 5 to 92 bar (336 K against 329 K at 5 bar and 160 W/m2; 20 bar and 300 W/m2 is the worst, 509 K against 475 K), on the Earth-like airs between 225 and 255 W/m2 within about 3 K, and 1 to 2 percent warm on the three giants at their measured effective temperatures. It is blind to the methane share of a hydrogen air: the model cools by 7 percent as methane goes from 0.5 to 4 percent of a Uranus-like air (condensation shallows the lapse rate), and the table by under 1, so a drawn ice giant with a few percent of methane is read about 5 percent warm. A further round with 3,000 designed one-bar hydrogen and methane airs did not teach it, and it moved the Earth-like airs 5 to 8 K cold (a fit of this size is sensitive to its seed there), so it was not shipped. The 500 to 1000 K row is dominated by a few dense steam and carbon dioxide airs at 10 to 200 bar; the table is weakest there.

The loss weights (an error floor of 0.003 on `ln(T_s/T_e)` and a relative term of 0.025) make the fit follow the cool airs the habitability and the life read; the hot dense airs (over 500 K: runaway greenhouses and
Venus-like airs) are the model's own weakest and the table's: a hot dense air is right to the first tens of kelvin, not to the degree. The weights are stated in `scripts/fit_greenhouse_table.py` (`NOISE_FLOOR`, `RELATIVE_ERROR`).

## The refit for the temperate wet branch (GEN.1304)

GEN.852 part 3 left the table 12 K cold on the branch an Earth-like world stands on: the Earth twin of the fixed point (1 bar nitrogen, 40 Pa of carbon dioxide, vapor at
0.77 of saturation, 239 W/m2 absorbed) sat at 276.0 K where the model puts it at 288.0 K, and the runaway sweep at 255 and 285 W/m2 read 282.8 K and 295.7 K against the model's 294.9 K
and 307.7 K, which is why `test_atmosphere_state.py` was red on main. The table sat under the model by 4.8 K at 400 Pa of vapor, 8.2 K at 1,200 Pa and 10.3 K at 2,000 Pa at 239 W/m2:
it underweighted the feedback of the vapor, and `test_greenhouse_table.py` did not see it because its Earth air (1,000 Pa of vapor) was one anchor in 22,000.

The work folder of the GEN.852 fit was not kept (hours of model time), so the table was not fitted from the beginning. `scripts/fit_greenhouse_table.py` gained the steps of a refit from the file:

- **`wet`** draws temperate nitrogen-oxygen airs with their vapor free (100 Pa to 3e4 Pa), 20 kPa to 400 kPa, gravity 4 to 16, 150 to 400 W/m2 absorbed, carbon dioxide from 10 ppm to 2 percent,
  a third of them near one bar and the Earth's gravity with 100 to 3,000 ppm (the airs the checks use), and a third with no oxygen or no argon. 10,000 airs, all solved by the model.
- **`solve --cool-below 340`** runs the model only on the airs the file puts under 340 K (about half a second to ten seconds an air), in a random order, so any prefix is a fair sample. 12,388 airs of the generator
  (seeds 1000 to 1300 and perturbed copies, 7,388 of them random and 5,000 drawn from the hazy and hydrogen-chloride-bearing ones, because the first 7,000 left those the weakest) were solved.
- **`fit --start --distill --hold anchors`** begins from the weights and the standardisation of the file, takes the model's answer for every solved air, takes the file's own answer for the hot airs
  the model was not asked about (over 340 K, minutes an air), and holds 2,500 copies of the airs the tests hold the table to (Mars, Venus from 5 to 92 bar, a 5 bar carbon dioxide air, Titan clear and
  hazy, a cold nitrogen-methane air, hydrogen-helium, the three giants; not solved, held at the old answers) so the correction does not move them. A temperate wet air that the generator made is not
  held at the old answer (it is in the fit only where the model solved it).

Result against the model: the Earth twin 287.7 K (model 288.0), 255 W/m2 at 294.5 K (294.9), 285 W/m2 at 307.0 K (307.7); at 239 W/m2 the air with 400, 1,200, 2,000 and 5,000 Pa of vapor reads 279.4, 287.3, 290.4 and
293.5 K (model 279.3, 287.3, 291.9, 294.0). On the 10,000 wet airs of the fit the median error is 0.5 K and the 99th percentile 3.4 K (the GEN.852 file: 3.8 K and 20 K, and 2.8 K cold on average).

The validation is 1,081 airs of the model's own answers from seeds 5000 to 5100 of today's generator, which the fit never saw (all of them solved by the model, hot and dense included):

| Surface temperature | Airs | Median error | 90th percentile | 99th percentile |
|---|---|---|---|---|
| 0 to 300 K | 773 | 0.7 K (0.45 %) | 10.5 K (6.1 %) | 88 K (35 %) |
| 300 to 500 K | 243 | 13.9 K (3.9 %) | 169 K (42 %) | 428 K (93 %) |
| 500 to 1000 K | 53 | 29.9 K (4.7 %) | 170 K (23 %) | 487 K (79 %) |
| over 1000 K | 12 | 40.2 K (3.6 %) | 132 K (11 %) | 1,564 K (66 %) |
| all | 1,081 | 2.0 K (0.97 %) | 45.7 K (12.6 %) | 271 K (69 %) |

The cool airs by class (`cool_by_class` in the file; the classes are read from the features: haze is any sunlight taken aloft, `other_infrared_gases` is 1 percent or more of gases that have only stand-in bands,
hydrogen chloride and fluoride among them):

| Cool airs (under 300 K) | Airs | Median error | 90th percentile | 99th percentile |
|---|---|---|---|---|
| clear | 262 | 0.6 K (0.40 %) | 5.3 K (2.6 %) | 17.8 K (11 %) |
| other infrared gases | 144 | 1.0 K (0.62 %) | 10.5 K (5.2 %) | 128 K (52 %) |
| haze | 367 | 0.6 K (0.48 %) | 16.1 K (9.0 %) | 92 K (35 %) |

The GEN.852 file on the same airs read, for the cool airs, 2.4 K (1.5 %) median and 32.7 K (16 %) at the 90th percentile; for the clear ones 1.6 K, 9.5 K and 39.7 K. **The accuracy the GEN.852 file stated
(0.37 % median, 3 % at the 90th percentile of the cool airs) was true of the airs the generator made when it was fitted. The generator now makes many airs it never saw**: hazy ones with up to 95 percent
of the sunlight taken aloft (47 percent of the cool airs of this set) and ones rich in hydrogen chloride and fluoride (19 percent held 1 percent or more of gases with only stand-in bands). The table is weaker there, and `test_the_file_states_an_accuracy_the_fit_reached`
now holds the clear-sky airs to the old figures and the other two classes to what they reached.

**Not changed, and open:** the airs over 300 K are held at the GEN.852 answers where the model was not asked (the hot dense airs: 13.9 K median and 169 K at the 90th percentile in 300 to 500 K), and the hazy and
hydrogen-chloride-bearing airs need several tens of thousands more model solves (the generator makes them at 5 to 10 seconds an air). A refit of those is a separate item.

## The refit for hazy and hot dense airs (GEN.1310)

GEN.1304 left two open: the airs the file put at 300 to 500 K (median 3.9 %, 90th percentile 42 %) were held at the GEN.852 answers where the model had not been asked, and the hazy airs read 9.0 % at the 90th percentile. GEN.1310 solved them with the radiative-convective model (`solve`, resumable, about 3 s an air per 4 cores): 6,012 hot dense airs the GEN.1304 file puts at 300 to 520 K, and 8,644 hazy ones, both drawn in a random order so any prefix is a fair sample. They were added to the GEN.1304 training set and the file was fitted again from the GEN.1304 weights (`fit --start --distill --hold anchors`, 12,000 iterations, 128-128-64). Hot and hazy airs the model did not finish in the time allowed were dropped from the fit; a first fit that kept them pinned to the old answer did less well on the 90th percentile of 300 to 500 K (33 % against 31 %) and of the cool hazy airs (7.4 % against 5.5 %), so the second fit leaves them out.

Validation (the same 1,081 airs of seeds 5000 to 5100; the GEN.1304 file in brackets):

| Surface temperature | Airs | Median error | 90th percentile | 99th percentile |
|---|---|---|---|---|
| 0 to 300 K | 773 | 0.6 K (0.41 %) [0.45 %] | 8.3 K (4.3 %) [6.1 %] | 25 K (13 %) [35 %] |
| 300 to 500 K | 243 | 7.5 K (2.3 %) [3.9 %] | 130 K (31 %) [42 %] | 426 K (93 %) [93 %] |
| 500 to 1000 K | 53 | 27 K (3.9 %) [4.7 %] | 246 K (27 %) [23 %] | 510 K (83 %) [79 %] |
| over 1000 K | 12 | 71 K (4.6 %) [3.6 %] | 136 K (10 %) [11 %] | 1,621 K (68 %) [66 %] |
| all | 1,081 | 1.4 K (0.69 %) [0.97 %] | 28.7 K (9.3 %) [12.6 %] | 296 K (66 %) [69 %] |

| Cool airs (under 300 K) | Airs | Median error | 90th percentile | 99th percentile |
|---|---|---|---|---|
| clear | 262 | 0.6 K (0.47 %) | 4.9 K (2.7 %) | 18 K (9 %) |
| other infrared gases | 144 | 0.7 K (0.40 %) | 5.3 K (2.6 %) [5.2 %] | 26 K (14 %) |
| haze | 367 | 0.5 K (0.37 %) | 12.3 K (5.5 %) [9.0 %] | 27 K (14 %) |

`test_greenhouse_table.py` passes unchanged: the anchors (Mars, Venus, 5 bar carbon dioxide, Titan, hydrogen-helium, the giants, the Earth grid) stay within their tolerances.

**Not done, and open:** the GEN.1189 ensemble of three seeds (members fitted from scratch, not from the shipped file) was not built; one longer fit shipped. The 300 to 500 K airs still have a 90th percentile of 31 % and a bias of 32 K cold: the worst are thick steam-rich and hazy airs the table overshoots (model 450 K, table 600 to 850 K), and the 99th percentile there did not move. The airs over 500 K have too few validation airs (65) to say more than "about as accurate as before". A model run that converges for the hot-hazy corner is the next step, not more fitting.

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

## Minor gases (GEN.1588)

The fit was made on airs whose minor gases are absent or at least a few parts in a hundred thousand of the air, so between the two the network invented a slope: ten parts per million of hydrogen cyanide added 76 K to Earth's air, a part in a hundred million of ozone 10 K, a trace of helium up to 100 K, and a trace of ammonia or hydrogen took the greenhouse away. `warming` now reads each gas other than water, carbon dioxide, methane, nitrogen and oxygen (`TRACE_EXEMPT`: the gases the fit saw at every share) that stands under `TRACE_FULL_SHARE` (1e-3 of the air) by its share over that share of what the table says it adds to the air without the minor gases, the optically thin limit, where the effect goes with the column. A gas at the share counts in full and a gas at zero not at all, so the answer is continuous where a gas crosses the share (only the small interaction between two traces steps). Registry row IR18 (tier T5, tweak L3, hand-chosen against the sweep). The better answer is a refit on airs with minor gases from 1e-9 up (`scripts/fit_greenhouse_table.py`); until then the table's own answer for a gas above 1e-3 stands, helium and hydrogen cyanide included, whose full-share effects are still larger than the physics gives.

## What moves

Every planet and moon with an air: the surface temperature, and with it the class, the habitability and the life draw. 28 percent of the airs the generator makes shift by more than 5 K against the
grey closure, the dense ones by 60 to 190 K (the closure was too cold); the worlds with a haze are cooler. The Earth twin of the fixed point stays at 288 K (287.7 K since GEN.1304; it read 276.0 K between GEN.852 part 3 and GEN.1304). Airs above 1000 bar (the pressure the table and the model stop at) rise past the table's top value by the adiabat, `ln(T_s / T_e)` gaining `(1 - 1/gamma) ln(P / P_top)` with the air's own heat-capacity ratio (GEN.1240, stand-in IR16): a convective air is as warm as its adiabat takes it from a radiating level the opacity sets and the column does not move, so a steam air of 1e5 bar reads 1,800 K where the held table read 700 K. It is not fitted to the model, which does not converge there, and the solver's own ceiling (2,000 K) still bounds it.

## Limits

- The table is as good as its training airs: an air of a kind the generator and the designed families never made is the network's extrapolation, held at the feature ranges' edge.
- The derivative of the temperature with the flux (the climate sensitivity of GEN.972 to GEN.974) is smooth but not exact: the table's slope error is a few hundredths of the slope.
- The line data are the ones of GEN.715; the stand-ins of the model (IR1 to IR13) are in the table with it.
