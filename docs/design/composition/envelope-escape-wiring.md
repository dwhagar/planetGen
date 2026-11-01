# The envelope escape of a close-in planet, wired to rows e1 and e2 (GEN.822)

Status: module and tuning, written on 2026-10-25 (as the lane reported) by the Citations lane. The module is `src/planetgen/physics/envelope_escape.py`; rows e1 and e2 of `scripts/atmosphere_recipes_data.py` call it through the recipe function `envelope_lost_fraction`; the valley survey of GEN.509 and GEN.722 (`radius_valley.py`) uses it too. It replaces the energy-limited loss and its efficiency knob of GEN.509 (`radius-valley-tuning.md`, sections 1 and 2, which this note supersedes). The mixture draw of 4c (GEN.411, GEN.422) reads these rows, so the weight of a sub-Neptune's bare-core and envelope options changes for every generated planet of the classes S, V, I and T: generator epoch 110, stage `fill`.

## 1. The answer in short

- Rows e1 (the planet is a bare core) and e2 (it keeps its envelope) now weigh the share of the starting envelope that the hydro-based mass-loss rate of [R2512] (`thermosphere.hydro_based_mass_loss_g_s`, GEN.738) has removed, integrated over the star's XUV history (`activity.coronal_xuv_at`) while the envelope thins and the planet shrinks with it (radius from [R2741], `bulk_structure.envelope_radius_earth`). The old row used `10 f_env0 (15 / v_esc)^2 pi R_core^3 E_xuv S / (G M^2)`, an energy-limited loss with a fitted efficiency; both the formula and the second knob are gone.
- **What the model says, and what had to be tuned.** Started at the 10 Myr of disc dispersal, the Kubyshkina rate is 35 to 60 times the energy-limited rate of 0.15 for a 4 Earth-mass planet at 0.09 AU (its escape parameter is 15, the boil-off regime), and it strips every core below about 15 Earth masses inside 30 days whatever the starting envelope (up to 0.5 of the mass). No valley comes out: the survey gives a slope of +0.07 to +0.14 and a valley at 2.1 to 2.3 Earth radii. The rate is inside the published grid (1 to 39 Earth masses, 1 to 10 Earth radii, 300 to 2,000 K, 1e26 to 5e30 erg/s) and reproduces the nine planets of the paper's Table 2, so the early loss is real in the model; but real sub-Neptunes of 4 to 6 Earth masses with a few percent of envelope exist at 10 days. The integration therefore starts at `ONSET_GY` = 1 Gyr, and what is lost earlier (the boil-off of the hot, inflated envelope, and the saturated first hundred million years) is part of the starting envelope `f_env0`: the number is what is left when the model begins. That is a tuned stand-in (tier T9, tweak level L3, registry row AM81), not a measurement of the onset.
- **The result.** With the onset at 1 Gyr and the core-mass width narrowed from 0.4 to 0.25 in ln (RV2 of the registry, section 6.10), the simulated Kepler surveys put both published valley numbers inside 2 sigma for every starting envelope `f_env0` from 0.04 to 0.2: slope -0.099 (target -0.11 plus or minus 0.02, z 0.56) and 1.84 Earth radii at 10 days (target 1.84 plus or minus 0.03, z -0.1 to 0.1) for 0.07, 0.1 and 0.2; at 0.04 the radius is 1.88 (z 1.2). The recipe bracket for `f_env0` is 0.05 to 0.2, log (it was 0.035 to 0.08). The old ridge `eta = 10 f_env0` is gone: there is one knob (`f_env0`, a free bracket) and one tuned onset.
- What the model gives for a planet: with a starting envelope of 10 percent, a 4 Earth-mass core has lost 35 percent of it at 10 times the Earth's insolation and 91 percent at 30 times; a 6 Earth-mass core 10 and 57 percent; a 10 Earth-mass core 2 and 13 percent, and 70 percent at 100 times. Nothing is lost at 1 times the Earth's insolation from a 5 Earth-mass core (under 5 percent). In the survey the share of planets under 1.8 Earth radii is 100 percent for periods of 2 to 8 days, 97 percent for 8 to 16, 60 percent for 16 to 32 and 9 percent for 32 to 64 days.

## 2. Sources, and what was read

| Source | Read | Used for |
|---|---|---|
| [R2512], already read for GEN.738 | the grid boundaries, Table 1 and Table 2, the comparison with the energy-limited rate | the rate and its range of validity; the statement that the energy-limited form is off by up to 1e8 at small escape parameters and by a factor of 10 to 20 at intermediate ones |
| [R3500], through the tuning target table (GEN.509) | the two numbers, -0.11 plus or minus 0.02 and 1.84 plus or minus 0.03 Earth radii at 10 d | the targets `e1-valley-slope` and `e1-valley-radius` |
| [R2741] | the radius law, through `bulk_structure` | the radius of the planet as its envelope thins |
| [R3797], through GEN.509 | the peak of the core mass, 3.72 Earth masses | the median of the survey's core masses (kept) |
| Recalled (tier 3) | gas discs last a few to ten Myr | the 10 Myr of the first try |

## 3. The model

- `envelope_escape.envelope_left_earth(core_mass, f_env0, insolation, host_mass, luminosity, age, teq)`: the envelope left, in Earth masses, 0 for a bare core. The step is the smaller of 8 percent of the age and 3 percent of the envelope over the rate, so a fast loss is followed; the integration stops when the envelope is under 0.2 percent of the planet's mass (AM82). The XUV flux at the planet is the star's coronal XUV at that age (`activity.coronal_xuv_at`, which integrates to 3.0e16 J/m2 at 1 AU over 5 Gyr, the number the old row held as a constant) over the distance, which comes from the insolation and the luminosity. The equilibrium temperature is the stored one, or 255 K times the fourth root of the insolation (AM84), held to the fit's 300 to 2,000 K.
- `lost_fraction` is `1 - left / start`; the recipe function `envelope_lost_fraction` is `lost_fraction`. Row e1's weight is `(1 - p(b4)) f_lost` and row e2's is `(1 - p(b4)) (1 - f_lost)`, as before in form. The rows still pass the planet's mass as the core mass, as the old row did.
- The survey (`radius_valley.planet`) draws the starting envelope log-normally around `f_env0` (width 0.6, held to 0.5 of the mass), integrates to 5 Gyr, and gives the planet the radius of its core or of core and envelope.

## 4. Tests

`src/tests/test_envelope_escape.py` (the share lost against insolation, core mass, age, envelope, the bare core, the equilibrium temperature, the agreement of the left envelope and the share), `src/tests/test_radius_valley.py` (both valley targets inside 2 sigma at `f_env0` 0.05, 0.1 and 0.2; a weaker or a stronger loss moves the valley out; the recipe bracket and the row formula) and `src/tests/test_atmosphere_recipes.py` (the e-row probabilities follow the lost share).

## 4a. The stellar-mass slopes (GEN.722)

The two targets added by GEN.722, the valley against insolation (alpha 0.069 +0.019 -0.023) and against stellar mass (beta -0.046 +0.125 -0.117; [R0354]), were run against the hydro-based loss: both pass inside 2 sigma at f_env0 0.05, 0.1 and 0.2 (8 surveys of 12000 planets each, `test_radius_valley.py`). The mechanism is still not chosen (RV1).

## 5. Limits and what is not done

- **The onset is a fit, not physics.** A start at 1 Gyr throws away the saturated phase, where most XUV is emitted. A physical treatment would be the boil-off and core-powered mass loss of the first hundred million years (which GEN.722 asks about) with the photoevaporation of the XUV rate after it. Until it exists, `ONSET_GY` and `f_env0` together stand for it, and the starting envelope is read as the one that is left at 1 Gyr.
- **The valley is matched, not explained.** Both numbers come out only with the onset and the core width tuned together (a sweep of the onset from 0.3 to 1 Gyr and the width from 0.25 to 0.55 found the pair; a smaller loss scale at 10 Myr, 0.01 to 0.3 of the rate, did not: slopes of -0.06 to -0.19 with valleys at 2.0 Earth radii or more). Other onsets would fit with other widths.
- **Core mass and envelope widths** are design stand-ins; the core width changed from 0.4 to 0.25 (RV2 of the registry, section 6.10); the envelope width 0.6 did not.
- **The stellar-mass slope of the valley (the e-link) and the choice of mechanism stay with GEN.722**; the share of stripped cores follows from the model and the survey (section 1), not from a fit to a published fraction.
- **The Kubyshkina fit is for a hydrogen-dominated envelope** and is held to its mass, radius and temperature grid; a planet colder than 300 K is given 300 K (it loses little in any case) and a planet of more than 10 Earth radii is outside the fit.
- **What changes in a generated planet.** The drawn option (bare core or envelope) and the air that follows. The rows still pass the planet's mass as the core mass, as before; whether a stored radius follows the option is not changed here.

## 6. Registry

Rows AM81 to AM84 of `docs/design/research/gap-filling-plan.md` (and the generated `stand-in-values.md`): the onset age (T9, L3), the bare-core cutoff, the starting-envelope cap, and the equilibrium-temperature factor. The core-mass width is row RV2 (section 6.10 of the plan, GEN.722), changed from 0.4 to 0.25 (T9, L3), and RV1 and RV4 are updated for the new loss.
