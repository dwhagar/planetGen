# Cloud albedo on the surface-air fixed point (GEN.974)

Status: design only. Nothing in the generator reads it yet. The build is lock work (generated output changes), waits for GEN.972 (the albedo already moves with the temperature there), and rides the same epoch if the train allows. It adds one pure module, `src/planetgen/physics/cloud_albedo.py`, and the wiring below.

## What is missing

The albedo of a wet body is drawn once (`ALBEDO_WET_RANGE`, 0.25 to 0.40) and held while the air is solved. GEN.972 lets the ice move it. The clouds do not move it: a world that warms keeps the cloud it was drawn with. Real clouds change with the climate, and that is the largest single spread in the sensitivity of the Earth (IPCC AR6, 2021, table 7.10, from memory, not opened: cloud feedback +0.42 W/m2 per K, very likely range -0.10 to +0.94).

The published range matters for the test the sweep holds (`test_climate_sweep.py`). At 288 K the fixed point responds with 0.40 K per W/m2 (TEST.509), its outgoing flux rising 2.5 W/m2 per K. With the ice-albedo gain (+0.35, GEN.972) and the cloud gain (+0.42) the response is 1 / (2.5 - 0.77) = 0.58, over the lower bound of 0.54. With the low ends of both published ranges (+0.10 and -0.10) it is 0.40: the test passes only if the cloud gain is near its central value, so the size of the gain, not just its sign, is the point of this item.

## What the sources say

* **Earth, the cold and moderate side.** Warming thins the low cloud and brightens less: the net low-cloud feedback of the CMIP models rose from 0.09 to 0.21 W/m2 per K between CMIP5 and CMIP6, almost entirely in the shortwave (Zelinka and others 2020, GRL 47, section 3.2, opened; the paper gives the total only in a figure, not in its text, and the abstract names the cause as falling extratropical low-cloud cover and albedo). The sign is a loss of reflection with warming, so the albedo falls as the temperature rises.
* **The hot side of a slow rotator.** A tidally locked ocean planet builds a thick substellar cloud deck as the light rises, and its albedo goes up with the flux: 0.42 at 1,000 W/m2, 0.47 at 1,400, 0.51 at 2,000 and 0.53 at 2,200 for an M-star spectrum ([R0914] L45, table 1 group 3, opened). A fast rotator in the same paper does the opposite: 0.21 at 1,600 W/m2 for the 6:1 case. The two signs belong to the circulation, not to the temperature, so the rotation is an input.
* **Water limits the cloud.** The aerosol model (`aerosols.py`, GEN.414) already puts a cloud where a condensable crosses saturation going up and bounds its depth by the vapour above the base (`cloud_depth_limit`). A world with no ocean has no water cloud (the sulfuric cloud of Venus is GEN.1019 and GEN.1020, wired elsewhere).

None of these gives a law for a generated body, so the module is a registered stand-in with the sources as its bounds (rows below).

## The law

The cloud term moves the Bond albedo of the drawn body by a bounded amount, in two regimes that add:

    A(T, S) = A_ice(T) + dA_amount(T) + dA_deck(S)

`A_ice` is the law of GEN.972 (the drawn albedo is its ice-free value). The two cloud terms:

1. **Amount, every wet body with an ocean:** `dA_amount = -beta (T - T_REF)`, clipped to plus or minus `DELTA_MAX`. Warmer than `T_REF` the low cloud thins and the albedo falls; cooler it thickens and rises. `beta` is set so the cloud feedback of the Earth at 288 K is the assessment's central +0.42 W/m2 per K: `beta = 0.42 x 4 / 1361 = 1.2e-3 per K` (the same arithmetic as the ice law's calibration to AR6 table 7.10). `T_REF = 288 K`, the Earth's mean, because the albedo range of a wet body is centred on the Earth's 0.30: a body at 288 K keeps exactly its drawn albedo, and a hot one is 0.04 darker at 320 K. `DELTA_MAX = 0.12`, the shortwave cloud effect of the Earth (about 47 W/m2 over 340, Hartmann's budget, from memory, not opened) less the share that cannot thin, so no body loses or gains more cloud than the Earth holds.
2. **Deck, a locked ocean body only:** `dA_deck = +gamma (S - S_REF)` for `S` above `S_REF`, with `gamma = (0.53 - 0.42) / (2200 - 1000) = 9.2e-5 per W/m2` from table 1 of Yang and others, `S_REF = 1000 W/m2`, clipped at `DELTA_DECK_MAX = 0.12`. A body that is not locked in 1:1 (spin from the tide model, GEN.842) has no deck; that is the sign of the fast rotators of the same paper, which the amount term already carries. The deck term takes the stellar spectrum of the Yang M-star runs as the one in hand and is a stand-in for the others.

Both terms vanish where there is no liquid ocean (`air.ocean_present` false) or no air; the amount term vanishes where the air has no water vapour to condense. A giant, a dry body and an airless body are unchanged, as in GEN.972.

## How it is wired

`Setting` gains `cloud_law` (the pure function of `cloud_albedo`, None for none) and `locked` (a bool the caller reads from the body's spin state; GEN.842 holds it). Each step of `surface_air.solve` reads `albedo = ice_law(T) + cloud_law(T, S, locked)` at the temperature of the step before, as GEN.972 does for the ice. Only the temperature is an unknown, so the loop gain is the sum of the two slopes: the amount term adds 0.0012 per K to the ice slope, small beside the ramp's, and the deck term has no temperature in it, so the iteration's damping stays as it is. `Air.albedo` carries the sum and the stored GEN.412 pair is the last step's, so the inversion reproduces the temperature.

## Checks to build

* The cloud feedback at 288 K for the Earth twin is `-beta x 1361 / 4 = +0.42 W/m2 per K`, inside AR6's -0.10 to +0.94 (a pure-module test, like the ice law's).
* The climate sweep of the Earth twin (TEST.509) responds above 0.54 once GEN.972 and GEN.974 are both in; the strict expected failure comes off with whichever of GEN.972, GEN.973 (closed) and GEN.974 lands last.
* The deck term reproduces the albedos of Yang's table 1 (0.42 at 1,000, 0.53 at 2,200 W/m2) at its two ends.
* A dry body, a giant and an airless body do not move (their temperatures and albedos are unchanged against the golden).
* The census (Bug 2) shows the albedo and temperature shift of wet bodies by temperature band; only bodies away from 288 K move, by at most `DELTA_MAX`.

## Stand-ins for the registry (GEN.725), with tier and tweak level

| id | what | tier | level | replaced by |
|---|---|---|---|---|
| CL1 | `beta` 1.2e-3 per K (the cloud feedback of +0.42 W/m2 per K at 288 K) | T6 educated guess from a published statistic (AR6, from memory) | L2 tuned to a published statistic | a cloud-resolving or at least a two-layer cloud model on the fixed point |
| CL2 | `T_REF` 288 K and `DELTA_MAX` 0.12 | T6 | L0 as envisioned for the reference, L2 for the bound | the same |
| CL3 | `gamma` 9.2e-5 per W/m2, `S_REF` 1,000 W/m2, `DELTA_DECK_MAX` 0.12 | T2 Confident Unverified (Yang and others 2013, table 1, M-star spectrum, one model) | L1 nudged inside the source's spread | the circulation (a 3D model, or its published scalings with rotation and flux) |
| CL4 | the cloud acts only with a liquid ocean and (for the deck) a 1:1 lock | T6 | L0 | an amount from the condensate depth of `aerosols.cloud_depth_limit` |

No row is at L3 or L4, so none needs a burn-down line.

## What this does not do

* No cloud longwave effect: the model's cloud is a reflector only. The AR6 figure is a net of shortwave and longwave, so sizing the albedo term to the net lets one term stand for both, which overstates the albedo change by the longwave share (unknown from the text opened; Zelinka's paper puts most of the CMIP6 increase in the shortwave). The longwave effect of clouds is part of the reason the model stands 6.5 K under the Earth's 288 K at the Earth's own flux (TEST.509), and belongs to a later item.
* No change to the cloud of the colour model: `aerosols.draw_aerosols` still draws `cloud_cover` and `tau_cloud` from the finished air for the sky (they do not feed back on the temperature here).
* No runaway: the water-window opacity that lets the outgoing flux plateau is GEN.1194, which may share the epoch. The deck term raises the albedo of a locked world with the flux, as the cloud does in Yang's runs, but the dayside's moist greenhouse is not modelled here.
