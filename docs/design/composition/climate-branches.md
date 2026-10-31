# Climate branches: the ice-albedo feedback and the two climates of one flux (GEN.1046, with GEN.972)

Status: design and a pure module (`src/planetgen/physics/climate_branches.py`), lock-free. Nothing in the generator reads it yet. The wiring is GEN.972 (the feedback inside the surface-air fixed point) and GEN.1051 (W20 of `docs/design/systems/weird-worlds.md`: both branches handed back, one taken by a history draw).

## What the feedback does

A world with an ocean that cools freezes, and ice and snow reflect more of the light (Bond albedo 0.5 to 0.8 against 0.30 for open sea and cloud), so it cools further (Budyko 1969, Sellers 1969). Written as a loop, with `response(A)` the mean surface temperature the rest of the physics gives when the albedo is held at `A`:

    T = response(albedo(T))

For one flux the loop can have three answers: a warm world, a snowball, and an unstable one between them. A small push away from the middle one grows; a push away from either end dies out. The loop is stable where its slope `d response(albedo(T)) / dT` is under 1.

That is the whole of the physics this module holds. Everything else (the greenhouse, the vapour, the cloud) stays in `response`, so the three feedbacks of GEN.972 to GEN.974 add their terms there and this module does not change.

## The module

* `IceAlbedoLaw`: the albedo as a function of the mean surface temperature. Ice-free (0.30) above 305 K, covered (0.62) below 240 K, and between them the ice cover grows as the cube of the way down the ramp. The ramp is wide because a global mean carries the spread from equator to poles, and cubic because the ice stays near the poles until the world is well into the cold and then runs to the tropics, which a one-dimensional ice-line model (North 1975) gives and which would replace it. The ends and the exponent are set by one number from outside: the surface albedo feedback of the Earth at 288 K, 0.35 W/m2 per K very likely 0.10 to 0.60 (IPCC AR6 2021, table 7.10). The law gives 0.34. The first version, a straight ramp from 263 to 283 K, gave no feedback at 288 K at all, so an Earth twin could never feel its own ice.
* `fixed_points(response, law)`: every root of `T - response(albedo(T))` between 150 and 700 K, found by a 1 K scan and refined with Brent's method, coldest first, each as a `Branch(temperature_k, albedo, ice_fraction, stable)`.
* `stable_branches`, `is_bistable`: the stable ones; two means a snowball and a warm world both exist.
* `bistable_window(response_at_flux, fluxes)`: the lowest and highest flux of a list at which a body is bistable. The caller builds `response_at_flux` from the body (the surface-air fixed point with its absorbed light scaled), so the window of an Earth twin can be compared with the 10 to 20 percent W20 asks for.

It reads no body and draws nothing. The history draw that takes one branch belongs to `oddities.pick_branch`, with its own sub-stream (weird-worlds section 4.4).

## What the grey toy shows, and what it does not

With a held greenhouse factor of 0.63 (the 288 K Earth) and the law above, the Earth's own flux holds one climate, a warm world at 287 K with 2 percent of the ramp iced. A flux 15 percent lower (1,150 W/m2, a faint young Sun) holds three: a snowball at 237 K, an unstable middle at 243 K and a warm world at 271 K. The toy is bistable from 1,070 to 1,210 W/m2, 10 percent of the Earth's flux, inside the 10 to 20 percent W20 asks for and lying under the present flux, as the snowball Earth of Hoffman (1998) needed a fainter Sun. A real snowball has almost no vapour, so its greenhouse falls and its branch is colder, and the warm branch holds more; with the vapour and cloud feedbacks inside `response` (GEN.973, GEN.974) the window moves by their gain. The ramp ends and the exponent (row CC16) are the knob the wiring tunes if the real response puts the window outside W20's range, so they carry tweak level L2 (tuned to one published statistic, the AR6 feedback).

## How the wiring uses it (GEN.972, GEN.1051)

1. GEN.972 puts `law.albedo(T)` into the fixed point of `surface_air.solve`: each step of the iteration reads the albedo at the temperature of the last step instead of the drawn one. Started from the no-air temperature of the drawn albedo, the iteration falls to the branch nearest its start. The drawn albedo then stays the ice-free value of the body (rock, ocean, cloud) and the ramp raises it with the ice.
2. GEN.1051 asks `fixed_points` for the branches of the body, hands both stable ones to `oddities.pick_branch`, which takes one by the history draw (the warm one the more often, by the share of the window the flux sits in), and stores the branch taken as the body's climate. A snowball also stores that it needs a build-up of carbon dioxide to leave (Hoffman 1998), which the habitability scoring reads.
3. The stored albedo and greenhouse factor are the pair of the taken branch, so the inversion of GEN.412 reproduces the temperature.

## The wiring of GEN.972 in detail

Decisions, so the build is not a second design.

**Which bodies.** Only a body whose `Setting` holds surface water (`water_kg > 0`) and is not a giant gets a law. Dry rock keeps its drawn albedo, and a giant its own range (CC13). An airless body, or one whose air is lost, has no ice feedback in the fixed point (its temperature is the no-air one at its drawn albedo, as now).

**What is drawn.** Nothing new. The albedo drawn today (`planets.calculate_atmospheric_conditions`, one `draw.uniform` of `climate_draws.albedo_range`) becomes the ice-free albedo of the law, and the ice-covered albedo is `max(ICE_COVERED_ALBEDO, ice_free)`: the draw count does not change, so no other body's stream moves. The two regimes of CC4 are merged: a body with surface water draws from `ALBEDO_WET_RANGE` always (0.25 to 0.40), and the frozen range (0.45 to 0.70) and `ALBEDO_FROZEN_BELOW_K` leave the draw, because the ramp now makes a cold wet world bright by its temperature, which is what the regime was standing in for. CC4 is rewritten, not added to.

**In the fixed point.** `Setting` gains `starlight_w_m2` and `albedo_law` (None for no feedback). Each step reads the albedo at the temperature of the step before (`law.albedo(temperature)`), gets its absorbed flux from it (`greenhouse.absorbed_flux_w_m2`, floored at `MIN_ABSORBED_FLUX_W_M2` as now), adds the giant's own heat as now, and asks the table. The iteration is the existing damped one; the albedo is a function of the temperature, so it adds no new unknown, only a slope: where the ramp is steep (20 K for 0.32 of albedo) the loop gain can reach 1 and the damped step must not overshoot. The convergence test is the temperature's (`TOLERANCE_K`); the albedo follows from it. `Air` gains `albedo` (None without a law), the albedo of the last step.

**What is stored.** `climate_draws.record(planet, air.albedo, air.greenhouse_factor)`: the albedo of the last step and the factor of the same step, so the GEN.412 pair reproduces the temperature exactly (`test_every_body_stores_the_draws_its_temperature_is_made_of`). The drawn ice-free albedo is not stored on its own; the ice fraction is the law applied to the stored temperature (`IceAlbedoLaw.ice_fraction`), which the hydrosphere and the page read.

**Which branch.** The iteration starts at the no-air temperature of the ice-free albedo, so it falls to the warm branch when one exists and to the snowball when the world has only that (GEN.1051 adds the second start and the draw). A body that finds no fixed point inside the step budget keeps the last step and says `converged=False` as now.

**What moves.** Every wet body whose temperature falls inside 263 to 283 K, or whose old regime was frozen: colder bodies get brighter (and so colder, to a snowball), warmer ones are unchanged. Dry bodies, giants and airless bodies do not move. The census (Bug 2) should show the share of wet worlds with ice cover and the share that flip.

**Tests.** The Earth twin of the climate sweep stays warm at its flux (288 K is above the ramp) and its response to a 10 percent cut of the flux rises (the ice-albedo term of the sensitivity; GEN.973 and GEN.974 add the others, and whichever of the three lands last removes the strict expected failure). A wet world at 0.9 of the twin's flux falls to a snowball. A dry twin does not change. Replays of the same body give the same branch.

**Order.** GEN.972 waits for GEN.1030 and GEN.1032 (the table it iterates is refitted there), then GEN.973, then GEN.974. It is lock work (generated output changes) and takes an epoch with the vapour feedback if the train allows, so the fixed point is retuned once.

## Invariants (weird-worlds section 5)

* Each branch is a fixed point of its own loop (`test_every_branch_is_a_fixed_point_of_its_own_loop`): the albedo it reports is the one that gives its temperature, so the energy balance closes on every branch.
* A world with no ice sensitivity has one branch; a world far inside the cold or the hot side has one.
* The coldest stable branch is never warmer than the warmest.

## Stand-ins

Rows CC15 (the two albedos) and CC16 (the ramp ends and the exponent) in `gap-filling-plan.md`, section 6.22. The replacement for both is an ice line in a one-dimensional energy balance, which would also give the partial ice cover of a body on the edge of freezing (the eyeball worlds of W18).
