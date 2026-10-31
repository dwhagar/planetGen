# Climate branches: the ice-albedo feedback and the two climates of one flux (GEN.1046, with GEN.972)

Status: design and a pure module (`src/planetgen/physics/climate_branches.py`), lock-free. Nothing in the generator reads it yet. The wiring is GEN.972 (the feedback inside the surface-air fixed point) and GEN.1051 (W20 of `docs/design/systems/weird-worlds.md`: both branches handed back, one taken by a history draw).

## What the feedback does

A world with an ocean that cools freezes, and ice and snow reflect more of the light (Bond albedo 0.5 to 0.8 against 0.30 for open sea and cloud), so it cools further (Budyko 1969, Sellers 1969). Written as a loop, with `response(A)` the mean surface temperature the rest of the physics gives when the albedo is held at `A`:

    T = response(albedo(T))

For one flux the loop can have three answers: a warm world, a snowball, and an unstable one between them. A small push away from the middle one grows; a push away from either end dies out. The loop is stable where its slope `d response(albedo(T)) / dT` is under 1.

That is the whole of the physics this module holds. Everything else (the greenhouse, the vapour, the cloud) stays in `response`, so the three feedbacks of GEN.972 to GEN.974 add their terms there and this module does not change.

## The module

* `IceAlbedoLaw`: the albedo as a function of the mean surface temperature. Ice-free (0.30) above 283 K, covered (0.62) below 263 K, a straight ramp between. The ramp is wide because a global mean carries the spread from equator to poles; a one-dimensional ice-line model (North 1975) would replace it.
* `fixed_points(response, law)`: every root of `T - response(albedo(T))` between 150 and 700 K, found by a 1 K scan and refined with Brent's method, coldest first, each as a `Branch(temperature_k, albedo, ice_fraction, stable)`.
* `stable_branches`, `is_bistable`: the stable ones; two means a snowball and a warm world both exist.
* `bistable_window(response_at_flux, fluxes)`: the lowest and highest flux of a list at which a body is bistable. The caller builds `response_at_flux` from the body (the surface-air fixed point with its absorbed light scaled), so the window of an Earth twin can be compared with the 10 to 20 percent W20 asks for.

It reads no body and draws nothing. The history draw that takes one branch belongs to `oddities.pick_branch`, with its own sub-stream (weird-worlds section 4.4).

## What the grey toy shows, and what it does not

With Earth's flux, a held greenhouse factor of 0.63 (the 288 K Earth) and the ramp above, the loop has a snowball at 247 K, an unstable middle at 277 K and a warm world at 288 K, and it is bistable from 1,280 to 1,750 W/m2 (35 percent of the Earth's flux). That is wider than the 10 to 20 percent W20 asks for, for a known reason: a real snowball has almost no vapour, so its greenhouse falls and its branch is colder; and a real warm branch has more, so the window moves. With the vapour and cloud feedbacks inside `response` (GEN.973, GEN.974) the window narrows or widens by their gain. The ramp ends (row CC16) are the knob the wiring tunes to W20's width, so they carry tweak level L3.

## How the wiring uses it (GEN.972, GEN.1051)

1. GEN.972 puts `law.albedo(T)` into the fixed point of `surface_air.solve`: each step of the iteration reads the albedo at the temperature of the last step instead of the drawn one. Started from the no-air temperature of the drawn albedo, the iteration falls to the branch nearest its start. The drawn albedo then stays the ice-free value of the body (rock, ocean, cloud) and the ramp raises it with the ice.
2. GEN.1051 asks `fixed_points` for the branches of the body, hands both stable ones to `oddities.pick_branch`, which takes one by the history draw (the warm one the more often, by the share of the window the flux sits in), and stores the branch taken as the body's climate. A snowball also stores that it needs a build-up of carbon dioxide to leave (Hoffman 1998), which the habitability scoring reads.
3. The stored albedo and greenhouse factor are the pair of the taken branch, so the inversion of GEN.412 reproduces the temperature.

## Invariants (weird-worlds section 5)

* Each branch is a fixed point of its own loop (`test_every_branch_is_a_fixed_point_of_its_own_loop`): the albedo it reports is the one that gives its temperature, so the energy balance closes on every branch.
* A world with no ice sensitivity has one branch; a world far inside the cold or the hot side has one.
* The coldest stable branch is never warmer than the warmest.

## Stand-ins

Rows CC15 (the two albedos) and CC16 (the ramp ends) in `gap-filling-plan.md`, section 6.22. The replacement for both is an ice line in a one-dimensional energy balance, which would also give the partial ice cover of a body on the edge of freezing (the eyeball worlds of W18).
