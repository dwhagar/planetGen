# The tide chain: rheology, spin torques and moon orbits (GEN.842, GEN.954 to GEN.960, GEN.778)

Status: partly built. Stage 1a (the rheology module, GEN.954) stage 2a (the torque-sum module, GEN.956) and stage 3a (the moon-orbit module, GEN.958) are built and called by nothing; every other stage is a plan. Plan and census: `/mnt/project-files/reports/gen842-tides-2026-10-29.md` (sections 4 and 9) and `/mnt/project-files/reports/gen842-spin-down-2026-10-29.md`. Owner: Tides 1. Boss, 8:16 PM EDT 2026-10-29: build all of it in 8.1; secondary sources are accepted by default; registered stand-ins are fine (GEN.725).

## What is built before this chain

The tide a star or planet raises in the body that orbits it is built (GEN.737 heat, GEN.731 despin and resonances, GEN.727 circularisation): `tidal_heating.response` gives each body one `k2` and one `Q`, and `tidal_spin.evolve` slows the spin by the host's tide alone. The chain adds the three things that leave out: a rheology that depends on temperature, melt and tide frequency; the other torques on a planet's spin (its moons', its air's); and the moons' own orbits under the planet's tide.

## The stages

One branch, a commit per stage, because the stages edit the same functions (`tidal_heating.response`, `tidal_spin.evolve`, `planets.generate_spin`, `planet.add_moons`). Each physics stage is a lock-free step (a pure module, tests, registry rows; merges freely) and a wiring step (changes generated output; rides train 5, spin and heat, after the star epoch GEN.932 and planet formation train 4).

| Stage | Lock-free | Wiring | What it gives |
|---|---|---|---|
| 1 rheology | GEN.954 `physics/tidal_rheology.py` (built) | GEN.955, with GEN.778 (fills `tidal_flux_w_m2`) | Q and k2 that follow the mantle's temperature and melt; retires the fitted `MOLTEN_K2_OVER_Q` and the chosen melt state |
| 2 spin torques | GEN.956 torque-sum module | GEN.957 (after the birth-spin retune GEN.953) | A planet's day from its star's tide, its moons' and its air's thermal tide, with an attractor choice |
| 3 moon orbits | GEN.958 moon-orbit module | GEN.959 | Outward drift, spiral-in of the moons of locked planets, rings from a fallen moon, the moons' tide as planet heat |
| ice | (rides GEN.954) | GEN.960 | Ice shell rheology for the hydrosphere and the icy branch of `tidal_heating` |

## GEN.954: the rheology module

`tidal_rheology.response(mass, radius, frequency, temperature, eta0, activation, melt)` returns `(k2, Q)`; `rock_response` and `ice_response` are the two ready-made forms. The compliance is Andrade (default) or Maxwell; the Love number is that of a homogeneous sphere with the complex rigidity (the correspondence principle), which reduces to `tidal_heating.rigid_k2` for a stiff mantle. Viscosity is the body's own drawn `eta0 exp(Ea / R T)` (`interior_heat.Interior`), lowered by melt; rigidity falls with melt. The module imports none of `tidal_heating`, `interior_heat`, `orbit_elements` or `tidal_spin`, which will import it.

Checked on (`src/tests/test_tidal_rheology.py`): the Moon (k2 0.022 against 0.0242 observed, Q 61 at 1,450 K against 38 to 60), the solid Earth (semi-diurnal Q a few hundred), Io (k2/Q a few hundredths at the heat-pipe melt share 0.3 with no fitted number, against 0.015 from its orbital acceleration, Lainey et al. 2009), a magma ocean (k2 to 3/2, no dissipation) and cold ice (Europa's 3e-4 to a factor of 2, a hundred times more near melting).

Census check (800 systems from seed 1000, 743 with planets, 8,822 rocky planets and moons; screening, the same caveat as the 10-29 report: the mantle temperature is the one the built model stores and the birth spin of a locked planet is taken as 3 h, of a locked moon 20 h; `/mnt/project-files/reports/tides1/census-q-of-t-2026-10-30.txt`): with Q from the rheology instead of 100, 18.9 percent of free rocky planets and 5.9 percent of free moons would be locked within their age, 0.1 percent of locked ones would stay free, 6.1 percent of planets and 3.4 percent of moons flip between a tide that is at least their own heat and one that is not, and the tide's power changes by more than a factor of 3 for 51 and 57 percent. The median Q comes out 155 (5th to 95th percentile 17 to 1,000). Without the cold ceiling (TR10) a quarter of planets had Q above 6e4 and magma-ocean planets never locked, so the ceiling is a finding as well as a stand-in.

Stand-ins: registry rows TR1 to TR10 (`docs/design/research/gap-filling-plan.md` section 6.23). The two numbers the Moon and Io tests do not pin are the Andrade timescale ratio (1) and the rigidity floor of a mush (1e-3).

## GEN.956: the torque-sum module

`physics/spin_torques.py` sums the torques on a planet's spin and runs the spin under them. A `Torques` holds the planet, its gravitational `Raiser`s (the star first, then moons), a `k2/Q` function (`constant_dissipation` for today's behaviour, `rheology_dissipation` for the GEN.954 rheology at the mantle's temperature) and an optional `ThermalTide`. Gravitational: `-F K_g b_g(sigma)` with `sigma = 2 (omega - n)`; thermal: the Debye torque of Leconte et al. 2015, `+F K_a q0 x / (1 + x^2)`, `x = (omega - n) / omega0`, with `q0` and `omega0` interpolated from that paper's Table 1 (`thermal_pressure_amplitude_pa`, `thermal_relaxation_rate`). `evolve` runs the one-dimensional flow `I d omega/dt = N` from the birth spin to the first attractor on its way: a resonance whose eccentricity window holds the orbit's, the synchronous state (a raiser's sign flip), a moon's corotation, or an asynchronous balance (`equilibria`: with constant Q at `x = (1 +- sqrt(1 - 4 c^2)) / (2 c)`, the outer stable, which is Leconte's equation 13). With one star and constant Q it reproduces `tidal_spin.evolve` exactly (a test over 60 random bodies).

Checked on (`src/tests/test_spin_torques.py`): the Earth-Moon day (2.0 ms per century from the lunar tide alone at k2/Q = 0.3/12 and the standard 3/2), the built rock Q (the Moon changes the Earth's day by under an hour in 4.5 Gy), Venus (thermal and gravitational torques the same size at Q = 100, opposite in sign), Mercury's 3:2 capture, the closed-form asynchronous equilibria, and a moon inside corotation that drags the planet's spin up. Stand-ins: TQ1 to TQ6 (plan section 6.24). One finding for GEN.957: with `tidal_rheology`'s cold ceiling the gravitational torque goes to a constant floor at low frequency, so the sign jump at synchronism returns for a cold planet and the synchronous state stays stable unless the thermal tide wins at the relaxation frequency, as in the paper's Andrade case.

## GEN.958: the moon-orbit module

`physics/moon_orbits.py` integrates the planet's spin and all its moons' distances together (LSODA, time in units of the age, state `ln a`). Each moon's tide puts the torque `N_i` of `spin_torques` on the planet and the opposite on its own orbital angular momentum, `d ln a / dt = -2 N_i / L_i` (signed `L`, so a retrograde moon always loses it); the star's tide and the air's thermal tide act on the spin. A moon reaching the fluid Roche limit (`rings.roche_limit_km`, for its own density) is removed: its mass is `ring_mass_kg` and its angular momentum goes to the planet's spin; a moon pushed past 0.4 of the Hill radius escapes. `planet_heat_j` is the energy the moons' tides dissipate in the planet, which GEN.959 reads as a forcing on `interior_heat`. Angular momentum is conserved to 1e-5 and the heat equals the loss of spin and orbital energy to 2 percent (tests).

Checked on (`src/tests/test_moon_orbit_evolution.py`): the Moon's recession (3.7 cm a year at k2/Q = 0.025 against 3.82 measured), Phobos (falls, 3 cm a year against 1.8), Triton (retrograde, falls in some Gy at Q = 1e4), a moon of a star-locked planet inside the synchronous radius (falls; the same moon round a fast Earth drifts out), a close moon of a slowly turning planet (falls into a ring, the spin rises) and several moons together. Stand-ins: TM1 to TM4 (plan section 6.25).

## What GEN.955 does with it

`tidal_heating.response` takes the mantle temperature, melt and frequency; `solid_power_w` and `molten_power_w` follow; the melt-limited state is found by balancing the tide against the heat pipe's loss in the box instead of chosen (TH9); `interior_heat.update_heat` fills `tidal_flux_w_m2` (GEN.778); `orbit_elements.circularisation_time_yr` takes the new Q. Rows TH3, TH4, TH5 and TH7 retire.
