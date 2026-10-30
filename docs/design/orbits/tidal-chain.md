# The tide chain: rheology, spin torques and moon orbits (GEN.842, GEN.954 to GEN.960, GEN.778)

Status: partly built. Stage 1a (the rheology module, GEN.954) is built and called by nothing; every other stage is a plan. Plan and census: `/mnt/project-files/reports/gen842-tides-2026-10-29.md` (sections 4 and 9) and `/mnt/project-files/reports/gen842-spin-down-2026-10-29.md`. Owner: Tides 1. Boss, 8:16 PM EDT 2026-10-29: build all of it in 8.1; secondary sources are accepted by default; registered stand-ins are fine (GEN.725).

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

## What GEN.955 does with it

`tidal_heating.response` takes the mantle temperature, melt and frequency; `solid_power_w` and `molten_power_w` follow; the melt-limited state is found by balancing the tide against the heat pipe's loss in the box instead of chosen (TH9); `interior_heat.update_heat` fills `tidal_flux_w_m2` (GEN.778); `orbit_elements.circularisation_time_yr` takes the new Q. Rows TH3, TH4, TH5 and TH7 retire.
