# A neutron star's or stellar black hole's age (GEN.553)

**Status (2026-10-21):** built (`physics/stellar_evolution.collapse_age_gy`, `tests/test_remnant_age.py`).

A neutron star or stellar-mass black hole used to draw its age uniform over 10 Gy
(`tuning.COMPACT_REMNANT_AGE_RANGE_GY`), whatever stars lived around it. It is now a remnant of its population: its age is
the time since the star that made it collapsed.

| Step | Rule |
|---|---|
| Progenitor | the IMF (`sample_imf_mass_sol`) over the kind's progenitor range, `tuning.NEUTRON_STAR_PROGENITOR_SOL` (8 to 30) or `BLACK_HOLE_PROGENITOR_SOL` (20 to 100), no lighter than the lightest star that can have collapsed inside the population's window (`youngest_progenitor_sol`) |
| Time to collapse | the end of the progenitor's last giant phase, `main_sequence_lifetime_gy * remnant_age_fraction` (`remnant_formation_time_gy`, 10 to 40 Myr for a neutron star's progenitor) |
| Birth age | as any star's: the population's window (`population_age_range_gy`, `age_window_gy`) and its age law (`galaxy.age_law`, GEN.546), no younger than the time to collapse |
| Remnant age | birth age less time to collapse, never under `COMPACT_REMNANT_AGE_RANGE_GY[0]` (10 kyr) |

The population is the one the sector fill gives the body from its place's mix (`FillContext.apply_population`), as it gives a
free-floating brown dwarf. A remnant made without a place has the whole disk's window, 0 to 10 Gy, so a standalone neutron star
is almost what it was. A bulge or halo remnant is 8 to 13 Gy old, which the uniform draw never gave.

## What changes downstream

The pulsar chance falls with age (`exp(-age / PULSAR_ACTIVE_LIFETIME_YR)`), so a young population's neutron stars pulse and an old
population's do not, and a magnetar is born only inside a magnetar's active life.

## What it does not do

- **The progenitor is not stored.** Only the age comes out of the draw; GEN.554 (derived counts) and the remnant row of the body
  draw keep the progenitor with the kind.
- **Intermediate-mass and supermassive black holes** keep the uniform age: no star collapsed to them (their own origin rows).
- **A remnant given an age** (the core of a supernova remnant, a scatter row's magnetar) keeps it.
- **The progenitor range overlap.** A star of 20 to 30 solar masses is in both ranges; the kind is the scatter's, the progenitor
  only gives the age.
- **A remnant of a multi-star tree** (`multistar._remnant`) already draws its progenitor inside the system's age; the same
  `youngest_progenitor_sol` serves both.
