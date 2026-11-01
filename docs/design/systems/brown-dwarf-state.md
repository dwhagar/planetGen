# A brown dwarf as a star: its type, temperature and radius (GEN.252, GEN.549)

**Status (2026-10-30):** built (`physics/brown_dwarf.py`, `tests/test_brown_dwarf.py`); GEN.549 replaced the power laws with a cooling track, and GEN.948 reads the track from the whole Baraffe table instead of eleven thinned rows.

A brown dwarf is too light to hold hydrogen fusion, so it only cools. It used to be a fixed 0.1 solar radius with no state
(a multi-star member) or a Jupiter-sized body with a log-uniform mass (the free-floating row). Now both have the type, the
temperature and the radius its mass and age give.

| Piece | Value | Source (read) |
|---|---|---|
| Cooling track | Teff and radius for every mass of the table (0.01 to 1.4 solar masses) at its 30 ages from 0.5 Myr to 10 Gy, linear in log mass and log age (`datasets/baraffe_isochrones.py`, `cooling_track`; GEN.948 replaced `tuning.BROWN_DWARF_TRACKS`, eleven thinned rows) | [R0254], the BHAC15 isochrones (`BHAC15_iso.2mass`) |
| Past a track's last age | a track that reaches 10 Gy holds its last values (the hydrogen-burning plateau); a shorter one cools at `Teff ~ t^-0.32` with the last radius | [R0606], eq. 2 [R0606] |
| Temperature | the track's, no hotter than the temperature of the hydrogen-burning limit at the age (`substellar_ceiling.cap_k`: 2,950 K at the peak, 2,050 K at 1 Gy, 1,630 K at 5 Gy; a body at or above 0.07 solar masses only the 2,950 K peak) | [R0254] section 3.3 and the BHAC15 table (GEN.494) |
| Radius | the track's, 0.09 to 0.16 solar radii old, larger when young | as above |
| Luminosity | `L = R^2 (Teff / Teff_sun)^4` of the two | Stefan-Boltzmann |
| Mass function | `dN/dM ~ M^-0.6` between 13 Jupiter masses and 0.08 solar masses (83.8 Jupiter masses, where the stars' mass function begins; GEN.936: it ended at 80 and left 0.0764 to 0.08 empty) | [R2374]: alpha 0.6 +- 0.1 from 525 L, T and Y dwarfs within 20 pc |
| Type from temperature | Teff of the type number (L0 = 0, T0 = 10, Y0 = 20) is `2237.5 - 144.96 x + 4.0301 x^2` to L8.75, `1437.9 - 18.309 x` to T4.75, `5141.3 - 368.65 x + 6.7301 x^2` beyond; inverted. L0 is 2,238 K, T0 1,255 K, Y0 461 K | the same paper, Table 13 (rms 134, 79, 79 K) |

Above L0 a brown dwarf is a late-M dwarf, named by the stellar scale (`stellar_evolution.spectral_letter_and_subclass`).
The paper's scatter is its measurement scatter, not a spread of the objects, so it is not drawn: a temperature has one type.

## Where it shows

- **A free-floating brown dwarf** (the `brown-dwarf` fill kind, 0.03 per pc^3, which is 0.21 per star against the census's
  "1 per 4.7") takes its mass from the census mass function, draws its age as a star does, and takes its radius and temperature from the track at that age. The age (GEN.553) is the age window of the population the sector fill gives it from its place's mix, spread by that population's star formation rate (`galaxy.age_law`); made without a place it has the whole disk's window, 0 to 10 Gy. Its page names its type
  ("a free-floating T5 brown dwarf", header "Rogue Brown Dwarf, Type T5"). The type is read from the stored effective
  temperature, so there is no new column. It is its own system; its companion (15% of brown dwarfs are in pairs) is
  GEN.275's.
- **A brown dwarf member of a multi-star tree** (`multistar._star_leaf`, any member under the hydrogen-burning limit of its place's [Fe/H], 0.074 solar masses at solar Z up to 0.092 at none (GEN.936)) has the type
  as its `state` (it was NULL) and the radius the track gives at the system's age (it was 0.1 solar radii). The members
  list on a scatter row's page shows "Brown dwarf (T5)".
- **The class page** `/classes/rogue-planet/brown-dwarf` lists the temperatures where L, T and Y begin.

The desert (a brown dwarf is kept inside 5 AU in 1% of draws) and the companion rows and rates stay as they were
(`multistar`, GEN.275).

## What it gives

Test targets from the census, none an input (`tests/test_brown_dwarf.py`): with a constant birth rate over 10 Gy, the
mass function and the track, the shares of brown dwarfs in the 150 K bins from 600 to 1,950 K are each within a factor
of 2 of Table 15's densities (2.80, 1.99, 1.72, 1.11, 1.95, 0.94, 0.81, 0.78, 0.50, 0.72 per 1e-3 pc^3 for the bins from
600 to 1,950 K, the test covering 600 to 1,050 and 1,500 to 1,950 K). The hot L dwarfs are the heaviest brown dwarfs, which
cool to a plateau near 2,000 K at the hydrogen-burning limit and stay there for the age of the disk.

## What it does not do

- **The table is derived, not downloaded.** The track is a thinned copy of the model's values committed with the code,
  which the 2026-10-17 ruling allows (derived values may be committed; the model files stay out of the repository).
  A build-time dataset step (GEN.548) could replace it with the full grid; the functions would not change.
- **0.072 to 0.0764 solar masses.** In the model these objects are stable stars at the bottom of the main sequence, yet
  they sit inside the 80 Jupiter mass brown dwarf sampling range, so they stay here as L and late-M dwarfs.
- **The 1,200 K bin** (the L/T transition) comes out at 0.44 of the census's share, below the factor of 2 band the test
  covers; the model has no clouds, so the L/T transition's slow cooling is not there.
- **Low gravity** lowers the type by up to 300 K (the census) and metallicity is solar.
- **The stellar mass function is not continued below 0.08 solar masses.** The census's slope joined to it at 0.08 would
  give 0.71 brown dwarfs per star against the observed 0.2; the free-floating rate already carries the observed one.
- **The hydrogen-burning edge** is a step: a 0.0764 solar mass brown dwarf (80 Jupiter masses, the top of the sampling range, not the burning limit) is a cooling body, a 0.08
  one a star of the main sequence relations. Decided (Boss, contradiction sweep card, 7:29 AM EDT 2026-10-21): the burning limit itself is `hydrogen_burning_limit_sol(z)` (0.074 solar masses at solar Z), one function that feeds GEN.110's fusion check and the brown dwarf split.

## Code

- `physics/brown_dwarf.py`: mass function, track, temperature, radius, luminosity, type.
- `tuning`: `BROWN_DWARF_*` with their sources.
- `generation/phenomena/rogue.py`: the free-floating brown dwarf's mass, radius, `spectral_type`, page text.
- `generation/multistar._star_leaf`; `web/scatter_members_view`; `web/lib/classref`.
- A change of the generator's output: epoch bumped (OPS.56).


## The boundary with the stars (GEN.936)

The stars' mass function (`imf`) begins at 0.08 solar masses and the brown dwarf one ends there: one number, `mass_boundaries.STAR_BROWN_DWARF_SPLIT_SOL`,
so no mass between them is missing and none is counted twice. Whether an object of a mass burns hydrogen is a question of its
place: `stellar_abundances.hydrogen_burning_limit_for_fe_h` (0.074 solar masses at solar Z, 0.092 at none, linear in Z between;
[R0606] parametrisation). A primary drawn under it is drawn again (a system is built around a star, and the
generator makes none around a brown dwarf), and a companion under it is a brown dwarf member of its tree. The main sequence from
0.07 to 0.1 solar masses is the Baraffe table's (`main_sequence`), scaled to meet MIST at 0.1.
