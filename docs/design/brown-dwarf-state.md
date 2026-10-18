# A brown dwarf as a star: its type, temperature and radius (GEN.252)

**Status (2026-10-18):** built (`physics/brown_dwarf.py`, `tests/test_brown_dwarf.py`).

A brown dwarf is too light to hold hydrogen fusion, so it only cools. It used to be a fixed 0.1 solar radius with no state
(a multi-star member) or a Jupiter-sized body with a log-uniform mass (the free-floating row). Now both have the type, the
temperature and the radius its mass and age give.

| Piece | Value | Source (read) |
|---|---|---|
| Luminosity | `L = 4e-5 Lsun (t / 1 Gy)^-1.3 (M / 0.05 Msun)^2.64`, the rogue planets' law | Burrows et al. 2001, Rev. Mod. Phys. 73, 719, eq. 1 (arXiv astro-ph/0103383) |
| Temperature | `Teff = 1550 K (t / 1 Gy)^-0.32 (M / 0.05 Msun)^0.83`, no hotter than 2,800 K (a young brown dwarf is larger and the old-age law overshoots) | the same paper, eq. 2; the cap is `tuning.ROGUE_MAX_EFFECTIVE_TEMPERATURE_K` |
| Radius | from the two: `R = sqrt(L / 4 pi sigma Teff^4)`, 0.088 solar radii at 0.05 solar masses, `M^-0.34`, `t^-0.01` (age left out) | derived; the paper's figure 3: older radii within 30% of Jupiter's |
| Mass function | `dN/dM ~ M^-0.6` between 13 and 80 Jupiter masses | Kirkpatrick et al. 2021, ApJS 253, 7 (arXiv 2011.11616): alpha 0.6 +- 0.1 from 525 L, T and Y dwarfs within 20 pc |
| Type from temperature | Teff of the type number (L0 = 0, T0 = 10, Y0 = 20) is `2237.5 - 144.96 x + 4.0301 x^2` to L8.75, `1437.9 - 18.309 x` to T4.75, `5141.3 - 368.65 x + 6.7301 x^2` beyond; inverted. L0 is 2,238 K, T0 1,255 K, Y0 461 K | the same paper, Table 13 (rms 134, 79, 79 K) |

Above L0 a brown dwarf is a late-M dwarf, named by the stellar scale (`stellar_evolution.spectral_letter_and_subclass`).
The paper's scatter is its measurement scatter, not a spread of the objects, so it is not drawn: a temperature has one type.

## Where it shows

- **A free-floating brown dwarf** (the `brown-dwarf` fill kind, 0.03 per pc^3, which is 0.21 per star against the census's
  "1 per 4.7") takes its mass from the census mass function and its radius from the laws. Its page names its type
  ("a free-floating T5 brown dwarf", header "Rogue Brown Dwarf, Type T5"). The type is read from the stored effective
  temperature, so there is no new column. It is its own system; its companion (15% of brown dwarfs are in pairs) is
  GEN.275's.
- **A brown dwarf member of a multi-star tree** (`multistar._star_leaf`, any member under 0.08 solar masses) has the type
  as its `state` (it was NULL) and the radius the laws give at the system's age (it was 0.1 solar radii). The members
  list on a scatter row's page shows "Brown dwarf (T5)".
- **The class page** `/classes/rogue-planet/brown-dwarf` lists the temperatures where L, T and Y begin.

The desert (a brown dwarf is kept inside 5 AU in 1% of draws) and the companion rows and rates stay as they were
(`multistar`, GEN.275).

## What it gives

Test targets from the census, none an input (`tests/test_brown_dwarf.py`): with a constant birth rate over 10 Gy, the
mass function and the laws, the shares of brown dwarfs in the 150 K bins from 600 to 1200 K are within a factor of 2
of Table 15's densities (2.80, 1.99, 1.72, 1.11 per 1e-3 pc^3 against 0.229, 0.208, 0.192, 0.143 of the sample).

## What it does not do

- **L dwarfs hotter than 1,500 K are 2 to 4 times too few** (`test_the_hot_l_dwarfs_reach_the_censuss_share`, a strict
  expected failure). The old-age power law is too cool for the heaviest brown dwarfs (a 0.07 solar mass one at 5 Gy is
  1,224 K here, so few reach L types hotter than 1,500 K, where the census counts 0.5 to 0.8 per 1e-3 pc^3 in every
  150 K bin). A track from Baraffe et al. 2003 or Phillips et al. 2020 (ATMO 2020), the grids the paper's own mass
  function fits use, would close it.
- **The L/T transition's slow cooling** (the density spike at 1,200 to 1,350 K) is only the fit's plateau; no clouds.
- **A young brown dwarf is bigger and cooler for its mass than the old-age laws say**; the 2,800 K cap stands in.
  Low gravity also lowers the type by up to 300 K (the census) and metallicity is solar.
- **The stellar mass function is not continued below 0.08 solar masses.** The census's slope joined to it at 0.08 would
  give 0.71 brown dwarfs per star against the observed 0.2; the free-floating rate already carries the observed one.
- **The hydrogen-burning edge** is a step: a 0.0764 solar mass brown dwarf (80 Jupiter masses) is a cooling body, a 0.08
  one a star of the main sequence relations.

## Code

- `physics/brown_dwarf.py`: mass function, temperature, radius, type.
- `tuning`: `BROWN_DWARF_*` with their sources, next to `ROGUE_BROWN_DWARF_*`.
- `generation/phenomena/rogue.py`: the free-floating brown dwarf's mass, radius, `spectral_type`, page text.
- `generation/multistar._star_leaf`; `web/scatter_members_view`; `web/lib/classref`.
- A change of the generator's output: epoch bumped (OPS.56).
