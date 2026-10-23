# Brown dwarf type against temperature (GEN.517)

Status: built 2026-10-23 by lane M (builder B9). Module `planetgen/physics/brown_dwarf_types.py`, data `physics/data/brown_dwarf_teff_map.json`, rows g1 to g5 of the atmosphere recipes (`scripts/atmosphere_recipes_data.py`), tests `src/tests/test_brown_dwarf_types.py`.

## 1. The finding

The type of a brown dwarf is not a weighted draw: it is a function of the effective temperature with a scatter (Kirkpatrick et al. 2021, Table 13, 525 L, T and Y dwarfs within 20 pc; READ). The rows g1 to g5 of part 09 are that map seen from the temperature, as soft windows with one softness of 100 K and edges that were design values (2200, 1300, 700, 500 K). The map gives the edges and the softness:

| Boundary | Spectral type | Temperature (K) | Scatter (K) | Design edge |
|---|---|---|---|---|
| late M and L | L0 | 2238 | 134 | 2200 |
| L and T | between L8.75 and T0 | 1273 | 79 | 1300 |
| T and late T (NH3-rich) | T8 | 686 | 79 | 700 |
| T and Y | Y0 | 460 | 79 | 500 |

The rows now use these edges and the scatter of each boundary (`phi((edge - teff_k)/scatter)` differences), so a dwarf at the temperature of a boundary is half one class and half the other, and the weights follow `brown_dwarf_types.type_probabilities`. The test pins the recipe to the module at twelve temperatures.

## 2. What it does not do

- The late-M rows (g1) have no fit in Table 13 (the map starts at L0); the 2238 K edge is the map's L0 and the 2800 K cap is the repo's maximum for a rogue body. The row stays partly verified with that limit.
- Row g4 is not a type: T8 and T9 are late T, NH3-rich, and the Y0 dwarfs start at 460 K (observed 371 to 510 K). The row stays partly verified; its edge at T8 is a design choice on the map.
- The map is for field dwarfs. Young low-gravity L dwarfs are up to 300 K cooler than field dwarfs of the same type (Filippazzo 2015), so a gravity term is needed below about 500 Myr; it is not modelled.
- The mass function (alpha 0.6 +/- 0.1) and the evolution tracks that give a dwarf its temperature are not part of this change; the census of Kirkpatrick Table 15 is not machine readable here (NOT REACHED), so no count-per-bin target exists.
- The recipes are not yet read by the generator, so this changes no generated value and takes no generator epoch.
