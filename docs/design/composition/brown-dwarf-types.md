# Brown dwarf type against temperature (GEN.517)

Status: built 2026-10-23 by lane M (builder B9). Module `planetgen/physics/brown_dwarf_types.py`, data `physics/data/brown_dwarf_teff_map.json`, rows g1 to g5 of the atmosphere recipes (`scripts/atmosphere_recipes_data.py`), tests `src/tests/test_brown_dwarf_types.py`.

## 1. The finding

The type of a brown dwarf is not a weighted draw: it is a function of the effective temperature with a scatter ([R2374: Table 13], 525 L, T and Y dwarfs within 20 pc; READ). The rows g1 to g5 of part 09 are that map seen from the temperature, as soft windows with one softness of 100 K and edges that were design values (2200, 1300, 700, 500 K). The map gives the edges and the softness:

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
- The map is for field dwarfs. Young low-gravity L dwarfs are up to 300 K cooler than field dwarfs of the same type [R1370], so a gravity term is needed below about 500 Myr; it is not modelled.
- The mass function (alpha 0.6 +/- 0.1) and the evolution tracks that give a dwarf its temperature are not part of this change; the census of Kirkpatrick Table 15 is not machine readable here (NOT REACHED), so no count-per-bin target exists.
- The recipes were not read by the generator when this was built (GEN.517); the type draw is wired in section 3.

## 3. Wired into the generator (GEN.760)

`brown_dwarf.spectral_type(temperature_k)` draws the type of every brown dwarf through this map: the temperature it reads on the map is the dwarf's own plus a Gaussian offset of the scatter of the piece its mean type is on (134 K among L dwarfs, 79 K among T and Y dwarfs), and the type is the map read backwards (`brown_dwarf_types.number_of_teff`, `drawn_number`) at that temperature. The L/T transition, where the map is nearly flat (18 K per subtype), spreads a type over several subtypes, as the census does; a dwarf that reads as L0 or hotter after the offset is a late-M dwarf, named by the stellar scale.

The draw is on a stream of the temperature alone, so a free brown dwarf's page, which reads the stored effective temperature, gives the type the generator drew, and a brown dwarf member of a multi-star system stores it with the member. The old inversion of the mean relation (`tuning.BROWN_DWARF_SPECTRAL_TYPE_FITS`, `brown_dwarf.type_temperature_k`, `spectral_number`) is gone: the data file is the one table. Source tier: Verified ([R2374: Table 13]); the offset is the paper's scatter (tweak scale L0); the field-dwarf map's missing gravity term for young dwarfs stays a gap (section 2). Generator epoch: the types of brown dwarfs change (GEN.760).

## 4. The hot end: the hydrogen-burning limit by age (GEN.494)

The two edges the map could not give are the hot ones. [R0254] (A&A 577, A42, section 3.3; READ) put the hydrogen-burning limit at 0.07 solar masses and its temperature at about 2900 K at 1 Myr and 2950 K at 10 Myr. The BHAC15 table held in the project files (`baraffe_isochrones`) gives the whole curve for the limit mass: 2897 K at 1 Myr, 2946 K at 10 Myr, 2794 K at 100 Myr, 2047 K at 1 Gyr, 1627 K at 5 Gyr and 1626 K at 10 Gyr. `physics/substellar_ceiling.py` reads it (`ceiling_k(age)`, `HOTTEST_K` 2950 K, `oldest_age_gyr(temperature)`).

What it settles: the late-M-like brown dwarf (row g1, 2238 K and up) is a young object. A brown dwarf is as hot as 2300 K only until about 0.55 Gyr (as hot as the L0 edge, 2238 K, until 0.64 Gyr), and at 1 Gyr the hottest is already near 2050 K, at 5 Gyr near 1630 K (the mid-L dwarfs). The 2800 K of the rogue-body clip is 150 K under the table's peak and is the same for every age. The row's hot edge becomes 2950 K at the wiring step, and the age dependence stays in the cooling track itself once the tracks are the Baraffe ones (GEN.948), which have the limit built in. Row g4 (T8 to Y0) is not a separate type; its edges are the map's own (T8 686 K, Y0 460 K, scatter 79 K), so it is a derivation from Table 13 and the "design choice" in its limit text goes when the wiring step edits the rows.

**Wired (GEN.494, step 2, train 7):** `brown_dwarf.effective_temperature_k`, `luminosity_sol` and `state` and `rogue_surface.giant_internal_flux_w_m2` cap a body below 0.07 solar masses at `substellar_ceiling.hottest_k(age)` and one at or above it at `HOTTEST_K`; `tuning.ROGUE_MAX_EFFECTIVE_TEMPERATURE_K` (2800 K) is gone, and the hot edge of row g1 is 2950 K. The cooling tracks (`baraffe_isochrones.cooling_track`, GEN.948) are BHAC15 already, so the change is mostly at the peak: a brown dwarf under about 0.1 Gy that the old clip held at 2800 K now keeps its track temperature up to 2950 K (0.065 solar masses at 10 My: 2924 K, was 2800 K). Under the GEN.948 tracks no body below the limit mass is above the limit mass's curve (the table itself is BHAC15), so the cap binds only at the peak. A cap asked for with the Baraffe table missing and told to disregard that gives the 2950 K peak at every age (recorded as a stand-in); the cooling track itself needs the table.

Not done: the limit moves with [Fe/H] (GEN.936); the young low-gravity L dwarfs of section 2 still need a gravity term. Source tier T1 (Verified), tweak L0 (registry SB1 and SB2).
