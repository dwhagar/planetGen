**Headline:** A brown dwarf now has the spectral type, temperature and radius its mass and age give (late M, L, T or Y), where a free-floating one had only a random size and a multi-star member had no state.

### Changed
- New `physics/brown_dwarf.py`: Burrows et al. 2001's two cooling laws (luminosity `t^-1.3 M^2.64`, temperature `t^-0.32 M^0.83`) fix the temperature and the radius, and Kirkpatrick et al. 2021's Table 13 (525 L, T and Y dwarfs within 20 pc) turns the temperature into a type: L0 at 2,238 K, T0 at 1,255 K, Y0 at 461 K, late M above. `docs/design/brown-dwarf-state.md` names every source and what is left out.
- A free-floating brown dwarf's mass is drawn from the census's `dN/dM ~ M^-0.6` between 13 and 80 Jupiter masses (it was log-uniform) and its radius follows the laws, 0.07 to 0.14 solar radii (it was Jupiter's times a random 0.75 to 1.1). Its page names its type ("a free-floating T5 brown dwarf", "Rogue Brown Dwarf, Type T5"), read from its stored temperature.
- A brown dwarf member of a multi-star tree has its type as the member's `state` (it was empty) and the radius the laws give at the system's age (it was 0.1 solar radii for all); the scatter row's members list reads "Brown dwarf (T5)". The brown dwarf class page gives the temperatures where L, T and Y begin.
- Temperatures of the hot L dwarfs are 2 to 4 times too few against the census (a strict expected failure in `tests/test_brown_dwarf.py`).

generation-output: changed

TODO items resolved: GEN.252 / Added: hot L dwarfs (above 1,500 K) come out 2 to 4 times too rare because the old-age power law is too cool for the heaviest brown dwarfs (a Baraffe or ATMO 2020 track would close it)
