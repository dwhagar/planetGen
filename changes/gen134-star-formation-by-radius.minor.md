**Headline:** The young stars of the disk now follow the galaxy's own star formation by radius: fewer in the outer disk and in the bar's gap, a ring at the bar's end.

### Added
- `galaxy/star_formation.py`: the share of a place's thin-disk stars that formed recently, derived from the galaxy's shape. The rate per area follows the gas to the power 1.4 (Kennicutt 1998) over a gas disk as long as the stellar one, the bar depletes the gas inside its length (a gap of 2.2 bulge scale radii), and the factor is normalised so the disk makes as many stars as before. Nothing is taken from the Milky Way: its profile is the test target in `tests/test_star_formation_profile.py`.
- The constants `STAR_FORMATION_LAW_EXPONENT`, `GAS_TO_STELLAR_SCALE_LENGTH_RATIO`, `BAR_GAP_RADIUS_RATIO`, `BAR_GAP_STAR_FORMATION_FLOOR` and `BAR_GAP_EDGE_POWER` in `tuning.py` (the last three are open defaults, named as such there and in the design note).

### Changed
- The young and intermediate populations carry the factor in `population_densities` and its array form, so sector profiles, the bright-star scatter and the nebula scatter follow it. No cell's total density changes. The sector profile of an inner cell shows fewer young stars (on an arm crest at 0.25 kpc: young 2% to under 1%, intermediate 15% to 4%) and the bar's core holds 9% to 11% of its stars under 5 Gyr.
- The object-first band bound widens the young and intermediate weights by the factor's true range over each band's radii (`recent_factor_range`).
- The magnetar rate's young share of all stars in a Milky Way-sized galaxy is 0.44% (was 0.48%): the bulge now holds fewer young stars, and `test_gen_magnetar_scatter` caught the old figure. A full scatter holds about 9% more magnetars per young star.
- `scripts/star_type_by_radius.py` runs again (it called the retired `--arm-amplitude`).
- `docs/design/star-types-by-galactic-radius.md` records the derivation, the sources and the result. Proposal 2 (local B stars at 0.04%) is not met and the note says why.

TODO items resolved: GEN.134 (proposals 1 and 3; 2 is left open and described in the note) / Added: the star-formation history of the disk from the gas law (the follow-up to proposal 2)
