**Headline:** Two moon tests that had stopped testing anything (one sampled systems with almost no moons, one compared a giant's moons with the radius the giant has shrunk to) now test what they say.

### Changed

- `src/tests/test_moon_history.py` (TEST.657): `_systems` made plain systems with `SystemConfig()`, whose `MOONS` is off by default. Since GEN.962 and GEN.1082 to GEN.1084 a rocky planet has a moon only by a giant impact (1 in 12) or capture, so the 16 sampled systems held 2 moons and `test_every_planet_with_moons_ends_with_moons_it_can_keep` failed on `seen_moons > 5`. `_systems` now sets `config.MOONS = True`.
- `src/tests/test_moon_disk_zone.py` (TEST.657): `test_the_disk_moons_of_generated_giants_sit_inside_the_zone` compared a disk moon's birth distance with 30 radii of the giant's final radius. A giant's radius is derived again at the star's final age (`_settle_to_final_age`) after its moons are drawn, so the final radius is a few percent smaller than the one the zone was made for. The test now records the zone each giant's moons were drawn in (`moon_disk_zone.zone_km`), asserts the birth distance lies inside that zone, and that the zone is 30 radii of the giant to within 6 percent.

generation-output: unchanged
