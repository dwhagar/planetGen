**Headline:** The test that saves a rogue planet and a moon no longer depends on seed 3's system having a moon.

### Fixed
- `test_a_rogue_planet_and_a_moon_of_it_are_saved_and_read_back_with_no_giant_to_share_an_orbit_with` drew a G2V system at seed 3 and took its first moon; the ring coefficient of GEN.1202 moved that system's draws and it has no moon any more. The test now tries seeds 3 up to 14 and takes the first system that has a moon. No generator change.
