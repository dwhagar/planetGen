**Headline:** The test that checks every generated planet is in a state the tide allows now samples systems with planets switched on, so it meets planets that are still turning freely.

### Changed
- `test_every_generated_planet_is_in_a_state_the_tide_allows` took the first 12 default systems; since the star changes they are all dwarf stars with planets inside the tidal-locking distance (25 planets, every one locked), so its "some planets are still free" assertion failed on main. It now uses the 12 seeds of the sheet audit with planets switched on, as that audit does. Nothing a generator draws or stores changes.
generation-output: unchanged
