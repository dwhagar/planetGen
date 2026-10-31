**Headline:** A System Map test that compared a planet's glow under two stars now compares the right thing, so main is green again.

### Changed
- **Glow test (TEST.528).** `test_a_planet_and_its_moon_glow_in_the_colour_of_their_air_under_their_own_star` compared the blue channel of the limb glow under a Sun-like star and under a red dwarf. The glow colour is scaled so its brightest channel is 255, so the blue channel was 255 under both stars and the test could not tell them apart. It now compares how much bluer than red each glow is (the Sun's air is clearly blue, a red dwarf's is whiter). No code changes.

No generator output changes.

generation-output: unchanged

TODO items resolved: TEST.528.
