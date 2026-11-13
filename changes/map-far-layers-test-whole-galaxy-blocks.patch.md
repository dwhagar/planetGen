**Headline:** The map's far-layer browser test holds the whole-galaxy block part to its own bound again.

### Fixed
- **A browser test was red on main since ADM.148 to ADM.151.** `test_layers_far_from_the_center_are_built_relative_to_origins_of_their_own` used to open the map in the "made" view, which drew no whole-galaxy blocks. With that view gone it opens the ordinary view, whose one part of block prisms (about the galaxy's middle, 15.5 kpc across in the fixture) was held to the 10 kpc a window's layers keep. The part is documented as spanning more (`static/camerarelative.js`), so the test now holds it to 50 kpc and every other layer to 10 kpc. No page code changes.
generation-output: unchanged
generation-affects: none
