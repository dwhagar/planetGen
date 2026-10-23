**Headline:** The test that clicks a star in a far sector waits for the star to hold still and for the panel to name it.

### Fixed
- **`test_a_star_in_a_sector_far_from_the_core_is_selectable` failed in three gates in a row (TEST.291).** The test worked out where the star was drawn once the map's centre had held still, clicked, and read the panel 200 ms later. On a loaded machine the star was then still moving: in a run that failed it was at (417, 299) when worked out and at (535, 310) a moment later, so the click landed beside it. The test now waits for the star's own place to hold for four looks, clicks it up to four times (working the place out again each time) and waits up to 3 s for the panel to name it. It failed four runs in four beside four busy CPU loops before and passed twelve in twelve after.

TODO items resolved: TEST.291. Added: none.
