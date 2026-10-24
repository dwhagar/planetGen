**Headline:** Two tests of scattered stars no longer assume the row's star is the system's primary or its first star.

### Fixed
- **`test_plan_scatter_stores_stars_and_fill_builds_their_systems` failed on clean main (TEST.338).** It joined the stored row to the built star with the role `primary` or `single`, but a system built around a scattered star is a tree (GEN.631): the row's star is member 0 of that tree and the `primary` and `secondary` roles are the top pair of its shape, so in 2 of 7 systems the row's star was a `member`. The test now finds the row's star in its system by the luminosity the row stored and checks its [Fe/H] there.
- **`test_a_pre_placed_hot_star_rolls_the_same_region_as_the_star_its_sector_builds` failed on clean main (TEST.340).** It built the row's star as a plain `StarSystem`, which draws a tree of its own, where the fill builds a row with no stored members as a single star (`stored_trees.SINGLE`). The same cause as TEST.311.
TODO items resolved: TEST.338, TEST.340. Added: none.
