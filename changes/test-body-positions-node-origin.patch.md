**Headline:** The scene-position test no longer fails when a drawn system has planets round the inner pair of a triple star.

### Fixed

- `test_the_scene_positions_at_the_epoch_match_the_stored_ones` read a planet's or comet's `around` straight out of the positions map, which holds no `node:N` entry (the center of mass of a tree node, which the scene's own code resolves inside `positions_at`), so a random draw of a system with planets round an inner pair raised `KeyError: 'node:1'`. The test now finds a node's center the way `positions_at` does. No assertion is loosened and nothing generated changes.

generation-output: unchanged

- `test_a_planet_made_from_a_belt_names_that_belt_on_its_sheet` read the reloaded system's planets through the removed `.planets`; it reads `primary.orbiters`.
