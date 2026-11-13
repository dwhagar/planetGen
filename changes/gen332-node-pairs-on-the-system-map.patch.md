**Headline:** The system page's map now shows the planets that circle a pair of stars deep inside a larger system: the pair is a marker in the whole-system scene, and clicking it opens a scene with only that pair's planets.

### Changed

- GEN.332, slice A follow-up: `web/maps/systemmap.py` (`_node_pair_hosts`) groups the planets that carry `host_node_id` by node and draws each pair as one body, as the generator's proxy does: a marker named for the pair (`Voranthis CD`, `names.bodies.node_pair_label`, which `StarSystem.node_label` now shares), with the stars' masses and lights added, the larger radius, the temperature weighted by light and the place of their center of mass, in the whole-system scene and a `node-<n>` scene of its own, reached through `_member_stars_scenes` like a further star. Before, those planets were drawn in the pair's own scene round the wrong center. The marker's id is `node-<n>`, which no page answers to, so the info panel offers no View or NAV links for it.
- The tree comes from the `orbits` and the stars' `parent_node_no` that `system_detail` already carries; a system stored before the tree (no star under the node) keeps its planets where they were.

generation-output: unchanged
generation-affects: none
