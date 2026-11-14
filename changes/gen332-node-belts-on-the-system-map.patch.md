**Headline:** The system page's map now draws a belt that circles an inner pair of stars in that pair's own scene, beside the planets of the pair.

### Changed

- GEN.332, slice C display: `web/maps/systemmap.py` (`_node_pair_hosts`) groups the belts that carry `host_node_id` (slice C stores a belt of an inner pair node's circumbinary disk with no star) with the planets of the same node, so the pair's marker and its `node-<n>` scene hold them. A node that has a belt and no planet gets its marker and scene as well. Before, such a belt was drawn in the whole-system scene round the system's own pair, at a distance that belongs to the inner pair. A system stored before the tree keeps its belts where they were.

generation-output: unchanged
generation-affects: none
