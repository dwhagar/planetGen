**Headline:** The planets round an inner pair of a larger system now come with the belts the chain leaves: a disk of small bodies round a close pair is a belt of its own, hosted at the pair rather than at one star.

### Changed

- GEN.332, slice C: the planets of an inner pair node (slice A, `StarSystem.node_hosts`) are followed by `_add_belts`, the same belt step a member's planets get (`belt_birth`, GEN.1068): the stretches the node's chain leaves unfinished inside its P-type window, never inside the pair's stability limit. A circumbinary disk is stored as an `asteroid_belts` row hosted at the node, as the Boss decided for the card on slice C, and not as a note on the pair: the table gains `host_node_id` (migration 0254, schema v254; `star_id` NULL, `parent_type` `node`, `parent_id` the node's number), `insert_asteroid_belt` and the edit save (`db/edits.py`) take it, and `load_star_system` reads a node's belts back with its planets, in one orbital-index space. The dust note on a star (`disk_kind`) is left alone for a node: the belt row is the disk.
- Reuse: the belt draw, its origin and history, the member list, the composition rows and the describe, map and scene code are the ones a star's belts already have; the node belt differs only in its host and its storage keys.

generation-output: changed
generation-affects: fill
