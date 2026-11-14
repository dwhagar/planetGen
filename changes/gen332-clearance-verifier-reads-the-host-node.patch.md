**Headline:** The overlap check now holds a planet or belt round an inner pair of a nested star tree to that pair's limits, so it no longer calls a correctly placed body an overlap.

### Changed
- **The independent clearance verifier reads `host_node_id` (GEN.332).** `physics/clearance.py` took a body with no star as going round the one pair of stars of the tree, which was right until GEN.332 put planets and belts about the inner pair node and the root pair node of a nested tree. `stored_rows` now reads `planets.host_node_id` and `asteroid_belts.host_node_id` and `_host_of` puts a body about the node it names, so it is held to that node's Holman-Wiegert limits and S-type limits. All 16 multi-star sweep systems of `test_a_generated_system_holds_no_overlap` that Train 27 reported (quadruple 4, 7, 8, 11, 16, triple 7, 9, 10, 11, 12, 16, 17, 19, 22, five 3, 5) are clear with it; the generator was right and the test's host was wrong.

generation-output: unchanged
