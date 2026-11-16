**Headline:** The warm-up stages for sectors, systems, phenomena and the Galaxy Map's nodes are now cut into units by the queue engine from the profile and the free resources, as the three galaxy stages already are.

### Changed
- `warm-sectors`, `warm-systems`, `warm-phenomena`, `warm-galaxy-nodes` and `warm-galaxy-deep-nodes` are flagged byte identical (PERF.260 already proves each stores the same cache bytes at two cuts, `tests/cut_neutral_recipes.py`), and their `expand` calls `splitter.cut` (OPS.243) instead of a typed `units.coarsen` size. `OBJECT_CHUNK` and `NODE_CHUNK` stay the suggested size the stage keeps when the profile has nothing usable, in at most `units.MAX_UNITS` units, so a fresh install cuts as it did. The stage plan records the engine's decision.

generation-affects: none
