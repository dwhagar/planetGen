**Headline:** Kill matrix over a whole galaxy run and the warm-up after it, with every engine transition named (OPS.246).

### Added

- A seeded chain test: the sectors, the sky charting and the settling of the paths, then the object index warm-up, with units killed at drawn commits along the whole chain, ends as an undisturbed chain ends (galaxy, stars, cells, spheres, and an index holding every object once). The gate runs 2 seeds with at most 4 kills each; `PLANETGEN_ENGINE_CHAIN_SEEDS` and `PLANETGEN_ENGINE_CHAIN_KILLS` raise both for the end-of-version hunt.
- A table of the engine's named transitions (claim, commit, end of a unit, hold, retry, cancel and finish, end of a stage, window write, launcher hand-over, template and kernel kills) with the test that kills each, and a test that every one of those tests exists.
