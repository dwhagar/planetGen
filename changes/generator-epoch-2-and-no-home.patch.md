**Headline:** The plan records that the galaxy keeps no home position, and retires the deterministic-math item.

### Changed
- GEN.135 (deterministic distances, radii and rotations through `planetgen.util.detmath`, generator epoch 2) and TEST.205 (the reserved-word scanner reads `ALGORITHM=` and `LOCK=` inside `ALTER TABLE` as options) are done in #1927 and leave the TODO list; TEST.77 notes that the epoch and golden-files rule is waived until it is built.
- Boss decided the galaxy stores no home or Sun position (GEN.368 row 11 goes; the Sun-radius uses become per-position calculations; the map starts where sectors are generated); GEN.388 is dropped, and the 9.4 plan, the roadmap and the galaxy types note lose their home position step.

TODO items resolved: GEN.135, TEST.205 / Added: none (GEN.388 dropped)
