**Headline:** A sky's mass rings are costed by what their walk does, so a sparse ring is no longer counted one visit to every cell of its sphere (GEN.1485).

### Changed
- **`sky_charting.cost_rings`.** The queue's sky units (`run_galaxy._sky_units`) now give each mass ring the visits `bright_stars.mass_sphere_walk_visits` counts over its drawn sphere: for each row, one visit to every cell of its window when a band of the ring is over `CANDIDATE_PATH_BOUND`, else one to every block of the window. It is the walk's own rule (`_sphere_batches`) applied to the rows with no star drawn, looked at on about 2,000 rows of the sphere and scaled. `sky_tiers.build` still costs a ring by its cells, since it has no galaxy to ask. The estimate of the time a sky unit takes drops for sparse rings (the 5 to 8 and 8 to 14 solar mass rings, which the audit found overstated); no cell is drawn differently and no star changes. Reused: `sphere_rows`, `bound_terms`, `bound_of` and `_row_chunks` of the walk itself, and the ring bands of `mass_sphere_batches` (now `_mass_paths`, shared by the draw and the cost) so the two cannot differ. Tests: the cost against an instrumented `mass_sphere_batches` for a sparse and a dense ring.

generation-output: unchanged
generation-affects: none
