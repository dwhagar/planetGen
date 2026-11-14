**Headline:** The star and phenomenon scatters, the nebula field and the map luminosity fill now report what they have done a piece at a time, so a long run shows real counts and can be stopped or paused between pieces.

### Changed

- OPS.198: `scatter-mass`, `scatter-luminosity`, `scatter-band`, `scatter-phenomena`, `draw-nebula-field` and `fill-map-luminosity` hand the engine `pieces`, `do` and `gather` (OPS.191) in place of a `unit` function with a loop of its own. The pieces are a shell (with its part beyond the outline for the phenomena), the far shells together, the clusters' and satellites' stars, a first-index slab of nebula cells, and a range of 20,000 ids of the map luminosity fill. The engine counts pieces done of pieces in all, honours a cancel and a pause between two pieces, and each piece commits through `ctx.commit()`. The second level counts the stars or phenomena placed and the columns of nebula cells (`ctx.step`, `ctx.advance`, OPS.192). The unit's result has the same fields as before, so the stage's followers, statistics and speed samples are unchanged.
- Removed with the old loops: `scatter_shells_unit`, `scatter_phenomena_unit`, `draw_nebula_field_unit`, `fill_map_luminosity_unit`, the shells and phenomena tasks, `_nebula_field_task`, `_ReportChannel` and `CHANNEL_INTERVAL_SECONDS`. The free-text phases of the nebula field and the map luminosity fill are gone with them.
- `engine.Context` gains `kept(name, make)`, so a unit reads the galaxy's shape and its components once for all its pieces.
- A map luminosity fill piece covers `MAP_LUMINOSITY_PIECE_IDS` (20,000) ids with its own statement and commit; the block light parts of its cells are put right in the same piece, and the values are those of one statement over the whole range.

generation-output: unchanged
generation-affects: none
