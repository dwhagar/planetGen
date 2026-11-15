**Headline:** The bright star backfill, the map's block totals refresh at the end of a run, the sector stages and the single system run now report what they have done a piece at a time, so a long run shows real counts and can be stopped or paused between pieces.

### Changed

- OPS.202 (second part): `scatter-backfill`, `refresh-run-blocks`, `generate-start-sector`, `generate-neighborhood`, `generate-sectors` and `generate-system` hand the engine `pieces`, `do` and `gather` (OPS.191) in place of a `unit` function with a loop of its own. The pieces are a chunk of `BACKFILL_CHUNK_SECTORS` sectors (the transaction the backfill already made of them) or one observer cell's sky tiers; the nebulae and then the sectors of the block totals; one sector of a sector stage (a stage run again leaves out the ones already saved); and the tries and then the saving of a system. The engine counts pieces done of pieces in all, honours a cancel and a pause between two pieces, and the second level counts the sectors drawn, the nebulae and blocks brought in, the cell visits of the sky tiers, and the tries of a system (`ctx.step`, `ctx.advance`, OPS.192). The unit's result has the same fields as before, so the stage's followers, statistics and speed samples are unchanged.
- The unit's galaxy shape, scatter limits, blocks' candidates and sky ledger are read once for all its pieces (`ctx.kept`), and `run_plan._draw_sector_masses` takes the last two, so the backfill draws what it drew before the engine cut it into pieces.
- Removed with the old loops: `backfill_bright_stars_unit`, `refresh_run_blocks_unit`, `generate_sectors_unit`, `system_unit` and `_sky_unit` (now the sky piece `_sky_cell`, which goes on after a lost MySQL connection as before).

generation-output: unchanged
generation-affects: none
