**Headline:** A job-engine unit that is stopped part way and run again starts after the pieces it finished, and gives the result an undisturbed unit gives.

### Added
- `ctx.commit()` stores the number of pieces the unit has finished and what each returned, in the transaction that holds the data and fences the attempt (`ledger.commit_item`; the table `engine_piece_results`, schema v262). The engine checkpoints after every piece.
- A command that says `resumes` (`CommandSpec`) starts a later attempt at `pieces(args, cursor)` and hands `gather` the stored results beside the new ones. Ten commands say so: the three star scatters, the phenomenon scatter, the backfill, the run's block refresh, the map luminosity fill, the nebula field draw, link-neighbors and settle-sector-paths. The others start from their first piece, as before, and a test lists each with why.
- `tests/test_checkpoints.py`: a unit killed at each of its pieces, run again, does the pieces left and gives the undisturbed result with no row written twice (a synthetic command, and the map luminosity fill on a real galaxy); every command with pieces is audited.

### Changed
- A new attempt of a job keeps the cursor of the one before it (`ledger.begin_attempt` no longer clears it).
- A piece's result goes through JSON before the unit keeps it, so a resumed unit and an undisturbed one hand `gather` the same kinds of value.
