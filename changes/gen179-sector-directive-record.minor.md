**Headline:** A sector now keeps the directive it was drawn under and how the draw came out, and its page and `check-db` show it.

### Changed

- A run given `--directive` minimums (GEN.96) used to throw the attempt record away once the winning draw was saved. The sector now keeps it: two new tables, `sector_directives` (the directive's digest, the status `met_naturally`, `met_after_k_attempts`, `met_forced` or `unmet`, and the draws made) and `sector_directive_minimums` (for each minimum asked: its key, what was wanted and what the kept sector got), written in the sector's own save transaction (migration 0275, a placeholder the merging train renumbers). A sector drawn with no directive stores nothing, and nothing is backfilled for sectors made before.
- `GET /api/sectors/<id>` gains `directive` (`digest`, `status`, `attempts`, `minimums` of `{key, wanted, got}`; `null` when none). The sector page shows a "Drawn under" row in its At a glance panel, for example "at least 2 habitable worlds (got 3): met after 4 draws".
- `check-db counts` gains a `directives` piece: a record whose status says `unmet` while every minimum was met, or the other way round, or that has no minimum rows, is reported as damage.
- The fingerprint leaves the two tables out: they record what was asked of a run, not what it made. This is also the record DB.17's parity repair reads when it draws a sector again.

generation-output: unchanged
