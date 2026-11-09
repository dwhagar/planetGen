**Headline:** Emptying the earlier phenomena and taking off their indexes are steps of their own in a scatter run, no longer a hook of the phenomena shells (OPS.123 slice 3, phenomena verb).

### Changed
- **`empty-phenomenon-scatter` and `drop-phenomenon-indexes` are their own commands and stages.** The `prepare` of `scatter-phenomena` is gone; its two effects are two one-unit stages, "Clear the earlier phenomena" (`phenomena-clear`) and "Take off the phenomenon indexes" (`phenomena-unindex`), before "Scatter the phenomena". The first reads the same deferred scatter settings as the shells, so a `scatter --phenomena` that is refused because phenomena are stored is refused before anything is emptied. A galaxy run counts two more stages (18). `finish_run.REWIND_WITH` chains a stopped phenomena scatter back through the indexes to the clearing (and the special phenomena back through the scatter), so a carried-on scatter empties the table again first instead of adding to the rows already written. Reused: `store.clear_phenomenon_scatter`, `store.drop_phenomenon_scatter_indexes`, the deferred `scatter-settings` parameter, the one-unit whole expansion of `clear-scatter` and the stage pattern of the nebula field and bright star verbs. `clear-scatter` itself was not reused: it refreshes the block lights and resets the scatter epoch.

generation-output: unchanged
generation-affects: none
