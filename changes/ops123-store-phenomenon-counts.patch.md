**Headline:** Storing how many phenomena of each class the scatter placed is a step of its own in the run, no longer the hidden end of the shells' pass (OPS.123 slice 2, phenomena scatter).

### Changed
- **`store-phenomenon-class-counts` is its own command and stage.** The `conclude` of `scatter-phenomena` is gone; its effect is a one-unit stage, "Store the phenomenon counts" (`phenomena-counts`), between the shells and the indexes. It takes what each shell unit found through a new deferred parameter, `wrapper.stage_results(stage, field)` (every unit's result of a stage, in unit order, from this run or from the run that did the stage when this one carries it on). A galaxy run counts one more stage (21), a phenomena scatter one more. Reused: `store.record_phenomenon_scatter_classes`, the one-unit whole expansion, the deferred parameter machinery of `opened_by` and `scatter-settings`, and the stage pattern of PR 6349.

generation-output: unchanged
generation-affects: none
