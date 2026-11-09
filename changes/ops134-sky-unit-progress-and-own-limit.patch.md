**Headline:** A sky unit's progress counter now moves with its work, and its time limit comes from its own predicted time, not the average of the sector units it shares a stage with (OPS.134, heartbeat part 1; the failure of the 11-08 full-galaxy run).

### Changed
- **The counter moves.** `sky_charting.chart` takes `progress(done, total)`, called after every batch and every tier, in the predicted cell visits (`sky_tiers` `visits`) of the tiers the ledger lacks (`owed_tiers`, now shared with `_sky_units`). Within a tier the counter follows the cells walked in a straight line to 90% of the prediction, then creeps toward the tier's end without reaching it (`tier_progress`), so a tier that runs over its prediction still moves the counter. `_sky_unit` passes `report` as `progress`; it used to call `report(steps_done[0], steps_done[0] + 1)` with a `steps_done` that never changed, so the unit read 0 of 1 for its whole run.
- **The unit has its own time limit.** `_sky_units` puts the unit's predicted seconds (visits x `BLOCK_VISIT_SECONDS` x `SKY_WRITE_FACTOR`) in the unit's arguments as `expected_seconds`, and `units.enqueue_units` gives a unit the limit of its own prediction (ten times, 15 minute floor, 6 hour ceiling) when that is longer than the stage's. The stage's expected seconds are what earlier runs of the stage took, mostly its small sector units: 188.5 s gave the 11-08 run a 1,885 s limit for a cell that takes 1,826 s, and RQ ended it. A measured cell now gets about 18,260 s. Any other command may name `expected_seconds` the same way.
- No cell is drawn differently and no star changes.

generation-output: unchanged
generation-affects: none
