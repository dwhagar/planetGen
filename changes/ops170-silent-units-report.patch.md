**Headline:** The eleven jobs of a generation run that used to show nothing between their start and their end (the galaxy skeleton, the clears and empties, the phenomena and nebula draws, the luminosity fill, the habitability and ledger passes) now say each phase they are in and move a bar.

### Changed
- **Phase lines and progress for the 11 silent units (OPS.170).** `plan_skeleton_unit` (nine phases), `clear_bright_stars_unit`, `empty_nebula_field_unit`, `empty_phenomenon_scatter_unit`, `drop_phenomenon_indexes_unit` and `clear_scatter_unit` say each phase with the time since the unit began (`Emptying the bright stars (0:03 elapsed).`), set the unit's phase for the heartbeat and report the steps done. `draw_nebula_field_unit` reports a slab at a time, `scatter_phenomena_unit` a run of shells at a time (with a ticker line, OPS.169), `fill_map_luminosity_unit` its three phases, and `rescore_habitability_unit` and `refill_effect_ledger_unit` every 100 systems with a ticker line every 30 seconds. The units that come many to a stage keep their own begin lines at debug, so the stage log does not grow by a line for each.
- **A bar for the stages of one unit.** `UnitBarFollower` (the skeleton, the clears, the empties and the index drop) draws a share-done bar moved by the unit's reports; it appears only if the unit runs past `PROGRESS_BAR_SECONDS`, and is full when the unit is back.
- No unit takes `report` and ignores it any more; `test_unit_reports.py` fails when one does.

generation-output: unchanged
generation-affects: none
