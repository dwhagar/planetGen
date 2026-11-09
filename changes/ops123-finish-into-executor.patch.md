**Headline:** The four scatter commands no longer carry a finish hook: the work queue's executor records their speeds and counts from what each unit reports (OPS.123 slice 1).

### Changed
- **A unit reports what it drew and wrote; the executor records it.** `scatter_shells_unit` and `scatter_phenomena_unit` now return `samples` (kind, density, seconds, systems, stars) and `row` (objects written, shells touched), and the plan carries `row_totals` (the shell count only the cut knows). `finish_stage` calls `stages.record_unit_samples` and `stages.stage_row_metrics` for every stage, so the `finish=` hooks of `scatter-mass`, `scatter-luminosity`, `scatter-phenomena` and `scatter-backfill` (and `run_plan.finish_scatter`, `finish_scatter_phenomena`, `_speed_stats`) are gone. The five warm-family finish hooks stay until slice 2. The stored speeds and the stage rows are the same.

generation-output: unchanged
generation-affects: none
