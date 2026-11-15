**Headline:** finish-run now decides which stages of a stopped run are finished from the run's ledger rows, so it no longer depends on Redis keeping the run's stage records.

### Changed

- OPS.238, start paths: `finish_run.stage_finished` reads the stage's end rows (`planrows.ended`): a stage is finished when its `release` row is done (its units, conclusion and statistics are done and the next opening was queued). A stage that opens with no unit to wait for (a stage a run carries over from the run before, or one cut to no unit) now marks its three end rows done where `wrapper.advance` passes it, so the rows say it ended as they do for any other stage.
