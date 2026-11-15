**Headline:** a run's whole sequence of stages is now stored as waiting rows in the ledger the moment the run starts, so the job view can say which stage each one waits for.

### Added

- OPS.238, third part (stored plan): `planrows.record_plan` writes one waiting `open` row per planned stage when `wrapper.run` starts. Each row's `needs` names the stage before it (`after:<stage>`) plus the artifacts that stage requires (`needs.required_stage_keys`), and its `makes` lists the artifacts it produces (`needs.makes_of`). When a stage opens, `planrows.record_cut` fills that row in with its plan and marks it done instead of adding a new row. The job view words the wait as "stage <key> to end".
