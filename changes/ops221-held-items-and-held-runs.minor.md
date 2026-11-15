**Headline:** A unit that fails for good no longer ends the run: it is held, the rest goes on, and a run held for a person can be carried on after the fix.

### Added

- OPS.221 (job engine stage 3, L7): a unit that fails for good (the error reached `units.unit_failed`, or a kill with no retry left) is written to its stage's held hash and its job row is `failed`; its queued retry finds it held and does nothing. Every other unit goes on. When every unit of the stage is done or held, the run takes the new state `held` (`runrecord.hold_run`): it frees the machine, shows as Held on the Queue page and in `GET /api/admin/work`, no worker takes a unit of it, and the ledger's run row says `held` with why. Retry held (`held.retry_held`; the Queue page, `planetgen queue retry-held RUN_ID`, `POST /api/admin/work/<id>/control` `{"action": "retry-held"}`) takes the machine back, sets each held job to `waiting` and queues it again with the arguments its row holds (`units.requeue_unit`, reused from `enqueue_units`). Cancel ends a held run. `planetgen queue held` lists them.
- A unit that has no ledger row, and a piece of an opening, a conclusion or an end of a stage, still fail the run as before until OPS.218's jobs carry the same state. The failure kind seam for OPS.220 is the one call `units.unit_failed` makes: everything that reaches it is lasting.

### Fixed

- A held run's starter ends with the reason and does not end the run (`runrecord.finish_run` leaves a held run alone unless it is the cancel of it).

generation-output: unchanged
generation-affects: none
