**Headline:** A unit can run as a job process the template forked, built and tested, not yet used.

### Added
- `queue/forkedjob.py` runs a unit in a forked job: its own Redis connection, the same `units.run_unit_on` the RQ job runs, and a say on the job's pipe when it ends. A stop request or a supervisor that went ends it between two pieces; a failure fails the run as RQ's failure callback does (`units.fail_unit`).
- The engine asks an optional `halt` before every piece (`engine.run(..., halt=)`), a local read.

### Changed
- `units.run_unit` and `units.unit_failed` are split so the RQ job and the forked job share the code that runs a unit and that fails a run (`run_unit_on`, `fail_unit`); nothing changes for the RQ path (OPS.226, the second slice).
