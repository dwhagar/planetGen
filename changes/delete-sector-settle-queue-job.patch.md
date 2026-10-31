**Headline:** The queued sector-settle job that no page or route reached any more is gone.

### Changed
- **`api_jobs.settle_sectors` and the test helper `queue_settle` are deleted (TEST.531).** Since PERF.149 no route generates a sector on a visit, so only a test helper queued this job. The tests that existed only for the queue are removed (a sector generated on the spot queues a settle job, a settle that can be neither queued nor run, the conftest rule that ran settle jobs inline); the behavior a run still has stays tested through `run_galaxy.settle_named_sectors`, which `make_sector` now calls directly. No generated output changes.

generation-output: unchanged

TODO items resolved: TEST.531.
