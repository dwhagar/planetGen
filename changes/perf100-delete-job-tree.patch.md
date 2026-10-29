**Headline:** The old MySQL job tree is gone: its code and its four tables are deleted now that every run lives on Redis.

### Changed
- **The work queue's MySQL machinery is deleted (PERF.100).** `queue/work.py` keeps only `worker_count`, `cpu_count`, `lower_priority` and `task_seed`; the job tree (`open_node`, `job_node`, `load_tree`, the control requests, the lease, the stale-lease and pruning code) and the update guard's `WorkRunSource` are removed. Control schema v21 drops `work_jobs`, `work_job_args`, `work_tasks` and `work_lease` the first time `update.sh` runs it, and a new control database never creates them. `planetgen.cli.reset_jobs` and `planetgen.cli.update_guard` no longer take database options, and reset-jobs no longer clears any control database rows, because there are none. The Queue page, the update guard and the pause flag already read Redis. No generated value, epoch or golden changes.

TODO items resolved: PERF.100. Added: none.
