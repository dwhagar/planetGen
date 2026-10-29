**Headline:** The notes and code comments that describe the work queue now describe the queue as it was built.

### Changed

- **The queue and scatter notes say what was built (DOC.20).** The architecture overview (its Generate flow and its table of the queue package), the benchmark report, the performance and time-left note, the queue audit, the resumable-runs note, the library-migration plan, the scatter queue study and the Redis review each say that a run lives only in Redis (a run record with stages, units that are RQ jobs on the bulk queue, the last unit of a stage queuing the next) and that the MySQL job tree, the lease and the task and subtask model are gone. The docstrings of the queue package, of the command registry and of the wrapper no longer say the massive-star scatter is the only stage on it or that stages are still inline, and the old generation queue's name says that nothing is queued on it. Documentation and comments only; no code behaviour changes.
