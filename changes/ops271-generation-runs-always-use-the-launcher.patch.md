**Headline:** The `queue.forked_jobs` setting is removed: a generation run that has a job ledger always hands its units to the supervisor's launcher.

### Changed
- `queue.forked_jobs` (and `PLANETGEN_FORKED_JOBS`) is gone, since the launcher path is the only one that is kept for a run with a ledger (OPS.271). A generation run with no ledger (the control database does not answer or is not migrated yet) and the other run kinds go to workers as before.
