### Added

- Every generation job is now a tree, with timing on every node (ADM.12):
  a Generate page job, its steps, each `generate.py` run, its phases
  (skeleton, bright stars, population), its work queues and their sectors
  or layers, each with its own state, start, end and duration, and
  totals, timings and an ETA added up from the children. One-worker runs
  are recorded too. Control schema v7: run `update.sh` (or `update.ps1`)
  after updating.
