**Headline:** A generation stage whose longest unit ran over five minutes is now reported in the run's log.

### Changed
- **A unit-length ceiling (PERF.155).** New setting `queue.unit_ceiling_seconds` (default 300, "Unit length ceiling" under Speed and disk): the longest one unit of a generation stage is meant to run. When a stage ends and its longest unit ran longer, the run's log says so, naming the stage and how many units it was cut into.
- Every queued command now says whether its stage is one unit by nature or how long its longest unit runs at the default galaxy size, and a test fails for a new command that says neither, or says more than the ceiling.
