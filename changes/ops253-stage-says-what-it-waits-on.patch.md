**Headline:** A stage that waits says what it waits on, in the run line and the log.

### Fixed

- OPS.253: a stage that finishes no unit says what it waits on. The run line of a stalled stage names what its units are doing ("waiting on progress (4 waiting on sky lock, 1 working)"), the job log says so every 30 minutes while no unit finishes, and the bright-star backfill's bar counts a sky cell by its own counter instead of standing still for the hours one cell takes.
