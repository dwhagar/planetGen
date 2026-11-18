**Headline:** the job engine is now tested on whole random runs with crashes anywhere.

### Added
- OPS.245: `tests/test_engine_random_runs.py` draws a run from a seed (two to four stages of units made of pieces that commit, stages written in windows by the spawner, units that fail for good so the run is held and retried, kills at commits, the spawner stopped after its rows went in or after its cursor moved) and checks that the run ends done, every piece of every unit is written exactly once, no unit began before the stage before it was released, every windowed item has one row and one finished attempt, and every end row is done. `PLANETGEN_ENGINE_SOAK_SEEDS=200` runs more seeds.

generation-output: unchanged (tests only; no golden, no epoch, no lock item)
