**Headline:** The merge train's gate, number checks and machine start-up are scripts in the repository, so a restart or a second machine no longer loses them.

### Changed

- `scripts/train_gate.py` lists the test files a train reaches, splits them over several machines by their last durations, runs them one at a time with a results file that a restart resumes from, and reports only the failures that `scripts/data/known_reds.json` (the failures present on a clean main) does not explain.
- `scripts/train.py` merges a train's heads in stacking order, shows the next migration, epoch, stand-in section and citation numbers, finds a number claimed by two heads or already taken on main, and runs the rebuild, golden and final checks.
- `scripts/lane_up.sh` starts Redis and MySQL after a restart and sets the test variables.

generation-output: unchanged
