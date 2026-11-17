**Headline:** A generation run whose units are forked jobs no longer starts bulk workers that sit idle, and its `--workers` count still limits how many units of a stage run at once.

### Changed
- A forked run's stages ask for no bulk workers (`wrapper.stage_ensure`), and the opening of a stage no longer asks the other bulk workers to exit (`_clear_the_floor_for_the_opening`): there are none (OPS.228 slice D2b-1).
- The stage's `open` row records the run's `--workers` count (`cap`), and the launcher hands out no more units of a stage at once than that; without `--workers` the memory gate alone decides, as before.
