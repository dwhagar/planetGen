**Headline:** The warm-up in a log leaves one progress line a minute at most, not one per chunk.

### Changed

- `planetgen warm` without a terminal (a log, `tee`, the background warm-up, `--follow` into a file) prints a stage's progress as `galaxy-nodes 16% (28 of 180 chunks), 1:56 elapsed, about 10:10 left; warm-up about 11:07 left (rough guess)`, and only when a minute and ten percentage points have both passed since the last one (OPS.99), so a stage of hours leaves a line at each tenth and a stage of minutes one a minute. The old line per chunk, with its nine facts, and the ASCII bar line after each stage are gone.
- A stage's closing line carries what the chunk lines used to: how many chunks ran at once and how long one took on its worker, and `step 2 of 5` with the whole warm-up's time left. A stage of one unit prints its own line only.
- With a terminal the `rich` bars are unchanged (update.sh shows them now that nothing pipes the warm-up through `tee`, OPS.96).

generation-output: unchanged
