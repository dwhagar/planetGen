**Headline:** A random start is picked by the start stage's cut from the run's seed, not by the starter's own draws.

### Changed
- **The random start of `planetgen sectors` and `galaxy` is picked when the start stage is cut (PERF.123, third slice).** The run draws a seed of its own; the starter picks each start from it for the estimate, and the start stage picks the same one again from the same seed and its position in the list (the starts before it keep it from their neighborhoods), so the stage's record carries the seed and no sector list. The run seed is not an option, so nothing a given set of options produces changes; `--debug` says the seed.

TODO items resolved: none (PERF.123 slice 3 of 4). Added: none.
