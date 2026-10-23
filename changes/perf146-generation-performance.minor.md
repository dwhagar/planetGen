**Headline:** `planetgen benchmark --suite` times every generation step on a galaxy of a chosen size, and the first round of fixes it found makes the massive-star scatter at the default size take minutes, not hours.

### Added
- **`planetgen benchmark --suite`** (PERF.146). Builds a scratch galaxy of a chosen size (`tiny`, `small`, `medium`, `large`), runs each generation step as its own measured command (wall time, machine CPU, peak memory, database questions and rows, time per star or cell, the stages the step reports) and prints a table ranked by share of the whole run; `--repeats` gives the run-to-run spread, `--steps` picks steps, `--step-timeout` stops a step that runs away and records it as not finished, `--output` writes the JSON. See `docs/cli.md`.
- **`docs/design/generation-performance.md`**: the measured cost of every step, a ranked list of what to speed up next, and the py-spy recipe the profiles were taken with.

### Changed
- **A queue worker builds the galaxy's component field and sampling tables once, before it forks a job, not once per unit.** The scatter units were each rebuilding the field from about 22,000 clusters and satellites; the massive-star scatter at the default size took over two hours and now takes about three minutes. The drawn values are unchanged.
- **The mean system mass at an age is about four times faster** (the remnant states and turn-off mass are worked out once per call, and the lifetime terms once per metallicity), with the same numbers.
- **The scatter's band tiers are looked up together** (one pass over the density arrays for all the places a component's stars land in) instead of one band at a time, bit for bit the same tiers; the supervisor's tier lookup for marked cells is vectorised the same way.
