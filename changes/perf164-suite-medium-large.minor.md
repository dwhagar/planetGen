**Headline:** The generation performance report now covers the medium and large galaxy, and the benchmark suite can carry on after a stopped run.

### Added
- **`planetgen benchmark --suite --resume SUFFIX`** carries on in the scratch databases of a run that was stopped or failed, with `--steps` naming the steps left. The suite now keeps those databases when a step did not finish, prints the suffix at the start and prints a line as each step ends, so a run killed part-way (a container restart) keeps what finished.

### Changed
- `docs/design/generation-performance.md` has the small, medium and large times side by side: the medium galaxy takes 46 minutes without the warm-up and 74 with it; the warm-up is the largest single cost and still fails on the search page answering 202.
