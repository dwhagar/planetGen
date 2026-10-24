**Headline:** Each update now says whether it changed the generated output.

### Changed
- **The update's version-key step reports it (OPS.15).** After it records the galaxy's row, `update.sh` prints `Generated output unchanged.`, or `Generated output CHANGED since the update of <time> UTC; the first case that differs is <case>.`, by comparing the generator battery's per-case digests (OPS.57) with the newest history row before it is replaced. The first update of a galaxy, and a galaxy whose last row recorded no battery, say they have nothing to compare. A new lock file or release alone does not count as a change of output. The live galaxy is not touched. `version_history.record` now returns the key and that sentence (`Recorded`).
