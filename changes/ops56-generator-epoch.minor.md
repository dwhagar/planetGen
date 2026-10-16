**Headline:** The generator now has an epoch counter, and a pull request that changes what a seed generates must bump it.

### Added
- `src/tests/golden/epochs.json` holds the generator epoch, `planetgen.util.generator_epoch.generator_epoch()` reads it, and `python -m tests.golden_update --reason "..."` bumps it (refusing under CI, when no golden file changed, and a second time in one change).
- `scripts/bump_version.py --check-pr` fails a PR unless the golden files changed, the epoch went up by one, and the note carries a `generation-output: changed` line, all three or none. Stamping the release writes it into the epoch's `first_release`.
