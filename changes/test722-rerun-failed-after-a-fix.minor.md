**Headline:** `scripts/test.py failed` reruns only the tests that failed and the tests a fix can affect, within a time limit, instead of the whole selection again.

### Added
- **`scripts/test.py failed [--base REF] [--limit MINUTES]`** (TEST.722): after a red gate and a fix, it runs the tests that failed in the last run plus the tests `changed` selects for the changes since the base, and keeps the added ones within `--limit` minutes of recorded test time (default 20) on the workers asked for; what does not fit is counted in the output and left to the end-of-version run. It uses two new pytest options in `src/tests/depmap.py`, `--also-run FILE` (node ids that stay whatever the map says) and `--time-limit SECONDS`. Documented in `docs/reference/testing.md`.

generation-output: unchanged
generation-affects: none
