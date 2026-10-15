### Added
- **A big-galaxy query budget tier: `scripts/test.py budget`** (TEST.133). It builds, once and keeps, a database of 2,000,000 systems and 150,000 sectors by cloning one real saved system, runs the Systems and Sectors page and count queries on it, and fails a page that reads far more rows than it shows (counted from the server's own read counters, so the result does not depend on the machine). It is opt-in and not part of the merge gate or CI. Known slow pages are expected failures that name the item that fixes them. A 3,000-system version of every case runs in the ordinary suite so the tier does not rot.

### Changed
- **`scripts/test.py` ends every run with one clear line: `PASS`, `FAIL` (with the test counts) or `INTERRUPTED`.** A run that stops without any of the three was killed. `docs/testing.md` says to start long runs in the background with a timeout of at least 90 minutes.
