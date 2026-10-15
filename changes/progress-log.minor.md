### Added
- **A progress log and a build id for every test run** (TEST.150, the first step of resuming a run). `scripts/test.py` prints the id of the build it is running and keeps `.test-progress/<id>.jsonl`, a line for every test that finished, written as it finishes. The id covers the code, uncommitted files, installed packages, Python, MySQL and Redis versions, the command and the day. Nothing reads the log yet; `--resume` is next.

### Fixed
- **Two map browser tests wait for what they check instead of a fixed pause** (TEST.126 group). The Color by test now waits for the legend to be visible after the reload, and the IndexedDB tiles test waits (up to 20 s) for the tiles to be stored; both failed once under load and passed alone.
