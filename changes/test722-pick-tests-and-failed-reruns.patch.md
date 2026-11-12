**Headline:** `scripts/pick_tests.py` picks the tests a PR needs, stops at a stated time limit, and re-runs only the files that failed.

### Added
- **`scripts/pick_tests.py` (TEST.722).** Without arguments it lists, in priority order, the test files the PR edits, then those whose name holds a changed module's name, then those that import a changed module (fewest imports of it first). `--run [--minutes N]` runs them in small batches and starts no batch after the limit (default 20 minutes), then names the files "left to the end-of-version run"; `--failed` (or `--failed-from LOG` for a pytest log) runs only the files that failed in the last run, so a fix is proved by its failed files and not by the whole set again. The pre-gate set (`scripts/data/pregate_set.txt`) is left out because `scripts/pregate.py` runs it. Reused: `pregate.changed_files` and `pregate._lines`; the gate's by-hand name matching becomes the second rule. `src/tests/test_pick_tests.py` pins the module names, the import forms, the order, the limit and the failed-file reading.

generation-output: unchanged
generation-affects: none
