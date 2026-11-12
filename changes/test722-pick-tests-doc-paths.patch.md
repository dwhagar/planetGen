**Headline:** `scripts/pick_tests.py` no longer names two example repo paths that do not exist, which turned `test_scripts_only_name_repo_paths_that_exist` red on main.

### Fixed
- **`scripts/pick_tests.py` (TEST.722).** Two docstring examples wrote made-up repository paths (`src/planetgen/a/b/c.py`, `src/tests/test_x.py`); the deploy-scripts test that checks every path a script names found them. They are reworded so the docstring names no path. No behaviour changed.

generation-output: unchanged
generation-affects: none
