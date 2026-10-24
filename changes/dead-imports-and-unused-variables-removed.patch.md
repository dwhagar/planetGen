**Headline:** 56 unused imports, 12 duplicate imports and 13 unused variables are gone from the code, and a loop variable no longer hides the `html` module in the TODO page builder.

### Fixed
- **Dead code that `ruff` (rule set F) and `vulture` reported (TEST.364).** Removed: 56 imports nothing used (for example `CAMERA_START_THETA_DEG` in `web/warmup.py`, `to_paragraph` in `generation/planet.py`), the physics imports `db/store.py` listed twice, and 13 local variables that were set and never read (for example `extents` and `e_value` in `generation/run_plan.py`, `by_node` in `physics/system_tree.py`). Before each import went, the source and the tests were searched for anything that imported it from that module or patched it by name: nothing did. `ruff check --select F` is clean on `src/planetgen` and `scripts`.
- **`scripts/build_todo_docs.py` reused the name `html` for a loop variable,** hiding the `html` module it imports; the loop now says `page_html`.
TODO items resolved: TEST.364. Added: none.
