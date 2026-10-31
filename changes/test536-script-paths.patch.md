**Headline:** A test that checks the deploy scripts name only files that exist no longer fails on a comment in the wiring guard.

### Changed

- `scripts/wiring_guard.py` names `planetgen/db/kind_tables.py`, `planetgen/db/models.py` and `planetgen/db/store.py` in full in its text; it wrote `db/kind_tables.py` and so on, which the deploy-script path test read as repo paths that do not exist (the files live under `src/planetgen/`). Behaviour is unchanged. generation-output: unchanged
