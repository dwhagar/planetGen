**Headline:** Only the TODO thread regenerates the global documents, and a stale one no longer fails a lane's tests.

### Changed
- **The global documents are TODO's (DOC.69; Boss, 2026-10-18, 11:24 PM EDT).** The architecture map, the plan pages and the docs index are built from the whole repository, so any lane's merge made them stale and every lane regenerated them. `python scripts/build_global_docs.py` now rebuilds all three in order (`--check` names the stale ones); the TODO thread runs it after a batch of merges, and a lane's pull request never contains a regenerated global file.
- **The tests that compare each with its source are marked `globaldocs` and left out of every run** unless `PLANETGEN_GLOBAL_DOCS=1` (`scripts/test.py <command> --global-docs` sets it): the architecture map's module list (two tests), the docs index (three), and the plan pages (one). The release stamp runs `scripts/test.py full --check --global-docs`. Tests of the pages' own rules (a description on every page, a link back to the index, what the generators produce) still run.
- **Two stale docstrings now name `planetgen warm`** (`admin/maintenance.py`, `cli/__init__.py`), and `docs/testing.md` and `docs/design/architecture.md` describe the rule.

TODO items resolved: DOC.69. Added: none.
