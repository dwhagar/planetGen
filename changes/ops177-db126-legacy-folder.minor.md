**Headline:** A new legacy/ folder holds the hand-run scripts that bring data made before 8.0 to 8.0, starting with the database bump.

### Changed
- **OPS.177, DB.126: `legacy/` (Boss's hard cut, 3:53 PM EST 2026-11-12: no generation or database structure before 8.0 is accepted by 8.0 and above).** `legacy/convert_to_8_0.py` converts the data of the last release before 8.0 (it runs the database bump; `config.json` needs nothing). `legacy/bump_database_to_8_0.py` takes a content database at schema v61 to v129 to v130, the oldest structure the live migrations accept, with the Alembic revisions v61 to v130 kept in `legacy/db_revisions/` (copies of the steps removed by the re-baseline, DB.98). Run by hand only; `update.sh`, the site and the tests never call them. Checked once by hand on a v67 schema fixture from before the re-baseline: 63 steps to v130, then `planetgen.cli.migrate` carried it to the current head. `scripts/american_spelling_sweep.py` skips `legacy/db_revisions/` as it skips the live migrations (stored words stay as they were). No generated output moves.

generation-affects: none
