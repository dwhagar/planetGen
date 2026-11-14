# legacy/

Scripts that bring data made by an older planetGen to the current one. Boss's hard cut (2026-11-12): **no
generation and no database structure before 8.0 is accepted by 8.0 and above**; the code refuses an old database
and points here.

Rules for everything in this folder:

- Run by hand, on purpose, by someone who has data they want to keep. `update.sh`, the site, the queue and the tests
  never run or import anything here, and no test covers it (Boss: the script only has to exist).
- One script per job, nothing else in it. A new job gets a new file and a line in the table below.
- Back the data up first.

| Script | What it does | Item |
| --- | --- | --- |
| `convert_to_8_0.py` | Converts the data of the last release before 8.0 to the 8.0 form (runs the database bump, notes that `config.json` needs nothing). | OPS.177 |
| `bump_database_to_8_0.py` | Brings a content database at schema v61 to v129 to v130, the oldest structure 8.0 and above accept, with the revisions in `db_revisions/` (copies of the steps `planetgen/db/migrations/` held before the re-baseline, DB.98). From v130 `update.sh` carries it on. | DB.126 |
| `db_revisions/` | The Alembic revisions v61 to v130 that the bump runs. Not part of the live migration chain. | DB.126 |
