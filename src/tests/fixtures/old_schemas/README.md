# Old galaxy schemas (TEST.8, TEST.9)

`schema_v<N>.sql.gz` is `src/planetgen/db/schema.sql` as `main` shipped it
while `SCHEMA_VERSION` was N, with its `--` comment lines and blank lines
removed. Only the Alembic baseline is kept: the older ones went with the
legacy migration steps (DB.11). `tests/test_db_old_schemas.py` loads each
into an empty database, migrates it, and compares the result with a new one.

When a schema change adds a revision, add the previous version's file here
(`git show <commit before the change>:src/planetgen/db/schema.sql`, comments
and blank lines removed), so the revision is tested against the schema it
starts from.

Commit each file came from:

- v61: 4c029efa (the DB.11 baseline, with the `alembic_version` table)
- v62: 3a2e936a (GEN.123 sector paths, before DB.13)
- v63: 70f7fa0
- v64: 3b9bb5f0 (GEN.125 facility velocity, before GEN.100)
- v65: 15d3f7a (GEN.100 phenomenon scatter, before GEN.106)
- v66: eebcfb3b (GEN.106 update clocks, before DB.7)
- v67: 500facd (DB.7 sector version, before GEN.104)
- v68: de8d3a7 (GEN.104 spin, before GEN.137)
- v69: 222288d (GEN.137 scatter epoch, before GEN.85)
- v70: 1dece90 (GEN.85 atmosphere species, before GEN.167)
- v71: 376d120 (GEN.167 phenomenon mass cut, before GEN.86)
