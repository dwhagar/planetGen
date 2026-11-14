# Old galaxy schemas (TEST.8, TEST.9)

`schema_v<N>.sql.gz` is `src/planetgen/db/schema.sql` as `main` shipped it
while `SCHEMA_VERSION` was N, with its `--` comment lines and blank lines
removed. Only the Alembic baseline (v130, DB.98) and the newest three versions below the current one are kept (TEST.729, Boss's
hard cut of 2026-11-12: no database structure before 8.0): the rest went, with the migration tests that loaded them.
`tests/test_db_old_schemas.py` loads each into an empty database, migrates it, and compares the result with a new one.
A database older than the baseline is brought up with `legacy/bump_database_to_8_0.py` (DB.126).

When a schema change adds a revision, add the previous version's file here
(`git show <commit before the change>:src/planetgen/db/schema.sql`, comments
and blank lines removed), so the revision is tested against the schema it
starts from, and delete the oldest file above the baseline when more than three are kept.

Commit each file came from:

- v130: b9bf305 (DB.85 galaxy containers) (before OPS.58)
- v243: after GEN.1497's migration 0243 (body history table), before GEN.1314 class fudge columns (0244)
- v244: after GEN.1314's migration 0244 (class fudge columns), before GEN.1425 cloud_chemistry column (0245)
- v245: after GEN.1425's migration 0245 (cloud_chemistry column), before GEN.1498 system_summary ledger_stamp column (0246)
- v246: after GEN.1498's migration 0246 (ledger_stamp column), before GEN.1535 obliquity range columns (0247)
- v247: after GEN.1535's migration 0247 (obliquity range columns), before GEN.1510 planets belt_order and belt_rank columns (0248)
- v248: after GEN.1510's migration 0248 (planets belt_order and belt_rank columns), before GEN.1567 ice cover column (0249)
- v249: after GEN.1567's migration 0249 (ice cover column), before GEN.558 rogue_planets.has_moons dropped (0250)
- v250: after GEN.558's migration 0250 (rogue_planets.has_moons dropped), before GEN.1267 drop of the nebula shape table and columns (0251)
- v251: after GEN.1267's migration 0251 (nebula shape dropped), before GEN.1201 shape_temperature_k column (0252)
- v252: after GEN.1201's migration 0252 (shape_temperature_k column), before OPS.180 chart ledger repair (0253, no schema change)
