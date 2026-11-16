# Old galaxy schemas (TEST.8, TEST.9)

`schema_v<N>.sql.gz` is `src/planetgen/db/schema.sql` as `main` shipped it
while `SCHEMA_VERSION` was N, with its `--` comment lines and blank lines
removed. Only the Alembic baseline (v130, DB.98) and the versions after it
are kept: the older ones went when the migration history was re-baselined. `tests/test_db_old_schemas.py` loads the baseline and the newest three
into an empty database, migrates each, and compares the result with a new one (TEST.729: the run from the baseline passes through every revision, and the fixtures hold no rows, so the
ones between add nothing). A fixture stays only while a test loads it by number (`load_old_schema(config, N)`): when such a test is retired, its fixture goes with it,
unless it is the baseline or one of the newest three. The commit lines below are for the fixtures that are left.

When a schema change adds a revision, add the previous version's file here
(`git show <commit before the change>:src/planetgen/db/schema.sql`, comments
and blank lines removed), so the revision is tested against the schema it
starts from.

Commit each file came from:

- v130: b9bf305 (DB.85 galaxy containers) (before OPS.58)
- v131: d10849e (OPS.58 generator epoch) (before GEN.458)
- v135: 9f96ca8 (GEN.413 ozone column) (before GEN.655)
- v136: 11e5784 (GEN.655 survival equipment, before GEN.656)
- v140: 2ae3043 (DB.73 parent link) (before GEN.459)
- v150: e229208 (GEN.649 life scores, before GEN.461)
- v165: a03f5260 (MAP.305 system interests, before MAP.321 uninhabitable worlds and asteroid fields)
- v173: 849970d3 (MAP.343 and later, before MAP.344 civilization technology)
- v186: 26009638e (PERF.211 phenomenon sort indexes, before GEN.768 hours outside in a day)
- v192: train 6 (GEN.1057 ring age, before GEN.1097 day and night state)
- v213: train 12 after GEN.1254's migration 0213 (scatter Object IDs), before GEN.1276 nebula composition (0214)
- v214: train 12 after GEN.1276's migration 0214 (nebula composition), before GEN.1277 nebula field composition (0215)
- v220: train 12 after TEST.602's migration 0220 (sector content leaf), before GEN.1256 Object ID names (0221)
- v242: after GEN.1484's migration 0242 (light_factor column), before GEN.1497 body history table (0243)
- v243: after GEN.1497's migration 0243 (body history table), before GEN.1314 class fudge columns (0244)
- v245: after GEN.1425's migration 0245 (cloud_chemistry column), before GEN.1498 system_summary ledger_stamp column (0246)
- v246: after GEN.1498's migration 0246 (ledger_stamp column), before GEN.1535 obliquity range columns (0247)
- v259: after GEN.1433's migration 0259 (flavor_key columns), before GEN.332 slice D's rogue_planets pair columns (0260) (train 27's schema.sql with its comment and blank lines removed)
- v260: after GEN.332 slice D's migration 0260 (rogue_planets pair columns), before OPS.215 engine_items (0261) (train 28's schema.sql with its comment and blank lines removed)
- v263: after OPS.219's migration 0263 (engine_piece_results), before DB.112's migration 0264 (sector_system_counts loses quadrant) (main's schema.sql with its comment and blank lines removed)
- v264: after DB.112's migration 0264 (sector_system_counts loses quadrant), before GEN.1141's migration 0265 (designations reason `deleted`) (main's schema.sql with its comment and blank lines removed)
- v265: after GEN.1141's migration 0265 (designations reason `deleted`), before GEN.1753's migration 0266 (body_history birth numbers, lost_moons) (main's schema.sql with its comment and blank lines removed)
- v266: after GEN.1753's migration 0266 (body_history birth numbers, lost_moons), before GEN.1231's migration 0267 (belt thrown inward) (main's schema.sql with its comment and blank lines removed)
- v267: after GEN.1231's migration 0267 (belt thrown inward), before DB.132 and DB.133's migration 0268 (stars birth mass and phase end, body_history held and forced eccentricity) (main's schema.sql with its comment and blank lines removed)
