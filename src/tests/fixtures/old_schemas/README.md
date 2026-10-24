# Old galaxy schemas (TEST.8, TEST.9)

`schema_v<N>.sql.gz` is `src/planetgen/db/schema.sql` as `main` shipped it
while `SCHEMA_VERSION` was N, with its `--` comment lines and blank lines
removed. Only the Alembic baseline (v130, DB.98) and the versions after it
are kept: the older ones went when the migration history was re-baselined. `tests/test_db_old_schemas.py` loads each
into an empty database, migrates it, and compares the result with a new one.

When a schema change adds a revision, add the previous version's file here
(`git show <commit before the change>:src/planetgen/db/schema.sql`, comments
and blank lines removed), so the revision is tested against the schema it
starts from.

Commit each file came from:

- v130: b9bf305 (DB.85 galaxy containers) (before OPS.58)
- v131: d10849e (OPS.58 generator epoch) (before GEN.458)
- v132: c42411611 (GEN.458 snow line) (before GEN.630)
- v133: 6f0754dda (GEN.630 star members) (before GEN.267)
- v134: e15c1bc (GEN.267 star count settings, before GEN.413)
- v135: 9f96ca8 (GEN.413 ozone column) (before GEN.655)
- v136: 11e5784 (GEN.655 survival equipment, before GEN.656)
- v137: ce50ec2 (GEN.656 regional temperatures, before GEN.414)
- v138: f65a213 (before MAP.143)
- v139: 3daf30ba3 (before DB.73)
- v140: 2ae3043 (DB.73 parent link) (before GEN.459)
- v141: 7cbd9d6 (GEN.459 radiogenic module, before GEN.412)
- v142: 77fdfbb96 (GEN.412 albedo and greenhouse, before GEN.666)
