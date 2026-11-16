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
- v132: c42411611 (GEN.458 snow line) (before GEN.630)
- v133: 6f0754dda (GEN.630 star members) (before GEN.267)
- v134: e15c1bc (GEN.267 star count settings, before GEN.413)
- v135: 9f96ca8 (GEN.413 ozone column) (before GEN.655)
- v136: 11e5784 (GEN.655 survival equipment, before GEN.656)
- v137: ce50ec2 (GEN.656 regional temperatures, before GEN.414)
- v139: 3daf30ba3 (before DB.73)
- v140: 2ae3043 (DB.73 parent link) (before GEN.459)
- v146: d68eabf (MAP.258 block ladder root) (before MAP.259 block lights)
- v150: e229208 (GEN.649 life scores, before GEN.461)
- v152: d0c4d6c (GEN.700 moon tide flux, before GEN.699)
- v160: cca3b2b03 (GEN.732 rings, before DB.107 region stats)
- v165: a03f5260 (MAP.305 system interests, before MAP.321 uninhabitable worlds and asteroid fields)
- v167: 766381158 (MAP.321 inhabited worlds, before DB.108 sector mass)
- v168: 8c67694d5 (DB.108 sector mass, before GEN.721 lid regime)
- v169: GEN.721 lid regime as built in this change (the v168 schema plus `lid_regime`, before GEN.767 day and night shares)
- v173: 849970d3 (MAP.343 and later, before MAP.344 civilization technology)
- v174: 95652091e (MAP.344 civilization technology, before GEN.889 life class switch)
- v177: aebe4cab (PERF.205 species sort indexes, before PERF.205 phenomenon name indexes)
- v178: 0d3e81011 (PERF.205 phenomenon name indexes, before GEN.908 brine seep)
- v179: 0d3e81011 (GEN.908 brine seep, before DB.70 moon origin and giant disk)
- v180: 50dfa10d9 (DB.70 moon origin and giant disk, before PERF.208 sector sort keys)
- v185: 452eb9f1 (GEN.846 star rotation period, before PERF.211 phenomenon sort indexes)
- v186: 26009638e (PERF.211 phenomenon sort indexes, before GEN.768 hours outside in a day)
- v187: train 5 (GEN.768 hours outside in a day, before GEN.1035 surface liquid; the train-5 merge before the liquid branch)
- v189: train 5 end (cac3113d5, GEN.1009 member record, before the train 6 ground intermediate)
- v190: train 6 (GEN.474 ground intermediate, before GEN.1034 ocean tides)
- v192: train 6 (GEN.1057 ring age, before GEN.1097 day and night state)
- v194: train 7 (before GEN.983 belt dose)
- v196: train 7 (GEN.1108 latitude bands, before GEN.1215 surface liquid solutes; built by Tides 1, the Train Manager may replace it)
- v204: train 9 (belt origin, before GEN.1218b surface liquid mix; the schema of main after train 9, built by Tides 1, the Train Manager may replace it)
- v210: train 10 (GEN.1232a belt history text, before DOC.103 American spelling; the train schema before the sweep)
- v213: train 12 after GEN.1254's migration 0213 (scatter Object IDs), before GEN.1276 nebula composition (0214)
- v214: train 12 after GEN.1276's migration 0214 (nebula composition), before GEN.1277 nebula field composition (0215)
- v216: train 12 after VIEW.25's migration 0216 (sky tier draws), before GEN.1265 temperature coverage (0217)
- v220: train 12 after TEST.602's migration 0220 (sector content leaf), before GEN.1256 Object ID names (0221)
- v221: train 13 after GEN.1256's migration 0221 (Object ID names, no column change), before GEN.1257 charted reason (0222)
- v231: train 18 after GEN.1364's migration 0231 (core composition), before GEN.1376 dust species (0232)
- v234: train 18 after GEN.1380's migration 0234 (ore deposit tables), before GEN.1377 and GEN.1420 isotope record tables (0235)
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
