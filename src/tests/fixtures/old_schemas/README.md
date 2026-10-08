# Old galaxy schemas (TEST.8, TEST.9)

`schema_v<N>.sql.gz` is `src/planetgen/db/schema.sql` as `main` last
shipped it while `SCHEMA_VERSION` was N, with its `--` comment lines and
blank lines removed. Versions that never reached `main` (10-13, 16) are
missing. `tests/test_db_old_schemas.py` loads each into an empty database,
migrates it, and compares the result with a fresh one.

Commit each file came from:

- v8: 72741e0
- v9: 25e3b6b
- v14: 2ed84be
- v15: 41c8775
- v17: 3d15324
- v18: 51ecf21
- v19: a7f431c
- v20: 87451a9
- v21: cef1069
- v22: ef59195
- v23: 5447db1
- v24: bea66ed
- v25: 9844b59
- v26: 1bb1bec
- v27: 2af6694
- v28: bc20f5c
- v29: fe4699c
- v30: faffaa7
- v31: 75e21be
- v32: c1e2aa4
- v33: 4a3ddde
- v34: ef3be19
- v35: 1ccc549
- v36: 317cb9d
- v37: 7f7d181
- v38: fd3bb7d
- v39: f5b7d54
- v40: 2e34748
- v41: 125449b
- v42: 0b8edd2
- v43: 329380f
- v44: 1bbd586
- v45: f8c456d
- v46: c0442d2
- v47: 65633a7
- v48: c414509
- v49: 5b77021
- v53: 85c0415
- v54: ba21569
- v55: 3284654
- v56: 7a46269
- v57: 5528f58

To add the schema a release is leaving behind (run from the repo root,
before bumping `SCHEMA_VERSION`):

    N=$(grep -m1 '^SCHEMA_VERSION = ' src/stellarObjects/store.py | awk '{print $3}')
    grep -v '^\s*--' src/planetgen/db/schema.sql | sed 's/\s\+--.*$//' | grep -v '^\s*$' \
        | gzip -9n > src/tests/fixtures/old_schemas/schema_v$N.sql.gz
