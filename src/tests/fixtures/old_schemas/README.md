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
- v72: df5df13 (GEN.86 activity and magnetism, before GEN.185)
- v73: 7c06976 (GEN.185 star scatter mass limit, before GEN.88)
- v74: ea2e00e (GEN.88 hydrosphere, before GEN.87)
- v75: 62f3842 (GEN.87 surface radiation, before GEN.89)
- v76: 3668b58b (GEN.89 habitability score, before GEN.187)
- v77: 9a023637 (GEN.187 scatter classes, before DB.20)
- v78: 3d44328c (object IDs, before GEN.187)
- v79: 5f09bea2 (GEN.187 mass backfill, before MAP.165)
- v80: dc510bb4 (MAP.165 scatter mass, before ADM.49)
- v81: ff1be90c (ADM.49 shape density, before PERF.68/70/74)
- v82: d1920845 (PERF.68/70/74 slow-read indexes, before DB.24)
- v83: fa95544a (DB.24 nebula_field, before DB.23)
- v84: 380bdc3e (DB.23 map luminosity, before PERF.92)
- v85: 04b69361 (PERF.92 object timing, before DB.25)
- v86: ac6330ca (ADM.54 Redis health, before DB.26)
- v87: 9c725da0 (DB.26 habitability columns, before the POP.9 columns)
- v88: 29f7a086 (DB.27 species tech columns, before DB.28)
- v89: c6c22061 (DB.28 facility types, before the MAP.151b ladder_step)
- v90: b73bb27 (MAP.151b block ladder, before DB.32 spherical sectors)
- v91: c8fec05c (DB.32 spherical sectors, before the MAP.205 block totals on shells)
- v92: 0980a8e (MAP.205 block totals on shells, before GEN.230/GEN.233/GEN.223 96-bit IDs)
- v93: 1936db3 (GEN.230/GEN.233/GEN.223 96-bit IDs, before GEN.234 halo populations)
- v94: a91a7156 (GEN.234 halo populations, before GEN.158 metallicity_fe_h)
- v95: 74ff7e3 (GEN.158 metallicity_fe_h, before GEN.153 magnetar pulsar type)
- v96: 5583b91d (GEN.153 magnetar pulsar type, before GEN.159 galaxy_components)
- v97: fd18454 (GEN.159 galaxy_components, before GEN.281 cell star ceiling and fill luminosity cap)
- v98: d320eba5 (GEN.281 cell star ceiling and fill luminosity cap, before GEN.280 nuclear structures)
- v99: 47705939 (GEN.280 nuclear structures, before GEN.281 cell star ceiling default and CHECK)
- v100: 21884c05 (GEN.364 cell star ceiling 5,000 and its CHECK, before DB.54 scatter limits by place)
- v101: a438b42f (DB.54 scatter limits by place, before DB.55 point columns)
- v102: the schema of migration 0102 (DB.55 point columns, before GEN.336 movement vector); built in one change with v103 and v104, so no commit of main carries it alone
- v103: the schema of migration 0103 (GEN.336 movement vector, before DB.53 orbit nodes); same change, no commit of main carries it alone
- v104: 3940ad23 (DB.55, GEN.336 and DB.53 together as #1952, before GEN.324 satellite stars)
- v105: 0690c541 (GEN.324 satellite stars, before DB.65 cell key edge)
- v106: 9c172019 (DB.65 cell key edge, before GEN.526 satellite young and intermediate populations)
- v107: 9a8f7671 (GEN.526 satellite young and intermediate populations, before DB.72 band_index INT)
- v108: 039042e (DB.72 band_index INT, before GEN.289 planet state vectors)
- v109: 6d0c7d3 (GEN.289 planet state vectors, before DB.50 life structure columns)
- v110: a8946a3 (DB.50 life structure columns, before DB.52 drops the prose; taken from the branch, as no commit of main carries it alone)
- v111: 2f34d931 (DB.52 drops the prose, before GEN.244 open cluster populations)
- v112: 422c052 (GEN.244 open cluster populations, before DB.79 object index)
- v113: b4d4e8e (DB.79 object index, before DB.81 shape pins)
- v114: 48e8dfe9 (DB.81 shape pins, before DB.44 sector storage checksum)
- v115: 45824355 (DB.44 sector storage checksum, before MAP.172 nebula block aggregates)
- v116: 836b3d8e (MAP.172 nebula block aggregates, before DB.66 family columns)
- v117: be4a2e3 (DB.66 galaxy family columns, before DB.56 orbital events)
- v118: ce94274 (DB.56 orbital events) (before DB.30)
