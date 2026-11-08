# Items: DB, Database

The life of every DB item that a commit title names, from the git history as of `8487fa43b`: first and last date, whether it is open in [docs/TODO.md](../TODO.md) now, and its commit titles in order (the pull request numbers are the GitHub pull requests). An item not open in TODO.md has shipped, been closed as a finding or been dropped; [dropped-and-closed.md](dropped-and-closed.md) lists the closings.

## DB.2

2026-10-01 to 2026-10-01; 6 commits; not open.

- 2026-10-01: Docs: file DB.2-4, bugs found by the Database tests thread
- 2026-10-01: Docs: move the DB.2-4 rows to the ID table
- 2026-10-01: TODO: Boss's decisions on MAP.81, DB.2, GEN.46 and API.6; file API.15 (API call log) in phase 1 (#336)
- 2026-10-02: DB.2: read asteroid field and interstellar comet composition rows on the page and in the API
- 2026-10-02: DB.2: the composition rows are read now; drop them from the never-read list
- 2026-10-01: TODO: retire DB.2 and DB.3 (done in PR #347) (#348)

## DB.3

2026-10-02 to 2026-10-01; 2 commits; not open.

- 2026-10-02: DB.3: resetDb keeps id_blocks so ids are never reused across a reset
- 2026-10-01: TODO: retire DB.2 and DB.3 (done in PR #347) (#348)

## DB.4

2026-10-02 to 2026-10-01; 2 commits; not open.

- 2026-10-02: DB.5, DB.4: schema creation under a per-database lock; detect version from table shape when schema_migrations is empty
- 2026-10-01: TODO: NAV.12 is Boss's routing guarantee (phase 0) split into NAV.35-37 and UX.35; version and arc size decided; retire DB.4 and DB.5 (PR #342); file TEST.78 (#345)

## DB.5

2026-10-01 to 2026-10-01; 3 commits; not open.

- 2026-10-01: Docs: file DB.5, schema creation race on first connections (bug)
- 2026-10-02: DB.5, DB.4: schema creation under a per-database lock; detect version from table shape when schema_migrations is empty
- 2026-10-01: TODO: NAV.12 is Boss's routing guarantee (phase 0) split into NAV.35-37 and UX.35; version and arc size decided; retire DB.4 and DB.5 (PR #342); file TEST.78 (#345)

## DB.6

2026-10-01 to 2026-10-02; 4 commits; not open.

- 2026-10-01: TODO: version + seed reproduce the same galaxy (GEN.55 chain: DB.6, OPS.10-12, GEN.56-58, TEST.77, ADM.17, API.16) (#340)
- 2026-10-02: TODO: retire GEN.39 (DB.6 and OPS.10 stand alone, schema v51); file TEST.86 (#383)
- 2026-10-02: DB.6: the galaxy records the version key that made it, and every run (schema v52) (#387)
- 2026-10-02: TODO: retire DB.6 (PR #387); DB.7, OPS.10, OPS.13 point at versionKey.version_key; schema v52 (#388)

## DB.7

2026-10-01 to 2026-10-09; 4 commits; not open.

- 2026-10-01: TODO: merge the reproducibility plan (DB.7, API.17, GEN.59; OPS.12 to phase 4) and revise NAV.12 from the hop-length study (NAV.34 islands bug) (#341)
- 2026-10-02: TODO: retire DB.6 (PR #387); DB.7, OPS.10, OPS.13 point at versionKey.version_key; schema v52 (#388)
- 2026-10-09: Each sector records its generating version; mixed-version warning (DB.7) (#813)
- 2026-10-09: Retire DB.7, NAV.10, GEN.108, ADM.18 and GEN.107 (PRs #813-#817) (#818)

## DB.8

2026-10-01 to 2026-10-09; 3 commits; not open.

- 2026-10-01: TODO: add DB.8 consistency check, DB.9 parity repair, ADM.18 galaxy settings JSON (#352)
- 2026-10-09: Add planetgen check-db and a Check the database job (DB.8) (#893)
- 2026-10-09: Retire DB.8 (PR #893) (#894)

## DB.9

2026-10-01 to 2026-10-28; 2 commits; not open.

- 2026-10-01: TODO: add DB.8 consistency check, DB.9 parity repair, ADM.18 galaxy settings JSON (#352)
- 2026-10-28: Retire DB.48 and DB.9, file GEN.847 and DB.109 (#3461)

## DB.10

2026-10-01 to 2026-10-01; 3 commits; not open.

- 2026-10-01: TODO: daily maintenance run and 18-slot JSON backups (GEN.61, OPS.16-18, ADM.19, DB.10); admin changes as pending deltas (#358)
- 2026-10-01: TODO: DB.10 and OPS.12's backup/pending parts low priority (Boss 02:31Z) (#359)
- 2026-10-01: TODO: ADM.20 merge-now button (phase 3, low priority); undo the low-priority marks on DB.10 and OPS.12 (#360)

## DB.11

2026-10-08 to 2026-10-09; 5 commits; not open.

- 2026-10-08: Use SQLAlchemy's QueuePool for database connections instead of DBUtils (DB.11 step 1) (#734)
- 2026-10-08: Add Alembic with a v61 baseline; revisions are schema versions (DB.11 step 2) (#744)
- 2026-10-08: Add SQLAlchemy table models generated from schema.sql, with a match test (DB.11 step 3) (#745)
- 2026-10-09: Remove the legacy migration steps and old-schema fixtures; Alembic carries every migration from v61 (DB.11) (#751)
- 2026-10-09: TODO: retire DB.11 (PR #751) and update the architecture flow for Alembic (#752)

## DB.13

2026-10-09 to 2026-10-09; 2 commits; not open.

- 2026-10-09: Store run arguments and job argv in columns; add missing indexes (DB.13) (#766)
- 2026-10-09: TODO: retire DB.13 (PR #766); phase 0 is complete, major note for 8.0 (#773)

## DB.14

2026-10-07 to 2026-10-08; 3 commits; not open.

- 2026-10-07: File sector and block coloring items DB.14 and MAP.128-MAP.132 (Boss 2026-10-08) (#527)
- 2026-10-08: Color sectors by star age, density and luminosity from raw sector stats (DB.14, MAP.128-130) (#556)
- 2026-10-08: Retire DB.14, MAP.129, MAP.130 (PR #556); narrow MAP.128 to the Sector Map (#557)

## DB.15

2026-10-09 to 2026-10-09; 3 commits; not open.

- 2026-10-09: Progress bars: one mechanism and the 15 second rule (UX.84, PERF.51); PERF.33 and DB.15 into Phase 1 (#898)
- 2026-10-09: Migration progress: step bar, batch bar, plain lines off a terminal, lock wait (DB.15) (#924)
- 2026-10-09: TODO: retire DB.15 (PR #924) (#925)

## DB.16

2026-10-22 to 2026-10-22; 1 commits; not open.

- 2026-10-22: DB.16: sectors record the generator epoch and run id instead of four version text columns (migration 0121)

## DB.17

2026-10-09 to 2026-10-09; 1 commits; open.

- 2026-10-09: Record Boss's answers on API.17, OPS.16, DB.17 and the seed field (#858)

## DB.18

2026-10-16 to 2026-10-16; 2 commits; not open.

- 2026-10-16: DB.18: migration helpers for slow DDL (instant columns, online indexes, batched updates) (#1901)
- 2026-10-16: TODO: retire DB.18 (#1901); browser-hang sighting on TEST.193/195; DB.65 revision 0105 (#1902)

## DB.19

2026-10-09 to 2026-10-09; 5 commits; not open.

- 2026-10-09: File the generation performance study: PERF.42 to PERF.48 and DB.19 (#836)
- 2026-10-09: Add the phenomenon scatter mass-cut note: recommended cut, error tables, sector-time fill, DB.19 effect (#837)
- 2026-10-09: File the phenomenon scatter mass cut (GEN.166-GEN.169), rework DB.19, fix the 1.6e8 row figures (#839)
- 2026-10-09: Retire GEN.166 to GEN.168; DB.19 becomes conditional (#867)
- 2026-10-09: Retire OPS.14; close DB.19 (#877)

## DB.20

2026-10-09 to 2026-10-10; 3 commits; not open.

- 2026-10-09: File the object ID items (GEN.170 to GEN.176, DB.20, API.23, TEST.110, DOC.5) (#874)
- 2026-10-10: DB.20, GEN.171: object IDs in the schema and the sector fill (schema v78) (#1055)
- 2026-10-10: Docs: retire DB.20 and GEN.171 (PR #1055) (#1056)

## DB.21

2026-10-09 to 2026-10-10; 3 commits; not open.

- 2026-10-09: File DB.21: deep pass for the database check (#897)
- 2026-10-10: Deep database check validates every star system, with a time estimate (DB.21) (#1082)
- 2026-10-10: TODO: retire DB.21 (PR #1082) (#1084)

## DB.22

2026-10-10 to 2026-10-10; 2 commits; not open.

- 2026-10-10: TODO: retire API.22, API.23 stage 1 done, file DB.22 (PR #1119) (#1120)
- 2026-10-10: TODO: retire DB.22, API.23 stage 2 done (PR #1127) (#1128)

## DB.23

2026-10-10 to 2026-10-14; 2 commits; not open.

- 2026-10-10: TODO: file map luminosity items (GEN.200, DB.23, GEN.201, MAP.171) and PERF.73, retire PERF.63 (PR #1174) (#1176)
- 2026-10-14: TODO: retire PERF.93, DB.23, GEN.201; file TEST.147 (#1366)

## DB.24

2026-10-10 to 2026-10-14; 3 commits; not open.

- 2026-10-10: TODO: file the 14 nebula visibility items (DB.24, MAP.172-182, GEN.202, TEST.135) (#1183)
- 2026-10-14: The nebula_field table: one row per kept class M and N field cloud (DB.24) (#1317)
- 2026-10-14: TODO: retire DB.24, OPS.43, OPS.44, ADM.53 (#1317, #1307, #1318) (#1319)

## DB.25

2026-10-14 to 2026-10-14; 2 commits; not open.

- 2026-10-14: TODO: file DB.25 and TEST.145; MAP.151 notes #1335 (#1336)
- 2026-10-14: TODO: retire DB.25 (#1376); schema v86; MAP.151b unblocked (#1377)

## DB.26

2026-10-14 to 2026-10-14; 5 commits; not open.

- 2026-10-14: Add DB.26: one migration for the habitability columns (#1400)
- 2026-10-14: Retire ADM.54 and add the induced moment to DB.26 (#1402)
- 2026-10-14: DB.26 adds the nebula columns; GEN.150 depends on it (#1403)
- 2026-10-14: Habitability, flux, induced-moment and nebula columns (DB.26, schema v87) (#1423)
- 2026-10-14: Retire DB.26 and OPS.50; add DB.27 and DB.28 migrations for POP.9 and POP.10 (#1426)

## DB.27

2026-10-14 to 2026-10-14; 3 commits; not open.

- 2026-10-14: Retire DB.26 and OPS.50; add DB.27 and DB.28 migrations for POP.9 and POP.10 (#1426)
- 2026-10-14: Species tech-level columns (DB.27, schema v88) (#1446)
- 2026-10-14: TODO: retire DB.27 (PR #1446) (#1449)

## DB.28

2026-10-14 to 2026-10-15; 4 commits; not open.

- 2026-10-14: Retire DB.26 and OPS.50; add DB.27 and DB.28 migrations for POP.9 and POP.10 (#1426)
- 2026-10-15: TODO: DB.28 and DB.29 sit with Phase 1 (#1465)
- 2026-10-15: Facility types and facility columns (DB.28, schema v89) (#1476)
- 2026-10-15: TODO: retire DB.28 (PR #1476) (#1477)

## DB.29

2026-10-14 to 2026-10-29; 7 commits; not open.

- 2026-10-14: Add DB.29 migrations refuse to run while jobs are active (#1435)
- 2026-10-14: Widen DB.29: update.sh refuses to start while jobs are active without --yes or confirmation (#1437)
- 2026-10-14: DB.29: systemd or cron runs and 30-second silence mean no (#1438)
- 2026-10-15: TODO: DB.28 and DB.29 sit with Phase 1 (#1465)
- 2026-10-28: Docs: retire TEST.211 (#3625), give DB.29 to Core 4 as its lock-free item (#3627)
- 2026-10-29: DB.29: update.sh and migrations refuse to start while jobs are active unless --yes or a typed yes (#3645)
- 2026-10-29: Docs: retire DB.29 (#3645), ADM.58 to Core 4, PERF.100 note (#3646)

## DB.30

2026-10-15 to 2026-10-22; 3 commits; not open.

- 2026-10-15: TODO: add DB.30, POP.12, POP.13 (species environment limits, could-live-here list, species editor) (#1463)
- 2026-10-22: DB.30: per-species environment limits table (migration 0119)
- 2026-10-22: docs: DB galaxy-schema items DB.30 to DB.63 retired (#2520 to #2530, Alembic head 0129), GEN.379 retired (#2531), DB.85 proposed for the DB lane (#2532)

## DB.31

2026-10-21 to 2026-10-21; 1 commits; not open.

- 2026-10-21: Retire DB.31 (umbrella, subitems done in #2206); fix the control version wording on DB.62 (#2393)

## DB.32

2026-10-15 to 2026-10-15; 6 commits; not open.

- 2026-10-15: DB.32: convert the galaxy in place, keeping every position (#1541)
- 2026-10-15: GEN.224, GEN.225, DB.32: no undefined space (#1547)
- 2026-10-15: TODO: record object ID renumbering in the DB.32 conversion (#1560)
- 2026-10-15: TODO: DB.32 merged into the switch branch (#1574) (#1575)
- 2026-10-15: Spherical sector coordinates: the switch to shells, bands and wedges (DB.32, GEN.227, UX.111, UX.112) (#1607)
- 2026-10-15: TODO: retire DB.32, GEN.226, GEN.227, UX.111, UX.112 (spherical switch #1607) (#1608)

## DB.33

2026-10-15 to 2026-10-22; 9 commits; not open.

- 2026-10-15: TODO: file the test failures from Quality's 4:56 PM EDT run (DB.33, TEST.163; notes on GEN.232, TEST.111, TEST.138) (#1652)
- 2026-10-15: TODO: record Quality's confirmations on DB.33, TEST.111 and GEN.232 (#1662)
- 2026-10-15: TODO: design index lists multi-star-systems.md; findings on TEST.111, DB.33, TEST.163, TEST.148; TEST.167 filed (#1696)
- 2026-10-15: TODO: retire TEST.158 (#1697) and DOC.21 (#1698); Quality findings on GEN.231, TEST.156, DOC.24, DB.33 (#1699)
- 2026-10-15: TODO: DB.33 confirmed again on clean main by Core 3 (#1730)
- 2026-10-15: TODO: GEN.158 model slice noted, DB.33 owner is Core 2, TEST.163 sighting, file DOC.53 (#1740)
- 2026-10-21: Skip and retire PERF.131 (Boss, 11:27 AM EDT); DB.33 owner Core 1; restore a catalogue row label that test_system_type_rates reads (#2315)
- 2026-10-22: DB.33: tests that keep a foreign key's default delete rule from reading as drift
- 2026-10-22: Retire DB.33 (#2436); file TEST.267; gate environment and grant notes corrected (#2437)

## DB.34

2026-10-15 to 2026-10-15; 1 commits; not open.

- 2026-10-15: TODO: retire GEN.223, GEN.230, GEN.233, DB.34, DB.35, GEN.260 (#1757), file TEST.170 (#1759)

## DB.35

2026-10-15 to 2026-10-15; 1 commits; not open.

- 2026-10-15: TODO: retire GEN.223, GEN.230, GEN.233, DB.34, DB.35, GEN.260 (#1757), file TEST.170 (#1759)

## DB.36

2026-10-22 to 2026-10-22; 1 commits; not open.

- 2026-10-22: DB.36: change plans, plan operations and the edit journal (migration 0122)

## DB.37

2026-10-22 to 2026-10-22; 1 commits; not open.

- 2026-10-22: DB.37: object locks table (migration 0123)

## DB.38

2026-10-22 to 2026-10-22; 1 commits; not open.

- 2026-10-22: DB.38: region density overrides table (migration 0124)

## DB.39

2026-10-22 to 2026-10-22; 1 commits; not open.

- 2026-10-22: DB.39: request comments and notifications tables (migration 0126)

## DB.40

2026-11-05 to 2026-11-05; 2 commits; not open.

- 2026-11-05: DB.40: the parity file of the DB repair is a backup part (settings.backup_parity) (#5173)
- 2026-11-05: docs: 8.2 headline, logging foundation OPS.104 and OPS.105, retire DB.40 and DB.41 (#5190)

## DB.41

2026-11-05 to 2026-11-05; 2 commits; not open.

- 2026-11-05: DB.41: put back a backed-up parity file and check it against the database (#5188)
- 2026-11-05: docs: 8.2 headline, logging foundation OPS.104 and OPS.105, retire DB.40 and DB.41 (#5190)

## DB.43

2026-10-22 to 2026-10-22; 2 commits; not open.

- 2026-10-22: Add the Reed-Solomon codec in the standard library (DB.43)
- 2026-10-22: docs: DB.43 retired (#2423) (#2489)

## DB.44

2026-10-22 to 2026-10-22; 2 commits; not open.

- 2026-10-22: DB.44: the sector storage checksum and the parity dirty flag (migration 0115) (#2420)
- 2026-10-22: Retire DB.44 (#2420); file DB.83 (galaxy_block_stats maintainer) and DB.84 (sectors FULLTEXT rebuilds); Alembic head 0115 (#2421)

## DB.45

2026-10-28 to 2026-11-04; 3 commits; not open.

- 2026-10-28: DB.45: storage checksum and the dirty flag on every writer (#3414)
- 2026-10-28: Docs: retire DB.45 (#3414), MAP.320 (#3415); file TEST.418, UX.233 (#3416)
- 2026-11-04: A rogue's deletion and an uncharted landing sector mark the storage checksum (DB.45)

## DB.46

2026-10-28 to 2026-10-28; 2 commits; not open.

- 2026-10-28: Add the sector parity file and the refresh-parity command (DB.46) (#3451)
- 2026-10-28: Retire DB.46, file TEST.423 (#3452)

## DB.47

2026-10-28 to 2026-10-28; 2 commits; not open.

- 2026-10-28: Add planetgen repair-db, the sector rebuild from the parity file (DB.47) (#3456)
- 2026-10-28: Retire DB.47 (#3457)

## DB.48

2026-10-28 to 2026-10-28; 2 commits; not open.

- 2026-10-28: Add the table rebuild steps of repair-db, salvage and reload (DB.48) (#3460)
- 2026-10-28: Retire DB.48 and DB.9, file GEN.847 and DB.109 (#3461)

## DB.49

2026-10-22 to 2026-10-22; 1 commits; not open.

- 2026-10-22: DB.49: leaf digest column and sector tree nodes (migration 0120)

## DB.50

2026-10-18 to 2026-10-18; 1 commits; not open.

- 2026-10-18: The description is one Document built from stored data (DB.50, DB.51, DB.52, GEN.218-221, UX.108, UX.119, UX.120, TEST.155) (#2206)

## DB.51

2026-10-18 to 2026-10-18; 1 commits; not open.

- 2026-10-18: The description is one Document built from stored data (DB.50, DB.51, DB.52, GEN.218-221, UX.108, UX.119, UX.120, TEST.155) (#2206)

## DB.52

2026-10-18 to 2026-10-18; 1 commits; not open.

- 2026-10-18: The description is one Document built from stored data (DB.50, DB.51, DB.52, GEN.218-221, UX.108, UX.119, UX.120, TEST.155) (#2206)

## DB.53

2026-10-17 to 2026-10-17; 4 commits; not open.

- 2026-10-17: Point columns: DB.55, GEN.336, GEN.337 and the orbit node table of DB.53 (migrations 0102 to 0104) (#1952)
- 2026-10-17: TODO: retire DB.55, GEN.336, GEN.337 (#1952); DB.53 first slice and TEST.128 sighting recorded (#1954)
- 2026-10-17: Docs: DB.53, GEN.289 and GEN.294 go to Blocker 1 (Boss, 3:40 PM EDT) (#2082)
- 2026-10-17: Docs: retire DB.53 and GEN.294 (done on main by #1952 and its subitems); clean stale prerequisites (#2083)

## DB.54

2026-10-16 to 2026-10-16; 4 commits; not open.

- 2026-10-16: Record migration revisions 0102 (DB.59) and 0103 (DB.54) (#1854)
- 2026-10-16: Retire DB.59 (already built, found by Core 2); DB.54 takes revision 0102 (#1880)
- 2026-10-16: DB.54: the scatter limits by place as three galaxy settings (#1887)
- 2026-10-16: Retire DB.54 (#1887); file TEST.197 and TEST.198; migration order 0101 to 0104 (#1888)

## DB.55

2026-10-17 to 2026-10-17; 2 commits; not open.

- 2026-10-17: Point columns: DB.55, GEN.336, GEN.337 and the orbit node table of DB.53 (migrations 0102 to 0104) (#1952)
- 2026-10-17: TODO: retire DB.55, GEN.336, GEN.337 (#1952); DB.53 first slice and TEST.128 sighting recorded (#1954)

## DB.56

2026-10-22 to 2026-10-22; 4 commits; not open.

- 2026-10-22: docs: TEST.268 retired (#2502), migration 0118 held by DB.56 (#2503)
- 2026-10-22: DB.56: the orbital_events and orbital_event_bodies tables (migration 0118) (#2505)
- 2026-10-22: docs: DB.56 (#2505) and the four closing audits (#2504) retired, GEN.622 filed, PERF.146 owned (#2506)
- 2026-10-22: DB.56: the orbital event columns declare their units in the column registry (#2510)

## DB.58

2026-10-22 to 2026-10-22; 2 commits; not open.

- 2026-10-22: DB.58: control schema v15, account role, email, status, single owner, account_tokens and invite_uses
- 2026-10-22: docs: DB.58 and DB.62 retired (#2536, #2537), TEST.280 filed, DB.85 gets 0130 and epoch 51 (#2538)

## DB.59

2026-10-16 to 2026-10-16; 2 commits; not open.

- 2026-10-16: Record migration revisions 0102 (DB.59) and 0103 (DB.54) (#1854)
- 2026-10-16: Retire DB.59 (already built, found by Core 2); DB.54 takes revision 0102 (#1880)

## DB.60

2026-10-22 to 2026-10-22; 1 commits; not open.

- 2026-10-22: DB.60: density multiplier used and override id on sectors (migration 0125)

## DB.61

2026-10-22 to 2026-10-22; 1 commits; not open.

- 2026-10-22: DB.61: unnamed clusters and the star-to-cluster link (migration 0127)

## DB.62

2026-10-21 to 2026-10-22; 3 commits; not open.

- 2026-10-21: Retire DB.31 (umbrella, subitems done in #2206); fix the control version wording on DB.62 (#2393)
- 2026-10-22: DB.62: control schema v16, the upload bookkeeping tables
- 2026-10-22: docs: DB.58 and DB.62 retired (#2536, #2537), TEST.280 filed, DB.85 gets 0130 and epoch 51 (#2538)

## DB.63

2026-10-22 to 2026-10-22; 4 commits; not open.

- 2026-10-22: DB.63: ships and stations table (migration 0129)
- 2026-10-22: docs: DB galaxy-schema items DB.30 to DB.63 retired (#2520 to #2530, Alembic head 0129), GEN.379 retired (#2531), DB.85 proposed for the DB lane (#2532)
- 2026-10-22: DB.63: the precision test and design table name vessels (floor 9,000 km like the other galaxy-frame tables)
- 2026-10-22: DB.63: vessels is a kind table in the uid group and the ids order

## DB.64

2026-10-22 to 2026-10-22; 1 commits; not open.

- 2026-10-22: DB.64: constellations table (migration 0128)

## DB.65

2026-10-16 to 2026-10-17; 4 commits; not open.

- 2026-10-16: TODO: file the GEN.372 audit's items (MAP.239/240, GEN.374-376, DB.65, TEST.200/201); GEN.109 precision note (#1899)
- 2026-10-16: TODO: retire DB.18 (#1901); browser-hang sighting on TEST.193/195; DB.65 revision 0105 (#1902)
- 2026-10-17: DB.65: the cell key holds the last cube at the edge of its box (#2025)
- 2026-10-17: TODO: DB.65 retired (#2025), migration head 0106 (#2026)

## DB.66

2026-10-17 to 2026-10-22; 7 commits; not open.

- 2026-10-17: Docs: GEN.533 (the bar draw) moves into 8.1 without DB.66; Core 1's scope for GEN.532 (#2054)
- 2026-10-21: Docs: galaxy types table approved by Boss (DB.66 and GEN.379 unblocked) (#2384)
- 2026-10-21: Docs: owner lines for DB.66 and GEN.379 (#2386)
- 2026-10-21: Docs: file DB.82 (VALUES() test failure), DB.66 and GEN.379 owned by Blocker 2 (#2387)
- 2026-10-22: DB.66: the galaxy family and component columns of galaxy_shape (migration 0117) (#2439)
- 2026-10-22: Retire DB.66 (#2439); Alembic head 0117; GEN.379 reads the NULL columns (#2440)
- 2026-10-22: GEN.379: the density as a sum of named components, component_masses in solar masses, the DB.66 columns read and written (#2531)

## DB.67

2026-10-28 to 2026-10-28; 3 commits; not open.

- 2026-10-28: Docs: owners for GEN.415, GEN.423, DB.67 (Core 4) and GEN.479 (Aux 9) (#3435)
- 2026-10-28: Hold the stored density, scale height and greenhouse to the stored mixture (GEN.423, DB.67) (#3474)
- 2026-10-28: Retire GEN.423 and DB.67 (#3475)

## DB.68

2026-10-24 to 2026-10-24; 2 commits; not open.

- 2026-10-24: GEN.414: haze, dust and cloud depths (DB.68) (#2856)
- 2026-10-24: Docs: GEN.414, DB.68 merged (#2856); epoch 66, migration 0138; file TEST.318; Q5 and OPS.67 decided; 832 of 1243 (#2858)

## DB.69

2026-10-24 to 2026-10-24; 2 commits; not open.

- 2026-10-24: GEN.429, GEN.428, DB.69, GEN.459: bulk structure kernels, bulk composition draw, the composition columns (schema 144, epoch 70); the radiogenic heat reads the body's own rock share and Al2O3
- 2026-10-24: Retire GEN.428, GEN.429, DB.69 (#2751); next free 0145 and epoch 71 (#2909)

## DB.70

2026-10-28 to 2026-10-30; 6 commits; not open.

- 2026-10-28: DB.70: store a moon's origin and its giant parent's circumplanetary disk
- 2026-10-29: DB.70: renumber to migration 0176 on main's v175
- 2026-10-29: Train 1: merge claude/project-thread-iku5yl (GEN.462, DB.70, TEST.427: sky colour over the stored ground mix, moon origin and giant disk stored, three red tests fixed)
- 2026-10-29: Train 1: DB.70 takes migration 0177 (GEN.908 holds 0176); v176 fixture and README line
- 2026-10-30: Train 1: merge main; migrations renumbered to 0178 (GEN.908) and 0179 (DB.70), registry sections 6.47 to 6.49, schema fixtures v177 and v178
- 2026-10-30: Train 1: merge main (migrations 0179 GEN.908 and 0180 DB.70; tests follow a renumbering by stem)

## DB.71

2026-10-17 to 2026-10-24; 3 commits; not open.

- 2026-10-17: TODO: atmosphere composition research merged into docs/design, GEN.431 retired, GEN.452-456 and DB.71 filed (#1970)
- 2026-10-24: DB.71, GEN.424: the air is stored gas by gas in planet_gases and moon_gases (schema 155, epoch 95) (#3110)
- 2026-10-24: docs: retire GEN.726, DB.71, GEN.424; captain-first pages chain UX.188 to UX.194 (#3112)

## DB.72

2026-10-17 to 2026-10-17; 4 commits; not open.

- 2026-10-17: Docs: file DB.72, the band_index SMALLINT overflow that stopped Boss's whole-sphere scatter; PERF.36 sighting (#2055)
- 2026-10-17: Docs: DB.72 is in progress with Bug Hunt (#2056)
- 2026-10-17: DB.72: band_index is INT in the three scatter tables (#2060)
- 2026-10-17: Docs: DB.72 merged (#2060); AD1 methane filled from free data (partly verified) (#2065)

## DB.73

2026-10-17 to 2026-10-24; 6 commits; not open.

- 2026-10-17: TODO: file the unified orbited-body chain under GEN.550 (GEN.551 to GEN.560, DB.73, UX.138, TEST.239, PERF.136), notes on GEN.315, 428, 142, 291, 346, 332 (#2170)
- 2026-10-17: TODO: file GEN.564 (refill job), record Core 1's refill and migration answers on GEN.550 and DB.73 (#2176)
- 2026-10-18: docs: GEN.565b re-scoped to motion in the core (Core 4), noted on GEN.565, GEN.584, DB.73 and DB.74 (#2290)
- 2026-10-22: Retire TEST.241 (#2454); assign GEN.569 to Core 4; record Core 1 order GEN.598, DB.76, DB.73, DB.74 (#2455)
- 2026-10-24: DB.73: one parent link for every kind, and moons under a free planet (#2863)
- 2026-10-24: Docs: GEN.693, GEN.694, DB.73 merged, retired; file GEN.696 to GEN.698, TEST.322; migration 0140 taken; 836 of 1264 (#2864)

## DB.74

2026-10-17 to 2026-10-22; 12 commits; open.

- 2026-10-17: TODO: mass object class decided (GEN.565, DB.74, TEST.240, GEN.566), stable object IDs (DB.75) (#2177)
- 2026-10-18: docs: retire GEN.404 (#2270), epoch 34 stamped, record Core 4's GEN.553 and GEN.565 plan and the DB.74 notes (#2271)
- 2026-10-18: docs: GEN.565b re-scoped to motion in the core (Core 4), noted on GEN.565, GEN.584, DB.73 and DB.74 (#2290)
- 2026-10-22: Retire TEST.241 (#2454); assign GEN.569 to Core 4; record Core 1 order GEN.598, DB.76, DB.73, DB.74 (#2455)
- 2026-10-22: docs: GEN.598 retired (#2478), epoch 49 taken, DB.74 plan (#2480)
- 2026-10-22: DB.74: one registry of kind tables replaces the per-site lists
- 2026-10-22: DB.74: one insert driven by a list of fields for every mass-object kind
- 2026-10-22: DB.74: one load for every kind, from the same field lists that write the rows
- 2026-10-22: DB.74: design note for the third slice
- 2026-10-22: DB.74: tests pass the star's node and offset context
- 2026-10-22: DB.74: readers keep the order of the lists the registry replaced
- 2026-10-22: docs: DB.77 retired (#2534), DB.74 slices 1 to 3 recorded (#2535)

## DB.75

2026-10-17 to 2026-10-22; 4 commits; open.

- 2026-10-17: TODO: mass object class decided (GEN.565, DB.74, TEST.240, GEN.566), stable object IDs (DB.75) (#2177)
- 2026-10-22: DB.75: a regenerated system's bodies keep the ID of the body that held their slot
- 2026-10-23: DB.75: the API regenerate test expects the slot's ID to be kept
- 2026-10-22: Docs: DB.75 first part merged (#2568), item stays open (#2569)

## DB.76

2026-10-18 to 2026-10-22; 5 commits; not open.

- 2026-10-18: TODO: file the unit plan (GEN.568 to GEN.575, UX.139, UX.140, DB.76, DB.77, DOC.70) and the design note from Boss's six cards (#2181)
- 2026-10-22: Retire TEST.241 (#2454); assign GEN.569 to Core 4; record Core 1 order GEN.598, DB.76, DB.73, DB.74 (#2455)
- 2026-10-22: DB.76: every numeric column declares its stored unit in one registry
- 2026-10-22: DB.76: the column registry names units of the unit registry
- 2026-10-22: docs: DB.76 retired (#2487), DB.86 filed, TEST.268 recurrence (#2488)

## DB.77

2026-10-18 to 2026-10-22; 5 commits; not open.

- 2026-10-18: TODO: file the unit plan (GEN.568 to GEN.575, UX.139, UX.140, DB.76, DB.77, DOC.70) and the design note from Boss's six cards (#2181)
- 2026-10-22: DB.77: the precision table, each length column's role, frame and reach
- 2026-10-22: DB.77: no TODO code in a string the code shows
- 2026-10-22: DB.77: classify the orbital-event columns in the precision table
- 2026-10-22: docs: DB.77 retired (#2534), DB.74 slices 1 to 3 recorded (#2535)

## DB.78

2026-10-18 to 2026-10-18; 1 commits; open.

- 2026-10-18: TODO: file the System plan (TEST.241, GEN.576 to GEN.580, DB.78, MAP.253, TEST.242, DOC.71), edits to GEN.552/GEN.315/GEN.264, design note from Boss's nine cards (#2185)

## DB.79

2026-10-21 to 2026-10-21; 3 commits; not open.

- 2026-10-21: TODO: retire PERF.139 and PERF.133 (#2316), record cache.max_mb notes, DB.79 and PERF.140 owner (#2317)
- 2026-10-21: DB.79: the object index table (migration 0113, schema v113) (#2361)
- 2026-10-21: TODO: retire DB.79 (#2361, Alembic 0113), file TEST.263 (test_every_check_rejects_a_bad_row fails on main), TEST.262 sighting (#2362)

## DB.80

2026-10-21 to 2026-10-22; 4 commits; not open.

- 2026-10-21: TODO: file GEN.597 (flux_variation above 1) and DB.80 (VALUES ROW on MariaDB) (#2334)
- 2026-10-21: Raise the documented database minimum to MySQL 8.0.19 and drop MariaDB from the install docs (DB.80) (#2336)
- 2026-10-22: DB.80: drop the MariaDB leftovers (MySQL 8.0.19 is the minimum) (#2406)
- 2026-10-22: Docs: UX.141 second part, retire DB.80, GEN.621 in the dials goal, load flake sightings (#2407)

## DB.81

2026-10-21 to 2026-10-22; 4 commits; not open.

- 2026-10-21: TODO: the galaxy shape dials chain (GEN.601 to GEN.611, DB.81, ADM.124, goal #2355): one lever table for the CLI and the Generate form, Shape selector unfilled (Boss, 6:41 to 6:49 PM EDT) (#2356)
- 2026-10-21: TODO: the Galaxy dials lane row and the owner lines for GEN.601 to GEN.611 and DB.81 (coordinator, 6:55 PM EDT) (#2357)
- 2026-10-22: GEN.603 to GEN.606, GEN.609 to GEN.611, DB.81: Structure, Scale, Halo and Edge dials, stored pins, arm contrast (#2403)
- 2026-10-22: Retire the galaxy dials delivered by #2403 (GEN.602 to 606, 609 to 611, DB.81); file GEN.621 (Auto default); 0114 and epoch 44 taken (#2405)

## DB.82

2026-10-21 to 2026-10-21; 3 commits; not open.

- 2026-10-21: Docs: file DB.82 (VALUES() test failure), DB.66 and GEN.379 owned by Blocker 2 (#2387)
- 2026-10-21: DB.82: the object index upsert takes its rows from a derived table, not the deprecated VALUES() (#2390)
- 2026-10-21: Record the GEN.537 plan configuration decisions and the region removal item GEN.620; retire DB.82 (#2391)

## DB.83

2026-10-22 to 2026-10-28; 6 commits; not open.

- 2026-10-22: Retire DB.44 (#2420); file DB.83 (galaxy_block_stats maintainer) and DB.84 (sectors FULLTEXT rebuilds); Alembic head 0115 (#2421)
- 2026-10-24: DB.83: a sector generated on a visit or again joins the block totals (#2960)
- 2026-10-24: Docs: DB.83 partial note, TEST.342, DB.100 to the pool (#2961)
- 2026-10-28: Retire GEN.844, release DB.83, note the red golden seed test (#3444)
- 2026-10-28: Keep the block totals true after an edit (DB.83) (#3476)
- 2026-10-28: Retire DB.83 (#3477)

## DB.84

2026-10-22 to 2026-10-24; 3 commits; not open.

- 2026-10-22: Retire DB.44 (#2420); file DB.83 (galaxy_block_stats maintainer) and DB.84 (sectors FULLTEXT rebuilds); Alembic head 0115 (#2421)
- 2026-10-24: DB.84: write down the rule for changing the sectors table (#2984)
- 2026-10-24: Retire GEN.649, GEN.653, DB.91 (#2930) and DB.84 (#2984); renumber epochs and migrations (#2985)

## DB.85

2026-10-22 to 2026-10-22; 6 commits; not open.

- 2026-10-22: Docs: retire MAP.139 (#2468) and GEN.617 (#2469); narrow MAP.254; file DB.85 (#2470)
- 2026-10-22: docs: DB galaxy-schema items DB.30 to DB.63 retired (#2520 to #2530, Alembic head 0129), GEN.379 retired (#2531), DB.85 proposed for the DB lane (#2532)
- 2026-10-22: docs: DB.58 and DB.62 retired (#2536, #2537), TEST.280 filed, DB.85 gets 0130 and epoch 51 (#2538)
- 2026-10-22: DB.85: globular clusters and superbubbles stored as containers (schema v130)
- 2026-10-23: DB.85: route snapshot gains the four container routes
- 2026-10-22: Docs: DB.85 (#2582) and GEN.498 (#2583) retired; full progress documentation update (versions, tree, plans, HTML) (#2584)

## DB.86

2026-10-22 to 2026-10-23; 2 commits; not open.

- 2026-10-22: docs: DB.76 retired (#2487), DB.86 filed, TEST.268 recurrence (#2488)
- 2026-10-23: Unit registry sweep and debris density: DB.86, GEN.587 group 1, GEN.659, UX.128 (#2801)

## DB.87

2026-10-22 to 2026-10-22; 1 commits; not open.

- 2026-10-22: docs: GEN.348, TEST.279, GEN.585, GEN.586 retired, DB.87 retired, DB.88 filed (#2543)

## DB.88

2026-10-22 to 2026-10-22; 2 commits; not open.

- 2026-10-22: docs: GEN.348, TEST.279, GEN.585, GEN.586 retired, DB.87 retired, DB.88 filed (#2543)
- 2026-10-22: docs: DB.88 and TEST.282 retired (#2545, #2547), PERF.146 findings filed as PERF.156 to PERF.162 and GEN.624 (#2548)

## DB.89

2026-10-23 to 2026-10-24; 4 commits; not open.

- 2026-10-23: Docs: giant interiors (GEN.641, GEN.642, DB.89), habitable = PHI-6 Blue or Green (GEN.643, GEN.644), lane orders (#2710)
- 2026-10-24: DB.89: a giant's interior is saved with the world (schema 154) (#3089)
- 2026-10-24: docs: retire DB.89, PERF.173; physics contradiction items; file TEST.357 (#3091)
- 2026-10-24: TEST.358: the property sheet shows the interior columns of DB.89 (#3123)

## DB.90

2026-10-23 to 2026-11-06; 3 commits; not open.

- 2026-10-23: Docs: conditions-based planet classes (GEN.645 to GEN.648, DB.90, UX.150, TEST.301, DOC.76), Boss's GEN.644 answer (#2711)
- 2026-10-31: File best-fit class items GEN.1205 to GEN.1208; notes on GEN.647, GEN.689, GEN.1144, GEN.1147, GEN.1160, DB.90, TEST.443, DOC.76 (#4464)
- 2026-11-06: Retire ten rows already built: GEN.690, 687, 685, 647, 689, DB.90, 863, 862, 1211, 1296 (#5373)

## DB.91

2026-10-23 to 2026-10-24; 2 commits; not open.

- 2026-10-23: Docs: the four life scores (GEN.649, DB.91, UX.151, UX.152, DOC.77), human operability removed (#2713)
- 2026-10-24: Retire GEN.649, GEN.653, DB.91 (#2930) and DB.84 (#2984); renumber epochs and migrations (#2985)

## DB.92

2026-10-23 to 2026-10-24; 3 commits; not open.

- 2026-10-23: Docs: GEN.648 ladder, GEN.639 decided, aquifers GEN.654 and DB.92 (Boss) (#2722)
- 2026-10-24: GEN.654, DB.92: aquifers, liquid water in the crust of a world with no liquid at its surface; stored depth and thickness; the Water rung, the solvent, the dose and the energy read them
- 2026-10-24: Docs: unused stored values chains (POP.40-45, GEN.804-811, UX.224-228, PERF.195-197, OPS.72-74), air fixed point findings GEN.812-818, bugs ADM.134 and TEST.387; retire GEN.654, DB.92, GEN.751 (#3223)

## DB.93

2026-10-23 to 2026-10-23; 3 commits; not open.

- 2026-10-23: Docs: survival equipment lane and items GEN.655 to GEN.658, DB.93, UX.153, UX.154, POP.15 (#2736)
- 2026-10-23: GEN.655, DB.93: the survival equipment line replaces equipment_tier (schema v136) (#2766)
- 2026-10-23: Docs: GEN.655, DB.93, POP.15 merged (#2766), retired; migration 0137 next; 806 of 1206 (#2819)

## DB.94

2026-10-23 to 2026-10-24; 3 commits; not open.

- 2026-10-23: Docs: retire GEN.622 (#2739), file UX.155, UX.156, DB.94, DB.95, GEN.659, correct the interstellar debris density (#2741)
- 2026-10-24: DB.94: moons of free planets in the object index, lookup, Search and their own page (+ TEST.322, TEST.325) (#2881)
- 2026-10-24: Retire DB.94, GEN.671, UX.163, TEST.322, TEST.325; drop DB.95; file TEST.327, TEST.328 (#2882)

## DB.95

2026-10-23 to 2026-10-24; 2 commits; not open.

- 2026-10-23: Docs: retire GEN.622 (#2739), file UX.155, UX.156, DB.94, DB.95, GEN.659, correct the interstellar debris density (#2741)
- 2026-10-24: Retire DB.94, GEN.671, UX.163, TEST.322, TEST.325; drop DB.95; file TEST.327, TEST.328 (#2882)

## DB.96

2026-10-23 to 2026-10-28; 4 commits; not open.

- 2026-10-23: Docs: GEN.450 merged (#2776, epoch 64), file GEN.669 residuals and the body-rendering items MAP.270 to MAP.278, DB.96, GEN.670, TEST.308, UX.158 (#2777)
- 2026-10-24: docs: file the four wiring items (GEN.759-762), trim GEN.737, assign MAP.287 and DB.96 (#3100)
- 2026-10-28: DB.96: store the atmosphere colour of planets and moons (schema v161, epoch 110) (#3393)
- 2026-10-28: Docs: retire DB.96 (#3393); epoch 112 and migration 0162 taken; file TEST.413 (#3394)

## DB.97

2026-10-23 to 2026-10-24; 5 commits; not open.

- 2026-10-23: Docs: file TEST.309 (merged #2771), UX.159, DB.97, DB.98 for the regression-test cleanup lane (#2780)
- 2026-10-24: Refuse a database with no recorded schema version (DB.97)
- 2026-10-24: Merge main into the DB.97 branch
- 2026-10-24: Merge main into the DB.97 branch
- 2026-10-24: Retire DB.97 (#2905); record the no-marker migration contract (#2906)

## DB.98

2026-10-23 to 2026-10-24; 5 commits; not open.

- 2026-10-23: Docs: file TEST.309 (merged #2771), UX.159, DB.97, DB.98 for the regression-test cleanup lane (#2780)
- 2026-10-24: Re-baseline the migration history at v130 (DB.98)
- 2026-10-24: Merge main into the DB.98 branch
- 2026-10-24: Merge main into the DB.98 branch
- 2026-10-24: Retire DB.98 (#2931); regenerate the architecture map (#2932)

## DB.99

2026-10-23 to 2026-11-04; 2 commits; not open.

- 2026-10-23: Docs: body-rendering items placed by version (Boss's plan, 10:26 PM EDT), build lane row, new items MAP.279 to MAP.285 and DB.99, TEST.308 dropped (#2782)
- 2026-11-04: Retire GEN.1147 (merged with train 8, #4655) and drop DB.99 (ground type and fractions already stored) (#4945)

## DB.100

2026-10-24 to 2026-10-24; 4 commits; not open.

- 2026-10-24: Docs: file OPS.69 (strip wiki, keep Markdown export) and DB.100 (store life and evolution data) (#2959)
- 2026-10-24: Docs: DB.83 partial note, TEST.342, DB.100 to the pool (#2961)
- 2026-10-24: DB.100: store the milestone schedule of a world's pace; the life text reads the row alone; an edit keeps the life a world has
- 2026-10-24: Docs: retire UX.125 and DB.100, file TEST.377 (#3184)

## DB.101

2026-10-24 to 2026-10-24; 1 commits; open.

- 2026-10-24: docs: wiring audit, file DB.101 and UX.187 (#3109)

## DB.102

2026-10-24 to 2026-10-25; 4 commits; not open.

- 2026-10-24: docs: file DB.102 (version key carries API, schema and epoch) (#3138)
- 2026-10-24: docs: file ADM.129 (admin buttons for habitability re-score and map luminosity fill); DB.102 order fix (#3139)
- 2026-10-25: The version key carries the API version, the schema and the epoch (DB.102) (#3323)
- 2026-10-25: Docs: file MAP.326 and MAP.327, MAP.325 darker background, pragmatic-tweak rule, retire DB.102 (#3324)

## DB.105

2026-10-24 to 2026-10-24; 1 commits; open.

- 2026-10-24: Docs: file the world population and trade routes tree (POP.29-POP.36, DB.105, MAP.309, MAP.310, NAV.63, TEST.382) (#3204)

## DB.106

2026-10-24 to 2026-10-24; 1 commits; open.

- 2026-10-24: Docs: retire TEST.375, drop TEST.374; file the resource locations tree (GEN.798-799, DB.106, MAP.312-313, UX.223, TEST.386, NAV.64, ADM.133, DOC.81, API.34) (#3219)

## DB.107

2026-10-28 to 2026-10-28; 3 commits; not open.

- 2026-10-28: Docs: retire GEN.825 and GEN.827 (#3375); migration 0161 goes to DB.107 (#3376)
- 2026-10-28: Stored regions: region_stats table, hooks, migration 0161 (DB.107) (#3383)
- 2026-10-28: Docs: retire DB.107, GEN.797, TEST.401, TEST.394, TEST.383; file TEST.410, TEST.411; GEN.328 queued (#3386)

## DB.108

2026-10-25 to 2026-10-28; 3 commits; not open.

- 2026-10-25: Docs: POI icon items MAP.320 to MAP.322 and DB.108, regions plan priority; retire TEST.370, TEST.390 (#3295)
- 2026-10-28: Stored sector, block and region mass (DB.108) (#3465)
- 2026-10-28: Retire DB.108 (#3466)

## DB.109

2026-10-28 to 2026-10-28; 1 commits; open.

- 2026-10-28: Retire DB.48 and DB.9, file GEN.847 and DB.109 (#3461)

## DB.112

2026-11-05 to 2026-11-05; 1 commits; open.

- 2026-11-05: docs: file UX.272 and DB.112, replace the quadrants with the charted region names (#5215)

## DB.113

2026-11-06 to 2026-11-06; 4 commits; not open.

- 2026-11-06: Docs: file DB.113, the massive-star scatter deadlock crash (stage 7, unit 11) (#5390)
- 2026-11-06: Docs: DB.113/114 dynamic retries; fudge table selection rule and cap tiers (#5395)
- 2026-11-06: DB.113: scatter batch writes retry after a MySQL deadlock (#5392)
- 2026-11-06: Retire DB.113 and PERF.245 after #5392 and #5399 merged (#5408)

## DB.114

2026-11-06 to 2026-11-07; 4 commits; not open.

- 2026-11-06: Docs: file DB.114, the backfill and sky-charting ID writes need a deadlock retry (#5393)
- 2026-11-06: Docs: DB.114 plan uses run_with_retries with the whole transaction as the attempt (#5396)
- 2026-11-07: DB.114: backfills, special phenomena and sky charting retry after a MySQL deadlock (#5504)
- 2026-11-07: docs: retire train 13 items and DB.114; file TEST.670 to TEST.673 (reds on main) (#5508)

## DB.115

2026-11-07 to 2026-11-07; 1 commits; open.

- 2026-11-07: File DB.115: find why MySQL dropped the connection in a long generation job (#5557)

## DB.116

2026-11-08 to 2026-11-08; 1 commits; open.

- 2026-11-08: Docs: file the twenty fictional materials parts for 8.2 (GEN.1408 to GEN.1419, DB.116/117, ADM.143/144, API.42, UX.293, TEST.699, DOC.110) (#5710)

## DB.118

2026-11-08 to 2026-11-08; 3 commits; not open.

- 2026-11-08: Docs: file Lane 05's red causes as items (UX.294, TEST.701-705, DB.118-120, DOC.111-112, GEN.1439, PERF.255) (#5718)
- 2026-11-08: Retire TEST.706, TEST.707, TEST.708, DB.118 (PR 5765) (#5766)
- 2026-11-08: Reds: DB.118, TEST.702, DOC.112 (formulas on the References pages) (#5767)

## DB.119

2026-11-08 to 2026-11-08; 1 commits; not open.

- 2026-11-08: Reds: UX.294, TEST.701, TEST.703, DB.119, DB.120, DOC.111 (#5755)

## DB.120

2026-11-08 to 2026-11-08; 1 commits; not open.

- 2026-11-08: Reds: UX.294, TEST.701, TEST.703, DB.119, DB.120, DOC.111 (#5755)

