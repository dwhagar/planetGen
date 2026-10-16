# Issue groups

Boss (2026-10-15): issues are groups of TODO items, the goals we are working toward. Each row is one GitHub issue and the items it covers. `scripts/build_todo_docs.py` reads this table and links every item to its issues in the tier plan, the TODO tree and the TODO reference; `docs/plan/execution-plan.html` links them too.

Rules:

- One issue per goal. An item belongs to the issue that groups it; an item may appear in more than one row.
- Issues Boss opened keep his text. The TODO lane adds a "TODO sync" comment that lists their items and does not edit his words.
- Findings from any thread that belong to a goal (a failed test, a bug noticed in passing, an oddity) go to the TODO thread, which attaches them to the goal's issue as a comment headed `Finding (date, source, kind)` with the test or object, the failure text, whether it fails on main, how often, and a reproducer. An item with no goal yet gets a new issue.
- Issue labels: `goal` marks a goal issue; one label per TODO ID (for example `GEN.250`) lets the issue list filter to an item's issue. Bugs also carry `bug`.
- An issue closes when every item in its row is done or dropped; the row stays so the reference page keeps its links.

| Issue | State | Goal | Items |
|---|---|---|---|
| [#1750](https://github.com/dwhagar/planetGen/issues/1750) | open | Site cache warm-up after an update, then rebuild only what a change invalidates | PERF.104 |
| [#1749](https://github.com/dwhagar/planetGen/issues/1749) | closed | Hypervelocity stars show no accessible info; button to the phenomenon page or the star system (MAP.211) | MAP.211 |
| [#1748](https://github.com/dwhagar/planetGen/issues/1748) | closed | Phenomena list: the Apply filters button disappears while the list loads (UX.117) | UX.117 |
| [#1785](https://github.com/dwhagar/planetGen/issues/1785) | open | test_gen_resume forced-scatter test failed once under gate load with no scattered system (TEST.171) | TEST.171 |
| [#1811](https://github.com/dwhagar/planetGen/issues/1811) | open | Sector page accessibility test failed once under gate load and passes alone (TEST.174) | TEST.174 |
| [#1808](https://github.com/dwhagar/planetGen/issues/1808) | open | A fill-after-failed-scatter test failed once under gate load with an empty ring (TEST.173) | TEST.173 |
| [#1802](https://github.com/dwhagar/planetGen/issues/1802) | open | Galaxy Map forward-to-current test failed once: a click on a block with 0 sectors charted did nothing (TEST.172) | TEST.172 |
| [#1758](https://github.com/dwhagar/planetGen/issues/1758) | open | test_math_check::test_whole_run_is_fast and the Galaxy Map breadcrumb folding test fail under full-suite load (TEST.170) | TEST.170 |
| [#1726](https://github.com/dwhagar/planetGen/issues/1726) | open | Interface update: every control works on the whole sphere | UX.115, TEST.168 |
| [#1694](https://github.com/dwhagar/planetGen/issues/1694) | open | Test failures under parallel load and on MySQL 8.4 | TEST.116, TEST.167, TEST.126, TEST.127, TEST.128, TEST.134, TEST.137, TEST.140, TEST.141, TEST.144 |
| [#1693](https://github.com/dwhagar/planetGen/issues/1693) | open | View from a planet: sky maps and constellations | VIEW.1, VIEW.2, VIEW.3, VIEW.4, VIEW.6, VIEW.7, VIEW.8, VIEW.9, VIEW.10 |
| [#1691](https://github.com/dwhagar/planetGen/issues/1691) | open | Documentation corrections and the docs index | DOC.4, DOC.21, DOC.24, DOC.52, DOC.53 |
| [#1690](https://github.com/dwhagar/planetGen/issues/1690) | open | Static help section for the web interface | DOC.6, DOC.7, DOC.8, DOC.9, DOC.10, DOC.11, DOC.12, DOC.13, DOC.14, DOC.15, DOC.16, DOC.18 |
| [#1689](https://github.com/dwhagar/planetGen/issues/1689) | open | One logging standard and a themed, componentised interface | UX.102, DOC.19, UX.43, UX.49, UX.82 |
| [#1688](https://github.com/dwhagar/planetGen/issues/1688) | open | Admin tools: editing, configuration and incomplete uploads | ADM.13, ADM.32, ADM.35, ADM.43, ADM.44, ADM.50, ADM.55, ADM.56, TEST.169, UX.3 |
| [#1687](https://github.com/dwhagar/planetGen/issues/1687) | open | Install, update and daily maintenance | OPS.16, OPS.17, OPS.29, OPS.30, OPS.35 |
| [#1686](https://github.com/dwhagar/planetGen/issues/1686) | open | Database repair and migration helpers | DB.9, DB.17, DB.18, DB.31 |
| [#1685](https://github.com/dwhagar/planetGen/issues/1685) | open | Performance guards, caches and benchmarks | PERF.20, PERF.41, PERF.47, PERF.48, PERF.75 |
| [#1684](https://github.com/dwhagar/planetGen/issues/1684) | open | API keys, call log and compatibility data | API.4, API.6, API.15 |
| [#1683](https://github.com/dwhagar/planetGen/issues/1683) | open | Remote generation: generate locally and upload through the API | API.3, API.5, API.7, API.8, API.10, API.11, API.12, API.13, API.14, API.17, API.18, API.19 |
| [#1682](https://github.com/dwhagar/planetGen/issues/1682) | open | Bookmarks and the user-level interface | USR.7, UX.45, UX.46, UX.47, UX.48, UX.92 |
| [#1681](https://github.com/dwhagar/planetGen/issues/1681) | open | User accounts, roles, invitations and password security | USR.1, USR.2, USR.3, USR.4, USR.5, USR.6, USR.8, USR.9, SEC.32 |
| [#1680](https://github.com/dwhagar/planetGen/issues/1680) | open | Courses that steer around gravity wells and keep-out spheres | NAV.6, NAV.25, NAV.26, NAV.27, NAV.28, NAV.51, NAV.54, NAV.55, NAV.56 |
| [#1679](https://github.com/dwhagar/planetGen/issues/1679) | open | Saved courses and waypoints (headline of 8.5) | NAV.4, NAV.5, NAV.17, NAV.18, NAV.19, NAV.21, NAV.22, NAV.23, NAV.36, NAV.39, NAV.49, NAV.58, NAV.59 |
| [#1678](https://github.com/dwhagar/planetGen/issues/1678) | open | Same seed, same data: determinism, fingerprint and generator version | GEN.55, GEN.57, GEN.135, GEN.136, OPS.15, OPS.28, OPS.37, DB.16, TEST.77 |
| [#1677](https://github.com/dwhagar/planetGen/issues/1677) | open | Later catalogue additions and compact-object details (low priority) | GEN.280, GEN.281, GEN.282, GEN.283, GEN.284, GEN.285, GEN.286, GEN.287, UX.116, GEN.244, GEN.245, GEN.246, GEN.254, GEN.255, GEN.256, GEN.257 |
| [#1676](https://github.com/dwhagar/planetGen/issues/1676) | open | One description layer for habitability, units and every object page | GEN.218, GEN.219, GEN.220, GEN.221, UX.32, UX.105, UX.107, UX.108, UX.109, DOC.45, TEST.155 |
| [#1675](https://github.com/dwhagar/planetGen/issues/1675) | open | Planet classes refactored around the habitability index | GEN.27, GEN.28, GEN.29, GEN.33, GEN.90, GEN.91, GEN.92, GEN.145, GEN.146 |
| [#1674](https://github.com/dwhagar/planetGen/issues/1674) | open | Fill planning: order, directives, density and sector weeding | GEN.9, GEN.40, GEN.42, GEN.102, GEN.103, GEN.134, GEN.144, GEN.179, GEN.180, GEN.181, PERF.35, PERF.36 |
| [#1673](https://github.com/dwhagar/planetGen/issues/1673) | open | Neighbour galaxies, satellites, wanderers and intergalactic frame | GEN.157, GEN.238, GEN.239, GEN.248, MAP.145, MAP.207, NAV.57 |
| [#1672](https://github.com/dwhagar/planetGen/issues/1672) | open | Nebulae on the Galaxy Map and in generation | MAP.171, MAP.172, MAP.173, MAP.174, MAP.176, MAP.177, MAP.178, MAP.179, MAP.180, MAP.181, MAP.183, GEN.93, GEN.95, GEN.99, GEN.202, TEST.135, GEN.295 |
| [#1671](https://github.com/dwhagar/planetGen/issues/1671) | open | Fly-through Galaxy Map and scale hand-offs | MAP.59, MAP.75, MAP.119, MAP.120, MAP.121, MAP.122, MAP.139, MAP.140, MAP.141, MAP.142, MAP.143, MAP.146, MAP.147, MAP.152, MAP.156, MAP.160, MAP.195, MAP.209, MAP.210 |
| [#1670](https://github.com/dwhagar/planetGen/issues/1670) | open | Gravity field and the gravity map | GEN.115, GEN.198, GEN.199, MAP.167, MAP.168, MAP.169, MAP.170, API.24, UX.94, PERF.72, TEST.131 |
| [#1669](https://github.com/dwhagar/planetGen/issues/1669) | open | Orbital update: solver, thresholds, collisions and the nightly ceiling | GEN.105, GEN.109, GEN.110, GEN.111, GEN.112, GEN.139, GEN.140, GEN.141, GEN.142, GEN.143, GEN.197, ADM.36, TEST.129, TEST.130, GEN.289, GEN.290, GEN.291, GEN.292, GEN.293, GEN.294 |
| [#1668](https://github.com/dwhagar/planetGen/issues/1668) | open | Nearby-object queries by cell key (stage C) | GEN.276, NAV.45, NAV.60, NAV.61, TEST.165, TEST.166, DOC.51 |
| [#1667](https://github.com/dwhagar/planetGen/issues/1667) | open | Multi-star systems and the scatter member tree (stage M) | GEN.261, GEN.263, GEN.264, GEN.267, GEN.268, GEN.269, GEN.270, GEN.271, GEN.272, GEN.273, GEN.274, GEN.275, GEN.129, GEN.130, TEST.162, TEST.164, DOC.50, UX.114 |
| [#1666](https://github.com/dwhagar/planetGen/issues/1666) | open | Phenomena registry, fill stages and the remaining star types | GEN.250, GEN.251, GEN.252, GEN.253, GEN.247, GEN.249, GEN.258, GEN.259, DOC.49 |
| [#1665](https://github.com/dwhagar/planetGen/issues/1665) | open | Whole-sphere reach: every cell out to 262 kpc has contents, from one rule set | GEN.236, GEN.240, GEN.231, API.26, API.27, UX.113, DOC.47, DOC.48, TEST.156, TEST.158, TEST.159, TEST.160, TEST.161 |
| [#1657](https://github.com/dwhagar/planetGen/issues/1657) | closed | GEN.232 Cap the sector-path integration for fast bodies so one hypervelocity star cannot run for min | GEN.232 |
| [#1656](https://github.com/dwhagar/planetGen/issues/1656) | open | TEST.111 test_ensure_sector_generated_creates_then_reuses_the_same_sector fails in a busy parallel r | TEST.111 |
| [#1655](https://github.com/dwhagar/planetGen/issues/1655) | open | TEST.138 test_open_map_menus_hold_no_overlap fails under load: [web.galaxy-1280] once, [web.sector-6 | TEST.138 |
| [#1654](https://github.com/dwhagar/planetGen/issues/1654) | closed | TEST.163 test_nav_ends_picked_from_bookmarks_on_every_page fails under load: #course-heading count i | TEST.163 |
| [#1653](https://github.com/dwhagar/planetGen/issues/1653) | open | DB.33 MySQL 8.0.46 reflects fk_facilities_type and fk_facility_field_values_field as RESTRICT, so ch | DB.33 |
| [#1565](https://github.com/dwhagar/planetGen/issues/1565) | open | Wide binary star systems must be stable: far orbits and Hill spheres must not overlap | GEN.228, GEN.229 |
| [#1537](https://github.com/dwhagar/planetGen/issues/1537) | open | Orbital update skips the galaxy position move of systems with no placed sector (GEN.223) |  |
| [#1532](https://github.com/dwhagar/planetGen/issues/1532) | open | Overview page: the phenomena count is wrong (UX.110) | UX.110 |
| [#1503](https://github.com/dwhagar/planetGen/issues/1503) | closed | Galaxy Map: regions under the galactic plane are hard to highlight and select (MAP.204) | MAP.204 |
| [#1500](https://github.com/dwhagar/planetGen/issues/1500) | closed | Generation output: nebula kinds and layers per slab; no progress for the phenomena indexes | GEN.222, PERF.101 |
| [#1464](https://github.com/dwhagar/planetGen/issues/1464) | open | List the species that could live in an environment; species environment limits and a species editor | POP.12, POP.13, DB.30, POP.14 |
| [#1454](https://github.com/dwhagar/planetGen/issues/1454) | closed | Habitability text: PHI-6 name, Air at zero pressure, wording of every line | GEN.214, GEN.215, GEN.216, GEN.217 |
| [#1436](https://github.com/dwhagar/planetGen/issues/1436) | open | update.sh and migrations refuse to start while jobs are active unless --yes or confirmed (DB.29) | DB.29 |
| [#1420](https://github.com/dwhagar/planetGen/issues/1420) | closed | Sector view: stars look too big compared to the Galaxy Map (MAP.203) | MAP.203 |
| [#1416](https://github.com/dwhagar/planetGen/issues/1416) | closed | Galaxy Map: a selected sector that is not filled should stay highlighted (MAP.201) | MAP.201 |
| [#1414](https://github.com/dwhagar/planetGen/issues/1414) | closed | Galaxy Map: unclear which block is selected, awkward scroll zoom, and no on-map zoom and move contro | MAP.198, MAP.199, MAP.200 |
| [#1406](https://github.com/dwhagar/planetGen/issues/1406) | closed | Mercury-like bodies never keep their magnetic field (GEN.213) | GEN.213 |
| [#1392](https://github.com/dwhagar/planetGen/issues/1392) | closed | Keep the pymysql login workaround and make it compatible with the upstream fix (OPS.50) | OPS.50 |
| [#1381](https://github.com/dwhagar/planetGen/issues/1381) | open | Investigate resuming a test run where it left off when the build is the same (TEST.148) | TEST.148 |
| [#1362](https://github.com/dwhagar/planetGen/issues/1362) | closed | update.sh fails at the database migrate step: 'NoneType' object has no attribute 'is_auth_switch_req | OPS.49 |
| [#1344](https://github.com/dwhagar/planetGen/issues/1344) | closed | update.sh stops at the database migrate step in the venv: cryptography missing (OPS.47) | OPS.47 |
| [#1337](https://github.com/dwhagar/planetGen/issues/1337) | closed | Installer: check the Apache mod_wsgi venv rules (OPS.46) | OPS.46 |
| [#1328](https://github.com/dwhagar/planetGen/issues/1328) | closed | Save generation timing per object and per object type | PERF.92 |
| [#1296](https://github.com/dwhagar/planetGen/issues/1296) | closed | Admin login returns 500: lock wait timeout on admin_users, and daemon timeouts | ADM.53 |
| [#1295](https://github.com/dwhagar/planetGen/issues/1295) | closed | Install planetGen into its own virtual environment | OPS.45 |
| [#1277](https://github.com/dwhagar/planetGen/issues/1277) | closed | Docs site build (GitHub Pages) fails on every release commit | OPS.44 |
| [#1270](https://github.com/dwhagar/planetGen/issues/1270) | closed | A fresh install gets markdown 3.11.1, which breaks the page renderer | OPS.43 |
| [#1248](https://github.com/dwhagar/planetGen/issues/1248) | open | Sector path integration on numba | PERF.88 |
| [#1247](https://github.com/dwhagar/planetGen/issues/1247) | closed | Drop astropy for a small constants table | GEN.211 |
| [#1246](https://github.com/dwhagar/planetGen/issues/1246) | closed | Require Python 3.12 and install numba | OPS.42 |
| [#1241](https://github.com/dwhagar/planetGen/issues/1241) | open | Faster math and fewer heavy libraries (numba, Python 3.12, constants table) | OPS.42, PERF.87, PERF.88, PERF.89, PERF.90, PERF.91, GEN.211, GEN.212 |
| [#1234](https://github.com/dwhagar/planetGen/issues/1234) | open | Nav Markers | MAP.132 |
| [#1233](https://github.com/dwhagar/planetGen/issues/1233) | open | Default View Galaxy | MAP.191 |
| [#1232](https://github.com/dwhagar/planetGen/issues/1232) | open | Restrict the Field | MAP.189, MAP.190, USR.11 |
| [#1231](https://github.com/dwhagar/planetGen/issues/1231) | open | See most recently changed positions | MAP.192 |
| [#1230](https://github.com/dwhagar/planetGen/issues/1230) | closed | Complete review and rebuild of testing system | TEST.139 |
| [#1229](https://github.com/dwhagar/planetGen/issues/1229) | open | Work queue on Redis: remaining steps and follow-ons | PERF.86, PERF.97, PERF.98, PERF.99, PERF.100, PERF.18, PERF.29, PERF.30, DOC.20 |
| [#1228](https://github.com/dwhagar/planetGen/issues/1228) | open | Star system editor | ADM.51, ADM.52 |
| [#1227](https://github.com/dwhagar/planetGen/issues/1227) | open | Show axis, tilt, and rotations | MAP.193 |
| [#1226](https://github.com/dwhagar/planetGen/issues/1226) | open | New feature:  User objects | USR.10, POP.11 |
| [#1225](https://github.com/dwhagar/planetGen/issues/1225) | closed | Phenomena scatter possible hang before start. | PERF.83 |
| [#1224](https://github.com/dwhagar/planetGen/issues/1224) | closed | Neighbor linking stalled out | PERF.84 |
| [#1223](https://github.com/dwhagar/planetGen/issues/1223) | closed | Sector View Hides Black Holes | MAP.185 |
| [#1222](https://github.com/dwhagar/planetGen/issues/1222) | open | Galactic Core should be a Sector | GEN.210 |
| [#1221](https://github.com/dwhagar/planetGen/issues/1221) | open | Add a 5th and 6th Dimention to PHI-4 (PHI-6) | GEN.207, UX.101 |
| [#1220](https://github.com/dwhagar/planetGen/issues/1220) | closed | Contradictory Planet Information | GEN.206 |
| [#1219](https://github.com/dwhagar/planetGen/issues/1219) | closed | Rogue Planet Gravity | UX.99 |
| [#1218](https://github.com/dwhagar/planetGen/issues/1218) | closed | Mark Rogue Planets Redundant | MAP.184 |
| [#1217](https://github.com/dwhagar/planetGen/issues/1217) | open | Habitability Notes in Planet List | UX.100, GEN.208 |
| [#1216](https://github.com/dwhagar/planetGen/issues/1216) | closed | Fix weird log wording | UX.98 |
| [#1215](https://github.com/dwhagar/planetGen/issues/1215) | closed | Binary Star Sector Display Error | MAP.186 |
| [#1214](https://github.com/dwhagar/planetGen/issues/1214) | closed | Generation Hang after Phenomena / Layer Count Problem | PERF.85 |
| [#1213](https://github.com/dwhagar/planetGen/issues/1213) | closed | Stars Vanish in new 'Window' algorithm | MAP.187 |
| [#1212](https://github.com/dwhagar/planetGen/issues/1212) | closed | Galactic Map Sector Boarders | MAP.202 |
| [#1211](https://github.com/dwhagar/planetGen/issues/1211) | closed | View artifacts in Galactic View | MAP.187 |
| [#1204](https://github.com/dwhagar/planetGen/issues/1204) | open | Nebulas and Globular Clusters | GEN.209, MAP.194, GEN.237, GEN.241, GEN.242, GEN.243, GEN.296, GEN.297, GEN.160, GEN.161, GEN.162, GEN.163, GEN.164 |
