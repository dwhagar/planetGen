# Rebuilding the test suite: faster, modular, and able to run only what a change touches (TEST.139)

Boss (2026-10-14 12:17Z): "rework, optimize, and rebuild our entire test suite from the ground up so that it runs faster and is modular (so we can if we choose only test the things that we changed in the last run)."

Status: measurement and plan, 2026-10-14. Stage 1 onward is not built yet. Evidence tags: [M] measured in this work (numbers below), [S] seen in the repo, [E] estimated from a measurement.

## Summary

- **The full suite is 15,622 tests and took 55 min 12 s on 4 cores** (MariaDB 10.11, Redis, Chromium, `pytest -n 4 --dist worksteal`, 2026-10-14 12:28Z). Adding up every test's own time gives 215 min, so the machine was busy the whole 55 minutes [M].
- **Where the time goes** [M]: the browser tests (571 tests, 5 files) are 42% of the test time; the 25 biggest files are 75%; creating and dropping a database for each database test is about 15%; fuzz and bug-hunt files are only 4%. The 2,206 tests marked `slow` are cheap; the expensive tests are not the ones marked that way.
- **The cheapest large win is not to run browser tests when nothing they cover changed**, and to run them as a separate phase with fewer workers (13 of the 14 failures in the baseline were browser timeouts under load, the open TEST.126). That, together with change-based selection, is what turns a typical lane run from about 55 minutes into a few.
- **Change-based selection works without any new dependency.** A pytest plugin built on the standard library (`sys.monitoring` for Python code, an audit hook for the data files a test reads) records which source files each test executes, at about zero overhead [M]. On a 3,796-test sample, one changed source file selects 9.5% of the tests on average (median 0.8%); with a safety margin of one import hop, 17.6% (median 8.5%). Reading the import statements alone, as most tools do, selects 55% of the test files for a typical change because the packages import each other everywhere [M]. So the recorded map is the right tool.
- **The full suite stays the merge gate** until the selection has been proven (below). Nothing here adds CI cost: no new job, no new matrix leg.
- Plan in five small stages (section 6). Stages 1 and 2 give the most.

## 1. The baseline

Machine: 4 cores, 15 GB, MariaDB 10.11 with the container's defaults, Redis 7, Python 3.13, Chromium 141 through Playwright 1.56. The suite as `docs/testing.md` says to run it, with `--durations=0` and a JUnit file for the per-test times.

| Run | Result |
|---|---|
| Full, `-n 4 --dist worksteal` | 15,548 passed, 14 failed, 60 skipped, 55 min 12 s |
| Same, database and browser unavailable (a first attempt with a wrong database user) | 11,763 passed, 3,791 skipped, 5 min 11 s |

The second row matters: **the part of the suite that needs neither MySQL nor a browser is 5 minutes**. Everything else is the database, the browser or generation.

The 14 failures: 13 are `Page.goto`/`click`/`screenshot` timeouts of 30 s in the WebGL browser tests running four at a time on four cores (TEST.126 describes this), and 1 was my own account lacking `GRANT OPTION` for a test that creates a restricted user (not a bug). Two environment findings for the TODO thread are in section 7.

### Time by file (sum of each test's setup, call and teardown) [M]

| File | Tests | Time | Share |
|---|---|---|---|
| test_web_browser_fixture_maps | 86 | 3,506 s | 27.1% |
| test_neighborhood_edge | 18 | 753 s | 5.8% |
| test_galaxy_gen | 105 | 628 s | 4.9% |
| test_web_browser_maps | 17 | 614 s | 4.8% |
| test_web_browser_layout | 137 | 475 s | 3.7% |
| test_web_browser_controls | 170 | 424 s | 3.3% |
| test_parallel_generation | 15 | 329 s | 2.5% |
| test_prevalence | 27 | 319 s | 2.5% |
| test_web_a11y | 138 | 290 s | 2.2% |
| test_api | 129 | 273 s | 2.1% |
| test_phenomenon_scatter | 76 | 260 s | 2.0% |
| test_fuzz_web_routes | 1,453 | 249 s | 1.9% |
| test_work_queue_failures | 22 | 239 s | 1.9% |

The top 5 files are 46% and the top 25 are 75% of the 12,920 s. The other 282 files are the remaining quarter.

| Group | Time | Share |
|---|---|---|
| Browser files (`test_web_browser*`, `test_web_a11y`) | 5,390 s | 42% |
| Fuzz files (`test_fuzz_*`) | 425 s | 3% |
| Bug-hunt files (`test_bughunt_*`) | 115 s | 1% |
| Everything else | 6,989 s | 54% |

A browser test spends most of its time drawing WebGL in software (SwiftShader), the same CPU the other workers need. The slowest test is 280 s (`test_galaxy_map_wheel_carries_the_zoom_into_a_system_and_back`); a typical fixture-map test is 40 s.

### Setup and teardown [M]

Of the 12,920 s, calls are 9,654 s, setup 1,919 s and teardown 1,322 s. For the database files outside the browser group alone: setup 789 s, teardown 1,137 s, call 4,249 s. The cost is in the fixture that gives every database test its own database: `CREATE DATABASE`, then `schema.sql` (122 tables, applied on the first connection), then `DROP DATABASE`. On an idle machine that is 0.04 s + 0.34 s + 0.28 s = 0.66 s per test; under load the drop doubles. Server settings help little (`innodb_flush_log_at_trx_commit=0` and `innodb_file_per_table=0` together took it to 0.33 s), because the cost is the DDL itself, so the fix is to create fewer databases (stage 3).

### The `db` marker misses most database tests [M]

`conftest.py` adds `db` to a test that names `mysql_config`. 1,465 tests carry it. With the server unreachable, 3,677 tests skipped for lack of one, so **2,210 tests need MySQL but are not marked `db`**: they reach the server through `_mysql_server_available` (`fuzz_db`, `site_db` and a few more in test_fuzz_web_routes, test_web_a11y, test_parallel_generation, test_db_collation and rich_galaxy_support). So `-m "not db"` is not the database-free slice that `docs/testing.md` says it is.

## 2. Principles

1. **One way to pick tests.** A person says what they want (`quick`, `changed`, `area web`, `full`) and a runner builds the pytest command. Plain `pytest` keeps working.
2. **A test is picked by what it depends on, not by where its file sits.** Moving files around does not make the suite modular; knowing what each test touches does.
3. **Over-select, never under-select.** Every shortcut in the selection errs toward running more. The full suite is the check on the shortcuts.
4. **No new dependency that is not light.** Everything below uses pytest, xdist, the standard library and git, all already here.
5. **Small stages.** Other lanes are editing tests all the time; each stage is a small PR and none rewrites a test body that a lane is editing.

## 3. Modular: three layers

### 3.1 Markers (what a test is)

Kept as they are (`db`, `slow`, `browser`, `generation`, `mathcheck`), with these changes:

- `db` also marks a test that uses `_mysql_server_available` (fixes the 2,210).
- New tier markers by cost, set automatically in `conftest.py` from the file name, as now: `browser` (571 tests), `heavy` (the files that take over 60 s of the baseline: a short list kept in `conftest.py`, about 20 files). `slow` stays what it is.
- New area markers set from the file name prefix, with no file moved: `area_web` (`test_web_*`, `test_fuzz_web_*`), `area_api` (`test_api*`), `area_db` (`test_db_*`), `area_generation` (`test_gen*`, `test_galaxy*`, `test_phenom*`, `test_scatter*` and so on), `area_physics` (`test_physics*`, `test_orbit*`, ...), `area_cli`, `area_names`, `area_map` (maps, tiles). The table lives in one place, `src/tests/areas.py`, and a test fails if a test file matches no area, so a new file cannot go unmapped.

### 3.2 Runner (how a person picks)

`scripts/test.py` (no dependency), thin over pytest:

| Command | Runs |
|---|---|
| `scripts/test.py quick` | everything not `db`, `browser` or `heavy`: about 5 minutes today, less after stage 1 |
| `scripts/test.py changed` | the tests selected from the change since the base (3.3) |
| `scripts/test.py area web db` | the tests carrying those area markers |
| `scripts/test.py browser` | the browser phase alone, with 2 workers |
| `scripts/test.py full` | everything: the merge gate. Runs the non-browser tests with all cores, then the browser tests with 2 workers |
| any of them `--record` | also records the map (3.3) |

The `full` command writes the start time and the previous duration to its first line, so a lane can paste it in its status line.

### 3.3 Change-based selection (what a change touches)

**Recording.** A pytest plugin (`src/tests/depmap.py`, on with `--record-deps`) tracks, for every test, the set of repository files it used:

- Python files: `sys.monitoring` (Python 3.12 and later) with a `PY_START` callback that stores the file of each function that starts and then disables itself for that function, and is switched on again at the start of each test. The overhead was not measurable (a 389-test run took 16.7 s with it and 18.6 s without).
- Data files (templates, static JS and CSS, `schema.sql`, fixtures, JSON): an `open` audit hook (`sys.addaudithook`, any Python 3.8+) records every path opened under the repository during a test.
- A test that starts a subprocess or forks (audit events `subprocess.Popen`, `os.fork`, `os.posix_spawn`) is recorded as opaque, since what the child runs is not seen: it depends on everything under `src/planetgen/`.
- The result is one JSON file per checkout (`.pytest-deps.json`, ignored by git), written when the run ends. Under xdist each worker sends its part to the controller, which merges them.
- Recording is only on a run that asks for it, and needs 3.12 or later. Selection itself reads the JSON and runs on any Python.

**Safety margins.** `sys.monitoring` sees functions that run, so a test that only reads a constant from module B, through a function in module A, would not list B. The selection is therefore widened by one import hop, read from the `import` statements of every file under `src/` when a changed file needs it (a file that imports a changed module counts as changed; a changed `__init__.py` reaches everything importing inside its package; test files count, so a test that only imports a constant is caught), and the plugin always treats these as changes that select everything: `conftest.py`, `src/tests/*support*.py`, `pytest.ini`, `setup.py`, `requirements*.lock`, anything under `src/tests/fixtures/`, and a test the map has never seen (that test alone is selected).

**Selection.** `--changed [BASE]` takes the files that differ from BASE (default: the merge base with `origin/main`, plus uncommitted and untracked files), and runs: every test in a changed test file, every test that is new to the map, every opaque test whose area changed, and every test whose recorded set (widened) contains a changed file. Non-code changes (docs, `changes/`, `TODO.md`) select nothing, except `docs/` files that a test reads, which the audit hook already finds.

**Keeping the map fresh.** A selected run re-records the tests it ran and updates their entries; the entries of tests not run are kept (they depend on none of the changes). `scripts/test.py full` records on every run (`--record-deps`), so the map stays current without anyone asking; a test that failed or skipped is not recorded and is therefore selected next time. A lane's first `changed` run with no map runs the full suite and records, so no one has to build a map by hand.

**Measured selectivity** (sample of 72 test files, 3,796 tests, browser and the heaviest files left out) [M]:

| Per changed source file (210 files) | Mean | Median | 90th percentile | Files hit by over half the tests |
|---|---|---|---|---|
| Functions that ran | 9.5% | 0.8% | 67% | 23 |
| Plus one import hop | 17.6% | 8.5% | 76% | 33 |
| Static `import` graph (not used) | 55% of test files | n/a | n/a | n/a |

A test ran 20.3 source files on average. The handful of files every test touches (`util/draw.py`, `physics/orbits.py`, `generation/config.py`, and `admin/throttle.py`, which an autouse fixture in `conftest.py` imports) are the widely shared core: a change there selects most of the suite, which is correct. A later refinement can record the fixture phase apart from the call phase, so an autouse fixture's imports stop counting for every test.

**Proving it before it is trusted.** The full suite stays the merge gate. For a number of merges (suggested: the next 20, or two weeks) each lane runs `changed` first and then the full suite, and the plugin reports any test that failed in the full run but was not selected; any such test is a hole in the map or the rules, and is fixed. Only when that list stays empty does Boss decide whether `changed` may replace the full run for small changes, with the scheduled daily full CI run as the backstop. Until then `changed` is a fast loop, not a gate.

Why not an existing tool: `pytest-testmon` does the same job with `coverage.py` (a C tracer, about 2 to 3 times slower while recording, one more dependency to pin) and does not see the data files or the subprocesses either. `sys.monitoring` plus the audit hook is smaller, and its overhead is nothing. If the plugin proves too much to maintain, testmon is the fallback.

## 4. Faster

In order of return per unit of risk. Times are from section 1.

1. **Browser phase** (42%). Run it apart from the rest, with 2 workers so SwiftShader is not starved, after the other tests; skip it entirely unless a change touches `web/`, templates, `static/`, `html/` or the shared map code (selection does this). Expected: the 13 load timeouts go away (TEST.126), and the all-core phase stops sharing CPUs with Chromium. A second step shares one page load across the several tests that load the same page (the fixture-map file loads a fresh page 86 times; its setup alone is 792 s).
2. **Database per worker, reset between tests** (about 10 to 13%: the 1,926 s of setup and teardown in the database files). A `mysql_config`-compatible fixture that takes a ready database from a per-worker pool, and on teardown deletes the rows from the tables the test filled (found with a `SELECT 1 ... LIMIT 1` per table, about 15 ms for the 122) instead of dropping and re-creating the schema. Tests that change the schema (old-schema fixtures, migrations, version detection) keep the current fixture. This is the one stage that touches test behaviour, so it comes with a leak check that compares a reset database to a fresh one table by table.
3. **The heavy generation files** (`test_neighborhood_edge` 753 s, `test_galaxy_gen` 628 s, `test_parallel_generation` 329 s, `test_prevalence` 319 s, `test_phenomenon_scatter` 260 s, `test_work_queue_failures` 239 s): about 2,500 s together. File by file, shrink the galaxy each one builds to the least that still shows the behaviour (the `sparse` neighbourhood case alone is 244 s), and share a built galaxy among the tests of a file when they only read it. These are the files the Foundations lanes edit, so each is done with its owner's word.
4. **Cheap wins**: `--dist worksteal` and the 600 s timeout are already on; set `-p no:cacheprovider` in the runner; collect once (collection is 5 s, run it with `--import-mode=importlib` only if it measures faster).

A rough result, adding the stages: the full suite from 55 to about 35 to 40 minutes (browser phase, database reset, three heavy files); a typical `changed` run to a few minutes, because most changes select under a fifth of the non-browser tests and no browser tests.

## 5. What stays and what does not move

- **Test files do not move.** A folder-per-area layout would make every other lane's open PR conflict, and the markers and the map give the same selection. If Boss wants the folders anyway, it is one mechanical PR at a quiet moment, last.
- **The merge gate** (full local suite, no waiting for CI) and the daily CI run are unchanged.
- **Flaky-test items** (TEST.111/116/125-132/137/138) stay with the Bugfixes lanes. This work removes one cause of them (the load) and does not edit their tests.

## 6. Stages

| Stage | Change | Touches | Risk |
|---|---|---|---|
| 1 | `db` marker fix, `heavy` and area markers (`areas.py`), `scripts/test.py` with `quick`, `area`, `browser`, `full` (two phases), docs/testing.md | conftest.py (marker hook only), new files | low |
| 2 | `depmap.py` plugin: `--record-deps`, `--changed`, the fixed select-everything rules, the proof report | new file plus two lines in conftest.py | low (opt-in) |
| 3 | Database per worker with reset, leak check | conftest.py fixture, new fixture module | medium |
| 4 | Browser page-load sharing; heavy generation files, one PR per file with its owner | test files | medium, per file |
| 5 | After the proof period: Boss decides whether `changed` may replace the full run for small changes | policy | n/a |

Each PR is merged on a green full run and reports the TODO items it resolves and adds.

## 7. For the TODO thread

- TEST.139 (this item) gets the stages in section 6 as sub-items.
- New bug: `markdown` 3.11.1 (installed by `pip install -e ".[test,api]"` today, since `setup.py` only says `markdown>=3.6`) removes the inline pattern `not_strong` that `web/lib/mdrender.py:34` names, so 68 tests fail in the first baseline run (74 when the database is on and the seven files that render Markdown run together: test_fuzz_text_rendering, test_web_system_phen, test_bughunt_render_edges, test_mdrender, test_web_population, test_systems, test_page_shell) with `ValueError: No item named "not_strong" exists`, and the site's Markdown pages would fail; `requirements.lock` pins 3.11, which works. Fix: cap the requirement (`markdown>=3.6,<3.11.1`) or make `mdrender.py` tolerate the missing name. (Confirmed: the same seven files give 466 passed on 3.11 and 74 failed on 3.11.1.)
- New: the `db` marker misses 2,210 tests (section 1); fixed in stage 1.
- `playwright` is unpinned (`>=1.49.0`) and a new release needs a newer Chromium than the one the cloud image carries (1.64 wants build 1248, the image has 1194). Lanes need `playwright==1.56.0` there. Worth a line in docs/testing.md.

## 8. Questions for Boss (with the default I will use)

1. **Run the browser tests as their own phase with 2 workers, and only when web or map code changed, in `changed` runs?** Default: yes. The `full` command still runs them.
2. **Keep test files where they are and group by markers?** Default: yes (section 5).
3. **Let `changed` replace the full run for small changes after the proof period?** Default: no change to the gate until the proof report is clean; I will ask again then.

## 9. Stage 1 result (2026-10-14, `scripts/test.py full`, started 13:47Z)

| Phase | Tests | Result | Time |
|---|---|---|---|
| Non-browser, 4 workers | 15,123 | 15,115 passed, 8 skipped, 0 failed | 21 min 33 s |
| Browser, 2 workers | 571 | 513 passed, 52 skipped, 6 failed | 34 min 43 s |
| Whole run | 15,694 | | 56 min 24 s |

- **The wall time did not shrink** (55 min 12 s before). The machine now does the same work without the contention: the browser phase on 2 workers spends 2,083 s of wall time where the browser files spent 5,390 s of summed test time before, but two cores sit idle during it. Stage 1 is about reliability and the fast loops, not the full-run time; that comes from stages 3 and 4 and from `changed` runs.
- **Browser failures went from 13 timeouts to 6**: four fail alone on main (TEST.125 three moon tests, TEST.142 the `after` argument), two are load failures that pass alone (`test_bookmarks_keep_a_course_pick_on_the_galaxy_map`, `test_galaxy_map_buttons`), confirmed by running them alone. Two workers is better than four but not enough; stage 4 gives the heavy WebGL tests state-based waits.
- Non-browser tests have no failures at all in this run (the 14 of the baseline were 13 browser timeouts and my own database grant).
- `quick` (not database, browser or heavy) is the fast loop; its time is printed by the runner.

## 10. The test fixes: decisions (Boss approved 2026-10-14 15:06Z)

Boss (14:01Z) moved every test-fix item to this rebuild and asked that none be implemented until it was decided which should be. Decided 2026-10-14; each item has a decision, a reason and a close condition. The TODO items keep the IDs; the owner of every one is the Test Suite Rebuild thread. The text of the decision list as posted is also kept in `/mnt/project-files/test-suite-rebuild/test-fix-decisions.md`.

### 10.1 Implemented as small fixes (group A)

These fail repeatably (alone, or on one database engine), not only under load, so they do not wait for the rebuild. One small PR each, in this order.

| Item | Decision | Reason | Closes when |
|---|---|---|---|
| TEST.125 three browser-map tests, "no generated system has a moon" | Implement | Failed in the 2026-10-14 baseline and fails alone and on main: the map fixture's galaxy has no moon since the generator's draws shifted. | The fixture always holds a planet with moons (a system generated last, with moons forced, so no earlier draw moves) and the three tests pass in a full browser phase. |
| TEST.142 fake `get_systems` lacks `after` | Implement | Failed in the baseline and alone; PERF.75's keyset paging added the argument. | The fake takes and honours `after`, and no other fake of a paged reader is left (checked with a grep). |
| TEST.143 install test fails from a git worktree | Implement | The test's stand-in Python imported `planetgen` from the machine's editable install, not from the checkout under test. Reproduced from a worktree. | The stand-in is pinned to the checkout with `PYTHONPATH`; the test passes from a worktree and from the main checkout. |
| TEST.129 orbit-update tests on MySQL 8.4 (`last_updated_at` rounds up) | Implement | MySQL 8.4 rounds fractional seconds where MariaDB truncates; may give a negative elapsed time in a real run. | The orbit update stores and compares a time both servers treat alike, checked on MySQL 8.4 and MariaDB. |
| TEST.130 scene vs stored positions for comets | Implement | A real disagreement for some random systems with comets (about 1e10 km). 400 seeded systems with comets and a second star did not reproduce it; systems with three stars are being tried next. | The cause is found and fixed, and the test draws a system that has comets (and the star layout that triggers it) every time instead of by luck. |

### 10.2 Folded into the rebuild, not fixed one by one (group B)

| Item | Decision | How the rebuild settles it | Closes when |
|---|---|---|---|
| TEST.126 WebGL tests time out under 4 workers | Fold into stage 1 and stage 4 | 13 of the baseline's 14 failures were this. Stage 1 runs the browser tests as a separate phase on 2 workers (6 failed in its first full run, 4 of them TEST.125 and 142; 2 pass alone). Stage 4 gives the heavy WebGL tests state-based waits. | 5 consecutive full runs show none of them failing. |
| TEST.138, TEST.144 browser menu and button tests fail once under load | Fold into TEST.126 | Same family. | With TEST.126. |
| TEST.111, 116, 127, 128, 132, 134, 137, 140, 141: failed once in a full run, passed alone, no traceback kept | Fold; no one-by-one hunt | None has a repro (TEST.111: none in 15 loaded runs; TEST.116: none in 83). Stage 1 keeps every failure's full report and re-runs it alone, saying "passes alone" or "FAILS ALONE". Stage 3 (one database per worker) removes the shared-database and create-and-drop-under-load causes. | 5 consecutive full runs after stage 3 show none of them; or a recurrence is fixed from its captured report. |
| TEST.116 specifically | Stays open for its first recurrence | The test pins `--star-type K2V` and a different stored type may mean the pin was not honoured on that draw, which would be a product bug, not test luck (inferred, not confirmed). | Its first recurrence is investigated with the captured seed. |

### 10.3 New work with an order (group C)

| Item | Decision | Closes when |
|---|---|---|
| TEST.133 big-galaxy query-budget test (EXPLAIN every page query on 2,000,000 systems) | Implement as its own opt-in tier (marker `budget`, `scripts/test.py budget`, run by hand before a PERF merge). Not in the merge gate and not in CI: building 2,000,000 systems takes about 2 minutes, and nothing here may raise CI cost. After stage 3. | The tier exists and passes on 2,000,000 systems. |
| TEST.135 nebula reads in that budget | Same tier, after TEST.133 and DB.24 (its stated prerequisite). | The nebula reads are in the budget and pass. |
| TEST.131 gravity evaluator tests | Stays with the gravity work: it needs GEN.198, which does not exist yet. Whoever builds GEN.198 writes these in the same PR as a normal test file, which the markers and the dependency map pick up with no change. | GEN.198 lands with its tests. |

### 10.4 Order

1. Stage 1 (merged): markers, runner, browser phase, failure capture.
2. Group A, one small PR each: TEST.125, 142, 143, then 129, 130.
3. Stage 2: change-based selection.
4. Stage 3: one database per worker; then re-check group B.
5. TEST.133, then TEST.135 when DB.24 lands.
6. Stage 4: browser page-load sharing and the heavy generation files, one PR each with the owning lane's word.

A lane that sees a test failure still reports it to the TODO thread; with the stage 1 failure capture the report is kept.


## 11. Stage 2 result: change-based selection (2026-10-14)

`src/tests/depmap.py` is loaded by `conftest.py` (`pytest_plugins`). It adds:

- `--record-deps`: writes `.pytest-deps.json`; on in `scripts/test.py full`.
- `--changed [BASE]` (`scripts/test.py changed [--base REF]`): runs the same two phases as `full`, but only the selected tests; prints how many of how many and the recorded seconds skipped. No map yet: runs everything and records one.
- `--deps-check [BASE]` (`scripts/test.py full --check`): runs everything, then prints `MISSED:` for each failed test `--changed` would not have chosen, and appends the run to `.pytest-deps-proof.json`.

Checked on a small file set: touching a file nothing used selected 0 of 43; touching `names/__init__.py` selected the 22 tests that import inside the package; touching `conftest.py` or the plugin selected all 43. A collect-only run no longer overwrites the map. The first full run with `--check` records the first complete map; the proof period (section 3.3) starts from the full runs after that, and `changed` stays a fast loop, not a gate, until Boss decides.

## 12. Stage 3 result: one database per worker (2026-10-14)

`src/tests/dbpool.py` and the new `mysql_config` (section 4, item 2). Measured on the 3,243 database tests (`-m "db and not browser" -n 4`), 16 min 49 s with the deep check on every reset (15 min 48 s without it); on four database files, 26 s became 14 s.

- What a reset does: empties the tables with rows (one `UNION ALL` probe, then `DELETE` with foreign keys off), sets an advanced `AUTO_INCREMENT` back to 1 (so a test still gets row id 1), puts back the baseline `schema_migrations` row, closes the pool, and drops the process's cached facts about the database name. The caches are found by name, not listed: every `_`-prefixed dict or set in a `planetgen` module is searched for a key naming the database (this found `_id_blocks`, `_facet_cache`, `_token_sizes`, the count cache, the queue connections); the Redis keys holding the name are deleted. A new cache keyed by database is covered with no change here.
- A test that changes the schema is noticed by comparing table, column, index, routine, trigger and constraint structure and each table's creation time with how the test began; the database is dropped and the next test gets a new one.
- Tests that need no schema at all are marked `fresh_database` and get the old fixture: `test_db_old_schemas.py`, `test_db_schema_version_detect.py`, `test_facility_velocity_migration.py`, `test_generation_size_estimate.py` (it measures table bytes, which a reused database has grown), the v63 migration test, and three "uninitialized database" tests. The first run with the reuse turned on had 55 failures; they were these and one test that read a count another test had left in Redis under the same database name.
- The proof: `PLANETGEN_TEST_DEEP_CHECK=1` compares each reset with a new database (every table definition, rows, counters). All 3,243 database tests passed with it on. `test_dbpool.py` repeats the check on a real sector save and covers the schema-change and counter cases.
- `PLANETGEN_TEST_FRESH_DB=1` returns every test to its own database.
- **Measured on a full run** (`scripts/test.py full --check`, started 18:47Z): the non-browser phase went from 21 min 40 s (stage 2 run) to 19 min 23 s, a saving of 2 min 17 s, not the 8 minutes the setup and teardown total suggested: the database files are dominated by the work the tests do, and the heaviest files run on four workers either way. The browser phase (39 to 40 min on 2 workers) is now two thirds of a full run, so it is the next lever (stage 4). The full run took 59 min 43 s (60 min 54 s before). One failure, `test_ratelimit_redis.py::test_two_apps_on_one_redis_share_the_request_counts`, which passed alone: its fixed one-minute request window can end between its requests (the limiter counts from the clock minute), a rare time flake; reported for the TODO list.

## 13. Policy: `changed` as the merge gate (Boss, 2026-10-14 19:44Z)

`scripts/test.py changed` may replace the full local suite as the merge gate for small changes and for the steps of a sequence of related items. A full run (`scripts/test.py full --check`) is for when a sequence is finished, not for every merge. The safeguards are in the selection (a change to the test setup, or no map yet, selects everything), in `--check` (every full run lists any failing test that `changed` would have missed; that is a bug to fix), and in the daily CI run on main. A PR description says which run it used.

## 14. TEST.133: the big-galaxy query budget tier (2026-10-14)

`scripts/test.py budget` (`src/tests/test_budget_queries.py`, `budget_support.py`). Opt-in, not in the merge gate or CI (decision list, TEST.133).

- The galaxy: `ensure_big_galaxy` clones one real saved single-star system into 2,000,000 systems with stars on 150,000 cloned sectors, with a real-galaxy mix of star types (three in four red dwarfs). The full-text indexes are dropped while filling and rebuilt after, then `ANALYZE`. The database is named for its size and the schema version and kept between runs.
- The measure is rows read by the server for one call (the session's `Handler_read*` counters before and after), not seconds, so the budget holds on any machine. A page of 50 may read 5,000 rows; a count may read about the table.
- Cases are the project's own query functions: the Systems and Sectors pages (several sorts, a star-type filter, offsets 0, 20,000 and a quarter of the table) and the counts. TEST.135 adds the nebula reads once DB.24 is in.
- First run on 200,000 systems found: sort by sector and by octant are now within budget (PERF.70 holds); the Sectors list sorted by name, systems or density, the Systems page filtered by star type, and deep offsets by name read far more than a page and are expected failures until their fixes land.
- `test_budget_support.py` builds 3,000 systems and runs every case in the ordinary run, so a schema or query change that breaks the tier is seen at once.
- **Measured on the full size** (2026-10-15, MariaDB 10.11): building 2,000,000 systems took 474 s (397 s of it the system rows, whose many indexes are the cost), 45 s for the stars, 22 s to rebuild the full-text indexes; the 13 cases then ran in about 60 s: 7 within budget, 6 expected failures. Later runs reuse the database.

## 15. A clear last line from `scripts/test.py` (Boss, 2026-10-15 00:09Z)

Boss's rule: a run must never be cut by the 30-minute background limit, and it must be certain when a run is done. Start every long run with `run_in_background` and a timeout of at least 90 minutes (7,200,000 ms), with its output and start time written to a file. `scripts/test.py` now passes pytest's output through and ends with exactly one of:

- `scripts/test.py <command>: PASS in <time>: <counts>`
- `scripts/test.py <command>: FAIL in <time>: <counts>` (the failed tests were already listed with their verdicts above it)
- `scripts/test.py <command>: INTERRUPTED (<signal>) in the <phase> phase after <time>; not finished, <counts> so far` on Ctrl-C or a SIGTERM; the child pytest is stopped first.

A run whose output ends with none of the three was killed (SIGKILL or the background limit, which cannot be caught) and counts as not run.

## 16. Flakes that were not load (2026-10-15)

The first capture with full reports (`/mnt/project-files/test-suite-rebuild/failure-records/2026-10-15-0103Z-changed.jsonl`) showed two "load" failures were random data. Both reproduce on an idle machine in a loop:

- TEST.149 `test_only_due_planets_and_moons_move_and_count`: 4 of 150 random systems have a slow, distant moon whose 100,000 km takes longer than a year, so it is correctly not due. Fixed in the test (expects the due count it computes, checks the rest did not move); 0 failures in 400.
- TEST.147 the nebula field table tests: the galaxy seed was random per database and 5 of 80 seeds had no class M or N cloud reaching the sampled sectors. Fixed with a fixed seed.

Lesson for the rest of the flake list: before calling a failure load, run it in a loop on an idle machine; a test that builds random data without a seed can fail at any load.
