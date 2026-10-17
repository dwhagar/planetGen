# Testing

The suite lives in `src/tests/` and runs with plain `pytest` from the
repository root (`pytest.ini` puts `.`, `src` and `src/html` on the path).

```sh
pip install -e ".[test,api]"
python -m nltk.downloader words
pytest -n auto
```

## Running in parallel

`-n auto` (pytest-xdist, part of the `test` extra) runs one worker process
per core. Plain `pytest` still runs everything serially in one process,
which is easier to read when debugging a single failure. On a 4-core
machine with the database tests on, the whole suite took 17 min 21 s
serially and 5 min 57 s with `-n 4` (2026-10-01).

Workers never share database state:

- each test process keeps its own uniquely named database for `mysql_config` (see "The database tests");
- each worker process gets its own control database
  (`PLANETGEN_CONTROL_DATABASE` is set to `pgtest_control_<worker>_<random>`
  by `conftest.py`, overriding any value in your environment, and dropped
  at the end of the run), so nothing in the suite touches
  `planetgen_control` or a real control schema;
- generation runs in-process (`PLANETGEN_WORKERS=1`) unless you set
  `PLANETGEN_WORKERS` yourself.

## Generation at several worker counts

A real run hands its sectors and bright-star layers to worker processes.
CI runs the tests marked `generation` (the files that run `planetgen` or
the work queue; `conftest.py` adds the marker) again with 2 and 4 workers
(TEST.74). To do the same locally:

```sh
PLANETGEN_WORKERS=2 python3 -m pytest -n auto -m generation
```

`monkeypatch.setattr` only changes the test's own process, never a spawned
worker. A test that patches something a generation run calls uses
`tests/worker_patches.py` instead: `patch_everywhere(monkeypatch, module,
name, "tests.my_test:factory", **params)` patches this process and every
worker started while the patch is in place (`$PLANETGEN_WORKER_HOOKS`).
The factory is a module-level function that builds the replacement, and
shared state such as a call count goes in a file (`worker_patches.bump`).
With more than one worker, the sectors or layers already handed out when
a fault hits can still finish, so a test checks "at least" where one
worker would stop exactly.

`--dist worksteal` (what CI uses) lets idle workers take queued tests from
busy ones, which helps when a few slow tests land on the same worker.

## The database tests

Tests that save and load systems, sectors or users need a MySQL or MariaDB
server. They skip (not fail) without one. Point them at a server with:

```sh
export PLANETGEN_TEST_MYSQL_HOST=127.0.0.1
export PLANETGEN_TEST_MYSQL_PORT=3306
export PLANETGEN_TEST_MYSQL_USER=root
export PLANETGEN_TEST_MYSQL_PASSWORD=...
```

Each one falls back to the matching `PLANETGEN_MYSQL_*` variable, then to
the built-in default (`127.0.0.1:3306`, user `planetgen`).

The user needs to create and drop databases. Each test process (each
`pytest -n` worker) keeps one uniquely named database with the schema in it,
and `mysql_config` hands it to one test at a time, so the database tests are
safe under `pytest -n auto` (see "Running in parallel" above). After each
test `src/tests/dbpool.py` puts the database back the way a new one is: the
tables the test filled are emptied, the auto-increment counters it advanced
go back to 1, the baseline `schema_migrations` row returns, and the
process's cached facts about that database name (id blocks, stored counts,
facet caches, what was put in Redis under the name) are dropped. This
replaced creating, filling with the schema and dropping a database for every
test, which cost 0.66 s each on an idle machine.

- A test that changes the schema leaves a database that cannot be reused;
  the pool notices (table, column, index, routine, trigger or constraint
  added, dropped or altered) and drops it, and the next test gets a new one.
- A test that needs a database with **no schema at all** (it checks what a
  first connection does, or builds an old schema and migrates it) takes the
  old fixture with `@pytest.mark.fresh_database`, or `pytestmark` for a
  whole file. So does a test that measures the bytes a fill adds, which a
  reused database has grown.
- `PLANETGEN_TEST_FRESH_DB=1` gives every test a new database again: the
  way to tell a problem in the pool from a problem in a test.
- `PLANETGEN_TEST_DEEP_CHECK=1` makes every reset compare each table's
  `SHOW CREATE TABLE`, the tables holding rows and the auto-increment
  counters with a new database, and fail the test that left a difference.
  Run it after changing the schema or adding a cache keyed by database.

## The math check runs first

`src/planetgen/physics/mathcheck.py` is a fixed list of checks that the
generator's math gives the right answers (TEST.63): known values from real
astronomy (the Sun, Earth's and Jupiter's orbits, the habitable zone, white
dwarf sizes, Holman & Wiegert's stability limits, the Kepler and Barker
equations), identities that hold for any input (unit round trips, the
constants agreeing with each other, energy conservation around a Kepler
orbit, the sector grid tiling every shell), and seeded draws from every
sampler checked against its intended shares with a chi-square test. Each
check names the function it calls, the expected value, the tolerance and
where the expected value comes from. It takes well under a second.

- Before any test runs, `conftest.py` runs the whole list once; if a check
  fails, the run stops there with the failed checks listed, since every
  other failure would be noise.
- `test_math_check.py` repeats each check as its own test (marker
  `mathcheck`), sorted to the front of the run.
- CI runs it as its own first job, `mathcheck`; every other job waits for it.
- The website runs it once when it starts and shows admins a warning on
  every page if it failed (the site keeps serving).
- Bulk generation runs it first and refuses to start if it fails:
  `planetgen check-math` by hand, and every bulk path (see
  [cli.md](cli.md#subcommands)); `update.sh` warns.

```sh
python -m planetgen.physics.mathcheck -v            # the report, exit 1 on failure
pytest -m mathcheck                               # just the math check tests
```

A check that fails after a deliberate change to a constant or formula
means the check's expected value, its source or the change needs a second
look; update the check only with a source for the new value.

## Picking a slice of the suite

`conftest.py` marks every test (markers registered in `pytest.ini`):

- `db`: uses a real database (the `mysql_config` fixture, directly or
  through another fixture);
- `slow`: the brute-force and seeded-sweep files (`test_fuzz_*`,
  `test_bughunt_*`);
- `browser`: the headless-browser checks (`test_web_a11y.py`,
  `test_web_browser_layout.py`, `test_web_browser_maps.py`,
  `test_web_browser_fixture_maps.py`);
- `mathcheck`: the math check gate (`test_math_check.py`,
  `test_web_math_check.py`), run first.

- `heavy`: the files that took over 100 s each in the 2026-10-14 baseline
  (`HEAVY_FILES` in `src/tests/areas.py`);
- `area_web`, `area_api`, `area_db`, `area_generation`, `area_physics`,
  `area_names`, `area_wiki`, `area_ops`: the part of planetGen a test file
  covers, by its file name (the rules are in `src/tests/areas.py`; no file is
  moved, and a name no rule matches is `area_other`).

```sh
pytest -n auto -m "not db and not slow"   # quick loop without a database
pytest -n auto -m db                      # just the database tests
pytest -n auto -m area_api                # one part of the program
```

### The slice picker, `scripts/test.py`

`scripts/test.py` builds those pytest commands for you (TEST.139,
[design/test-suite-rebuild.md](design/test-suite-rebuild.md)):

```sh
scripts/test.py quick           # not database, browser or heavy tests
scripts/test.py area web db     # the tests of those areas
scripts/test.py browser         # the headless-browser tests alone
scripts/test.py full            # everything: the merge gate
scripts/test.py full -- -x      # anything after -- goes to pytest
```

Every command but `quick` and `browser` runs two phases: the tests that need
no browser on every core, then the browser tests on 2 workers
(`--browser-workers`). Software WebGL needs the CPU the other workers would
have, and running both together is what timed the browser tests out
(TEST.126). The first line printed is the start time and how long the same
command took last time on this checkout (kept in `.test-runs.json`, not
tracked by git); the last line is how long it took now.

**Failures.** The runner keeps every failing test's full report (traceback,
worker, the machine's load average) in `.test-failures.jsonl`, and when the run
ends it runs each failed test again alone, on one process. The verdict is
"passes alone" (it failed only under load or beside other tests) or "FAILS
ALONE" (a real failure). A lane quoting a flake in a TODO item should paste
that line and the record from `.test-failures.jsonl` instead of "failed once,
passed alone". `--no-rerun` skips the second step.

A pull request needs a green `scripts/test.py changed` (below) or, at the end of
a sequence of related items, a green `scripts/test.py full` (or `pytest -n auto`
as before).

### Telling when a run is done

`scripts/test.py` always ends with one line: `PASS in <time>: <counts>`, `FAIL in <time>: <counts>`, or `INTERRUPTED ...`
(Ctrl-C or a SIGTERM; it stops pytest first). Output that ends with none of these was killed. A full run takes about
35 to 50 minutes (34 min 32 s on a 4-core machine, 2026-10-15), so start it in the background with a timeout of at least
90 minutes (7,200,000 ms; the default 30 minutes cuts it), with its output and start time written to a file, and read the last line of that file rather than
looking for processes.

### A test that hangs

Every test has a time limit of 900 s (`timeout` in `pytest.ini`, pytest-timeout, method `signal`; the slowest
test recorded took 479 s under load). A test that runs past it is stopped and fails with "Timeout (>900.0s) from
pytest-timeout", which stops one test and not the worker or the run. The runner lists such tests by name under
"TIMED OUT" before its last line and counts them there. A test that needs longer says so with
`@pytest.mark.timeout(seconds)` (the budget tier gives itself an hour for building its galaxy). A helper thread
in a test that waits for something must still do what lets the test go on if the wait fails
(`test_work_queue_admin.py` `_ResumeLater`), and raise the failure from the test.

### The progress log and the build id

Every `scripts/test.py` run works out the id of its build (printed on the first line: `build 1a2b3c4d5e6f7a8b`) and
keeps `.test-progress/<id>.jsonl`: a header with the parts the id is made of and a manifest of the files (path and
hash, for comparing two runs' trees), then one line for every test that finished (node id, outcome, seconds, worker)
and, for a passed test that recorded its dependencies, a `deps` line with the files it used, written as the test
finishes, so a killed run loses only the tests that were running. `--resume` reads it (below). The id changes when the
git tree, an uncommitted or untracked file (git-ignored files do not count), an installed package, the Python,
MySQL or Redis version, the command and its pytest arguments, the commit `changed` compares with, a `PLANETGEN_*`
or `PYTEST_*` variable or the UTC date changes, and not otherwise. A new run of the same build starts its log
again; logs older than a day are removed.

### Resuming a run

After a run was killed, hung or stopped (its last line says `INTERRUPTED`, or there is no last line), run the same
command again with `--resume`:

    scripts/test.py full --resume

The tests that passed in the earlier run are skipped (`--skip-passed`, a list of node ids) and the rest run. A test
is skipped when it passed in a run of the same build, or (TEST.208) in an earlier run, of any command, that kept the
same environment (installed packages, Python, MySQL and Redis versions, `PLANETGEN_*` and `PYTEST_*` variables, UTC
day) and none of the files that differ between that run's tree and this one is a file the test used. Which files a
test used is what the dependency map records for `changed` (its own test file, every Python file it ran, every file it
opened, and for a test that starts a subprocess all of `src/planetgen/`, `src/html/` and `scripts/`); it is written to
the progress log as a `deps` line when the test passes, so a cancelled run keeps it. The rules of `changed` decide:
a changed test file, a changed file the test used or that imports one it used, or a change to the test setup
(`conftest.py`, the `*_support.py` modules, fixtures, `pytest.ini`, the lock files) reruns it, the last one all of
them. So after a gate was cancelled and main was merged in, only the tests the merge touched run again. A pass that
has no `deps` line (the run was killed between the two lines, or `quick`, `area` or `browser` ran it, which record
no dependencies) is reused only by a resume of the same build.

It prints `RESUMING build <id>: the N tests that passed earlier are skipped (K in the run started <time>[ of build
<id>]; ...), the rest run now` first, and its last line ends `(RESUMED build <id>: N more passed earlier, K in the run
started <time> ...)`, so the verdict can be read as: the tests run now passed, and N more passed before, in those
runs. The reused passes are copied into the new build's log, so a second resume counts them too. A resumed run's
duration is not stored as the command's usual time. A different environment (a package, Python, a server, a
`PLANETGEN_*` variable or the day) reuses nothing: the run says `RESUMING build <id>: nothing from the earlier runs
can be reused (<why>); every test runs` and runs every test. When no earlier run of the command left a log it runs
nothing, exits with 2 and says so. Resuming is never automatic and not part of the merge gate, and a PR that was
merged on a resumed run says which build and runs it continued. Tests that failed, errored, were running when the run
stopped or never started run again; a test that passed earlier is not run again, so a test that only passes after
another test ran before it is not checked again either.

#### What a run clears first, and what a resumed run records

Before pytest starts (every command, a resumed run too) `scripts/test.py` clears what earlier runs left on the test
Redis (`src/tests/redis_cleanup.py`). A queued job leaves its job record, result and a place in the queue's finished,
failed or canceled registry, none of which expire, so the server gained a few hundred keys a run. The clean-up uses
only `PLANETGEN_TEST_REDIS_URL` (never `PLANETGEN_REDIS_URL`) and removes only keys a test makes: `planetgen-test:*`,
the API test queues `planetgen-api-test-*`, the run queues `planetgen-run-*` and queues `test-*`, with the job
records those queues list. It does not flush the database and leaves the rate-limit keys (`LIMITS:*`, which expire on
their own and are shared by every running test). A run queue of the application is named like a test's, so if any
worker has sent a heartbeat in the last 90 seconds the clean-up does nothing and says so, and a queue that lists a
job made in the last 15 minutes is left whole (another run may be queueing for workers that have not yet started: its
test would end with `WorkerDied: the task's job ended canceled`, TEST.203). It prints
`removed N keys earlier test runs left on the test Redis` when it cleared something.

A resumed run with `--record-deps` or `--check` keeps the dependency-map entries of the tests it skipped and
replaces those of the tests that ran. Its `--check` output says it was resumed and that the check covers only the
tests that ran now, and the record in `.pytest-deps-proof.json` carries `"resumed": true`.

### Only the tests a change touches, `scripts/test.py changed`

`scripts/test.py full` also records, for every test that passed, which
repository files it used (Python files, templates, static scripts, `schema.sql`,
anything it opened; a test that starts a subprocess is marked as depending on
all of `src/planetgen/`). The record is `.pytest-deps.json` (not tracked by git;
recording needs Python 3.12 or later). Then:

```sh
scripts/test.py changed             # tests the changes since origin/main can affect
scripts/test.py changed --base HEAD # only the uncommitted changes
scripts/test.py full --check        # full run; lists failures `changed` would have missed
```

It runs every test in a changed test file, every test that used (or whose
used file imports) a changed file, every test that starts a subprocess when
`src/planetgen/`, `src/html/` or `scripts/` changed, and every test the record
has not seen or saw fail. A change to `conftest.py`, `pytest.ini`, a support
module or the test fixtures runs everything, and so does a missing record
(which that run then writes). A green `changed` run is the merge gate for small changes and for the steps
of a sequence of related items (Boss, 2026-10-14); the whole suite
(`scripts/test.py full --check`) is run when a sequence is finished, and the
daily CI run on main backstops both. A PR description says which one it used.
`scripts/test.py full --check` is also how the selection is proved: it
should never print a `MISSED:` line, and if it does, tell the test-suite thread.

### Time limit and knowing a run is done

Start every test run with a timeout of at least 90 minutes (Boss,
2026-10-15), so a slow machine or a full run never hits the 30 minute
default and dies half way. A run that finishes early is fine; the limit only
has to be long enough. Know a run is done from its own completion: the
notice the run sends when it exits, or its final summary line (`scripts/test.py`
prints one). Do not poll with `pgrep`, `ps` or `sleep` loops; a process that
looks gone may only be between phases. Resuming a stopped run from where it
left off, when the build is the same, is being investigated (see TEST.148 in
the TODO).

## Setting up a cloud session

A fresh container needs `pip install -e ".[test,api,browser]"`, the NLTK
corpus (`python -m nltk.downloader words`), MariaDB or MySQL with a user that
has `ALL ... WITH GRANT OPTION` (one test creates a restricted user), Redis
(`redis-server --daemonize yes`, then `export PLANETGEN_TEST_REDIS_URL=redis://127.0.0.1:6379/0`; the tests that
need it skip without that variable), and a Playwright whose Chromium matches the one on the machine: the image
carries Chromium 1194, so use `pip install playwright==1.56.0` (a newer
Playwright asks for a build the image does not have and every browser test
skips). `markdown` must be 3.11 as in `requirements.lock`; 3.11.1 fails 74
tests (see TEST.139 in the TODO).

To have all of that ready in every new session, paste
[`examples/cloud-environment/setup-planetgen-env.sh`](../examples/cloud-environment/setup-planetgen-env.sh) into the
cloud environment's Setup script field (it installs MySQL 8, Apache 2, Redis and every dependency, and leaves a
`planetgen-env-up` command that starts the services each session); its
[README](../examples/cloud-environment/README.md) says what each part does.

## Database engines

Local runs here use MariaDB 10.11. CI's `test` job has five legs, each with its own Redis container: the whole
suite on Python 3.12 against MySQL 8.0, MySQL 8.4 and MariaDB 11.4, then the `generation` tests alone with 2 workers
(MySQL 8.4, Python 3.12) and with 4 workers (MariaDB 11.4, Python 3.13). MySQL 8 is stricter about reserved words and
`GROUP BY` than MariaDB, so a query that works locally can still fail there.

## Kinds of test

- **Unit and behaviour tests** (`test_<area>.py`): one module or page each,
  with fixed inputs and exact expectations.
- **Bug-hunt tests** (`test_bughunt_*.py`, support in
  `bughunt_support.py`): seeded, exhaustive sweeps for the recurring bug
  categories (placement overlaps, unit conversions, NaN physics, name
  exhaustion, rendering breakage). Plain asserts are hard invariants;
  tests marked `tier2` (registered in `pytest.ini`) check softer,
  statistical properties but still run every time.
- **Edge-case tests** (`test_edge_*.py`): the admin scripts (`jobRunner`,
  `resetDb`, `migrateDb`, `updateOrbits`) run in-process against empty,
  broken and unreachable databases, cancelled jobs and malformed input.
- **Brute-force (fuzz) tests** (`test_fuzz_*.py`): property-based tests
  written with [Hypothesis](https://hypothesis.readthedocs.io/). Each one
  states a rule that must hold for every input (a generated planet is
  inside its class's ranges, a sector cell contains every point it
  samples, a web route never answers 500, a CLI option never prints a
  traceback) and lets Hypothesis search for a counterexample, including
  NaN, infinity, negative zero, huge and tiny numbers, empty and hostile
  strings. They cover:
  - `test_fuzz_galaxy_grid.py`: the shell/band/wedge sector grid,
    neighbours, designations, density, skeleton and map tiles.
  - `test_fuzz_system_generation.py`: every star type and config option,
    binaries, remnants, moons and belts.
  - `test_fuzz_planet_inputs.py`: every way of pinning a planet's class,
    radius and mass.
  - `test_fuzz_sector_placement.py`: sector filling, spacing, explicit
    positions and save/load round trips.
  - `test_fuzz_names.py`: name generation, validation and uniqueness.
  - `test_fuzz_text_rendering.py`: wikitext and Markdown output.
  - `test_fuzz_utils_and_config.py`: physics helpers, formatting, config
    loading, logging and progress files.
  - `test_fuzz_cli.py`: every `planetgen` command and option.
  - `test_fuzz_web_routes.py` and `test_fuzz_api_auth.py`: every web page
    and API endpoint, logins, CSRF and API keys.

- **Browser checks** (`test_web_a11y.py`): every page in headless
  Chromium at phone and desktop widths, in light and dark, against
  axe-core and the layout rules described in
  [`html-interface.md`](html-interface.md#flask-pages) ("Browser checks").
  They need `pip install -e ".[browser]"` and `python -m playwright install
  chromium` plus the MySQL test server, and skip without them. CI runs
  them in their own `browser-a11y` job.
- **More browser checks** (`test_web_browser_layout.py`,
  `test_web_browser_maps.py`), with the same needs and database: no two
  controls overlap or run off the screen on any page at 390, 600, 820 and
  1280 px in both themes (with folded sections open too); every map
  button changes the view (the 3D maps draw with SwiftShader WebGL, and
  a canvas screenshot that doesn't change is a dead button); the System
  Map's selection and measuring; and the Galaxy Map's drill-down by
  clicks, with Back, Forward and the free camera.
- **Map behaviour with no database** (`test_web_browser_fixture_maps.py`,
  TEST.70): the Galaxy Map and the Sector Map served by the real Flask
  views over fixture data (`map_site_support.py`: a real density shape,
  a few generated sectors near the core, one sector with stars across the
  luminosity range, a nebula, a neutron star and rogue planets), so it
  needs Playwright and Chromium but not MySQL. It pins what the shared
  map code (MAP.61) will move: picking (a click on a point of light, a
  drag that picks nothing), hover, keys, Back and Forward and the URL,
  bookmarks and their keys, and the scale line.
- **Page-script tests under node** (`js/*.test.mjs`, run by
  `test_js_unit.py`, skipped without node): `js/fakedom.mjs` is a small
  stand-in for the browser's DOM (elements, selectors, events, fake
  timers, history, fetch), enough to run the Galaxy Map's stage view and
  buttons, the Sector Map (with a renderer that draws nothing), the
  diagram's zoom, the Generate job panel and the facility form without a
  browser. `test_js_unit.py` hands them what they compare against from
  Python (the density shape, sector designations, the panels' buttons)
  in `PLANETGEN_JS_FIXTURES`. To run one by hand, from `src/`:
  `PLANETGEN_JS_FIXTURES="$(python -m tests.test_js_unit)" node --test tests/js/<name>.test.mjs`.

- **Census audits** (`tests/audit/census.py`, TEST.186): the one harness for tests that fill seeded sectors in each
  region class and compare what the fills made with an expectation. A `RegionCase` names a position in a galaxy shape;
  `run_census` fills its seeded sectors (the seeds are fixed by the region name and the sector number, so every run
  fills the same ones) and counts them by stage and kind from `sector.stage_counts`; `compare` tests each expected
  count within four Poisson sigma, applies the rule that a kind expected 20 or more times appears at least once, flags a
  kind found where its expectation is zero, and `format_report` prints the one report. An audit adds its expectations
  and checks and keeps no counters of its own.

A known bug that is not fixed yet is kept as a test marked
`xfail(strict=True)` whose reason names the bug and a repro. When someone
fixes it, the test starts passing, strict mode fails the run, and the
marker has to come off: the fix can't land without its regression test.

## Opt-in tests

These skip unless you turn them on:

- `test_wikiclient_wikijs_integration.py` and
  `test_wikiclient_mediawiki_integration.py` talk to a real wiki: set
  `PLANETGEN_TEST_WIKIJS_BASE_URL` and `PLANETGEN_TEST_WIKIJS_TOKEN`, or
  `PLANETGEN_TEST_MEDIAWIKI_BASE_URL`, `PLANETGEN_TEST_MEDIAWIKI_USERNAME`
  and `PLANETGEN_TEST_MEDIAWIKI_PASSWORD`.
- `src/tests/test_budget_queries.py` is the big-galaxy query budget (TEST.133): `scripts/test.py budget`. It
  builds (once; kept between runs) a database of 2,000,000 systems on 150,000 sectors by cloning one real saved
  system, about eight minutes and 3 GB, named for its size and the schema version so a new schema builds it again.
  Each case runs one page or list query and counts the rows the server read (the session's `Handler_read*`
  counters, which do not depend on how fast the machine is); it fails when a page reads far more rows than it
  shows. It is not in the merge gate or CI: run it by hand before a performance merge. `PLANETGEN_BUDGET_SYSTEMS`
  and `PLANETGEN_BUDGET_SECTORS` make a smaller galaxy for writing a case, `PLANETGEN_BUDGET_REBUILD=1` builds
  again. A case that is known to read too much is an expected failure naming the item that will fix it; when the
  fix lands the case fails as an unexpected pass, and its mark is removed. `test_budget_support.py` runs every
  case on a 3,000-system galaxy in the ordinary run, so the tier cannot rot unseen.
- `test_sector_generation_perf.py` is a timing benchmark:
  `PLANETGEN_RUN_PERF_BENCHMARK=1 pytest -s src/tests/test_sector_generation_perf.py`.

## Fuzz profiles

`src/tests/fuzz_support.py` sets how hard Hypothesis searches:

| Profile | Examples per test | Randomised | Used by |
|---|---|---|---|
| `ci` (default) | 60 | No: the same inputs every run | every `pytest` run and CI |
| `deep` | 3000 | Yes: new inputs each run | the weekly Deep fuzz workflow |

Choose one with `PLANETGEN_FUZZ_PROFILE`, and override the count with
`PLANETGEN_FUZZ_EXAMPLES`:

```sh
PLANETGEN_FUZZ_PROFILE=deep pytest src/tests/test_fuzz_*.py
PLANETGEN_FUZZ_PROFILE=deep PLANETGEN_FUZZ_EXAMPLES=20000 pytest src/tests/test_fuzz_galaxy_grid.py
```

The `ci` profile is derandomised so a pull request never goes red because a
new random input happened to be drawn. The deep run is where new inputs are
explored.

## When a fuzz test fails

Hypothesis shrinks the failing input to the smallest one it can find and
prints it, together with a `@reproduce_failure(...)` line. To fix it:

1. Copy the printed input into an `@example(...)` on the test (or a
   `parametrize` case), so it runs on every future run.
2. Fix the code, not the test, unless the rule itself was wrong.

## When CI runs

`.github/workflows/ci.yml` (the math check, the five test legs described under "Database engines", browser and
accessibility, the JS tests, the installer checks, the dependency audit) never runs on a push or a pull request. It runs:

- **On a schedule**, on main (every day at 07:00 UTC; the `cron` line in
  the file). The same scheduled run starts CI on every remote branch named
  `dev-*` (the `trigger-dev-branches` job).
- **By hand**: **Actions > CI > Run workflow**, on any branch except
  `claude/*` and `claude-*`. Those branches never run the test jobs: the
  run is created and every job is skipped.

A newer run on the same branch cancels the older one that is still going;
cancelled runs stay in the Actions history. The post-merge
`stamp-version.yml` and the PR check `release-note.yml` still run on their
own. Run the suite locally before merging (see above).

## The Deep fuzz workflow

`.github/workflows/deep-fuzz.yml` runs every fuzz file under the `deep`
profile each Monday, and can be started by hand from the Actions tab
(**Deep fuzz** > **Run workflow**, optionally with a number of examples).
A red deep run is a real bug report: its log carries the shrunk input.

## CI runners

CI can run on your own computers. Which jobs run where, what each machine
needs, the security settings and troubleshooting are in
[`ci-runners.md`](ci-runners.md).
