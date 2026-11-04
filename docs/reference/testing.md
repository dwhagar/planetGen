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

Tests that save and load systems, sectors or users need a MySQL
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

## The math check

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
- `update.sh` runs it after every update and warns when it fails;
  `planetgen check-math` runs it by hand, and the Admin page's "Check the
  math" button runs it on request (PERF.170). Generation does not run it.

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
  `area_names`, `area_ops`: the part of planetGen a test file
  covers, by its file name (the rules are in `src/tests/areas.py`; no file is
  moved, and a name no rule matches is `area_other`).

```sh
pytest -n auto -m "not db and not slow"   # quick loop without a database
pytest -n auto -m db                      # just the database tests
pytest -n auto -m area_api                # one part of the program
```

### The slice picker, `scripts/test.py`

`scripts/test.py` builds those pytest commands for you (TEST.139,
[design/operations/test-suite-rebuild.md](../design/operations/test-suite-rebuild.md)):

```sh
scripts/test.py quick           # not database, browser or heavy tests
scripts/test.py area web db     # the tests of those areas
scripts/test.py browser         # the headless-browser tests alone
scripts/test.py full            # everything: the end-of-version run, not the gate of an ordinary merge
scripts/test.py env             # only check (and repair) MySQL, Redis, the test variables and packages
scripts/test.py full -- -x      # anything after -- goes to pytest
```

Every command but `quick` and `browser` runs three phases: the tests that need
no browser and no real queue worker on every core, then the `real_queue` tests
(they start real worker processes, so two at once kill each other's workers)
on one worker, then the browser tests on 2 workers (`--browser-workers`).
Software WebGL needs the CPU the other workers would have, and running both
together is what timed the browser tests out (TEST.126). The two parallel
phases are also capped by what the machine has free when they start (TEST.344):
one worker per 1000 MB of available memory (`--worker-mb`; browser workers 2500
MB, `--browser-worker-mb`) and, in the first phase, no more than the cores the
current load leaves idle, never fewer than one; the phase line says when it
capped, and `--no-worker-cap` turns the cap off. The first line printed is the start time and how long the same
command took last time on this checkout (kept in `.test-runs.json`, not
tracked by git); the last line is how long it took now.

**Failures.** The runner keeps every failing test's full report (traceback,
worker, the machine's load average) in `.test-failures.jsonl`, and when the run
ends it runs each failed test again alone, on one process. The verdict is
"passes alone" (it failed only under load or beside other tests) or "FAILS
ALONE" (a real failure). A lane quoting a flake in a TODO item should paste
that line and the record from `.test-failures.jsonl` instead of "failed once,
passed alone". `--no-rerun` skips the second step.

A pull request needs a green `scripts/test.py changed` (below), finished with
`--resume` if the run was interrupted, and no more: the whole suite
(`scripts/test.py full`, `full --check`) is the end-of-version bug hunt
(TEST.261, SEC.33, GEN.301), not the gate of an ordinary merge, and not run
again on a merge candidate after main moved (run `changed --resume`, which
reruns only what main's change touched). The PR description says which run it used.

### Containers: the check before a run and the watchdog (TEST.376)

A lane machine is a container that restarts, and MySQL and Redis go down with it. A
database test whose server is down is skipped, not failed, and a run without
`PLANETGEN_TEST_MYSQL_USER` logs in as `planetgen` and fails with
`DatabaseLoginError`, so a broken container used to look like a green or a red
run. `scripts/test.py` (`src/tests/lane_services.py`) now looks first and watches
while it runs:

- **Before the run**, every command (`scripts/test.py env` does only this, and
  `--no-preflight` skips it): the `PLANETGEN_TEST_*` variables are loaded from
  `/etc/profile.d/planetgen-test.sh` when they are missing; MySQL and Redis are
  started when they do not answer (the commands `planetgen-env-up` runs, or
  `PLANETGEN_TEST_MYSQL_START` / `PLANETGEN_TEST_REDIS_START`); the MySQL binary
  logs are purged when under 3 GB is free; `argon2-cffi` must be installed and,
  for the browser tests, `playwright` must be 1.56.0 (the image's Chromium; do
  not run `playwright install`). What it cannot repair is printed as
  `PREFLIGHT FAILED: ...` and nothing runs (exit 2). The same check runs again
  at the start of each phase.
- **During the run**, a watchdog asks MySQL and Redis every 15 seconds. A server
  that is gone is restarted at once, the phase's pytest is stopped and started
  again without the tests that already passed (the progress log, as `--resume`
  does), up to 5 times a phase, and the closing line says
  `(N service restart(s): ...)`. A server that cannot be brought back ends the run
  with exit 3, `FAIL`, and the same command with `--resume` continues it once the
  services are up. A server that still has its process (`mysqld_safe` restarts a
  crashed `mysqld` itself, crash recovery answers late) is waited for, not
  started over, and a MySQL that only holds its port must miss three answers in a
  row before it counts as dead.

### Queue workers a test left behind (TEST.393)

A test that starts real queue workers (`workers.ensure_workers`) starts them in a session of
their own, so they outlive it; one that was stopped part way never tells them to stop, and each
keeps about 300 MB (fifty filled a 16 GB machine and MySQL was killed). Every test worker's
Redis is a database of its own, 1 to 15, and a real worker is started with `--url` naming it, so
`src/tests/stray_workers.py` can tell a test's worker from the machine's (database 0). The
conftest ends the workers on its own xdist worker's database when its session starts and when it
ends, and `scripts/test.py` ends every one still running when a run is over and says how many
and how old the oldest was. A test that starts real workers should still stop them; the count in
that line is the way to find the test that did not.

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
nothing, exits with 2 and says so. Resuming is never automatic; it is how an interrupted `changed` run is finished, and a PR that was
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
(which that run then writes). A green `changed` run is the merge gate (Boss, 2026-10-14, restated 2026-10-24: every
merge tests only what changed or would change); the whole suite
(`scripts/test.py full --check`) is the end-of-version bug hunt. A PR description says which one it used.
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

## The prevalence census (TEST.430)

A standing seeded report of how often the generator makes what it makes. Part A, the class census
(`src/tests/audit/prevalence.py`, `src/tests/test_prevalence_census.py`), draws planets and moons of six star types (M2V, K5V,
G2V, F5V, A0V and a K0III giant) in three zones (hot, ecosphere, cold) and checks, per cell, the share of each planet and moon
class, the share of each overall habitability colour (Blue, Green, Yellow, Red) and the share of bodies each life ceiling
allows (complex life, microbial only, none) against `src/tests/data/prevalence_expectations.json`.

- **Ranges and their kind.** Every range has a `min`, a `max`, a `kind` (`observed` rate, `norm` researched, or `stand-in` with
  its tier and tweak level per docs/design/research/stand-in-values.md) and a `source`. The class census has no observed class
  rates, so its ranges are stand-ins taken from the prevalence before the class cut-over: the reference commit, the draws per
  cell and the seeds are recorded in the file's `reference`.
- **Four flags.** A share outside its range (`out-of-range`), a class missing where expected (`missing`, a range with a `min`
  above zero and a share of zero), a class present where the class table bars it (`barred`), and a Blue or Green share outside
  its range (the colour measures). A planet draw that runs past 45 s is a fifth, `hang`.
- **Known anomalies.** A range the generator is known to break carries the item that will fix it (`anomaly`), and the test fails
  only on a flag with no such item. The report lists the known ones apart, and names a tag that no longer trips (stale: remove it).
- **The report.** `class-census-latest.md`, written to `$PLANETGEN_REPORT_DIR` or `.test-reports/`.
- **Re-baselining.** `python -m tests.audit.prevalence_baseline measure` draws the records at the code being measured, and
  `python -m tests.audit.prevalence_baseline build RECORDS --commit SHA --count N [--tag-current RECORDS --anomaly ITEM]` rewrites
  the file: the ranges are the reference's shares plus or minus four standard errors and 0.04, so they start wide and tighten as
  the classes settle. Never widen a range to hide a new fault.
- **When it runs.** It is a `slow` test that reads the generation code, so `scripts/test.py changed` selects it whenever
  generation changes. About 8 minutes on 4 cores.

Part B, the system census (TEST.431, `src/tests/audit/prevalence_systems.py`, `src/tests/test_prevalence_systems_census.py`), draws
100 full systems of each census star type and checks what only a whole system decides: planets, belts and comets per system, moons
per planet, the share with no planet, the share that are pairs, the share with a habitable-class world or a Blue or Green one, and the
share where life reaches each milestone (abiogenesis to a technological civilization), against
`src/tests/data/prevalence_systems_expectations.json`. Part E, planet occurrence (TEST.395,
`src/tests/audit/prevalence_occurrence.py`), reads the same systems and checks giants, hot Jupiters, cold giants, 1 to 20
Earth-mass planets inside 0.5 AU, Earth-mass planets in the ecosphere and the super-Earth to Neptune ratio per host group against
observed rates (`src/tests/data/prevalence_occurrence_expectations.json`, kind `observed`, recalled and unverified until the papers
are read). One live test draws the systems once for both; their reports are `system-census-latest.md` and
`occurrence-census-latest.md`. A system that raises is the `error` flag, as one over 180 s is `hang`.

Re-baselining B and E: `python -m tests.audit.prevalence_systems measure` and `build`, `python -m tests.audit.prevalence_occurrence`.

Part C, the star census (TEST.432, `src/tests/audit/prevalence_stars.py`, `src/tests/test_prevalence_stars_census.py`), draws 40,000
stars from the population model at the default galaxy's reference solar position and 2,500 planetless systems around each of M2V,
K5V, G2V, F5V and A0V. It checks the star type mix (M, K, G, F, A, B, O, white dwarfs, evolved stars) as shares of all stars against
the RECONS 10 parsec census and Mamajek's number densities, and the share of systems with a companion star (and, for G2V, with three
or more) against [R1175, R3646], and [R3108]
(`src/tests/data/prevalence_stars_expectations.json`, kind `observed` or `norm`). The file holds the norm range, never a generator
run: where two sources differ the widest range stands, and the sampling slack of the run (four binomial sigma at its own draw count)
is added when it is compared. A known fault carries its item as in the other parts. The report is `star-census-latest.md`; about 15
seconds on 4 cores. `python -m tests.audit.prevalence_stars` prints it without the test.

The star mass function census (TEST.454, `src/tests/audit/mass_curve.py`, `src/tests/test_mass_function_curve.py`; item 1 of the star mass
function chain, GEN.922) measures the same population bin by bin in mass against the published disk curves held as data in
`src/tests/data/mass_function_curve.json` ([R0705]'s log-normal fits for single stars and for systems, [R2498]'s broken power law,
all recalled, tier T2, until a lane reads the papers). It draws 60,000 systems at the reference solar position and counts two samples below
0.8 solar masses, where every star born is still alive: the primaries (held to the system function) and every star of the systems' member
trees (held to the single-star functions), and reports the generator's primaries above 0.8 solar masses by phase (main sequence, subgiant,
giant, supergiant, white dwarf) and the white dwarfs' mean mass (a norm). A bin's range is the lowest and highest share the sample's curves
give it, widened by the sample's relative tolerance and by four binomial sigma of the run. A bin a later item of the chain fixes carries
that item's number. The report is `mass-curve-latest.md`; about 20 seconds on 4 cores. `python -m tests.audit.mass_curve` prints it.

Part D, the phenomena census (TEST.433, `src/tests/audit/prevalence_phenomena.py`, `src/tests/test_prevalence_phenomena_census.py`),
first cut: the galaxy-wide totals of the classes the phenomenon scatter places (stellar-mass and intermediate-mass black holes, neutron
stars, magnetars, planetary nebulae, supernova remnants, and Wolf-Rayet stars), against `src/tests/data/prevalence_phenomena_expectations.json`. A total is the
scatter's own expectation (`shell_expected_by_class`) summed over every shell of the Milky Way-shaped reference galaxy with no mass cut, so
it is the same on every run and a range is the researched range itself (a catalogue count is a floor). The Wolf-Rayet total is the
population model's share of Wolf-Rayet stars among its living stars of 8 solar masses or more (120,000 draws) times the galaxy's stars
by population, with the range widened by four Poisson sigma of the draws. About a minute on 4 cores; `python -m tests.audit.prevalence_phenomena`
prints the report without the test. Giant molecular clouds, planetary nebulae and supernova remnant counts per star, and the globular
cluster count already have tests of their own (`test_nebula_rates.py`, `test_cluster_count_law.py`); the second cut adds seven measures that are not scatter totals: the rates the sectors' fill draws from, read where the fill reads them
(`tuning.phenomenon_rate_per_star`, `ROGUE_PLANET_MASS_BINS` and `PHENOMENON_RATE_SCALE`, so a dialled rate shows): terrestrial rogue planets
per star and per pc^3, Jupiter-mass rogue planets per star, brown dwarfs per pc^3, runaway stars per star and interstellar debris per pc^3, each
against its researched range; and the open clusters a galaxy of the reference shape draws (`open_clusters.generate_open_clusters`, about 17,500,
against the 1,867 Gaia clusters as a floor and the generator's cap of 100,000). These are rates, not realized sector draws. The nebula classes
(emission and reflection nebulae around hot stars) are the remainder of TEST.433.

Part F, the body sweeps (`src/tests/audit/sweep.py`): one shared set of generated bodies that the consistency sweeps read. It builds a
planet of every class and zone pair the class table allows round an M2V, a G2V and an F5V host, with the moons each drew, and adds the
planets and moons of 24 seeded generated systems (the history the class-by-class builds skip); the draws are seeded from the host, so a
repeat gives the same bodies. A rule is a `Rule`: `check(entry)` lists the violations, `breaker(entry)` returns a copy of the body with
exactly that quantity broken, and `prove_rules_can_fail` checks the sweep reports every broken copy, so a rule that cannot fail does not
pass. The first sweep is the atmosphere and temperature (TEST.209, `tests/audit/sweep_atmosphere.py`, `tests/test_sweep_atmosphere.py`):
11 rules (the partial pressures add up; the mean molar mass, the air density and the scale height follow from them; the temperature is the
Stefan-Boltzmann value on the stored albedo, greenhouse quantity and the starlight at the orbit; the ozone follows the oxygen and the star;
the aerosol depths follow their rules; every gas in the text is in the table). The second is the bulk composition (TEST.210,
`tests/audit/sweep_bulk.py`, `tests/test_sweep_bulk.py`): 8 rules (the five shares make one; the details are shares and the envelope
details exist exactly where there is an envelope; the density and gravity follow from the mass and radius; the radius is the structure
model's, to a millionth for a planet and 4 percent for a moon, with bodies that have an envelope left out; the body is at least as dense
as its materials at zero pressure; the phrase names exactly the materials the shares hold). About a minute each on 4 cores. A sweep
that ever fails on a real fault is filed as a bug and tagged in the test with its item until it is fixed.

The third is the possibility tables (TEST.214, `tests/audit/sweep_tables.py`, `tests/test_sweep_tables.py`). It walks every table of
`possibility_table.REGISTRY` (it imports the modules that use the loader, so a table added later is covered without a change to the test)
with 60 seeded subjects for each class a table names, plus two that satisfy each row's own conditions, and checks that: every row has a
source and a state (no licence assertion; licences are in `DATA-LICENSE.md`); every named class reaches a row; the gated weights of every
subject add up to 1 within 1e-12; no subject fails for a reason but an empty gate; each row's drawn fractions add up to 1 and name only
components of the shared materials table; `Table.draw` picks each row as often as its probability says (five binomial sigma); every row is
given a probability for every class it names; the open rows are the ones listed with a note in `tests/data/possibility_sweep_notes.json`
(a new open row needs a note there); and every `Table` that a fresh process holds after the imports is the registered one. The same file
lists, each with a note, the tables whose rows are independent flags and not exclusive odds (an all-zero answer is allowed there) and the
rows that no subject reaches today. Subjects take their numbers from physical ranges and agree with one another where the generator's do
(the escape speed follows the mass, the equilibrium temperature the insolation).

The fourth is the classes from conditions (TEST.301, `tests/audit/sweep_classes.py`, `tests/test_sweep_classes.py`). On the same bodies it
checks that every planet and moon has a class letter of the rule file and of its own kind, that `classify` is total and pure (the same
twice and from the stored columns, and it leaves the body alone), and that the stored letter is the one the conditions give, re-derived from
the rule rows. It also generates 12 seeded systems with `+habitable_world` and counts the requests answered `aim` or `fallback`. Two
results are known anomalies listed in `tests/data/class_sweep_expectations.json`, so they are reported in `class-sweep-latest.md` and
`habitable-request-latest.md` and do not fail: the stored letter is still drawn before the body is built (GEN.646), and the air solve of 4c
leaves no free oxygen so most requests answer unmet (GEN.855). When an anomaly's tag is removed the rule fails on any miss, and a tag on a
rule that holds everywhere is flagged stale in the report. The reclassification after an edit is left for 8.6.


Part G, the life census (TEST.435, `src/tests/audit/prevalence_life.py`, `src/tests/test_prevalence_life_census.py`): 240 seeded whole
systems, every terrestrial planet and moon with a life chemistry put on the highest rung it reached (origin, producers, complex cells,
multicellular, animals, sentient), and the shares of life worlds checked against Boss's targets for the life class system
(`docs/design/habitability/life-classes.md` 5.2: the majority stop below animals, 15 to 35 percent reach animals, at most 2 percent are
sentient and fewer than one in twenty animal worlds) in `src/tests/data/prevalence_life_expectations.json`, widened by four binomial sigma of
the run's life worlds. The shadow life class of each life world (`life_classes.best`, GEN.886) is counted in the report and not checked, as the
design says per-class bands are measured and not forced; the class store (GEN.889) will let the census read the stored class. The report is
`life-census-latest.md`; about a minute on 4 cores; `python -m tests.audit.prevalence_life [--host-systems N]` prints it without the test (the live test draws 300 host systems each; with only about 30 oxygenic worlds in a run the four-sigma slack makes the animal band wide, so a baseline run takes several times more). The report also lists the highest rung reached by each stored life class (epoch 120). Known anomaly: GEN.855
(no biogenic oxygen, so no world reaches animals).

Part H, the class-blind life test (TEST.452, `src/tests/audit/sweep_life.py`, `src/tests/test_sweep_life.py`; Boss 2:46 PM EDT 10-29, GEN.915):
life depends only on the conditions, never on the class letter. Every terrestrial body of the shared sweep is given each of the 20 class
letters with its stored conditions unchanged, and the life model (`life_classes.score_classes`, its best class) and the habitability index
must answer the same on every letter; a failure names the body, the letter and the life class or index field that moved. The life code is
scanned (AST) for every read of the class letter or the class tables, and each one must be listed with its reason in
`src/tests/data/class_blind_life_notes.json` (today the class description text, the planet or orbital flavour line, and the value copy that
the index does not read); an unlisted read fails, a listed place that no longer reads it is reported stale. Both checks are proved able to
fail by breaking the life model, the index and a planted module. The report is `class-blind-life-latest.md` (with the life-capable bodies of
each letter). The draws that run before the life model and key on the letter (water fraction, albedo range, porosity, mass-radius branch)
change the conditions and not the way they are judged; they are GEN.915 and GEN.646. The second layer (epoch 120) is on the generator: every planet and moon of the shared sweep (M2V, G2V, F5V; the class is derived from the
conditions since epoch 121, so it cannot be forced on a bare `Planet`) goes through the life step (`generation.life.apply_life_data`) on a copy; a
planet the life model supports must be given life and one it does not support must not, the class drawn must be a supported one, and the same
conditions put under every other class letter must give the same class, stage and schedule on the same stream. A failure names the letter and the
class. The report `class-blind-life-draws-latest.md` lists, for each derived letter, the planets, the supported ones and the ones given life.

## Setting up a cloud session

A fresh container needs `pip install -e ".[test,api,browser]"`, the NLTK
corpus (`python -m nltk.downloader words`), MySQL with a user that
has `ALL ... WITH GRANT OPTION` (one test creates a restricted user), Redis
(`redis-server --daemonize yes`, then `export PLANETGEN_TEST_REDIS_URL=redis://127.0.0.1:6379/0`; the tests that
need it skip without that variable), and a Playwright whose Chromium matches the one on the machine: the image
carries Chromium 1194, so use `pip install playwright==1.56.0` (a newer
Playwright asks for a build the image does not have and every browser test
skips).

To have all of that ready in every new session, paste
[`examples/cloud-environment/setup-planetgen-env.sh`](../../examples/cloud-environment/setup-planetgen-env.sh) into the
cloud environment's Setup script field (it installs MySQL 8, Apache 2, Redis and every dependency, and leaves a
`planetgen-env-up` command that starts the services each session); its
[README](../../examples/cloud-environment/README.md) says what each part does.

## Database engines

Tests run on MySQL 8.0.19 or later; MariaDB is not supported and a run against it fails by design. CI's `test` job
has four legs, each with its own Redis container: the whole suite on Python 3.12 against MySQL 8.0 and MySQL 8.4,
then the `generation` tests alone with 2 workers (MySQL 8.4, Python 3.12) and with 4 workers (MySQL 8.0,
Python 3.13).

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
  - `test_fuzz_text_rendering.py`: the Document renderers (text, Markdown, HTML).
  - `test_fuzz_utils_and_config.py`: physics helpers, formatting, config
    loading, logging and progress files.
  - `test_fuzz_cli.py`: every `planetgen` command and option.
  - `test_fuzz_web_routes.py` and `test_fuzz_api_auth.py`: every web page
    and API endpoint, logins, CSRF and API keys.

- **Guard tests** (`test_no_*`): read the source and fail when a rule is
  broken, with a named list of exceptions each carrying its reason.
  `test_no_own_log_formatting.py` (TEST.183) fails when code under
  `planetgen/` writes or logs a line built with its own timestamp or
  divider instead of the logging layer; an exception goes in its
  `ALLOWED` list with the reason it is not a log line (the form and the
  rules are in [`logging.md`](logging.md)).

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
  The representation audit (`tests/audit/representation.py`, TEST.160) is the first audit on it: eight region classes of
  the Milky Way's shape (the bulge, which holds the bar; an arm and the space between arms; the thick disk; the halo; a
  Magellanic-type satellite; a globular cluster; the intergalactic floor), filled with every per-star rate dialed up to
  0.3 so that a kind the place's gate forbids has a chance to appear and the disk, where it is allowed, shows it does.
  A fill holds a handful of systems, so the kinds that need a free spot of their own (a neutron star, a black hole, a
  brown dwarf, a planetary nebula's central star) are audited from above only. The whole-sphere tests of TEST.184 and
  TEST.185 (`test_census_nebula_wholesphere.py`, `test_census_component_wholesphere.py`) read the same regions.
  `python -m tests.audit.representation [report.txt] [--sectors N]` prints the report and exits 1 when a count breaks a
  rule.

A known bug that is not fixed yet is kept as a test marked
`xfail(strict=True)` whose reason names the bug and a repro. When someone
fixes it, the test starts passing, strict mode fails the run, and the
marker has to come off: the fix can't land without its regression test.

## Opt-in tests

These skip unless you turn them on:

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
- The `globaldocs` tests (DOC.69) compare a generated global document with its source: the architecture map's module list
  (`test_architecture_map.py`), the docs index (`test_docs_index.py`) and the plan pages (`test_todo_docs.py`). Only the TODO
  thread regenerates those files (`python scripts/build_global_docs.py`, after a batch of merges; `--check` lists the stale
  ones), so a lane's pull request never carries one and a stale one never fails a lane's gate: the tests skip unless
  `PLANETGEN_GLOBAL_DOCS=1`, which `scripts/test.py <command> --global-docs` sets. The release stamp runs
  `scripts/test.py full --check --global-docs`.
- `test_golden_seed.py::test_the_pinned_galaxy_is_made_at_one_worker_and_at_four` builds the pinned galaxy at four workers and
  compares it with the one-worker galaxy (TEST.407). The ordinary run builds the galaxy once, at one worker, and
  compares it with `galaxy.json`, which is also all that `python -m tests.golden_update` builds; the second build,
  about as long again, is the end-of-version run's: `scripts/test.py full` (and `full --check`) set
  `PLANETGEN_VERSION_GATE=1`, or set it by hand to run the comparison alone.
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

`.github/workflows/ci.yml` (the math check, the four test legs described under "Database engines", browser and
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
