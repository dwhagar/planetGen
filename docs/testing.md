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

- each test that uses `mysql_config` gets its own uniquely named database;
- each worker process gets its own control database
  (`PLANETGEN_CONTROL_DATABASE` is set to `pgtest_control_<worker>_<random>`
  by `conftest.py`, overriding any value in your environment, and dropped
  at the end of the run), so nothing in the suite touches
  `planetgen_control` or a real control schema;
- generation runs in-process (`PLANETGEN_WORKERS=1`), so a test's
  patched functions are the ones that run.

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

The user needs to create and drop databases: each test gets its own
uniquely named, throwaway database, so the database tests are safe under
`pytest -n auto` (see "Running in parallel" above).

## Picking a slice of the suite

`conftest.py` marks every test (markers registered in `pytest.ini`):

- `db`: uses a real database (the `mysql_config` fixture, directly or
  through another fixture);
- `slow`: the brute-force and seeded-sweep files (`test_fuzz_*`,
  `test_bughunt_*`);
- `browser`: the headless-browser checks (`test_web_a11y.py`).

```sh
pytest -n auto -m "not db and not slow"   # quick loop: about a minute on 4 cores
pytest -n auto -m db                      # just the database tests
```

A pull request still needs the whole suite green.

## Database engines

Local runs here use MariaDB 10.11. CI's `test` job runs the whole suite
three times: Python 3.9 on MySQL 8.0, and Python 3.12 on MySQL 8.4 and on
MariaDB 11.4. MySQL 8 is stricter about reserved words and `GROUP BY`
than MariaDB, so a query that works locally can still fail there.

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
  - `test_fuzz_galaxy_grid.py`: the ring/layer/slot sector grid,
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
  - `test_fuzz_cli.py`: every `generate.py` command and option.
  - `test_fuzz_web_routes.py` and `test_fuzz_api_auth.py`: every web page
    and API endpoint, logins, CSRF and API keys.

- **Browser checks** (`test_web_a11y.py`): every page in headless
  Chromium at phone and desktop widths, in light and dark, against
  axe-core and the layout rules described in
  [`html-interface.md`](html-interface.md#flask-pages) ("Browser checks").
  They need `pip install -e ".[browser]"` and `python -m playwright install
  chromium` plus the MySQL test server, and skip without them. CI runs
  them in their own `browser-a11y` job.

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

## The Deep fuzz workflow

`.github/workflows/deep-fuzz.yml` runs every fuzz file under the `deep`
profile each Monday, and can be started by hand from the Actions tab
(**Deep fuzz** > **Run workflow**, optionally with a number of examples).
A red deep run is a real bug report: its log carries the shrunk input.

## CI runners

Every workflow job picks its runner from a repository variable (Settings >
Secrets and variables > Actions > Variables). Each holds JSON: a label
list such as `["self-hosted", "Linux", "X64"]`, or a hosted runner name in
quotes such as `"ubuntu-latest"`. Unset, a job uses the default below.

| Variable | Jobs | Default |
|---|---|---|
| `RUNNER_LINUX` | tests, browser checks, dependency audit, deep fuzz, release note, version stamp | `["self-hosted", "Linux"]` |
| `RUNNER_WINDOWS` | Generate page jobs on Windows | `["self-hosted", "Windows"]` |
| `RUNNER_MACOS_INSTALLERS` | `install.sh` on macOS | `"macos-latest"` (GitHub-hosted) |
| `RUNNER_WINDOWS_INSTALLERS` | `install.ps1` on Windows | `"windows-latest"` (GitHub-hosted) |

The installer jobs stay on GitHub's throwaway machines by default because
they install planetGen as a system service (launchd daemons, scheduled
tasks, a server on port 8000) with sudo or admin rights, which would stay
behind on a machine of your own.

A self-hosted Linux runner needs Docker (the MySQL 8 service container) and
passwordless `sudo` for `playwright install --with-deps`. The MySQL
container gets a free host port, so several jobs (or a MySQL of your own
on 3306) can share a machine. Each job installs into its own virtualenv
under `RUNNER_TEMP`, so packages don't carry over between jobs.
`actions/setup-python` downloads the Python versions the jobs ask for
(3.9 and 3.12) into the runner's tool cache on first use.
