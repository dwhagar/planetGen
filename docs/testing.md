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
- generation runs in-process (`PLANETGEN_WORKERS=1`) unless you set
  `PLANETGEN_WORKERS` yourself.

## Generation at several worker counts

A real run hands its sectors and bright-star layers to worker processes.
CI runs the tests marked `generation` (the files that run `generate.py` or
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

The user needs to create and drop databases: each test gets its own
uniquely named, throwaway database, so the database tests are safe under
`pytest -n auto` (see "Running in parallel" above).

## The math check runs first

`src/planetgen/physics/mathcheck.py` is a fixed list of checks that the
generator's math gives the right answers (TEST.63): known values from real
astronomy (the Sun, Earth's and Jupiter's orbits, the habitable zone, white
dwarf sizes, Holman & Wiegert's stability limits, the Kepler and Barker
equations), identities that hold for any input (unit round trips, the
constants agreeing with each other, energy conservation around a Kepler
orbit, the sector grid tiling each ring), and seeded draws from every
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
  `generate.py check-math` by hand, and every bulk path (see
  [cli.md](cli.md#subcommands)); `update.sh`/`update.ps1` warn.

```sh
cd src && python -m planetgen.physics.mathcheck -v   # the report, exit 1 on failure
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

CI can run on your own computers. Which jobs run where, what each machine
needs, the security settings and troubleshooting are in
[`ci-runners.md`](ci-runners.md).
