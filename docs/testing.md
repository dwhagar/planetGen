# Testing

The suite lives in `src/tests/` and runs with plain `pytest` from the
repository root (`pytest.ini` puts `.`, `src` and `src/html` on the path).

```sh
pip install -e ".[test,api]"
python -m nltk.downloader words
pytest
```

## The database tests

Tests that save and load systems, sectors or users need a MySQL or MariaDB
server. They skip (not fail) without one. Point them at a server with:

```sh
export PLANETGEN_TEST_MYSQL_HOST=127.0.0.1
export PLANETGEN_TEST_MYSQL_PORT=3306
export PLANETGEN_TEST_MYSQL_USER=root
export PLANETGEN_TEST_MYSQL_PASSWORD=...
```

The user needs to create and drop databases: each test gets its own
uniquely named, throwaway database. Run the suite serially (no `pytest -n`)
when these are set; a few API tests share server-wide state and can trip
over each other under parallel workers.

## Kinds of test

- **Unit and behaviour tests** (`test_<area>.py`): one module or page each,
  with fixed inputs and exact expectations.
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

A known bug that is not fixed yet is kept as a test marked
`xfail(strict=True)` whose reason names the bug and a repro. When someone
fixes it, the test starts passing, strict mode fails the run, and the
marker has to come off: the fix can't land without its regression test.

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
