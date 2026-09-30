### Added
- **Brute-force (property-based) tests across the whole codebase.** New
  `test_fuzz_*.py` files use Hypothesis to hunt for inputs that break the
  galaxy grid, system and planet generation, sector placement, names, text
  rendering, the physics helpers, config loading, every `generate.py`
  command, and every web page and API endpoint. New `test_edge_admin_scripts.py`
  covers the job runner, `resetDb`, `migrateDb` and `updateOrbits` against
  empty, broken and unreachable databases. Every normal `pytest` run uses
  a fixed, repeatable `ci` profile; a weekly **Deep fuzz** workflow runs
  3000 random examples per test. See `docs/testing.md`.

### Fixed
- **Bugs the new tests found.** Among them: a point a sector cell sampled
  could count as outside that cell; negative ring or out-of-range slot
  designations were accepted; the galaxy density could exceed 1; zero,
  NaN or infinite map radii picked the wrong tile level; a planet pinned
  by mass alone crashed; star types with trailing junk (`G2Vjunk`) were
  accepted; NaN or infinite values in the physics helpers came back as
  NaN instead of an error; NaN sector edges and positions poisoned
  placement; a sector's grid cell was lost on save and load; `generate.py`
  accepted NaN or infinite numbers, bad star types, broken
  `--system-file`s and out-of-range ports, and printed tracebacks for
  database and file errors; `generate.py plan` accepted shapes that could
  not be built; a huge `--radius-pc` tried to enumerate far more sectors
  than the galaxy has; decoration words and offensive words hidden by an
  apostrophe or space slipped into generated names; and the debug log did
  not redact every password or token option.
- **Web and API inputs that crashed or misbehaved.** A non-ASCII CSRF
  token, a query string that is not valid UTF-8, a page number past
  2^63, NaN or infinite search bounds and radii, a Unicode digit in a
  sector address, a huge galaxy cell coordinate, a non-JSON request body,
  non-string logins or key labels, and over-long names, URLs, usernames
  and key labels now get a clear 400 instead of a 500 or a traceback.
  `/api/health?db=` with an unknown database answers 404, and
  `/api/databases` lists a schema it can't open instead of failing.
- **A flaky job runner test.** It checked for the job lock in the instant
  between the runner writing its final state and releasing the lock.
