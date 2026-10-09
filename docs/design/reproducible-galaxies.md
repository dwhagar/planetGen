# Reproducible galaxies

How a galaxy will be rebuilt from a 128-bit seed and a version key, how
admin changes and the passage of time are kept on top of that, and how a
damaged database is found and repaired. Recorded 2026-10-02 from Boss's
decisions of 01:34Z to 02:31Z and the seed math report (dependency tree
thread). This is the design note OPS.11 asks for.

**Scope change (Boss, 2026-10-09 20:42Z).** "Remove the idea of us letting
the USER regenerate an entire galaxy from the seed and version, we'll use
that internally, but no need for it to go anywhere else." What stays, as
internal mechanism only: the 128-bit galaxy seed, per-unit seeds
(GEN.39, GEN.56), GEN.57 (a sector's contents depend only on the seed, the
version and its address, so parallel workers, fills, backfills, the
phenomenon scatter and its sector fill, and settle agree), the version
stamp that rejects stale data (DB.7, OPS.13, OPS.14, OPS.28), the
fingerprint and golden-seed tests (GEN.58, GEN.135, GEN.136, TEST.77), and
the seed-based sector repair (DB.9, DB.17). What is dropped: `generate.py
reproduce` (OPS.12), the pages and routes that show the seed and version
(ADM.17, API.16), the net-difference JSON of admin changes (GEN.59), the
daily merge, the 18 backup slots and their dashboard (GEN.61, OPS.18,
ADM.19, ADM.20), and repair from the JSON deltas (DB.10). OPS.16 now runs
the positional update only. Sections below that describe the dropped
parts are kept as history and are not planned.

**Status (2026-10-09, checked against main at 27220d8): partly built.**
GEN.39 (per-unit seeds), GEN.56 (every draw through `util/draw.py`), DB.6
(seed, version key, run history) and OPS.10 (the log line) are in, each
described under "As built" in section 3 or 5. Everything else here is
planned: the fingerprint, the golden test, the update history, the
settings file, the maintenance run, the check and repair, and `reproduce`.
Each piece names its TODO item and phase; `docs/TODO.md` holds each item's
full text and `docs/plan/` the order.

Measurements of what breaks reproducibility across Python versions, operating
systems and CPUs, and the design built on them (deterministic-math helpers,
fingerprint encoder, generator epoch, net-diff delta format, golden-test
structure), are in [generation-determinism.md](generation-determinism.md)
(research, 2026-10-09). Where that note recommends something different from
the TODO text, this note says so in the section concerned.

The check, the parity codec and file, the stale-record rules, how MariaDB
behaves on a damaged page, and the long-migration helpers are in
[db-check-and-parity-repair.md](db-check-and-parity-repair.md) (research,
2026-10-09); section 9 below summarises and corrects the TODO text.

## 1. Decisions already taken

Boss's decisions, as recorded; sections 2 onward say which parts are built and
which are still recommendations.

> Ok use a 128 bit value and store the seed in the database, and put it
> in the log at the top of any generation, also populate the TODO upward
> from here to eventually build a system that a version number and a
> seed value would reproduce the same galaxy by the end of the phases.
> (01:40Z)

Later the same night: the version packed in hex (01:46Z, 01:53Z); two
more build digits and the Python version, OS and architecture in the
key (02:08Z); each update records the key, keeping the last 10 (02:08Z);
a consistency check, parity repair and a settings JSON file (02:13Z);
admin changes as a net difference, regenerate seeds and the word list in
the file (02:20Z); the file name (02:21Z); the daily maintenance run and
18 backups (02:28Z). GEN.39 starts from a wiped galaxy (01:46Z: "I'm
going to nuke the galaxy anyway"), so nothing here migrates or supports
an unseeded galaxy.

Later decisions that touch this design: every generation draw reaches the
unit's stream through a `ContextVar` instead of an `rng` argument on every
function (2026-10-09 07:02Z); IDs and the naming key name the interstellar
objects while stars and sectors keep word-salad names, and planets, moons
and belts keep the "<system> I" pattern (2026-10-08 03:57Z and 04:00Z,
[object-ids.md](object-ids.md)).

## 2. What "the same galaxy" means

Every generated object has the same address, position, properties and
name. Database ids, timestamps and population data rebuilt later are not
compared. Nearest-system links are rebuilt from positions, so they are
left out too. Admin changes are compared after replaying the settings
file (section 7).

A seed reproduces a galaxy only on the exact setup that made it: the
same PlanetGen release, Python version, OS, architecture and word list.
New code turns the same seed into a different galaxy and no formula
converts between versions, so an older galaxy is rebuilt by running its
own release. The guarantee is exact equality on the same setup; across
Python versions or platforms it does not hold today (Python 3.9 and 3.10 and
later give different sectors from one seed, because `math.hypot` changed;
[generation-determinism.md](generation-determinism.md) section 2), and
the recommended fix makes integers, strings and structure exact and floats
equal to 9 significant digits there. A seed cannot be worked out for a
galaxy that already exists: one made before GEN.39 came from true
randomness, and finding a seed for a given galaxy means running SHA-256
backwards (about 2^128 tries).

## 3. The seed (GEN.39, DB.6, phase 0)

- One 128-bit galaxy seed, stored in `galaxy_shape` as `BINARY(16)`,
  shown as 32 uppercase hex digits, written once when the galaxy is
  first planned. `planetgen plan --seed <32 hex>` sets it; without it
  one is drawn at random.
- Each unit of work gets its own seed:

  ```
  unit seed = SHA-256(galaxy seed || "kind:" address)
  ```

  for example `sector:12/3/0`. Units are a sector fill, a bright-star
  layer, a backfill block, a phenomenon and a population pass. The full
  256-bit digest seeds a `random.Random` that is passed down
  explicitly, so a sector's numbers depend only on the galaxy seed and
  its address, never on the order sectors run in or the worker count.
- The version is not mixed into the hash: a release that changes one
  formula changes only what that formula touches.
- The bright-star scatter's 63-bit seeds come from the unit's stream
  (`draw.getrandbits(63)`), so no step throws bits away.
- Odds, from the seed math report: two of a billion random 128-bit seeds
  match with a chance of about 1.5 x 10^-21; two of 12 billion sectors
  sharing a unit seed, about 10^-57.

GEN.39 builds after PERF.21 in the parallel path thread. Its tests: one
seed gives the same sectors at 1, 2 and N workers, and two seeds that
differ only in their high 64 bits give different output.

As built (GEN.39): `planetgen/galaxy/seed.py` holds the seed helpers.
A sector fill (its save included) and the bright-star scatter, bands
and backfill blocks run on their unit seed; the sector seeds the
module-level `random` stream for the length of the unit
(`galaxySeed.seeded`, which puts the stream back afterwards) rather than
passing a `random.Random` down, which is GEN.56's change. `secrets` and
`SystemRandom` are gone from generation (`utils.reseed_rng` is removed,
`spaceSector._rng` reads the module stream). A save retried after a
deadlock replays the same draws. The `phenomenon` and `system` subcommands and the run's own
choices (a random start's address) still draw from the run's stream.

As built (GEN.56): every generation draw goes through
`planetgen/util/draw.py`. A `draw.Stream` keeps a `random.Random` for its
`random()` alone and builds `uniform`, `randint`, `randrange`, `choice`,
`choices`, `shuffle`, `sample`, `gauss` and `getrandbits` on it, so the
numbers don't change with Python's own recipes. `seeded` binds the unit's
stream for the length of the unit with a `ContextVar` (Boss chose this
over an `rng` argument on every function, 2026-10-09 07:02Z): generators
call `draw.uniform(...)` and get the running unit's numbers, each thread
its own, with no shared stream and no lock. Draws outside a unit come
from the run stream (`draw.set_run_seed`, from `planetgen`'s run seed or a
work queue task's seed). A pre-placed bright star's system binds its
stored seed; a regenerated sector binds a fresh one. The debug log's
per-draw tracing wraps `draw`'s functions. `test_reproducible_draws.py`
scans the generation packages (`galaxy`, `generation`, `names`,
`physics`, `population`, `util`) for the stdlib `random`, `secrets`,
`uuid`, `os.urandom`, `SystemRandom`, `numpy.random` and clock seeds, pins
the helpers' recipes, and generates one sector under two hash seeds and
two locales. Population streams still key on row ids (GEN.57).

As built (DB.6): `planetgen/galaxy/version_key.py` computes the key
(`version_key`) and the versions stored beside it (`current`). The galaxy
records them in `galaxy_shape` whenever its seed is written (schema
v52), and every `planetgen` run that changes the galaxy writes a
`generation_runs` row with its command line, its own run seed, the
galaxy seed, the key and its outcome.

## 4. The version key (DB.6, phase 0)

22 uppercase hex digits, no separators:

| Part | Digits | Values |
|---|---|---|
| PlanetGen MAJOR, REVISION, BUILD | 4 + 4 + 6 | MAJOR.REVISION.BUILD, each packed (a plain sum collides: 7.127.352 and 7.128.351 both sum to 486) |
| Python major, minor, micro | 2 + 2 + 2 | One digit each would overflow at 3.16, and micro releases pass 15 |
| OS | 1 | 0 Linux, 1 Windows, 2 other Unix and BSD, 3 macOS, F unknown |
| Architecture (`platform.machine()`) | 1 | 0 x86-64, 1 ARM64 (Apple Silicon included), 2 32-bit x86, 3 32-bit ARM, 4 RISC-V 64, F unknown |

PlanetGen 7.127.352 on Python 3.12.3, Linux, x86-64 is
`0007007F000160030C0300`. One helper computes it from the running code,
so the log line, the update history and the mismatch warning agree. It
is stored with the galaxy (DB.6) and with each sector (DB.7, phase 1),
next to the full version string, Python version and platform.

Inputs outside the key that can change a galaxy are hashed (SHA-256) and
kept with each update's history row (section 6): nltk's `words` corpus
files (downloaded apart from the pinned packages), `offensive_words.txt`
and any name lists, and `requirements.lock`. The word list itself goes in
the settings file (section 7).

The OS, architecture and Python digits do not cover everything that moves
a float. Generation also uses numpy and scipy (`galaxy/nebula_shape.py`,
`physics/sector_path.py`, `physics/kepler.py`), scikit-image (marching cubes)
and astropy (`physics/constants.py`), the lock pins a different version of
each per Python leg, and the C maths library and BLAS kernel depend on the
OS and CPU in ways the digits do not name. The lock hash covers the
versions; [generation-determinism.md](generation-determinism.md) section
2.3 lists the specific fixes (freeze astropy constants, no BLAS in
generation).

**The key follows releases, not output (recommendation).** MAJOR.REVISION.BUILD
moves on nearly every release (each patch or minor note bumps REVISION and
BUILD is a sum of TODO counters), so the key alone would flag almost every
update as a change. Beside the key, the history row, the settings file and
the log line carry a `generator_epoch` (an integer bumped only by a PR that
changes the golden digest, enforced by `bump_version.py --check-pr`) and a
battery digest (a hash of a canned set of tiny generations run in the
running environment, OPS.15). A galaxy is reproducible here when its epoch
matches and the battery digest now equals the stored one. DB.7 then stores
`generator_epoch` and a run id per sector instead of four text columns, the
full key text living on the `generation_runs` row. Detail:
[generation-determinism.md](generation-determinism.md) section 5.

## 5. Every draw seeded, every sector independent (phase 1)

- **The log line (OPS.10, phase 0).** Every `planetgen` subcommand,
  web or API job and work queue run first writes, at normal level, the
  galaxy seed, the version with its key, and the run's command, for
  example `Galaxy seed 3f2a...c901, PlanetGen 7.127.352
  (0007007F000160030C0300), run: sector 12 3 0`.
  As built (OPS.10): `versionKey.run_line` writes it, first thing after
  `planetgen` sets up logging (a `plan --seed` names the seed given,
  otherwise the stored one, `none yet` before the first plan), and
  `planetgen.web.job_runner` writes it at the top of a job's `output.log` with the
  job's title. A work queue always runs inside a `planetgen` run, so
  its line is that run's.
- **The run history (DB.6, phase 0).** A `generation_runs` table, one
  row per run that changes the galaxy (command and options, version,
  start and end, outcome). A galaxy is built by a series of commands
  (plan, sectors, scatter, backfill), not by the seed alone, so this is
  what makes a rebuild possible.
- **Every draw from the derived seeds (GEN.56, built).** Before GEN.56,
  `secrets.choice` and `secrets.randbelow` picked planet classes, life,
  moons and name syllables, `spaceSector._rng` was a `secrets.SystemRandom`
  and `utils.reseed_rng()` reseeded the global `random` from `secrets` in 13
  generator modules; all of that now draws from the unit's stream (the
  "As built (GEN.56)" paragraph in section 3). Python promises the same
  sequence only for `random()` itself, so generators draw through the
  in-house helpers built on `random()` alone. Measured on Python 3.9 to 3.14,
  `draw.Stream` gives identical output (stdlib `choice`, `shuffle` and the
  others also agreed there, so the rule is insurance, not a response to
  observed breakage; only `paretovariate` changed), with one exception:
  `Stream.gauss` uses `log` and `cos` from the C maths library
  ([generation-determinism.md](generation-determinism.md) section 2).
  Nothing may depend on
  set or string-keyed dict iteration order (`PYTHONHASHSEED`) or locale
  sorting. A test scans the generation modules for `secrets`,
  `SystemRandom`, `os.urandom`, `uuid4`, `time`, unseeded `random.seed()`
  and direct calls to the unsafe helpers (with an allowlist for login,
  CSRF, API keys and job ids), and runs a small generation under two
  hash seeds and two locales.
- **Order independence (GEN.57).** Names, seeds and skips key on
  addresses and the galaxy seed, not on database ids or arrival order:
  in a name collision the sector with the lower address keeps the name
  and the other is renamed from its own seeded stream (with GEN.46);
  population seeds stop keying on row ids; the backfill stops depending
  on which sectors are already filled (with GEN.44).
- **Deterministic-math helpers (proposed item, no ID yet).** Measured:
  `math.hypot` and `math.dist` give different results on Python 3.9 and
  3.10 and later, and `sum()` of three or more floats differs between 3.11
  and 3.12, so one seed gives different sectors on the 3.9 CI leg and the
  others. Generation code uses `planetgen/util/detmath.py` (`hypot`, `dist`,
  `fsum`) instead, `test_reproducible_draws.py`'s scan flags the builtins,
  astropy-derived constants become literals, and `@` (BLAS) leaves
  `nebula_shape.py` and `sector_path.py`. The change alters stored values, so it
  bumps the epoch. It must land before TEST.77 or that test fails on 3.9.
  Rules and measurements: [generation-determinism.md](generation-determinism.md)
  sections 2 and 3.
- **Fingerprint (GEN.58).** `planetgen fingerprint` prints a canonical
  SHA-256 per sector and for a region, in address order with canonical
  number formatting, skipping ids and timestamps, either as first
  generated or with the settings file applied.
  As built (GEN.58, `planetgen/db/fingerprint.py`): each row is one
  line of canonical JSON (floats in shortest round-trip form, bytes in
  hex), foreign keys replaced by what they point at (a sector's address,
  an object's `uid`, a system configuration's digest); a sector's digest
  is the SHA-256 of its lines sorted, so save order doesn't count; the
  region's is over its sectors' lines in address order, with the plan's
  first for the whole galaxy. Clocks (GEN.106), nearest systems, the
  location text, sector paths and stats, and population data are left
  out. It reads the galaxy as stored; the settings-file layer waits for
  GEN.59.
  Not built, recommended in
  [generation-determinism.md](generation-determinism.md) section 4: round
  floats to 9 significant digits before hashing, so 1-ulp differences in
  `exp`, `sin` and `atan2` between machines do not show; add separate
  structural and float digests over named sections to localise a failure;
  save a leaf digest per sector with a (ring, layer) node table and ring
  and galaxy roots (DB.9). The built encoder hashes floats exactly, so run
  TEST.77 on every CI leg before deciding whether rounding is needed.
- **Golden test (TEST.77).** A fixed seed builds a small galaxy (plan, a
  few sectors, a scatter, a backfill) at 1 and 4 workers on every CI
  Python leg, and its fingerprint must match the one pinned in
  `src/tests/golden/`. A PR that changes generation output updates the
  pinned fingerprint and says so in its `changes/` note;
  `bump_version.py --check` checks the two go together. It also watches
  for the one known risk: a maths function (`exp`, `pow`) differing in
  the last digit between machines and tipping a value over a threshold.
  Recommended changes to the TODO text: the golden is pinned per
  `generator_epoch`, one digest for all legs, not per release with the
  version it was made on (REVISION moves every release, and per-leg
  goldens would hide cross-version drift); a database-free tier A (the
  OPS.15 battery, 3 to 5 s) also runs on the Windows and macOS CI jobs,
  tier B is the 1 against 4 workers galaxy, tier C the hash-seed and
  locale probe; one small sector is stored as full text so a failure can
  be read as a diff; a `golden_update` tool rewrites the goldens and bumps
  the epoch; `bump_version.py --check-pr` requires golden change, epoch
  bump and a `generation-output: changed` line in the note together.
  Details: [generation-determinism.md](generation-determinism.md)
  section 7.
- **Mixed versions (DB.7, built).** Each sector records the key, release,
  Python and platform that generated it (schema v67, `sectors.version_key`
  and friends). Extending a galaxy whose sectors came from a different
  version warns first, naming what differs ("PlanetGen 7.381.1 now,
  7.379.678 when generated"), since a mixed-version galaxy reproduces only
  sector by sector, each on its own version: `planetgen galaxy` prints it,
  and the Generate page shows it (`GET /api/galaxy/shape`'s
  `version_warning`). Sectors from before v67 have no version and are not
  compared.
  Not built, recommended: a `generator_epoch` (SMALLINT) and a
  `generation_runs` id per sector in place of repeated text columns, which at
  12 billion sectors repeat the same strings.

## 6. Updates and the version history (phase 1, then 2)

- **History (OPS.13, built).** After the databases are migrated,
  `update.sh` and `update.ps1` run `planetgen.cli.version_history`, which
  adds one row per planned galaxy to the control database's
  `version_key_history` (control schema v11): the galaxy seed (unchanged by
  an update), the key, the release, the SHA-256 of `requirements.lock` and
  the date. The corpus and name-list hashes are not recorded (the name
  codec replaced them). Only the last 10 rows per galaxy are kept;
  `planetgen versions` lists them. A failure to record only warns.
  Not built, recommended: the row also holds
  the `generator_epoch`, the battery digest and an environment record
  (libc, numpy, astropy, scipy, scikit-image versions). A row is added only
  when the key or a hash differs from the newest row (an update that
  changes nothing only refreshes that row's last-seen date), so repeated
  "already up to date" runs do not push the older distinct versions out.
  The lock file is pinned to LF in `.gitattributes` (or hashed as
  normalised text) so a Windows checkout hashes the same as a Linux one.
  *The seed question.* Boss's wording (2026-10-02 02:08Z) was that each
  update makes "a sweet to recalculate the seed value based on the current
  code version". The seed is not changed: it identifies the galaxy whose
  rows exist, and mixing the version or epoch into it would re-roll every
  sector at each epoch bump, against section 3's rule that a release
  changes only what its formula touches. The sweep recomputes and records
  what depends on the code (key, epoch, dependency versions and hashes,
  battery digest) beside the unchanged seed
  ([generation-determinism.md](generation-determinism.md) section 5.3).
- **Mismatch warning (OPS.14, built).** `version_check.galaxy_warning` is the
  one check: it names the sectors made by another version (DB.7) and how the
  running release, Python, platform, version key and `requirements.lock` hash
  differ from those in the galaxy's settings file (ADM.18). `planetgen
  galaxy`, `planetgen fingerprint` and the Generate page (`/api/galaxy/shape`'s
  `version_warning`) show it; a galaxy planned before the settings file has
  only the sector part. The reproduce report (GEN.61) will print it too. As
  designed: One check compares the running
  key and hashes with the galaxy's (and each sector's) and names every
  field that differs ("Python 3.12.3 now, 3.11.9 when generated").
  `planetgen`, the Generate page, the fingerprint output and the
  reproduce report all print it. Recommended severities: an epoch
  difference (output differs, and `epochs.json` names the release that
  introduced it); the epoch equal but the battery digest different (the
  environment or a library moved the output, naming the first differing
  battery case); everything else (release, Python micro, OS, architecture)
  as notes.
- **Did this update change output? (OPS.15, phase 2).** The update
  fingerprints a small fixed region from the galaxy seed under the new
  code, compares it with the previous history row's, stores it, and says
  "generated output unchanged" or names the sectors that differ. The
  live galaxy is not touched. Recommended: the compared fingerprint is the
  battery digest of [generation-determinism.md](generation-determinism.md)
  section 5.4, a pure database-free function over 6 to 8 tiny canned
  cases (3 to 5 s), compared together with the epoch, so a new Python or
  library that moves a value shows up even when no code changed.

## 7. The settings file (phase 1, then 2)

**Written at creation (ADM.18, phase 1).** When a galaxy is created a
JSON file is written holding only what reproduction needs:

- every creation setting, defaults included (all `plan` and new-galaxy
  options: disk scale length and height, bulge radius, max ring,
  bright-star floor, prevalence settings and the rest);
- the 128-bit seed and the version key with its parts spelled out;
- the corpus and lock hashes;
- the word list: the filtered list the name generator actually uses
  (nltk's `words`, about 2.5 MB raw, after `offensive_words.txt` and any
  name lists), gzip-compressed and base64-encoded, with its SHA-256. A
  rebuild reads names from this list, not from whatever nltk has
  installed;
- the naming key and the codec version ([object-ids.md](object-ids.md),
  GEN.70), for the objects the codec names;
- the `generator_epoch`, the fingerprint format version and a
  `content_sha256` of the canonical body, so the read-back check after a
  write is one comparison;
- the admin changes and the positional-update epoch (below).

The name is `<32-hex seed>-<22-hex key>-<YYYYMMDD>-<HHMMSS>Z.json` in UTC,
with no colons so it is valid on Windows, for example
`9F3A07C2E81B44D5A1C06E7B3D2F9081-0007007F000160030C0300-20261002-022133Z.json`.
It lives in the site's data directory (path in `config.json`). The Admin
dashboard offers the current file as a download at any time, with the
last 10 history rows beside it.

**As built (ADM.18).** `galaxy/settings_file.py` writes the file at the end of
`planetgen plan` (`run_plan.build_skeleton`), into `galaxy-settings` inside
the Generate page's jobs directory (`PLANETGEN_SETTINGS_DIR` overrides it).
It holds the `plan` options (`PLAN_SETTINGS`) plus the resolved edge, outer
ring and normalisation, the seed, the version key with its parts, the naming
key, the `requirements.lock` hash and both word lists (the dictionary and the
offensive list) as gzip + base64 with their SHA-256. The prevalence options
(ADM.45) are options of the generate runs, not the plan, so they are not in
the file yet; they join it when that item lands. A plan that changes nothing
writes nothing; a plan with other settings writes a new file and the old one
stays as a dated backup (the newest file for the seed is current). The Admin
dashboard's "Galaxy settings" panel lists the files and offers each as a
download (`GET /api/admin/galaxy-settings[/<name>]`). The admin-change
difference and the epoch below are GEN.59 and later.

**Admin changes as a net difference (GEN.59, phase 1).** The file keeps
what seed + key would not produce, not a history:

- each changed object's final state only: changed fields with their
  current values, deleted objects as tombstones;
- an edit that puts a value back to the original drops its entry;
- a regenerate draws a fresh random 128-bit seed for that object, which
  replaces its derived seed, and clears earlier entries for the object
  and its children; edits after it are recorded on top. (GEN.56 already
  removed the old `random.seed()` from `api/edits.py`. What is left is that
  `admin/edits.py: regenerate_planet`, `regenerate_moon` and `regenerate_belt`
  draw from the run stream with no stored seed, so GEN.59 adds a seed
  argument to them, binds `draw.bound(seed)` around the draws and stores
  the seed.)
- objects are named by a stable address path (sector ring, layer and
  slot, then system, body and moon by generated index), never by
  database id, since ids change on every rebuild. The index must be a
  stored generation index that is never renumbered, not a row rank:
  `store.assign_uids` and `galaxy/uid.py` rank an object among its
  parent's rows by row id, which can shift when an admin deletes one and
  the sector is saved again;
- the positional-update epoch: when `planetgen.cli.orbits` last moved the
  systems, so a rebuild reaches the same positions.

Admin changes go into a pending-delta table in the control database as
they happen. The file is not rewritten per change.

*Delta format (recommended).* One entry per changed object, keyed by
address path: `set` holds each changed field as `{was, now}` (`was` is
the generated value captured before the first edit, so a field put back
to its original drops out without regenerating); `regen` holds the new
128-bit seed; `deleted` marks a tombstone; `added` marks an admin-created
object, addressed with a random `+hex` index so it cannot collide with a
generated one. A regenerate or delete clears every entry under its path
(prefix operations), later edits go on top, and the merge folds pending
rows in time order, so the result is the same however they are grouped. A
custom form was chosen over RFC 6902 (it silently patched the wrong
planet once a later generator version inserted a system) and RFC 7396
(arrays are replaced whole and a field cannot be set to null).
Floats in the file are the exact repr, not rounded; a changed epoch
means the file is refused, since there is no automatic migration. Format,
rules and the demonstration: [generation-determinism.md](generation-determinism.md)
section 6.

**The daily merge (GEN.61, phase 2).** Boss (02:28Z): "a JSON file is
only changed with the deltas at the end of the day." A merge step reads
the pending deltas and the newest file, applies the net-difference rules,
records the epoch and writes a new file under the new date-time; the old
one stays as a backup, so the newest file is always the current one.
Pending rows are cleared only after the new file is written and read
back. With nothing pending and no epoch change it writes nothing.

Changes made since the last merge live only in the database until the
next one; the parity file (section 9) protects them in between.

**Merge now (ADM.20, phase 3, low priority).** Boss (02:31Z): "Let's
put that part of Phase 3, low priority." An admin-only button on the
Admin dashboard runs the delta merge at once, under the same lock and
rules as the daily run, and writes a new seed-key-date-time file. The
newest file of a day holds that day's daily backup slot, so the night
run's file replaces a merge-now file from earlier the same day. If the
daily run holds the lock, or a Generate job is running, the button says so
and does nothing. The merge runs as a job or thread, not inside the web
request, as the same account that runs maintenance (lock file and JSON
ownership). Needs OPS.16, GEN.61 and OPS.18.

## 8. The daily maintenance run (phase 2)

- **The script (OPS.16).** The logic lives once, in
  `planetgen.cli.maintenance`; `scripts/maintenance.sh` (Linux and macOS)
  and `scripts/maintenance.ps1` (Windows) are thin launchers that run it as
  the web user. It runs once a day: the positional update
  (`planetgen.cli.orbits`, for every galaxy database), then the delta
  merge (GEN.61), then the backup rotation (OPS.18). A lock file stops two
  runs overlapping, and the run is skipped (exit 75) while a Generate job
  is running; each step logs and any failure exits non-zero. `--if-due`
  makes the run a no-op when the last success is under about 20 hours old,
  which enforces "once every 24 hours" however often the scheduler fires.
  Redis is optional here, so an unreachable Redis is a warning. Optionally
  it runs OPS.15's fingerprint check too.
- **The schedule (OPS.17).** Install and update set it up by calling the
  installers that already exist in `examples/maintenance/`: a systemd
  timer (cron.d where there is no systemd) on Linux, launchd on macOS, Task
  Scheduler on Windows. The daily run replaces the monthly orbit timers,
  tasks and plists (`planetgen-orbits@`, `org.planetgen.orbits.*`), which
  are removed when it is installed, since it runs orbits itself. Update
  keeps an existing schedule and adds a missing one; a
  `maintenance.schedule` key in `config.json` (`auto` or `off`) is the off
  switch that survives updates. The deployment docs say how to change the
  time. Lands after OPS.7, OPS.8 and OPS.13, which change the same scripts.
- **18 backups (OPS.18).** After each merge the files are kept by
  grandfather-father-son rotation with calendar periods taken from the UTC
  time in each file name: the newest file of each of the last 7 days, 4 ISO
  weeks, 6 months and 1 year that have a file, a file kept by a finer tier
  not counting again. That is exactly 18 distinct files once history
  exists; the yearly file is the last file of a past calendar year, 188 to
  553 days old, not exactly 365. Older files are deleted, the current one
  is always kept, and only parsed names of the same seed are touched.
- **Listing them (ADM.19).** The Admin dashboard lists all 18 with date,
  slot and version key, each downloadable. The slot comes from the same
  selection function as the rotation, not from a stored label.

The reload step, the per-OS schedulers, the lock, the rotation function and
its tests are in [ops-scheduling-and-rotation.md](ops-scheduling-and-rotation.md).

## 9. Damage: check and repair

- **Check (DB.8, phase 0).** `planetgen check-db` and an Admin
  dashboard button (run as a job) check the galaxy and control
  databases without changing anything: the Alembic revision (`alembic_version`,
  `schema_migrations` and `store.SCHEMA_VERSION` agree, and the live shape
  matches `db/models.py` through `compare_metadata`), `CHECK TABLE ... QUICK`
  per table, orphan rows (foreign-key queries generated from
  `information_schema`, plus hand-written ones for the polymorphic
  `object_table`/`object_id` columns), ids against `id_blocks`, sector
  addresses in bounds, stored values inside plausibility ranges (a DOUBLE
  column cannot hold NaN or infinity, so a NaN check can only run before a
  save), every system passing `validation.check_star_system`, and, once they
  exist, the per-sector stats (GEN.44, PERF.11), version keys (DB.6, DB.7)
  and per-sector checksums (DB.9). The name-registry check goes with the
  registries (GEN.71). One pass or fail line per check; non-zero exit on
  damage, a different code for "could not check"; `--sector` and `--region`
  limit it. Queries, costs and lock behaviour:
  [db-check-and-parity-repair.md](db-check-and-parity-repair.md) section 6.
- **Repair (DB.9, phase 1).** Each sector gets a storage checksum: SHA-256
  over a lossless, id-preserving dump of its stored rows, read back after
  the write. It is not the GEN.58 fingerprint, which is rounded to 9 digits,
  skips ids and describes the sector as generated, so it would miss a
  flipped low bit and would not match an edited sector. Neither the checksum
  nor the parity covers the columns the daily orbit update rewrites (phases,
  positions, velocities, reflex offsets, a placed object's sector), which
  would make every record stale each day. A parity file outside the database
  (path in `config.json`) holds Reed-Solomon parity over groups of sector
  exports, so one damaged sector per group (two with the recommended
  G = 32, m = 2) can be rebuilt; groups are sectors of similar export size
  far apart in id (about 7 percent overhead), it is updated by delta as
  sectors are saved or edited, and each record carries its own checksum.
  `planetgen repair-db` finds damage with DB.8's check and rebuilds from
  parity; where parity can't and the key matches the running code (OPS.14),
  it regenerates the sector from its seed and replays its edits. It checks
  again and lists anything it could not repair. A damaged clustered-index
  page cannot be read or deleted through SQL, so for that case repair
  rebuilds the affected table from the surviving rows and the rebuilt
  exports. Codec, layout, stale-record table and the page-damage behaviour:
  [db-check-and-parity-repair.md](db-check-and-parity-repair.md) sections 4
  and 5.
- **Repair with the settings file (DB.10, phase 3).** Where repair
  regenerates a sector, it applies the newest file's diff and epoch plus
  any pending deltas, so changes since the last merge survive; if the
  newest file is damaged it falls back to the next backup and says so.

The seed and the parity file cover different risks: the seed rebuilds
the galaxy as generated, with admin changes replayed on top; the parity
file repairs damage to the database as it was at the last save, with no
rebuild (the orbit-updated columns excepted: how repair restores them is
open, [db-check-and-parity-repair.md](db-check-and-parity-repair.md)
section 5.2).

## 10. The end state: `planetgen reproduce` (OPS.12, phase 3+)

`planetgen reproduce --seed X --version Y` rebuilds a galaxy or region
into a fresh database from the seed, the run history and the settings
file, then compares fingerprints with the live galaxy (or a given one)
and lists any sector that differs. It refuses, naming the release to
check out, when the running release isn't Y, and prints OPS.14's
comparison of keys and hashes. Without the settings file it rebuilds the
galaxy as first generated; with it, as it is now. It reads the newest
file plus any pending deltas, or any of the 18 backups (`--as-of DATE`).
No automatic migration of old galaxies to a new release's output.
Recommended: it refuses on an epoch mismatch (naming the first release of
that epoch from `epochs.json`) or a battery-digest mismatch rather than on a
release mismatch, and rebuilds per address, which GEN.57's order
independence allows, rather than replaying commands where it can.

Seed and version also reach the web and API: the Generate page shows
the seed (with a copy button), the version and the run history; its
new-galaxy form takes no seed (Boss, 2026-10-09 20:52Z); a remote run with the same seed and
release produces exactly what the server would, checked by fingerprint
(API.17, phase 3). Anything that draws new randomness later (GEN.47,
GEN.42, PERF.18, API.12, API.13) uses the derived seeds and keeps
TEST.77 green.

## 10a. Changes from the 2026-10-07 plan

- IDs and a naming key (object-ids.md) name only the interstellar objects
  and constellations; stars and sectors keep word-salad names (Boss,
  2026-10-08 03:57Z and 04:00Z). So the settings file stores the naming key
  and codec version in addition to the filtered word list, the update
  history keeps the corpus and name-list hashes along with the lock-file
  hashes, and the name-collision rule in section 5 stays for stars and
  sectors.
- Schema changes are Alembic migrations, and the consistency check reads
  Alembic's revision.
- Generation runs as RQ jobs on Redis; results must not depend on the
  worker count or the order jobs finish.

## 11. Order

| Phase | Items | Needs |
|---|---|---|
| 0 | GEN.39 (per-unit seeds), DB.6 (seed, key, run history), OPS.10 (log line) | PERF.21; then in that order |
| 0 | DB.8 (consistency check) | none; later checks switch on as DB.6, DB.7, GEN.44 and PERF.11 land |
| 1 | Deterministic-math helpers (proposed item) | GEN.56; before TEST.77 |
| 1 | OPS.11 (this note), GEN.56, GEN.57, DB.7, GEN.58, TEST.77 | GEN.39 and DB.6; GEN.57 also GEN.46 and GEN.44 |
| 1 | OPS.13, OPS.14 | DB.6, OPS.7, OPS.8 |
| 1 | ADM.18, GEN.59 | DB.6, DB.7, OPS.13; GEN.59 also GEN.56, GEN.58 |
| 1 | DB.9 (parity repair) | DB.8, GEN.39, GEN.57, GEN.44, GEN.58, OPS.14 (recommended: the parity half needs only DB.8, and the regenerate-from-seed fallback needs the rest) |
| 2 | GEN.61, OPS.18, OPS.16, OPS.17, ADM.19 | GEN.59, ADM.18; OPS.17 also OPS.13 |
| 2 | OPS.15, ADM.17, API.16 | OPS.13 and GEN.58; DB.6; DB.6 and API.5 |
| 3 | API.17, DB.10, ADM.20 (merge now, low priority) | API.12, API.13, GEN.57, GEN.58; DB.9, GEN.61, OPS.18; OPS.16, GEN.61, OPS.18 |
| 3+ | OPS.12, closing GEN.55 | everything above it |

See also [ops-scheduling-and-rotation.md](ops-scheduling-and-rotation.md) (the daily maintenance run and the backup rotation, OPS.16 to OPS.18).
