# Reproducible galaxies

How a galaxy will be rebuilt from a 128-bit seed and a version key, how
admin changes and the passage of time are kept on top of that, and how a
damaged database is found and repaired. Recorded 2026-10-02 from Boss's
decisions of 01:34Z to 02:31Z and the seed math report (dependency tree
thread). This is the design note OPS.11 asks for.

**Status (2026-10-02, checked against main after PR #359, 7.132.433):
planned, nothing in this note is built.** Today every generation draw
comes from the operating system's random source, so no seed reproduces a
galaxy (GEN.39). Each piece below names its TODO item and phase;
`docs/TODO.md` holds each item's full text and `docs/plan/` the order.

## 1. What Boss asked for

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
own release. A seed cannot be worked out for a galaxy that already
exists: today's galaxy came from true randomness, and finding a seed for
a given galaxy means running SHA-256 backwards (about 2^128 tries).

## 3. The seed (GEN.39, DB.6, phase 0)

- One 128-bit galaxy seed, stored in `galaxy_shape` as `BINARY(16)`,
  shown as 32 uppercase hex digits, written once when the galaxy is
  first planned. `generate.py plan --seed <32 hex>` sets it; without it
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
- The bright-star scatter's 63-bit seed (today
  `random.SystemRandom().getrandbits(63)` in `generate.py`) is derived
  from the galaxy seed too, so no step throws bits away.
- Odds, from the seed math report: two of a billion random 128-bit seeds
  match with a chance of about 1.5 x 10^-21; two of 12 billion sectors
  sharing a unit seed, about 10^-57.

GEN.39 builds after PERF.21 in the parallel path thread. Its tests: one
seed gives the same sectors at 1, 2 and N workers, and two seeds that
differ only in their high 64 bits give different output.

As built (GEN.39): `stellarObjects/galaxySeed.py` holds the seed helpers.
A sector fill (its save included) and the bright-star scatter, bands
and backfill blocks run on their unit seed; the sector seeds the
module-level `random` stream for the length of the unit
(`galaxySeed.seeded`, which puts the stream back afterwards) rather than
passing a `random.Random` down, which is GEN.56's change. `secrets` and
`SystemRandom` are gone from generation (`utils.reseed_rng` is removed,
`spaceSector._rng` reads the module stream). A save retried after a
deadlock replays the same draws. The `phenomenon` and `system` subcommands and the run's own
choices (a random start's address) still draw from the run's stream.

As built (DB.6): `stellarObjects/versionKey.py` computes the key
(`version_key`) and the versions stored beside it (`current`). The galaxy
records them in `galaxy_shape` whenever its seed is written (schema
v52), and every `generate.py` run that changes the galaxy writes a
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
the settings file (section 7). No numpy or other numeric library is
used, so the OS, architecture and Python parts cover the maths library.

## 5. Every draw seeded, every sector independent (phase 1)

- **The log line (OPS.10, phase 0).** Every `generate.py` subcommand,
  web or API job and work queue run first writes, at normal level, the
  galaxy seed, the version with its key, and the run's command, for
  example `Galaxy seed 3f2a...c901, PlanetGen 7.127.352
  (0007007F000160030C0300), run: sector 12 3 0`. Today `generate.py`
  logs its `secrets.randbits(128)` run seed at debug level only.
  As built (OPS.10): `versionKey.run_line` writes it, first thing after
  `generate.py` sets up logging (a `plan --seed` names the seed given,
  otherwise the stored one, `none yet` before the first plan), and
  `jobRunner.py` writes it at the top of a job's `output.log` with the
  job's title. A work queue always runs inside a `generate.py` run, so
  its line is that run's.
- **The run history (DB.6, phase 0).** A `generation_runs` table, one
  row per run that changes the galaxy (command and options, version,
  start and end, outcome). A galaxy is built by a series of commands
  (plan, sectors, scatter, backfill), not by the seed alone, so this is
  what makes a rebuild possible.
- **Every draw from the derived seeds (GEN.56).** Today `secrets.choice`
  and `secrets.randbelow` pick planet classes, life, moons and name
  syllables (`planetPhysics.py`, `planetLife.py`, `planetData.py`,
  `utils.py`); `spaceSector._rng` is a `secrets.SystemRandom`; and
  `utils.reseed_rng()` reseeds the global `random` from `secrets` in 13
  generator modules. All of that moves to the unit's `random.Random`.
  Python promises the same sequence only for `random()` itself, not for
  `choice`, `uniform`, `gauss` or `shuffle`, so generators draw through
  small in-house helpers built on `random()` alone. Nothing may depend on
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
- **Fingerprint (GEN.58).** `generate.py fingerprint` prints a canonical
  SHA-256 per sector and for a region, in address order with canonical
  number formatting, skipping ids and timestamps, either as first
  generated or with the settings file applied.
- **Golden test (TEST.77).** A fixed seed builds a small galaxy (plan, a
  few sectors, a scatter, a backfill) at 1 and 4 workers on every CI
  Python leg, and its fingerprint must match the one pinned in
  `src/tests/golden/`. A PR that changes generation output updates the
  pinned fingerprint and says so in its `changes/` note;
  `bump_version.py --check` checks the two go together. It also watches
  for the one known risk: a maths function (`exp`, `pow`) differing in
  the last digit between machines and tipping a value over a threshold.
- **Mixed versions (DB.7).** Extending a galaxy with a different release
  warns first, since a mixed-version galaxy reproduces only sector by
  sector, each on its own version.

## 6. Updates and the version history (phase 1, then 2)

- **History (OPS.13, phase 1).** After updating the code, `update.sh`
  and `update.ps1` compute the running key and add one row per galaxy
  to a control-database table: the galaxy seed (unchanged by an update),
  the key, the corpus and lock hashes, and the date. Only the last 10
  rows per galaxy are kept; `generate.py` can list them. Lands after
  OPS.7 and OPS.8, which change the same scripts.
- **Mismatch warning (OPS.14, phase 1).** One check compares the running
  key and hashes with the galaxy's (and each sector's) and names every
  field that differs ("Python 3.12.3 now, 3.11.9 when generated").
  `generate.py`, the Generate page, the fingerprint output and the
  reproduce report all print it.
- **Did this update change output? (OPS.15, phase 2).** The update
  fingerprints a small fixed region from the galaxy seed under the new
  code, compares it with the previous history row's, stores it, and says
  "generated output unchanged" or names the sectors that differ. The
  live galaxy is not touched.

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
- the admin changes and the positional-update epoch (below).

The name is `<32-hex seed>-<22-hex key>-<YYYYMMDD>-<HHMMSS>Z.json` in UTC,
with no colons so it is valid on Windows, for example
`9F3A07C2E81B44D5A1C06E7B3D2F9081-0007007F000160030C0300-20261002-022133Z.json`.
It lives in the site's data directory (path in `config.json`). The Admin
dashboard offers the current file as a download at any time, with the
last 10 history rows beside it.

**Admin changes as a net difference (GEN.59, phase 1).** The file keeps
what seed + key would not produce, not a history:

- each changed object's final state only: changed fields with their
  current values, deleted objects as tombstones;
- an edit that puts a value back to the original drops its entry;
- a regenerate draws a fresh random 128-bit seed for that object, which
  replaces its derived seed, and clears earlier entries for the object
  and its children; edits after it are recorded on top (this replaces
  today's `random.seed()` in `api/edits.py`);
- objects are named by a stable address path (sector ring, layer and
  slot, then system, body and moon by generated index), never by
  database id, since ids change on every rebuild;
- the positional-update epoch: when `updateOrbits.py` last moved the
  systems, so a rebuild reaches the same positions.

Admin changes go into a pending-delta table in the control database as
they happen. The file is not rewritten per change.

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
rules as the daily run, and writes a new seed-key-date-time file that
counts toward that day's daily backup slot. If the daily run holds the
lock, the button says so and does nothing. Needs OPS.16, GEN.61 and
OPS.18.

## 8. The daily maintenance run (phase 2)

- **The script (OPS.16).** `scripts/maintenance.sh` (Linux and macOS)
  and `scripts/maintenance.ps1` (Windows) run once a day: the positional
  update (`updateOrbits.py`), then the delta merge (GEN.61), then the
  backup rotation (OPS.18). A lock stops two runs overlapping, each step
  logs, and any failure exits non-zero. Optionally it runs OPS.15's
  fingerprint check too.
- **The schedule (OPS.17).** Install and update set it up: a systemd
  timer or cron entry on Linux, launchd on macOS, Task Scheduler on
  Windows. Update keeps an existing schedule and adds a missing one; the
  deployment docs say how to change the time or turn it off. Lands after
  OPS.7, OPS.8 and OPS.13, which change the same scripts.
- **18 backups (OPS.18).** After each merge the files are kept by
  grandfather-father-son rotation: the newest 7 daily, then 4 weekly, 6
  monthly and 1 yearly, 18 in all; older files are deleted and the
  current one is always kept.
- **Listing them (ADM.19).** The Admin dashboard lists all 18 with date,
  slot and version key, each downloadable.

## 9. Damage: check and repair

- **Check (DB.8, phase 0).** `generate.py check-db` and an Admin
  dashboard button (run as a job) check the galaxy and control databases
  without changing anything: schema against the recorded migration
  level, orphan rows, ids against `id_blocks`, name registries against
  names in use, sector addresses in bounds, no NaN or infinite values,
  every system passing `validation.check_star_system`, and, once they
  exist, the per-sector stats (GEN.44, PERF.11) and version keys (DB.6,
  DB.7). One pass or fail line per check; non-zero exit on damage;
  `--sector` and `--region` limit it.
- **Repair (DB.9, phase 1).** Each sector gets a content checksum (the
  hash of its GEN.58 fingerprint). A parity file outside the database
  (path in `config.json`) holds Reed-Solomon parity over groups of
  sector exports, so one damaged sector per group can be rebuilt; it is
  updated as sectors are saved or edited and carries its own checksum.
  `generate.py repair-db` finds damage with DB.8's check and rebuilds
  from parity; where parity can't and the key matches the running code
  (OPS.14), it regenerates the sector from its seed and replays its
  edits. It checks again and lists anything it could not repair.
- **Repair with the settings file (DB.10, phase 3).** Where repair
  regenerates a sector, it applies the newest file's diff and epoch plus
  any pending deltas, so changes since the last merge survive; if the
  newest file is damaged it falls back to the next backup and says so.

The seed and the parity file cover different risks: the seed rebuilds
the galaxy as generated, with admin changes replayed on top; the parity
file repairs damage to the database as it is now, with no rebuild.

## 10. The end state: `generate.py reproduce` (OPS.12, phase 3+)

`generate.py reproduce --seed X --version Y` rebuilds a galaxy or region
into a fresh database from the seed, the run history and the settings
file, then compares fingerprints with the live galaxy (or a given one)
and lists any sector that differs. It refuses, naming the release to
check out, when the running release isn't Y, and prints OPS.14's
comparison of keys and hashes. Without the settings file it rebuilds the
galaxy as first generated; with it, as it is now. It reads the newest
file plus any pending deltas, or any of the 18 backups (`--as-of DATE`).
No automatic migration of old galaxies to a new release's output.

Seed and version also reach the web and API: the Generate page shows
the seed (with a copy button), the version and the run history, and its
new-galaxy form takes an optional seed (ADM.17, phase 2); an API route
returns the same (API.16, phase 2); a remote run with the same seed and
release produces exactly what the server would, checked by fingerprint
(API.17, phase 3). Anything that draws new randomness later (GEN.47,
GEN.42, PERF.18, API.12, API.13) uses the derived seeds and keeps
TEST.77 green.

## 11. Order

| Phase | Items | Needs |
|---|---|---|
| 0 | GEN.39 (per-unit seeds), DB.6 (seed, key, run history), OPS.10 (log line) | PERF.21; then in that order |
| 0 | DB.8 (consistency check) | none; later checks switch on as DB.6, DB.7, GEN.44 and PERF.11 land |
| 1 | OPS.11 (this note), GEN.56, GEN.57, DB.7, GEN.58, TEST.77 | GEN.39 and DB.6; GEN.57 also GEN.46 and GEN.44 |
| 1 | OPS.13, OPS.14 | DB.6, OPS.7, OPS.8 |
| 1 | ADM.18, GEN.59 | DB.6, DB.7, OPS.13; GEN.59 also GEN.56, GEN.58 |
| 1 | DB.9 (parity repair) | DB.8, GEN.39, GEN.57, GEN.44, GEN.58, OPS.14 |
| 2 | GEN.61, OPS.18, OPS.16, OPS.17, ADM.19 | GEN.59, ADM.18; OPS.17 also OPS.13 |
| 2 | OPS.15, ADM.17, API.16 | OPS.13 and GEN.58; DB.6; DB.6 and API.5 |
| 3 | API.17, DB.10, ADM.20 (merge now, low priority) | API.12, API.13, GEN.57, GEN.58; DB.9, GEN.61, OPS.18; OPS.16, GEN.61, OPS.18 |
| 3+ | OPS.12, closing GEN.55 | everything above it |
