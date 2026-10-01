# planetGen

**Version:** 7.47.0 &middot; [Changelog](CHANGELOG.md) &middot; [Repository](https://github.com/dwhagar/planetGen) &middot; License: [CC0 1.0 Universal](LICENSE.md)

A procedural planet and star system generator, designed for the Molten Aether FFRP game. The output is designed to be easily copied and pasted into the wiki.

## Features

*   **Star generation**: Stars are generated with physically-grounded mass, radius, temperature, and luminosity, either fully at random (weighted by realistic galactic spectral-class prevalence) or pinned to a specific spectral type and Yerkes luminosity class (e.g. `G2V`) via `--star-type`. `+large_star` biases generation toward hotter, more massive stars.
*   **Binary star systems**: `+binary_system` generates a binary pair — a primary and secondary star orbiting each other — in one of the two real binary-star configurations, selected by `+wide_binary`/`-wide_binary` (random if omitted): a **P-type (close/circumbinary)** pair, represented as a single unified effective star (combined mass, luminosity, and habitable zone) for the purposes of planet placement, while still reporting each star's individual properties in the output; or an **S-type (wide)** pair, separated by tens to thousands of AU, where each star keeps its own identity and hosts its own independently-generated planets, with each star's maximum stable orbit limited by the companion's gravity (Holman & Wiegert 1999) and a mutual-Hill-radius safety check (Gladman 1993) preventing the two stars' own planetary disks from encroaching on each other.
*   **Planet and moon generation**: Planets are drawn from 25 distinct planet classes (terrestrial and gas giant), each with its own composition, atmosphere, and valid orbital zones (hot/ecosphere/cold). Planets can generate their own moons, with orbital placement, atmospheric conditions, and surface gravity calculated per body. Orbital spacing is validated against each object's Hill sphere to keep the system physically plausible.
*   **Asteroid belts**: Belts can appear between planets (or be forced/forbidden with `+asteroid_belt`/`-asteroid_belt`), each with a randomly generated density and mineral/gem composition.
*   **Explicit orbit/slot specification**: The number of orbital slots (planets and asteroid belts combined) can be pinned exactly with `--num-orbits`, and a `--system-file` JSON specification can dictate the exact contents of any specific orbital slot — whether it's a planet or an asteroid belt, the planet's class, and how many moons it has — leaving unspecified slots to normal random generation.
*   **Stellar age and lifespan modeling**: Star age and lifespan are derived from spectral and Yerkes class data, then adjusted (via `--age young|old`) or extended as needed so any habitable planets' life stages remain consistent with how long the star has existed and how much longer it has left. A star's age is always capped at the actual age of the universe (~13.8 billion years), even for classes (like M dwarfs) whose theoretical lifespan runs to trillions of years.
*   **Life chemistry and evolutionary timelines**: Habitable planets are evaluated against the star's spectral class to determine which of several life chemistries (e.g. Chlorophyll a, Melanin, Retinal) are viable, each with its own evolutionary pace (fast/normal/slow). A speculative evolutionary timeline (abiogenesis through technological civilization) is generated for habitable worlds, influenced by `+intelligent_life` / `-intelligent_life`.
*   **Flavor text**: Randomly-selected descriptive "sensor readings" flavor text can be appended to systems and planets, with limits and chances controllable via `--flavor-chance-system`, `--flavor-chance-planet`, and `--max-planet-flavor`.
*   **Dual output formatting**: Every generated system can be rendered as either MediaWiki wikitext templates (default, ready to paste into the wiki) or Markdown (`--markdown`).
*   **Sector generation**: `generate.py sector` generates a whole sector of independently-random star systems in one pass, with an optional guaranteed minimum number of habitable systems (`--min-habitable`) — see [Sector Generation](#sector-generation) below. Each system placed in a sector records a `location`: its sector's name plus distance (in light-years) to its 3 nearest neighboring systems.
*   **Exotic stellar phenomena**: `generate.py phenomenon` generates a single black hole, neutron star, nebula, supernova remnant, rogue planet, interstellar comet, standalone asteroid field, or quasar on demand — kept separate from normal system generation odds (`generate.py system`/`sector` never produce one) — see [Exotic Phenomena Generation](#exotic-phenomena-generation) below. Every standalone phenomenon except a quasar still orbits the galactic center like a lone star does, advanced over time the same way by `updateOrbits.py`; a quasar is the galaxy's own active nucleus and sits at the center itself.

## Setup

This project uses the `nltk` library to generate phonetically pleasing names. Install it (and the rest of the project's dependencies) with:

```bash
pip install .
```

run from the root of the project. This installs `nltk` (the 'words' corpus it needs for name generation is then downloaded automatically, the first time it's needed, by whichever user first imports `stellarObjects`), plus `pymysql`/`DBUtils` for database persistence. (Do not run `python setup.py install` directly -- it's deprecated by setuptools itself and, on some systems, an old system-installed copy of a dependency setuptools vendors internally can shadow the working bundled one and crash the install; a normal `pip install` avoids this by building in an isolated environment. See [`CHANGELOG.md`](CHANGELOG.md) for the incident this was fixed after.)

Every generation run saves to a MySQL database (see [`docs/database-schema.md`](docs/database-schema.md)) — you'll need a MySQL server (8.0.16+, for `CHECK` constraint support) reachable from wherever you run `generate.py`, with a database/user already created. Point it at your database with the `$PLANETGEN_MYSQL_HOST`/`$PLANETGEN_MYSQL_PORT`/`$PLANETGEN_MYSQL_USER`/`$PLANETGEN_MYSQL_PASSWORD`/`$PLANETGEN_MYSQL_DATABASE` environment variables (or the equivalent `--mysql-*` flags every subcommand accepts) — tables are created automatically on first connection.

Deploying the web interface to a server is a separate, more involved process — see [Web Interface](#web-interface) below and the [deployment guides](docs/deployment/README.md).

## Usage

`generate.py` is the single command-line entry point for every generator in this project — one script, with a subcommand per generation scale:

```bash
python generate.py system [options]      # one star system
python generate.py sector [options]      # one or more independent sectors
python generate.py galaxy [options]      # many sectors placed as one galaxy
python generate.py plan [options]        # the galaxy's density skeleton
python generate.py phenomenon [options]  # one exotic stellar phenomenon
```

Run `python generate.py <command> --help` for that command's own full option list. Every subcommand saves what it generates to the database.

To generate a new star system:

```bash
python generate.py system [options]
```

To generate a whole sector of star systems at once, see [Sector Generation](#sector-generation) below.

### Options

Most generation options use a `+name`/`-name` tri-state syntax: `+name` forces that feature to be present, `-name` forces it to be absent, and leaving it off the command line leaves it up to chance.

*   `--version`: Prints the program's version, repository URL, and license, then exits immediately.
*   `+habitable_world` / `-habitable_world`: Force / forbid the generation of a habitable world in the system.
*   `+asteroid_belt` / `-asteroid_belt`: Force / forbid the generation of an asteroid belt in the system.
*   `+large_star` / `-large_star`: Force / forbid the generation of a large star.
*   `+moons` / `-moons`: Force / forbid moons on the system's planets.
*   `+max_planets` / `-max_planets`: Force the system to the maximum, or the minimum, number of orbital objects it can support.
*   `+intelligent_life` / `-intelligent_life`: Ensure a planet with intelligent life is (or is not) generated. Either implies `+habitable_world`.
*   `+binary_system` / `-binary_system`: Force / forbid a binary star system.
*   `+wide_binary` / `-wide_binary`: Force an S-type (wide) binary, or force a P-type (close) one instead. Only meaningful together with `+binary_system`; if omitted, one of the two is chosen at random.
*   `+planets` / `-planets`: Ensure the system has at least one planet or asteroid belt, or none at all (star only).
*   `--system-file`, `-f <path>`: Load a system generation specification from a JSON file (see below). Any of the options above, or the value options below, given on the command line override the corresponding value from the file.
*   `--num-orbits <int>`: Force an exact number of orbital slots (planets and asteroid belts combined) to be generated.
*   `--markdown`, `-m`: Formats the stored write-up in Markdown instead of the default wikitext.
*   `--output`, `-o <file>`: Writes the system's page (wikitext, or Markdown with `--markdown`) to `file` instead of saving the system to the database; `-` writes it to stdout. The admin site's one-off system page (`/admin/generate/system`) does the same from the browser.
*   `--star-type <type>`: Force the generation of a specific star type (e.g., G2V).
*   `--name <name>`: Specifies a name for the star system, overriding the default random generation. Everything in the system is named from it: planets are numbered in orbit order (`Sol I`, `Sol II`, ...), moons add a letter (`Sol IIIa`), and a binary's two stars each put their own word after the system name (`Sol Kelmoor`, `Sol Ostra`), with no A/B letters. A wide binary numbers each star's planets after that star (`Sol Kelmoor I`). Asteroid belts aren't numbered.
*   `--age <young|old>`: Specifies the age of the star system (young or old).
*   `--flavor-chance-system <float>`: Overrides the default system-level flavor text chance (0.0 to 1.0).
*   `--flavor-chance-planet <float>`: Overrides the default planet-level flavor text chance (0.0 to 1.0).
*   `--max-planet-flavor`: Sets the maximum flavor text total for planets to 99.
*   `--debug [file]`: Logs every choice the generator makes, and why, with timestamps, to the console. If `file` is given, also mirrors that output to `file`. Available on every subcommand. For a permanent, far more detailed log of everything (every random roll, SQL statement and web request too), set `"debug": true` in `config.json` instead; it writes to `/var/log/planetgen.log` (see [`docs/config.md`](docs/config.md)).
*   `--quiet` / `--silent`: Suppresses all output except errors. Available on every subcommand. Combined with `--debug`, a filename is required (there would otherwise be nowhere for debug output to go).

**Note on Incompatible Options:**

*   `-planets` cannot be combined with `+moons`, `+max_planets`, or `+habitable_world`.
*   `--star-type` cannot be combined with `+large_star`.
*   `+intelligent_life`/`-intelligent_life` cannot be combined with `-habitable_world`.
*   `+habitable_world` and `+asteroid_belt` together cannot be combined with `-large_star` (both objects require the room a large star provides).
*   `--num-orbits` cannot be combined with `-planets`.
*   `--num-orbits` must be from 0 to 500 (the generator's own ceiling on objects in a system).
*   `--flavor-chance-system` must be a float between 0.0 and 1.0.
*   `--flavor-chance-planet` must be a float between 0.0 and 1.0.

### System specification files

`--system-file`/`-f` loads a JSON file describing exactly what to generate. Any key may be omitted, in which case that aspect is generated normally. The tri-state options above become `true`/`false`/omitted JSON keys, e.g.:

```json
{
  "star_type": "G2V",
  "name": "Sol",
  "age": "young",
  "num_orbits": 5,
  "habitable_world": true,
  "asteroid_belt": true,
  "slots": [
    {"type": "planet", "planet_class": "M", "moons": 1},
    {"type": "asteroid_belt"},
    null,
    {"type": "planet", "planet_class": "J", "moons": 4},
    null
  ]
}
```

`slots` is an optional, per-orbit list: each entry is either `null` (generate that slot normally) or an object with a required `"type"` (`"planet"` or `"asteroid_belt"`) and, for planets, an optional `"planet_class"` (e.g. `"M"`) and/or `"moons"` (an exact moon count, `0` for none). The list doesn't need to cover every orbit — slots past the end of the list are generated normally too.

### Sector Generation

`generate.py sector` generates a whole sector of independently-random star systems in one pass, reusing `generate.py system`'s own generation logic for each one:

```bash
python generate.py sector [options]
```

Most of `system`'s options work here too, but apply *uniformly* to every system in the sector — `+asteroid_belt` guarantees a belt in every system, `--star-type G2V` makes every star in the sector a G2V, and so on. `--system-file`, `--num-orbits`, and `system`'s own per-system `--name` aren't offered here, since those describe one specific, hand-crafted system rather than a sector of varied ones; use `generate.py system --system-file` directly for that.

Sector-specific options:

*   `--num-systems <int>`: How many star systems the sector contains. Defaults to 10.
*   `--name <name>`, `-n <name>`: Hard-sets the sector's own name, overriding the default random two-word name (e.g. `"Voranthis Kelmoor"`) generated the same phoneme-salad way as star names.
*   `--min-habitable <int>`: Guarantees at least this many systems in the sector have a habitable world, chosen randomly among them — without forcing *every* system to have one the way a uniform `+habitable_world` would. Extra systems can still turn out habitable by chance on top of this minimum. Cannot exceed `--num-systems`, and cannot be combined with a uniform `-habitable_world`.
*   `--mysql-host <host>`, `--mysql-port <port>`, `--mysql-user <user>`, `--mysql-password <password>`, `--mysql-database <database>`: Where the generated sector is saved. Each defaults to the matching `$PLANETGEN_MYSQL_*` environment variable, or a built-in default (`127.0.0.1:3306`, user/database `planetgen`) -- see [`docs/database-schema.md`](docs/database-schema.md).

Each run saves the whole generated sector — every system, star, planet, moon, and asteroid belt, plus any exotic phenomena and a rendered copy of the wiki page in both wikitext and Markdown — to the MySQL database described in [`docs/database-schema.md`](docs/database-schema.md), printing only a short status line and summary per saved sector to the console (see "Logging" below for `--debug` if you want to see more).

### Sector-Level Exotic Phenomena

Every generated sector also seeds a realistically sparse population of `generate.py phenomenon`'s own seven exotic phenomenon types (black holes, neutron stars, nebulae, supernova remnants, rogue planets, interstellar comets, standalone asteroid fields) — sampled independently per type via a Poisson draw whose mean is a real (or, where flagged, deliberately conservative) astrophysical rate per star system, scaled by however many systems the sector actually ended up with (see `program_constants.PHENOMENON_RATE_PER_STAR_SYSTEM` for where each rate comes from and its citations). At this generator's own sector scale, most of these rates are low enough that a typical sector shows none at all — which is realistic; real space this size is usually devoid of black holes, neutron stars, and visible nebulae, exactly as it's usually devoid of Alpha-Centauri-close star systems. Separately, a galaxy's nucleus is active 10% of the time (`program_constants.QUASAR_ACTIVE_NUCLEUS_CHANCE`): when it is, the core sector at ring 0, layer 0, slot 0 gets a quasar at the galactic center, so a galaxy never has more than one.

Black holes and neutron stars are real, stellar-mass gravitating bodies, so they're placed the same Hill-sphere-aware way every star system itself already is: never within a neighboring star system's or another compact remnant's own Hill sphere (see "Minimum separation (Hill spheres)" in `stellarObjects/spaceSector.py`'s own module docstring). Every other phenomenon type has no comparable gravitational footprint at this scale and is placed at a random point in the sector's cell instead.

## Galaxy Generation

`generate.py galaxy` generates and persists many sectors as one galaxy, placed
in real galaxy-frame 3D space (see
[`docs/design/galaxy-coordinate-system.md`](docs/design/galaxy-coordinate-system.md)):

```bash
python generate.py galaxy --ring I [--layer J] [--slot K [--radius-pc R]] [options]
python generate.py galaxy --ring I --slot K --column [options]
python generate.py galaxy --ring I --shell [--limit N | --yes] [options]
python generate.py galaxy --center-sector ID --radius-pc R [options]
python generate.py galaxy [options]
```

It reuses `generate.py sector`'s own generation logic per sector, so a sector
generated here and one generated by `generate.py sector` directly are built by
the same code path — the only difference is a real galaxy-frame position
(and, in turn, the same sparse exotic-phenomena population every sector
gets — see "Sector-Level Exotic Phenomena" above).

Sectors sit on a cylindrical grid with one standard sector size, 4 parsecs
(about 13 ly): the galaxy is a stack of flat layers 4 pc tall (layer 0
centered on the galactic plane), each cut into rings 4 pc wide around the
galactic axis, and each ring cut into wedge-shaped slots about 4 pc across.
A ring has the same slots on every layer, so sectors line up in vertical
columns. Each layer reaches out only as far as the galaxy still expects at
least one star per sector. `--ring I` generates
one whole layer of a ring (layer 0 unless `--layer J` is given), and adding
`--slot K` generates just that one sector (and, with `--radius-pc R`, its
neighborhood within R parsecs too). `--ring I --slot K --column` generates
that slot through every layer the galaxy reaches at that ring, and
`--ring I --shell` generates the whole ring through every layer, a
cylindrical shell usually thousands of sectors large, so it needs `--limit`
or `--yes`.

Run with neither `--ring` nor `--center-sector` (i.e. no arguments at
all), `generate.py galaxy` picks a uniformly random (by volume), not-yet-
occupied sector address inside the galaxy's planned outline (every layer
out to its stored edge, from the top of the galaxy to the bottom),
generates it, and then generates every not-yet-generated sector within
100 ly of it too, in every direction — a whole small starmap around a
fresh, randomly chosen starting point in one run. `--max-ring` bounds how
far out the random starting address can land (defaults to the galaxy's
own edge from `generate.py plan`), `--radius-pc` overrides the default 100 ly
neighborhood radius, and `--min-start-density` requires the randomly
chosen starting sector's own real density to be at least that many times
local (e.g. `--min-start-density 1.0` for at least as dense as the
galaxy's own real local density) before accepting it, retrying otherwise
— useful for skipping past the galaxy's own vast, sparse outskirts to
start somewhere with more to look at. A high threshold combined with a
large `--max-ring` can take many retries to satisfy, since a
volume-weighted random draw favors the sparser outskirts to begin with.

Each of these has an upper bound, so a typo can't start a run that never
ends: `--radius-pc` at most 200 (about 650 ly), `--ring` and `--max-ring`
at most 100,000, and `--limit` at most the slot count of ring 100,000.
The web Generate page and the API check the same bounds
(`src/stellarObjects/generationLimits.py`).

Most of the galaxy is never actually visited or generated; `generate.py plan`
builds a small, cheap-to-recompute density "skeleton" (one singleton shape
row plus one row per layer, from the top of the galaxy to the bottom, naming
the last ring that layer reaches) that `generate.py galaxy` consults to decide, per
address, whether anything exists there at all before generating it lazily
on demand. It also stores each ring's column bound (the highest and lowest
layer that ring reaches). Every `generate.py galaxy` mode checks its
address against this outline before generating anything, even when
`--density` or `--num-systems` is given, so nothing is ever placed outside
the galaxy: an address outside it is refused with the reason, and a
neighborhood near the edge simply leaves out the sectors past it.
`generate.py galaxy` refuses to run until `generate.py plan` has been run.

After the outline, `generate.py plan` also places every star of 500 solar
luminosities or more across the whole galaxy, before any sector is
filled (about 60 million in a Milky Way, roughly 20 minutes of drawing
plus the database load, about 10 GB of rows). `--bright-star-min-luminosity
100` goes down to 100 solar luminosities instead: about 220 million stars,
roughly an hour and a quarter, and about 35 GB. Each one is a finished star at a fixed point in
its sector, stored in `bright_stars`, so the Galaxy Map can show the
bright stars tracing the spiral arms right away. Filling a sector later
builds a full system around each of its bright stars first, then draws
the rest of its systems from dimmer stars only, so its expected total is
unchanged. Every system's age comes from the stellar population mix where
its sector sits (young stars crowd the arms near the plane; the bulge is
old). `--bright-star-min-luminosity` changes the threshold (100 or more),
`--no-bright-stars` skips the step, and `--bright-stars-only` re-scatters
on the stored outline. The scatter refuses a galaxy whose sectors are
already filled (they would never get their bright stars) unless `--force`
is given, which leaves those sectors out.

## Exotic Phenomena Generation

`generate.py phenomenon` generates a single exotic stellar phenomenon on demand,
independent of any one sector (see "Sector-Level Exotic Phenomena" above
for the population every generated sector gets automatically):

```bash
python generate.py phenomenon --type {black-hole,neutron-star,nebula,supernova-remnant,rogue-planet,comet,asteroid-field,quasar} [options]
```

Omitting `--type` picks uniformly at random among the first seven; a quasar is only made when asked for by name, and `--sector-id` only accepts a ring-0, layer-0 (galactic core) sector for one, placing it at the galactic center. `--anchor-system`
(black hole/neutron star only) builds a full star system around the
compact remnant instead of describing it standalone — real pulsar planets
exist (PSR B1257+12) — reusing all of `generate.py system`'s own orbit-placement
logic; since a compact remnant's near-zero luminosity naturally collapses
the disk-physics planet-count estimate toward zero (matching the real
rarity of confirmed planets around black holes/neutron stars), pass
`--num-orbits` to force orbiting bodies. `--sector-id` (any type except an
anchored one) links the generated phenomenon to an already galaxy-placed
sector — for nebula/asteroid-field/black-hole/neutron-star, it also
computes a real galaxy-frame position near that sector, instead of
leaving it unplaced; other types are linked only (no placement columns of
their own) — see `docs/database-schema.md`'s "v18"/"v21" notes.
`--markdown`, `--debug`/`--quiet`/`--silent`, and the
`--mysql-*` connection options all work the same way they do on
`generate.py system`.

## Other Tools

- **`src/queryDb.py`** — a read-only search CLI over an existing database:
  `sectors` (list every sector), `systems --star-type ... --sector-id ...`
  (filtered system list), and `near <system_id> --radius <ly>` (systems
  within a radius of another, same sector).
- **`src/updateOrbits.py`** — advances every planet's/moon's live orbital
  position based on real elapsed time since the last run, including each
  binary pair's own barycentric wobble (both stars now visibly orbit
  their common center, not one fixed while the other circles it) and a
  planet's/star's own small "reflex offset" from the combined pull of its
  moons/planets; meant to be run periodically (e.g. via cron, "once a
  month or so"), not on every generation run. See [`docs/database-schema.md`](docs/database-schema.md#orbit_simulation_state)
  for a cron example, or [`examples/maintenance/`](examples/maintenance/)
  for a systemd timer to run it the Ubuntu/Debian-native way.
- **`src/migrateDb.py`** — brings an existing database's schema up to the
  version this checkout expects, applying any migration steps in between
  (a no-op if it's already current). Run automatically by
  `install.sh`/`update.sh` on every deploy.
- **`src/migrateSqliteToMysql.py`** — one-time import of a pre-MySQL-port
  SQLite database into MySQL.
- **`src/resetDb.py`** — wipes every generated sector/system/galaxy row
  (`TRUNCATE`, not `DROP`) so the database is empty and ready for a fresh
  galaxy, leaving the schema itself and the separate control (admin) schema
  untouched. Destructive and unrecoverable — prompts for the database name
  to be typed back before doing anything (`--yes` skips this, for
  scripted use only); `--dry-run` lists what would be wiped without
  touching it.

## Web Interface

One Flask app serves a JSON API ([`src/html/api/`](docs/api.md), needs
`pymysql`/`DBUtils` and MySQL access) and, on top of it, the web pages
([`src/html/web/`](docs/html-interface.md)), which fetch everything they
show from that API in-process. Drill into sectors and star systems, view
sector/system/galaxy visualizations (an interactive 3D "Sector Map", a
scaled "System Map" orbit diagram, and a galaxy-scale "Galaxy Map"), copy
the rendered wikitext/Markdown page saved for each system, or jump
straight to an object via the `/search` page's faceted search: click-to-filter
tag buttons (object type, star spectral/luminosity class, planet
class/body type, supported life chemistry) built from only the values
actually present in the chosen database, plus an autocompleting name
search across sectors, systems, stars, and planets/moons -- a search box
is always in the page header. The site shows the one database configured
in `config.json`.

**Deploying:** [`docs/deployment/`](docs/deployment/README.md) compares the
supported platforms and has a guide and example configs for each: Apache2
with mod_wsgi (the reference setup), nginx or Caddy with gunicorn on
Linux, IIS, Caddy or Apache with waitress on Windows (plus WSL2), and
Homebrew nginx with gunicorn on macOS. On Debian or Ubuntu, `sudo
./install.sh` does the whole install except the web server's site file,
and `sudo ./update.sh` pulls and applies updates (plain `git pull` isn't
enough on its own). Every deployment-level setting
(MySQL connection details, rate limits, site name/base URL, and more) can
be set once in a `config.json` file at the repo root, instead of (or
alongside) the `PLANETGEN_*` environment variables, kept outside the
served `src/html/` tree -- see [`config.md`](docs/config.md). See
[`docs/html-interface.md`](docs/html-interface.md) for how the interface works and how to
deploy or test it locally.

### Additional Information

This tool is designed as a personal tool for the Molten Aether FFRP game. Everything it generates is saved to the database in both wikitext and Markdown, designed to be simply cut and paste from the web interface into the wiki. See https://wiki.moltenaether.com for wiki and game information.

Using commands to force a habitable world and an asteroid belt will automatically force a large star to ensure there is room for both objects. The most common stars are small dwarf stars which make a smaller star system. Forcing a large star as well as the maximum number of planets will cause the generated system to be very large with an extremely high number of planets. Do not assume that just because it is generated here, it is accurate or possible, such large systems may require editing as some worlds may end up saying they are several thousand AU's from the central star.

## Testing

`pip install -e ".[test,api]"` and then `pytest` runs the whole suite,
including the brute-force (Hypothesis) tests that try to break generation,
placement, the CLI and every web page. Database tests need a MySQL server
and skip without one. See [`docs/testing.md`](docs/testing.md) for the
database settings, the fuzz profiles and the weekly deep fuzz run.

## License

This project is licensed under the [CC0 1.0 Universal](LICENSE.md) license.