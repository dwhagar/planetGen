# Command-line reference

Every option of `generate.py`, the one entry point for every generator.
For what planetGen is and a short tour, see the [README](../README.md);
for installing it, [INSTALL.md](../INSTALL.md). On Windows run it with
the venv's Python (`C:\srv\planetgen-venv\Scripts\python.exe`); on
Linux and macOS the installer also adds a `planetgen` command that runs
it. `python3 generate.py <command> --help` prints the same options.
## Subcommands

`generate.py` is the single command-line entry point for every generator in this project — one script, with a subcommand per generation scale:

```bash
python3 generate.py system [options]      # one star system
python3 generate.py sector [options]      # one or more independent sectors
python3 generate.py galaxy [options]      # many sectors placed as one galaxy
python3 generate.py plan [options]        # the galaxy's density skeleton
python3 generate.py phenomenon [options]  # one exotic stellar phenomenon
```

Run `python3 generate.py <command> --help` for that command's own full option list. Every subcommand saves what it generates to the database.

To generate a new star system:

```bash
python3 generate.py system [options]
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
*   `--debug [file]`: Logs every choice the generator makes, and why, with timestamps, to the console. If `file` is given, also mirrors that output to `file`. Available on every subcommand. For a permanent, far more detailed log of everything (every random roll, SQL statement and web request too), set `"debug": true` in `config.json` instead; it writes to `/var/log/planetgen.log` (see [`docs/config.md`](config.md)).
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
python3 generate.py sector [options]
```

Most of `system`'s options work here too, but apply *uniformly* to every system in the sector — `+asteroid_belt` guarantees a belt in every system, `--star-type G2V` makes every star in the sector a G2V, and so on. `--system-file`, `--num-orbits`, and `system`'s own per-system `--name` aren't offered here, since those describe one specific, hand-crafted system rather than a sector of varied ones; use `generate.py system --system-file` directly for that.

Sector-specific options:

*   `--num-systems <int>`: How many star systems the sector contains. Defaults to 10.
*   `--name <name>`, `-n <name>`: Hard-sets the sector's own name, overriding the default random two-word name (e.g. `"Voranthis Kelmoor"`) generated the same phoneme-salad way as star names.
*   `--min-habitable <int>`: Guarantees at least this many systems in the sector have a habitable world, chosen randomly among them — without forcing *every* system to have one the way a uniform `+habitable_world` would. Extra systems can still turn out habitable by chance on top of this minimum. Cannot exceed `--num-systems`, and cannot be combined with a uniform `-habitable_world`.
*   `--mysql-host <host>`, `--mysql-port <port>`, `--mysql-user <user>`, `--mysql-password <password>`, `--mysql-database <database>`: Where the generated sector is saved. Each defaults to the matching `$PLANETGEN_MYSQL_*` environment variable, or a built-in default (`127.0.0.1:3306`, user/database `planetgen`) -- see [`database-schema.md`](database-schema.md).

Each run saves the whole generated sector — every system, star, planet, moon, and asteroid belt, plus any exotic phenomena — to the MySQL database described in [`database-schema.md`](database-schema.md), printing only a short status line and summary per saved sector to the console (`--debug` shows more). The wiki page is rendered from those rows whenever it's viewed, not stored.

### Sector-Level Exotic Phenomena

Every generated sector also seeds a realistically sparse population of exotic phenomena: black holes, neutron stars, nebulae (including molecular clouds and planetary nebulae), supernova remnants, rogue planets and free-floating brown dwarfs, and interstellar comets. Each type is drawn per star from the research densities in [`design/interstellar-object-rates.md`](design/interstellar-object-rates.md) (`program_constants.PHENOMENON_DENSITY_PC3`, with a `PHENOMENON_RATE_SCALE` dial per type), so most sectors hold none of the rarer kinds, as real space this size usually doesn't. Isolated asteroid fields are not generated in sectors (they disperse); `generate.py phenomenon --type asteroid-field` still makes one on request. Nebulae come with the stars that make them: every O star and half the B0-B2 stars sit in an H II region, and every planetary nebula has its own new hot white dwarf at the center. A galaxy's nucleus is active 10% of the time (`program_constants.QUASAR_ACTIVE_NUCLEUS_CHANCE`): then the core sector at ring 0, layer 0, slot 0 gets a quasar at the galactic center; otherwise it gets a quiescent supermassive black hole.

Black holes and neutron stars are real, stellar-mass gravitating bodies, so they're placed the same Hill-sphere-aware way every star system itself already is: never within a neighboring star system's or another compact remnant's own Hill sphere (see "Minimum separation (Hill spheres)" in `stellarObjects/spaceSector.py`'s own module docstring). Every other phenomenon type has no comparable gravitational footprint at this scale and is placed at a random point in the sector's cell instead.

## Galaxy Generation

`generate.py galaxy` generates and persists many sectors as one galaxy, placed
in real galaxy-frame 3D space (see
[`docs/design/galaxy-coordinate-system.md`](design/galaxy-coordinate-system.md)):

```bash
python3 generate.py galaxy --ring I [--layer J] [--slot K [--radius-pc R]] [options]
python3 generate.py galaxy --ring I --slot K --column [options]
python3 generate.py galaxy --ring I --shell [--limit N | --yes] [options]
python3 generate.py galaxy --block M.I.S.SLAB [--block-layer J] [options]
python3 generate.py galaxy --center-sector ID --radius-pc R [options]
python3 generate.py galaxy [options]
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
`--block M.I.S.SLAB` generates one Galaxy Map drill-down block (size M of
243, 27 or 3, in ring I, wedge S and slab SLAB; see
[`design/galaxy-drilldown-navigation.md`](design/galaxy-drilldown-navigation.md)),
and `--block-layer J` narrows it to one of the block's layers. Blocks past
the large-ring threshold need `--limit` or `--yes` too.

Run with neither `--ring` nor `--center-sector` (i.e. no arguments at
all), `generate.py galaxy` picks a uniformly random (by volume), not-yet-
occupied sector address inside the galaxy's planned outline (every layer
out to its stored edge, from the top of the galaxy to the bottom),
generates it, and then generates every not-yet-generated sector within
12 pc (about 39 ly) of it too, in every direction — a whole small starmap
around a fresh, randomly chosen starting point in one run. Every sector
any `galaxy` mode generates (and every sector the map generates on a
visit) first gets the bright stars around it: each sector block (3x3x3
sectors) within 100 ly is given every star from 100 L_sun up to what was
already placed there, once per block, leaving filled sectors alone
(GEN.23). `--max-ring` bounds how
far out the random starting address can land (defaults to the galaxy's
own edge from `generate.py plan`), `--radius-pc` overrides the default 12 pc
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

### Parallel generation

`generate.py sector` (with `--num-sectors`) and every `generate.py galaxy`
mode fill several sectors at once, each in its own worker process: by
default 80% of the machine's cores, one fewer when MySQL runs on the same
machine (3 workers on 4 cores without a local MySQL, 2 with one). Workers
run at a lower priority (`nice` 10 on Linux and macOS, below normal on
Windows), so the web site and anything else on the machine come first.
Each worker generates a whole sector and saves it in one transaction;
the run prints each sector as it's saved, so sectors can finish out of
order. `--workers N` sets the count (`PLANETGEN_WORKERS` does the same for
every run), and `--workers 1` generates one sector at a time in the run's
own process, as before. On a 4-core machine with MySQL local, 60 sectors
of ring 2000 took 16.5 s on one worker and 9.8 s on the default two.

Only one run's workers use the machine at a time: a run started while
another is generating (from the command line or the Generate page) says
it's waiting and starts when the other finishes. The control database
keeps the lease and one row per run and per sector (`work_jobs`,
`work_tasks`, see [`database-schema.md`](database-schema.md)); without
it (before `update.sh` has created the tables) runs don't wait for each
other. A run stopped with Ctrl-C, or Cancel on the Generate page, lets
the sectors already being saved finish and queues nothing more; a run
that was killed outright frees the lease after 30 seconds.

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
python3 generate.py phenomenon --type {black-hole,neutron-star,nebula,supernova-remnant,rogue-planet,comet,asteroid-field,quasar} [options]
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
their own) — see `database-schema.md`'s "v18"/"v21" notes.
`--markdown`, `--debug`/`--quiet`/`--silent`, and the
`--mysql-*` connection options all work the same way they do on
`generate.py system`.

## Population

`generate.py population` is an optional pass over what is already
stored: it names the dominant species of every world with complex life,
gives technological civilizations an age and an era (Industrial through
Elder), founds one polity per spacefaring species, and works out which
systems each polity owns (reach up to 100 ly from its capital). It is
off by default and never runs on its own.

```bash
python3 generate.py population                     # scan what isn't scanned yet
python3 generate.py population --rescan            # forget everything and start over (new names, ages and borders)
python3 generate.py population --territories-only  # only recompute who owns which system
```

`generate.py sector` and `generate.py galaxy` run it after saving when
given `--population`. The install and update scripts offer it (y/N,
default N after 30 seconds; `POPULATION=1` or `-Population` runs it
without asking). The website hides species, polities and the Galaxy
Map's Territories overlay until it has made something to show
(`GET /api/population`). The `--mysql-*`, `--debug` and `--quiet`
options work as on the other subcommands. Design and decisions:
[`design/population-and-politics.md`](design/population-and-politics.md).
