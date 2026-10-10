# Command-line reference

Every option of `planetgen`, the one entry point for every generator.
For what planetGen is and a short tour, see the [README](../README.md);
for installing it, [INSTALL.md](../INSTALL.md). The installer adds a
`planetgen` command that runs it. `planetgen <command> --help` prints the same options.
## Subcommands

`planetgen` is the single command-line entry point for every generator in this project — one command, with a subcommand per generation scale:

```bash
planetgen system [options]      # one star system
planetgen sector [options]      # one or more independent sectors
planetgen galaxy [options]      # many sectors placed as one galaxy
planetgen plan [options]        # the galaxy's density skeleton
planetgen phenomenon [options]  # one exotic stellar phenomenon
planetgen population [options]  # species, civilizations and territories
planetgen check-math [-v]       # the math check bulk runs start with
planetgen fingerprint [options] # a digest of the generated content, to compare two builds
```

Run `planetgen <command> --help` for that command's own full option list. Every subcommand except `check-math` and `fingerprint` saves what it generates to the database.

**The math check comes first.** `check-math` runs `planetgen/physics/mathcheck.py` (known answers from real astronomy, identities, and sampler distributions; see [testing.md](testing.md#the-math-check-runs-first)) and exits 1 if a check fails; `-v` lists every check. Every bulk run (`galaxy`, `plan`, `population`, and `sector --num-sectors` above 1) runs it first and refuses to start if a check fails, naming the failed checks and writing nothing, not even the activity log line. So do the Generate page's jobs (their first step, "Check the math") and the Sector page's "generate the neighbourhood" button. One system, one sector or one phenomenon is not gated. `update.sh` runs it after updating and warn (and skip the population pass) if it fails.

To generate a new star system:

```bash
planetgen system [options]
```

To generate a whole sector of star systems at once, see [Sector Generation](#sector-generation) below.

### Options

Most `system` options use a `+name`/`-name` tri-state syntax (these forcing options are for a single system only; `sector` and `galaxy` refuse them): `+name` forces that feature to be present, `-name` forces it to be absent, and leaving it off the command line leaves it up to chance.

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
*   `--name <name>`: Specifies a name for the star system, overriding the default random generation. Everything in the system is named from it: planets are numbered in orbit order (`Sol I`, `Sol II`, ...), moons add a letter (`Sol IIIa`), a close binary's two stars are the system's A and B (`Sol A`, `Sol B`, with planets named for the system), and a wide binary's two stars are two words sharing the system name's first word (`Sol Kelmoor`, `Sol Pikkita`; the second star's word sounds like a diminutive). A wide pair's primary numbers its planets after that shared word (`Sol I`) and its secondary after its own word (`Pikkita I`). Asteroid belts aren't numbered.
*   `--age <young|old>`: Specifies the age of the star system (young or old).
*   `--flavor-chance-system <float>`: Overrides the default system-level flavor text chance (0.0 to 1.0).
*   `--flavor-chance-planet <float>`: Overrides the default planet-level flavor text chance (0.0 to 1.0).
*   `--max-planet-flavor`: Sets the maximum flavor text total for planets to 99.
*   `--debug [file]`: Logs every choice the generator makes, and why, with timestamps, to the console. If `file` is given, also mirrors that output to `file`. Available on every subcommand. For a permanent, far more detailed log of everything (every random roll, SQL statement and web request too), set `"debug": true` in `config.json` instead; it writes to `/var/log/planetgen.log` (see [`docs/config.md`](config.md)).
*   `--quiet` / `--silent`: Suppresses all output except errors. Available on every subcommand. Combined with `--debug`, a filename is required (there would otherwise be nowhere for debug output to go).

**Note on Incompatible Options:**

*   `-planets` cannot be combined with `+moons`, `+max_planets`, `+habitable_world`, or `+asteroid_belt` (on the `system` command this also covers the same keys, `num_orbits` and `slots`, in a `--system-file`).
*   `--star-type` cannot be combined with `+large_star`.
*   `+intelligent_life`/`-intelligent_life` cannot be combined with `-habitable_world`.
*   `+habitable_world` and `+asteroid_belt` together cannot be combined with `-large_star` (both objects require the room a large star provides).
*   `--num-orbits` cannot be combined with `-planets`.
*   `--num-orbits` must be from 0 to 500 (the generator's own ceiling on objects in a system).
*   `--flavor-chance-system` must be a float between 0.0 and 1.0.
*   `--flavor-chance-planet` must be a float between 0.0 and 1.0.

A forced body some stars can't host (a habitable world around a hot O or B
star, for example) is tried on up to five whole systems. If none of them has
it, `system` exits with an error and saves nothing.

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

`planetgen sector` generates a whole sector of independently-random star systems in one pass, reusing `planetgen system`'s own generation logic for each one:

```bash
planetgen sector [options]
```

Some of `system`'s options work here too, and apply *uniformly* to every system in the sector: `--star-type G2V` makes every star in the sector a G2V, and so on. The `+name`/`-name` forcing options are for a single system only: `sector` and `galaxy` refuse them with an error naming the option (a saved or queued command line that still has one gets the same message). `--system-file`, `--num-orbits`, and `system`'s own per-system `--name` aren't offered here either, since those describe one specific, hand-crafted system rather than a sector of varied ones; use `planetgen system` for that.

Sector-specific options:

*   `--num-systems <int>`: How many star systems the sector contains. Defaults to 10.
*   `--name <name>`, `-n <name>`: Hard-sets the sector's own name, overriding the default random two-word name (e.g. `"Voranthis Kelmoor"`) generated the same phoneme-salad way as star names.
*   `--min-habitable <int>`: Guarantees at least this many systems in the sector have a habitable world, chosen randomly among them — without forcing *every* system to have one. Extra systems can still turn out habitable by chance on top of this minimum. Cannot exceed `--num-systems`.
*   `--prevalence <feature>=<percent>`: Makes a feature more or less common across every system the run generates, as a percentage of its normal chance: `--prevalence comets=+50` gives 1.5 times as many systems with comets, `--prevalence binary_system=-100` none at all. Repeat it for more features (also on `galaxy`). The features are the forcing options' names: `habitable_world`, `asteroid_belt`, `comets`, `large_star`, `moons` (each planet's chance of moons), `max_planets` (the chance a star gets the most orbits it can hold), `intelligent_life` (the population pass's chance of a civilization), `binary_system`, `wide_binary` and `planets`. Each system stores the setting with its config.
*   `--mysql-host <host>`, `--mysql-port <port>`, `--mysql-user <user>`, `--mysql-password <password>`, `--mysql-database <database>`: Where the generated sector is saved. Each defaults to the matching `$PLANETGEN_MYSQL_*` environment variable, or a built-in default (`127.0.0.1:3306`, user/database `planetgen`) -- see [`database-schema.md`](database-schema.md).

Each run saves the whole generated sector — every system, star, planet, moon, and asteroid belt, plus any exotic phenomena — to the MySQL database described in [`database-schema.md`](database-schema.md), printing only a short status line and summary per saved sector to the console (`--debug` shows more). The wiki page is rendered from those rows whenever it's viewed, not stored.

### Sector-Level Exotic Phenomena

Every generated sector also seeds a realistically sparse population of exotic phenomena: black holes, neutron stars, nebulae (including molecular clouds and planetary nebulae), supernova remnants, rogue planets and free-floating brown dwarfs, and interstellar comets. Each type is drawn per star from the research densities in [`design/interstellar-object-rates.md`](design/interstellar-object-rates.md) (`program_constants.PHENOMENON_DENSITY_PC3`, with a `PHENOMENON_RATE_SCALE` dial per type), so most sectors hold none of the rarer kinds, as real space this size usually doesn't. Isolated asteroid fields are not generated in sectors (they disperse); `planetgen phenomenon --type asteroid-field` still makes one on request. In a galaxy, molecular clouds (classes M-Q) belong to the galaxy rather than to one sector: they are drawn per 50 pc cell from the galaxy seed, more on the spiral arms and near the plane, and a sector gets every cloud that reaches it, stored once whichever sector it reaches is generated first (GEN.47). On an arm at the Sun's distance about one sector in eight sits in one, between the arms about one in twenty-five. Nebulae come with the stars that make them: every O star and half the B0-B2 stars sit in an H II region, and every planetary nebula has its own new hot white dwarf at the center. A galaxy's nucleus is active 10% of the time (`program_constants.QUASAR_ACTIVE_NUCLEUS_CHANCE`): then the core sector at ring 0, layer 0, slot 0 gets a quasar at the galactic center; otherwise it gets a quiescent supermassive black hole.

Black holes and neutron stars are real, stellar-mass gravitating bodies, so they're placed the same Hill-sphere-aware way every star system itself already is: never within a neighboring star system's or another compact remnant's own Hill sphere (see "Minimum separation (Hill spheres)" in `planetgen/galaxy/sector.py`'s own module docstring). Every other phenomenon type has no comparable gravitational footprint at this scale and is placed at a random point in the sector's cell instead.

## Galaxy Generation

`planetgen galaxy` generates and persists many sectors as one galaxy, placed
in real galaxy-frame 3D space (see
[`docs/design/galaxy-coordinate-system.md`](design/galaxy-coordinate-system.md)):

```bash
planetgen galaxy --ring I [--layer J] [--slot K [--radius-pc R]] [options]
planetgen galaxy --ring I --slot K --column [options]
planetgen galaxy --ring I --shell [--limit N | --yes] [options]
planetgen galaxy --block M.I.S.SLAB [--block-layer J] [options]
planetgen galaxy --center-sector ID --radius-pc R [options]
planetgen galaxy [options]
```

It reuses `planetgen sector`'s own generation logic per sector, so a sector
generated here and one generated by `planetgen sector` directly are built by
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
all), `planetgen galaxy` picks a uniformly random (by volume), not-yet-
occupied sector address inside the galaxy's planned outline (every layer
out to its stored edge, from the top of the galaxy to the bottom),
generates it, and then generates every not-yet-generated sector within
12 pc (about 39 ly) of it too, in every direction — a whole small starmap
around a fresh, randomly chosen starting point in one run. Every sector
any `galaxy` mode generates (and every sector the map generates on a
visit) first gets the bright stars around it: each sector within 100 ly
is given every star from its distance tier's floor up to what was already
placed there, once per sector, leaving filled sectors alone (GEN.23,
GEN.44). `--max-ring` bounds how
far out the random starting address can land (defaults to the galaxy's
own edge from `planetgen plan`), `--radius-pc` overrides the default 12 pc
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
(`src/planetgen/generation/limits.py`).

### Parallel generation

`planetgen sector` (with `--num-sectors`) and every `planetgen galaxy`
mode fill several sectors at once, each in its own worker process: by
default 80% of the machine's cores, one fewer when MySQL runs on the same
machine (3 workers on 4 cores without a local MySQL, 2 with one). Workers
run at a lower priority (`nice` 10), so the web site and anything else on the machine come first.
Each worker generates a whole sector and saves it in one transaction;
the run prints each sector as it's saved, so sectors can finish out of
order. `--workers N` sets the count (`PLANETGEN_WORKERS` does the same for
every run), and `--workers 1` generates one sector at a time in the run's
own process, as before. The workers are RQ workers on the Redis server
at `redis.url` ([config.md](config.md)), started by the run for its own
tasks and gone when it finishes. Without a Redis server the run says so
and generates one sector at a time. A task whose worker dies (killed, out
of memory) runs once more on a fresh worker before the run fails. On a 4-core machine with MySQL local, 60 sectors
of ring 2000 took 16.5 s on one worker and 9.8 s on the default two.

`planetgen plan` draws its bright stars the same way, one layer of the
galaxy per task, densest layers first (`--workers` works there too).
Each layer draws from its own random stream, so the same seed places the
same stars on any number of workers. That seed is drawn at random and
stored (`galaxy_shape.bright_star_seed`); there is no `--seed` option
yet, and the rest of generation is not reproducible yet (GEN.39, see
[Planned commands](#planned-commands)).

Every progress bar shows the time elapsed and the time remaining until
it's done. The remaining time comes from a decaying average of how many
sectors (or layers) finished per second, weighted toward the last minute
or so, so it follows the run's current speed and doesn't jump about when
several workers finish at once. The Generate page shows the same
estimate ("about 4 m 10 s left").

The bright-star bar ("Bright stars (12 of 1,271 layers)") shows a share
done rather than a count: before the scatter starts, every layer gets an
expected star count from the same density model the scatter draws from
(a quick sample of its rings, a second or two for a whole galaxy), plus
five stars' worth per ring walked, so the nearly empty layers at the top
and bottom of the disk count for little and the dense middle for a lot.
Layers in progress report their stars a few times a second and the bar
moves with them, so the time left follows the work left rather than the
layers left. While layers take longer than 30 seconds each, a second bar
appears under it with the stars of the layers being drawn, done of their
estimate, and its own time left ("Layer 0: stars"); it goes again once
layers finish faster than one every 20 seconds. Sector fill keeps its one
bar. The Generate page shows both lines too.

Only one run's workers use the machine at a time: a run started while
another is generating (from the command line or the Generate page) says
it's waiting and starts when the other finishes. The control database
keeps the lease and one row per run and per sector (`work_jobs`,
`work_tasks`, see [`database-schema.md`](database-schema.md)); without
it (before `update.sh` has created the tables) runs don't wait for each
other. A run stopped with Ctrl-C, or Cancel on the Generate page, stops
its workers mid-sector (each sector is saved in one transaction, so an
unfinished one is rolled back, never half written), queues nothing more
and is recorded as cancelled, with the lease freed at once; a run that
was killed outright frees the lease after 30 seconds. A worker that
dies fails the run; losing the control database mid-run doesn't stop
it (only its rows stop being updated).

Every run is also recorded as a job tree (control schema v7): the run
at the top, its phases (the skeleton, the bright stars, the population
pass) under it, each work queue under those, and the queue's sectors or
layers as its leaves, each with its own start, end and duration. A run
the Generate page started hangs under that page job and its step. This
happens with one worker too (without the lease), and is skipped without
the control database.

### Size and time estimates

Before a bulk run writes anything (`planetgen galaxy` in every mode,
`planetgen sector --num-sectors`), it works out how big and how long
it will be and prints it:

```
Estimate for ring 12 layer 0: About 1.6 GB and 10 m 59 s for 4,289 sectors (~23,926 star
systems, ~31,198 stars, 2 workers). 30 GB free of 271 GB.
```

- **Systems** come from each sector's expected density (the galaxy's
  skeleton, or `--density`/`--num-systems`).
- **Size** is those systems times this galaxy's own bytes per system,
  measured from its tables after every run (about 62 KB on a test
  galaxy, over half of it moons), plus 10%.
- **Time** adds up each sector's expected systems times this server's
  measured seconds per system at that density, divided by the workers.
  The speed is kept per density bucket (two per decade from 0.01) as a
  decaying average of every sector ever filled, in the control
  database (`generation_stats`), so dense sectors aren't priced like
  sparse ones. Until a bucket has been measured the nearest measured
  one is used, and before anything has been measured a default of
  0.2 s per system (the estimate says so).

On a terminal, a run of more than one sector then asks `Generate these
N sectors? [y/N]` (`--yes` skips the question). Off a terminal (the
Generate page's jobs, scripts and cron) it prints the estimate and
carries on; the Generate page shows the estimate and asks first.

A run is **refused**, with nothing written, when it would take more than
a quarter of the database disk or leave less than 5 GB free; the
message says how much it needs. The disk is the one holding MySQL's
data directory, asked of the server itself (`SELECT @@datadir`) and never
assumed to be the boot drive; the estimate names the path and the mount
measured, with symlinks and bind mounts resolved. A
server on another machine is measured only if that directory is visible
here or the server reports its own disks (MariaDB); otherwise the disk
shows as "unknown" with the reason and nothing is refused for space. `--yes` never overrides a refusal;
generate fewer sectors (`--limit`, a smaller radius) instead.

`--estimate-only` prints the estimate (and any refusal) and stops
without writing anything, ending with one `ESTIMATE {json}` line. In
random-start mode it is the estimate for one random start, so the real
run's start (drawn again) differs.


Most of the galaxy is never actually visited or generated; `planetgen plan`
builds a small, cheap-to-recompute density "skeleton" (one singleton shape
row plus one row per layer, from the top of the galaxy to the bottom, naming
the last ring that layer reaches) that `planetgen galaxy` consults to decide, per
address, whether anything exists there at all before generating it lazily
on demand. It also stores each ring's column bound (the highest and lowest
layer that ring reaches). Every `planetgen galaxy` mode checks its
address against this outline before generating anything, even when
`--density` or `--num-systems` is given, so nothing is ever placed outside
the galaxy: an address outside it is refused with the reason, and a
neighborhood near the edge simply leaves out the sectors past it.
`planetgen galaxy` refuses to run until `planetgen plan` has been run.

**The galaxy seed.** The first `planetgen plan` stores the galaxy's
128-bit seed, shown as 32 hex digits (`Galaxy seed 3F2A...C901 (drawn at
random)`). `--seed <32 hex digits>` sets it instead; a later plan keeps
the stored seed, and a different `--seed` is refused once any sector
exists. Each sector, the bright-star scatter, each band added with
`--bright-stars-down-to` and each backfill block draws from its own seed,
the SHA-256 of the galaxy seed and its address (`sector:12/3/0`), so one
seed fills a sector the same way at any `--workers` count and in any run
that reaches it, on the same PlanetGen release. Two systems that draw the
same name are still told apart in save order, and the bright-star
backfill still depends on which sectors are filled first (both phase 1,
GEN.57). The galaxy also records the version key of the code that made
it (22 hex digits: release, Python version, OS and architecture), and
every run that changes it adds a row to `generation_runs` (command line,
the run's own seed, version key, start, end and outcome; DB.6). Every
run's first line names the galaxy seed, the release with its key, and
the command line without the `--mysql-*` and `--debug` options, for
example `Galaxy seed 3F2A...C901, PlanetGen 7.127.352
(0007007F000160030C0300), run: galaxy --ring 3` (`none yet` before the
first plan; OPS.10). It goes to the console, the `--debug` file and the
debug log, and a Generate page job's log starts with the same line for
the job before its first step. Design:
[`design/reproducible-galaxies.md`](design/reproducible-galaxies.md).

After the outline, `planetgen plan` also places every star of 500 solar
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
on the stored outline. `--redo-scatters mass luminosity phenomena` (GEN.196) redoes only the scatters named, each with this run's settings: `mass` clears and rewrites the stars born at the mass limit or more (`--phenomenon-min-mass`), `luminosity` the lighter stars at `--bright-star-min-luminosity`, and `phenomena` the neutron stars, black holes and the rest (`--compact-min-mass`). A new mass limit redoes both star passes. The scatter refuses a galaxy whose sectors are
already filled (they would never get their bright stars) unless `--force`
is given, which leaves those sectors out.

The scatter can go down in stages. `--bright-stars-down-to N` keeps every
bright star already placed and adds only those from `N` up to (not
including) the galaxy's star-fill level, the threshold already scattered
(500 by default), then lowers that level to `N`. So a quick test galaxy
scattered at 500 can later go down to 100 without redrawing anything
brighter:

```bash
planetgen plan --bright-stars-down-to 100
```

Sectors already filled get none of the new stars: their own systems were
drawn below the old level, so they already hold stars that bright. Asking
for a level at or above the current one does nothing and says so. The
Generate page shows the current level and offers the same step as "Add a
dimmer layer of bright stars".

## Exotic Phenomena Generation

`planetgen phenomenon` generates a single exotic stellar phenomenon on demand,
independent of any one sector (see "Sector-Level Exotic Phenomena" above
for the population every generated sector gets automatically):

```bash
planetgen phenomenon --type {black-hole,neutron-star,nebula,supernova-remnant,rogue-planet,comet,asteroid-field,quasar} [options]
```

Omitting `--type` picks uniformly at random among the first seven; a quasar is only made when asked for by name, and `--sector-id` only accepts a ring-0, layer-0 (galactic core) sector for one, placing it at the galactic center. `--anchor-system`
(black hole/neutron star only) builds a full star system around the
compact remnant instead of describing it standalone — real pulsar planets
exist (PSR B1257+12) — reusing all of `planetgen system`'s own orbit-placement
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
`planetgen system`.

## Population

`planetgen population` is an optional pass over what is already
stored: it names the dominant species of every world with complex life,
gives technological civilizations an age and an era (Industrial through
Elder), founds one polity per spacefaring species, and works out which
systems each polity owns (reach up to 100 ly from its capital). It is
off by default and never runs on its own.

```bash
planetgen population                     # scan what isn't scanned yet
planetgen population --rescan            # forget everything and start over (new names, ages and borders)
planetgen population --territories-only  # only recompute who owns which system
```

`planetgen sector` and `planetgen galaxy` run it after saving when
given `--population`. The install scripts offer it (y/N,
default N after 30 seconds; `POPULATION=1` or `-Population` runs it
without asking). The website hides species, polities and the Galaxy
Map's Territories overlay until it has made something to show
(`GET /api/population`). The `--mysql-*`, `--debug` and `--quiet`
options work as on the other subcommands. Design and decisions:
[`design/population-and-politics.md`](design/population-and-politics.md).

## Fingerprint

`planetgen fingerprint` prints a canonical SHA-256 digest of each
sector's generated content, one line a sector (`ring layer slot
digest`, unplaced sectors last), then the plan's digest and the
region's. Two builds of one galaxy are the same when their digests
match; a sector whose line differs is where they part. It writes
nothing and isn't logged as a run.

```bash
planetgen fingerprint                      # the whole galaxy, plan included
planetgen fingerprint --ring 3 --ring 4    # every sector on rings 3 and 4
planetgen fingerprint --sector 12 3 0      # one sector (repeatable)
```

It compares what [`design/reproducible-galaxies.md`](design/reproducible-galaxies.md)
section 2 defines: every object's address, position, properties and
name, with their child rows and the cell's bright stars and phenomenon
scatter. Row ids, timestamps, each row's update clock, the nearest
systems and the location text written from them, sector paths and
stats, population data and bookkeeping are left out; a foreign key
counts as what it points at (a sector's address, an object's unique
ID). It reads the galaxy as stored; comparing a rebuild with the
stored galaxy is internal testing only.

## Planned commands

Not built yet. Each names its TODO item and phase; the design is in
[`design/reproducible-galaxies.md`](design/reproducible-galaxies.md).

| Command or option | Item | Phase | What it will do |
|---|---|---|---|
| `planetgen check-db [--sector S] [--region R]` | DB.8 | 0 | Check the galaxy and control databases without changing anything (schema, orphans, ids, names, values, counts) and exit non-zero when damage is found; also a button on the Admin dashboard. |
| Version-key history listing | OPS.13 | 1 | List the last 10 version keys recorded for a galaxy by `update.sh`. |
| `planetgen repair-db` | DB.9 | 1 | Rebuild damaged sectors from the parity file, or regenerate them from their seed when the version key matches, then check again. |
| `--strict` | GEN.81 | 0 | Today's refusals (density, qualify, size, no room) become warnings and the run goes ahead; `--strict` keeps the old stop for scripts. |
| `--resume` | PERF.30 | 1 | Finish the runs an interrupted fill left, from the step each reached. |
| Layer, ring and column ranges; a radial cylinder; N random neighborhoods | ADM.29, ADM.30, GEN.97 | 1 | New fill shapes, also on the Generate page. |
| `--directive` | GEN.96 | 1 | Generation directives for a sector (density, at least N stars of a type, at least N habitable worlds). |

Work runs as RQ jobs on Redis from phase 0 (PERF.24), so `--workers` sets the RQ worker count.

Phase 2 also adds a daily maintenance run, `scripts/maintenance.sh`, set up as a
scheduled job by install and update (OPS.16, OPS.17): the positional
update (`planetgen.cli.orbits`). (The daily settings-file merge and the
18 backup slots, GEN.61 and OPS.18, were dropped on 2026-10-09.)
