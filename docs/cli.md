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
planetgen scatter --KIND        # one scatter on the plan: --mass, --luminosity, --phenomena, --nebula-field, --band
planetgen clear-scatter --KIND  # empty one scatter, so scatter can draw it again: --mass, --luminosity, --phenomena, --nebula-field
planetgen fill-map-luminosity   # fill the map luminosity of the scattered phenomena
planetgen rescore-habitability  # re-derive every stored habitability score
planetgen clean-nebula-field    # count, or with --apply drop, the field clouds the keep rates now reject
planetgen phenomenon [options]  # one exotic stellar phenomenon
planetgen population [options]  # species, civilizations and territories
planetgen check-math [-v]       # the math check bulk runs start with
planetgen fingerprint [options] # a digest of the generated content, to compare two builds
planetgen check-db [options]    # look for damage in the database (writes nothing)
planetgen versions              # the version keys the updates recorded
planetgen benchmark [options]   # time where generation spends its time, in scratch databases
planetgen warm [FAMILY|PHASE ...]     # build the site's caches ahead of the first visit
planetgen commands              # every command the work queue can run
planetgen command-run KIND ARGS # one unit of one of them, alone
```

Run `planetgen <command> --help` for that command's own full option list. `system`, `sector`, `galaxy`, `plan`, `phenomenon` and `population` save what they generate to the database. `check-math`, `fingerprint`, `check-db` and `versions` write nothing, and `benchmark` writes only to scratch databases of its own.

**The math check comes first.** `check-math` runs `planetgen/physics/mathcheck.py` (known answers from real astronomy, identities, and sampler distributions; see [testing.md](testing.md#the-math-check-runs-first)) and exits 1 if a check fails; `-v` / `--verbose` lists every check. Every bulk run (`galaxy`, `plan`, `population`, and `sector --num-sectors` above 1) runs it first and refuses to start if a check fails, naming the failed checks and writing nothing, not even the activity log line. So do the Generate page's jobs (their first step, "Check the math") and the Sector page's "generate the neighbourhood" button. The same runs also need the downloaded data files (second step, "Check the data"; see [datasets](#datasets)). One system, one sector or one phenomenon is not gated. `update.sh` runs it after updating and warn (and skip the population pass) if it fails.

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
*   `+comets` / `-comets`: Force / forbid comets in the system.
*   `+large_star` / `-large_star`: Force / forbid the generation of a large star.
*   `+moons` / `-moons`: Force / forbid moons on the system's planets.
*   `+max_planets` / `-max_planets`: Force the system to the maximum, or the minimum, number of orbital objects it can support.
*   `+intelligent_life` / `-intelligent_life`: Ensure a planet with intelligent life is (or is not) generated. Either implies `+habitable_world`.
*   `+binary_system` / `-binary_system`: Force / forbid a binary star system.
*   `+wide_binary` / `-wide_binary`: Force an S-type (wide) binary, or force a P-type (close) one instead. Only meaningful together with `+binary_system`; if omitted, one of the two is chosen at random.
*   `+planets` / `-planets`: Ensure the system has at least one planet or asteroid belt, or none at all (star only).
*   `--system-file`, `-f <path>`: Load a system generation specification from a JSON file (see below). Any of the options above, or the value options below, given on the command line override the corresponding value from the file.
*   `--num-orbits <int>`: Force an exact number of orbital slots (planets and asteroid belts combined) to be generated.
*   `--format <text|markdown|wikitext>`: How `--output` writes the system: plain text for the terminal (the default), Markdown or MediaWiki text (the same description the site shows; UX.106).
*   `--output`, `-o <file>`: Writes the system's description (in `--format`) to `file` instead of saving the system to the database; `-` writes it to stdout. The admin site's one-off system page (`/admin/generate/system`) does the same from the browser.
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
*   `--density <x>`: Scales the number of systems by a multiplier on the real local stellar density (1.0 is a realistic sector this size, 2.0 twice as dense). The count is drawn from a Poisson distribution, so it still varies from sector to sector. Cannot be combined with `--num-systems`.
*   `--min-habitable <int>`: Guarantees at least this many systems in the sector have a habitable world, chosen randomly among them — without forcing *every* system to have one. Extra systems can still turn out habitable by chance on top of this minimum. Cannot exceed `--num-systems`.
*   `--directive <minimum>`, `--directive-attempts <n>`: Redraws each sector until it holds at least this much: `systems>=N`, `habitable>=N` or `type:X>=N` (X is one of O B A F G K M wd bh ns). Repeat `--directive` for more than one minimum. After `--directive-attempts` draws (200 by default) the closest sector is kept and the run reports which minimum it missed.
*   `--star-mix <single> <close> <wide>`: The share of systems, in percent, with one star, a close binary and a wide pair; the three must total 100. It replaces `--prevalence` for `binary_system` and `wide_binary`, which it cannot be combined with.
*   `--num-sectors <int>`, `--yes`, `--estimate-only`, `--workers <int>`, `--population`, `--strict`: A run of more than one sector works as a `galaxy` run does (see [Size and time estimates](#size-and-time-estimates) and [Parallel generation](#parallel-generation)); `--strict` turns the warnings the console would otherwise go ahead after into errors that stop the run with nothing generated.
*   `--prevalence <feature>=<percent>`: Makes a feature more or less common across every system the run generates, as a percentage of its normal chance: `--prevalence comets=+50` gives 1.5 times as many systems with comets, `--prevalence binary_system=-100` none at all. Repeat it for more features (also on `galaxy`). The features are the forcing options' names: `habitable_world`, `asteroid_belt`, `comets`, `large_star`, `moons` (each planet's chance of moons), `max_planets` (the chance a star gets the most orbits it can hold), `intelligent_life` (the population pass's chance of a civilization), `binary_system`, `wide_binary` and `planets`. Each system stores the setting with its config.
*   `--mysql-host <host>`, `--mysql-port <port>`, `--mysql-user <user>`, `--mysql-password <password>`, `--mysql-database <database>`: Where the generated sector is saved. Each defaults to the matching `$PLANETGEN_MYSQL_*` environment variable, or a built-in default (`127.0.0.1:3306`, user/database `planetgen`) -- see [`database-schema.md`](database-schema.md).

Each run saves the whole generated sector — every system, star, planet, moon, and asteroid belt, plus any exotic phenomena — to the MySQL database described in [`database-schema.md`](database-schema.md), printing only a short status line and summary per saved sector to the console (`--debug` shows more). The wiki page is rendered from those rows whenever it's viewed, not stored.

### Sector-Level Exotic Phenomena

Every generated sector also seeds a realistically sparse population of exotic phenomena: black holes, neutron stars, nebulae (including molecular clouds and planetary nebulae), supernova remnants, rogue planets and free-floating brown dwarfs, and interstellar comets. Each type is drawn per star from the research densities in [`design/interstellar-object-rates.md`](design/interstellar-object-rates.md) (`program_constants.PHENOMENON_DENSITY_PC3`, with a `PHENOMENON_RATE_SCALE` dial per type), so most sectors hold none of the rarer kinds, as real space this size usually doesn't. Isolated asteroid fields are not generated in sectors (they disperse); `planetgen phenomenon --type asteroid-field` still makes one on request. In a galaxy, molecular clouds (classes M-Q) belong to the galaxy rather than to one sector: they are drawn per 50 pc cell from the galaxy seed, more on the spiral arms and near the plane, and a sector gets every cloud that reaches it, stored once whichever sector it reaches is generated first (GEN.47). On an arm at the Sun's distance about one sector in eight sits in one, between the arms about one in twenty-five. Nebulae come with the stars that make them: every O star and half the B0-B2 stars sit in an H II region, and every planetary nebula has its own new hot white dwarf at the center. A galaxy's nucleus is active 10% of the time (`program_constants.QUASAR_ACTIVE_NUCLEUS_CHANCE`): then the core ball (shell 0) gets a quasar at the galactic center; otherwise it gets a quiescent supermassive black hole.

Black holes and neutron stars are real, stellar-mass gravitating bodies, so they're placed the same Hill-sphere-aware way every star system itself already is: never within a neighboring star system's or another compact remnant's own Hill sphere (see "Minimum separation (Hill spheres)" in `planetgen/galaxy/sector.py`'s own module docstring). Every other phenomenon type has no comparable gravitational footprint at this scale and is placed at a random point in the sector's cell instead.

## Galaxy Generation

`planetgen galaxy` generates and persists many sectors as one galaxy, placed
in real galaxy-frame 3D space (see
[`docs/design/spherical-sector-coordinates.md`](design/spherical-sector-coordinates.md)):

```bash
planetgen galaxy --shell I [--band J [--wedge K [--radius-pc R]]] [--limit N | --yes] [options]
planetgen galaxy [--shells A:B] [--bands C:D] [--wedges E:F] [--limit N | --yes] [options]
planetgen galaxy --block M.I.S.SLAB [--limit N | --yes] [options]   # refused for now (see below)
planetgen galaxy --center-sector ID --radius-pc R [options]
planetgen galaxy --course FROM TO [--course-border] [--yes] [options]
planetgen galaxy [--radius-pc R] [--max-shell I] [--neighborhoods N] [options]   # random start
```

It reuses `planetgen sector`'s own generation logic per sector, so a sector
generated here and one generated by `planetgen sector` directly are built by
the same code path — the only difference is a real galaxy-frame position
(and, in turn, the same sparse exotic-phenomena population every sector
gets — see "Sector-Level Exotic Phenomena" above).

Sectors sit on a spherical grid with one standard sector size, 4 parsecs
(about 13 ly). The core is a ball, shell 0; every other shell is a
hollow sphere 4 pc deep around it, numbered outward from 1. A shell is
cut into **bands** (horizontal strips of latitude, signed: band 0 is the
equatorial band that holds the galactic plane, positive bands are north,
negative south, and the two end bands are the circular polar caps), and each band is
cut into **wedges** (vertical slices from pole to pole). A sector is one
(shell, band, wedge) address. Every band of every shell has 6 times m wedges,
where m has no prime factor above 7 and is picked so each sector is close to
4 pc across. The grid covers all of space out to its reach (262 kpc by
default); the galaxy's outline (below) says where the disk and bulge expect
a system, which is where random starts, the maps' planned sectors and a
neighborhood around an outline sector look: every shell out to the last one
whose equatorial band could still hold a star, and in each shell the bands
the galaxy still expects a star in. A sector that is asked for outside it is
generated from the profile of its place (`galaxy/profile.py`).

`--shell I` generates every not-yet-generated sector in shell I that the
outline holds; adding `--band J` narrows it to one band, and adding
`--wedge K` too generates exactly that one sector (and, with
`--radius-pc R`, its neighborhood within R parsecs as well). A single
address is generated wherever the grid's reach holds it, inside the outline or
out in the halo (its density comes from the sector's profile, so a halo sector
is usually saved empty); one beyond the reach, or in a band its shell does not
hold, stops with exit status 1, and one that was already generated reports its
sector ID. A range named in a span (`--shells`, `--bands`) is taken as named
anywhere in the reach; a span with no range means the outline's cells. A
`--center-sector` neighbourhood around a halo sector takes every cell of the
reach within the radius and says how many sectors that is, and how many
systems they expect, before it runs. A whole shell is usually
thousands of sectors, so past 2,000 sectors a run needs `--limit` (generate
only the first N) or `--yes`.
The span form picks sectors by range instead: `--shells FIRST:LAST`,
`--bands FIRST:LAST` (negative south of the equator, so write
`--bands=-2:2`) and `--wedges FIRST:LAST` (only with exactly one shell and
one band; a range whose first is past its last wraps through wedge 0, as in
`50:5`). Each of the three defaults to everything the outline holds, so
`--bands` alone fills those bands out to the galaxy's edge.
`--block M.I.S.SLAB` was the old drill-down's block key (size M, ring I, wedge
S, slab SLAB) and is refused: the Galaxy Map's blocks are now keyed
`m.shell.band.wedge` (`docs/design/spherical-block-ladder.md`), and the
command does not take those yet. Fill a shell, a band or a span instead.
`--center-sector ID --radius-pc R` generates every not-yet-generated sector
within R parsecs of an already placed sector's own center, the sector at the
center first and the rest outward (GEN.101): the sectors are laid on a lattice
of one sector edge around the center, cut into blocks of four sectors a side,
and the blocks come nearest first, in Hilbert order among blocks the same
distance away and inside a block, so the same neighbourhood always fills in the
same order. A ball cannot always be walked step by step between neighbouring
sectors, so the run says how many of its steps are not (a step of more than two
and a half sector edges), with their mean and longest length: about 97% of
them are between neighbours in a 50 pc ball.
`--course FROM TO` (NAV.48) generates the uncharted sectors that block the NAV
course between two objects (a system's object ID, or `<kind>:<ID>` for a body
or phenomenon): the cells of the route's hops through unknown space, inside the
galaxy's outline, not every cell on the straight line. `--course-border` adds
the uncharted cells that share a face with them. Past 5,000 sectors it needs
`--yes`; the size and disk checks are the same as for any run. The NAV page
(`/nav/chart`, admins) starts it as a Generate page job.

Run with none of `--shell`, `--shells`, `--bands`, `--block`, `--course` and
`--center-sector` (random-start mode, no arguments at all), `planetgen galaxy` picks a uniformly random (by volume),
not-yet-occupied sector address inside the galaxy's planned outline,
generates it, and then generates every not-yet-generated sector within
12 pc (about 39 ly) of it too, in every direction — a whole small starmap
around a fresh, randomly chosen starting point in one run. Every sector
any `galaxy` mode generates (and every sector the map generates on a
visit) first gets the bright stars around it: each sector within 100 ly
is given every star from its distance tier's floor up to what was already
placed there, once per sector, leaving filled sectors alone (GEN.23,
GEN.44). `--max-shell` bounds how
far out the random starting address can land (defaults to anywhere inside the
galaxy's stored outline), `--radius-pc` overrides the default 12 pc
neighborhood radius, and `--min-start-density` requires the randomly
chosen starting sector's own real density to be at least that many times
local (e.g. `--min-start-density 1.0` for at least as dense as the
galaxy's own real local density) before accepting it, retrying otherwise
— useful for skipping past the galaxy's own vast, sparse outskirts to
start somewhere with more to look at. A high threshold combined with a
large `--max-shell` can take many retries to satisfy, since a
volume-weighted random draw favors the sparser outskirts to begin with.
`--neighborhoods N` generates N random neighborhoods instead of one; each
start lies inside the galaxy with its whole neighborhood and at least twice the
radius from every other start, so they never overlap. `--neighborhood-gamma G`
biases the starts toward dense space (a start is kept with probability
min(1, relative density)^G; the default 0 keeps every one, so starts are
uniform by volume), and `--avoid-filled-space` drops any start whose
neighborhood touches a sector that is already generated.

Each of these has an upper bound, so a typo can't start a run that never
ends: `--radius-pc` at most 200 (about 650 ly), `--shell`, `--shells` and
`--max-shell` at most 100,000, at most 100 `--neighborhoods`, and `--limit` at most the
number of sectors in shell 100,000.
The web Generate page and the API check the same bounds
(`src/planetgen/generation/limits.py`).

After the sectors, a run saves the path every star system, rogue planet and
comet takes through its sector, for the sectors it created and the ones around
them; `--no-settle` skips that step (`planetgen orbits` saves the paths later).
Once every sector of the run is generated, the bright stars around the run's
edge are backfilled: four steps outward around the generated sectors (a face
apart, no diagonals) get every star born from 1 solar mass up in the nearest
step, 2 in the next, then 5, then 8, and a sector already generated never
gets any. `--backfill-from none` skips it (`edge` is the default).
`--then-scatter` first scatters the bright stars galaxy-wide, as
`planetgen scatter --mass` and `scatter --luminosity` do, and then the phenomena,
as `scatter --phenomena` does, leaving out every sector already generated; a new
galaxy uses it so the scatter never draws stars for sectors the run fills
anyway. With it, `--bright-star-min-luminosity` (9000 solar luminosities by
default), `--phenomenon-min-mass` and `--compact-min-mass` set the
scatter, as on `scatter`, which is described below.

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
of one ring of the old cylindrical grid took 16.5 s on one worker and 9.8 s on the default two (the grid has since changed; the speed per sector has not been measured again).

`planetgen plan` draws its bright stars the same way, one shell of the
galaxy per task, densest shells first (`--workers` works there too).
Each shell draws from its own random stream, so the same seed places the
same stars on any number of workers. The seed is the galaxy seed that
`planetgen plan --seed <32 hex digits>` sets (see "The galaxy seed" under
[Size and time estimates](#size-and-time-estimates)); each shell's stream
comes from it and the shell's address.

Every progress bar shows the time elapsed and the time remaining until
it's done. The remaining time comes from a decaying average of how many
sectors (or shells) finished per second, weighted toward the last minute
or so, so it follows the run's current speed and doesn't jump about when
several workers finish at once. The Generate page shows the same
estimate ("about 4 m 10 s left").

A `galaxy` or `plan` run that does more than one stage also draws one "Whole job (stage 3 of 7)" bar at the top
of the display (PERF.55). It counts the time since the run began against the stored time of all its stages (the
same stage times the Generate page's overall bar uses), so it shows the time left across every phase; with no
stored time for a stage it shows only the elapsed time. A run that outlasts its estimate keeps the bar just short
of full. When the run ends it stores two more rows in the stage history, with the run's settings: "Whole job"
(the seconds of the whole run, the stages it ran and the shells the scatter stages visited, those that placed
nothing included, per second; `generation_stats.job_shells_per_second`) and, when scatter stages ran, "All
scatters" (their operations, shells, objects and seconds, per shell and per operation). Each stage's own row also
carries its seconds per shell. A scatter with no stored time of its own is estimated from the "All scatters"
figure, so every bar and estimate reads the one stored-time function (`generation_stats.stage_seconds`).

The bright-star bar ("Bright stars (12 of 1,271 shells)") shows a share
done rather than a count: before the scatter starts, every shell gets an
expected star count from the same density model the scatter draws from
(a quick sample of its bands, a second or two for a whole galaxy), plus
a small allowance per band, so the nearly empty outer shells count for
little and the dense inner ones for a lot.
Shells in progress report their stars a few times a second and the bar
moves with them, so the time left follows the work left rather than the
shells left. While shells take longer than 30 seconds each, a second bar
appears under it with the stars of the shells being drawn, done of their
estimate, and its own time left ("Shell 0: stars"); it goes again once
shells finish faster than one every 20 seconds. Sector fill keeps its one
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
shells as its leaves, each with its own start, end and duration. A run
the Generate page started hangs under that page job and its step. This
happens with one worker too (without the lease), and is skipped without
the control database.

### Size and time estimates

Before a bulk run writes anything (`planetgen galaxy` in every mode,
`planetgen sector --num-sectors`), it works out how big and how long
it will be and prints it:

```
Estimate for shell 12: About <size> and <time> for <N> sectors (~<n> star systems,
~<n> stars, 2 workers). <free> GB free of <total> GB.
```

(The figures depend on the galaxy, the density and this server; the line names
the run, `shell 12 band 0` for a single band, and so on.)

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
- **A run of more than 20,000 sectors** is estimated from 200 of them, spread
  evenly through the run, and the line says so ("worked out from 200 sampled
  sectors"); its sectors are resolved as they are queued, never all held in
  memory. A radial fill (`--center-sector` and `--radius-pc`) asked for by the
  Generate page's estimate, or the sector page's neighborhood action, stops at
  2 million sectors in the ball and says to choose a smaller radius; the command
  line has no such limit.

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
row plus one row per shell, from the core outward, naming the highest band
that shell reaches; the outline is symmetric about the plane) that `planetgen galaxy` consults to decide, per
address, whether anything exists there at all before generating it lazily
on demand (the `galaxy_shell` table in
[`database-schema.md`](database-schema.md)). Every `planetgen galaxy` mode checks its
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
`scatter --band` and each backfill block draws from its own seed,
the SHA-256 of the galaxy seed and its address (`sector:12/3/0` is shell 12, band 3, wedge 0), so one
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
(0007007F000160030C0300), run: galaxy --shell 3` (`none yet` before the
first plan; OPS.10). It goes to the console, the `--debug` file and the
debug log, and a Generate page job's log starts with the same line for
the job before its first step. Design:
[`design/reproducible-galaxies.md`](design/reproducible-galaxies.md).

`planetgen plan` draws the shape of the galaxy from its seed (GEN.401), and each part can be
set as an option: `--disk-scale-length-pc` and `--disk-scale-height-pc` for the thin disk,
`--bulge-scale-radius-pc` and `--bulge-amplitude` for the bar bulge, `--arm-count` and
`--pitch-angle-deg` for the spiral arms. An option left out is drawn: the galaxy's stellar mass
comes from the seed (`galaxy/mass_chain.py`), the disk length from the size-mass relation of
late-type galaxies, the height from the length, the bulge's share and size from the galaxy's
stage and the arm count and pitch angle from the arm statistics
([`design/intergalactic-and-halo-contents.md`](design/intergalactic-and-halo-contents.md),
section 5.13, rows 4 to 6). The plan prints what was drawn. `--arm-density` (1.4) and
`--interarm-density` (0.6) set the arms' contrast, `--nuclear-amplitude` (1: the nuclear star
cluster and nuclear stellar disc at the centre, scaled to the bulge; 0 is none).
The galaxy holds the systems its stellar mass makes (no point of it is pinned to the local density, GEN.404). `--max-shell` caps how far out the
shells are walked (100,000 by default) and `--seed <32 hex digits>` sets the
galaxy seed. `--workers` is how many shells of bright stars are drawn at
once.

`planetgen scatter` places the galaxy's stars and phenomena on the stored outline, before any sector is
filled (PERF.121). It draws exactly one kind of scatter per run, so each kind is a command of its own that a job
can queue, retry or leave out:

- `--mass` places every star born at the stellar mass limit or more, whatever its luminosity.
  `--phenomenon-min-mass` is one of 8, 10, 12, 14, 16, 18 or 20 solar masses (14 by default).
- `--luminosity` places the stars born lighter than the mass limit that are at least
  `--bright-star-min-luminosity` bright: from 2500 to 4,000,000 solar luminosities, 9000 by default (the Generate
  page offers presets). A lower level places far more stars and takes longer, and the run prints an estimate first.
  It works at the mass limit the mass pass stored, and stops with "run 'planetgen scatter --mass' first" when
  there is none.
- `--phenomena` places the neutron stars, black holes and the rest. `--compact-min-mass` is the lowest mass of a
  neutron star or black hole placed galaxy-wide: one of 1, 2, 4, 6, or `star` for the stellar limit, which is the
  default. A lower limit stores far more rows, and `scatter --help` lists the sizes. The stars are left alone and the
  central black hole or quasar is always rebuilt.
- `--nebula-field` draws the `nebula_field` table of molecular clouds again from the galaxy seed.
- `--band --down-to N` adds a band of dimmer stars (described below).

Each star is a finished star at a fixed point in its sector, stored in `bright_stars`, so the Galaxy Map can show the
bright stars tracing the spiral arms right away. Filling a sector later
builds a full system around each of its bright stars first, then draws
the rest of its systems from dimmer stars only, so its expected total is
unchanged. Every system's age comes from the stellar population mix where
its sector sits (young stars crowd the arms near the plane; the bulge is
old). Every option belongs to one kind and is refused with another. `planetgen plan` takes none of them, and
a plan alone places no stars (and empties what an earlier plan's scatters left); `galaxy --then-scatter` runs the
kinds as steps of one run.

A scatter draws only where none of its kind is stored (Boss's F7): it stops and names the clear to run first, so
nothing is replaced by accident. `planetgen clear-scatter` empties one kind and nothing else: `--mass` (refused while
the bright stars are stored, which are split from the massive stars at the mass limit, so clear `--luminosity`
first), `--luminosity` (the band goes with it), `--phenomena` or `--nebula-field`. To redo a scatter at a new
setting, clear it and draw it again; to change the mass limit, clear `--luminosity` and `--mass`, then draw
`--mass` and `--luminosity`. Clearing a kind with nothing stored says so and changes nothing.


`--cell-star-ceiling` and `--fill-luminosity-cap` (GEN.281, `plan`) set, for the whole galaxy, how deep a crowded cell is filled: a cell that expects more systems than the ceiling (default 5,000) is filled only with the stars at or above the cap (default 100 solar luminosities, raised for a cell those stars would still fill past the ceiling, never under the brightest white dwarf), and the rest is reported as unresolved on its sector page and in the API (`unresolved`). Ordinary space never reaches the ceiling. Both are stored with the galaxy and apply to sectors generated afterwards. `--core-min-luminosity`, `--core-mass-limit` and `--halo-min-luminosity` (DB.54, `plan`) store the scatter limits by place: the core's luminosity floor and mass limit (unset, they are the baseline's times 1.1, 9,900 and 15.4 at the defaults) and the halo's luminosity floor (default 500, from 50 to 4,000,000; the halo has no mass limit of its own). The tiers that read them are GEN.333's. The satellites and clusters (`galaxy_components`) are drawn only by a full `plan`, so a galaxy planned before GEN.159 and GEN.237 has an empty component table until a full plan runs. The measured cost of a full scatter is in [generation-performance-study.md](design/generation-performance-study.md).

`planetgen clean-nebula-field` (GEN.205) brings an existing galaxy's stored molecular clouds in line with the current nebula keep rates (`tuning.NEBULA_FIELD_KEEP`): without `--apply` it counts, per class, the stored field clouds the rule now rejects and changes nothing; `--apply` deletes them (everything inside them loses its `inside_nebula_id`, and containment is refreshed in the sectors they reached). A nebula that was not made by the field is left alone. `planetgen scatter --nebula-field` then draws the `nebula_field` table again.

`planetgen fill-map-luminosity` (DB.23, GEN.201) fills `phenomenon_scatter.map_luminosity_sol`, the luminosity a star of a row's mass would have on the Galaxy Map, for a galaxy scattered before schema v84: id-range batches of 200,000 rows, each its own short statement and commit, resumable (it only fills rows still empty). A new scatter writes the value with each row. `planetgen rescore-habitability` (GEN.217) re-derives the stored habitability numbers of every planet and moon from its stored conditions: the flux swing of a body in a double-star system (loading that system, as generation does), the PHI-6 index, the six domains and their tiers, the equipment tier and the note. It is for a galaxy generated before a change to the rules (GEN.148, GEN.207, GEN.215), so the planet chip, the search page and the facet counts agree with the description, which is worked out from the same conditions whenever a page is drawn. Nothing else is touched; each batch of 500 systems is its own transaction and the result is the same however often it is run, so a stopped run is started again. A gas giant and a body saved before schema v76 (no score) are left alone. Both take only `--workers` and the connection options.

The scatter always leaves the sectors that are already filled out (they
would never get their bright stars).

The scatter can go down in stages. `scatter --band --down-to N` keeps every
bright star already placed and adds only those from `N` up to (not
including) the galaxy's star-fill level, the threshold already scattered
(9000 by default), then lowers that level to `N`. So a quick test galaxy
scattered at 9000 can later go down to 5000 without redrawing anything
brighter:

```bash
planetgen scatter --band --down-to 5000
```

Sectors already filled get none of the new stars: their own systems were
drawn below the old level, so they already hold stars that bright. Asking
for a level at or above the current one does nothing and says so. The
Generate page shows the current level and offers the same step as "Add a
dimmer step of bright stars".

## Exotic Phenomena Generation

`planetgen phenomenon` generates a single exotic stellar phenomenon on demand,
independent of any one sector (see "Sector-Level Exotic Phenomena" above
for the population every generated sector gets automatically):

```bash
planetgen phenomenon --type {black-hole,neutron-star,nebula,supernova-remnant,rogue-planet,comet,asteroid-field,quasar} [options]
```

Omitting `--type` picks uniformly at random among the first seven; a quasar is only made when asked for by name, and `--sector-id` only accepts a shell-0 (galactic core) sector for one, placing it at the galactic center; there is one per galaxy. `--anchor-system`
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
`--name <name>` forces the phenomenon's name. `--no-settle` (with `--sector-id`) skips saving the sector paths of that sector and the sectors around it, which a black hole or neutron star can bend (`planetgen orbits` saves them later). `--debug`/`--quiet`/`--silent`, and the
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
sector's generated content, one line a sector (`shell band wedge
digest`, unplaced sectors last), then the plan's digest and the
region's. Two builds of one galaxy are the same when their digests
match; a sector whose line differs is where they part. It writes
nothing and isn't logged as a run.

```bash
planetgen fingerprint                      # the whole galaxy, plan included
planetgen fingerprint --shell 3 --shell 4  # every sector in shells 3 and 4
planetgen fingerprint --sector 12 3 0      # one sector: shell 12, band 3, wedge 0 (repeatable)
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

## Check the database

`planetgen check-db` looks for damage and changes nothing: the schema revision, table health (`CHECK TABLE`), rows without parents, ids, impossible values, sector counts and version keys. It ends with one pass or fail line per check, exits 1 when it found damage and 2 when a check could not run (DB.8, DB.21, PERF.126). The Generate page runs the same checks.

Each check is a command of its own (`check-revision`, `check-tables`, `check-tables-extended`, `check-orphans`, `check-ids`, `check-values`, `check-counts`, `check-version-keys`, `check-systems`, `check-object-index`; `planetgen commands` lists them). `check-db` only builds the list of them and the wrapper queues them one after the other as a check run on the Redis queue, so it needs Redis (and a worker running as part of it) like every queue run, never holds up a generation run and shows no banner to the site's visitors. Name the checks to run only those. `object-index` is not in the default list: it compares the object index with the tables (PERF.140) and warns when the index has not been built yet.

```bash
planetgen check-db                              # the default checks, the whole database
planetgen check-db ids orphans                  # only those checks, in the order the list keeps
planetgen check-db --shell 12                   # values, counts and version keys of one shell's sectors (repeatable)
planetgen check-db --sector 12 3 0              # the same for one sector: shell 12, band 3, wedge 0 (repeatable)
planetgen check-db systems --estimate-only      # print the time the systems check would take, then stop
planetgen check-db tables-extended systems --yes   # the slow checks: CHECK TABLE ... EXTENDED, then every star system
planetgen command-run check-ids '{"database": {...}}'   # one check alone
```

*   Without a name the list is `revision tables orphans ids values counts version-keys`; `tables-extended` (`CHECK TABLE ... EXTENDED`) and `systems` (validates every star system) are the slow ones and run only when named.
*   `systems` is slow, so its time is estimated first and a yes is asked for; `--yes` goes ahead after the estimate without asking (needed off a terminal), and `--estimate-only` prints the estimate and stops (the Generate page uses it).
*   `--wide-binaries` (with `systems`): validates only the wide binaries (S-type pairs), to count the pairs whose orbits overlap (GEN.228).
*   `--shell` and `--sector` scope `values`, `counts`, `version-keys` and `systems`; the other checks always look at the whole database.

## Version keys

`planetgen versions` lists the version key each update recorded for the galaxy database, newest first (OPS.13): the release, the Python version, the operating system and the architecture that made the galaxy, in the 22 hex digits explained under [the galaxy seed](#size-and-time-estimates). It takes only `--mysql-*`, `--debug` and `--quiet`.

## Benchmark

`planetgen benchmark` times a small galaxy from the plan to the finished fill in scratch databases and prints the time and database work of every stage, at each worker count asked, with an optional profile of the fill (PERF.31). It makes a scratch database and a scratch control database, runs the real `plan` and `galaxy` commands against them in child processes with a fixed seed, and drops both afterwards unless `--keep` is given; nothing is written to the real database.

```bash
planetgen benchmark                          # 2 sectors in shell 12, one worker
planetgen benchmark --sectors 20 --workers 1,2,4 --profile --output bench.json
```

*   `--sectors N`: the sectors the fill makes (2 by default). `--shell I`: the shell the fill starts in (12 by default; a lower shell is denser). `--max-shell I`: the plan's outermost shell (12 by default, a small galaxy). `--seed HEX`: the galaxy seed (32 hex digits).
*   `--workers N[,N...]`: the worker counts to run at, one run each (1 by default).
*   `--profile`: also profiles one single-worker fill and lists its costliest functions. `--keep`: leaves the scratch databases in place. `--output FILE`: writes the results to `FILE` as JSON.

## Warm the caches

`planetgen warm` builds the site's caches so the first visitor after a release does not wait for them (MAP.134, PERF.104, PERF.127; design in [cache-warm-up.md](design/cache-warm-up.md)). `update.sh` runs it in the background after every update; run it by hand after clearing a cache, as the web server's user so Apache can read what it writes. Each cache family is a command of its own, listed by `planetgen commands` and runnable alone with `planetgen command-run KIND ARGS`: `warm-object-index` (the object index: one row per object of every kind the site lists, written when it changed since the last run; built whole after `planetgen orbits` or a new epoch, because those move objects between sectors without touching their rows' `modified_at`; PERF.140), `warm-galaxy-map` (the opening tiles and stage, the shape, the sector and phenomena lists), `warm-galaxy-tiles` (the Galaxy Map's tiles: the first `tile_cache.warm_levels` levels of the tile octree in full, cut into chunks of 128 tiles that the workers build; PERF.130), `warm-galaxy-deep-tiles` (the zoom run in to the galaxy's centre, down to level `tile_cache.warm_zoom_levels`; PERF.141), `warm-lists` (the first page of every list and its counts), `warm-sectors`, `warm-systems` and `warm-phenomena` (every page of that kind, nearest the galaxy centre first, up to `page_cache.warm_max_objects` of each) and `warm-search`. `planetgen warm` only builds the list of them (every family in that order, or the families named) and hands it to the wrapper, which queues them one after the other on the Redis queue as a maintenance run: it does not wait for a generation run, holds none up and shows no banner to visitors. A list of every family ends with `end-map-rebuild-maintenance`, which opens the site if the sector-grid switch closed it (UX.112). It waits for a running update to finish, prints one line per family (how many pages it warmed, how many answers it built, how many were already cached and how many the cap left out), and exits 1 naming the command that failed, or when Redis does not answer or another warm-up is running; there is no silent skip.

```bash
planetgen warm [FAMILY|PHASE ...] [--cap N]
```

The warm-up runs in five phases (PERF.141): `opening` (`object-index`, `galaxy-map`, `lists`: the object index, the opening views and the first page of each list), `built-levels` (`galaxy-tiles`), `deep-tiles` (`galaxy-deep-tiles`), `objects` (`sectors`, `systems`, `phenomena`) and `search`. After each phase it prints one line: what it built, what was already stored, what the cache budget (`cache.max_mb`) or the cap left out, or `unchanged, nothing built` when the data under it has not changed since the last run. After a generation run, `planetgen orbits`, an edit or any other change that moves a data version tag it also starts by itself in the background (`planetgen warm --on-change DATABASE`, which waits `cache.warm_settle_seconds` for more changes and then warms once; off with `cache.warm_on_change`, never during an update). `FAMILY` is any of the families above and `PHASE` any of the phase names, a phase standing for all its families (the list is those, in order, whatever order they are named in); `--cap` overrides `page_cache.warm_max_objects`.

## Helper tools

Besides `planetgen`, a few modules run on their own with `python3 -m`. Run them as the account that owns the install (the web server's user for `planetgen warm`, so Apache can read what it writes). `orbits`, `query`, `migrate`, `reset`, `forget_devices` and `reset_two_factor` take the same five database options: `--mysql-host`, `--mysql-port`, `--mysql-user`, `--mysql-password` and `--mysql-database`, each defaulting to its `PLANETGEN_MYSQL_*` variable (see [config.md](config.md)). `worker`, `config`, `lockouts` and `lift_lockout` take no database options. `orbits` and `query` also take `--version`.

| Command | What it is for |
|---|---|
| `python3 -m planetgen.cli.orbits` | Moves every object that is due along its orbit. |
| `python3 -m planetgen.cli.worker` | Runs one queue worker. |
| `python3 -m planetgen.cli.query` | Lists what is stored. |
| `python3 -m planetgen.cli.migrate` | Brings the database up to the current schema. |
| `python3 -m planetgen.cli.datasets` | Downloads the static data the generator reads (isochrone grids and the like) or says where each set stands. |
| `python3 -m planetgen.cli.reset` | Empties the galaxy. |
| `python3 -m planetgen.cli.lockouts` | Lists the login lockouts. |
| `python3 -m planetgen.cli.lift_lockout` | Lifts a login lockout. |
| `python3 -m planetgen.cli.forget_devices` | Revokes an admin's trusted-device cookies. |
| `python3 -m planetgen.cli.reset_two_factor` | Turns off an admin's two-factor sign-in. |
| `python3 -m planetgen.cli.config` | Checks, documents and edits the settings. |
| `python3 -m planetgen.cli.backup` | Takes, lists, verifies and restores backup archives. |

### orbits

The positional update. Each planet, moon, star, binary pair, comet, facility and standalone phenomenon keeps its own clock (`epoch_unix`) and a time it next needs storing (`next_update_due`); a run moves only what is due, from its own epoch to now by the database server's clock, then refiles anything that drifted into another sector, refreshes containment (nebulae and remnants), nearest systems and location text, and recomputes the sector paths. A run prints how many years passed since the last one, what it moved and how many objects changed sector. It never moves anything backwards: if the server's clock went back it says so, moves nothing and restarts the clock. There is no `update-orbits` command; this is it. Run it from cron or by hand, as often as you like (a run with nothing due changes nothing):

```bash
python3 -m planetgen.cli.orbits
```

### worker

One RQ worker, at the lowered priority generation workers run at. The site and the generation runs start these themselves (`planetgen.queue.workers.ensure_workers`); run one by hand only to add capacity or to watch one. Redis is read from `--url`, else `PLANETGEN_REDIS_URL`, else `redis.url` in `config.json`.

*   `QUEUE ...`: the queues to serve, first one first (default: the interactive queue). A worker serves one kind of queue, the web queues or a run's own, never both.
*   `--burst`: exit once the queues are empty.
*   `--idle-seconds <int>`: exit after this long without work (default 60; ignored with `--burst`).
*   `--name <name>`: the worker's name (default: RQ picks one).
*   `--url <redis-url>`: the Redis server.

```bash
python3 -m planetgen.cli.worker --burst
```

### query

Lists what is stored, read-only. The subcommands are `sectors` (every sector with its size and system count), `systems`, `near`, `planets` and `moons`.

*   `systems [--star-type TYPE] [--sector-id ID]`
*   `near PLACE --distance PC [--kinds a,b] [--limit N] [--offset N]`: everything within a distance (at most 50 pc) of an object (`system:<ID>`, `planet:<ID>`, `nebula:<ID>` ...; a bare three-part ID is a system) or of a galaxy-frame point `x,y,z` in parsecs.
*   `planets` and `moons`: `[--class C] [--min-radius-km N] [--max-radius-km N] [--sector-id ID] [--system-id ID]`

```bash
python3 -m planetgen.cli.query sectors
python3 -m planetgen.cli.query planets --class M --sector-id <ID>
```

### migrate

Applies any schema migration the database is missing and does nothing if it is current. `install.sh` and `update.sh` run it, so you rarely need to. It shows a progress bar (a line per step when output is not a terminal). `--status` prints `<current version> <target version> <pending steps> <database>` and changes nothing, which is how `update.sh` decides whether to ask about the database. A database older than schema v61 is refused.

```bash
python3 -m planetgen.cli.migrate --status
```

### datasets

Downloads the static data files the generator reads (GEN.548) into a cache folder outside the checkout, `/var/cache/planetgen/datasets` (`PLANETGEN_DATA_DIR` moves it). `install.sh` and `update.sh` run `fetch` once at setup, as root, and make the folder readable for Apache's user; nothing downloads while the site or the generator runs. Each file is requested with the standard `User-Agent`, must have the size and sha256 the catalogue states, and is kept by an atomic rename, so an interrupted run leaves no half file. A dataset whose host does not answer, or whose file does not match, stays missing; `fetch` then prints the reason and exits 1, which the scripts report as a warning and not as a failure. `status` fetches nothing, says which datasets are present and exits 1 while one is missing. **Bulk generation needs them (GEN.590):** `galaxy`, `plan`, `population` and `sector --num-sectors` above 1 refuse to start while a dataset is missing, naming it and writing nothing, as with the math check; so do the Sector page's neighbourhood button and the Generate page's jobs, whose second step, "Check the data", is `status` (before the reset, so nothing is wiped for a run that cannot start).

```bash
python3 -m planetgen.cli.datasets fetch
python3 -m planetgen.cli.datasets status
```

### reset

Empties every generated table of the galaxy (`TRUNCATE`, so the schema stays), keeping `schema_migrations` and the id counters (`id_blocks`), so the next galaxy's ids carry on past the old one's. It never touches the control database (admin accounts, sessions, API keys, audit log). **It cannot be undone.** By default it asks you to type the database's name back; `--yes` skips that for scripts. Look first with `--dry-run`, which lists each table and its row count and changes nothing.

```bash
python3 -m planetgen.cli.reset --dry-run
```

### lockouts

For an admin locked out of the site. It lists what is locked and does nothing else; the lockouts live in Redis, on the rate-limit storage (`ratelimit.storage_uri`, else `redis.url`).

```bash
python3 -m planetgen.cli.lockouts
```

### lift_lockout

Lifts login lockouts, on the same storage. One of these is needed:

*   `--ip ADDRESS`: lift one address (an IPv6 address by its /64).
*   `--user NAME`: lift one username.
*   `--all`: lift every lockout.

Each lift is written to the activity log as `DB lockout.lift`.

```bash
python3 -m planetgen.cli.lift_lockout --ip 203.0.113.5
```

### forget_devices

`forget_devices USER`: that admin's browsers lose their trusted-device cookies, so the next sign-in asks for the two-factor code again. Written to the activity log as `DB devices.revoke`.

### reset_two_factor

`reset_two_factor USER`: turns off that admin's two-factor sign-in, for a lost phone and lost recovery codes. The admin signs in with the password and can set it up again. Written to the activity log as `DB totp.reset`.

### config

`config.json` is the one file that holds all program configuration, your own copy of `config.json.example`. `check` validates it the way the app reads it (an unknown key is an error) and exits 1 naming every field that does not fit. `docs` rewrites the generated parts, the option table in [config.md](config.md), `config.json.example` and `config.schema.json`, from the settings model; a test fails when they differ from what it would write, so run it after changing an option.

`get [OPTION]` prints the value that applies and where it comes from (default, `config.json` or an environment variable; a secret shows only `(set)` or `(empty)`). `set OPTION=VALUE ...` and `unset OPTION ...` change `config.json` through the same writer as the Admin settings page ([config.md](config.md#saving-settings)): every value is checked first and one bad value saves nothing, the old file is backed up, and each change gets an audit row (when the control database is reachable) and an activity-log line. The command prints what needs a restart. A value is read as JSON for a number, `true` or a list and as text otherwise (`set site_name="My galaxy" cache.max_mb=300`); `OPTION=-` reads it from standard input, which keeps a secret out of the shell history. `restore ARCHIVE` puts the settings of a backup archive back (`planetgen.cli.backup list` shows the archives), backing up the current file first. Every option may be set here, the database, the secrets and the Owner-only ones included, because whoever can run the command owns the server; the Admin page cannot change those. If `config.json` does not exist the command says so and creates nothing.

```bash
python3 -m planetgen.cli.config check
python3 -m planetgen.cli.config docs
python3 -m planetgen.cli.config get cache.max_mb
python3 -m planetgen.cli.config set cache.max_mb=300 site_name="My galaxy"
python3 -m planetgen.cli.config unset cache.max_mb
python3 -m planetgen.cli.config restore 20261016T193012Z.tar.bz2
```

### backup

Backup archives (ADM.87, [config.md](config.md#backups)): one `.tar.bz2` per backup, holding a part per kind of data (`settings`, a copy of `config.json`, now; the control database and the DB repair parity data in 9.0). Run it as root or the web user (the archives are mode 0600).

```bash
python3 -m planetgen.cli.backup take                    # every part (--part NAME for some)
python3 -m planetgen.cli.backup list                    # newest first, with each archive's parts
python3 -m planetgen.cli.backup verify ARCHIVE          # check every hash
python3 -m planetgen.cli.backup restore ARCHIVE         # the whole archive, or --part NAME for some
python3 -m planetgen.cli.backup parts                   # the parts this build can back up
```

A restore verifies the archive, backs up what it is about to replace in a safety archive, skips and reports a part this build has no provider for, and writes an audit row. `settings.backups_kept` of the archives are kept.

## Planned commands

Not built yet. Each names its TODO item and phase; the design is in
[`design/reproducible-galaxies.md`](design/reproducible-galaxies.md).

| Command or option | Item | Phase | What it will do |
|---|---|---|---|
| `planetgen repair-db` | DB.9 | 1 | Rebuild damaged sectors from the parity file, or regenerate them from their seed when the version key matches, then check again. |
| `--resume` | PERF.30 | 1 | Finish the runs an interrupted fill left, from the step each reached. |

Work runs as RQ jobs on Redis from phase 0 (PERF.24), so `--workers` sets the RQ worker count.

Phase 2 also adds a daily maintenance run, `scripts/maintenance.sh`, set up as a
scheduled job by install and update (OPS.16, OPS.17): the positional
update (`planetgen.cli.orbits`). (The daily settings-file merge and the
18 backup slots, GEN.61 and OPS.18, were dropped on 2026-10-09.)
