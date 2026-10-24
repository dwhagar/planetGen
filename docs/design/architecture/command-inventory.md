# Every command does one thing: the command inventory

Status: plan (Boss, 2026-10-17, 4:19 PM EDT: "redesign ALL the commands like that so we can build whatever sequences we want"; the principle is Boss's of 4:18 PM EDT, in [redis-and-the-work-queue.md](../operations/redis-and-the-work-queue.md) section 5.12). Read from the code on main on 2026-10-17 (4:30 PM EDT); the items that reshape each command are named in the last column.

## The principle

Modular commands: each command does one and only one thing and does it exceptionally well, and never worries about anything but what it does. No command knows about another stage, about what ran before it or what runs after it, or about the order. A list of steps is built from what the user asked for and handed to the wrapper (PERF.118), which alone handles the sequence; the order lives only in the wrapper's step list, kept as data, so any sequence can be built and changed without touching a command. Each command has its own parameters, set by the UX when it runs the command through the API and by Boss on the command line.

Naming (decided by Boss, 5:14 to 5:16 PM EDT on 2026-10-18, cards F1 to F7): one root command per kind of work and exactly one kind chosen per run, as in Boss's example `scatter --mass`. A run of `scatter` takes one of `--mass`, `--luminosity`, `--phenomena`, `--nebula-field`, `--band`, and only that kind's parameters. The map-luminosity fill, the habitability re-score and the nebula-field clean-up are one root command each, named like the registered commands: `fill-map-luminosity`, `rescore-habitability` and `clean-nebula-field` (with `--apply`). Redoing a scatter is a separate `clear-scatter` command, not a `--redo` flag, and `--redo-scatters` is removed with no shim. Every start path (the Generate page, the API jobs, `planetgen sector` and `planetgen system`) builds a step list and hands it to the wrapper: the Generate page builds the list and one wrapper run follows it (`job_runner` only starts it detached), each API job is a small run of its own, and reset, check-math and the like are one-unit registry commands.

## What exists today (`planetgen` and `python3 -m planetgen.cli.*`)

| Command | What it does today | One thing? | Reshape |
|---|---|---|---|
| `planetgen plan` | The skeleton and nothing else (PERF.121, done) | Yes | none |
| `planetgen scatter` | One kind of scatter per run: `--mass`, `--luminosity`, `--phenomena`, `--nebula-field` or `--band --down-to N`, with only that kind's parameters; refuses to draw over stored rows (PERF.121, done) | Yes | none |
| `planetgen clear-scatter` | Empties one kind of scatter's stored rows: `--mass`, `--luminosity`, `--phenomena` or `--nebula-field`, run before `scatter` to redo it (PERF.121, Boss's F7, done) | Yes | none |
| `planetgen fill-map-luminosity`, `rescore-habitability`, `clean-nebula-field` | One fill or clean-up each, with no setting of a scatter (PERF.122, done) | Yes | none |
| `planetgen galaxy` | Sector generation in nine modes (random start with neighborhood, shell, block, center sector, band, wedge, course), the neighbour link, the three scatters (`--then-scatter`), the backfill, the population pass (`--population`) and the settle | No: eight jobs in one run | PERF.123 |
| `planetgen sector` | Generates sectors; `--population` also runs the population pass | Almost: population leaves it | PERF.124 |
| `planetgen phenomenon` | Generates one phenomenon, then settles its sector unless `--no-settle` | Almost: settle leaves it | PERF.124 |
| `planetgen system` | Generates one system and prints or writes it | Yes | none (registered by PERF.118) |
| `planetgen population` | The population pass and the territories pass (`--territories-only`, `--rescan`) | No: two passes | PERF.125 |
| `planetgen check-db` | Builds the list of the nine `check-*` commands (one per check: revision, tables, extended tables, orphans, ids, values, counts, version keys, systems) the arguments name, and hands it to the wrapper as a check run | Yes: each check is its own command; this only builds the list (PERF.126, done) | PERF.126 |
| `planetgen check-math`, `versions`, `fingerprint` | The math check; list the version history; print the content digests | Yes | none |
| `planetgen benchmark` | Runs `plan`, then `galaxy` over a small shell, in child processes and prints the timings | One purpose, but a fixed sequence of its own | PERF.128 |
| `cli.migrate`, `cli.reset`, `cli.worker`, `cli.query`, `cli.maintenance`, `cli.update_guard`, `cli.config`, `cli.backup` | Bring the schema up; wipe the generated rows; run a queue worker; list stored rows; switch maintenance mode; close the site for an update; check and edit settings; take, list, verify and restore archives | Yes (one verb each) | none |
| `planetgen warm` (was `cli.warm_map`) | Builds the list of the seven `warm-*` commands (one per cache family) and the closing `end-map-rebuild-maintenance`, and hands it to the wrapper as a maintenance run | Yes: each family is its own command; this only builds the list (PERF.127, done) | PERF.127 |
| `cli.lockouts`, `cli.version_history` | List or lift login lockouts by flags; record the version key, or `--list` it (which `planetgen versions` already does) | Yes, done (OPS.60, 2026-10-17, 4:36 PM EDT): `cli.lockouts` only lists, `lift_lockout.py` lifts, `cli.version_history` only records | OPS.60 |
| `cli.render_parity`, `cli.dedupe` | One-off checks for a schema (v29) and a naming rule long since migrated | Done, removed (OPS.61, 2026-10-17, 4:36 PM EDT): `render_parity.py` and `dedupe.py` are deleted | OPS.61 |
| `cli.orbits` | Advances every planet, moon, star and binary by the real time since the last run | Replaced by the nightly orbital update of GEN.105 (8.3), which follows this principle | note on GEN.105 |
| `update.sh`, `install.sh` | Shell sequences of the commands above | They are wrappers by nature and stay shell scripts that call the commands | none |
| Generate page, API jobs, one-off pages | Start runs | Each must build a step list for the wrapper | PERF.119 |

`scripts/*.py` are development tools (the test runner, the docs builder, the version bump, the lock files), not product commands, and are out of this scope.
