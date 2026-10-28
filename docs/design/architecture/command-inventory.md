# Every command does one thing: the command inventory

Status: plan (Boss, 2026-10-17, 4:19 PM EDT: "redesign ALL the commands like that so we can build whatever sequences we want"; the principle is Boss's of 4:18 PM EDT, in [redis-and-the-work-queue.md](../operations/redis-and-the-work-queue.md) section 5.12). Read from the code on main on 2026-10-17 (4:30 PM EDT); the items that reshape each command are named in the last column.

## The principle

Modular commands: each command does one and only one thing and does it exceptionally well, and never worries about anything but what it does. No command knows about another stage, about what ran before it or what runs after it, or about the order. A list of steps is built from what the user asked for and handed to the wrapper (PERF.118), which alone handles the sequence; the order lives only in the wrapper's step list, kept as data, so any sequence can be built and changed without touching a command. Each command has its own parameters, set by the UX when it runs the command through the API and by Boss on the command line.

Naming (decided by Boss, 5:14 to 5:16 PM EDT on 2026-10-18, cards F1 to F7): one root command per kind of work and exactly one kind chosen per run, as in Boss's example `scatter --mass`. A run of `scatter` takes one of `--mass`, `--luminosity`, `--phenomena`, `--nebula-field`, `--band`, and only that kind's parameters. The map-luminosity fill, the habitability re-score and the nebula-field clean-up are one root command each, named like the registered commands: `fill-map-luminosity`, `rescore-habitability` and `clean-nebula-field` (with `--apply`). Redoing a scatter is a separate `clear-scatter` command, not a `--redo` flag, and `--redo-scatters` is removed with no shim. Every start path (the Generate page, the API jobs, `planetgen sector` and `planetgen system`) builds a step list and hands it to the wrapper: the Generate page builds the list and one wrapper run follows it (`job_runner` only starts it detached), each API job is a small run of its own, and reset, check-math and the like are one-unit registry commands.

## Rules for a new command (PERF.120)

A command is added only if it follows all of these. `src/tests/test_command_rules.py` holds the structural ones for every kind in `queue/commands.REGISTRY`.

1. **One job.** It does one thing and nothing else, and it is named for it (`scatter-mass`, `run-territories`). If a change makes a command do a second job, the second job is a command of its own.
2. **Its own parameters.** It takes only the settings its job needs, as plain JSON in the step's `params`; it never reads another stage's settings.
3. **No knowledge of the order.** It does not know what ran before it, what runs after it or where it sits in the list; it does not call another command. The order lives only in the step list `commands.build_steps` makes, which is data.
4. **A registered kind.** An entry in `commands.REGISTRY` with a `unit(params, report)` that returns a JSON-able result and an `expand(params)` that cuts its stage into units (one unit, with `whole=True`, when the job is one by nature; otherwise `unit_seconds` says how long the longest unit runs). There is no inline command and no second way to run it.
5. **Runnable alone.** `planetgen command-run KIND ARGS` runs one unit from its JSON, and the command-line verb (`planetgen <verb>`) is a run of that one command through the wrapper, with the same rows stored as a run inside a longer list.
6. **A test that runs it alone,** with the arguments it takes and the result it gives.
7. **A docs line:** its kind name in this inventory or in the document that describes its job, and the command-line reference if it has a verb.

## What exists today (`planetgen` and `python3 -m planetgen.cli.*`)

| Command | What it does today | One thing? | Reshape |
|---|---|---|---|
| `planetgen plan` | The skeleton and nothing else (PERF.121, done) | Yes | none |
| `planetgen scatter` | One kind of scatter per run: `--mass`, `--luminosity`, `--phenomena`, `--nebula-field` or `--band --down-to N`, with only that kind's parameters; refuses to draw over stored rows (PERF.121, done) | Yes | none |
| `planetgen clear-scatter` | Empties one kind of scatter's stored rows: `--mass`, `--luminosity`, `--phenomena` or `--nebula-field`, run before `scatter` to redo it (PERF.121, Boss's F7, done) | Yes | none |
| `planetgen fill-map-luminosity`, `rescore-habitability`, `clean-nebula-field` | One fill or clean-up each, with no setting of a scatter (PERF.122, done) | Yes | none |
| `planetgen galaxy` | Sector generation in nine modes (random start with neighborhood, shell, block, center sector, band, wedge, course), the neighbour link, the three scatters, the backfill, the population pass and the settle, with no option to leave one out (PERF.123) | No: eight jobs in one run | PERF.123 |
| `planetgen sector` | Generates sectors and nothing else (`--population` is gone, PERF.124, done) | Yes | none |
| `planetgen phenomenon` | Generates one phenomenon and saves it (`--no-settle` is gone; the settle is `planetgen settle --sector-id ID`, PERF.124, done) | Yes | none |
| `planetgen system` | Generates one system and prints or writes it | Yes | none (registered by PERF.118) |
| `planetgen population` | The population pass: species, civilizations and polities (`--rescan`); `--territories-only` is gone (PERF.125, done) | Yes | none |
| `planetgen territories` | The territories pass: which polity owns which system (PERF.125, done) | Yes | none |
| `planetgen bring-to-epoch` | Brings a stored galaxy to the running generator epoch: generates again the sectors whose fill is behind, makes the stale species again, runs the population pass (GEN.764) | One purpose, a fixed sequence of its own like `benchmark` | none |
| `planetgen refresh-parity` | Makes the sector parity file current for the whole database: marked and never-sealed sectors, sealed sectors the file lacks, members whose sector is gone; `--compact` drops the file's dead space (DB.46) | One purpose, a serial sweep with its own locking of the parity file, no generation | none |
| `planetgen repair-db` | Rebuilds the sectors that no longer hash to their stored checksum from the parity file (`--shell`, `--sector` limit it, `--dry-run` writes nothing); makes a damaged parity file again from the database (DB.47) | One purpose, a serial pass over the sectors in scope, no generation; writes only the covered columns of the damaged sectors' rows | none |
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
