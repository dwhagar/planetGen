**Headline:** Every log line says where it is with one of three words, Stage, Part or Phase, and every line names its component from one table.

### Changed

- A count is worded in one place (OPS.103): `Stage N of M` for the top-level steps of an update, an install or a generate job, `Part N of M` for the runs inside one (the update's two warm-up runs), and `Phase N of M` for a named group inside a run (the warm-up's phases). The warm-up's families are named, and its progress counts them as `family 2 of 2`. `chunk` is only ever a count of pieces of work. `update.sh` words the warm-up runs `part 1 of 2` and `part 2 of 2` (they were `stage 1` and `stage 2`, which collided with the update's own stages), the queue page and the cluster-field messages use the same wording, and the update guard's failed step reads `stage 4 of 9, checking the generator's math`, set by `log_stage` itself.
- The component of a record is one of the names in `logrecord.COMPONENTS` and nothing else (`update`, `install`, `warm`, `tiles`, `orbits`, `generate`, `queue`, `web`, `job`, `activity`, `command`, `lib`, `logs`, `permissions`, `cache-dir`, `swap`, `python-deps`). `Record` refuses another name, `logrecord.component_of` finds the one a logger or a source file belongs to, and the shell scripts take their list from the same table. The queue run page's records now say `queue` (they said `run`), and a library's own log lines say `lib`.
- `test_no_own_log_formatting.py` (TEST.183) fails on a hand-made count in a Python module or a shell script, and on a component named as a literal that is not in the table.
- `docs/design/operations/logging-standard.md` has a "Counting" section and a "Components" table.

generation-output: unchanged
