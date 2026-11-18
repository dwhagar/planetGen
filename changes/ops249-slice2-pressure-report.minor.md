**Headline:** `python3 -m planetgen.cli.pressure report` says what the memory controller's record of the host's pressure holds and the edges it gives, and the record now keeps how slow the database was each minute.

### Added
- **OPS.249 slice 2.** `planetgen.cli.pressure report` prints the number of minutes recorded and their span, the memory kills no pressed reading came before, the holds that were not needed (against what the learning needs before it moves an edge), the minutes that carry a probe of the database and the edges the record gives now; `--csv FILE` also writes every recorded minute, and the kills beside it, so a gate keeps its record with its results.
- Each minute of `pressure_minutes` also holds the slowest probe of the database in that minute and the calm figure it was judged against (control schema v26, two new columns added by the update), the record the factors of the database signal are to be learned from.

### Changed
- `docs/reference/cli.md` no longer lists the `worker` command among the modules without database options.
