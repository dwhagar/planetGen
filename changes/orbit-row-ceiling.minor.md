**Headline:** The orbital update has a ceiling on a night's work and warns with the backlog it leaves.

### Changed

- `planetgen.cli.orbits` moves at most `--row-ceiling` rows a run (default 10,000,000, 0 for no ceiling), counted across every table. When more are due it moves the oldest of each table up to a share of the ceiling, prints a warning naming the backlog and where it sits, and leaves the rest due for the next run, which takes them first.
- A due row that is left waiting gets its epoch written first, so no elapsed time is lost.

TODO items resolved: PERF.103. Added: none.
