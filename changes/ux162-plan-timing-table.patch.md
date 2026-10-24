**Headline:** The plan no longer prints its table of where the time went unless asked.

### Changed
- **The "Where the plan's time went" table is printed only with `--timing` or `--debug` (UX.162).** Every galaxy plan ended with a step table made for finding slow steps (PERF.163). A plain plan's text is back to what it was before; `planetgen plan --timing` prints the same table.
