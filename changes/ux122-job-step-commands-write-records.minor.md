**Headline:** A generation job's log now ends each stage with a finish line, says why a failed run failed, and the data, orbit and reset steps write their lines into the job's records file.

### Changed
- **Stage finish and failure lines (UX.122).** A generation run writes `=== Stage N of M finished after 31 s ===` when a stage ends and `=== Stage N of M failed after 31 s: reason ===` when it fails, on the command line and in the job's records file, the same wording the job runner uses.
- **Data, orbit and reset commands (UX.122).** `planetgen datasets`, `planetgen orbits` and `planetgen reset` report through the logging layer, so what they say in the Generate page's log is also in the job's `records.log`; warnings and errors still go to stderr and the text is unchanged. The interactive prompts and the listing commands (`query`, `config`, `backup`, `queue`, `maintenance`, `accounts`) print as before.
- **Activity log (UX.122).** Its line keeps its form on purpose, because the fail2ban filters match it; the logging standard records the decision.
