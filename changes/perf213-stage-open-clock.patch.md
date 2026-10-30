**Headline:** The link, backfill, settle and block totals stages of a generation run now work on the sectors saved since the sectors stage opened, not since the run was built.

### Changed

- **Sectors of the run, by the clock of their stage (PERF.213).** A galaxy or sectors run used to take "the sectors saved since the run started" from the database's clock read when the run was built, before it waited for the machine and before the sectors stage opened. The link, bright-star backfill, sector path settle and block totals stages now take the clock when the first sectors stage opens (a step parameter worked out at that moment, so a stage of the same run can take it from another), and a run carried on with finish-run keeps the clock the stopped run took. The same sectors are linked, backfilled and settled, so nothing stored changes; a run that waited in the queue no longer takes in sectors another run saved meanwhile.
