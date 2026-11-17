**Headline:** A unit the follower finds stuck is stopped through the supervisor (the job process killed, its attempt closed `stuck`), and the reset command and the Stats page's worker count read the ledger and the job processes, not RQ.

### Changed
- **Stuck units of a ledgered run (OPS.234, second slice).** `unitstuck.Stuck` asks the supervisor to stop the unit (`controls`, the new `stop_unit` request): the job process it forked is killed and its attempt is closed `stuck` with the reason and the stack the follower had written (a file the attempt's row names). The launcher settles the end of the process and hands the unit out again, and the follower clears the stale beat when it sees the row move on; a unit stuck a second time fails the run as before. Before, the follower looked for an RQ job a forked unit does not have, so the first stuck unit failed its run at the next poll. A dead beat of a unit whose row is `ready` again is a retry waiting, not a dead unit. Runs with no rows in the ledger keep the RQ path until OPS.271. The solo hold (OPS.145) is RQ's: the launcher's hand-outs do not read it.
- `planetgen reset-jobs` has the supervisor stop every job process it forked and ends the ledger's rows (a line says how many), and reads the queues and the workers of an older run as Redis keys, not through RQ.
- The Stats page's worker count is the job processes the ledger says run and whose process exists (no worker registers in Redis any more); a control database that does not answer says "workers unknown".

generation-output: unchanged
generation-affects: none
