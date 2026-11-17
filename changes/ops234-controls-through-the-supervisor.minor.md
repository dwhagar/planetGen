**Headline:** Pause, Stop All, cancel and clear are requests to the supervisor, which stops the job processes it forked itself and holds its hand-outs while the queues are paused.

### Added
- **Controls through the supervisor (OPS.234, first slice).** A request list in Redis and a reply list per request (`queue/controls.py`: `ask`), and a supervisor duty (`controls.Controls`, added where `cli/supervisor.py` builds its duties, before the launcher) that takes the requests each turn: `pause`, `resume`, `stop_jobs` (a run's jobs) and `stop_all`. For the ledger it cancels the rows of the jobs not started, asks each job process it forked to stop between two pieces through the pipe it holds (the pieces done so far kept), kills one still alive after 10 s, closes its running attempt `lost` with the kind `canceled`, and ends the run's row. Stop All's answer comes when the last job process is gone, so the process scan that follows it, which stays as the second line, has none of the engine's to find. The scan could not find a forked job in any case: its command line is the template's, not a worker's.
- **A pause holds the hand-outs.** While the queues are paused by hand the launcher only settles the jobs that end and the rows a dead supervisor left; it forks nothing new until Resume. A solo hold for one restarted unit is not a pause of the hand-outs. Before, the launcher did not read the pause at all.

### Changed
- `runcontrol.cancel_run` and `clear_run` leave the ledger's half to the supervisor (`controls.stop_run`): a request when one is alive, the rows ended directly when none is (no job process of the ledger can be running then). The Redis record and the RQ units are as before.
- `stopall.stop_all` asks the supervisor to stop every job process and reports them (`jobs` in the summary, and a line in `planetgen queue stop-all`); the Queue page's Pause and Resume go through the supervisor when one is alive, and write the flag themselves when none is.

generation-output: unchanged
generation-affects: none
