**Headline:** A Generate page job is now a run of the job ledger, so it survives a web server restart and Cancel works with or without the supervisor.

### Changed
- **OPS.244 slice B: the Generate page's job is rows of the ledger, not a folder of files and a private queue.** The job is a run of kind `page` with one row for the job and one row for each step (`queue/pagerun.py`). The launcher hands a step out like any other row, a forked job process runs it (`forkedjob.run_page_step`, `web/job_runner.run_step`) and follows the step's own process, and the page reads the job's status, step, start, end and error from the rows. `job.json`, `state.json`, `runner.pid`, the `cancel` file, the job's private RQ queue and its burst worker are gone; `output.log`, `records.log`, `progress.json` and `memory.json` stay in the job folder. A step is never run again by the launcher after a restart or a kill; the page's Retry starts the rest.
- **Cancel is a request to the supervisor** (`controls.stop_engine_run`), or the rows are ended in the page's own process when no supervisor is alive; a step whose runner is gone is stopped from the process group named in `step.pid`.
- **The supervisor stays awake while a page job runs** (`supervisor.pages_key`).
- **The Queue page no longer lists a Generate page job among the beating jobs**; its state is the ledger's.
- **The engine design note says where the old queue stands** (`docs/design/operations/job-engine.md`, sections 2 and 7).

generation-output: unchanged
generation-affects: none
