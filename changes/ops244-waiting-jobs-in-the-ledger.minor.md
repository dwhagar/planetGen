**Headline:** A Generate job that waits for the running one is now rows of the job ledger, so the waiting line survives a web server restart and two processes start it once.

### Changed
- **OPS.244 slice C: a waiting Generate job is a waiting run of the ledger, not a `queued-<ticket>.json` file and a job entry in the run line.** `jobs.queue_next` writes the run with every row waiting (`pagerun.queue`); `jobs.promote_queued` claims it with one conditional update (`pagerun.claim_queued`), starts it as the job (`jobs.start_job(..., queued=...)`) and puts it back at its place with the reason when it could not be started; Discard deletes the run.
- **The run line holds runs only** (`queue/runline.py`: the job type, its fields and `claim` and `restore` are gone). The strip and the Generate page show both in the order they arrived, and the cap `queue.max_waiting` counts both.
- **The engine design note says where the waiting line stands** (`docs/design/operations/job-engine.md`, stage 5 and section 7).

generation-output: unchanged
generation-affects: none
