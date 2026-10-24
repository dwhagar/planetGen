**Headline:** The Queue page lists the runs in Redis with their stages and time left, cancel stops a run at once, and a stuck run can be cleared.

### Changed
- **The Queue page shows runs, not job trees (PERF.99).** `/admin/queue` lists the runs in Redis newest first, each with its stage, how far it is, the time left and whether it has stopped answering. `/admin/queue/<run id>` shows one run's stages and the tail of its log, with cancel, clear and delete. The old job tree pages and the lease clear are gone; API version 6 (`/api/admin/work` returns runs, `/log` is new).
- **Cancel is instant (PERF.99).** Cancelling a run takes its waiting units off the queue and sends the stop command to the started ones in one pass (0.1 s for 4,000 units), and the machine is free at once. A run that has stopped answering, or whose follower process is gone, can be cleared from the page.
- **Run logs are kept (PERF.99).** Each run's event stream is written to `runs/<run id>.log` in the jobs directory when it ends, the last 20 kept, with a download link on the run's page.
- **The Generate page lists a job's runs (PERF.99).** The job page has a Runs panel linking each run and its log, and Cancel job cancels the job's live runs at once instead of waiting for the step to notice.
