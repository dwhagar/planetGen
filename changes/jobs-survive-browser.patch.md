### Changed

- "Generate the neighborhood" on the Sector Map now starts a background job (followed on the Generate page and Admin, Queue) instead of running inside the page request, so closing the browser no longer matters (ADM.11).

### Fixed

- Two admins starting a job at the same moment could both start one: the job lock was briefly empty and the second caller cleared it as stale. The lock now appears with the job id in it, and only one caller clears a stale lock (TEST.40).
- A job id drawn twice in the same second gets a fresh one, and pruning old jobs never removes the running one (TEST.41).
