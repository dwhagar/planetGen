### Added
- The Work Queue's Jobs list shows a "Time left" column for running jobs, from the job's measured pace or the estimate a Generate page job publishes; it stays blank when none can be told.
- A "Clear finished jobs" button on the Work Queue page (and `POST /api/admin/work/clear-finished`) deletes every finished job's record in one step; running jobs stay.

### Changed
- The Jobs list is minimal: job, status, started, duration, time left, progress and its controls. The Generate page's recent-jobs list drops its "By" column, is paged like every list, and points to the Work Queue for every job, including ones started from the command line.
