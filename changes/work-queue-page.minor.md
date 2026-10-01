### Added

- An admin Queue page (`/admin/queue`, in the settings menu) to view and
  manage generation jobs (ADM.10): workers active and the server's load
  as "x / x / x" (CPU percent over 1, 5 and 15 minutes on Windows), the
  job trees with timing, progress and ETA on every node, and Pause,
  Resume, Cancel, Retry and Delete on any job, part of a job or failed
  sector. A paused job finishes its running tasks and stands by without
  holding up other jobs; "Pause the queue" stops every job from starting
  or taking another task until it is resumed. Each control is confirmed
  first and written to the activity log.
