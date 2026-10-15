### Added
- Two shared work queues, `planetgen-bulk` and `planetgen-interactive`, each served by its own workers. Workers wait a minute for work before they exit, so a run's queue can be empty for a moment between stages. The API's jobs and the admin's long operations use the interactive queue, so a long run never holds them up.
- An admin can set the bulk worker count from the Queue page, or pick "Ludicrous Speed" (95% of the cores). Lowering the count asks the extra workers to finish their unit and exit.
- A running job can be stopped at once, not only a waiting one.

### Changed
- A worker killed with `kill -9` is replaced the next time work is queued or a job's status is read.
