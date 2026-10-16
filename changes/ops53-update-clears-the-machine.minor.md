**Headline:** Before an update changes anything it lists the running and queued work and stops unless you agree to cancel it, and nothing new can start while it runs.

### Changed
- **`update.sh` lists the work that runs first.** A Generate job, a `planetgen` run from the command line, jobs on the shared queues, an orbital update and a Galaxy Map warm-up are printed with their kind, id, start time and what they hold. On a terminal the update asks before cancelling them; with no terminal (a timer, cron) it refuses with a non-zero exit and changes nothing unless you pass `--cancel-running`. A migration run by hand (`python -m planetgen.cli.migrate`) does the same when it has something to migrate.
- **Cancelling goes through each job's own cancel path and the update waits for it.** Queued jobs are removed, running ones are stopped so they end as cancelled, not failed, and the update starts only when none is left. A job that does not stop within 10 minutes is reported and the update stops before it has changed anything.
- **Nothing new can start while the update runs.** The Generate page and the API refuse to queue work and say that the site is being updated, a command-line run that would open a job tree is refused, and the queue's workers exit with status 75. The Galaxy Map warm-up and the orbital update register themselves so a later update can see them.

TODO items resolved: OPS.53 / Added: none
