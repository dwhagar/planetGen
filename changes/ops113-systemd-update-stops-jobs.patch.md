**Headline:** An update started by the Update button or the monthly timer now stops the running jobs and goes on, instead of refusing.

### Changed
- **OPS.113: a systemd-started update failed whenever the warm-up or another job was running.** Started by systemd, `update.sh` now cancels everything it lists (the warm-up, queued jobs, tile and orbit work) without `--yes`, leaving each command-line run a stop reason that says an update started by systemd needed the machine. The one thing it never cancels is generation (a Generate job, a generation run, a job on the generation queue): with that running it changes nothing and the next firing or press tries again. Runs by hand and from cron keep the question and `--yes` rules.
