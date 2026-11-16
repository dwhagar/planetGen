**Headline:** With `queue.forked_jobs` on, the supervisor no longer waits a minute for a template to build, takes back rows a dead supervisor left running, and ends a unit whose run is gone.

### Changed
- The supervisor builds a lane's warm template in the background (`TemplateSet(background=True)`); a job whose template is not up yet stays `ready` with the reason `the template ... is being built` and is handed out at a later turn.
- A row left `running` by a supervisor that ended is settled when its attempt's process is gone (or no attempt was taken): handed out once more, then failed. Rows whose process is alive are left to the unit.
- A unit that skipped itself (its run is no longer live, or the unit was already done) has its row ended `canceled`, not left `running`.
- The template treats a reset channel as a closed one and exits.
