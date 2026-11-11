**Headline:** After the stuck action stops a unit and queues it again, every other task is held so the restarted unit runs alone, the Queue page says so and who holds it, and the hold lifts by itself when the unit completes or fails (OPS.145, heartbeat part 8).

### Added
- **A solo hold on the one pause flag.** `queue/pause.py` keeps its flag (`planetgen:queue:paused`) and gives it two fields, `reason` and `solo` (the unit whose job is let through). `pause.hold_solo` is called by `unitstuck.Stuck` when it queues the stopped unit again, at the front of its queue (`redisqueue.requeue_job(..., at_front=True)`); `pause.blocks` is what the workers ask. A worker that takes any other job while the hold is on puts it back at the end of its queue and waits (`cli/worker.py`), so nothing new starts beside the unit; units already running finish or stop at their next batch. A pause made by hand is never replaced by a hold, and the follower never lifts a pause it did not make.
- **The hold lifts by itself.** The follower lifts it when the unit is complete, when the run fails (including the unit's second stuck verdict) or ends, and `pause.state` drops it when the unit's job has ended or is gone (a follower that died). An admin can lift it from the Queue page: it reads "The queue is held by run ...: solo run of unit N (stage S) after a restart", with a Lift the hold button (`queue-resume`). `GET /api/admin/work` status gains `paused_reason` and `solo_unit`.
- **The memory gate gives the restarted unit the host.** `memorygate.Claim.solo`: the unit asks for the whole memory budget of the lanes that are not the site's (`memorygate.solo_mb`: the host's memory less its own reserve and the places kept for the interactive and long lanes), waits until no claim of the bulk or maintenance lanes is held, and no bulk or maintenance claim is admitted while it is held. The site's lanes are not held.
- **The scheduler treats the solo run as a turn that does not rotate.** The solo ticket goes before every other (`scheduler.SOLO_WEIGHT`), holds the bulk lanes back while it waits, and its time is not charged to its kind (`memorygate.release`), so the kinds that waited are not behind it when the hold lifts and keep the waiting they did.

### Changed
- `pause.state()` returns `reason` and `solo` too (`None` for a pause by hand).
- No star, planet or stored number changes.

generation-output: unchanged
generation-affects: none
