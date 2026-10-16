**Headline:** The orbital update does nothing while an update is running, and a cancel never leaves positions half moved.

### Changed
- **The orbital update refuses to run during an update.** `python -m planetgen.cli.orbits` (the nightly timer, a Generate job's step, a queued run) checks whether an update is running, or has left the site closed after failing, before it touches the database. If so it moves nothing, leaves every epoch, due time and the clock as they were, says why and exits with status 75 so the next timer firing simply tries again.
- **A cancel finishes the update of the positions first.** Once the first stage has begun, a stop request (the update's cancel, Ctrl+C) lets the run finish moving, refiling, refreshing containment, nearest systems and sector paths and recording the clock, then exits, so stored positions are never left half moved. A cancel before anything moved stops at once, and a second Ctrl+C stops it regardless.

TODO items resolved: OPS.54 / Added: none
