**Headline:** The supervisor's hand-out turn is built and tested, used by nothing yet.

### Added
- `queue/launcher.py` (OPS.228, slice A): each turn closes the attempt of every job process that ended (done, failed, or lost with the signal named), reads the `ready` rows of the ledger as tickets, asks `dispatch.dispatch` which go, writes each reason a job waits onto its row, writes the attempt (and its number to the galaxy's `engine_items`) and forks the job from its lane's template. The bulk lane moves to it in the cutover slice, a run at a time.
