**Headline:** One controller now watches the host's memory: on lasting pressure the supervisor stops admitting bulk and upkeep jobs, and if the host stays starved it asks the youngest bulk job to stop at its next piece and puts it back, with no retry used.

### Changed
- **One memory controller (OPS.231, D6 of the job engine).** `queue/pressure.py` reads the kernel's memory pressure (PSI) or, where there is none, the host's available memory against the gate's reserve, and runs a hold, yield, cooldown and release cycle in the supervisor. While held, the launcher and the light jobs' duty admit only the site's lanes and say why the rest waits. While starved, it asks the youngest bulk job (a command that resumes before one that does not) to stop at its next piece through the job pipe (`Launcher.yield_youngest`). The row goes back to *ready*; the attempt is stored as `lost` with `error_kind="stopped"` (`ledger.STOPPED`) and does not count against the retry limit. The edges are stand-ins that the job engine's H4 learns.
- **The Queue page** words a yielded attempt as a stop on request, not as a loss (`jobview.attempt_words`).
- A unit that is stopped by request ends its attempt with `error_kind="stopped"` (`units._do_unit`, a one-line change).

### Removed
- **The adaptive worker limit** (`queue/adaptive.py`, OPS.148): the lowering and raising by waits at the gate, its start limit per stage kind, `memorygate._adapt_to_wait`, the `adaptive:settled` profile entry and the log of recent changes.
- **The lend rule and the per-lane places** in `queue/workers.py`: `shed_bulk`, `staying_bulk_workers`, `maintenance_share`, `reserved_places` and the lent-worker keys; `wrapper` no longer asks the bulk workers to go before a stage opens. Each of the lane's workers is held by the host cap and the stage's own cap alone, and memory is the controller's.

generation-affects: none
