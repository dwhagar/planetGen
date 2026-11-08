**Headline:** A warm-up or another bulk command can no longer block the site: waiting jobs are now ordered by a fair-share scheduler, with places kept for the site.

### Added
- **A fair-share scheduler on the memory gate (OPS.130, `queue/scheduler.py`).** A job waiting for room is a ticket, ordered by its class (interactive 100, node prefetch 90, long 80, a run's units 40, the maintenance lane 30), plus a rise for each second it waited, less a penalty for the share of recent system time its kind has used. A job that has waited 600 s outranks a new job of any class, so nothing starves. A lower job goes ahead of a higher one that has no room (the host is not left idle), unless the higher one has waited 120 s, when the room is kept for it. The scheduler reads the same load and memory as the gate.
- **Cores are kept for the site's places** (half a core a place, at most a quarter of the host) against the bulk and maintenance lanes, as memory already was.

### Changed
- Every worker and `memorygate.hold` use the scheduler. `PLANETGEN_MEMORY_GATE=0` still turns the gate off.
- Generated output is unchanged.
