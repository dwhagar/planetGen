**Headline:** The block light totals stage is now held to the workers the host's memory holds, and a stage's cut into units runs with the other workers gone, so a full-galaxy run is less likely to be killed for memory at the end of the block lights.

### Changed
- **The block lights have a worker memory (PERF.226).** `refresh-block-lights` is exempt from the unit-length ceiling (`whole`) but is cut into one unit per run of shells (90 on the full galaxy), and it had no figure in `workers.STAGE_WORKER_MEMORY_MB`, so the run's 15 workers were all kept (15,128 MB of a 15,992 MB host when the run of 2026-11-08 was killed). It now has the lean figure (700 MB), 11 workers on that host.
- **The workers beside the opener leave before a stage is cut (PERF.226).** The unit that completes a stage cuts the next one in its own work horse; the other workers only waited, each holding what its last unit left. The cut now asks them to finish and exit first (`workers.shed_bulk`), and the run's `ensure` does not start workers while a cut is in progress; it starts the next stage's workers once the cut is done.

Output unchanged.
