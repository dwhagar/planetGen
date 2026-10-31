**Headline:** Every stage with units now has the workers its memory allows, sector generation is capped at what a sector worker really takes, and a run's report names the largest worker of each stage.

### Changed
- **Sector stages are limited by memory (PERF.226).** Measured on a populated database (a two-shell galaxy, the starting sector and 12 of its 107 neighbours, 4 workers): one worker of the sector stages holds up to 1,487 MB resident. `workers.SECTOR_WORKER_MEMORY_MB` (1,600) sizes `generate-start-sector`, `generate-neighborhood`, `generate-sectors` and `generate-unplaced-sector`; they had the ordinary 300 MB guess, so 36 workers could take 54 GB of a 16 GB host. `stage_worker_cap` now covers every command with units: the measured table (`STAGE_WORKER_MEMORY_MB`, renamed from `MEMORY_HEAVY_KINDS`), the unmeasured list (`UNMEASURED_KINDS`, 900 MB each until measured), and nothing for a command that is one unit or inline. A test holds every command to one of the three.
- **The run report names the memory.** Each unit records the most resident memory its process held; the stage keeps the largest. The run's log says it beside the stage's figure (a warning when it is over the figure), the stage's finish line carries "(largest worker N MB)", and the admin queue page has a "Largest worker" column.

No generator output changes (the work seed does not depend on the worker count).

generation-output: unchanged

TODO items resolved: none (PERF.226 stays open: the unmeasured stages and the worker footprint itself).
