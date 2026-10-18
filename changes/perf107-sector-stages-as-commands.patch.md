**Headline:** The sector stages of `planetgen galaxy` run on the work queue as commands.

### Changed
- **`generate-start-sector`, `generate-neighborhood` and `generate-sectors` are commands (PERF.107).** The starting sector, the neighborhood and every other sector `galaxy` generates (random start, shell, band, wedge, span, course) are cut into static chunks of consecutive sectors, about 4,000 at most, queued through the wrapper and run by the workers; the console still reports each sector as it is saved, and the galaxy holds the same stars, planets and moons at any worker count, as it did before. A failing sector now stops the run with the unit's error, and the sectors already saved stay saved, so a rerun resumes.
- **`link-neighbours` and `settle-sector-paths` are commands (PERF.107).** The neighbour link and the sector paths run as one unit each on the sectors the run created, found from the database's clock at the start of the run.
- **The galaxy job tree no longer lists a queue of sector tasks.** The sectors are units of the run record on Redis. `docs/design/redis-and-the-work-queue.md` describes them.

TODO items resolved: none (the sector and link and settle stages of PERF.107; the backfill is PERF.108). Added: none.
