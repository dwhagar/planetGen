### Fixed

- The size estimate before a bulk run is closer to what the run adds (PERF.26). Fills measured on MariaDB add 37 to 52 KB per star system, so the default before anything is measured is now 48 KB instead of 70 KB. The size measured after a run leaves out the galaxy-wide tables a plan or the bright-star scatter writes, and on MySQL 8 it reads the tables' current sizes instead of a copy cached for a day.
