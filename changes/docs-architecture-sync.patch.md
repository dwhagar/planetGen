**Headline:** The architecture document describes the code as it is now, module by module.

### Changed
- `docs/design/architecture.md` has a code map built from the real `src/planetgen/` packages (generation, physics, galaxy, population, names, describe, util, queue, db, admin, cli, wiki, api, web) and the current static files and scripts, in place of the old flat `stellarObjects` map.
- The generating flow now follows the shell, band and wedge grid, the two scatter passes and the steps after a run; the serving flow describes the Galaxy Map as it is while its block ladder is rebuilt on shells; the database, login, job and release flows use today's module names, and the CI row reflects the scheduled, manual-only workflow.
