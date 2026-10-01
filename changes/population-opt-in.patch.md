### Changed

The population pass (species, civilizations, territories) is now optional and off by default. `generate.py sector` and `generate.py galaxy` run it after saving only with the new `--population` flag (this replaces `--no-population`). `install.sh` and `update.sh` (and `install.ps1`/`update.ps1`) ask whether to run it after the database step, defaulting to No after 30 seconds and skipping it with no terminal; `POPULATION=1` (`-Population` on Windows) runs it without asking. `generate.py population` still runs it by hand.

### Added

`GET /api/population` (and `population.population_status`) says whether a population pass has run and whether any species, polity or owned system exists, so the pages and the Galaxy Map's Territories button can hide themselves when there is no population data.
