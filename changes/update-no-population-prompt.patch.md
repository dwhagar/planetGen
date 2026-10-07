### Changed

- `update.sh` and `update.ps1` no longer ask about or run the population pass, so an update that just wiped the database doesn't offer to fill it; the closing message says how to run `generate.py population` by hand. `POPULATION=1` and `-Population` are gone from the update; the installers keep their prompt (OPS.7).
