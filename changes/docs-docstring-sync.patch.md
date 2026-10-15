**Headline:** The code's own docs name the modules and the shell, band and wedge grid as they are now.

### Changed
- About 270 mentions in docstrings and comments of modules that no longer exist by that name (`queryDb`, `_db`, `sectorGen`, `brightStars`, `adminAuth` and the like) now name the real module, and the scatter's "layers" and "rings" are shells and bands. Only comments and docstrings change, plus the sector grid comment in `schema.sql`; nothing behaves differently.
- The `planetgen.galaxy` package docstring lists its modules and says what each of the map's modules does on the shell grid.
