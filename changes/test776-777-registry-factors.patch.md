**Headline:** two tests that were red on main pass again, and the browser's clock formatter now reads its time factors from the unit registry.

### Fixed

- TEST.776: `tests/test_photochemical_loss.py` read `pl.YEAR_S`, which the module no longer has, so the file did not collect. It now reads the Julian year from `planetgen.physics.unit_registry`, with no constant of its own.
- TEST.777: `html/static/unitladder.js` (the clock formatter of UX.319) declared the year, the 86400 and the 3600 factors itself, which `test_no_browser_script_declares_a_unit_factor` forbids. It reads the year, day, hour and minute from `unitregistry.js`; the text it gives is the same (the factors are the same doubles).

generation-output: unchanged
