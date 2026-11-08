**Headline:** a script that measures how far each class quantity moves when one stand-in is moved across its registered range.

### Changed

- `scripts/stand_in_sensitivity.py` (GEN.1311, part one of the class fudge table): reads the registry's constant and knobs cells, overrides one module constant at a time (low end, high end) on a fixed sample of systems, and records how far gravity, mean temperature, pressure, gas shares, ocean cover, the regional, day and night temperatures and the wet bulb move. `rows` lists which registry rows can be overridden and why the others cannot; `screen` runs the pass and writes a report. Nothing it does changes a generated world.
