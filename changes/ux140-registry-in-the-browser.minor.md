**Headline:** The browser reads its unit factors from one generated copy of the unit registry.

### Added
- **The unit registry in the browser** (UX.140). `scripts/build_unit_registry.py` writes `html/static/unitregistry.js` from the Python registry: every unit's kind, symbol, names, factor and offset, every kind's base unit and bounds, and `convert`, `toBase` and `fromBase` with the same arithmetic. A test fails when the file is stale, compares every conversion with Python through Node, and fails when a script declares a unit factor.

### Changed
- `orbitpositions.js`, `systemview3d.js`, `galaxysystem.js`, `orbitclock.js` and `systemreach.js` no longer declare `AU_KM`, `KM_PER_PC` or the seconds in a year, and the nominal 50 AU system of the map uses the exact AU per parsec instead of the rounded 206264.806.
