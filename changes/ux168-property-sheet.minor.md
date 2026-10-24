**Headline:** Every star, planet, moon, asteroid belt and comet now has one property table served at `GET /api/objects/<ref>/sheet`, with a test that fails when a stored column is missing from it.

### Added

- **The property sheet (UX.168).** `planetgen/describe/sheet.py` lists, for each kind, the stored properties with their labels, groups (Identity, Orbit, Physical, Rotation, Composition, Atmosphere, Climate, Water, Magnetic field, Radiation, Habitability, Life, and a star's Zones, Stellar activity, Galactic orbit and Debris disk) and units. `GET /api/objects/<ref>/sheet` returns it for a star, planet, moon, belt or comet: every value is the description's own quantity (so a length, a mass, a temperature or a pressure reads in the reader's unit system) or a number with its unit's symbol, and a property the row holds no value for says "not stored".
- `tests/test_property_sheet.py` checks every column of the `stars`, `planets`, `moons`, `asteroid_belts` and `comets` tables is a sheet row or is on the registry's `NOT_SHOWN` list with its reason, so a column a later change adds cannot be forgotten.
