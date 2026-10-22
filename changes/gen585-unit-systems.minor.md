**Headline:** The unit registry now has unit systems, and a physical constant's number in AU, year and solar mass is worked out from its dimension.

### Added
- **Unit systems** (GEN.585): `unit_registry.UnitSystem` is four registry units (length, mass, time, temperature);
  `si` and `au-yr-msun` are declared. `constant_in("G", "au-yr-msun")` (also c, sigma, k_B, hbar), `in_system` and
  `system_factor` work a number out from its dimension and the system's base units, so a typed shortcut such as 4 pi
  squared can be generated and checked. The computed G for the solar mass in the Julian year is 39.4769, 3.8e-5 below
  the typed 39.4784 of `kepler.AU3_PER_YR2_PER_SOLAR_MASS` (the Gaussian year's); the difference is stated in the
  module header and pinned by a test, and GEN.587's sweep decides which stands.
- Kinds force and mass flow, with newtons, dynes, kg/s, g/s and solar and Earth masses per year; a test that every kind
  has a five-integer dimension and that the derived kinds the code uses have theirs.
- The column registry (DB.76) lost its local `km3`: the unit registry now holds it with the area and volume kinds.
