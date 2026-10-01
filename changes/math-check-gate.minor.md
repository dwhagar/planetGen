### Added
- **A math check that runs first (TEST.63 to TEST.67).** A new module,
  `stellarObjects/mathCheck.py`, checks the generator's math against 49
  fixed answers before anything trusts it: known values from real
  astronomy (the Sun's luminosity, lifetime and temperature, Earth's and
  Jupiter's orbits, Earth's Hill sphere, the habitable zone and snow line,
  white dwarf sizes, the Sun's Schwarzschild radius and galactic orbit,
  Holman & Wiegert's stability limits, the Kepler and Barker equations),
  identities that hold for any input (unit conversions, constants that
  agree with each other, energy conservation around an orbit, the sector
  grid), and seeded draws from every sampler (star masses, ages, sector
  counts, planet sizes and classes) checked against their intended shares.
  Each check says where its expected value comes from. It takes under a
  second. The test suite runs it before any test and stops if it fails,
  CI runs it as its own first job, and the website runs it at startup and
  shows admins a warning if it fails. `python -m stellarObjects.mathCheck -v`
  prints the report.

### Fixed
- **The speed of light disagreed with the light-year.** `SPEED_OF_LIGHT_M_S`
  was 2.998e8 m/s while the light-year used the exact 299,792,458 m/s; it
  is now exact too, so black hole and quasar event horizons come out
  0.005% larger.
