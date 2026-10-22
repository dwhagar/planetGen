**Headline:** Every unit the project uses is now one row in a single registry, with its kind, exact factor, symbol and allowed range.

### Added
- **`planetgen.physics.unit_registry`** (GEN.568). One row per unit: an ASCII key, a display symbol, a name, the kind it
  measures, the factor to the kind's SI base unit and the offset of its zero (only the temperature scales have one).
  Kinds carry a five-exponent dimension and the values they allow (a mass, duration, power or pressure is at least 0, a
  speed is below c, everything is finite). A temperature difference is its own kind, so 10 K warmer is 10 degrees
  Celsius warmer. Compound units (km/s, AU/yr, pc/Myr) are declared by their parts and tested against them. Relations
  between kinds are data: a mass and its gravitational parameter, a radius and its diameter, a period and its
  frequency, a rotation period and its angular rate. Angular rate and angular momentum are new kinds.
- Tests hold every factor against its definition written out in full, each ladder rung of `util/ladders.py` against a
  registry unit of the same factor, round trips, and the bounds. No existing file changes: the duplicate names in
  `constants.py` stay until GEN.571 sweeps them.
