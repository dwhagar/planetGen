**Headline:** Every planet and moon now says what its day and night are like (an ordinary cycle, tidally locked, a spin-orbit resonance, hardly rotating, polar, or no solid surface), how many hours of night it has and how many solar days make a year.

### Added
- **Day-and-night state (GEN.1097).** A new `physics/day_night` module names the state of a body's day and night from its spin, its orbit and its polar circle: Earth and Mars are an ordinary cycle, Mercury a spin-orbit resonance (the star returns every two years, 0.5 days in a year), a world that turns once an orbit is tidally locked, Uranus is polar, a giant has no solid surface. The state, the hours of night in a solar day at the equator and the solar days in one orbit are stored on every planet and moon (`day_night_state`, `night_hours`, `days_per_year`, schema v188, migration 0188) and shown on the property sheet. Belt bodies, stars and rogue planets wait for their own storage items. The registry has a new row (DN1).

### Changed
- `physics/daylight` returns the three new fields with the ones it already had; a gas giant is `no solid surface` whatever its spin.

generation-output: changed
