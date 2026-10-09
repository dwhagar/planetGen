### Added

- Every rotating object stores a spin axis and an axial tilt (GEN.104, schema v68): stars, planets, moons, comets, rogue planets, interstellar comets, and standalone black holes and neutron stars. Stars, comets, rogue planets and interstellar comets also store a rotation period. Cool stars spin by gyrochronology, hot stars by a log-normal speed held under breakup, small bodies never faster than the 2.2-hour spin barrier, and black hole spin follows Beta(1.4, 3.6).

### Changed

- Planets now tidally lock to their star by the same rule moons use; a locked body turns once per orbit, upright. Seeds give different systems than before.
