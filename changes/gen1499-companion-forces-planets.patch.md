**Headline:** a companion star now pushes the eccentricity of the planets of the star it orbits toward its own periapsis, so a planet in a binary settles on a forced eccentricity instead of ignoring the other star.

### Changed

- `physics/secular_chain.py` (GEN.1499): `evolve` takes `apsidal_extra_rad_yr` and `forcing`, the turn a companion adds to each planet's periapsis and the push `sum B_c exp(i varpi_c)` on its eccentricity vector (`dz/dt = i (A z - f)`); the planets then circle the forced vector `A^-1 f`, or, with the tide of GEN.1459, settle on it. The forced vector matches `belt_forcing.forced_vector` for a lone planet (tested).
- `physics/planets.py`, `generation/system.py`: `secular_exchange` is given the star's companions (`radiation_field.companions_of`, the same companion list and periapsis draw the belt members' forcing uses, `belt_forcing.companion_rates`, `belt_forcing.periapsis_rad`), and a lone planet beside a companion is exchanged too. A companion past `belt_forcing.ALPHA_LIMIT` of a planet forces nothing.
- Registry XD4 added, XD3 trimmed to the moon-moon part (6.244); the design table row for settled eccentricity is updated.
- Census, 50 binary systems (46 wide, 4 close; 550 planets): the distribution barely moves (mean eccentricity 0.0431 to 0.0434, 90th percentile 0.1071 to 0.1023) while 450 of the 550 individual values change, 192 of them by more than 0.01, since a forced vector and a new damping shift the phase of a cycle that has run for gigayears. No system changes its planet count.

generation-affects: fill
