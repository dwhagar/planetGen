**Headline:** a planet's eccentricity now settles where its neighbours' pumping and its own tide's wearing-down balance, instead of the tide simply capping what the neighbours gave it.

### Changed

- `physics/secular_chain.py` (GEN.1459): `evolve` takes `damping_per_yr`, each planet's tidal decay rate, and follows `dz/dt = (i A - Gamma) z` for the eccentricity vectors (the matrix exponential, checked against `scipy.linalg.expm`); the inclinations are not damped.
- `physics/planets.py`: `secular_exchange` starts a planet from its birth eccentricity (`_orbit_birth`) with its cold circularisation rate (`orbit_elements.circularisation_time_yr`) as the damping, over the star's age. `settle_orbit` no longer wears a damped planet from its birth again: it applies only the harder tide of a strong one (the heat softening the mantle), as the share it adds to the cold decay. A planet with no birth value, a lone planet, one outside the theory, and every captured moon keep the old rule.
- Registry XD1 to XD3 (6.244); the design table row for settled eccentricity is updated.
- Census, 60 systems (538 planets, 470 damped): the stored eccentricity of a planet changes in 358 of them, in both directions (180 up, 178 down), because the secular cycle has run for gigayears and any change of the start or the damping moves its phase at the star's age; the distribution is what moves: mean 0.0457 to 0.0417, median 0.0225 to 0.0171, for planets inside 0.3 AU (352) mean 0.0305 to 0.0237. No system changes its planet count.

generation-affects: fill
