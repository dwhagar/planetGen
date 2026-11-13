**Headline:** a belt's mass now pulls on the planets of its star, turning their orbits as it would, so the eccentricity and tilt a planet has at the star's age depend on the belts beside it.

### Changed

- `physics/secular_chain.py` (GEN.1464, secular side): a belt is cut into eight equal-width rings (`belt_rings`, mass spread over the annulus at one surface density) that add their pull to each planet's apsidal and nodal rate. The rings keep a flat circular orbit and are not moved; a ring within 10 percent of a planet's axis is left out. `matrices`, `precession_rates` and `evolve` take the rings (`fixed_axes_au`, `fixed_masses_kg`), and `planets.secular_exchange` passes every belt of the star (`planets.belt_rings_of`). `belt_forcing.planet_rate` takes them too. Registry BR1 to BR3 (6.244); SC2 no longer says "no belt".
- Census, 150 systems: 88 have a belt with mass and in 84 of them a planet's orbit at the star's age moves (median change in eccentricity 0.027, 90th percentile 0.077, at most 0.13; tilt median 1.9 degrees, at most 6.6). A small belt moves a planet as much as a heavy one, because a small change of the turning rate shifts the phase of a slow cycle that has run for gigayears. Systems with no belt do not change at all. The exchange's outcome word changes in 8 of 150 star systems (4 become 'outside the theory' or leave it).
- The tide terms of belts, asteroids and comets (part 2 of GEN.1464) wait for GEN.1446 on main.

generation-affects: fill
