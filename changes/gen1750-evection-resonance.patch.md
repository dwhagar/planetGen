**Headline:** A moon that the planet's bulge brings into the evection resonance (its orbit turning once a year of the star's) is now caught and held there with a high eccentricity, or passes, by how fast the resonance sweeps it.

### Changed

- GEN.1750: `physics/evection.py` (new): the resonance's distance from the planet's spin (`kozai_bodies.bulge_precession_rad_s`, one shared bulge precession), the sweep ratio `rho` of the resonance's coupling against its rate across the moon, and the capture rule (caught at `rho` 1.58 or under, from the likely and unlikely rates of the first analysis, 0.54 and 5.4 for the Earth).
- `physics/moon_orbits.py` (`evolve`): each prograde moon gets an evection event; at the crossing the moon is caught or passes, and a caught moon is held at the resonance, which moves in with the planet's spin as its tide still brakes the planet (`Fate.held`). `physics/moon_history.py` gives a caught moon the eccentricity it raises (`min(0.5, 1 - 2 R / a)`, kept by `settle_orbit` as a floor) and makes its tide again. Stored values change for the moons the resonance sweeps: their distance, eccentricity and tidal heat.
- Stand-ins EV1 to EV3 are in GEN.725 (section 6.302); the exit from the resonance is not modelled.

generation-output: changed
generation-affects: fill
