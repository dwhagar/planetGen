**Headline:** Two moons whose periods closed on a 2:1, 3:2 or 4:3 on the way out are now held in that resonance, and moons near such a ratio force an eccentricity on each other that the planet's tide does not wear away, so an Io and Europa pair is flexed by its neighbour.

### Changed

- GEN.1449 (slice 2): `physics/moon_resonances.py` (new): `capture` holds a pair of moons that swept a first-order resonance on the way out (born above its period ratio, ended below it) at `p/q (1 + 0.002)`, sharing the pair's orbital angular momentum by the resonance condition so the planet's spin is unchanged; a chain of pairs (Io, Europa, Ganymede) is held together; a pair whose heavy outer moon would hold the inner one inside its birth orbit is not captured. `give_forced_eccentricities` gives every moon the eccentricity its neighbours near a first-order ratio force on it (a moon already more eccentric keeps its own).
- `physics/resonant_chains.py` (reused): `pair_forced_eccentricities` (two bodies are enough) beside `forced_eccentricities` (chains of three), both through one `_forced`; for the Galilean masses and periods it gives Io 0.0051, Europa 0.0082, Ganymede 0.0007 against the observed 0.0041, 0.0094, 0.0015.
- `physics/moon_history.py`: the capture runs on the survivors after the tide's history and before the clearance pass. `generation/system.py`: `_moon_resonances` gives the forced eccentricity once the orbits are final, before the tide wears one down, and the moon's tide is made again by the final-age refresh. `physics/planets.py`: `settle_orbit` keeps the resonance's eccentricity as a floor. Stored values change for a planet with two or more moons near a first-order ratio: the moons' distances, eccentricities and tidal heat, and what follows from them.
- Stand-ins HP3 to HP5 are in GEN.725 (section 6.300).

generation-output: changed
generation-affects: fill
