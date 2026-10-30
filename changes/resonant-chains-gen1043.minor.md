**Headline:** The generator now has the physics of resonant chains, the rows of planets whose orbits lock in step and tug each other's orbits oval.

### Added
- **`physics/resonant_chains.py` (GEN.1043).** From a system's orbital periods it finds the chains of three or more planets whose neighbouring periods sit near 2:1, 3:2, 4:3 or 5:3 (TRAPPIST-1 comes out as one chain of seven, Kepler-223 of four, Jupiter's three big moons of three, and the Solar System as none), and it gives the eccentricity each neighbour forces on a planet, from the resonance's disturbing function: Io 0.005 against the observed 0.004, TRAPPIST-1's inner chain a few thousandths. It hands back the orbit in the form `tidal_heating` already flexes for moons, so an eccentric chain planet is heated the way Io is. Nothing puts a chain on a world yet; that is the weird-worlds wiring (GEN.1049).
- **Stand-in rows RC1 to RC4** in the gap-filling plan (section 6.92) and the stand-in list.
