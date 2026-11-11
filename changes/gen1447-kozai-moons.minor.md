**Headline:** A captured moon inclined to the Laplace plane at its distance now has its eccentricity pumped by the star (Kozai-Lidov), and a moon the cycle drives into its planet is lost; moons that formed in place are not touched, and nothing is stored yet.

### Added

- `physics/kozai_bodies.py`: the orbit-averaged quadrupole cycle of a body that has a distant third one, integrated (not fitted) for the maximum eccentricity, the period and the state at the system's age; the quenches that stop it (the hierarchy, the inclination window of 39.2 to 140.8 degrees, a Kozai time longer than the age, general relativity, a faster precession from the siblings or the planet's bulge) each name their cause; the octupole cap of 0.98; and the Laplace plane of a planet (`laplace_tilt_rad`, `moon_free_inclination`) so a moon is measured against the plane it precesses about and not against the planet's orbit.
- `StarSystem._kozai_pumping` runs after the secular exchange of the settling to the final age: a surviving moon takes the eccentricity of its place in the cycle at the star's age (which the tide then wears down from there), a lost one is removed from its planet, and the star's pull on every moon is noted as the effect ledger's secular rate.
- Stand-in rows KZ1 to KZ6 (GEN.725, section 6.262).

### Changed

- Over 60 seeds, 3 of 167 moons screened cycle and 2 are lost. The moons that formed in a disk, an impact or a ring keep the inclination to the equator they were drawn with as their free one, so none of them cycles; the planets (waiting for GEN.1500) and the belts and comets (GEN.1446) are not in this slice.

generation-output: changed
generation-affects: fill
