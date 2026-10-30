**Headline:** A moon orbit module follows a planet's moons drifting out, spiralling in or falling into a ring under the planet's tide, ready for the tide chain to use.

### Added
- **`physics/moon_orbits.py` (GEN.958, stage 3a of the tide chain).** A pure module that integrates a planet's spin and all its moons' orbits together under the tide each moon raises: a moon of a fast-spinning planet drifts outward (the Moon at 3.7 cm a year), a moon of a slow or tide-locked planet or a retrograde moon (Phobos, Triton) spirals in, and one that reaches the Roche limit becomes ring material and gives its angular momentum to the planet. It also returns the heat the moons' tides put into the planet. Angular momentum is conserved. Nothing calls it yet; GEN.959 wires it into the moons, the rings and the planet's heat.
- **`spin_torques.Torques.softness`:** an optional smoothing of the constant-Q torque's sign flip, used by the integration.
- **Stand-in rows TM1 to TM4** in the gap-filling plan (section 6.25) and the stand-in list.
