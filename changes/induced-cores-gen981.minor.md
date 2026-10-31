**Headline:** A moon's iron core and a rocky planet's core and ocean now answer the changing field around them with an induced magnetic moment, so a planet with no dynamo of its own still carries a small one from its star's field.

### Added
- **Induced moment of cores and of a planet in its star's field (GEN.981).** A metal core answers a changing field as a conducting sphere (the exact response, 1 for a core many skin depths across; conductivity 1e6 S/m, solid or liquid). A moon's core adds to its ocean's moment in its parent's turning field, past the screening of the ocean. A rocky planet is excited by its star's field at its orbit (5 nT at 1 AU scaled by the root of the wind's pressure, half of it perpendicular to the spin, at the planet's synodic rotation); its core and its ocean answer, and the moment is stored in `induced_moment_a_m2`, which was filled for moons only (a dry Earth-mass planet at 1 AU gets about 5e17 A m2, a dynamo planet 8e22; a close planet of a flaring M dwarf gets more). The property sheet now shows the induced moment on planets. Four stand-in rows (IN1 to IN4) are registered.

### Changed
- Generator epoch bump (numbered at merge): the induced moment of every rocky planet (new) and of every moon with an iron core (adds the core) changes; nothing else reads it yet.
generation-output: changed
