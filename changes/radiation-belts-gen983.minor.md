**Headline:** A moon that orbits inside the radiation belts of its parent planet's magnetic field now takes their dose at its surface, and the habitability radiation score and the survival gear read it.

### Added
- **Radiation belts reach the ground of a moon (GEN.983).** The trapped-particle belts of a planet (`physics/radiation_belts.py`: edges from the magnetopause and the spin, peaks, intensity scaled from Earth's and Jupiter's) are read when a moon is finished: the dose at its orbit, in sieverts by the weight of the electrons and the protons, from the half of the sky its body does not block, cut by its own air with a stopping column of 2 g/cm2. It is stored as `dose_belt_msv_yr` (planets and moons; 0 for a planet, which has its belts above its air, none for a giant) and added to `surface_dose_msv_yr`, so the habitability radiation score, the stored score recomputed under a squeezed heliosphere, the radiation terms on the page, the property sheet and the survival equipment plan (the belt part is cut by the shelter's thickness) all carry it. A moon of a giant, or of a rocky planet with a strong field, inside the belts takes hundreds to thousands of sieverts a year; one outside them takes nothing. Three stand-in rows (RB9 to RB11) are registered; the belt rows of section 6.15 are now read.

### Changed
- Generator epoch bump (numbered at merge): the surface dose and habitability score of every moon inside its parent's belts change.
- Migration (numbered at merge, 0195 provisional): `dose_belt_msv_yr` on `planets` and `moons`.
generation-output: changed
