**Headline:** The galaxy now has one pure module for its pull and its rotation curve, built from the galaxy's own masses.

### Changed
- **The galactic potential (GEN.338).** `planetgen/galaxy/potential.py` gives the potential, the pull (finite at the centre) and the circular speed of a generated galaxy from Hernquist spheres (the bulge, the stellar halo, the nuclear cluster), Miyamoto-Nagai disks (thin, thick and nuclear) and an NFW dark halo, with the central black hole softened at 1 pc. Every mass and scale comes from the galaxy's density components and its mass chain, so the rotation curve is a result, not a constant; a component of mass zero is legal. The matching of the density model to those three profiles is registered as stand-in rows GP1 to GP3. Nothing reads it yet: the acceleration kernel (GEN.342), the star stepper (GEN.340) and the velocity seeding (GEN.339) do, so generation output does not change. The dark halo's parameters moved into `enclosed_mass.dark_halo` so the enclosed mass and the potential use the same halo.

generation-output: unchanged
