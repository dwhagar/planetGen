**Headline:** A black hole with no accretion disk no longer carries planets or belts that its mass alone used to give it.

### Fixed
- **GEN.1648: planets around a black hole with no disk.** `StarSystem._generate_planets` sent a millisecond pulsar and a planetary nebula's pair to their own rules and every other remnant to the generic architecture, which places an inner chain from the mass alone (16 of 40 disk-less stellar holes got up to 7 planets at 0.1 to 0.6 AU, 2 of those also a belt). No disk survives the progenitor and no planet of a black hole is known, so `remnant_orbits.is_diskless_black_hole` now returns no planets and no belts for such a hole unless the run asks for planets (`PLANETS` true). A black hole with a disk (0.1 percent) and every other remnant keep their rules. The test that pinned seed 0 now covers 40 seeds.
generation-output: changed
generation-affects: plan, fill
