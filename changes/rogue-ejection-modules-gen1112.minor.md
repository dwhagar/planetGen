**Headline:** The generator can now say what a system throws out and what a thrown-out planet takes with it, from physics instead of rates: masses, air, water, moons, and the planets a birth cluster strips.

### Added
- **`physics/rogue_mass_function.py` (GEN.1112, RP2).** The masses of the bodies a system ejects follow the supply of a collisional cascade (Dohnanyi's 0.83 per decade of mass, up to 1.0), cut at the mass of the planet that throws them out. Its slope sits inside the microlensing measurement of free-floating planets (0.96, range 0.69 to 1.43) without being fitted to it.
- **`physics/rogue_volatiles.py` (GEN.1112, RP3).** The hydrogen a thrown-out planet keeps (the star's XUV over the time it spent lit), the pressure at the bottom of that envelope, whether it is a gas giant (the generator's own giant radius), the water from where it formed against the snow line, and whether it once held an air.
- **`physics/rogue_moons.py` (GEN.1112, RP4).** Whether a moon leaves with its planet (the moon's orbital time against the closest encounter, with the survival shares of Huang and Lai 2026) and the heat it raises there on an eccentric orbit (Debes and Sigurdsson 2007).
- **`physics/cluster_unbinding.py` (GEN.1112, RP7).** The share of planets of a given orbit a birth cluster has stripped from their stars, from the cluster's density and speed history (Spurzem et al. 2009).
- **Stand-in rows RM1 to RM2, EV1 to EV2, MO1 to MO4 and CU1 to CU3** in the gap-filling plan (sections 6.110 to 6.113), the stand-in list, the burn-down register (MO2) and the wiring register. Nothing reads the modules yet; the rogue rows (GEN.1111) wire them and remove the tuning defaults they replace.

### Changed
- The comments on the rogue planet mass bins in `tuning.py` now state the rates the bins actually give (5.8 per star in all, not 6.5).

generation-output: unchanged
