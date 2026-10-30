**Headline:** The speed a rogue planet leaves its star with is now taken from published scattering simulations instead of a guess, and a calculation for how often its moons stay bound exists; no world uses either yet.

### Changed
- **Rogue ejection speed (GEN.1113).** `physics/rogue_kick` now gives the speed of an ejected planet from the scatterer's mass and orbit, as the two N-body papers read in full (Huang and Lai 2025, Kramer and Profumo 2026) agree: about 1.2 km/s for a Jupiter at 5.2 AU, 3.5 km/s for ten Jupiter masses, falling with distance, independent of the star and of the planet thrown out. The earlier energy-step guess (median 6 km/s for a Jupiter) is gone. The median over giants of 0.3 to 10 Jupiter masses at 1 to 30 AU stays inside the 0.7 to 10.7 km/s of Sayed and others. The tail that eccentric giants and three-planet systems add is not modelled yet.

### Added
- **Moon survival (GEN.1113).** A moon of an ejected planet stays bound if it lies inside half the Hill sphere the closest approach leaves it (Huang and Lai Eq. 4); the closest approach is drawn from their simulated median and spread, and the share of ejections that keep a moon at a given distance is given. Registry rows RK1 and RK3 to RK8 (section 6.87) record the numbers, and the report `docs/design/research/rogue-ejection-kick.md` the reading. Nothing reads any of this yet (GEN.1111).

generation-output: unchanged
