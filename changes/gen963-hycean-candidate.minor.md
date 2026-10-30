**Headline:** The code can now say which small planets fall in the Hycean box, the sub-Neptunes that could hold a liquid ocean under a hydrogen envelope.

### Changed
- **The Hycean box (GEN.963).** `physics/hycean_candidate.py` flags a gas dwarf, ice giant or high-gravity world of 1 to 10 Earth masses, 1.4 to 2.6 Earth radii, an equilibrium temperature of 150 to 500 K and a density of 1.5 to 3 g/cc. The audit found the state `hycean` and the envelope survival (`envelope_escape`, the retention library) already in the code; only the box was missing. No class and no stored value moves, so no generated world changes.
- **Stand-in recorded.** Row HC1 in the stand-in registry (GEN.725): the four ranges, as the item quotes them (Madhusudhan et al. 2021 was not read).
