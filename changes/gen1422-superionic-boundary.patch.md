**Headline:** The deep ice layer of an ice giant is superionic, an ionic fluid or an ice by a measured table of pressure and temperature, not by a single pressure.

### Changed
- **Superionic water and ammonia come from a boundary table (GEN.1422).** `physics/superionic.py` reads `physics/data/superionic_boundary.json` (floors and ceilings of water from the diamond-cell work [R5168] and the shock work [R3081], ammonia and the metallic limits from [R5169]). The ice mantle of an ice giant takes its state, and its water becomes superionic ice, from the pressure and the temperature of the structure model's adiabat at its bottom; the single 100 GPa line is gone, and the weird-worlds line G4 is built for the superionic part. Stand-ins SI1 to SI4 are in section 6.233. This changes ice-giant layer states and materials, so it takes a generator epoch (the Train Manager assigns it).

