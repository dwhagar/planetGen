**Headline:** The rogue-planet physics modules lose the functions and constants nothing read.

### Changed
- **GEN.1111 register leftovers wired or deleted.** The seven functions that the rogue physics wrote ahead of their callers and that no part of the plan reads (the donor offset and along-orbit speed of `rogue_drift`; the moon survival draw and probability and the relative speed of `rogue_kick`; the moon survival draw and probability of `rogue_moons`), the survival curve helper of `rogue_moons`, and the unread moon loss fraction (stand-in row RK5) are removed with their tests and register lines. `rogue_moons.closest_approach_hill` stays, it is read by the starless host. Generation output is unchanged.
