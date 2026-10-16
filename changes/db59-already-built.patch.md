**Headline:** The TODO list retires DB.59, whose nebula and remnant age columns turned out to exist already.

### Changed
- **DB.59 is done: already built.** Core 2 found that `nebulae` has had `age_myr`, `q_total` and `hosts` since v87 and the remnants table `age_years` since v38, saved by the store and covered by the HII and column round-trip tests. GEN.357 can start now, and GEN.358 and GEN.95 no longer list DB.59. Migration revision 0102 is free: DB.54 takes 0102 after DB.55's 0101, and Core 1 fixes the order if DB.54 merges first.

TODO items resolved: DB.59 / Added: none
