**Headline:** The property sheet says why a close pair's stars have no critical orbit, and the comet period rule is confirmed.

### Changed
- **GEN.789: the critical orbit radius (`wide_binary_a_crit_km`) now applies only where the generator sets it.** It is set for the stars of a wide pair and for the members of a multiple system; a single star has none, and neither do the two stars of a close pair, whose planets orbit the pair out to its own limit. The sheet rule for it gained that last case with its reason, the fill audit reads each star's system, and the open-gap entry is gone (in 120 systems, no star is empty where the column applies).
- **GEN.790: the comet period columns (`mean_anomaly_deg`, `orbital_period_years`, `period_class`) were already ruled by orbit type.** A comet on a bound orbit owes all three and a parabolic or hyperbolic one has none; in 120 systems no bound comet was empty. A test now states the rule.
- Nothing the generator makes changes.
