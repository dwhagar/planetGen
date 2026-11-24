**Headline:** The star mass function's centre is nudged from 0.079 to 0.092 solar masses (inside the paper's own uncertainty), so the primaries of 0.6 to 0.8 solar masses are in their census band and the star mass function census has no known anomaly left.

### Changed

- `tuning.IMF_BELL_CENTER_SOL` (stand-in MF1, tier T2, tweak level L1) is 0.092, not the published 0.079. Chabrier (2003), Eq. 17, gives the centre as 0.079 (+0.021, -0.016) and the width as 0.69 (+0.05, -0.01); the primaries of 0.6 to 0.8 solar masses were 7.1 percent of those of 0.08 to 0.8 in the 60,000-system census, under the 7.3 percent floor of the system function's band, and neither the companion counts (x0.5 to x3 moves the bin 6.8 to 7.3 percent) nor the companion mass ratio law (GEN.951) move it, so the single-star curve is what is nudged. More stars of 0.2 to 0.8 solar masses, fewer of 0.08 to 0.2, in every galaxy; the census shows all bins in band.
- `tuning.PREVALENCE_DRAW_SHARES["multi_star"]` follows the primary mass function (0.0046 to 0.0049).
- GEN.922 closes with this note: the star mass function census (TEST.454) has no known anomaly (the GEN.951 and GEN.996 tags are dropped), and the guards that held the generator to the published centre now hold it to the registered 0.092 and to the paper's range 0.063 to 0.100.
