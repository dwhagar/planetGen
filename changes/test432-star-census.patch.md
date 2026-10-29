**Headline:** The prevalence suite now checks the stars the generator makes against the local star census and the multiplicity surveys.

### Added
- **The star census (TEST.432).** `src/tests/audit/prevalence_stars.py` draws 40,000 stars at the default galaxy's reference solar position and 2,500 planetless systems around each of M2V, K5V, G2V, F5V and A0V. It checks the share of M, K, G, F, A, B and O stars, white dwarfs and evolved stars against the RECONS 10 parsec census and Mamajek's number densities, and the share of systems with a companion star (G2V also with three or more) against Duchene and Kraus 2013, Raghavan 2010, De Rosa 2014 and Moe and Di Stefano 2017, in `src/tests/data/prevalence_stars_expectations.json`. The ranges are the researched norms, widest where two sources differ; each run adds four binomial sigma of its own sampling. The first draw found three faults (GEN.897 too many K dwarfs, GEN.898 too few A and B stars, GEN.899 too few companions of M dwarfs). A range the generator breaks carries the item that will fix it and prints in `star-census-latest.md` instead of failing.

TODO items resolved: TEST.432. Added: none.
