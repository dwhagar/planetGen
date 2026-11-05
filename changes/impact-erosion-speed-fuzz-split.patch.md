**Headline:** Generating a system around a luminous star is faster, and two small slips that the old time limits were hiding are fixed.

### Changed
- **The impact erosion sum is cheaper (PERF.229).** The mass-weighted erosion of a planet's air by the impactors that hit it recomputed the same impactor-density factors for every grid step. They are now worked out once per call, and the sum is taken in the same order, so the result is identical to the last digit (a test pins it against the sum of the single efficiencies).
- **The star-type fuzz test is split by spectral letter (PERF.229).** The heavy luminosity classes (Ia, Iab, Ib, II, III) each ran the whole letter range in one test, so the slowest cells hit the time limit and hid what failed behind them. Each letter now runs on its own, with its reason recorded in the test.
- **A belt's "stirred" and "moved through" flags are plain true or false (PERF.229).** They could come out as numpy booleans, which the stored system could not serialize.
- **A belt cut back keeps its middle inside its limits (PERF.229).** When a belt was cut back to the orbits a star holds, its distance could stay past the new outer limit (a luminous giant's belt had limits 497 to 575 AU and a distance of 643 AU). The distance now stops at the new limit.

generation-output: changed

TODO items resolved: PERF.229.
