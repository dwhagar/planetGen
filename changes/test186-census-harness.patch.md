**Headline:** The test suite has one census harness that fills seeded sectors per region class and checks the counts of each kind against their expectation.

### Added
- **The census audit harness (`src/tests/audit/census.py`).** It fills a fixed set of seeded sectors for each region class, counts what each fill made by stage and kind, tests each count against its expectation within four Poisson sigma, applies the rule that a kind expected 20 or more times appears at least once, flags a kind found where its gate expects none, and writes one report. The representation, fill-stage, multiplicity, catalogue-type and whole-sphere audits add their checks to it and keep no counters of their own.

TODO items resolved: TEST.186. Added: none.
