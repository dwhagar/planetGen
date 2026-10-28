**Headline:** The planet prevalence checks draw their sample once and share it, and the draw itself is faster.

### Changed
- **The checks of `test_prevalence.py` share one drawn sample, and drawing a system is faster (TEST.392).** Each check drew its own 2,000 systems; they now draw each seed's sample once and share it, with the same sample size and the same limits. Under the draw, the band insolation is worked out for every latitude at once (it differs from the loop's by at most 9e-16, the last bit, far inside the nine digits the golden fingerprints keep; the golden battery digests still match), and a suit's shield answers repeated limits from a table.

TODO items resolved: TEST.392. Added: none.
