**Headline:** The generator's close-binary shares are now checked against the published ranges by primary mass and [Fe/H], and a bin outside its range fails the test run.

### Added
- **The close-binary check (TEST.196).** `tests/audit/close_binary.py` draws primaries of six masses at seven [Fe/H] with the generator's own companion draw and compares the share with a companion inside 10 AU with Moe, Kratter and Badenes 2019 (solar-type ratios to [Fe/H] 0), Moe and Di Stefano 2017 (mid B at least 16%, O between 23% and 78%) and the no-trend finding for O and early B stars, plus the outer share's flatness in [Fe/H]. Absolute levels are printed, only ratios fail. `tests/test_close_binary_check.py` runs it.
