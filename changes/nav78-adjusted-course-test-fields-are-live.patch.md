**Headline:** The tests of the adjusted course place a field that is still a field, so the course bends around it again.

### Fixed
- `tests/test_nav_adjusted.py` made its asteroid field 20 light-years wide. A field wider than a light year is retired (GEN.1107) and is no obstacle, so the course rightly did not avoid it and three tests failed. The field is now 0.9 of the retirement radius, and every assertion is kept.
