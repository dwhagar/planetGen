**Headline:** The nearby page tests build their little galaxy in minutes instead of hours on a slow machine.

### Changed

- TEST.475: `test_web_nearby.py` fills its two sectors without the deepest floor of the first fill and without the bright-star bands of a capped cell (the page tests do not read either), so the file runs in about six minutes where four of its eight tests timed out in setup at 900 s on a slow machine. Its hand-placed bright star carries the rotation period column GEN.846 added (the row had one value too few).
