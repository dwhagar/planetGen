**Headline:** The test table of the Milky Way's satellites no longer gives Leo I a rising star-formation history and so gas it does not have.

### Changed
- **TEST.508.** The Leo I row of `tests/reference_satellites.py` had a mean age of 5 Gy, which the satellite history reads as a rising rate: it formed stars now and held 1.4e6 solar masses of gas, where Leo I has under 3e4 of HI (Grcevich and Putman 2009). The row is 9 Gy now, so it forms none; the age is the table's own choice (its docstring says so). `test_satellite_history` lists Leo I and Leo II among the satellites that form no stars. Test data only; nothing a generator draws or stores changes.
generation-output: unchanged
