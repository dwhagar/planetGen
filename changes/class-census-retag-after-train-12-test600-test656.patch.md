**Headline:** The class census's known-fault tags match the run after train 12: 69 tags that no longer trip are gone, and the nine new flags are tagged with the reason each one is not a new fault.

### Changed

- `src/tests/data/prevalence_expectations.json` (TEST.600, TEST.656): the 69 stale tags the census printed are removed (and the GEN.1226 note, which no range uses now). The new flags: the A0V moon Class B shares (5 of 24 at the ecosphere, 2 of 11 hot) are impact moons of 2 to 11 Earth mass planets and follow the planet Class B share of the same cells, so they join GEN.728. The G2V cold moon Class D share (42 of 47, 0.894 against 0.806) joins GEN.1160. The K5V ecosphere unlimited share (1 of 8) and the M2V cold microbial and none shares (7 of 78) are TEST.656: since GEN.1082 to GEN.1084 put captured moons in the irregular store, which the census does not read, a rocky planet's moon list holds only impact and binary moons (24 moons in 300 planets at A0V ecosphere against about 660 in the reference), so those cells hold 1 to 78 bodies against ranges measured on 500 to 700. The M2V hot moon cell is empty by physics: of its 300 planets (0.03 to 0.12 AU, days of 4 to 29 h) 15 were born with a moon, the tide brought 12 down onto the planet (`fallen_moons`, GEN.959) and the other 3 were lost, so none remains. The two draws over 45 s are PERF.246.
- `src/tests/audit/prevalence.py` (TEST.656): a cell with no bodies takes the `anomaly` of its spec, so an empty cell that is right by physics can carry its reason.
- `src/tests/test_prevalence_census.py`: a test of that, and a test that every cell-level anomaly names an item that has a reason.

generation-output: unchanged
