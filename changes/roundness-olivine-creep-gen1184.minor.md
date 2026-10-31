**Headline:** The roundness engine now knows that rock that took its shape hot is much weaker than cold rock, using the olivine creep laws, and the hot rubble-pile strength has a physical basis.

### Changed
- **GEN.1184.** `physics/roundness.py` gains an `olivine-creep-law` temperature rule for rock: with a shape-set temperature, rock follows dry olivine dislocation creep between the Hirth and Kohlstedt 2003 and Korenaga and Karato 2008 columns (by the body's strength draw), never above the 16.5 MPa cold rock strength, so rock stays 16.5 MPa up to about 1100 to 1150 K and falls to 2 to 5 MPa at 1300 K. The temperature rules are now per material in `physics/data/roundness_strengths.json` (`ice`, `rock`). The hot rubble strength (RN9, 2.5 MPa) is re-tiered to L1 as olivine creep near 1330 K. Nothing reads the engine yet, so generated worlds do not change.
generation-output: unchanged
