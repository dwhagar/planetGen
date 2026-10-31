**Headline:** The strengths used to decide whether a body is round now follow the newest strength research: ice depends on temperature, ice and rock mix by a measured law, and metal has read data behind it.

### Changed
- **GEN.1158.** The stand-in strengths in `physics/data/roundness_strengths.json` are updated from `reports/planetary-material-strength-2026-10-31.md`: ice is 1 MPa by default and, when a shape-set temperature is given, follows the creep law of ice (about 34 MPa at 100 K, 3.6 at 140 K, 0.9 at 160 K, 0.09 at 200 K); a mix of ice and rock has the strength of its ice times exp(2 x the solids' volume fraction), replacing the threshold rule; metal is 100 MPa (range 50 to 250); porosity stays off. Nothing reads the engine yet, so generated worlds do not change.
generation-output: unchanged
