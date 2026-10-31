**Headline:** The Milky Way satellite table the plausibility tests use is now checked against McConnachie (2012).

### Changed
- **GEN.366.** The 11 rows of the reference table (`tests/reference_satellites.py`, test data only; the generator reads none of it) are checked against the published tables 2, 3, 4 and 5: the distances of the Large and Small Magellanic Clouds, Sagittarius and Fornax (51, 64, 26 and 147 kpc), the stellar masses of the Large Magellanic Cloud (1.5e9, was 2.5e9), the Small Magellanic Cloud (4.6e8) and Sagittarius (2.1e7, was 1.0e8), the half-light radius of Ursa Minor (181 pc, was 405) and the mean [Fe/H] of every row are corrected, and a new test pins every value to the paper. Generated worlds do not change.
generation-output: unchanged
