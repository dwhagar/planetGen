**Headline:** The generated sky is now held against the real one: setup downloads the HYG star catalogue as a test reference, and a developer-run calibration counts the stars a generated sky shows at magnitude 6.5.

### Changed

- A new setup download, the HYG star database v4.1 (`datasets/hyg_stars.py`, listed in the README): 119,626 stars, 34 MB, only a test reads it, so bulk generation does not require it (`Dataset.required` False). It lists 8,920 stars at V 6.5 or brighter and 25,791 at 7.5.
- `tests/test_sky_hyg_calibration.py` (VIEW.30): plans a galaxy, scatters its bright stars, charts the sky of a place where the density is 1 and holds the dust-free count between 15 and 40 percent of the catalogue's (measured 24.5 percent: 2,184 stars). Run with `PLANETGEN_SKY_CALIBRATION=1` (about 40 minutes); skipped otherwise and when the catalogue is absent.
- The measurement and what it shows (a fifth of the real stars from absolute magnitude -3 to 4, 24 times too many brighter than -6 without dust, a default dust law that leaves nothing beyond 217 parsecs) are in reports/system-sky-hyg-calibration-2026-11-05.md.
