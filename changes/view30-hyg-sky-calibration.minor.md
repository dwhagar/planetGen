**Headline:** The generated sky is now held against the real one: setup downloads the HYG star catalogue as a test reference, and a developer-run calibration counts the stars a generated sky shows at magnitude 6.5.

### Changed

- A new setup download, the HYG star database v4.1 (`datasets/hyg_stars.py`, listed in the README): 119,626 stars, 34 MB, only a test reads it, so bulk generation does not require it (`Dataset.required` False). It lists 8,920 stars at V 6.5 or brighter and 25,791 at 7.5.
- `tests/test_sky_hyg_calibration.py` (VIEW.30): plans a galaxy, scatters its bright stars, charts the sky of a cell and holds the dust-free count of stars at magnitude 6.5 to a window of the catalogue's. Two cases, each run on request with `PLANETGEN_SKY_CALIBRATION` (`small`, about 25 minutes; `milky-way`, about 2 hours 15 minutes; `all`), skipped otherwise and when the catalogue is absent. The Milky Way-sized case measured 5,733 stars, 64 percent of the catalogue's 8,920 (band 50 to 80 percent); the small galaxy the default seed draws, 24.5 percent (band 15 to 40).
- The measurement and what it shows (a third of the real stars from absolute magnitude 0 to 4 in the Milky Way-sized galaxy, many too many brighter than -6 without dust, a default dust law that cuts the small galaxy's sky to a quarter) are in reports/system-sky-hyg-calibration-2026-11-05.md.
