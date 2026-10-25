**Headline:** The map browser tests have a dense fixture galaxy to fly into, ready for the zoom, label, pick and icon work.

### Added
- **A fixture galaxy with a dense region (TEST.368).** `tests/map_site_support.py` gains `dense_site` (a local map server with about three thousand seeded stars packed round the middle of the first fixture sector, a few named, six rogue planets, a neutron star, a black hole and five habitable-world systems, no scatter) and `dense_region_query()`; `tests/test_web_browser_fixture_maps.py` gains `open_dense_region()` and a test that zooms into the region and finds it on screen.
