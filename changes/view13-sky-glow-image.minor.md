**Headline:** A system's sky now has a glow image of the faint stars it cannot resolve, served with the sky for the page's band layer.

### Added

- `galaxy/sky_band.py` turns the glow band (the cube map of VIEW.27) into the equirectangular image `skyband.js` samples (VIEW.13): 256 by 128 texels in the sky's galactic axes (row 0 at the north pole, column 0 at longitude -180), three linear channels a texel, each the light in the page's unit `10^(-0.4 (SB - 22))`, tinted by a 4500 K blackbody (stand-in GLW1, 6.174), as little-endian half floats (196,608 bytes).
- `GET /system/<id>/sky/band?limit=&key=` (the page's route) and `GET /api/systems/<id>/sky/band?limit=` build that image for the 50 parsec cell of the system through the galaxy's density and dust. The sky's `meta.band` names it (`key`, `width`, `height`, `format`, `href`). The web layer keeps the bytes as `b<key>.bin` in the tile cache's generation, so the same cell and limit are worked out once, and `skycache.clear()` removes them with the skies. The key is a hash of the cell, the galaxy shape, the limit and what turns the band into the image. The API route is additive, so the API version stays at 12.
- The wire format is written down in `docs/design/interface/system-sky-wire-format.md` (section 3.2).

### Changed

- `skypack.SKY_VERSION` is 2, so every stored sky is rebuilt once with the `band` field in its `meta`.

generation-output: unchanged
