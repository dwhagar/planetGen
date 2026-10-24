**Headline:** A world's stored UV index now comes from the light that reaches its ground, so ozone, air, haze and the star's own ultraviolet move it.

### Changed
- **GEN.445: the UV index is derived, not a proxy.** `radiation.generate_dose` stores `uv_surface_index` as the DNA-damage (Setlow) dose of the two-stream ground spectrum of `physics/body_radiation.py` (the body's air, ozone column, haze, dust and cloud, lit by its star's own spectrum and UV-C share) divided by the same dose at a reference Earth, so Earth is still 1 and the limits of the protection rules stay in "x Earth". A body with no air reads the star's unfiltered light; a body with no star has no index. The ratio-to-the-Sun index and its ozone exponent (`uv_ratio`, `uv_index`, `ozone_ratio` and their constants) are removed; a test fails if one returns.
- **Generator epoch.** The stored UV index and what reads it (the UV line of the survival gear, the habitability) change for every world with a star; no draw moves. The golden record is regenerated for it.
generation-output: changed
