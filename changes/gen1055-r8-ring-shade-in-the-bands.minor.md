**Headline:** A tilted rocky planet with a broad ring has cooler bands behind the ring and warmer bands at its poles, because the ring now takes sunlight off the latitudes it shades.

### Changed
- **The ring's shadow in the seven latitude bands (GEN.1055, R8).** `latitude_bands.ring_shade` reads the ring (`ring_inner_km`, `ring_outer_km`, `ring_optical_depth` and the radius) and, season by season, the share of the day's sunlight that reaches each latitude through it (`ring_darkening.shadow_transmission`). A band keeps the insolation-weighted mean of that share, and its temperature follows. A world with no tilt, no ring, or a spin-orbit lock is unchanged. The world's own mean temperature is not lowered (it is solved before the ring is drawn), so the shade moves heat between bands: the equator band of a Saturn-like ring on a 27 degree tilt is about 8 K cooler, the poles about 10 K warmer. A moon's ring does not shade its bands, because a moon's bands are made before its ring is drawn.
- The climate values of a body (`habitability_world.values_from_body`) carry the three ring columns.

No migration (the band columns already exist).

generation-output: changed

Epoch affects: fill. The band temperatures of ringed rocky planets change; no new draw.

TODO items resolved: GEN.1055 part R8 (GEN.1055 is complete with this and R6).
