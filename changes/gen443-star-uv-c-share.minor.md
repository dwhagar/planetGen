**Headline:** Every star now carries its own ultraviolet-C share, and a planet or moon's ground light and UV index can be read from what it stores.

### Added
- **GEN.443: the ground light of a stored body.** `physics/body_radiation.py` turns a planet's or moon's record (its air, its star and distance, its ozone column and its haze, dust and cloud) into the band spectrum at its ground with the two-stream model, and gives its UV index and DNA dose. Migration 0145 adds `stars.uv_c_share` (schema v145), the share of the star's blackbody UV-C it gives off, drawn once from the star's own stream (log-uniform 0.1 to 1); a row saved before it keeps NULL and reads as the blackbody until the galaxy is regenerated.

### Changed
- **Generator epoch.** Each ordinary star carries one more stored value, drawn on its own stream so no other draw moves; the golden record is regenerated for it.
generation-output: changed
