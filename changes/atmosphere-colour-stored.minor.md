**Headline:** Every planet and moon with air now stores the colour of that air, ready for the maps to paint with.

### Changed
- **The atmosphere colour is data (DB.96, schema v162).** `planets` and `moons` keep `atmosphere_colour`, the colour of the body's air as one sRGB integer, worked out from its air, its star, its ozone and its stored haze, dust and cloud (`physics/body_colour.py`; a gas or ice giant takes the colour of its air at the reference level). A body with no air stores nothing. The cloud cover it used is the `cloud_cover` the body already kept. The ground under the sky is still the labelled stand-in reflectance. Old galaxies keep an empty column until they are regenerated.

generation-output: changed
