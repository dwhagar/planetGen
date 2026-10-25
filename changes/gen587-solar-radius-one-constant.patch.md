**Headline:** The tuning surveys read the Sun's radius from the one shared constant instead of four typed copies.

### Changed
- **Four copies of the solar radius are gone (GEN.587, a further file group).** The Kepler recovery, radius valley, mass-radius
  sample and JWST selection modules now take `SOLAR_RADIUS_M` from `physics/constants.py`. The values are bit-identical, so no
  generated output changes and no epoch is taken. The map script keeps its own literal on purpose (it runs without the package).
