**Headline:** The other stars of a system are suns of their own in a planet's sky, so a planet's seasons and the temperature of every place follow where each star stands.

### Changed
- **Seasons with a second sun (GEN.1458).** `physics/sky_suns.py` gives each planet of a multi-star system its other suns: each other star's share of the mean light and its direction in the planet's equinox frame, worked from the stored spin axis, the orbit's normal and the star's place in the orbit tree. The year's mean light at a latitude, the latitude bands and the cell-by-cell day of the temperature series add each sun at its own declination and right ascension, so a planet with a bright companion in its summer sky has a milder winter and a pole under a companion stays warm. A planet stores the suns as `sky_suns` (new column on `planets`, `moons` and `rogue_planets`, migration 0236; moons and rogue planets store none yet), so a re-scored row gives the same temperatures. A lone star gives the same numbers as before. The direction is the host-centered one and fixed through the year; moons and ring shadows read the host's sun alone (registry rows SS1 to SS3, section 6.246). New epoch, `--affects fill`.

generation-output: changed

generation-affects: fill
