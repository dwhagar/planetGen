**Headline:** The dust a planet takes in from its system's belts now follows a thin, drifting disk, so a Sun-like system gives the Earth about 1e7 kg a year, the measured order, where the old law gave 45 kg a year.

### Fixed
- **The belt and comet dust that reaches a planet was spread over a sphere (GEN.1543).** `physics/influx.py` spread a belt's dust as an isotropic stream through the sphere of the planet's orbit, so a Sun with a main-belt sized reservoir gave the Earth about 45 kg a year against the 1e7 kg a year the Earth takes. The dust is a disk that Poynting-Robertson drag carries inward: the planet now takes its focused cross-section times the grains' speed past it, over the annulus the stream fills (thickness 2 a sin i) times the drift speed `2 beta G M / (c a)` (`influx.capture_share`, `influx.radiation_ratio`, `influx.drift_speed_m_s`). The call takes the planet's mass and the star's luminosity too (`planet_influx_kg_per_yr` and the sheet row read them). New registered stand-ins 6.247 IF7 to IF10 (grain size, density, dust eccentricity and inclination); IF1 (the share of a belt's loss that is dust) is 0.03, tuned to the measured Earth influx at those values (tweak L2), and IF3 is the new geometry. The row is worked out when a sheet is made, so no stored value moves, and nothing the air reads changes (`late_accretion` and `late_supply` read only `influx.reach` and the reservoirs). The V8 dust supply of GEN.1542 is not built yet; it will read this influx when it is.

generation-affects: none

generation-output: unchanged
