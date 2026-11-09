**Headline:** A planet or moon with more than one bright sun in its sky has a longer day: the share of the time some sun is up is the union of the arcs each one lights.

### Changed
- **Day and night with every sun (GEN.1461).** `daylight.compute` now takes every star that puts at least civil-twilight light (0.036 W/m2) on a body as a sun of the day, the host included, and merges the arcs of the turning ground each one lights (`daylight.lit_arcs_share`, `sky_day_shares`). The suns stand where the body's place on its circle of mean flux and the star's offset from the host put them (from `sky_suns.placed`), so a close pair reads as nearly one sun and a wide pair as two: `day_share_equator`, `day_share_mean`, `daylight_hours` and the night and day/night state follow, on planets and on moons (a moon keeps its parent's shadow on top). A lone star, a host with only faint companions and a body without a stored spin axis read as before. The polar circle and the solar day stay the host's (registry rows DS1 to DS4, section 6.243). New epoch, `--affects fill`.

generation-output: changed
