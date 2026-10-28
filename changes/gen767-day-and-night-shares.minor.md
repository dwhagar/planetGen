**Headline:** Every planet and moon now knows how much of its time the star is up, how long its day is and how many hours of daylight it has, and a moon loses the time it spends in its planet's shadow.

### Changed
- **Day and night shares (GEN.767).** Worked out from the rotation, the orbit, the tilt and the star's size, with no draws: the share of the time the star is up at the equator and over the whole surface, the length of the solar day, the hours of daylight at the equator, the polar circle, and for a moon the share of its time in its planet's shadow (both day shares lose it). A world that turns once an orbit has no day, only a day side and a night side. Earth gets 50.15 percent at the equator and 12.04 hours of daylight; Io spends 5.3 percent of its orbit in Jupiter's shadow. The body's page shows them.
- **Stored (schema 170).** `planets` and `moons` gain `day_share_equator`, `day_share_mean`, `parent_shadow_share`, `solar_day_hours`, `daylight_hours` and `polar_circle_deg` (NULL for a row made before this release).
- Generator epoch bump (numbered at merge): new stored columns on every generated planet and moon.

generation-output: changed
TODO items resolved: GEN.767. Added: none.
