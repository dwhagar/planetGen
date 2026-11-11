**Headline:** A rogue planet, and the moon of one, no longer has a day, daylight or latitude bands from a star it does not have: it is dark all year.

### Changed
- **No sun, no day (GEN.1496).** A body under a star-less host (a rogue planet, or the moon of one) has the empty set of suns in `daylight.compute`: `day_share_equator`, `day_share_mean` and `daylight_hours` are 0, the polar circle, the solar day, the night hours and the days a year are empty, and the state word is the new `no sun` (`day_night.NO_SUN`; a giant still says `no solid surface`). `latitude_bands.compute` gives no bands for such a body (`habitability_world.values_from_body` now says whether the star is star-less), so the rogue planets of a galaxy store no `band1_*` to `band7_*` and the page shows no band table for them. Before, they stored a day share of 0.5 and half a rotation of daylight, and a rocky one stored bands in the state `polar day and night`. Surface temperatures are unchanged. New epoch, `--affects fill`.

generation-output: changed
generation-affects: fill
