**Headline:** Every temperature on a world's page now shows the day and night range of its own place.

### Changed
- **GEN.848: a day and night swing is stored for every place that has a temperature.** Until now only the equator had one, so the polar, northern, southern and mean temperatures read "no day and night range stored for this place". Seven new columns on planets and moons (`temp_polar_swing_k`, `temp_north_swing_k`, `temp_south_swing_k`, `temp_mean_swing_k`, `temp_day_side_swing_k`, `temp_terminator_swing_k`, `temp_night_side_swing_k`; migration, schema v178) hold them. A band's swing is the equator's cut by the cosine of its latitude (the daily variation of the sunlight, which the ground and air carry linearly), the mean swing is the area-weighted mean of the bands, and a tidally locked world, which turns once an orbit and has no day, stores 0 for its three zones and says so on its page. The regions object of the API gains the same swings by name (`polar_swing`, `mean_swing`, ...). The climate sentence already states the equator's low and high and agrees with the page; it is unchanged.
- **Migration 0182 adds the columns empty;** a galaxy made earlier gets them when it is regenerated.

generation-output: changed
