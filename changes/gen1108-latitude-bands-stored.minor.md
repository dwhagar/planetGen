**Headline:** A planet or moon now keeps its climate in seven bands of latitude, and its page shows the seven side by side.

### Changed
- **Seven bands stored (GEN.1108 part B, schema v196).** `planets` and `moons` hold, for each of the seven bands of `physics/latitude_bands.py` (the north pole to the south pole, or the point under the star to the point away from it for a world that turns once an orbit), the annual mean temperature, the day and night swing, the hours of day and of night, the longest polar day and night of the year and a state word (`band1_mean_k` to `band7_state`, 49 columns), and `temp_mean_swing_k` for the whole surface. Each is worked out from the stored temperature, spin, orbit and air with no new draw. The columns `temp_polar_k`, `temp_north_k`, `temp_equator_k`, `temp_south_k`, `temp_swing_k`, `temp_day_side_k`, `temp_terminator_k`, `temp_night_side_k` and the six swings by place are gone, with no compatibility shim; the migration adds the new columns empty and drops the old ones.
- **The sheet page shows the bands.** A planet's or moon's page has a "Climate by latitude" table with the seven bands as columns and the mean, the swing, the hours of day and night, the longest day and night and the state as rows. It replaces the old rows by place and the note "no day and night range stored for this place". The mean temperature row keeps its range over a day.
- **Text and hours read the bands.** The Temperature row's regional sentence names the equator, the poles and (for a locked world) the day side, the terminator and the night side from the bands, and the gear check runs over all seven. The hours a person can be outside read the equator band's swing, or a locked world's terminator ring.
- **API version 12.** `regions` of a planet's or moon's `equipment` is `{bands, mean_swing_k}` with the seven bands, and the sheet of a planet or moon has `bands`.

generation-output: changed
api-version: changed
