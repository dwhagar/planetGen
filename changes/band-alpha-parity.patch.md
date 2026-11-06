**Headline:** The latitude bands of a new rocky world are worked out with the size of its star in the sky, as the re-score of a stored world already did.

### Changed

- `habitability_world.values_from_body` did not hand `day_share_mean` (the mean share of the time the star is up, from which the star's angular radius follows) to the latitude bands, so generation computed every band as if the star were a point, while `rescore-habitability` and the API read the stored share and used the real angular radius. The two disagreed on the day and night hours and the temperatures of the bands of every rocky planet and moon with a rotation, and `test_habitability_fill` failed in two cases on main (a wide-pair planet's `band1_day_hours` 13.09 h stored against 12.65 h re-scored in one seed). Generation now passes the share, so the stored bands and the re-scored ones agree.

generation-output: changed
