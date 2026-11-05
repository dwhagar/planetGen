**Headline:** The seven latitude bands of a rocky world now take the size of its star in the sky at generation, as a re-score always did, so a re-score changes no stored number.

### Changed
- `habitability_world.values_from_body` did not pass the body's `day_share_mean` (the mean share of the time the star is up, which holds the star's angular radius) to `compute`, so generation worked the bands (`bandN_day_hours`, `bandN_night_hours` and the longest day and night) for a point star, while a re-score read the stored column and worked them for the star's disc. A 'wide' system's first band had a day of 13.60 hours at generation and 14.21 hours on a re-score. `values_from_body` now carries `day_share_mean`, so generation and a re-score agree on every `SCORE_FIELDS` column.
- Only the bands that depend on the star's angular size move, and only for a world that turns faster than it orbits; the draws are untouched.

### Added
- `tests/test_habitability_fill.py::test_a_re_score_reproduces_every_stored_number_of_systems_of_any_kind`: 14 saved systems, every score and band column compared after a re-score (it failed on the first band's day length before the fix). The two older fill tests that were red on main pass again.

generation-output: changed
TODO items resolved: GEN.1283 (the epoch and golden ride the train; GEN.1265 in PR #5053 carries the same one-line change in `values_from_body`)
