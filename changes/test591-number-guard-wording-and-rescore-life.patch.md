**Headline:** A re-score of a system now rates a body on every life class that lives on it, as generation does, so the stored life scores no longer drop after a re-score.

### Fixed
- `store.with_life_inputs` (the input of a habitability re-score, `refresh_habitability` and the fill command) now carries the body's `life_living` classes from `planet_life` and `moon_life`, the dominant class first. It left them out, so a body with two life classes was rated on its dominant class alone and a re-score gave lower `life_microbial` and `life_fauna` than generation stored (for example 2 and 0 against 3 and 1). This is the cause of the life half of the red `test_habitability_fill` parity tests (TEST.591).
- `test_number_guard` expects the message the code gives, "no unit, scale, or noun" (Oxford comma, rule 15).

### Notes
- TEST.591 stays open after this merges. The other half of the `test_habitability_fill` parity failure (`band1_day_hours`, because `values_from_body` lacks `day_share_mean`) changes what generation stores, so it is a lock change and rides GEN.1265 (#5053), which carries that line. The two parity tests pass only with both fixes.
