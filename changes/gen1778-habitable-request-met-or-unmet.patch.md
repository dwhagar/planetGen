**Headline:** A request for a habitable world now says "unmet" when its one temperate ocean world has no life, instead of "aim" over a count of 0.

### Fixed

- GEN.1778: `test_full_system_generates_with_belt_and_habitable_world_across_every_star_type[K5IV]` failed under xdist with "HABITABLE_WORLD=True produced no habitable world and did not say so". The cause is not shared state: about 1 K5IV system in 75 is a 50 Myr star (younger than the first step of every life class) whose one world meets the request's conditions (301.6 K, 56 kPa, a surface ocean over 71 percent), so the request read "aim" while no world had life and `habitable_class_count` was 0. The request is now answered "aim" only when a body meets the conditions and has life on it (`habitable_request.count_living`), else "unmet", so a finished system always holds a habitable world or says it does not. No draw changes: the odds the request raises, the bodies and the life pass are as they were; only the report answer of that case moves.
- `test_full_matrix` now binds a run seed of its own for each star type, so a failure replays alone or under xdist. The assertion is as it was, and no seed was picked for passing.

generation-output: unchanged
