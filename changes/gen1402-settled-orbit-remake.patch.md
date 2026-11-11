**Headline:** A rocky planet or moon whose orbit wore down after its air was made now gets its air, climate, dose and spin made again on the settled orbit, even when its star's age barely moved.

### Changed

- GEN.1402: `refinish_rocky` made a rocky body again only when the star's age rose past `AGE_SETTLE_TOLERANCE`; a planet whose eccentricity `settle_orbit` lowered past `ECCENTRICITY_SETTLE_TOLERANCE` in a young-enough system kept an air, temperature, dose and spectrum made on the older, more eccentric orbit. The body now notes the eccentricity its air was made on (`_air_eccentricity`, kept, not stored) and `needs_remake` is true on a moved age or a moved orbit, for a rocky planet and for its moons (their host's eccentricity); a gas giant is made again on a moved orbit only (its temperature reads the light of the orbit; its age does not move its end stages). Reused: `settle_orbit`, `age_moved`, `remake_end_stages`, `refinish_rocky`; the radiation sweep test no longer needs its "air on the light now" stopgap and checks every body.

generation-output: changed
generation-affects: fill
