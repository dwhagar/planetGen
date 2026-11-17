**Headline:** The members of a belt whose star has a tilted companion star now carry the inclination and node that the companion's secular pull and the planets give them; a lone star's belt is unchanged.

### Changed
- **A companion star tilts its star's belt members (GEN.1492, second part, slice 1).** `AsteroidBelt._force_by_companions` now turns each member's drawn orbit normal about the forced axis of `physics/belt_forcing.py` (`forced_normal`: the planets' plane weighted by the star's nodal rate on the member, plus each companion's orbit plane weighted by its secular rate), to a phase drawn from one stream of the belt, and keeps the angle to that axis (`tilted_orbit`). A bare belt follows the companion's plane, a belt beside a giant stays near the planets' plane. A lone star, a companion with no known plane, a coplanar companion, a member beyond `ALPHA_LIMIT` of the companion and a forced axis past 0.5 rad from the planets' plane keep their drawn inclination and node exactly. Reused, not rebuilt: slice 0's `belt_forcing`, `secular_chain.matrices`, `radiation_field.companions_of`, `spin.orbit_normal`. The three wiring-register lines are gone as the functions are now read. Census of 300 generated pairs (490 belts of stars with a companion): the tilt moves 127,724 of 143,510 members in 456 belts, and the mean inclination of a tilted belt goes from 22.3 to 24 degrees (`tests/audit/belt_forced_members_census.py`). Registry rows FE3 and FE5 of section 6.256 are unchanged. No migration.

generation-affects: fill

generation-output: changed
