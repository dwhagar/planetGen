**Headline:** A planet or moon loaded from the database now starts from the numbers it was born with when it is moved or edited, so its moons and its draws come out as they did when it was made.

### Changed

- GEN.1753: `body_history` keeps where each body formed (`formation_distance_au`, the key of its end-stage draw streams and the line its water's origin is read against), a moon's birth orbit (`birth_distance_au`) and its place in the planet's born moons (`born_index`); a new table `lost_moons` keeps the moons a planet was born with that its history took away (mass, birth orbit, density, direction). Migration 0266 (schema v266) adds them; nothing is backfilled.
- `store.load_star_system` gives each restored body its formation distance and each restored planet its born moons (`moon_history.restore_born`), the lost ones as `LostMoon` stand-ins that take part in the history and are never a body of the system. A moved planet or a star change then runs its moon history from the start it was born with and not from the orbits the moons have drifted to, and the end-stage draws of a moved body no longer change because its distance did.
- An admin edit writes the birth numbers again from the bodies the system has (`store.rewrite_body_history`): a body the edit made again has its own, no longer its predecessor's row.
- A planet saved before v266 behaves as before (its moon history starts from the moons it has).

generation-output: unchanged
generation-affects: none
