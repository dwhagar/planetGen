**Headline:** Moons and the planets of a wide pair are held to every neighbour they have, not to the last or the outermost one alone.

### Changed

- A moon now clears every moon inside it by that moon's own five Hill radii (`moon_clearance_au`), where only the moon just inside it was checked: a massive moon two places in keeps its room even when a light moon sits between. This holds in the generated lunar system, in the stabiliser of an edited system, in the validator, and in the history of the moons' tidal drift.
- The planets of the two stars of a wide pair are held to each other for every pair of bodies (`cross_star_violation`), by the criterion that was already there (Gladman's mutual Hill factor at the stars' closest approach), where only the outermost body of each star was checked. A heavy planet further in, that needs more room than a light outer one, can now be trimmed.
- The planet-planet spacing already held every inner planet and belt (the clearance ledger), so nothing changes there. The circular Hill radius stays for `min_orbit_distance`, registered as SP1.
- On 150 systems (69 wide pairs) and 202 further wide-pair builds no system changes and the pinned golden battery is unchanged: the rules now hold where no generated system broke them, so this is a guard and not a change of output.

generation-affects: none
