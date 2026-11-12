**Headline:** The pull of a distant companion star or a heavy outer planet on a planet's orbit can now be worked out in one place; nothing reads it yet.

### Changed

- GEN.1447 (slice 0, planet and sibling-planet cases): new `physics/kozai_planets.py`, pure. `screen_planets` runs the Kozai-Lidov screen of `kozai_bodies` for each planet of a star against every perturber past three times its axis: a companion star or dark member with the plane of its orbit, and an outer planet at least as heavy as the planet. The other precession a cycle has to beat is `secular_chain.precession_rates` (siblings and relativity). Each planet gets the outcome of its strongest cycling pull (the eccentricity it reaches, where it stands at the star's age, whether the periapsis falls inside the star), the perturbers that cycle it and the ledger rate of every perturber. Stand-ins KZP1 and KZP2 (tier T6, tweak L1) are in section 6.284 of the gap-filling plan.
- Nothing in the generator reads it yet (wiring register): applying the outcome to the planets is a lock item that waits for the stored companion planes of GEN.1500.

generation-output: unchanged
