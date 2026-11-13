**Headline:** A planet whose thermal and gravitational tides balance away from synchronism can no longer spin faster than 85 percent of its breakup rate.

### Fixed
- **The asynchronous balance is capped (GEN.1576).** `planets.generate_spin` took the rate of a balance of the two tides (`spin_torques` kind "asynchronous") as it came, so a puffy planet (3.5 Earth masses at 9 Earth radii, 0.15 AU from a G2 star) turned once in 6.2 hours where its surface would leave the ground at 24.4. The rate is now held at `spin.BREAKUP_FRACTION` of the breakup rate of the planet's mass and radius, the limit every other spin already obeys. Only a planet whose balance was past that limit changes. The red `test_orbital_motion` rotation-period test came from this (its sample is unseeded, so it failed now and then).

generation-output: changed
generation-affects: fill
