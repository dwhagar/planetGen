**Headline:** A belt no longer keeps a world-sized body: every belt member heavier than the isolation mass of its stretch is a planet, and the belt's largest-body figure names the largest member it still holds.

### Changed
- **Belt members over the oligarch mass become planets (GEN.1561).** `belt_planets.promote` still takes the ten heaviest round members of a belt (`tuning.BELT_PLANETS_PER_BELT`), and now also every round member heavier than the isolation mass of the stretch the belt lies in (`physics.belt_origin.isolation_cap_earth`, median 0.08 Earth masses), whatever the count. The belt keeps its mass in the members below the cap; belt masses and member counts are as they were drawn (a Train 25 sample of 79 belts: unchanged).
- **The belt's largest body is the largest member it holds.** `AsteroidBelt.largest_body_km` was the drawn body's diameter, which after GEN.1510 is a planet; it is now the largest remaining member's diameter. In the 79-belt sample the median fell from 2,391 to 1,236 km, belts with a largest body over Pluto's diameter (2,377 km) from 40 to 6, and over Mars's (6,779 km) from 2 to 0; none is over Earth's radius any more (4 were). Real belts: Ceres 940 km, Pluto 2,377 km across ([R3538, R3539], R3539; planetesimals born small, Krivov and Wyatt 2020, R5253).
generation-output: changed
generation-affects: fill
