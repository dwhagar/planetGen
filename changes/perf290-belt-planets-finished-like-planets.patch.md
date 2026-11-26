**Headline:** A belt member made a planet is finished the way every planet is, so a young star's belts no longer take minutes to build.

### Changed
- A belt member promoted to a planet is built bare and finished with the late accretion the system gives it (`Planet.finish_climate`, `belt_planets.finish_planet`), the way every other planet is, so its air is made once instead of with the law's draw and again after the belts are final. A system of a 4 Myr star with seven belts (seventy promoted bodies) builds in 113 s instead of 174 s, one with four belts (forty) in 54 s instead of 68 s.
- Every promoted body is drawn again from its own stream in the order the planets are (motion, air, spin, moons): its spin, tilt and day change, and its moons may; its class, mass, orbit and makeup do not. No other body of the system moves except in the last digit of a tide sum (1e-16 relative).
- Census (seeds 2, 3 and 4, 120 promoted bodies): the class of none of them changes; the spin, tilt and day of every one does.

generation-output: changed
generation-affects: fill
