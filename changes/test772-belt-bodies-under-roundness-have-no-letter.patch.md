**Headline:** The on-demand belt body tests look for a body that is round enough to have a class letter, as the generator does.

### Changed
- TEST.772: `test_belt_ondemand` took the body of one rank and wanted a class letter; a body whose planet comes out under a roundness ratio of 1 is an Asteroid, an object type with no letter (GEN.1160), and the belt's status (median strength) can call a body just over the round mass round while its planet draws its own strength. The two tests now take the first body with a letter among the next ranks and hold every body without one to a roundness ratio under 1; the on-demand planet is still the planet `belt_planets.planet_of` makes.
