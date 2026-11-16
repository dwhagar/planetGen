**Headline:** A comet's period now lies in the period range of its own family.

### Fixed
- GEN.1612: the comet draw never read its family's period range, so 37 percent of elliptical comets (47 percent of the Jupiter-family, 37 percent of the Halley-type) had a period outside it. The eccentricity is now drawn from the part of the family's range that gives a period in the family's range at the comet's perihelion, and a perihelion that cannot give one is drawn again; a comet whose star's ceiling leaves no such orbit is named for the family whose range holds its period (`period_class_of`, which extends GEN.1608 from the 200-year line to every family). Over 1 500 comets none is now outside, and the families come out at their table weights (0.50, 0.31 and 0.19 of 0.5, 0.3 and 0.2, where the long-period family had been 0.14).

generation-output: changed
generation-affects: fill
