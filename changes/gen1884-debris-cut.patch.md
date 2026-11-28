**Headline:** The share of stars with a cold debris disk now follows a survey's flux limit at the distance of its sample, so F and G stars land in the survey's range.

### Fixed
- **GEN.1884: the stars with a cold debris disk had grown to 40 percent of G dwarfs against the surveys' 14.** Since GEN.1797 planned the planets' spacing, the outer belts start nearer and are heavier, so more of them passed the cut.
- **GEN.1892: a cut that rose 4 times from A to K stars is now a physical law.** A survey finds dust by its flux, and its sample is volume limited by spectral type (the farthest F, G and K star of the Herschel DEBRIS sample is 24, 21 and 16 pc away, [R4081]), so the dust it can find around a brighter star is a smaller share of its light. `debris_dust` replaces the cube of the star's temperature (row DU5) by the distance of the sample, 21 pc times the star's luminosity to the 0.164 (fitted to the three published distances), and the cut DU1 becomes one flux (dust luminosity over distance squared, 1.9e-7 in solar luminosities per parsec squared) tuned on a sweep of 100 stars each of A1V, F5V, G2V and K3V: A 13, F 28, G 17, K 7 percent (before: 4, 26, 40, 54). F and G are inside their surveys' bands; K is 2 points and A 7 under (A: GEN.1895).

generation-affects: fill
generation-output: changed
