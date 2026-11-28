**Headline:** The share of Sun-like stars that hold a cold debris disk is back inside the survey's range.

### Fixed
- **GEN.1884: the stars with a cold debris disk had grown to 40 percent of G dwarfs against the surveys' 14.** Since GEN.1797 planned the planets' spacing, the outer belts start nearer and are heavier, so more of them passed the cut `debris_dust.DETECTABLE_FRACTIONAL_LUMINOSITY` (stand-in DU1). The cut is re-tuned in place from 5e-5 to 7.8e-5 on the row's own sweep of 100 stars each of A1V, F5V, G2V and K3V: A 0, F 7, G 17, K 26 percent (before: 4, 26, 40, 54). G sits in its band; A, F and K do not, recorded as a known gap in the row, since no one cut fits all four (the cut each type needs climbs from hot to cool stars, GEN.1892). The new cut is about eight times the surveys' own limit, so it is a cut on the mass of the generator's belts, not an instrument limit.

generation-affects: fill
generation-output: changed
