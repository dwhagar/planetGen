**Headline:** A and F stars now hold cold debris disks about as often as the surveys find them (17 and 20 percent of A1V and F5V stars, against 0 and 8 percent), because a wide belt's dust is the sum of its rings and a survey can find a fainter disk around a hotter star (GEN.1278).

### Changed

- `debris_dust.belt_fractional_luminosity` adds the dust of 32 rings of a belt, each with the belt's eccentricity and the share of the belt's mass the nebula's solids gave it, in place of one ring at the belt's middle radius (stand-in DU4).
- The fractional luminosity a survey can find falls as the cube of the star's temperature (`detectable_fractional_luminosity`, stand-in DU5); at the Sun's temperature it is 5e-5 (DU1 re-tuned from 3e-5). Shares of 100 stars each: A1V 17, F5V 20, G2V 17, K3V 18 percent (surveys 24, 24, 14, 13).
