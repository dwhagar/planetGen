**Headline:** Close-in planets no longer pile up at 12 Earth masses; they thin out toward 20, and a system holds fewer excess warm sub-Neptunes.

### Changed

- GEN.1861: the tuned 12 Earth-mass limit on a solid core inside the snow line (`PLANET_INNER_CORE_MAX_EARTH`, tweak level 3) is gone. A core that finds no gas now stops near an ice giant's mass, 20 Earth masses, inside the snow line as well as beyond it (`GIANT_NO_GAS_CORE_MAX_EARTH`, row PM8, the existing sourced limit). The old limit pressed 0.083 planets a star into 11.5 to 12 Earth masses and 0.16 a star into 10 to 12; those now spread from 10 to 20. Generated values change (the mass, radius and density of the heaviest inner planets, and the planet count of the systems they crowd); no epoch bump.

generation-affects: fill
