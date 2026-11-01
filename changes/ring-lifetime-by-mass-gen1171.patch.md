**Headline:** A ring lasts as long as its mass and its planet allow: light rings are worn away by meteoroids in thousands to millions of years, so far fewer planets keep one.

### Changed
- **GEN.1171.** `ring_origin.spread` no longer drops every ring by 1/e per half billion years. A ring now loses mass by viscous spreading (quicker the denser the ring, scaled from Saturn's rings, 2 Gyr e-fold at their mass) and by meteoroid erosion (Saturn's ring flux scaled by distance from the star and the planet's gravitational focusing, over the ring annulus from 1.24 radii to the Roche limit), `dM/dt = -k M^2 - b`, solved in closed form. A 1e13 kg ring of a fallen moon lasts about 1e5 years, a Saturn-mass ring most of the age of the Solar System. In a 60-system sweep rocky planets with a ring fall from 63 of 416 to 10 and ringed giants from 19 of 53 to 5; the median ring mass goes from 1.5e10 kg to 1e18. Stand-ins RO2 (changed) and RO9 to RO13 (new) in the registry; the old `RING_LIFETIME_YR` is gone. Dust rings fed by small moons (Jupiter's) are not modelled, a follow-up.

generation-output: changed

TODO items advanced: GEN.1171.
