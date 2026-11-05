**Headline:** Giants now keep faint dust rings fed by small inner moons, so most giants have a ring again instead of one in fifty.

### Changed
- **GEN.1233.** The ring lifetime law of GEN.1171 correctly kills a ring that nothing feeds, but dust rings that small moons keep supplying were not modelled. A giant now has a Poisson (mean 6) number of small moonlets inside 3 planet radii, drawn at ring time and not stored as moons (the generator makes no moon under a few hundred km), each 5 to 100 km. The meteoroid flux on their cross-sections, times a kept-dust yield of 10 and a dust life of 1e4 years, gives a steady ring mass, origin `moon_erosion`, added to the rings of fallen moons and comets. Jupiter's four inner moons come to 2.4e7 kg, inside the 1e7 to 1e9 kg of its ring. Over 60 systems the giants with a ring go from 1 of 52 to 39 of 52 (rocky planets are unchanged at 10 of 415); the median ring mass of ringed bodies is 3e8 kg. The draws change the unit's random stream after the ring step, so later values of a giant's unit change too. Stand-in rows RO14 to RO18 in the registry.

generation-output: changed

TODO items advanced: GEN.1233.
