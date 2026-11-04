**Headline:** Every belt now keeps its biggest bodies one by one: where each lies, how its orbit is shaped, what it is made of and whether it counts as a planet, a dwarf planet or a minor body.

### Changed
- **Belts store their members (GEN.1070, schema v205, migration 0205).** A new `belt_members` table holds, for each belt, up to 300 of its heaviest bodies, ranked from the heaviest: mass, radius, distance, eccentricity, inclination, the three orbit angles, the taxon (spectral class), bulk density, the rock, metal, water and ice shares, the IAU status word and its Soter mu, and a scatter tag that stays empty until the scatter step tags members. The masses are the ones settled with the belt (the bodies above 1e19 kg, up to 90 percent of the belt's mass), so they add up to no more than the belt weighs. A belt read back from the database lists the same members as one just generated, and a belt cut back to the orbits its star holds draws its members again for the shorter belt.
- **Where each member lies.** The distance of a member is drawn from where the belt's solids are (more bodies where the disc held more mass), and its eccentricity and inclination are drawn so a wide belt is more stirred than a narrow one (the rms eccentricity is the belt's width over its middle, with an inclination of the same size in radians). The angles are uniform. The draw is keyed on the belt and its star, so the heaviest bodies keep their rank and their place when the belt is cut back or its limit changes.
- The belt's member list is not yet shown on the belt page or in the API (that comes with the belt page items), so the API version does not change. The four stand-ins of the placement are registered in the stand-in list with their tier and tweak level (MR1 to MR3).
- Wiring register: lines for modules and numbers the belt generator now reads are removed, and the lines for the roundness and resonance code that is still waiting now name the open items that will read it.

generation-output: changed

TODO items resolved: none.
