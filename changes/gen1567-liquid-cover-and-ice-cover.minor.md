**Headline:** A world's ice is now counted apart from its sea, and a deep ocean no longer drowns the land at a few kilometres: the share of a planet's surface that is liquid ocean no longer reads 100 percent on a frozen world, so the ocean and desert classes and the weathering land see the sea and the land that are really there.

### Changed

- GEN.1567: `ocean_fraction` is the share of the surface under liquid water (a surface ocean, a brine seep, a hycean sea, or the sea under an ice lid) and no longer counts ice; the new `ice_cover_fraction` (migration 0244, `planets`, `moons` and `rogue_planets`, DOUBLE, NULL for a body with no ground) holds the share the ice covers (an ice lid, or a body that is ice through). `land_fraction` is the rest of the surface after the liquid sea, so a body that is ice through (`ice`) now stores 0 sea, its ice, and 1.0 land (the ice lies on land) where it stored 1.0 ocean and 0 land (census bug 9: 1,877 bodies).
- Readers that asked whether a surface is wet at all keep that: the ground shares (`surface_composition`) take the sea for a liquid world and the ice for a frozen one, the haze draw (`aerosols`) and the cosmogenic cascade (`isotope_record`) count either. The class fits read the liquid sea, as the water-land plan asks. The sheet shows "Surface under liquid ocean" and "Surface under ice".
- Planets that had a frozen surface change their stored sea, land and ground-dose air flag (`radiation`: an ice world with an air now has land to take radon from).
- GEN.1573: the ocean crust of the cover law sinks under the water it carries (Cowan and Abbot 2014, Eq. 10): its mean depth is the water layer over its share of the area where that is deeper than Earth's 4.5 km, up to 11.4 km divided by the gravity in Earth's, so the continents stay dry through about 7 km of water at Earth's gravity where the old law covered them at 4 to 5 km. Earth's 2.7 km layer still covers 0.71. Icy bodies keep their strength-limited relief. Cover changes only for bodies with more than Earth's water: those with roughly 3 to 7 km of liquid water (at 1 g) move from nearly or fully covered to about 0.71.

generation-output: changed
generation-affects: fill
