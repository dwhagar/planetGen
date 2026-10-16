**Headline:** A metal-poor star's disk makes far fewer planets, and the sector profile scales each rare kind by the place it is in.

### Added
- **Kind factors in the sector profile (GEN.236).** `SectorProfile.kind_factors` holds, per kind, the factor by which the place scales the kind's rate per star: neutron stars, black holes, planetary nebulae and brown dwarfs 1 everywhere; supernova remnants the gas fraction (0 in the halo and in clusters, and not zeroed in a superbubble); molecular clouds the cold gas fraction times the young share against the solar circle's; the nebulae grown around hot stars the cold gas fraction; rogue planets and comets the planet occurrence at the place's [Fe/H]. A sector fill multiplies each rolled kind's mean by its factor, and the phenomenon scatter weights supernova remnants by the gas (`nebula_field.gas_factor_array`).
- **The planet cull by metallicity.** Planet occurrence is unchanged at and above [Fe/H] -0.5 and falls ten-fold per dex below it (`tuning.PLANET_OCCURRENCE_FLOOR_DEX_BELOW_SOLAR`); gas giants fall faster (the giant-planet form of `giant_planet_factor`). A halo G dwarf now carries about 1.6 planets and 0.05 gas giants on average, where it carried 8.5 and 1.7, and no halo or cluster star holds more objects than its metallicity allows (`tuning.PLANET_OCCURRENCE_REFERENCE_OBJECTS`).

### Changed
- **A close binary's merged star now carries its primary's [Fe/H]**, which the planet cull and the disk's solids read (it had none, so a close pair counted as solar).
- **Fewer supernova remnants where there is little gas.** The scatter's remnants followed the stars; they now also follow the molecular gas, so the halo and the thick disk away from the plane hold almost none.
- **A seeded galaxy fills differently from before** in the places the factors reach (remnants, clouds, nebulae, rogue planets, comets, the planets of metal-poor stars).

TODO items resolved: GEN.236. Added: none.
