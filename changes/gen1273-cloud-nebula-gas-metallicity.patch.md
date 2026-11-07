**Headline:** A nebula inside the Large or Small Magellanic Cloud now has its gas metallicity from the Cloud's own [Fe/H], so it is more metal-poor than a nebula at the solar circle.

### Changed

- **A Cloud nebula takes the Cloud's metallicity (GEN.1273).** `nebula_columns._gas_oh` never passed a component's [Fe/H] to `nebula_composition.gas_oh`, so every nebula, including the ones GEN.367 places in the Clouds, used the Milky Way gas gradient at its galactic radius (a Cloud nebula, 50 kpc out, got the gradient's far-disk value regardless of the Cloud). It now asks the bound galaxy's components for the satellite whose gas holds the nebula's place (`ComponentField.gas_fe_h`, the satellite with the largest gas tracer there) and draws the gas [O/H] inside the widest range of the two methods above that [Fe/H] (`SATELLITE_GAS_OFFSETS`, registry row NBC3, already a GEN.725 row with tier T2 and tweak level L2): the Large Cloud's nebulae come out at [O/H] -0.34 to -0.14, the Small Cloud's at -0.66 to -0.28. Nebulae of the Milky Way, nebulae with no center, and galaxies with no satellites are unchanged. Supernova remnants keep the gradient (a remnant stores no center yet). The stored composition of a Cloud nebula changes, so this is a lock item; the golden galaxy has no Cloud, so its digests do not move.

generation-output: changed
