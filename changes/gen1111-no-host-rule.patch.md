**Headline:** A planet with no star can now be built through the one planet chain: it keeps the air its own gravity and heat allow, and it is warmed by its own heat instead of sitting at the 2.725 K floor.

### Changed
- **The air of a body with no light (GEN.1111).** The atmosphere table's shoreline rows (d1 to d5) measure a body against the starlight and the luminosity of its host; a body with no host has neither, so the star cannot strip its air and it sits far inside the shoreline (`NO_LIGHT_DELTA`, a sentinel, not a tuned value: any large negative number gives the same row). The star's wind pressure at a body with no star is zero, so it has no magnetopause to hold off a wind that is not there (`magnetism.magnetopause_rp` gives `None`).
- **A body with no star is warmed by its own heat.** The temperature of a rocky planet whose host is starless, with or without air, adds the heat of its interior and its tide (`internal_heat_flux_w_m2` and `tidal_flux_w_m2`) to the flux it radiates, as a giant's already did. Under a star that heat is below a part in a thousand of the absorbed flux for an ordinary rocky planet and is left out, so no lit planet changes. (A strongly tide-heated moon is the case where it is not that small; that is a separate question, not changed here.)

generation-output: unchanged

TODO items resolved: none (GEN.1111 stays open).
