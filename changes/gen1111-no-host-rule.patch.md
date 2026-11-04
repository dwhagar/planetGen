**Headline:** A planet with no star now keeps the air its own gravity and heat allow, and every rocky planet's temperature includes its own heat as well as starlight.

### Changed
- **The air of a body with no light (GEN.1111).** The atmosphere table's shoreline rows (d1 to d5) measure a body against the starlight and the luminosity of its host; a body with no host has neither, so the star cannot strip its air and it sits far inside the shoreline (`NO_LIGHT_DELTA`, a sentinel, not a tuned value: any large negative number gives the same row). The star's wind pressure at a body with no star is zero, so it has no magnetopause to hold off a wind that is not there (`magnetism.magnetopause_rp` gives `None`). These are what an ejected planet needs to be built through the one planet chain; nothing here changes a lit planet.
- **The planet's own heat warms it.** A rocky planet's temperature, with or without air, now adds the heat of its interior and its tide (`internal_heat_flux_w_m2` and `tidal_flux_w_m2`) to the flux it radiates, as a giant's already did. Under a star the shift is a few hundredths of a kelvin for an Earth-like heat flow, and a body far from any star settles at the temperature its own heat gives rather than the 2.725 K floor. This moves every rocky planet's temperature slightly, so the generator output changes (the epoch is taken at the train, --affects fill).

generation-output: changed

TODO items resolved: none (GEN.1111 stays open).
