**Headline:** A world with water vapour in its air is no longer too hot: the vapour thins with height, as it does on the Earth, and the greenhouse sums over that profile.

### Changed
- **The water vapour of the greenhouse thins with height (GEN.812).** `physics/infrared_opacity.py` took the surface's vapour pressure as the mixing ratio all the way up, so an Earth twin whose air held 0.77 of saturation at its own temperature came out at 319 K, not 288 K. It now thins the vapour with height (the humidity falls linearly with pressure, Manabe and Wetherald 1967; the temperature follows the troposphere), sums the water continuum over the layers and counts the H2O bands and far wing for the pressure each layer sits at. The Earth twin lands at 288 K, and the surface temperature of every body with water vapour in its air changes.
- **A grey cloud depth for the longwave effect of clouds (GEN.812).** The module is a clear-sky greenhouse, but a body's albedo is its all-sky one, so one grey depth (0.35 at the Earth's surface humidity, scaled with the humidity, nil without liquid water) stands for the cloud. It is a stand-in with registry rows GH1 to GH4 in `docs/design/research/stand-in-values.md`.
- **`heat_capacity_ratio` moved** from `atmosphere_state` to `mixture`, where the vapour profile and the closure both read it.

generation-output: changed
