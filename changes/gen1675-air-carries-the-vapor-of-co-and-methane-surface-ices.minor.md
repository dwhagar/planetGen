**Headline:** A cold body's air now carries the carbon monoxide and methane its surface ices give off, at the saturation pressure of the ice at the surface temperature.
generation-output: changed
generation-affects: fill

### Changed
- **Surface ice vapor in the air (GEN.1675).** `surface_air.solve` treated only water as a reservoir the air stands in contact with, so a cold body with CO or methane ice held an air with none of that vapor (or none at all). The carbon monoxide and methane ices of the body (`ice_co_mass_fraction`, `ice_ch4_mass_fraction` times its mass) are now reservoirs of `freeze_out` at full humidity (`ICE_HUMIDITY` 1.0, equilibrium with a solid), present only below the triple point of the ice and used up when small. An empty column over a humid reservoir fills to the saturation pressure, a body whose air is only ice vapor is no longer returned airless early, and a gas standing over its own ice is not also held in the clathrate. It reuses the GEN.1564 saturation table and freeze-out; no second path. Nitrogen ice is not included.
