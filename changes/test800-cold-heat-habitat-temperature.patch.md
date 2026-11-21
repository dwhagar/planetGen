**Headline:** The cold and heat adaptation traits of a life world now read the temperature of the water the life lives in, so a frozen world with warm crust no longer draws a cold trait.

### Changed

- TEST.800: the cold and heat stress of `adaptation_draw.stresses` reads `life_supply.habitat_temperature_k` (the surface under an open ocean or on a dry surface, the melting point of ice under an ice lid, the rise with depth in an aquifer) and not the mean surface temperature. In a 40-system sample the old reading gave a cold trait to 95 of 143 life worlds (94 of 137 chemical-class worlds, whose crust is warm under a frozen surface); the new one gives it to 13 of 144. The radiation, sulfide and carbon dioxide traits and their draws are unchanged.
- The life census test lists all five stresses (cold and heat came with GEN.1725) and the `any_trait` band is set again from a fresh census.

generation-output: changed
generation-affects: fill
