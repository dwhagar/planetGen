**Headline:** The speed at which an ejected neutron star streams out of a galaxy now follows its natal kick and the galaxy's escape speed.

### Changed
- An ejected neutron star's far speed is what its kick leaves above the escape speed at its birth radius (`galaxy.natal_kick.ejecta_speeds_kms`), kept as 16 quantiles in the galaxy's structure (`Structure.neutron_star_ejecta_speeds_kms`), instead of a speed uniform from 300 to 500 km/s. For a galaxy of the Milky Way's structure the mean is 386 km/s; a deeper well or a slower kick gives slower ones, and a galaxy that ejects none has no ejecta population. The population's density and outer radius take the mean speed (`ejecta_field.mean_speed_kms`); ejected black holes keep the 300 to 500 km/s range, a stand-in until they have a kick model. `docs/design/natal-kick.md` records the model.

generation-output: changed

TODO items resolved: GEN.554 (the ejected neutron stars' speed; GEN.554 stays open for the neutron star mass function and the placement of kicked remnants) / Added: none
