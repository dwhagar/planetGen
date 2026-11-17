**Headline:** A lone neutron star's surface temperature is now the cooling curve of its age, the same curve a neutron star in an orbit tree reads, so an old one is a few thousand kelvin and not a drawn 50,000 to 3,000,000.

### Changed

- GEN.1529: `NeutronStar` read its surface temperature from a range of 5e4 to 3e6 K at every age (`tuning.NEUTRON_STAR_SURFACE_TEMPERATURE_RANGE_K`, gone); it is now `dark_light.surface_temperature_k` of its age since the collapse, so a lone neutron star and one in a tree agree (4.6 Gy: about 4,600 K; 0.9 Gy: about 10,300 K). Its luminosity, X-rays and habitable zone follow. The plausibility check holds the temperature to the curve (the record carries `age_gy`). The snow line of its disk is read at the light the star had at the end of its neutrino era (`_formation_luminosity_w`, about 1e26 W, stand-in DL8), because the light of an old neutron star today put the snow line at 1e-5 AU and left a system that asks for two bodies with one; the line stays where it was, about an astronomical unit.

generation-output: changed
generation-affects: scatter, fill
