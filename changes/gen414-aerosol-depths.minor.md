**Headline:** Every planet and moon with air now has a cloud, haze and dust depth drawn from its own air and temperature, stored for the radiation model.

### Added
- **GEN.414 and DB.68: haze, dust and cloud depths.** `physics/aerosols.py` draws, for each planet and moon with air, the in-cloud optical depth and cover (the condensate that a dry adiabat to the skin temperature can put in the air: water, ammonia, methane or carbon dioxide, droplets 1 to 20 micron, depth log-uniform from 0.1 to that limit), the methane haze (the CH4/CO2 gate of 0.1 to 0.2, the diffusion-limited supply times a yield of 0.01 to 1, settled through the air column; a Titan-like body comes out between 0.02 and 3) and the dust of a dry ground (the wind ratio log-uniform on 0.3 to 3, depth 0.5 to 1.5). `layers` turns them into the grey cloud, the tholin haze and the Wolff dust spectra the band model of GEN.443 takes. Migration 0138 adds `tau_haze`, `tau_dust`, `tau_cloud`, `cloud_cover` and `aerosol_draw` to `planets` and `moons` (schema v138); a row saved before it keeps NULL until the galaxy is regenerated.

### Changed
- **Generator epoch.** The bodies now carry five more stored values, drawn on their own stream so no other draw moves; the golden record is regenerated for it.
generation-output: changed
