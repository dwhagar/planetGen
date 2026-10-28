**Headline:** A rocky world now has one tectonic regime, plate tectonics or a single rigid plate, and its interior heat, its magnetic field and its ground radiation all follow it.

### Changed
- **One regime, kept with the body (GEN.721).** Whether a rocky planet or moon has a mobile lid (plate tectonics) was decided twice by separate coin flips on the same chance, once by its magnetic field and once by its ground dose. It is now drawn once, after the hydrosphere, at the same chance (a liquid surface ocean and 0.5 to 5 Earth masses at 0.6, heavier at 0.3, any other body stagnant), stored as `lid_regime`, and shown on the body's page as "Tectonic regime". The dynamo's lifetime, the ground dose's felsic crust and the heat law all read it.
- **The heat of a mobile lid.** `interior_heat` has a mobile-lid law next to the stagnant one: the plates carry the whole drop from the mantle to the surface, so at one mantle temperature a mobile lid loses about ten times the heat. An Earth twin ends at a mantle temperature of about 1,430 K with a mobile lid and 1,950 K with a stagnant one (the Earth's potential temperature is near 1,600 K). The mobile exponent (1/3) and critical Rayleigh number are stand-ins, labelled in the design note.
- **Stored (schema 169).** `planets` and `moons` gain `lid_regime` (NULL for a gas giant and for a row made before this release, which keeps the heat flux of a stagnant lid).
- Generator epoch bump (numbered at merge): the regime is drawn in the place of two other draws, so the run's later draws move, and every rocky body's dynamo, ground dose and heat flux can change.

generation-output: changed
TODO items resolved: GEN.721. Added: none.
