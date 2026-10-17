**Headline:** The globular clusters and the satellite galaxies now place their own bright and massive stars and their black holes, neutron stars and planetary nebulae.

### Added
- **The component pass (GEN.324).** A component's systems are its mass over the mean system mass; each pass draws one Poisson count from them times the share of its population the pass keeps, puts every candidate where the component's profile says (a radius with that share of the mass inside it, a direction uniform on the sphere), and keeps it by its tier's share. The same is done for the black holes, neutron stars and planetary nebulae by their rates per system. This is what a sector's fill already takes off the component's cells, so no pre-placed star is drawn twice. Each component has its own random stream.

### Changed
- `bright_stars.population` takes `satellite` (migration 0105); a galaxy scattered before this has no stars in its clusters or satellites: `planetgen plan --redo-scatters luminosity mass phenomena` draws them.

generation-output: changed

TODO items resolved: GEN.324. Added: none.
