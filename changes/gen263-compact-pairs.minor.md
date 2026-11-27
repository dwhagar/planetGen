**Headline:** Bound pairs of compact objects (two neutron stars, two black holes, a black hole with a neutron star) are now rows of the scatter whatever their mass, so they can be seen on the map before their sector is filled.

### Added

- Two scatter classes, `neutron-star` `pair` and `black-hole` `stellar-pair` (GEN.263, multi-star-systems.md 4.5): about 1.7e5 neutron star pairs and about 2e5 black hole pairs in the default galaxy, drawn below the mass cut where the classes' own rows are not scattered. Each row stores a member tree with a compact companion (`scatter_companions._pair_companion`); the Galaxy Map and the Phenomena table list them under their own descriptors.

### Changed

- A sector's fill leaves out the neutron stars and black holes that are members of a promoted pair (`phenomenon_scatter.fill_share`, about 0.04 percent of the neutron stars and 0.45 percent of the black holes below the cut), so the objects of a class still add up to its population.

generation-output: changed
generation-affects: fill
