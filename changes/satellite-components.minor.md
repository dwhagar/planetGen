**Headline:** The Magellanic Clouds, Sagittarius and eight dwarf galaxies are in the galaxy, at their real places.

### Added
- **Satellite galaxies (GEN.237).** `planetgen plan` places the Milky Way's eleven satellites (the Large and Small Magellanic Clouds, the Sagittarius dwarf, Fornax, Sculptor, Carina, Draco, Ursa Minor, Sextans, Leo I and Leo II) as `satellite` components, each as an exponential or Plummer profile cut at six half-mass radii, at the direction and distance the sky gives them (the Large Cloud 49.6 kpc away and holding about 6.2e9 systems). The sector gate counts their density, so the outline reaches them and a cell at the centre of the Large Cloud is a satellite place, not empty halo. Their sizes, masses and distances are recalled from McConnachie (2012) and need checking.
- **A `satellite` stellar population.** The sector profile's mix, region label, ages (0.5 to 13 Gyr), dispersion and rotation now have it, and its [Fe/H] is the one of the satellites that reach the place.
- **[Fe/H] from the component.** A star in a cluster or satellite draws its [Fe/H] from the peaks of the place (`metallicity.sample_fe_h(peaks=)`, `SystemConfig.FE_H_PEAKS`), which the fill sets; before, filling a cluster's cell raised an error.

### Changed
- **One field for the clusters and the satellites.** `ClusterField` is `ComponentField` (`galaxy/component_field.py`) and `GalaxyShape.clusters` is `GalaxyShape.components`; every component needs a `tidal_radius_pc`, and `king_w0` only for a King profile.
- **A seeded galaxy's outline is different from before**, grown by the Clouds, Sagittarius and the dwarfs; plan the galaxy again.

TODO items resolved: GEN.237. Added: none.
