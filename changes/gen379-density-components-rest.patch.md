**Headline:** The galaxy's density is now the sum of named components, with each component's mass in solar masses, and the stored family and component numbers are read back.

### Changed
- **The density model is one list of components (GEN.379).** `relative_density`, `population_densities` and the two array forms take each component's value from one place (`galaxy.density.component_values` for a point, `component_arrays` for arrays of points): the thin disk with its arms, the thick disk, the bulge or bar, the stellar halo and the nuclear cluster and disc. A component with amplitude 0 is exactly 0. Every density is the same to the last digit as before (a fixture of six shapes at 24 points checks it) and the array forms are as fast.
- **`component_masses(shape)` now returns each component's mass in solar masses**, by name, from the closed forms; the six sum to the galaxy's stellar mass. The numerical integral it used to be is `integrated_masses`, which the math check and the tests use.
- **The stored family and component numbers are read and written (DB.66).** `get_galaxy_shape` puts `family`, the bar's presence, axis ratios and angle, the thick disk's amplitude, the spheroid's numbers and the clump share into `GalaxyShape`; a NULL bar or thick disk number means the structure draw's own value, a number pins it. `save_galaxy_shape` writes them, replacing what an earlier plan stored. A stored spheroid amplitude above 0 is refused until the spheroid component is built.

generation-output: unchanged
