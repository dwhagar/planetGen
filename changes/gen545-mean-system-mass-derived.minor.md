**Headline:** A galaxy's number of systems follows from a mean system mass worked out from the generator's own star mix, not a recalled 0.4 solar masses.

### Changed
- **The mean mass of a system is derived (`generation/system_mass.py`, `galaxy/system_mass.py`, `density.mean_system_mass_solar`).** `tuning.MEAN_SYSTEM_MASS_SOLAR` (0.4, "recalled, within a factor of 1.5") is gone. For stars formed at one age the mean is the placed star's present mass from the initial mass function (its birth mass until its last giant phase ends, then its white dwarf, neutron star or black hole) plus the other stars of its system: the `star-count` table's mean number by the placed star's class with the 8 to 12 tail (`multistar.expected_star_count`) times the mean present mass of a star born with 0.1 to 1 of the primary's mass. A population's mean averages the ages by its star formation (the disk's history; the other windows uniformly), and a galaxy's or cell's averages the populations by their systems. The reference galaxy comes out at 0.483 solar masses a system (young 0.71, intermediate 0.54, old 0.48, bulge 0.47, halo 0.47: a placed star of 0.35 and 1.42 stars a system), so it holds 1.03e11 systems where the constant gave 1.25e11, and the density at the reference radius is 0.62 of the pinned local density (0.750 after GEN.544).
  - `galaxy.density.build_galaxy_shape` divides the stellar mass by it; `galaxy.resolution` (the mass of a cell's unresolved systems), `generation.place_terms.from_profile` (a place's mass density) and the math check read the cell's or galaxy's own mean.
  - The initial mass function integrals and the cohort moments moved from `generation/place_terms.py` to `generation/system_mass.py`; the perturber mass is unchanged. The empty-cell fallback of `mass_weighted_perturber_mass_solar` is the disk's mean system mass.
  - `galaxy.disk_history.history_of` ignores `k_norm`, so the history computed while the density is normalized is the one every later read of the shape uses.
  - `multistar.MAX_MEMBERS`, `multistar.COMPANION_MASS_FRACTION` and `multistar.expected_star_count` name what `draw_member_tree` already did.
  - `CLUSTER_MEAN_SYSTEM_MASS_SOL` (the clusters' and satellites' 0.4) is a separate constant and stays.
- The generator epoch is 46.
generation-output: changed

### Added
- `tests/test_gen_system_mass.py`: the expected star count against the generator's own draw, the placed star's mean against a seeded initial mass function draw, the other stars against real member trees, and the whole against seeded systems at three ages.

TODO items resolved: GEN.545. Added: none.
