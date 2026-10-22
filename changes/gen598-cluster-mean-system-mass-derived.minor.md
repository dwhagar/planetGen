**Headline:** The systems of a cluster or a satellite galaxy follow from a mean system mass worked out from its own age, not a recalled 0.4 solar masses.

### Changed
- **A component's mean system mass is derived (`Component.mean_system_mass_solar`).** `tuning.CLUSTER_MEAN_SYSTEM_MASS_SOL` (0.4, "recalled, good to about a factor 1.5") is gone. A globular cluster, a young open cluster or a satellite galaxy holds its mass over the mean mass of a system whose stars formed `age_gy` ago (`generation.system_mass.mean_system_mass_at_age_solar`, the way GEN.545 does the galaxy's): 0.47 solar masses at 11 Gy, 0.49 for the Large Cloud, about 0.9 for a cluster a few Myr old. A cluster of 2e5 solar masses at 11 Gy holds 4.3e5 systems (it held 5.0e5), the Large Cloud 5.1e9 (6.2e9), so every seed's clusters and satellites hold about a sixth fewer systems and place their stars at that density.
  - `Component.expected_systems`, `generation.component_pass.systems_of` and `galaxy.component_field.ComponentField` read it; `physics.king.expected_systems` takes the mean as an argument.
  - `generation.open_clusters.ambient_density_msun_pc3` multiplies by the galaxy's own mean system mass (`density.mean_system_mass_solar`) instead of the cluster constant, because it is the density of the field stars around the clusters.
  - `docs/design/globular-clusters.md` and `intergalactic-and-halo-contents.md` carry the new exemplar and satellite counts.
- The generator epoch is 49.
generation-output: changed

### Added
- `tests/test_gen_component_system_mass.py`: a component's mean against seeded systems at four ages (the generator's own samplers end to end), the systems as the mass over the mean, and that no constant holds it.

TODO items resolved: GEN.598. Added: none.
