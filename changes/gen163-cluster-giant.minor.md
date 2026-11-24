**Headline:** A millisecond pulsar of a globular cluster can now hold a captured giant like PSR B1620-26 b, and every scattered object records which cluster it was drawn in.

### Added

- `phenomenon_scatter.cluster_id` (GEN.163, schema v275): the `galaxy_components.id` of the globular cluster a scattered object was drawn in. The component pass writes it; every other scatter writes NULL (field objects, satellites, rows made before v275). A built neutron star carries it (`NeutronStar.cluster_id`).
- The captured-giant origin of a pulsar's planets is built, for a pulsar of a cluster only: one planet of 1 to 4 Jupiter masses at 20 to 50 AU, both log-uniform (`tuning.PULSAR_GIANT_MASS_JUPITER_RANGE`, `PULSAR_GIANT_DISTANCE_AU_RANGE`), made by the ordinary planet classes (a gas giant, with the moons of its class). A field pulsar that draws the origin stays bare, as the design rule says.

### Changed

- The bare-pulsar stand-in for the captured giant is gone (`PULSAR_PLANET_BUILT_ORIGINS` now holds `disk` and `captured-giant`); the ablated companion stays unbuilt until the pair builder (GEN.555).

generation-output: changed
generation-affects: fill
