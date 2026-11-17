**Headline:** A star's own disk time now decides when its planets are thrown out, when its cluster starts to unbind them, how long the pebbles drift, and whether a close companion has already ended its disk.

### Changed

- GEN.1551 slice 2: `scatter.DISK_LIFETIME_YR` (3 Myr, EJ2) and `cluster_unbinding.DISK_LIFETIME_MYR` (CU4) were one disk time for every star; slice 1 gave each star its own. `disk_truncation.disk_time_yr(star)` is that time (shorter for a star with a companion inside 47 AU, the same draw). The planet architecture records its losses at it (`Lost.disk_time_yr`, `planet_architecture.build(disk_time_yr=...)`), `disk_ejection` throws a disk's bodies out then, `cluster_unbinding.unbinding_epoch_myr` counts only the cluster encounters after it (the cluster's history now starts at its formation), and `disk_drift.disk_end_myr(mass_solar)` ends the pebble flux at it. A young star with a close companion that the shorter time has already taken past class II has its disk taken off once the system attaches the companions (`StarSystem._lose_disks_to_close_companions`; no draw). A top-up rogue with no star keeps the 3 Myr reference. Pebble drift reads the lone-star time (it holds no companions).

generation-output: changed
generation-affects: scatter, fill, population
