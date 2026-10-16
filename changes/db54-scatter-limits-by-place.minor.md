**Headline:** The scatter limits by place are three settings of the galaxy: a core luminosity floor and mass limit and a halo luminosity floor.

### Added
- **`planetgen plan --core-min-luminosity`, `--core-mass-limit` and `--halo-min-luminosity` (DB.54).** Stored on `galaxy_shape` (migration 0101, schema v101): the core's luminosity floor and mass limit, unset until given and then the baseline's times 1.1 (9,900 and 15.4 at the defaults), and the halo's luminosity floor (default 500, from 50 to 4,000,000; the halo has no mass limit of its own). A floor under 50 is refused with the reason, and `plan` says what the three are. Nothing reads them yet: the tiers that choose a limit by place are GEN.333.

TODO items resolved: DB.54. Added: none.
