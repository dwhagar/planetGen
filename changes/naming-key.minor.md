### Added
- The galaxy has a naming key (8 hex digits) in the control database, drawn from the galaxy seed when `planetgen plan` makes a new seed and changeable by an admin on the Stats page or with `POST /api/admin/naming-key` (control schema v9, GEN.70). Run `update.sh` to add the table.
