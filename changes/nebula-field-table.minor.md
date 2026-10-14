### Added
- `nebula_field` (schema v83): one row for every kept class M and N cloud of the galaxy's molecular cloud field, written by `planetgen plan` (and by `plan --redo-scatters phenomena`) from the galaxy seed, so the Galaxy Map can find a nebula no generated sector holds yet. A sector's save marks the clouds it stores, with their object IDs. Classes P and Q stay in the seed (DB.24).

### Changed
- The migration to v83 only creates the empty table; the next plan fills it.
