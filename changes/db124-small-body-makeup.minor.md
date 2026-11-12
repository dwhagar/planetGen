**Headline:** A belt or comet read back from the database keeps its drawn makeup (porosity, density, albedo and the material shares).

### Changed
- **DB.124 (bug).** `asteroid_belts` and `comets` gain a `makeup` column (JSON text, migration 0247, placeholder number for Lane 02 to renumber) that holds the drawn `SmallBodyMakeup`; `store.BELT_FIELDS` and `COMET_FIELDS` write and read it, and the Details sheet notes it as not shown. Before, a loaded belt or comet had the words and the chondrite group but an empty makeup, so a loaded body could not say its porosity, density or albedo. A body saved before this migration reads back with no makeup, as before.
The drawn values are the same, but the stored rows of belts and comets gain a column, so this rides a train with its migration.
generation-output: changed
