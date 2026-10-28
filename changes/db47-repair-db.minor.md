**Headline:** `planetgen repair-db` rebuilds damaged sectors from the parity file, so a flipped value or a lost row in a sector can be mended instead of regenerating the sector.

### Changed
- **`planetgen repair-db [--shell N] [--sector SHELL BAND WEDGE] [--dry-run]` (DB.47).** Finds the sectors whose rows no longer hash to the checksum sealed for them and rebuilds them from their group in the parity file (up to two damaged sectors of a group at once with the default settings). The rebuilt block must hash to the stored checksum before anything is written, the rows are updated, inserted or deleted in one transaction, and the sector is checked again before the commit. The columns the orbit update rewrites are left alone on every row that still exists. A sector the parity cannot mend (a record behind the database as well, a sector in no group, more damage than parity blocks) is listed with the reason and left as it is, and the command exits 1; a table the server cannot read is DB.48's case. A parity file that is itself damaged is made again from the database. `--dry-run` writes nothing.
- A sealed sector that no longer hashes to its checksum is never taken into the parity file (`protect_sector`, so `refresh-parity` and the save hook leave it for `repair-db`): the sweep would otherwise have recorded the damage as the truth.

TODO items resolved: DB.47. Added: none.
