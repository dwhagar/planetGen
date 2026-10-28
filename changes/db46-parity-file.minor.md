**Headline:** An optional parity file beside the database can now rebuild a damaged sector: it is kept current as sectors are saved, edited and deleted, and `planetgen refresh-parity` builds it for a galaxy made before.

### Changed
- **The parity file (DB.46).** With the new setting `parity.enabled` on (off by default; `parity.dir`, `parity.group_size` 32 and `parity.parity_blocks` 2 are the others), a log of Reed-Solomon parity records is kept outside the database for groups of up to 32 sectors of about the same size that sit far apart in sector id, so that two damaged sectors of a group can be rebuilt from the rest (the rebuild is DB.47) at about 7 percent of the sectors' size in disk space. Every save, edit and deletion commits first and then updates the parity (`Connection.after_commit`); a sector whose parity is behind is reported as unprotected by `check-db`, never as damaged, and a record behind the database is made again, not treated as damage. A sector whose group has a damaged member is left alone.
- **`planetgen refresh-parity [--compact] [--limit N]`** makes the file current for the whole database (marked and never-sealed sectors, sealed sectors the file lacks, members whose sector is gone) and `--compact` drops the dead space an append-only log leaves.
- `sector_export` writes a sector as bytes (every covered row, ids included, in the checksum's own encoding) and reads them back, so a rebuilt block can be turned into rows with their ids; `storage_checksum.checksum_of` takes rows from either.
- The design note (`docs/design/operations/db-check-and-parity-repair.md`, 5.10) records what was built and where it differs from the sketch.

TODO items resolved: DB.46. Added: none.
