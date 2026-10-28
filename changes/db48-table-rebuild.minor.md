**Headline:** A table whose clustered index has a page the server cannot read can now be rebuilt in two steps with `planetgen repair-db --salvage` and `--reload`, and the rows the page took come back from the parity file.

### Changed
- **`planetgen repair-db --salvage TABLE` and `--reload TABLE` (DB.48).** With the server started with `innodb_force_recovery = 1` and the website and generation stopped, `--salvage` dumps every readable row of the table to a file (exactly, nothing rounded), then drops the table and makes it again empty. After a normal restart `--reload` inserts the rows back with their ids, checks the count, and runs the sector repair, which rebuilds the sectors that lost rows from the parity file. It refuses, changing nothing, when the table cannot be read and the server is not in force recovery, when a dump is already waiting, or when a normal server counts other rows than were dumped. An interrupted reload is finished by running it again. `--restore-sectors` (with `--reload sectors`) also puts back sector rows the parity file holds and the table lacks. `repair-db` now names the unreadable tables when it lists unreadable sectors. Only the tables a sector's checksum covers can be rebuilt this way.
- `repair-db` can rebuild a sector whose own `sectors` row is gone (it puts the row back first).

TODO items resolved: DB.48. Added: none.
