**Headline:** Resetting the galaxy names the table it is emptying in the beat, so a slow TRUNCATE shows what it waits on.

### Changed
- **Reset-galaxy steps by table (OPS.206).** Each `TRUNCATE` of the wipe goes through `queue/longsql.run` with the table, its row count and its place in the list ("Emptying scatter_masses (41,000,000 rows), table 31 of 120").

generation-output: unchanged
generation-affects: none
