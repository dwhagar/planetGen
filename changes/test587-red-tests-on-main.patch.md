**Headline:** Four more tests that were red on the main branch pass again, with no change to what the generator makes.

### Fixed
- `test_sectors_table_rejects_duplicate_address` (TEST.587) asks the table directly: `insert_sector` has filled the row kept at an occupied address since GEN.1127, so it no longer raises; a new test holds that behavior.
- The order test of the bright-stars backfill uses its own seed: with the old one the 200 cells held 4 stars where 9 are expected, too few to tell one order from another.
- The burn-down register has lines for the sky and glow stand-ins that were added without them (GLW1, SKT1, SKT4 to SKT6, SKD1, SKD2, SKD4, SKD5); the SKT1 line named VIEW.41, which is closed.
