**Headline:** A database older than the oldest structure 8.0 accepts is refused with a message that points to legacy/bump_database_to_8_0.py.

### Changed
- **DB.127: Boss's hard cut (3:53 PM EST 2026-11-12): no database structure before 8.0 is accepted.** The live chain already starts at the v130 baseline (2026-10-24, inside 8.0), so no revision is deleted or renumbered. `SchemaTooOldError` now tells the operator to run `legacy/bump_database_to_8_0.py` (DB.126) to keep the data, or to drop the database and let `update.sh` make a new one, instead of pointing to an earlier checkout. The old-schema fixtures and the migration tests below the current few were retired by TEST.729, TEST.731 and TEST.733. No migration, no epoch, no generated output moves.

generation-affects: none
