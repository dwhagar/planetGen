**Headline:** The sector-delete test now matches the rule that a sector holding a filed facility cannot be deleted until the facility is removed.

### Fixed
- **test_deleting_a_sector_with_facilities was stale.** Since GEN.1127 a facility filed in a sector stops the sector's delete with a 409 (`facilities.sector_id` is `ON DELETE RESTRICT`); the test now checks the refusal, deletes the facility, then the sector.
