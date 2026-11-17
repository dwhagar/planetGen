**Headline:** The cache test expects a deleted system's ID to answer 410 Gone, as the tombstone rule says.

### Fixed
- `test_cache_invalidation::test_deleting_a_system_is_never_served_again` read 404 for the page of a deleted system; the ID of a deleted object answers 410 (GEN.1141), as `test_api::test_update_and_delete_system` already expects (TEST.781). Test only.

generation-output: unchanged
