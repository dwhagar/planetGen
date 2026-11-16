**Headline:** the API test for deleting a system expects what the API answers, a 410 for the deleted system's ID.

### Fixed

- TEST.781: `test_update_and_delete_system` expected a 404 for a system after `DELETE /api/systems/<id>`. Since GEN.1141 a deleted object's ID is kept (a tombstone), so the API answers 410 Gone for it (`designation-history.md` invariant 7; `api/ids`, `test_api_retired_ids`, `test_admin_edits` say the same). The API was right and the test was left from before; it now expects the 410 and keeps the 404 for an ID nobody ever had.

generation-output: unchanged
