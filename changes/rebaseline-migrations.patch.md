**Headline:** The migration history now starts at schema v130, so a database older than that has to be dropped and regenerated.

### Removed
- The Alembic revisions 0061 to 0129 and the old `galaxy_containers` revision 0130, replaced by one empty baseline at v130; the revisions from 0131 on are unchanged.
- The 69 old-schema fixtures below v130 and the 36 per-revision migration tests (and parts of 4 other test files) that loaded them; the fixtures from v130 up stay (the missing v144 one is added), and `test_db_old_schemas` still upgrades each of them to the current shape.

### Changed
- A database below v130 is refused with `SchemaTooOldError` ("older than v130"). On Boss's test server, `python3 -m planetgen.cli.migrate --status` shows the version; at v130 or above `update.sh` upgrades it as before.
