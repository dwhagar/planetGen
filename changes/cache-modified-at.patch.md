### Changed
- **Editing a sector no longer throws away the whole Galaxy Map cache.**
  The tile cache used to go stale on any change to the database. It now
  asks the new `GET /api/galaxy/changes`, which reads the schema-v27
  `modified_at` columns plus new sector and system ids, which cube tiles
  changed, and deletes just those, on the server's disk and in each
  visitor's browser. A rename refetches the dozen tiles holding that
  sector. Deleting a sector, re-planning the galaxy or a new release
  still refreshes everything, since a deleted row leaves nothing to
  locate its tiles by.
