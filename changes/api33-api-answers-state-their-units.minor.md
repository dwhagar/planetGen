**Headline:** The sizes in the backup list now say they are in bytes, and a test checks that every number the API returns names its unit or says what it counts.

### Changed
- **Backup sizes are `size_bytes` in the API.** `GET /api/admin/backups` and `GET /api/admin/backups/<name>` return `size_bytes` for an archive and each of its parts where they returned `size`. The Backups pages read the new name and look the same. The API version is now 11.
- **A test now reads real API answers for numbers with no unit.** Every number in the answers for sectors, systems (with their planets and moons), search and the galaxy routes must end in a unit, be a stored column with a unit that says what it is, or be listed with what it counts or names. A new field with a bare number fails it.

api-version: changed
