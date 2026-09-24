### Added
- **Admin stats page.** Logged-in admins get a new Stats page (sidenav
  and a link on the Admin page) showing server health (API version,
  uptime, load, memory, MySQL version/uptime/connections, galaxy tile
  cache usage and free disk) and stats about the current database: exact
  sector and system counts, size on disk, schema version, when rows were
  last created or modified, and per-table row estimates and sizes.
- **Names made unique.** The same page counts every name the uniqueness
  rules had to decorate (Alpha/Beta..., Little..., ...Kin) and lists each
  one with links to every sector and system carrying it; planets and
  moons link to their system.
- New admin-only endpoints `GET /api/admin/stats` and
  `GET /api/admin/duplicate-names` (see `docs/api.md`).
