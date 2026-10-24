**Headline:** The admin status page and the API health check now show the generator version, the one number that says which rules made a galaxy.

### Added
- **The generator version (OPS.37).** `GET /api/health` and the database figures of `GET /api/admin/stats` carry `generator_version`, the current generator epoch (OPS.28), and the admin Stats page shows it as "Generator version" beside the schema row. It is the same counter as the generator epoch, not a second one; an added field, not a breaking change.
