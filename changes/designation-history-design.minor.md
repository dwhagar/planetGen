**Headline:** The design for objects that change name and ID while keeping their past designations is written.

### Added
- **Designation history design note (GEN.1121).** `docs/design/architecture/designation-history.md`: one `designations` table of retired names and IDs (a retired ID is reserved forever, a retired name is kept in history and may move by the collision rule), lookup of an old ID redirecting to the current object (a write against a retired ID refused with 409, a deleted object 410), bodies, belt members and captured bodies following a parent's ID change, a "formerly" line on the object page and an API field, the split into parts and three questions for Boss. Nothing is built.

No generator output changes.

generation-output: unchanged

TODO items resolved: none (GEN.1121 stays open for the build).
