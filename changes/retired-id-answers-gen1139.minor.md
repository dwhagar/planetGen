**Headline:** An old link to an object that has changed its ID now leads to the object: a read is redirected, a write is refused with the new ID, and a deleted object says it is gone.

### Added
- **Retired IDs answer (GEN.1139, part c).** `api/ids.py` looks a printed ID that names nothing live up in `designations`. A `GET` or `HEAD` of it is a `301` to the same URL under the object's current ID (with `X-Designation-Formerly`), a `PATCH`, `POST` or `DELETE` against it is a `409` that names the current ID and changes nothing, and an object that was deleted answers `410`; an ID nobody held is still a `404`. The pages follow: a page route asked for a retired ID answers a `301` to its own URL under the current ID, or a `410` page. Nothing changes an ID yet except `designation_history.renumber` (called by tests until the edit flow of 8.6), so no link is redirected in a galaxy that has not used it.

No generator output changes.

generation-output: unchanged

TODO items resolved: none (TODO closes GEN.1139 when the branch merges)
