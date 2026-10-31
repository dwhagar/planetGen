**Headline:** An ID an object gave up still works: the API sends a read to the object's current ID, refuses a write with the current ID named, and says gone for a deleted object.

### Added
- **Answers for retired IDs (GEN.1139, first slice: the API).** When a printed ID names nothing now but an object once had it (the `designations` table of GEN.1138), `GET` and `HEAD` answer `301` with `Location` set to the same URL under the current ID (path and query string kept) and `X-Designation-Formerly`; `PATCH`, `POST` and `DELETE` answer `409` and change nothing; the retired ID of a deleted object answers `410` to both. The JSON body names the retired ID, the current one, and the object's `former_designations`. An ID nobody had, or one retired from an object of another kind, is still `404`. `designations.status_of` is the rule, `api/ids.retired_answer` the lookup (one indexed read, only when the live lookup missed), and `ApiError` can now carry response headers.
- **Documented** in the Errors section of `docs/reference/api.md`. Additive: the API version does not change.

### Changed
- Nothing else. The page routes (a 301 to the current page URL and an error page for a deleted object) are the next slice of GEN.1139, so the item stays open.

No generator output changes.

generation-output: unchanged

TODO items resolved: none (GEN.1139 stays open for the page routes).
