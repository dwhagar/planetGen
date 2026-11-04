**Headline:** An ID an object gave up still works: the API sends a read to the object's current ID, refuses a write with the current ID named, and says gone for a deleted object.

### Added
- **Answers for retired IDs (GEN.1139).** When a printed ID names nothing now but an object once had it (the `designations` table of GEN.1138), `GET` and `HEAD` answer `301` with `Location` set to the same URL under the current ID (path and query string kept) and `X-Designation-Formerly`; `PATCH`, `POST` and `DELETE` answer `409` and change nothing; the retired ID of a deleted object answers `410` to both. The JSON body names the retired ID, the current one, and the object's `former_designations`. An ID nobody had, or one retired from an object of another kind, is still `404`. `designations.status_of` is the rule, `api/ids.retired_answer` the lookup (one indexed read, only when the live lookup missed), and `ApiError` can now carry response headers. A page asked for under a retired ID is a 301 to the same page under the current one, and the page of a deleted object is a 410 (`apiclient.MovedError` and `GoneError`, both kinds of `NotFoundError`; the HTTP transport no longer follows a redirect by itself).
- **Documented** in the Errors section of `docs/reference/api.md`. Additive: the API version does not change.

### Changed
- Nothing else.

No generator output changes.

generation-output: unchanged

TODO items resolved: GEN.1139.
