### Changed
- A saved sector's rows now get their unique ID in the INSERT instead of being selected back and updated afterwards, which saves a query pass per sector (PERF.44). The IDs are the same.

### Fixed
- The admin page's creation-settings download answers an anonymous or non-admin caller with a plain 403, like the other file and data views, instead of a redirect (ADM.18).
