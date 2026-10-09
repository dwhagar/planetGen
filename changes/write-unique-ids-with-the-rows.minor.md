### Changed
- A saved sector's rows now get their unique ID in the INSERT instead of being selected back and updated afterwards, which saves a query pass per sector (PERF.44). The IDs are the same.
