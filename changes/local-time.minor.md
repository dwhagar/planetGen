### Changed
- **Times show in the viewer's own time zone.** Pages write every time
  as UTC in a `<time>` element (labelled "UTC", so they read correctly
  without script), and the new `static/localtime.js` rewrites each in the
  browser's zone with its abbreviation: API key created/last used/revoked
  times, the Stats page's activity times, and a Generate job's start
  time. The database connection's session zone is now pinned to UTC, so
  `TIMESTAMP` columns read back the same whatever the server's own zone
  is, and the API's key times and the stats times end in `Z`.
