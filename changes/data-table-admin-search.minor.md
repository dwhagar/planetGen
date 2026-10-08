### Added

- The admin API keys list, the duplicate-names list on the stats page, the
  job queue's Jobs list and every Search result panel are data tables like
  the others (UX.41): they scroll through every row and keep the 50-row pages
  for scripts-off visitors. API keys sort by label, dates and status and
  filter by status; Revoke keeps working from the scrolled table.
- `GET /api/search` takes `panels` (comma-separated panel names) to run just
  those result panels.
- Rows come from `/table/api-keys`, `/table/duplicate-names`,
  `/table/queue-jobs` (admins only) and `/table/search-<panel>`.

### Changed

- Times in the API keys, duplicate-names and Jobs tables read in UTC (they no
  longer switch to the viewer's time zone).
- Revoking an API key returns to the API keys list, not to a numbered page.
