**Headline:** A data table whose rows are slow to build shows a progress bar instead of stalling.

### Changed
- **The tables' scroll and sort wait on a page job (UX.136).** When a table's data route (`/table/<name>`) would take long, it answers `202` with the job (`pending`, `job`, `percent`, `status_url`) as the pages do; `static/datatable.js` shows the same bar and estimate above the rows (the rows already loaded stay), polls `/api/page-jobs/<id>`, and asks again when the job has succeeded. A failed or cancelled job shows why with Try again (which replaces it); Stop waiting ends the polling and Check again picks the job up. Browser tests cover a slow scroll and a slow sort, in both themes and at phone width.
