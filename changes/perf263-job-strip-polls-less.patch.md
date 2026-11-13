**Headline:** A logged-in admin page no longer asks the server who is logged in every 3 seconds: the job strip asks once a minute when nothing runs or waits.

### Changed
- **The admin job strip polls by need (PERF.263).** Each ask of `/admin/generate/status` looks the admin up (`/api/auth/me`), which filled the warm-up log with a line every 3 seconds on any open admin page. The strip now asks every 3 seconds only while a job runs or waits, once a minute when none does, at once when the tab is shown or focused, and not again after a 401 or 403.

generation-output: unchanged
generation-affects: none
