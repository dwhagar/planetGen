**Headline:** The admin Update button: the Owner presses Update on a page and update.sh runs as root, with no root shell.

### Added
- **SEC.43: `/admin/update`, and `GET`/`POST /api/admin/update`.** The page shows whether an update is waiting or running and the end of `update.log`, and reloads itself every ten seconds while one runs. The Owner types the password again (the same ten-minute window as the Owner's settings) and presses Update. The web process only creates the empty `update.request` file in the state folder; the root-owned `planetgen-update.path` does the rest (see the units note), so nothing the browser sends reaches the command. A press while an update runs runs after it. Not set up (no unit) the page says so.
- The Update page is in the admin rail's search and the Overview list. The API route is documented in `docs/reference/api.md`; the route snapshot moved (additive, no version bump).
