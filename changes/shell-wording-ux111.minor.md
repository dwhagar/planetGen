### Added
- Maintenance mode: while it is on, every page and API call answers "being updated" (a page at `/maintenance`, a JSON 503 from `/api`) to anyone but a signed-in admin or an admin's API key; sign-in and static files stay reachable. It is one file (`maintenance.json` beside `settings.json`), switched from the admin page or `python3 -m planetgen.cli.maintenance on|off|status`. The site is open in its normal state: `planetgen.cli.migrate` turns it on before it converts a database to the spherical sector grid, and the map rebuild turns it off with `maintenance.turn_off_for(maintenance.SECTOR_SWITCH)`.

### Changed
- The API, the web pages and the Generate page speak the spherical sector grid: sectors are addressed by shell, band and wedge (`/api/galaxy/cell`, `/api/galaxy/bright-stars` and `/api/galaxy/uncharted` take `?shell=&band=&wedge=`, `/api/galaxy/shape` returns `shells`), the Generate page's modes are random start, a whole shell, around a sector, one address and a span of shells, bands or wedges, and the sector page and nav chart show the new address.
- The Galaxy Map redirects to the maintenance page while maintenance mode is on, and shows a notice in place of the map when it is off but the new map has not landed; the Sector Map shows the same notice.

### Removed
- The Generate page's whole-ring, column and Galaxy Map block modes, and the radius-in-sectors cylinder fields.

TODO items resolved: UX.111, UX.112 / Added: none
