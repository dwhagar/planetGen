### Changed
- **The sector page moved to `/sector/<id>`.** It is a Flask page now, with
  the new header and breadcrumbs, bookmarkable Contents pages
  (`?contents_page=N`) and no database in the URL. The Sector Map's info
  panel buttons and its no-JavaScript list are plain links. The admin
  forms (wiki upload, generate neighborhood) carry a CSRF token and, once
  they succeed, redirect back to the page with a message, so reloading
  never repeats them. `sector.py` answers 301 to the new address.
- **The NAV page moved to `/nav`, and every step is a GET URL.** Endpoints
  read as `<kind>:<id>`: `/nav?from=system:12&to=nebula:3`; the pickers
  use `from_sector`/`to_sector`. The NAV Map's points and the route's
  stops are plain links, and a "Reverse course" link swaps the endpoints.
  `nav.py` answers 301 to the new address, translating its old
  parameters, and the old `from_id`/`from_kind`/`from_type` style
  redirects to the new one.
