### Added
- **The home page is now served by the Flask app, with a new header.**
  `/` shows every sector and standalone system (each table paged on its
  own with `?sectors_page=N`/`?standalone_page=N`), and `/sectors` and
  `/systems` show one table each. These are plain, bookmarkable GET URLs
  with no database name in them (the database comes from `config.json`'s
  `mysql.database`). The new header has the sections (Galaxy, Sectors,
  Systems, Phenomena, Nav) with the current one marked, a search box,
  Login or Admin/Stats/Logout and the theme button. On phones these fold
  into a Menu that works without JavaScript. There is also a "Skip to
  content" link and breadcrumbs. The pages are Jinja2 templates in
  `src/html/web/` (autoescaped), served without a process start or an
  HTTP call back to the API. The other pages still run as CGI and move
  over in later releases; see `docs/html-interface.md`, "Flask pages".
- `secret_key` in `config.json` (or `PLANETGEN_SECRET_KEY`), used to sign
  CSRF tokens for the Flask pages' forms.

### Changed
- `index.py` and `browse.py` now answer 301 to `/`, keeping their page
  numbers, so old links and bookmarks still work.
- The example Apache vhost mounts the Flask app at `/` instead of `/api`,
  serves `/static/` with `Alias` and runs the remaining CGI pages through
  `ScriptAliasMatch`. **Existing servers need their vhost updated**: see
  `docs/apache-deployment.md`, "Updating an existing server".
- The API's app-wide default rate limit no longer counts the calls pages
  make in-process. Login and write limits still apply.
- Unknown URLs outside `/api` get an HTML 404 page instead of JSON.
