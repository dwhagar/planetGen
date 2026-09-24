# planetGen Web Interface

A small web interface (in [`../src/html/`](../src/html/)) for
browsing the MySQL databases described in
[`database-schema.md`](database-schema.md) -- pick a database (schema),
drill into its sectors and star systems, and view (or copy) the
rendered wikitext/Markdown page saved for each one. Sectors and systems
also get their own interactive visualizations: a draggable 3D "Sector
Map" per sector (translucent clouds for any nearby nebula/asteroid field,
glowing points for any nearby black hole/neutron star),
a galaxy-scale "Galaxy Map" (Quadrant/Ring drill-down) for sectors
actually placed in the galaxy (small dots for every galaxy-placed
phenomenon -- nebula, asteroid field, black hole, or neutron star), and a
true-position "System Map" per system, plotting
every star/planet/moon/belt at its actual angle and a shared log-scaled
distance from whatever it orbits.

This is a read-only browser, and the frontend half of `TODO.md`'s Phase 5
web application: every page here is a thin server-rendered client of the
Flask API in [`../src/html/api/`](api.md), never touching the database
directly itself, though the search page (`/search`, see "Flask pages"
below) does provide a faceted/name search (built on `GET /api/search`, same as every
other page). It exists so a generated galaxy can be looked at from a
browser today, on nothing more than Apache2, a system Python 3, and the
API process (see "Locating the database" below).

## How it works

Every page here is a plain [CGI](https://en.wikipedia.org/wiki/Common_Gateway_Interface)
script (`#!/usr/bin/env python3`, executed fresh by Apache on every
request) that fetches its data from `GET /api/...` (`../src/html/lib/apiclient.py`,
stdlib `urllib` only -- no Flask/FastAPI, no database driver, nothing
beyond the standard library) and renders it as HTML. Every navigational
link between these pages (`db`, a sector/system/phenomenon id, a search
filter, and the like) is a same-origin `POST` form instead of a plain
`<a href="page.py?...">` (`../src/html/lib/fmt.py`'s `post_link`, styled
to look exactly like a normal link -- `.link-btn` in `static/style.css`),
so none of it shows up in the browser's own address bar the way a GET
link's target URL always would; `lib/page.py`'s `nav_params`/
`nav_multi_params` are what a page reads these back with. This makes
every page un-bookmarkable and un-shareable by URL (a deliberate
trade-off -- see `lib/page.py`'s own module docstring) and, for the small
number of markers that still need a real, clickable `<a>` because a
`<form>` can't nest inside an SVG shape (the Galaxy Map's/Sector Map's/
NAV Map's own plotted dots), requires JavaScript (`static/navform.js`) --
every one of those markers also has a plain, no-JS-required row in a
table below its own map, so a marker click is never the only way to reach
something. `nltk` (from
generation) is still needed for a handful of display-formatting constants
a few pages import from `stellarObjects` (spectral-class colors, planet-
class descriptions, unit formatters) -- importing any `stellarObjects`
submodule runs that package's own `__init__.py`, which pulls in `names.py`
and its NLTK corpus dependency, even though nothing under `html/` talks
to MySQL anymore. That keeps `html/`'s own deployment to "copy the files,
set permissions, enable a vhost, point it at a running API" -- see
[`api.md`](api.md#deploying-behind-apache-mod_wsgi) for the API process
itself, which is what actually needs `pymysql`/`DBUtils` and a database
account now.

| File | Purpose |
|---|---|
| `../src/html/index.py` | Shim: 301 to `/`, the Flask-served home page (see "Flask pages" below). |
| `../src/html/browse.py` | Shim: 301 to `/` (keeping `sectors_page`/`standalone_page`); the sector and standalone-system lists are Flask pages now (`/`, `/sectors`, `/systems`, see "Flask pages" below). |
| `../src/html/galaxy.py` | "Galaxy Map": a real perspective-camera WebGL scene (`lib/galaxymap3d.py` + `static/galaxymap3d.js`, three.js) a visitor can rotate/dolly/click through, rather than only ever viewing the galaxy from directly above a flat projection (a former SVG-based version had exactly that limitation, and its markers grew relative to the view as you zoomed in with no camera to shrink them the opposite way -- replaced, not kept alongside). Only fetches the zoomed-all-the-way-out starting view itself (`GET /api/galaxy/shape`, then that view's cube tiles through `lib/tilecache.py`); every later view, as the camera moves, is fetched by the page's own client-side JS directly from `galaxy_tiles.py`, never through this handler again. Left-click zooms in on whatever was clicked (a placed sector, a real not-yet-generated address, or empty space), right-click zooms out -- both by a *logarithmic* step (big jumps zoomed out over the whole galaxy, fine ones once close to a single sector), not a flat factor. Clicking a not-yet-generated address shows its designation and a copyable `generate.py galaxy --shell K --slot N` command (see that flag's own docstring in `generate.py`) instead of navigating anywhere. Below the map, still keeps its own two data tables (`GET /api/galaxy/sectors`), independent of the map itself: every sector actually placed in the galaxy grouped into four azimuthal Quadrants (I-IV) and concentric Rings (`shell_index` bands, `lib/galaxymap.py`) -- the full-galaxy default view is a per-Quadrant summary table; a `quadrant` of `I\|II\|III\|IV` swaps it for a full, distance-sorted sector list for that Quadrant (no longer cropping the map itself -- the 3D map has its own free-flying camera). Sectors never placed in the galaxy (the common case for a plain `sectorGen.py` run) stay out of this page, visible only via `browse.py`'s flat sector table. |
| `../src/html/galaxy_tiles.py` | The 3D Galaxy Map's tile endpoint, fetched by `static/galaxymap3d.js` as its camera moves: returns the requested fixed cubes of space (`tiles=level/ix/iy/iz,...`) and optionally one density cloud (`density=level/ix/iy/iz`), JSON via `page.run_json`. Goes through `lib/tilecache.py`, so only tiles not already cached on disk reach the API (`GET /api/galaxy/tiles`). A malformed request is a 400. Not meant to be navigated to directly. |
| `../src/html/sector.py` | One sector's contents -- its systems and nearby phenomena in one table, nearest the sector's center first -- plus an interactive 3D "Sector Map" of the sector's cube -- drag to rotate, scroll (or the +/- buttons) to zoom, one dot per placed system (two, overlapping, for a binary) sized by star radius and colored by spectral type/brightness, plus a translucent spheroid cloud for every nebula/asteroid field/supernova remnant (or a small glowing point for every black hole/neutron star/rogue planet/interstellar comet) whose real galaxy-frame sphere reaches into this sector's cube, or that was generated as part of it (`queryDb.phenomena_near_sector`, part of `GET /api/sectors/<id>`); clicking a dot, cloud, or point fills an info panel with a link to `system.py` for a star system, via `lib/starmap.py` + `static/sectormap.js`. Shows a "View on Wiki" link once `sectors.wiki_url` is set; an admin session additionally gets an "Upload to Wiki" form (`POST /api/sectors/<id>/wiki`) while it isn't. |
| `../src/html/system.py` | One system rendered natively as an expandable "System" list: the overview (a binary pair's data, the system summary), then its stars, planets (moons nested under each; a 'wide' pair's bodies nested under their own star), asteroid belts and comets. Each row shows compact stats (class, terrestrial/gas giant, Habitable yes/no, Inhabited yes/no, moon count) and opens onto that body's own generated description (`GET /api/systems/<id>/sections`). Wikitext and Markdown buttons (`code=wikitext|markdown`) show the full generated page (`GET /api/systems/<id>/text`) in a code box with a Copy button (`static/copycode.js`) -- nothing is stored, both are rendered from the database on demand. Links to the wiki page(s) once uploaded (`wikijs_url`/`mediawiki_url`, opening in a new tab). Also shows the stars/bodies tables. Also renders the interactive "System Map" panel (`lib/systemmap.py` + `static/systemmap.js`) whenever the system has stars -- a true top-down plot of every body's *real* position (real angle from `position_x/y_km`, log-scaled distance from whatever it orbits) with one clickable marker per body, colored by planet class and sized by `radius_km`; a binary pair's two stars (and, for a 'wide' pair, each star's own independent planets) are placed at their real mass-weighted offsets from the system's shared barycenter (`binary_mutual_position_x/y/z_km`). Clicking a planet with moons zooms into that planet's own moon system, and any body with `life_chemical` set gets a small green "supports life" badge (surfaced in the info panel's "Life Chemistry" field too). Links to `nav.py` (a "Navigate from here" button) whenever the system is assigned to a sector. An admin session gets an "Upload to Wiki" form (`POST /api/systems/<id>/wiki`) offering whichever backend(s) are configured and not yet uploaded to. |
| `../src/html/nav.py` | Course, distance, and optimal route between two systems -- calls `GET /api/nav` (`queryDb.nav_between`; see [`api.md`](api.md#nav) for the availability rules and course convention). Reachable either via `system.py`'s "Navigate from here" button or the sidenav's own "Nav" link with no origin known yet, which shows a two-step sector-then-system picker (`GET /api/sectors` then that sector's own systems) to choose one. Without a destination yet, shows a dropdown of the origin's own sector-mates plus, when that sector has a galaxy placement, the same two-step sector-then-system picker (scoped to `to_sector`) for a cross-sector destination (there's no bounded way to offer every galaxy-placed system in one dropdown). With a destination, shows the direct course (distance/azimuth/altitude/warp-1-3-6-9 travel times), the NAV Map (a flat, top-down SVG plot of the origin/destination/route hops -- `html/lib/navmap.py`), and the optimal route via adjacent systems, each hop linking to `system.py`. |
| `../src/html/phenomena.py` | Flat, paginated list of every exotic phenomenon (nebula/asteroid field/black hole/neutron star/supernova remnant/rogue planet/interstellar comet) across every sector, regardless of galaxy placement -- `GET /api/phenomena`. Each row links to `phenomenon.py`. |
| `../src/html/phenomenon.py` | One phenomenon's detail/info page -- `GET /api/phenomena/<type>/<id>`, this project's first per-phenomenon page (previously a phenomenon had no page of its own, only a hover tooltip on the Sector Map). Reached from `phenomena.py`'s listing, or directly from a Sector Map marker (`lib/starmap.py`) -- the Galaxy Map (`galaxy.py`) doesn't plot phenomena of its own. |
| `../src/html/search.py` | Shim: 301 to `/search` (the Flask-served search page, see "Flask pages" below), carrying every search parameter (name fields, size ranges, repeated tag facets, `<panel>_page`) from the old GET query or POST body and dropping `db`. |
| `../src/html/lib/pagination.py` | The site's one pager, used under every paged table (Browse's two tables, Phenomena, a sector's Contents table, a Galaxy Map Quadrant's sector list, each Search result panel, the admin API key list and the admin stats page's duplicate-names list): a "Showing X-Y of Z" summary, then First/Prev, numbered pages and Next/Last, 50 rows a page. Each table has its own page parameter (e.g. `sectors_page`), posted like every other link here, and changing a Search filter starts its results back at page 1. |
| `../src/html/login.py` | Admin login form -- `POST /api/auth/login`, relaying the session cookie it sets back to the browser. Redirects to `changecreds.py` (still on the seeded `admin`/`password` default) or `admin.py` on success; re-renders the form with an inline error for a wrong username/password. See [`api.md`](api.md#authentication). |
| `../src/html/changecreds.py` | Change the logged-in admin's username/password (`POST /api/auth/change-credentials`) -- reached both by `login.py`'s forced redirect (default credentials) and voluntarily (an already-"fresh" admin rotating credentials). Always requires the current password. |
| `../src/html/admin.py` | Protected admin landing page: API key list/create/revoke (`GET`/`POST /api/auth/api-keys`, `DELETE /api/auth/api-keys/<id>`) -- keys are used as `Authorization: Bearer` credentials against the write endpoints (sector/system create/update/delete), not exercised as a web form here. The one exception is a small form to manually set/clear a sector's `wiki_url` (`PATCH /api/sectors/<id>`, by database name + sector id, via the session cookie). Redirects to `login.py`/`changecreds.py` if not authenticated / still on default credentials. |
| `../src/html/adminstats.py` | Admin-only server health and database stats (`GET /api/admin/stats`, `GET /api/admin/duplicate-names`): API/MySQL health, galaxy tile cache usage, exact sector/system counts, schema version, last created/modified times, per-table sizes, and every name made unique (Alpha/Beta..., Little..., ...Kin) with links to each sector and system (planets and moons link to their system). Same login gate as `admin.py`; linked from it and from the sidenav's Stats item. |
| `../src/html/logout.py` | Ends the current admin session (`POST /api/auth/logout`, called server-side) and redirects to `index.py`. A plain sidenav link, not a form -- nothing to confirm, no request body needed. |
| `../src/html/lib/apiclient.py` | The `GET /api/...` HTTP client (stdlib `urllib` only) every page above calls instead of querying MySQL directly -- one typed wrapper function per read endpoint, plus `auth_*` wrappers (`login.py`/`changecreds.py`/`admin.py`) supporting POST/DELETE, a request body, and `Cookie`/`Set-Cookie` relay, and `NotFoundError`/`ApiError` (`lib/page.py`'s `run` turns these into a 404/502 page; `ApiError.status_code` lets `auth_me`/the admin pages branch on a 401 without string-matching). `PLANETGEN_API_BASE_URL` (default `http://127.0.0.1/api`) is where it looks for the API. Not web-accessible. |
| `../src/html/lib/fmt.py` | HTML-escaping and small formatting helpers (`esc`, `linkify_location`, `format_density`, and `static_url`, the versioned `static/` URL every page uses) with nothing to do with fetching data -- what's left of the old `dbutil.py` once its database-access functions moved into `apiclient.py`/the API itself. Also `post_link` (a same-effect, no-JS-required replacement for `<a href="page.py?...">` that posts its params as hidden fields instead -- see "How it works" above) and `data_nav_params` (its JS-required counterpart for a marker embedded in an SVG map, since a `<form>` can't nest inside one -- paired with `static/navform.js`). Not web-accessible. |
| `../src/html/lib/page.py` | Shared CGI response/HTML-shell helpers -- `query_params`/`form_params`/`form_multi_params` (GET/POST parsing) and `nav_params`/`nav_multi_params` (POST body when present, else the GET query string -- what every page reads a `post_link`-followed link's params back with), `incoming_cookie_header` (relays the browser's own `Cookie` header to the API, unparsed), `send_headers`/`render`/`redirect` (all three accept `set_cookie_headers` to relay the API's own `Set-Cookie` back), the sidenav's Search/Galaxy/Sectors/Systems/Nav/Phenomena/Login/Admin/Logout items (via `apiclient.auth_me` and `fmt.post_link`) plus the theme button, `head_html` (the shared `<head>`), and `SECURITY_HEADERS` (the one source of the pages' CSP and other security headers). See "The page shell" below. Not web-accessible. |
| `../src/html/lib/mdconvert.py` | A small, purpose-built Markdown-to-HTML converter for the narrow Markdown subset `StarSystem.__str__` actually generates (headers, pipe tables, paragraphs, `<sup>` exponents) -- not a general-purpose parser. Not web-accessible. |
| `../src/html/lib/galaxymap.py` | Quadrant/Ring classification (`sector_quadrant`/`sector_ring`/`ring_bounds_ly`) shared by `galaxy.py`'s own data tables and `sector.py`/`browse.py`'s "Quadrant N" links back into it -- what's left after this module's former flat-SVG map rendering was superseded by a real 3D scene (`lib/galaxymap3d.py`). Not web-accessible. |
| `../src/html/lib/tilecache.py` | The web layer's on-disk cache of 3D Galaxy Map tiles: `fetch_tiles` serves each requested tile from `<cache dir>/<db>/<generation>/` when present and asks `GET /api/galaxy/tiles` only for the rest. Every 60 s at most it asks `GET /api/galaxy/changes` which tiles changed (from the rows' `modified_at`) and deletes only those, or starts a new generation when the answer is `full`. The changed keys are passed on to the browser's cache too. Size-capped (`tile_cache.max_mb`, oldest files pruned first); location from `PLANETGEN_TILE_CACHE_DIR`/`tile_cache.dir` (see `config.md`). Fails open: an unwritable directory or bad file just means an API call. Not web-accessible. |
| `../src/html/lib/galaxymap3d.py` | Builds `galaxy.py`'s "Galaxy Map" panel: the canvas/controls/info-panel markup, plus the one starting JSON payload (zoom-range numbers, tile settings, and the zoomed-all-the-way-out view's cube tiles, which `initial_tile_request` picks the same way the client does for every later view) `static/galaxymap3d.js` reads on first paint -- every later payload, as the camera moves, is fetched by that script directly and never passes through this module. `view_radius_bounds` (the zoom floor/ceiling: a couple of sector-widths up to this galaxy's own real outer edge) is a pure function `galaxy.py` also calls directly, before this module's own panel-rendering function runs. Not web-accessible. |
| `../src/html/lib/systemmap.py` | Builds `system.py`'s "System Map" panel: a fixed-size, square true-position plot (real angle from `position_x/y_km`, a shared log-radial scale from anchor -- star, or barycenter for a merged/'close' binary pair -- to body) with one colored/sized marker per body, plus a full ring (not a directional band) for each asteroid belt -- distinct from `starmap.py`'s draggable 3D cube, since this only ever needs a flat top-down projection. A once-only pairwise-repulsion pass (`_relax_markers`) nudges apart any two markers real placement happened to put too close together, real position first, decluttering only where needed. Clicking a planet with moons swaps to a "zoom into this planet's moons" scene (`static/systemmap.js` toggles which `<svg>` scene is visible); clicking anything else fills the info side panel, same click-for-info pattern as `starmap.py`/`sectormap.js`. Adds a small green badge to any body with `life_chemical` set. Also hands over each planet/moon marker's resolved class color (`_class_color`) plus its `atmosphere`/`composition`/`surface_temperature_k` (whether it has a real atmosphere at all, not just the description text) as `data-*`, plus a star's own spectral-type color (`_star_color`, same `data-color` attribute), for `static/systemmap.js`'s own per-marker live sphere rendering. Not web-accessible. |
| `../src/html/lib/tabledisplay.py` | Computes the same "Star Data"/"Planet Data" display strings once baked into the database's now-removed `table_*`/`binary_table_*` columns, but on demand from the raw numeric columns `system.py` already has -- reuses `stellarObjects.utils`'s formatters directly. Not web-accessible. |
| `../src/html/lib/starmap.py` | Builds `sector.py`'s "Sector Map" panel: computes every position/size/color/label the map needs (a wedge or fallback-cube outline, one entry per star -- billboarded, radius from `radius_km` square-root scaled against the Sun, color from `star_type`'s spectral letter (`SPECTRAL_CLASS_COLORS`) shaded by `luminosity_w` and nudged by where `temperature_k` falls in that spectral class's range, so "White Giant" reads white and "Blue Giant" reads blue regardless of temperature -- and one per nearby standalone phenomenon, sized by `radius_ly` (always 0 for the point-like types) and positioned directly in the galaxy frame, no rotation needed unlike a star system's sector-local position -- see `queryDb.phenomena_near_sector`) and serializes it as a `<script type="application/json">` block; `static/sectormap.js` is what actually renders it, this module builds no HTML scene of its own. Not web-accessible. |
| `../src/html/lib/navmap.py` | Builds `nav.py`'s "NAV Map" panel: a flat, static, top-down SVG plot of the galactic X-Y plane -- origin and destination as labeled points, a dashed line for the direct course, and (when one was found) a solid polyline through the optimal route's intermediate hops. Auto-scaled to whatever points it's given (no fixed sector size to normalize against), with one uniform light-years-per-pixel ratio on both axes so azimuth angles aren't visually distorted, plus a `+X` compass tick and a scale-bar legend. Deliberately blind to altitude/z, same as the flat SVG phenomenon Diagram panel (`lib/phenomenonmap.py`) -- the course panel's own Altitude figure already covers that axis. Not web-accessible. |
| `../src/html/static/style.css` | Shared stylesheet (CSS custom properties, light/dark via `prefers-color-scheme` or an explicit `data-theme` on `<html>`, card-style panels, phone layout under 40rem), served directly (not through CGI). |
| `../src/html/static/theme.js` | Loaded on every page, blocking, before `style.css` (`lib/page.py`'s `head_html`): applies the saved light/dark/system theme before the first paint and drives the side nav's theme button. See "The page shell" below. |
| `../src/html/static/favicon.svg` | The site icon (a small ringed planet), linked from every page's `<head>`. |
| `../src/html/static/vendor/` | Vendored third-party JS -- currently just `three.module.min.js` (three.js, bundled+minified from the `three` npm package), used by `sectormap.js`/`systemmap.js`. Vendored rather than loaded from a CDN so `lib/page.py`'s `Content-Security-Policy: default-src 'self'` needs no exception; see this directory's own `THIRD_PARTY_NOTICES.txt` for the license and how to rebuild it from a newer release. Served directly, same as `style.css`. |
| `../src/html/static/sectormap.js` | Renders `lib/starmap.py`'s `#starmap-data` JSON as a real WebGL scene (three.js, vendored at `static/vendor/` -- see that directory's `THIRD_PARTY_NOTICES.txt`): a perspective camera, GPU-billboarded sprites (always face the camera by construction, no manual per-frame counter-rotation needed) for stars/nebulae/asteroid fields/black holes/neutron stars, and a wireframe outline for the sector's wedge or fallback cube. Hand-rolled drag-to-rotate/scroll-or-button-to-zoom/click-for-info (a raycast against the sprites, not DOM hit-testing) plus a hidden, keyboard/screen-reader-focusable button list (a canvas has no focusable children of its own the way the old per-star `<div role="button">`s were) so every star/cloud stays reachable without a mouse. Clicking (or activating a fallback-list button) fills the info side panel and, for a "View system"/"View phenomenon" link, navigates via `data-nav-target`/`data-nav-params` (`static/navform.js`, whose global click handler also intercepts this dynamically-built link, same as any other in `html/`) instead of a plain `href`, so `db`/an id never shows up in the address bar. Served directly, same as `style.css`. |
| `../src/html/static/galaxymap3d.js` | Renders `lib/galaxymap3d.py`'s `#galaxymap3d-data` JSON as a real WebGL scene (three.js, vendored at `static/vendor/`) -- unlike `sectormap.js`'s camera (always orbiting a fixed origin, every star baked into one page load), this camera's own orbit target moves freely through the galaxy, so most of what it draws is fetched live from `galaxy_tiles.py` (debounced, on every camera move), one fixed cube of space at a time, rather than server-rendered once. Tiles are kept in memory and in the browser's `localStorage`, keyed by the database's content stamp, so panning back over seen space or reloading the page doesn't refetch them. Drag to rotate, scroll to zoom; left-click zooms in on whatever's under the cursor (a placed sector, a real address, or empty space -- via a raycast against a depth sphere at the camera's own current distance), right-click zooms out -- both by a logarithmic step (`clickZoomFactor`, interpolated between `lib/galaxymap3d.py`'s own `clickZoomFactorMin`/`Max` by the camera's current distance in log space) rather than a flat factor, so a handful of clicks still crosses the whole galaxy while a click near one sector stays fine enough not to overshoot it. Placed/planned sectors are diffed against what's already drawn (`syncTier`) rather than rebuilt on every fetch, so already-visible content doesn't flicker as the camera moves; the illustrative density tier is a single `THREE.Points` cloud instead of individual sprites. A real perspective camera means dot size naturally scales with distance for free -- no manual per-zoom counter-scaling needed the way a flat SVG `viewBox` zoom would. Clicking a not-yet-generated address's dot shows its designation/predicted count and a "Copy CLI command" button (`generate.py galaxy --shell K --slot N`, clipboard API with a selectable-field fallback). Served directly, same as `style.css`. |
| `../src/html/static/systemmap.js` | Toggles which System Map `<svg>` scene is visible (the whole-system view, or one per planet's own moon system) and fills the info side panel -- including "Atmosphere"/"Surface composition"/"Surface temperature"/"Life Chemistry" fields -- from a clicked marker's `data-*` attributes. Every visible scene's star/planet/moon markers also get their own live-rendered 3D sphere on `#sysmap-spheres-canvas` (the same vendored three.js build `sectormap.js` uses): one shared WebGL context, redrawn each frame via a scissored sub-viewport per marker (never one `<canvas>`/context per body -- browsers cap concurrent WebGL contexts), each sized and positioned to exactly cover that marker's own `<circle>` and colored by its `data-color`, banded with a tilted ring for a gas giant (`data-bodytype`), and wrapped in a fresnel-glow atmosphere shell (tinted by `data-surfacetemp`) when `data-hasatmosphere` is set -- the one genuinely 3D layer on this otherwise flat-SVG page, drawn behind the SVG so each marker's own stroke/label/life-badge still shows on top. An appearance layer only (no position data, and no cross-page link of its own to navigate) -- a marker whose sphere renders keeps its flat circle's fill transparent (`sysmap-sphere-active`) but never changes its actual plotted position. Served directly, same as `style.css`. |
| `../src/html/static/copycode.js` | The system page's Copy button: copies the generated Wikitext/Markdown out of the code box named by the button's `data-copy-target` (clipboard API, falling back to a selection copy on a plain-HTTP deployment). A separate file because the Content-Security-Policy allows no inline script; the box stays selectable by hand without it. Served directly, same as `style.css`. |
| `../src/html/static/navform.js` | Loaded on every page (`lib/page.py`'s `render`). Lets a `data-nav-target`/`data-nav-params` element (a Galaxy Map/Sector Map/NAV Map marker, or the Sector Map's own dynamically-built info-panel link -- see `lib/fmt.py`'s `data_nav_params`) navigate by posting a throwaway hidden form, since a `<form>` can't nest inside the SVG/canvas those are drawn in -- everywhere else in `html/`, `lib/fmt.py`'s `post_link` builds a real, no-JS-required `<form>` instead. Served directly, same as `style.css`. |

This project recommends (but doesn't enforce in code) pointing the API
(not `html/` itself anymore -- see "How it works" above) at a MySQL
account with `SELECT`-only grants, so it can't write to a database even
if a query were buggy -- see `queryDb.py`'s module docstring for the same
convention, and [`api.md`](api.md) for where that account is configured.
Database and system names pulled from generated data are HTML-escaped
before being placed in a page; a requested `?db=` schema name is
validated API-side against the actual, prefix-filtered schema listing
(exact match only, `stellarObjects._db.resolve_database`), which is what
prevents it from being used to select a schema this deployment never
meant to expose. `mdconvert` escapes every block in full before emitting
any markup, then narrowly re-enables only the one legitimate raw-HTML
pattern generated content ever contains (`<sup>...</sup>`) -- so a
mischievous `--name`/`--star-type` value can't inject live HTML into a
rendered page.

## The page shell

Every page goes through `lib/page.py`'s `render`, which writes the
security headers (`SECURITY_HEADERS`: a `Content-Security-Policy` of
`default-src 'self'; base-uri 'self'; form-action 'self'; frame-ancestors
'none'; object-src 'none'`, plus `X-Frame-Options`, `X-Content-Type-
Options` and `Referrer-Policy`) and then the shared `<head>` from
`head_html`:

- `<meta name="viewport" content="width=device-width, initial-scale=1">`,
  so phones use the narrow layout (`@media (max-width: 40rem)` in
  `style.css`: the side rail becomes a wrapping row of links at the top,
  and wide tables scroll inside their own `.table-scroll` box instead of
  widening the page).
- A `<meta name="description">` (a site-wide default, or `render(...,
  description=...)`), `<meta name="color-scheme">`, and two `<meta
  name="theme-color">` tags, one per OS colour scheme.
- The favicon, `static/favicon.svg`.
- `static/theme.js`, a small blocking script placed before the
  stylesheet: it applies the visitor's saved theme before the first
  paint (see below).
- `static/style.css` and `static/navform.js`.

Every `static/` URL, here and in the pages' own `<script>` tags, comes
from `lib/fmt.py`'s `static_url(name)`, which appends `?v=<package
version>` (read from `src/stellarObjects/_version.py`, which the release
Action stamps). A release therefore changes every static URL, and Apache
can let browsers cache them for a year (see
[`apache-deployment.md`](apache-deployment.md#static-files-compression-and-security-headers)).
The ES modules that import siblings (`sectormap.js`, `systemmap.js`,
`galaxymap3d.js` -> `bodyRendering.js`, `vendor/three.module.min.js`) use
`await import(...)` with their own `?v=` (from `import.meta.url`), so each
module has exactly one URL per page and is never loaded twice.

**Theme switch.** The last item in the side nav is a "Theme: System /
Light / Dark" button. `static/theme.js` reveals it (it renders `hidden`,
so a browser without JavaScript never sees a dead button), cycles
system -> light -> dark, sets `data-theme` on `<html>` (none for
"system", so `style.css` follows `prefers-color-scheme`), points the
theme-color tags at the chosen theme, and remembers the choice in
`localStorage` under `planetgen-theme` (guarded, so private windows still
work for the current page). `style.css` has matching
`:root[data-theme="light"]` and `:root[data-theme="dark"]` blocks; the
dark tokens appear twice (OS dark and explicit dark) and must be kept
identical. The map canvases read their colours when they start, so they
follow a theme change on the next page view.

## Flask pages (the pages are moving off CGI)

The site is moving from one CGI script per page to HTML pages served by
the same Flask app as the API (`../src/html/web/`, registered by
`api/app.py`'s `create_app`). Pages move one at a time; until a page
moves, its CGI script keeps working exactly as described above. Moved so
far:

| URL | Replaces | Shows |
|---|---|---|
| `/` | `index.py`, `browse.py` | Every sector and every standalone system, each table paged on its own (`?sectors_page=N`, `?standalone_page=N`). |
| `/sectors` | `browse.py#sectors` | The sectors table alone. |
| `/systems` | `browse.py#standalone-systems` | The standalone systems table alone. |
| `/search` | `search.py` | Faceted search (see below). |

`index.py` and `browse.py` are now CGI shims that answer `301 Moved
Permanently` to `/`, carrying `sectors_page`/`standalone_page` from the
old GET query or POST body (`lib/page.py`'s `moved_permanently`). The
cleanup PR removes them.

**The search page** (`/search`, `web/searchpage.py` + `templates/search.html`,
data from `GET /api/search` in-process) takes every filter as a GET
parameter, so a search is a bookmarkable URL:

| Parameter | Meaning |
|---|---|
| `q` | The main name search, and what the header search box sends. Searches sector, system, star, planet and moon names at once. With an Object Type tag active it searches only those object types (not sectors and systems). |
| `sector_q`, `system_q`, `star_q`, `planet_q`, `moon_q` | A name search for one kind of object; replaces `q` for that kind. Under "Search by object and size", each with a `<datalist>` of existing names. |
| `star_min_radius_km`, `star_max_radius_km` (and `planet_`, `moon_`) | A size range; a bound that isn't a number is ignored. |
| `type`, `spectral`, `luminosity`, `class`, `body`, `life`, `moon_class`, `moon_body`, `moon_life`, `density` | Tag facets, repeated once per selected value (`spectral=G&spectral=K`). |
| `sectors_page`, `systems_page`, `stars_page`, `planets_page`, `moons_page`, `belts_page` | Each result panel's page; each pager keeps the search and the other panels' pages. |

Tags, "remove filter" chips and "Clear all" are plain links to the same
search with one thing changed. Results come right under the form, the tag
browser below them. The search form sends every field, filled or not, so
a URL with empty (or unknown) parameters is redirected (302) to the same
search without them. `search.py` is a CGI shim that 301s to `/search`,
carrying every one of these parameters from the old GET query or POST
body (repeated facets included) and dropping `db`.

What changes for a visitor:

- **Plain GET links and bookmarkable URLs.** Links between moved pages
  are ordinary `<a href>`s (no hidden POST forms), so Back, reload,
  open-in-new-tab, bookmarks and sharing all work.
- **No database in the URL.** The pages show one database, taken from
  config (`mysql.database` / `PLANETGEN_MYSQL_DATABASE`), never from the
  request. Links from a moved page to a page still on CGI are plain GET
  links such as `/sector.py?db=planetgen&id=5` (the CGI pages read a GET
  query through `nav_params()`); that is the one place `db` still shows,
  until each page moves.
- **New header** instead of the side rail: site name, the sections
  (Galaxy, Sectors, Systems, Phenomena, Nav) with `aria-current="page"`
  on the current one, a search box, and Login or Admin/Stats/Logout plus
  the theme button. Below 56rem the sections, search and account links
  fold into a native `<details>` "Menu" (works without JavaScript). A
  "Skip to content" link comes first, and every page but the home page
  has a breadcrumb trail.
- **Faster**: no Python process start per page, and no HTTP call from
  the page back to the API (see below). Visitors without a session
  cookie cost no login lookup at all.

### How a Flask page is built

```
src/html/web/
  __init__.py     blueprint `web`, template globals, init_app()
  views.py        the routes
  helpers.py      db_name, page_url, crumb, render_page, trusted_html,
                  pager, current_admin, LEGACY_PAGES
  transport.py    in-process transport for lib/apiclient.py
  csrf.py         CSRF tokens for POST forms
  errors.py       HTML 404/502/500 pages
  templates/      base.html + one template per page (+ partials/)
```

A page is a route plus a template:

```python
# web/views.py
@bp.route("/phenomena")
def phenomena():
    envelope, page = fetch_page(
        lambda limit, offset: apiclient.get_phenomena(db_name(), limit=limit, offset=offset),
        parse_page(request.args.get("page")),
    )
    return render_page(
        "phenomena.html", title="Phenomena", section="phenomena",
        breadcrumbs=[crumb("Phenomena")],
        rows=envelope["items"],
        pager=pager("page", page, envelope["total"], anchor="list", label="Phenomenon pages"),
    )
```

```jinja
{# web/templates/phenomena.html #}
{% extends "base.html" %}
{% block content %}
<section class="panel" id="list">
  ... {{ row.name }} ...   {# escaped automatically #}
  {{ pager }}              {# Markup from lib/pagination.py #}
</section>
{% endblock %}
```

The helpers (all in `web/helpers.py`):

- `render_page(template, title=, section=None, breadcrumbs=(),
  description=None, status=200, **context)`: renders a template that
  extends `base.html`. `section` is one of `SECTIONS` (`galaxy`,
  `sectors`, `systems`, `phenomena`, `nav`) and gets `aria-current`.
- `crumb(label, name=None, **params)`: one breadcrumb. "Home" is added
  in front automatically; the last crumb (no `name`) is the current page.
- `page_url(name, **params)`: the URL of any page by endpoint name,
  moved or not. Moved pages go through `url_for("web.<name>")`; pages
  still on CGI through `LEGACY_PAGES`.
- `db_name()`: the one database. Never read `db` from the request.
- `trusted_html(html)`: passes HTML built by a `lib/` renderer (the star
  and system maps, `fmt.format_density`, ...) through unescaped. Those
  renderers escape their own inputs; never wrap request or database text
  in it directly.
- `pager(page_param, page, total, anchor=None, label="Pages", keep=None)`:
  the shared pager (`lib/pagination.render_pagination`, GET mode) linking
  back to the current path with `?<page_param>=N`; `keep` carries other
  query parameters (another table's page number).
- `current_admin()`: the logged-in admin or `None`, one lookup per
  request (cached on `g`). Also available in templates.
- In templates: `static_url("x.js")` (`/static/x.js?v=<version>`, from
  `fmt.static_url`), `page_url(...)`, `current_admin()`,
  `csrf_field()`, `site_name`, and blocks `head` (extra `<script>`/
  `<link>`, e.g. a map module), `heading`, `subhead` and `content`.

**Moving a page** (what each page PR does): add its route to `views.py`
under the endpoint name its `LEGACY_PAGES` entry uses (`sector`,
`system`, `galaxy`, ...), add its template, delete the `LEGACY_PAGES`
entry (every link to the page then switches to the new URL; a test fails
if a name is both a route and a legacy entry), and turn the old `.py`
script into a shim calling `page.moved_permanently(new_path, keep)`.

**Data in-process.** Views call the same `lib/apiclient.py` functions as
the CGI pages (`get_sector(db, id)`, `auth_me(cookie_header)`, ...).
Inside a Flask request, `web/transport.py` dispatches each call straight
through the app's own `/api` routes (same validation, auth and JSON as
over HTTP, in a fresh app context so database connections open and
close per call) instead of making an HTTP request to itself. Outside a
Flask request (a CGI page, a script) the functions use HTTP as before.
These in-process calls are exempt from the API's app-wide default rate
limit (a page view is not API abuse); a route's own limit (login,
writes) still applies to the visitor's address. Page views themselves
are not rate-limited, as the CGI pages never were.

**Forms.** Any POST/PUT/PATCH/DELETE to a path outside `/api` must carry
the CSRF token: put `{{ csrf_field() }}` inside the `<form>`. The token
is an HMAC (keyed with `secret_key` from `config.json`, see
[`config.md`](config.md)) of a random value in the `pg_csrf` cookie
(HttpOnly, SameSite=Strict); a missing or wrong token gets a 400 page
and the view never runs.

**Errors.** An `apiclient.NotFoundError` becomes a 404 page, an
`apiclient.ApiError` a 502 page (without the API's detail), anything
else a 500 page saying only "An unexpected error occurred."; the
traceback goes to Apache's error log and the debug log, never the page.
Unknown URLs get the HTML 404 page; `/api/...` keeps its JSON errors.

**Headers.** Every HTML response carries `lib/page.py`'s
`SECURITY_HEADERS` (the same CSP as the CGI pages), set in
`api/app.py`; JSON keeps `default-src 'none'`.

## Locating the database (and the API)

Every page here needs the planetGen API (`../src/html/api/`, see
[`api.md`](api.md)) reachable to work at all now -- set
`PLANETGEN_API_BASE_URL` (or `config.json`'s `api_base_url`; default
`http://127.0.0.1/api`, i.e. the same host this CGI script itself runs on)
if it's deployed somewhere else, e.g. `http://127.0.0.1:5000/api` for
`python src/html/wsgi.py`'s own local dev server. The *database*
server/account, and which schemas the picker (`index.py`/`?db=`) offers,
are entirely the API's own configuration now
(`PLANETGEN_MYSQL_HOST`/`_PORT`/`_USER`/`_PASSWORD`/`_DATABASE_PREFIX`, or
`config.json`'s `mysql` section, via
`stellarObjects._db.MySQLConfig`/`list_databases`) -- see
[`api.md`](api.md#running-locally) for those; `html/` itself no longer
reads any `PLANETGEN_MYSQL_*` variable at all.

Separately from all of the above, a `config.json`
file at the repo root (a sibling of `../src/html/`, not a file inside `../src/html/`
itself) holds every deployment-level setting in one place -- MySQL
connection details, the site's own `site_name`/`base_url`, and more --
edited once per deployment rather than passed through the vhost config --
see [`config.md`](config.md) for the full field list
and how it relates to the `PLANETGEN_*` environment variables.

## Deploying

1. Copy the repo (or at least `../src/html/`, `src/`,
   `install.sh`, `update.sh`, `setup.py`, and `examples/apache/`) to the
   server, e.g. `/var/lib/planetGen/`. Cloning it there as a git checkout
   (rather than copying a tarball) is what makes `update.sh` possible
   later. A MySQL server (8.0.16+) reachable from this host, with a
   database and account already created, is a separate prerequisite --
   see [`database-schema.md`](database-schema.md).
2. From that directory, run `sudo ./install.sh` -- installs the Python
   package, brings the configured MySQL database's schema up to date
   (a no-op if it's already current -- see
   [`database-schema.md`](database-schema.md)'s "Versioning"),
   pre-fetches the NLTK `words` corpus into a shared world-readable
   location (so it works under Apache's `www-data`, not just whatever
   user happens to run the CLI tools), makes the CGI scripts executable,
   enables Apache's CGI, headers and deflate modules, and sets `../src/html/` ownership for
   Apache via `examples/apache/set-permissions.sh`. See
   [`apache-deployment.md`](apache-deployment.md) for what each step does
   and how to re-run pieces of it individually.
3. `install.sh` prints one remaining manual step: copy
   `examples/apache/planetgen.conf.example` to
   `/etc/apache2/sites-available/planetgen.conf`, edit it (at minimum,
   `ServerName`), then `sudo a2ensite planetgen && sudo systemctl reload
   apache2`. This is deliberately not automated -- the vhost's
   ServerName/TLS/logging are your call, not something to silently create
   or overwrite.

### Updating an existing deployment

Run `sudo ./update.sh` instead of pulling manually. `git pull` on its own
isn't enough -- pulling a changed file rewrites it with whatever mode is
tracked in the repo, silently undoing any executable bit `install.sh`
previously fixed. `update.sh` pulls (refusing to run over uncommitted
local changes, and failing loudly rather than merging if history has
diverged) and then re-runs `install.sh`, so permissions are guaranteed
correct again afterward. Every `install.sh` step is idempotent (the
schema migration and corpus fetch both skip themselves if already
current/present, `chmod +x`/`a2enmod`/permission-setting are all safe to
repeat), and `install.sh` itself
remains safe to run directly any time you want to re-apply everything
without pulling first.

## Local testing without Apache

Every script is a normal CGI program: it reads `QUERY_STRING` (or, for a
`POST` -- what every in-app navigational link now sends, see "How it
works" above -- `CONTENT_LENGTH` bytes of `application/x-www-form-
urlencoded` body from stdin) from the environment and writes an HTTP
response (status + headers + body) to stdout. `nav_params`/
`nav_multi_params` fall back to `QUERY_STRING` whenever the request isn't
a POST at all, so a plain `QUERY_STRING`-only smoke test still works
exactly as before. That makes these scripts runnable directly for a quick
smoke test without standing up Apache at all -- start the API separately
first (see [`api.md`](api.md#running-locally)), then, from the repo root:

```bash
PLANETGEN_API_BASE_URL=http://127.0.0.1:5000/api QUERY_STRING="db=planetgen&id=1" python3 src/html/sector.py
```

The Flask pages need no CGI at all: `python3 src/html/wsgi.py` serves
them (and `/static/`) at `http://127.0.0.1:5000/` alongside the API. Set
`admin_cookie_insecure` in `config.json` to log in over plain HTTP
there.
