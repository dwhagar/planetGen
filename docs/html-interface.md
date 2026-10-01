# planetGen Web Interface

A small web interface (in [`../src/html/`](../src/html/)) for
browsing the MySQL databases described in
[`database-schema.md`](database-schema.md) -- pick a database (schema),
drill into its sectors and star systems, and view (or copy) the
rendered wikitext/Markdown page saved for each one. Sectors and systems
also get their own interactive visualizations: a draggable 3D "Sector
Map" per sector (translucent clouds for any nearby nebula/asteroid field,
glowing points for any nearby black hole/neutron star),
a galaxy-scale "Galaxy Map" (Quadrant/Zone drill-down) for sectors
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
browser today, on nothing more than a web server (Apache2 with mod_wsgi,
or nginx, Caddy or IIS in front of gunicorn or waitress) and Python 3
(see "Locating the database" and "Deploying" below).

## How it works

Every page is served by the same Flask app as the JSON API
(`../src/html/wsgi.py`; the pages live in `../src/html/web/`, see "Flask
pages" below). A page fetches its data through `../src/html/lib/apiclient.py`,
which inside the app dispatches straight through the API's own routes
(no HTTP round trip), and renders it with a Jinja2 template. Links are
plain GET `<a href>`s, so every page is bookmarkable. `nltk` (from
generation) is still needed for a handful of display-formatting constants
a few pages import from `stellarObjects` (spectral-class colors, planet-
class descriptions, unit formatters) -- importing any `stellarObjects`
submodule runs that package's own `__init__.py`, which pulls in `names.py`
and its NLTK corpus dependency. See
[`api.md`](api.md#deploying) for the process
itself, which needs `pymysql`/`DBUtils` and a database account.

| File | Purpose |
|---|---|
| `../src/html/wsgi.py` | The WSGI entry point (loaded by mod_wsgi, gunicorn or waitress): the Flask app serving the API under `/api` and every page. |
| `../src/html/web/` | The pages: routes, templates and helpers (see "Flask pages" below). `web/old_urls.py` answers the old `/<name>.py` CGI URLs with a 301 to the page that replaced them. |
| `../src/html/lib/pagination.py` | The site's one pager, used under every paged table (Browse's two tables, Phenomena, a sector's Contents table, a Galaxy Map Quadrant's sector list, each Search result panel, the admin API key list and the admin stats page's duplicate-names list): a "Showing X-Y of Z" summary, then First/Prev, numbered pages and Next/Last, 50 rows a page. Each table has its own page parameter (e.g. `sectors_page`), a plain GET link, and changing a Search filter starts its results back at page 1. |
| `../src/html/lib/apiclient.py` | The API client every page calls instead of querying MySQL directly -- one typed wrapper function per read endpoint, plus `auth_*` wrappers (the admin pages) supporting POST/DELETE, a request body, and `Cookie`/`Set-Cookie` relay, and `NotFoundError`/`ApiError` (`web/errors.py` turns these into a 404/502 page; `ApiError.status_code` lets `auth_me`/the admin pages branch on a 401 without string-matching). Inside the Flask app it runs in-process (`web/transport.py`); elsewhere it uses HTTP to `PLANETGEN_API_BASE_URL` (default `http://127.0.0.1/api`). Not web-accessible. |
| `../src/html/lib/fmt.py` | HTML-escaping and small formatting helpers (`esc`, `linkify_location`/`nearest_neighbors_location` -- which take a `system_url(system_id)` hook for the neighbour links -- `format_density`, and `static_url`, the versioned `static/` URL every page uses) with nothing to do with fetching data. Not web-accessible. |
| `../src/html/lib/mdconvert.py` | A small, purpose-built Markdown-to-HTML converter for the narrow Markdown subset `StarSystem.__str__` actually generates (headers, pipe tables, paragraphs, `<sup>` exponents) -- not a general-purpose parser. Not web-accessible. |
| `../src/html/lib/galaxymap.py` | Quadrant/Zone classification (`sector_quadrant`/`sector_zone`/`zone_bounds_ly`) shared by the `/galaxy` page's data tables and the sector/sectors pages' "Quadrant N" links back into it -- what's left after this module's former flat-SVG map rendering was superseded by a real 3D scene (`lib/galaxymap3d.py`). Not web-accessible. |
| `../src/html/lib/tilecache.py` | The web layer's on-disk cache (written by the WSGI process serving `/galaxy` and `/galaxy/tiles`) of 3D Galaxy Map tiles: `fetch_tiles` serves each requested tile from `<cache dir>/<db>/<generation>/` when present and asks `GET /api/galaxy/tiles` only for the rest. Every 60 s at most it asks `GET /api/galaxy/changes` which tiles changed (from the rows' `modified_at`) and deletes only those, or starts a new generation when the answer is `full`. The changed keys are passed on to the browser's cache too. Size-capped (`tile_cache.max_mb`, oldest files pruned first); location from `PLANETGEN_TILE_CACHE_DIR`/`tile_cache.dir` (see `config.md`). Fails open: an unwritable directory or bad file just means an API call. Not web-accessible. |
| `../src/html/lib/galaxymap3d.py` | Builds the `/galaxy` page's "Galaxy Map" panel: the canvas/controls/info-panel markup, plus the one starting JSON payload (zoom-range numbers, tile settings, and the zoomed-all-the-way-out view's cube tiles, which `initial_tile_request` picks the same way the client does for every later view) `static/galaxymap3d.js` reads on first paint -- every later payload, as the camera moves, is fetched by that script directly and never passes through this module. `view_radius_bounds` (the zoom floor/ceiling: a couple of sector-widths up to this galaxy's own real outer edge) is a pure function the `/galaxy` view also calls directly, before this module's own panel-rendering function runs. The panel's JSON names the tile URL (`fetchPath`, `/galaxy/tiles`) and a sector-page URL template (`sectorUrl`, from `page_url`) for the info panel's real "View sector" `<a href>`; it carries no database name except `storageKey`, which only namespaces the browser's `localStorage`. Not web-accessible. |
| `../src/html/lib/systemmap.py` | Builds the system page's "System Map" panel: a fixed-size, square true-position plot (real angle from `position_x/y_km`, a shared log-radial scale from anchor -- star, or barycenter for a merged/'close' binary pair -- to body) with one colored/sized marker per body, plus a full ring (not a directional band) for each asteroid belt -- distinct from `starmap.py`'s draggable 3D cube, since this only ever needs a flat top-down projection. A once-only pairwise-repulsion pass (`_relax_markers`) nudges apart any two markers real placement happened to put too close together, real position first, decluttering only where needed. Clicking a planet with moons swaps to a "zoom into this planet's moons" scene (`static/systemmap.js` toggles which `<svg>` scene is visible); clicking anything else fills the info side panel, same click-for-info pattern as `starmap.py`/`sectormap.js`. Adds a small green badge to any body with `life_chemical` set. Also hands over each planet/moon marker's resolved class color (`_class_color`) plus its `atmosphere`/`composition`/`surface_temperature_k` (whether it has a real atmosphere at all, not just the description text) as `data-*`, plus a star's own spectral-type color (`_star_color`, same `data-color` attribute), for `static/systemmap.js`'s own per-marker live sphere rendering. Not web-accessible. |
| `../src/html/lib/systempage.py` | The system page's Python-built HTML (moved from the old `system.py`): the expandable body list (`system_list_html`) and the Stars/Planets & Moons/Asteroid Belts/Comets tables (`stars_html`, `bodies_html`). Escapes every database value itself; `web/system_pages.py` passes the result through `trusted_html`. Not web-accessible. |
| `../src/html/lib/tabledisplay.py` | Computes the same "Star Data"/"Planet Data" display strings once baked into the database's now-removed `table_*`/`binary_table_*` columns, but on demand from the raw numeric columns the system page already has -- reuses `stellarObjects.utils`'s formatters directly. Not web-accessible. |
| `../src/html/lib/starmap.py` | Builds the sector page's "Sector Map" panel (`render_map_panel(link_url, ...)`, where `link_url` is `web.helpers.page_url`; every entry carries a plain `href`): computes every position/size/color/label the map needs (a wedge or fallback-cube outline, one entry per star -- billboarded, radius from `radius_km` square-root scaled against the Sun, color from `star_type`'s spectral letter (`SPECTRAL_CLASS_COLORS`) shaded by `luminosity_w` and nudged by where `temperature_k` falls in that spectral class's range, so "White Giant" reads white and "Blue Giant" reads blue regardless of temperature -- and one per nearby standalone phenomenon, sized by `radius_ly` (always 0 for the point-like types) and positioned directly in the galaxy frame, no rotation needed unlike a star system's sector-local position -- see `queryDb.phenomena_near_sector`) and serializes it as a `<script type="application/json">` block; `static/sectormap.js` is what actually renders it, this module builds no HTML scene of its own. Not web-accessible. |
| `../src/html/lib/navmap.py` | Builds the NAV page's "NAV Map" panel (`render_nav_map_panel(link_url, ...)`; each point is an SVG `<a href>`): a flat, static, top-down SVG plot of the galactic X-Y plane -- origin and destination as labeled points, a dashed line for the direct course, and (when one was found) a solid polyline through the optimal route's intermediate hops. Auto-scaled to whatever points it's given (no fixed sector size to normalize against), with one uniform light-years-per-pixel ratio on both axes so bearings aren't visually distorted, plus a compass arrow along the origin's bearing 000 (toward the frame's center) and a scale-bar legend. Deliberately blind to altitude/z, same as the flat SVG phenomenon Diagram panel (`lib/phenomenonmap.py`) -- the course's mark already covers that axis. Not web-accessible. |
| `../src/html/static/style.css` | Shared stylesheet (CSS custom properties, light/dark via `prefers-color-scheme` or an explicit `data-theme` on `<html>`, card-style panels, phone layout under 40rem), served directly by the web server. |
| `../src/html/static/theme.js` | Loaded on every page, blocking, before `style.css` (`web/templates/base.html`): applies the saved light/dark/system theme before the first paint and drives the header's theme button. See "The page shell" below. |
| `../src/html/static/favicon.svg` | The site icon (a small ringed planet), linked from every page's `<head>`. |
| `../src/html/static/vendor/` | Vendored third-party JS -- currently just `three.module.min.js` (three.js, bundled+minified from the `three` npm package), used by `sectormap.js`/`systemmap.js`. Vendored rather than loaded from a CDN so the pages' `Content-Security-Policy: default-src 'self'` needs no exception; see this directory's own `THIRD_PARTY_NOTICES.txt` for the license and how to rebuild it from a newer release. Served directly, same as `style.css`. |
| `../src/html/static/sectormap.js` | Renders `lib/starmap.py`'s `#starmap-data` JSON as a real WebGL scene (three.js, vendored at `static/vendor/` -- see that directory's `THIRD_PARTY_NOTICES.txt`): a perspective camera, GPU-billboarded sprites (always face the camera by construction, no manual per-frame counter-rotation needed) for stars/nebulae/asteroid fields/black holes/neutron stars, and a wireframe outline for the sector's wedge or fallback cube. Hand-rolled drag-to-rotate/scroll-or-button-to-zoom/click-for-info (a raycast against the sprites, not DOM hit-testing) plus a hidden, keyboard/screen-reader-focusable button list (a canvas has no focusable children of its own the way the old per-star `<div role="button">`s were) so every star/cloud stays reachable without a mouse. Clicking (or activating a fallback-list button) fills the info side panel, with a plain "View system"/"View phenomenon" link. Served directly, same as `style.css`. |
| `../src/html/static/generatebuttons.js` | The admin Generate buttons for a sector that isn't generated yet (it alone, its neighborhood, its column, or its whole shell after a confirm), shared by the Sector Map (`sectormap.js`, a neighbor) and the Galaxy Map (`galaxymap3d.js`, a sector cell). Each is a plain POST form to the Generate page with the CSRF token from the server's `generate` target (`web.helpers.generate_target`), which the pages only send a logged-in admin. Served directly, same as `style.css`. |
| `../src/html/static/galaxymap3d.js` | Renders `lib/galaxymap3d.py`'s `#galaxymap3d-data` JSON as a real WebGL scene (three.js, vendored at `static/vendor/`) -- unlike `sectormap.js`'s camera (always orbiting a fixed origin, every star baked into one page load), this camera's own orbit target moves freely through the galaxy, so most of what it draws is fetched live from `/galaxy/tiles` (plain GET, no database in the URL) (debounced, on every camera move), one fixed cube of space at a time, rather than server-rendered once. Tiles are kept in memory and in the browser's `localStorage`, keyed by the database's content stamp (under `planetgen:tile:<db>:`, the same keys as before the page moved), so panning back over seen space or reloading the page doesn't refetch them. Drag to rotate, scroll or the buttons to zoom; a click centers on the block under the cursor (or the sector at that spot in empty space) and shows its info, and a double-click also zooms in -- by a logarithmic step (`clickZoomFactor`, interpolated between `lib/galaxymap3d.py`'s own `clickZoomFactorMin`/`Max` by the camera's current distance in log space) rather than a flat factor, so a handful of clicks still crosses the whole galaxy while a click near one sector stays fine enough not to overshoot it. Everything is drawn as one solid of blocks the page computes itself (`static/galaxyprisms.js`, built into GPU-ready arrays by `static/galaxyblocks.js` in a Web Worker so zooming never stalls the page; the first frame, and a browser where the worker can't start, build on the page): each block is a power-of-3 cube of whole sectors, the smallest at least `blockMinPx` (4) pixels across at the focus, colored by predicted density, and the Slice button cuts it at the focus layer. Unfilled space is translucent (only its surface is drawn); blocks holding generated sectors (counted from each tile's `filled` summary) are drawn wherever they are, grow more solid and warmer with their filled share, and are fully solid once every sector in them is generated. At one sector per block a generated sector takes its real density's color and its panel links to its sector page; clicking along a ray prefers the nearest block holding generated sectors, and centers on their mean position, so double-clicking zooms toward them. A block's panel gives its ring, layer and slot ranges, exact sector count and generated count; an unfilled sector's gives its address and designation, plus, for a logged-in admin, the Sector Map's Generate buttons (`static/generatebuttons.js`: this sector, its neighborhood, its column, or its whole shell after a confirm), which post to the Generate page. Nebulae and supernova remnants are drawn as soft translucent spheres their real size (each tile lists the ones reaching into it, `queryDb.galaxy_clouds_in_box`), hidden when only a pixel or two across and faded out as the camera nears them; clicking one while it is small enough to aim at shows its type, class and radius and links to its phenomenon page. Zoom steps glide over 160 ms instead of jumping (they jump with `prefers-reduced-motion`), a change of block size crossfades, built views are kept so zooming back is instant, and while idle the page prepares the views and fetches the tiles one zoom step either way. In-plane lines along the sector grid's master wedges (3 from the core, doubling outward, each zone shown once its lines are far enough apart on screen) are toggled by the Wedges button, and a scale readout gives sectors, pc and ly. Served directly, same as `style.css`. |
| `../src/html/static/systemmap.js` | Toggles which System Map `<svg>` scene is visible (the whole-system view, or one per planet's own moon system) and fills the info side panel -- including "Atmosphere"/"Surface composition"/"Surface temperature"/"Life Chemistry" fields -- from a clicked marker's `data-*` attributes. Every visible scene's star/planet/moon markers also get their own live-rendered 3D sphere on `#sysmap-spheres-canvas` (the same vendored three.js build `sectormap.js` uses): one shared WebGL context, redrawn each frame via a scissored sub-viewport per marker (never one `<canvas>`/context per body -- browsers cap concurrent WebGL contexts), each sized and positioned to exactly cover that marker's own `<circle>` and colored by its `data-color`, banded with a tilted ring for a gas giant (`data-bodytype`), and wrapped in a fresnel-glow atmosphere shell (tinted by `data-surfacetemp`) when `data-hasatmosphere` is set -- the one genuinely 3D layer on this otherwise flat-SVG page, drawn behind the SVG so each marker's own stroke/label/life-badge still shows on top. An appearance layer only (no position data, and no cross-page link of its own to navigate) -- a marker whose sphere renders keeps its flat circle's fill transparent (`sysmap-sphere-active`) but never changes its actual plotted position. Served directly, same as `style.css`. |
| `../src/html/static/copycode.js` | The system page's Copy button: copies the generated Wikitext/Markdown out of the code box named by the button's `data-copy-target` (clipboard API, falling back to a selection copy on a plain-HTTP deployment). A separate file because the Content-Security-Policy allows no inline script; the box stays selectable by hand without it. Served directly, same as `style.css`. |

**Why the maps use three.js.** The Sector, System and Galaxy Maps all
draw with three.js (r186), vendored as one minified file. That choice
was checked again on 2026-09-30 for the Galaxy Map:

- Babylon.js ships several megabytes, and deck.gl needs a bundler. This
  project has no build step, and the pages' `default-src 'self'` policy
  favors one vendored file.
- regl or raw WebGPU would mean writing picking, sprites and lighting by
  hand.
- The Galaxy Map's slow part was the JavaScript that lists blocks, not
  the drawing. That work now runs in a Web Worker (`galaxyblocks.js`).

If the map outgrows WebGL, three's own `InstancedMesh`, `BatchedMesh` and
`WebGPURenderer` are the next steps, before any other library.

This project recommends (but doesn't enforce in code) pointing the app at a MySQL
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

Every page extends `web/templates/base.html`, and every HTML response
carries the security headers (`web.SECURITY_HEADERS`, added in
`api/app.py`: a `Content-Security-Policy` of `default-src 'self';
base-uri 'self'; form-action 'self'; frame-ancestors 'none'; object-src
'none'`, plus `X-Frame-Options`, `X-Content-Type-Options` and
`Referrer-Policy`), and a request that came in over HTTPS also gets
`Strict-Transport-Security: max-age=31536000` (`web.STRICT_TRANSPORT_SECURITY`;
never on plain HTTP, so a local `python src/html/wsgi.py` still works).
The shared `<head>` has:

- `<meta name="viewport" content="width=device-width, initial-scale=1">`,
  so phones use the narrow layout (`@media (max-width: 40rem)` in
  `style.css`: wide tables scroll inside their own `.table-scroll` box
  instead of widening the page).
- A `<meta name="description">` (a site-wide default, or `render_page(...,
  description=...)`), `<meta name="color-scheme">`, and two `<meta
  name="theme-color">` tags, one per OS colour scheme.
- The favicon, `static/favicon.svg`.
- `static/theme.js`, a small blocking script placed before the
  stylesheet: it applies the visitor's saved theme before the first
  paint (see below).
- `static/style.css`.

Every `static/` URL, here and in the pages' own `<script>` tags, comes
from `lib/fmt.py`'s `static_url(name)`, which appends `?v=<package
version>` (read from `src/stellarObjects/_version.py`, which the release
Action stamps). A release therefore changes every static URL, and the
web server can let browsers cache them for a year (see
[`deployment/apache.md`](deployment/apache.md#static-files-compression-and-security-headers);
every other guide in [`deployment/`](deployment/README.md) sets the same
rule).
The ES modules that import siblings (`sectormap.js`, `systemmap.js`,
`galaxymap3d.js` -> `bodyRendering.js`, `galaxyprisms.js`,
`galaxyblocks.js`, `vendor/three.module.min.js`) use `await import(...)`
with their own `?v=` (from `import.meta.url`), so each module has exactly
one URL per page and is never loaded twice. The Galaxy Map's worker is
started from `galaxyblocks.js` with the same `?v=`.

**Theme switch.** The header has a "Theme: System / Light / Dark"
button. `static/theme.js` reveals it (it renders `hidden`,
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

## Flask pages

Every page is served by the same Flask app as the API
(`../src/html/web/`, registered by `api/app.py`'s `create_app`). The
pages used to be one CGI script each; the "Replaces" column names the
old script, whose URL now answers with a 301 to the new page (see "Old
URLs" below).

| URL | Replaces | Shows |
|---|---|---|
| `/` | `index.py`, `browse.py` | Every sector and every standalone system, each table paged on its own (`?sectors_page=N`, `?standalone_page=N`). |
| `/sectors` | `browse.py#sectors` | The sectors table alone. |
| `/systems` | `browse.py#standalone-systems` | The standalone systems table alone. |
| `/galaxy` | `galaxy.py` | The 3D Galaxy Map plus the Quadrant summary; `?quadrant=I\|II\|III\|IV` lists that Quadrant's placed sectors nearest the core first (`?page=N`). The info panel's "View sector" is a plain link. |
| `/galaxy/stage?at=...` | `galaxy_views.py` | JSON for the map's drill-down: one stage's generated counts (`at=m.ring.wedge.slab`, or none for the galaxy; see `/api/galaxy/stage`). Through `lib/tilecache.py`'s disk cache (`fetch_stage`), where a sector change deletes only its chain's stages; 400 on a malformed key, 502 on an API failure, both `{"error": ...}`. Shares the `galaxy_tiles` rate limit. `Cache-Control: no-store`. |
| `/galaxy/locate?q=...` | `galaxy_views.py` | JSON for the map's address bar: sectors and star systems named like `q`, each with its sector address (see `/api/galaxy/locate`). A blank `q` asks the API nothing; 502 on an API failure, as `{"error": ...}`. Shares the `search` rate limit. |
| `/galaxy/tiles?tiles=...` | `galaxy_tiles.py` | JSON for the map's script: the requested cube tiles (`tiles=level/ix/iy/iz,...`) and, given the browser cache's `stamp`, the changed tiles since. Through `lib/tilecache.py`'s disk cache; 400 on a malformed request, 502 on an API failure, both `{"error": ...}`. `Cache-Control: no-store`. |
| `/system/<id>` | `system.py` | One star system: badges, "Navigate from/to here" (systems in a sector), nearest-neighbour location links, the System Map, the expandable body list, `?code=wikitext\|markdown` code views with a Copy button, the Stars/Planets/Belts/Comets tables, and for an admin the "Upload to Wiki" form (see below). |
| `/phenomena` | `phenomena.py` | Every exotic phenomenon, paged with `?page=N`. |
| `/phenomenon/<type>/<id>` | `phenomenon.py` | One phenomenon's data table and AU-scale diagram, with "Navigate from/to here". `<type>` is one of `nebula`, `asteroid_field`, `black_hole`, `neutron_star`, `supernova_remnant`, `rogue_planet`, `interstellar_comet`, `quasar`; anything else is a 404. |
| `/sector/<id>` | `sector.py` | One sector: badges, the 3D Sector Map, and its Contents table (systems and nearby phenomena, nearest the center first, `?contents_page=N`); admin forms (wiki upload, generate neighborhood). |
| `/nav` | `nav.py` | The NAV route planner; see "The NAV page's URLs" below. |
| `/search` | `search.py` | Faceted search (see below). |
| `/login` | `login.py` | The admin login form (`?next=<local path>` to return to afterwards). |
| `/logout` | `logout.py` | `GET` asks to confirm and changes nothing; the button `POST`s to end the session. |
| `/account` | `changecreds.py` | Change the admin username and password. |
| `/admin` | `admin.py` | API keys (list, create, revoke; `?keys_page=N`) and a sector's manual wiki link. |
| `/admin/stats` | `adminstats.py` | Server health and database stats, and every name made unique (`?names_page=N`). |
| `/admin/generate` | (new) | Admins only: generate, plan or reset the galaxy from the browser (see below). |
| `/admin/generate/system` | (new) | Admins only: one star system with every `generate.py system` option, shown as Markdown or wikitext and never saved (see below). |

**System page wiki upload.** For an admin (`current_admin()`), the
system page asks `GET /api/wiki-config` and offers an "Upload to Wiki"
form for each configured backend not yet uploaded to. The form POSTs to
`/system/<id>` with `{{ csrf_field() }}`; the view calls
`apiclient.upload_system_to_wiki` with the visitor's cookie and always
answers `303 See Other` back to `/system/<id>?wiki=<outcome>#wiki-upload`
(POST-redirect-GET, so a reload never re-posts). `<outcome>` is one of
`uploaded`, `exists` (409), `invalid` (400 or an unknown backend),
`forbidden` (401/403, e.g. default credentials not yet changed),
`unconfigured` (501) or `failed`; the page shows a fixed message for
each (and only to an admin), never text from the query string or the
API. A POST from a visitor who is not an admin gets a 403 page.

**Old URLs.** `web/old_urls.py` answers a GET of an old CGI page's URL
(`/<name>.py`) with `301 Moved Permanently` to the page that replaced it,
so old links and bookmarks still land somewhere sensible. The query
string is carried over minus `db` and empty values (the new pages kept
the old parameter names: `sectors_page`, `quadrant`, `page`, the search
filters, ...). `/sector.py?id=N` goes to `/sector/N` (keeping
`contents_page`), `/system.py?id=N` to `/system/N` (keeping `code`),
`/phenomenon.py?type=T&id=N` to `/phenomenon/T/N`, each falling back to
its list page without a valid id; `/nav.py`'s old `from=12&from_kind=
phenomenon&from_type=nebula` becomes `from=nebula:12`. `galaxy3d.py` and
`galaxy_view.py` (older Galaxy Map scripts) go to `/galaxy`. Any other
`.py` name, `wsgi.py` included, gets the 404 page. An old form POSTed to
a `.py` URL is not replayed (it fails the CSRF check like any stale
form).

The Galaxy Map's disk tile cache (`lib/tilecache.py`) is written by the
WSGI daemon. A tile cache location or size set only with `SetEnv
PLANETGEN_TILE_CACHE_DIR`/`_MAX_MB` in the vhost does not reach it: put it
in `config.json`'s `tile_cache` (or the app server's own environment) instead.

### The sector page

`/sector/<id>` (`web/sector_page.py`, `templates/sector.html`) shows the
sector's badges (cube edge, counts, a link to its Galaxy Map quadrant),
the interactive Sector Map (`lib/starmap.py` data rendered by
`static/sectormap.js`) and one Contents table of its systems and the
phenomena near it, nearest the center first, 50 per page
(`?contents_page=N`). Every map entry carries a plain `href`: the info
panel's "View system/phenomenon/sector" button and the `<noscript>` list
are ordinary links.

Admin actions are POST forms to the same URL with `{{ csrf_field() }}`
and an `action` field: `upload_wiki` (any logged-in admin, while the
sector has no wiki page; offers the configured backends) and
`generate_neighborhood` (an admin with current credentials, on a
galaxy-placed sector). A successful action flashes its message and
answers `303 See Other` back to `/sector/<id>` (POST-redirect-GET, so a
reload never repeats it); a failed one re-renders the page with the error
next to its form. A POST without an admin session changes nothing and
redirects back.

### The NAV page's URLs

`/nav` (`web/nav_page.py`, `templates/nav.html`) is all GET, so every
step can be bookmarked or shared:

| URL | Shows |
|---|---|
| `/nav` | Choose a starting sector. |
| `/nav?from_sector=5` | Choose a starting system in sector 5. |
| `/nav?from=system:12` | Destination pickers: the other systems in its sector, and (when that sector is galaxy-placed) a sector-then-system picker for another sector (`&to_sector=9`). |
| `/nav?from=system:12&to=system:40` | The direct course (distance, "bearing mark mark" and its frame, warp and fold travel times), the NAV Map and the optimal route. |
| `/nav?from=nebula:3&to=system:40` | The same with a phenomenon endpoint. |
| `/nav?to=black_hole:3` | "Navigate to here": the origin picker, carrying `to` along. |

An endpoint is `<kind>:<id>`: `system:<id>`, or a phenomenon's type and
id (`nebula`, `asteroid_field`, `black_hole`, `neutron_star`,
`supernova_remnant`, `rogue_planet`, `interstellar_comet`, `quasar`). A bare number
means a system. A malformed or unknown endpoint is a 404; two endpoints
that exist but cannot be navigated between (`GET /api/nav` answers 400)
get a message on the page. Build links with `web.nav_page.nav_url(origin,
destination)` and `endpoint(kind, id)`, or in a template
`page_url('nav', **{'from': 'system:12'})`. The older parameters
(`from_id`/`to_id`, or a numeric `from`/`to` with `from_kind`/`to_kind` =
`phenomenon` and `from_type`/`to_type`) answer 301 to the canonical URL.
A phenomenon is never listed in the pickers; it becomes an endpoint only
through a link from its own page.

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
search without them.

For a visitor:

- **Plain GET links and bookmarkable URLs.** Links between pages are
  ordinary `<a href>`s, so Back, reload, open-in-new-tab, bookmarks and
  sharing all work.
- **No database in the URL.** The pages show one database, taken from
  config (`mysql.database` / `PLANETGEN_MYSQL_DATABASE`), never from the
  request.
- **A header** with the site name, the sections
  (Galaxy, Sectors, Systems, Phenomena, Nav) with `aria-current="page"`
  on the current one, a search box, and Login or Admin/Stats/Generate/Logout
  plus the theme button. Below 56rem (92rem for a logged-in admin, whose
  header has more links) the sections, search and account links fold into
  a native `<details>` "Menu" (works without JavaScript). A
  "Skip to content" link comes first, and every page but the home page
  has a breadcrumb trail.
- **Fast**: no Python process start per page, and no HTTP call from the
  page back to the API (see below). Visitors without a session cookie
  cost no login lookup at all.

**The Generate page** (`/admin/generate`, `web/generate_page.py` +
`templates/generate.html`) runs the galaxy tools an admin used to run in a
terminal on the server, as background jobs:

| Form | Runs |
|---|---|
| New galaxy | `src/resetDb.py --yes`, then `generate.py plan`, then `generate.py galaxy` around a random start. |
| Generate sectors | `generate.py galaxy` in any of its modes: around a random start, a whole ring at one layer (`--ring --layer`, with `--limit`, or `--yes` for a very large one), around a sector (`--center-sector --radius-pc`), one address (`--ring --layer --slot`, with an optional neighborhood radius), a column (`--ring --slot --column`), or a shell (`--ring --shell`, marked not recommended). The Sector Map's Generate buttons on an unfilled neighbor post straight to this form. |
| Plan the galaxy | `generate.py plan` with the galaxy shape fields. |
| Reset | `src/resetDb.py --yes`. |

The number fields have upper bounds, the same ones `generate.py` checks
(`src/stellarObjects/generationLimits.py`): a radius of at most 200 pc,
rings up to 100,000, and at most 500 orbital slots on the one-off system
page.

New galaxy and Reset delete every generated row, so both need the
database name typed back. Every job writes the database this site shows,
passed to the child as `PLANETGEN_MYSQL_*` environment variables (so the
MySQL account needs the generator's grants, including `DROP` for
`TRUNCATE`). One job runs at a time; the page shows its step, a progress
bar (from `generate.py`'s `PLANETGEN_PROGRESS_FILE`, see
`stellarObjects/progressFile.py`), elapsed time and live output
(`static/generatejobs.js` polls `/admin/generate/status`), with a Cancel
button. The last jobs are listed with their full output at
`/admin/generate/jobs/<id>`.

**The one-off system page** (`/admin/generate/system`,
`web/system_page.py` + `templates/generate_system.html`, linked from the
Generate page) offers every `generate.py system` option: the ten
force/forbid choices, name, star type, age, orbital slots, the flavor
overrides, a pasted `--system-file` JSON, Markdown or wikitext, and the
`--debug` narration. It runs `generate.py system --output FILE` in a
temporary directory and waits for it (a system takes about a second), so
nothing touches the database. The result shows in a code box with Copy
and Download buttons (Download posts the text back to
`/admin/generate/system/download`, which returns it as a `.md`/`.wiki`
file), plus a rendered preview for Markdown.

A job is started as `python3 src/jobRunner.py <job dir>` in its own
session (on Windows, a detached process in its own process group,
broken away from the server's job object where the server allows it),
so it outlives the request and a graceful reload (Apache's, or
gunicorn's). A full stop or restart of the service under systemd (which
stops everything in the service's cgroup), or an IIS app pool recycle
that doesn't allow breakaway, does stop it; the page then shows it as
interrupted. Cancel writes a `cancel` file into the job's directory; the
runner sees it within a quarter second and stops the running step's
whole process tree (`os.killpg` on POSIX, `taskkill /T /F` on Windows).
Liveness comes from `/proc` on Linux, `os.kill(pid, 0)` on other POSIX
systems, and `OpenProcess`/`GetExitCodeProcess` on Windows.
Jobs live under `jobs.dir` (`docs/config.md`). A visitor who isn't a
logged-in admin is sent to the login page, and POSTs and status requests
without an admin session get a 403.

### How a Flask page is built

```
src/html/web/
  __init__.py     blueprint `web`, template globals, init_app()
  views.py        the routes
  helpers.py      db_name, page_url, crumb, render_page, trusted_html,
                  pager, current_admin
  transport.py    in-process transport for lib/apiclient.py
  csrf.py         CSRF tokens for POST forms
  errors.py       HTML 404/502/500 pages
  old_urls.py     301s from the old /<name>.py URLs
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
- `page_url(name, **params)`: the URL of any page by endpoint name
  (`url_for("web.<name>")`).
- `db_name()`: the one database. Never read `db` from the request.
- `trusted_html(html)`: passes HTML built by a `lib/` renderer (the star
  and system maps, `fmt.format_density`, ...) through unescaped. Those
  renderers escape their own inputs; never wrap request or database text
  in it directly.
- `pager(page_param, page, total, anchor=None, label="Pages", keep=None)`:
  the shared pager (`lib/pagination.render_pagination`) linking
  back to the current path with `?<page_param>=N`; `keep` carries other
  query parameters (another table's page number).
- `current_admin()`: the logged-in admin or `None`, one lookup per
  request (cached on `g`). Also available in templates.
- In templates: `static_url("x.js")` (`/static/x.js?v=<version>`, from
  `fmt.static_url`), `page_url(...)`, `current_admin()`,
  `csrf_field()`, `site_name`, and blocks `head` (extra `<script>`/
  `<link>`, e.g. a map module), `heading`, `subhead` and `content`.

**Data in-process.** Views call `lib/apiclient.py` functions
(`get_sector(db, id)`, `auth_me(cookie_header)`, ...).
Inside a Flask request, `web/transport.py` dispatches each call straight
through the app's own `/api` routes (same validation, auth and JSON as
over HTTP, in a fresh app context so database connections open and
close per call) instead of making an HTTP request to itself. Outside a
Flask request (a script) the functions use HTTP.
These in-process calls are exempt from the API's app-wide default rate
limit (a page view is not API abuse); a route's own limit (login,
writes) still applies to the visitor's address. Page views themselves
have their own per-IP limits; see "Rate limits" below.

**Forms.** Any POST/PUT/PATCH/DELETE to a path outside `/api` must carry
the CSRF token: put `{{ csrf_field() }}` inside the `<form>`. The token
is an HMAC (keyed with `secret_key` from `config.json`, see
[`config.md`](config.md)) of a random value in the `pg_csrf` cookie
(HttpOnly, SameSite=Strict) together with a hash of the admin session
cookie (empty when not logged in), so a token only works for the login
session it was rendered for; a missing or wrong token gets a 400 page
and the view never runs. Logging in, logging out or changing credentials
changes the session, so a form left open from before then fails once
and works after a reload.

**The admin pages** (`web/admin_pages.py`). `/login`, `/account` and `/logout` call `apiclient.auth_login`,
`auth_change_credentials` and `auth_logout` in-process, and the
`Set-Cookie` headers the API's `/api/auth/*` routes answer with are
copied onto the page's response unchanged, so the session cookie keeps
the API's HttpOnly/Secure/SameSite=Strict attributes, and `POST
/api/auth/login`'s own per-address limit (10 a minute) still applies to
the visitor (a limited login shows "Too many login attempts" with a
429). The rules:

- Every form that changes something (log in, log out, change
  credentials, create or revoke an API key, set a sector's wiki link) is
  a POST with `{{ csrf_field() }}`, answered with a 303 to a GET page. A
  result the next page has to show (a new API key, shown once; a
  message; an error) travels in `pg_flash`, a signed (`secret_key`),
  HttpOnly, SameSite=Strict cookie with `Path=/admin` (so a new key's
  value is never sent with requests to the rest of the site) that the
  next `/admin` view reads and deletes, never in the URL. The one exception is a wrong password on
  `/login` or `/account`: nothing changed, so the form is shown again
  directly with the error and the username typed.
- `GET /logout` never logs out: it shows a "Log out" button.
- A visitor who isn't logged in is sent to `/login?next=<the page>`, and
  an admin still on the seeded first login (random password printed once
  by the installer, see [`api.md`](api.md#the-first-admin-login)) to
  `/account?next=...`.
  `next` is followed only when it is a local path (starts with a single
  `/`, no scheme, host, backslash, whitespace or control character);
  anything else goes to `/admin`.
- Admin responses are `Cache-Control: no-store`.
- `/admin/stats` and the wiki-link form use the site's configured
  database (no database picker or field any more).

**Errors.** An `apiclient.NotFoundError` becomes a 404 page, an
`apiclient.ApiError` a 502 page (without the API's detail), anything
else a 500 page saying only "An unexpected error occurred."; the
traceback goes to the app server's error log (Apache's, under mod_wsgi)
and the debug log, never the page.
Unknown URLs get the HTML 404 page; `/api/...` keeps its JSON errors.

**Headers.** Every HTML response carries `web.SECURITY_HEADERS` (see
"The page shell" above), set in `api/app.py`; JSON keeps `default-src
'none'`.

**Browser checks.** `src/tests/test_web_a11y.py` loads every GET route
of the `web` blueprint (read from the app's `url_map`, so a new page is
covered automatically) in headless Chromium at 390px and 1280px,
in the light and dark color schemes, against a small database generated
by `generate.py`. Each page must have no serious or critical axe-core
violations of the WCAG 2.1 A/AA rules, no horizontal page scroll (wide
tables scroll inside `.table-scroll`), no console errors, failed
requests or CSP violations, a skip link that is the first Tab stop, and
an `aria-current="page"` marker; at 390px the header Menu is opened and
checked too. A route parameter needs a sample value in the test's
`sample_params` (the test says so when one is missing); a page that
redirects anonymous visitors is checked as a logged-in admin (the old
`/<name>.py` routes only redirect, so they are skipped). It needs
`pip install -e ".[browser]"` plus `python -m playwright install
chromium` (or `PLAYWRIGHT_BROWSERS_PATH` pointing at an existing
Chromium) and the MySQL test server, and skips without them; CI runs it
in its own `browser-a11y` job. `PLANETGEN_A11Y_SCREENSHOTS=<dir>` saves
a screenshot of every page checked. axe-core is vendored in
`src/tests/vendor/axe-core/` (MPL-2.0); to update it, copy `axe.min.js`
and the license files from the npm package.

### Rate limits

Every page view does real work on one of a fixed number of WSGI threads,
so the pages are rate-limited per client IP (`api/limiter.py`'s
`page_limit`, set from `config.json`'s `ratelimit.pages`, see
[`config.md`](config.md)):

| Name | Covers | Default |
|---|---|---|
| `search` | `/search` | 30 per minute |
| `galaxy` | `/galaxy` (the Galaxy Map page) | 60 per minute |
| `galaxy_tiles` | `/galaxy/tiles` (fetched as the map's camera moves) | 600 per minute |
| `health` | `/api/health` | 60 per minute |
| `other` | every other page, all counted together | 300 per minute |

An empty value turns that limit off. These replace the API's default
limits (`ratelimit.default`) for the pages, and the API calls a page makes
in-process are not counted against either. Over a limit, a page answers
`429` with the site's HTML error page ("Too many requests"), while
`/galaxy/tiles` (read by the map's script) and everything under `/api/`
answer JSON; both carry `Retry-After`. With more than one WSGI process,
`ratelimit.storage_uri` needs a shared backend for the counts to add up
(see [`api.md`](api.md#rate-limiting)).

## Locating the database (and the API)

The pages reach the API in-process, so nothing needs pointing at it for
a normal deployment. `PLANETGEN_API_BASE_URL` (or `config.json`'s
`api_base_url`; default `http://127.0.0.1/api`) only matters to
`lib/apiclient.py` used outside the app (a script). The *database*
server/account is the app's own configuration
(`PLANETGEN_MYSQL_HOST`/`_PORT`/`_USER`/`_PASSWORD`/`_DATABASE`, or
`config.json`'s `mysql` section, via `stellarObjects._db.MySQLConfig`) --
see [`api.md`](api.md#running-locally) for those.

Separately from all of the above, a `config.json`
file at the repo root (a sibling of `../src/html/`, not a file inside `../src/html/`
itself) holds every deployment-level setting in one place -- MySQL
connection details, the site's own `site_name`/`base_url`, and more --
edited once per deployment rather than passed through the vhost config --
see [`config.md`](config.md) for the full field list
and how it relates to the `PLANETGEN_*` environment variables.

## Deploying

[`deployment/README.md`](deployment/README.md) compares every supported
platform (Apache, nginx or Caddy on Linux; IIS, Caddy or Apache on
Windows; macOS; a VPS) and links a guide for each. The steps below are
the reference setup, Apache2 with mod_wsgi on Debian or Ubuntu
([`deployment/apache.md`](deployment/apache.md)).

1. Copy the repo (or at least `../src/html/`, `src/`,
   `install.sh`, `update.sh`, `setup.py`, and `examples/apache/`) to the
   server, e.g. `/var/lib/planetGen/`. Cloning it there as a git checkout
   (rather than copying a tarball) is what makes `update.sh` possible
   later. A MySQL server (8.0.16+) reachable from this host, with a
   database and account already created, is a separate prerequisite --
   see [`database-schema.md`](database-schema.md).
2. From that directory, run `sudo ./install.sh` -- installs the Python
   package (with pip, or on an externally managed Python such as Ubuntu
   24.04+'s, from apt packages, with system-wide pip only for anything apt
   lacks or ships too old; see
   [`deployment/apache.md`](deployment/apache.md#managed-python)), brings the configured MySQL database's schema up to date
   (a no-op if it's already current -- see
   [`database-schema.md`](database-schema.md)'s "Versioning"),
   pre-fetches the NLTK `words` corpus into a shared world-readable
   location (so it works under Apache's `www-data`, not just whatever
   user happens to run the CLI tools), makes the shell scripts executable,
   enables Apache's wsgi, headers and deflate modules, and makes `../src/html/`
   `root:<apache group>` (readable, not writable, by Apache) and `config.json`
   mode 640 via `examples/apache/set-permissions.sh`. See
   [`deployment/apache.md`](deployment/apache.md) for what each step does
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
diverged) and then checks everything the site needs without
reinstalling anything that's already there: the executable bits and
permissions, each Python library (installing only one that's missing,
too old or broken, see
[`deployment/apache.md`](deployment/apache.md#managed-python)), the NLTK
corpus, the schema migration, Apache's modules, the cache/jobs
directories and the debug log, and finally that the web app imports as
Apache's user. `install.sh` remains safe to run directly any time you
want a full reinstall without pulling first.

## Local testing without a web server

`python3 src/html/wsgi.py` serves the pages (and `/static/`) at
`http://127.0.0.1:5000/` alongside the API. Set `admin_cookie_insecure`
in `config.json` to log in over plain HTTP there.
