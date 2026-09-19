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
directly itself, though [`../src/html/search.py`](#how-it-works) does
provide a faceted/name search (built on `GET /api/search`, same as every
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
| `../src/html/index.py` | Renders `browse.py`'s own content directly, in place (`browse.handler`, called in-process -- see that page's module docstring for why not a redirect), for the first schema `GET /api/databases` returns. |
| `../src/html/browse.py` | A chosen database's sectors and standalone systems, linking each galaxy-placed sector's row into `galaxy.py`. |
| `../src/html/galaxy.py` | Galaxy-scale "Galaxy Map": every sector actually placed in the galaxy (`GET /api/galaxy/sectors`), grouped into four azimuthal Quadrants (I-IV) and concentric Rings (`shell_index` bands) -- the full-galaxy default view is a per-Quadrant summary table; a `quadrant` of `I\|II\|III\|IV` drills into a full, distance-sorted sector list for that Quadrant. Also plots every galaxy-placed nebula/asteroid field/black hole/neutron star (`GET /api/galaxy/phenomena`) as a small fixed-size dot -- a nebula/asteroid field's own real physical extent is instead shown on the Sector Map of any sector it reaches into (a black hole/neutron star stays point-like there too; see `sector.py` below). Rendering/geometry lives in `lib/galaxymap.py`. Sectors never placed in the galaxy (the common case for a plain `sectorGen.py` run) stay out of this page, visible only via `browse.py`'s flat sector table. |
| `../src/html/sector.py` | One sector's systems (name, quadrant, star type), plus an interactive 3D "Sector Map" of the sector's cube -- drag to rotate, scroll (or the +/- buttons) to zoom, one dot per placed system (two, overlapping, for a binary) sized by star radius and colored by spectral type/brightness, plus a translucent spheroid cloud for every nebula/asteroid field (or a small glowing point for every black hole/neutron star) whose real galaxy-frame sphere reaches into this sector's cube (`queryDb.phenomena_near_sector`, part of `GET /api/sectors/<id>`); clicking a dot, cloud, or point fills an info panel with a link to `system.py` for a star system, via `lib/starmap.py` + `static/sectormap.js`. Shows a "View on Wiki" link once `sectors.wiki_url` is set; an admin session additionally gets an "Upload to Wiki" form (`POST /api/sectors/<id>/wiki`) while it isn't. |
| `../src/html/system.py` | One system's stars/planets/moons/belts, plus its description -- rendered as HTML from `markdown_content` by default (`view=rendered`), with the original raw wikitext/Markdown source (`view=source`/`format=...`) still available for copy-pasting into a wiki, *unless* this system already has a wiki page (`wikijs_url`/`mediawiki_url`), in which case the whole Description section becomes a link to it instead (opening in a new tab). Also renders the interactive "System Map" panel (`lib/systemmap.py` + `static/systemmap.js`) whenever the system has stars -- a true top-down plot of every body's *real* position (real angle from `position_x/y_km`, log-scaled distance from whatever it orbits) with one clickable marker per body, colored by planet class and sized by `radius_km`; a binary pair's two stars (and, for a 'wide' pair, each star's own independent planets) are placed at their real mass-weighted offsets from the system's shared barycenter (`binary_mutual_position_x/y/z_km`). Clicking a planet with moons zooms into that planet's own moon system, and any body with `life_chemical` set gets a small green "supports life" badge (surfaced in the info panel's "Life Chemistry" field too). Links to `nav.py` (a "Navigate from here" button) whenever the system is assigned to a sector. An admin session gets an "Upload to Wiki" form (`POST /api/systems/<id>/wiki`) offering whichever backend(s) are configured and not yet uploaded to. |
| `../src/html/nav.py` | Course, distance, and optimal route between two systems -- calls `GET /api/nav` (`queryDb.nav_between`; see [`api.md`](api.md#nav) for the availability rules and course convention). Reachable either via `system.py`'s "Navigate from here" button or the sidenav's own "Nav" link with no origin known yet, which shows a two-step sector-then-system picker (`GET /api/sectors` then that sector's own systems) to choose one. Without a destination yet, shows a dropdown of the origin's own sector-mates plus, when that sector has a galaxy placement, the same two-step sector-then-system picker (scoped to `to_sector`) for a cross-sector destination (there's no bounded way to offer every galaxy-placed system in one dropdown). With a destination, shows the direct course (distance/azimuth/altitude/warp-1-3-6-9 travel times), the NAV Map (a flat, top-down SVG plot of the origin/destination/route hops -- `html/lib/navmap.py`), and the optimal route via adjacent systems, each hop linking to `system.py`. |
| `../src/html/phenomena.py` | Flat, paginated list of every exotic phenomenon (nebula/asteroid field/black hole/neutron star) across every sector, regardless of galaxy placement -- `GET /api/phenomena`. Each row links to `phenomenon.py`. |
| `../src/html/phenomenon.py` | One phenomenon's detail/info page -- `GET /api/phenomena/<type>/<id>`, this project's first per-phenomenon page (previously a phenomenon had no page of its own, only a hover tooltip on the Sector/Galaxy Map). Reached from `phenomena.py`'s listing, or directly from a Sector Map/Galaxy Map marker (`lib/starmap.py`/`lib/galaxymap.py`). |
| `../src/html/search.py` | Faceted search: click-to-filter tag buttons for object type, star spectral/luminosity class, and planet class/body type/supported life chemistry -- with a separate, identically-shaped set of tags for moons, since planets and moons live in their own tables (schema v2) and a "Class D" tag only ever means one or the other -- built only from values actually present in the chosen database (`GET /api/search`, which owns the query logic; this page just renders it). Plus a name search (with HTML5 `<datalist>` autocomplete, no JavaScript) across sectors, star systems, stars, planets, and moons. Asteroid belts have no name of their own, so they're reachable only via the "Asteroid Belt" object-type tag. Every tag toggle/filter-removal/search submit is a POST form carrying every other currently active filter forward, same as every other page's own navigation (see "How it works" above). |
| `../src/html/login.py` | Admin login form -- `POST /api/auth/login`, relaying the session cookie it sets back to the browser. Redirects to `changecreds.py` (still on the seeded `admin`/`password` default) or `admin.py` on success; re-renders the form with an inline error for a wrong username/password. See [`api.md`](api.md#authentication). |
| `../src/html/changecreds.py` | Change the logged-in admin's username/password (`POST /api/auth/change-credentials`) -- reached both by `login.py`'s forced redirect (default credentials) and voluntarily (an already-"fresh" admin rotating credentials). Always requires the current password. |
| `../src/html/admin.py` | Protected admin landing page: API key list/create/revoke (`GET`/`POST /api/auth/api-keys`, `DELETE /api/auth/api-keys/<id>`) -- keys are used as `Authorization: Bearer` credentials against the write endpoints (sector/system create/update/delete), not exercised as a web form here. The one exception is a small form to manually set/clear a sector's `wiki_url` (`PATCH /api/sectors/<id>`, by database name + sector id, via the session cookie). Redirects to `login.py`/`changecreds.py` if not authenticated / still on default credentials. |
| `../src/html/logout.py` | Ends the current admin session (`POST /api/auth/logout`, called server-side) and redirects to `index.py`. A plain sidenav link, not a form -- nothing to confirm, no request body needed. |
| `../src/html/lib/apiclient.py` | The `GET /api/...` HTTP client (stdlib `urllib` only) every page above calls instead of querying MySQL directly -- one typed wrapper function per read endpoint, plus `auth_*` wrappers (`login.py`/`changecreds.py`/`admin.py`) supporting POST/DELETE, a request body, and `Cookie`/`Set-Cookie` relay, and `NotFoundError`/`ApiError` (`lib/page.py`'s `run` turns these into a 404/502 page; `ApiError.status_code` lets `auth_me`/the admin pages branch on a 401 without string-matching). `PLANETGEN_API_BASE_URL` (default `http://127.0.0.1/api`) is where it looks for the API. Not web-accessible. |
| `../src/html/lib/fmt.py` | HTML-escaping and small formatting helpers (`esc`, `linkify_location`, `format_density`) with nothing to do with fetching data -- what's left of the old `dbutil.py` once its database-access functions moved into `apiclient.py`/the API itself. Also `post_link` (a same-effect, no-JS-required replacement for `<a href="page.py?...">` that posts its params as hidden fields instead -- see "How it works" above) and `data_nav_params` (its JS-required counterpart for a marker embedded in an SVG map, since a `<form>` can't nest inside one -- paired with `static/navform.js`). Not web-accessible. |
| `../src/html/lib/page.py` | Shared CGI response/HTML-shell helpers -- `query_params`/`form_params`/`form_multi_params` (GET/POST parsing) and `nav_params`/`nav_multi_params` (POST body when present, else the GET query string -- what every page reads a `post_link`-followed link's params back with), `incoming_cookie_header` (relays the browser's own `Cookie` header to the API, unparsed), `send_headers`/`render`/`redirect` (all three accept `set_cookie_headers` to relay the API's own `Set-Cookie` back), and the sidenav's Search/Galaxy/Sectors/Systems/Nav/Phenomena/Login/Admin/Logout items (via `apiclient.auth_me` and `fmt.post_link`). Not web-accessible. |
| `../src/html/lib/mdconvert.py` | A small, purpose-built Markdown-to-HTML converter for the narrow Markdown subset `StarSystem.__str__` actually generates (headers, pipe tables, paragraphs, `<sup>` exponents) -- not a general-purpose parser. Not web-accessible. |
| `../src/html/lib/galaxymap.py` | Builds `galaxy.py`'s Quadrant/Ring map geometry and rendering -- Quadrant/Ring classification (`sector_quadrant`/`sector_ring`), the per-Quadrant/full-galaxy panel markup, and every galaxy-placed nebula/asteroid field/black hole/neutron star's small dot (`_phenomenon_elements`) `galaxy.py` supplies data around. Not web-accessible. |
| `../src/html/lib/systemmap.py` | Builds `system.py`'s "System Map" panel: a fixed-size, square true-position plot (real angle from `position_x/y_km`, a shared log-radial scale from anchor -- star, or barycenter for a merged/'close' binary pair -- to body) with one colored/sized marker per body, plus a full ring (not a directional band) for each asteroid belt -- distinct from `starmap.py`'s draggable 3D cube, since this only ever needs a flat top-down projection. A once-only pairwise-repulsion pass (`_relax_markers`) nudges apart any two markers real placement happened to put too close together, real position first, decluttering only where needed. Clicking a planet with moons swaps to a "zoom into this planet's moons" scene (`static/systemmap.js` toggles which `<svg>` scene is visible); clicking anything else fills the info side panel, same click-for-info pattern as `starmap.py`/`sectormap.js`. Adds a small green badge to any body with `life_chemical` set. Also hands over each planet/moon marker's resolved class color (`_class_color`) plus its `atmosphere`/`composition`/`surface_temperature_k` (whether it has a real atmosphere at all, not just the description text) as `data-*`, for `static/systemmap.js`'s own 3D body-preview sphere. Not web-accessible. |
| `../src/html/lib/tabledisplay.py` | Computes the same "Star Data"/"Planet Data" display strings once baked into the database's now-removed `table_*`/`binary_table_*` columns, but on demand from the raw numeric columns `system.py` already has -- reuses `stellarObjects.utils`'s formatters directly. Not web-accessible. |
| `../src/html/lib/starmap.py` | Builds `sector.py`'s "Sector Map" panel: computes every position/size/color/label the map needs (a wedge or fallback-cube outline, one entry per star -- billboarded, radius from `radius_km` square-root scaled against the Sun, color from `star_type`'s spectral letter (`SPECTRAL_CLASS_COLORS`) shaded by `luminosity_w` and nudged by where `temperature_k` falls in that spectral class's range, so "White Giant" reads white and "Blue Giant" reads blue regardless of temperature -- and one per nearby nebula/asteroid field/black hole/neutron star, sized by `radius_ly` (always 0 for a black hole/neutron star) and positioned directly in the galaxy frame, no rotation needed unlike a star system's sector-local position -- see `queryDb.phenomena_near_sector`) and serializes it as a `<script type="application/json">` block; `static/sectormap.js` is what actually renders it, this module builds no HTML scene of its own. Not web-accessible. |
| `../src/html/lib/navmap.py` | Builds `nav.py`'s "NAV Map" panel: a flat, static, top-down SVG plot of the galactic X-Y plane -- origin and destination as labeled points, a dashed line for the direct course, and (when one was found) a solid polyline through the optimal route's intermediate hops. Auto-scaled to whatever points it's given (no fixed sector size to normalize against), with one uniform light-years-per-pixel ratio on both axes so azimuth angles aren't visually distorted, plus a `+X` compass tick and a scale-bar legend. Deliberately blind to altitude/z, same as `galaxymap.py`'s Quadrant view -- the course panel's own Altitude figure already covers that axis. Not web-accessible. |
| `../src/html/static/style.css` | Shared stylesheet (CSS custom properties, light/dark via `prefers-color-scheme`, card-style panels), served directly (not through CGI). |
| `../src/html/static/vendor/` | Vendored third-party JS -- currently just `three.module.min.js` (three.js, bundled+minified from the `three` npm package), used by `sectormap.js`/`systemmap.js`. Vendored rather than loaded from a CDN so `lib/page.py`'s `Content-Security-Policy: default-src 'self'` needs no exception; see this directory's own `THIRD_PARTY_NOTICES.txt` for the license and how to rebuild it from a newer release. Served directly, same as `style.css`. |
| `../src/html/static/sectormap.js` | Renders `lib/starmap.py`'s `#starmap-data` JSON as a real WebGL scene (three.js, vendored at `static/vendor/` -- see that directory's `THIRD_PARTY_NOTICES.txt`): a perspective camera, GPU-billboarded sprites (always face the camera by construction, no manual per-frame counter-rotation needed) for stars/nebulae/asteroid fields/black holes/neutron stars, and a wireframe outline for the sector's wedge or fallback cube. Hand-rolled drag-to-rotate/scroll-or-button-to-zoom/click-for-info (a raycast against the sprites, not DOM hit-testing) plus a hidden, keyboard/screen-reader-focusable button list (a canvas has no focusable children of its own the way the old per-star `<div role="button">`s were) so every star/cloud stays reachable without a mouse. Clicking (or activating a fallback-list button) fills the info side panel and, for a "View system"/"View phenomenon" link, navigates via `data-nav-target`/`data-nav-params` (`static/navform.js`, whose global click handler also intercepts this dynamically-built link, same as any other in `html/`) instead of a plain `href`, so `db`/an id never shows up in the address bar. Served directly, same as `style.css`. |
| `../src/html/static/systemmap.js` | Toggles which System Map `<svg>` scene is visible (the whole-system view, or one per planet's own moon system) and fills the info side panel -- including "Atmosphere"/"Surface composition"/"Surface temperature"/"Life Chemistry" fields -- from a clicked marker's `data-*` attributes. A planet/moon selection also redraws `#sysmap-preview`: a small rotating shaded three.js sphere (the same vendored build `sectormap.js` uses) colored by the body's own `data-color`, banded with a tilted ring for a gas giant (`data-bodytype`), and wrapped in a fresnel-glow atmosphere shell (tinted by `data-surfacetemp`) when `data-hasatmosphere` is set -- the one genuinely 3D element on this otherwise flat-SVG page, an appearance preview only (no position data, and no cross-page link of its own to navigate). Served directly, same as `style.css`. |
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
   enables Apache's CGI module, and sets `../src/html/` ownership for
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
PLANETGEN_API_BASE_URL=http://127.0.0.1:5000/api QUERY_STRING="db=planetgen" python3 html/browse.py
```
