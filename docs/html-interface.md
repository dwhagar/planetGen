# planetGen Web Interface

A small web interface (in [`../src/html/`](../src/html/)) for
browsing the MySQL databases described in
[`database-schema.md`](database-schema.md) -- pick a database (schema),
drill into its sectors and star systems, and view (or copy) the
rendered wikitext/Markdown page saved for each one.

This is a read-only browser, and the frontend half of `TODO.md`'s Phase 5
web application: every page here is a thin server-rendered client of the
Flask API in [`../src/api/`](api.md), never touching the database
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
beyond the standard library) and renders it as HTML. `nltk` (from
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
| `../src/html/index.py` | Lists every matching MySQL schema on the configured server (`GET /api/databases`). |
| `../src/html/browse.py` | A chosen database's sectors and standalone systems. |
| `../src/html/sector.py` | One sector's systems (name, quadrant, star type), plus an interactive 3D "Sector Map" of the sector's cube -- drag to rotate, scroll (or the +/- buttons) to zoom, one dot per placed system (two, overlapping, for a binary) sized by star radius and colored by spectral type/brightness; clicking a dot fills an info panel (star type, temperature, quadrant, location) with a link to `system.py`, via `lib/starmap.py` + `static/sectormap.js`. |
| `../src/html/system.py` | One system's stars/planets/moons/belts, plus its description -- rendered as HTML from `markdown_content` by default (`?view=rendered`), with the original raw wikitext/Markdown source (`?view=source&format=...`) still available for copy-pasting into a wiki. Links to `nav.py` (a "Navigate from here" button) whenever the system is assigned to a sector. |
| `../src/html/nav.py` | Course, distance, and optimal route between two systems -- calls `GET /api/nav` (`queryDb.nav_between`; see [`api.md`](api.md#nav) for the availability rules and course convention). Without a destination yet, shows a dropdown of the origin's own sector-mates plus, when that sector has a galaxy placement, a plain destination-system-id field for a cross-sector route (there's no bounded way to offer every galaxy-placed system in a dropdown). With a destination, shows the direct course (distance/azimuth/altitude/warp-1-3-6-9 travel times), the NAV Map (a flat, top-down SVG plot of the origin/destination/route hops -- `html/lib/navmap.py`), and the optimal route via adjacent systems, each hop linking to `system.py`. |
| `../src/html/search.py` | Faceted search: click-to-filter tag buttons for object type, star spectral/luminosity class, and planet class/body type/supported life chemistry -- with a separate, identically-shaped set of tags for moons, since planets and moons live in their own tables (schema v2) and a "Class D" tag only ever means one or the other -- built only from values actually present in the chosen database (`GET /api/search`, which owns the query logic; this page just renders it). Plus a name search (with HTML5 `<datalist>` autocomplete, no JavaScript) across sectors, star systems, stars, planets, and moons. Asteroid belts have no name of their own, so they're reachable only via the "Asteroid Belt" object-type tag. |
| `../src/html/lib/apiclient.py` | The `GET /api/...` HTTP client (stdlib `urllib` only) every page above calls instead of querying MySQL directly -- one typed wrapper function per endpoint, plus `NotFoundError`/`ApiError` (`lib/page.py`'s `run` turns these into a 404/502 page). `PLANETGEN_API_BASE_URL` (default `http://127.0.0.1/api`) is where it looks for the API. Not web-accessible. |
| `../src/html/lib/fmt.py` | HTML-escaping and small formatting helpers (`esc`, `linkify_location`, `format_density`) with nothing to do with fetching data -- what's left of the old `dbutil.py` once its database-access functions moved into `apiclient.py`/the API itself. Not web-accessible. |
| `../src/html/lib/page.py` | Shared CGI response/HTML-shell helpers. Not web-accessible. |
| `../src/html/lib/mdconvert.py` | A small, purpose-built Markdown-to-HTML converter for the narrow Markdown subset `StarSystem.__str__` actually generates (headers, pipe tables, paragraphs, `<sup>` exponents) -- not a general-purpose parser. Not web-accessible. |
| `../src/html/lib/starmap.py` | Builds `sector.py`'s "Sector Map" panel: a real CSS 3D scene (`perspective`-free/orthographic, `transform-style: preserve-3d`) of the sector's cube -- a 6-bordered-`<div>`-face wireframe plus one billboarded `<div>` per star, positioned via plain layout (`left`/`top`) for x/y and `transform: translateZ()` for z. Star dot radius comes from `radius_km` (square-root scaled against the Sun); dot color from `star_type`'s spectral letter (`SPECTRAL_CLASS_COLORS`) shaded by `luminosity_w` and nudged by where `temperature_k` falls in that spectral class's range -- not a raw Kelvin-to-RGB blackbody map, so "White Giant" reads white and "Blue Giant" reads blue regardless of temperature. Not web-accessible. |
| `../src/html/lib/navmap.py` | Builds `nav.py`'s "NAV Map" panel: a flat, static, top-down SVG plot of the galactic X-Y plane -- origin and destination as labeled points, a dashed line for the direct course, and (when one was found) a solid polyline through the optimal route's intermediate hops. Auto-scaled to whatever points it's given (no fixed sector size to normalize against), with one uniform light-years-per-pixel ratio on both axes so azimuth angles aren't visually distorted, plus a `+X` compass tick and a scale-bar legend. Deliberately blind to altitude/z, same as `galaxymap.py`'s Quadrant view -- the course panel's own Altitude figure already covers that axis. Not web-accessible. |
| `../src/html/static/style.css` | Shared stylesheet (CSS custom properties, light/dark via `prefers-color-scheme`, card-style panels), served directly (not through CGI). |
| `../src/html/static/sectormap.js` | Drag-to-rotate, scroll/button-to-zoom, and click/keyboard-for-info behavior for the 3D sector map -- tracks two rotation angles and a zoom factor fed to `#starmap-scene`'s CSS transform (the browser's own compositor does the actual 3D projection and occlusion), counter-rotates each star dot every frame so it stays facing the camera instead of going edge-on, and resolves which dot was clicked by geometry (`getBoundingClientRect`) rather than relying on native hit-testing through the rotated 3D stack. Clicking (or Enter/Space on a focused dot) fills the info side panel from the dot's `data-*` attributes (never `innerHTML`) instead of navigating straight to `system.py`, so activating a dot shows details first. Served directly, same as `style.css`. |

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

Every page here needs the planetGen API (`../src/api/`, see
[`api.md`](api.md)) reachable to work at all now -- set
`PLANETGEN_API_BASE_URL` (default `http://127.0.0.1/api`, i.e. the same
host this CGI script itself runs on) if it's deployed somewhere else,
e.g. `http://127.0.0.1:5000/api` for `python src/wsgi.py`'s own local dev
server. The *database* server/account, and which schemas the picker
(`index.py`/`?db=`) offers, are entirely the API's own configuration now
(`PLANETGEN_MYSQL_HOST`/`_PORT`/`_USER`/`_PASSWORD`/`_DATABASE_PREFIX`,
`stellarObjects._db.MySQLConfig`/`list_databases`) -- see
[`api.md`](api.md#running-locally) for those; `html/` itself no longer
reads any `PLANETGEN_MYSQL_*` variable at all.

Separately from all of the above, a `webconfig.json`
file at the repo root (a sibling of `../src/html/`, not a file inside `../src/html/`
itself) holds site-level settings such as `site_name` and `base_url`,
edited once per deployment rather than passed through the vhost config --
see [`webconfig.md`](webconfig.md) for the full field list
and how it relates to `PLANETGEN_DEBUG`.

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

Every script is a normal CGI program: it reads `QUERY_STRING` from the
environment and writes an HTTP response (status + headers + body) to
stdout. That makes them runnable directly for a quick smoke test without
standing up Apache at all -- start the API separately first (see
[`api.md`](api.md#running-locally)), then, from the repo root:

```bash
PLANETGEN_API_BASE_URL=http://127.0.0.1:5000/api QUERY_STRING="db=planetgen" python3 html/browse.py
```
