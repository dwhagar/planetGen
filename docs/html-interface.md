# planetGen Web Interface

A small web interface (in [`../src/html/`](../src/html/)) for
browsing the MySQL databases described in
[`database-schema.md`](database-schema.md) -- pick a database (schema),
drill into its sectors and star systems, and view (or copy) the
rendered wikitext/Markdown page saved for each one.

This is a read-only browser, not the Phase 5 web application described in
[`TODO.md`](TODO.md) -- there's no backend framework and no writing to
the database, though [`../src/html/search.py`](#how-it-works) does provide a
faceted/name search. It exists so a generated galaxy can be looked at from
a browser today, on nothing more than Apache2 and a system Python 3.

## How it works

Every page here is a plain [CGI](https://en.wikipedia.org/wiki/Common_Gateway_Interface)
script (`#!/usr/bin/env python3`, executed fresh by Apache on every
request) -- no Flask/FastAPI, just the `planetGen` package's own runtime
dependencies (`pymysql`/`DBUtils` for the MySQL connection, `nltk` from
generation). That keeps deployment to "copy the files, `pip install .`,
set permissions, enable a vhost" with no separate web framework to
install or run.

| File | Purpose |
|---|---|
| `../src/html/index.py` | Lists every matching MySQL schema on the configured server. |
| `../src/html/browse.py` | A chosen database's sectors and standalone systems. |
| `../src/html/sector.py` | One sector's systems (name, quadrant, star type), plus an interactive 3D "Sector Map" of the sector's cube -- drag to rotate, scroll (or the +/- buttons) to zoom, one dot per placed system (two, overlapping, for a binary) sized by star radius and colored by spectral type/brightness; clicking a dot fills an info panel (star type, temperature, quadrant, location) with a link to `system.py`, via `lib/starmap.py` + `static/sectormap.js`. |
| `../src/html/system.py` | One system's stars/planets/moons/belts, plus its description -- rendered as HTML from `markdown_content` by default (`?view=rendered`), with the original raw wikitext/Markdown source (`?view=source&format=...`) still available for copy-pasting into a wiki. Links to `nav.py` (a "Navigate from here" button) whenever the system is assigned to a sector. |
| `../src/html/nav.py` | Course, distance, and optimal route between two systems -- reuses `queryDb.nav_between` (the same function behind `GET /api/nav`; see [`api.md`](api.md#nav) for the availability rules and course convention). Without a destination yet, shows a dropdown of the origin's own sector-mates plus, when that sector has a galaxy placement, a plain destination-system-id field for a cross-sector route (there's no bounded way to offer every galaxy-placed system in a dropdown). With a destination, shows the direct course (distance/azimuth/altitude/warp-1-3-6-9 travel times), the NAV Map (a flat, top-down SVG plot of the origin/destination/route hops -- `html/lib/navmap.py`), and the optimal route via adjacent systems, each hop linking to `system.py`. |
| `../src/html/search.py` | Faceted search: click-to-filter tag buttons for object type, star spectral/luminosity class, and planet class/body type/supported life chemistry -- with a separate, identically-shaped set of tags for moons, since planets and moons live in their own tables (schema v2) and a "Class D" tag only ever means one or the other -- built only from values actually present in the chosen database. Plus a name search (with HTML5 `<datalist>` autocomplete, no JavaScript) across sectors, star systems, stars, planets, and moons. Asteroid belts have no name of their own, so they're reachable only via the "Asteroid Belt" object-type tag. |
| `../src/html/lib/dbutil.py` | Read-only database access and HTML-escaping helpers. Not web-accessible -- see the Apache config note below. |
| `../src/html/lib/page.py` | Shared CGI response/HTML-shell helpers. Not web-accessible. |
| `../src/html/lib/mdconvert.py` | A small, purpose-built Markdown-to-HTML converter for the narrow Markdown subset `StarSystem.__str__` actually generates (headers, pipe tables, paragraphs, `<sup>` exponents) -- not a general-purpose parser. Not web-accessible. |
| `../src/html/lib/starmap.py` | Builds `sector.py`'s "Sector Map" panel: a real CSS 3D scene (`perspective`-free/orthographic, `transform-style: preserve-3d`) of the sector's cube -- a 6-bordered-`<div>`-face wireframe plus one billboarded `<div>` per star, positioned via plain layout (`left`/`top`) for x/y and `transform: translateZ()` for z. Star dot radius comes from `radius_km` (square-root scaled against the Sun); dot color from `star_type`'s spectral letter (`SPECTRAL_CLASS_COLORS`) shaded by `luminosity_w` and nudged by where `temperature_k` falls in that spectral class's range -- not a raw Kelvin-to-RGB blackbody map, so "White Giant" reads white and "Blue Giant" reads blue regardless of temperature. Not web-accessible. |
| `../src/html/lib/navmap.py` | Builds `nav.py`'s "NAV Map" panel: a flat, static, top-down SVG plot of the galactic X-Y plane -- origin and destination as labeled points, a dashed line for the direct course, and (when one was found) a solid polyline through the optimal route's intermediate hops. Auto-scaled to whatever points it's given (no fixed sector size to normalize against), with one uniform light-years-per-pixel ratio on both axes so azimuth angles aren't visually distorted, plus a `+X` compass tick and a scale-bar legend. Deliberately blind to altitude/z, same as `galaxymap.py`'s Quadrant view -- the course panel's own Altitude figure already covers that axis. Not web-accessible. |
| `../src/html/static/style.css` | Shared stylesheet (CSS custom properties, light/dark via `prefers-color-scheme`, card-style panels), served directly (not through CGI). |
| `../src/html/static/sectormap.js` | Drag-to-rotate, scroll/button-to-zoom, and click/keyboard-for-info behavior for the 3D sector map -- tracks two rotation angles and a zoom factor fed to `#starmap-scene`'s CSS transform (the browser's own compositor does the actual 3D projection and occlusion), counter-rotates each star dot every frame so it stays facing the camera instead of going edge-on, and resolves which dot was clicked by geometry (`getBoundingClientRect`) rather than relying on native hit-testing through the rotated 3D stack. Clicking (or Enter/Space on a focused dot) fills the info side panel from the dot's `data-*` attributes (never `innerHTML`) instead of navigating straight to `system.py`, so activating a dot shows details first. Served directly, same as `style.css`. |

This project recommends (but doesn't enforce in code) pointing this web
interface at a MySQL account with `SELECT`-only grants, so these scripts
can't write to a database even if a query were buggy -- see
`queryDb.py`'s module docstring for the same convention. Database and
system names pulled from generated data are HTML-escaped before being
placed in a page; a requested `?db=` schema name is validated against the
actual, prefix-filtered schema listing (exact match only), which is what
prevents it from being used to select a schema this deployment never
meant to expose. `mdconvert` escapes every block in full before emitting
any markup, then narrowly re-enables only the one legitimate raw-HTML
pattern generated content ever contains (`<sup>...</sup>`) -- so a
mischievous `--name`/`--star-type` value can't inject live HTML into a
rendered page.

## Locating the database

By default, each script connects to the same MySQL server/account every
other tool in this project defaults to (`stellarObjects._db.MySQLConfig`,
`127.0.0.1:3306`, user/database `planetgen`), and the picker
(`index.py`/`?db=`) offers every schema on that server whose name starts
with `planetgen` (see `../src/html/lib/dbutil.py`). Set the
`PLANETGEN_MYSQL_HOST`/`PLANETGEN_MYSQL_PORT`/`PLANETGEN_MYSQL_USER`/
`PLANETGEN_MYSQL_PASSWORD`/`PLANETGEN_MYSQL_DATABASE_PREFIX` environment
variables (e.g. via `SetEnv` in the Apache vhost) to point somewhere else
or restrict/widen which schemas are offered.

Separately from these Apache-set environment variables, a `webconfig.json`
file at the repo root (a sibling of `../src/html/`, not a file inside `../src/html/`
itself) holds site-level settings such as `site_name` and `base_url`,
edited once per deployment rather than passed through the vhost config --
see [`webconfig.md`](webconfig.md) for the full field list
and how it relates to `PLANETGEN_MYSQL_*`/`PLANETGEN_DEBUG`.

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
standing up Apache at all, from the repo root:

```bash
PLANETGEN_MYSQL_HOST=127.0.0.1 QUERY_STRING="db=planetgen" python3 html/browse.py
```
