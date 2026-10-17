{% raw %}
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
Flask API in [`../src/planetgen/api/`](api.md), never touching the database
directly itself, though the search page (`/search`, see "Flask pages"
below) does provide a faceted/name search (built on `GET /api/search`, same as every
other page). It exists so a generated galaxy can be looked at from a
browser today, on nothing more than a web server (Apache2 with mod_wsgi,
or nginx or Caddy in front of gunicorn), Python 3 and a
MySQL server (see "Locating the database" and "Deploying" below).

## Wording: "uncharted"

Boss (2026-10-10 07:26Z): anything not yet generated is called **uncharted**
in text shown to the user (labels, tooltips, map legends, empty states,
messages, API text meant for people), never "unbuilt", "not generated" or
the like. The exception is the Generate system and the admin panels and
controls, which keep the plain technical words ("generated", "not
generated", "queued", and so on). New items and pages follow this rule;
The sweep of the existing wording is done (PR #1049).

## Habitability chips

Each rocky planet's and moon's description also carries the full PHI-6 explanation (UX.91, `physics/habitability_explain.py`): the score and equipment, then one paragraph per colour factor with its colour, this world's input values and the input that sets it.

Planet and moon rows carry habitability chips (GEN.89, design in
`docs/design/habitability-index.md`): an equipment label (Ideal, formerly Shirtsleeve;
Mask, Mask and scrubber, Environment suit, Full life support), a colour chip
that comes from the six PHI-6 domains (pressure, temperature, air, water,
radiation, suitability; Blue, Green, Yellow, Red) under the gates of GEN.207, and the Habitable,
Habitable moon and Inhabited chips. Planets around pulsars, neutron stars
and black holes follow a host rule. The Equipment search has a facet for
the labels. The explanation page is /classes/habitability under Classes (UX.90, PR
#1065), linked from the chips and the Equipment filter; planet and moon
rows carry the per-domain colours and the stored microbial, complex-life
and human-operability scores. UX.91 added the full per-factor explanation
to the planet and moon description.

## How it works

Every page is served by the same Flask app as the JSON API
(`../src/html/wsgi.py`; the pages live in `../src/planetgen/web/`, see "Flask
pages" below). A page fetches its data through `../src/planetgen/web/lib/apiclient.py`,
which inside the app dispatches straight through the API's own routes
(no HTTP round trip), and renders it with a Jinja2 template. Links are
plain GET `<a href>`s, so every page is bookmarkable. `nltk` (from
generation) is still needed for a handful of display-formatting constants
a few pages import from the `planetgen` generator modules (spectral-class
colors, planet-class descriptions, unit formatters), some of which pull
in `planetgen.names.wordlists` and its NLTK corpus dependency. See
[`api.md`](api.md#deploying) for the process
itself, which needs `pymysql` and a database account.

| File | Purpose |
|---|---|
| `../src/html/wsgi.py` | The WSGI entry point (loaded by mod_wsgi or gunicorn): the Flask app serving the API under `/api` and every page. |
| `../src/planetgen/web/` | The pages: routes, templates and helpers (see "Flask pages" below). `web/old_urls.py` answers the old `/<name>.py` CGI URLs with a 301 to the page that replaced them. |
| `../src/planetgen/web/lib/pagination.py` | The site's one pager, used under every paged table (Browse's two tables, Phenomena, a sector's Contents table, a Galaxy Map Quadrant's sector list, each Search result panel, the admin API key list and the admin stats page's duplicate-names list): a "Showing X-Y of Z" summary, then First/Prev, numbered pages and Next/Last, 50 rows a page. Each table has its own page parameter (e.g. `sectors_page`), a plain GET link, and changing a Search filter starts its results back at page 1. |
| `../src/planetgen/web/lib/apiclient.py` | The API client every page calls instead of querying MySQL directly -- one typed wrapper function per read endpoint, plus `auth_*` wrappers (the admin pages) supporting POST/DELETE, a request body, and `Cookie`/`Set-Cookie` relay, and `NotFoundError`/`ApiError` (`web/errors.py` turns these into a 404/502 page; `ApiError.status_code` lets `auth_me`/the admin pages branch on a 401 without string-matching). Inside the Flask app it runs in-process (`web/transport.py`); elsewhere it uses HTTP to `PLANETGEN_API_BASE_URL` (default `http://127.0.0.1/api`). Not web-accessible. |
| `../src/planetgen/web/lib/fmt.py` | HTML-escaping and small formatting helpers (`esc`, `linkify_location`/`nearest_neighbors_location` -- which take a `system_url(system_id)` hook for the neighbour links -- `format_density`, and `static_url`, the versioned `static/` URL every page uses) with nothing to do with fetching data. Not web-accessible. |
| `../src/planetgen/web/lib/mdrender.py` | A small, purpose-built Markdown-to-HTML converter for the narrow Markdown subset `StarSystem.__str__` actually generates (headers, pipe tables, paragraphs, `<sup>` exponents) -- not a general-purpose parser. Not web-accessible. |
| `../src/planetgen/web/maps/galaxymap.py` | Quadrant/Zone classification (`sector_quadrant`/`sector_zone`/`zone_bounds_ly`) shared by the `/galaxy` page's data tables and the sector/sectors pages' "Quadrant N" links back into it -- what's left after this module's former flat-SVG map rendering was superseded by a real 3D scene (`planetgen/web/maps/galaxymap3d.py`). Not web-accessible. |
| `../src/planetgen/web/lib/modulegraph.py` | The modules a map script imports, directly or not (`module_graph('galaxymap3d.js')`), which `galaxy.html` and `sector.html` turn into `<link rel="modulepreload">` hints (MAP.161) so the browser fetches the tree in parallel. |
| `../src/planetgen/web/lib/tilecache.py` | The web layer's on-disk cache (written by the WSGI process serving `/galaxy` and `/galaxy/tiles`) of 3D Galaxy Map tiles: `fetch_tiles` serves each requested tile from `<cache dir>/<db>/<generation>/` when present and asks `GET /api/galaxy/tiles` only for the rest. Every 60 s at most it asks `GET /api/galaxy/changes` which tiles changed (from the rows' `modified_at`) and deletes only those, or starts a new generation when the answer is `full`. The changed keys are passed on to the browser's cache too. Size-capped (`tile_cache.max_mb`, oldest files pruned first); location from `PLANETGEN_TILE_CACHE_DIR`/`tile_cache.dir` (see `config.md`). Fails open: an unwritable directory or bad file just means an API call. Not web-accessible. |
| `../src/planetgen/web/lib/pagecache.py` | The pages' cache of the API's public GET answers (`apiclient._request` asks it first; `web.init_app` gives each app its own copy in memory and the answers are shared through Redis, so the warm-up can fill them for every worker). An answer is stored with the data tags it read and served only while those tags' counters (`db/dataversion.py`) are unchanged; an API write, a command-line run and the galaxy's content stamp (checked every `page_cache.stamp_seconds`) bump the tags they touch; capped by age and size (see `config.md`, [cache-warm-up.md](design/cache-warm-up.md)). Not web-accessible. |
| `../src/planetgen/web/lib/slowpages.py` | The cost rule for slow page requests and what a page sees of them (PERF.105): a route that expects a slow answer runs inside `slowpages.allow_pending()`; a cache miss the rule calls slow raises `slowpages.Pending` carrying the job (`queue/pagejobs.py`) that builds the answer into the page cache, which the page polls at `GET /api/page-jobs/<id>`. Redis down builds it in the request and `take_notice()` says so (see [slow-page-answers.md](design/slow-page-answers.md)). Not web-accessible. |
| `../src/planetgen/web/lib/pagetags.py` | What each cached answer reads (`read_tags`) and what each API write changes (`change_for_write`), the two tables behind the page cache's invalidation. Not web-accessible. |
| `../src/planetgen/web/maps/galaxymap3d.py` | Builds the `/galaxy` page's "Galaxy Map" panel: the canvas/controls/info-panel markup, plus the one starting JSON payload (zoom-range numbers, tile settings, and the zoomed-all-the-way-out view's cube tiles, which `initial_tile_request` picks the same way the client does for every later view) `static/galaxymap3d.js` reads on first paint -- every later payload, as the camera moves, is fetched by that script directly and never passes through this module. `view_radius_bounds` (the zoom floor/ceiling: a couple of sector-widths up to this galaxy's own real outer edge) is a pure function the `/galaxy` view also calls directly, before this module's own panel-rendering function runs. The panel's JSON names the tile URL (`fetchPath`, `/galaxy/tiles`) and a sector-page URL template (`sectorUrl`, from `page_url`) for the info panel's real "View sector" `<a href>`; it carries no database name except `storageKey`, which only namespaces the browser's `localStorage`. Not web-accessible. |
| `../src/planetgen/web/maps/systemmap.py` | Builds the system page's "System Map" panel: a fixed-size, square true-position plot (real angle from `position_x/y_km`, a shared log-radial scale from anchor -- star, or barycenter for a merged/'close' binary pair -- to body) with one colored/sized marker per body, plus a full ring (not a directional band) for each asteroid belt -- distinct from `starmap.py`'s draggable 3D cube, since this only ever needs a flat top-down projection. A once-only pairwise-repulsion pass (`_relax_markers`) nudges apart any two markers real placement happened to put too close together, real position first, decluttering only where needed. Clicking a planet with moons swaps to a "zoom into this planet's moons" scene (`static/systemmap.js` toggles which `<svg>` scene is visible); clicking anything else fills the info side panel, same click-for-info pattern as `starmap.py`/`sectormap.js`. Adds a small green badge to any body with `life_chemical` set. Also hands over each planet/moon marker's resolved class color (`class_color`) plus its `atmosphere`/`composition`/`surface_temperature_k` (whether it has a real atmosphere at all, not just the description text) as `data-*`, plus a star's own spectral-type color (`star_color`, same `data-color` attribute), for `static/systemmap.js`'s own per-marker live sphere rendering. Not web-accessible. |
| `../src/planetgen/web/lib/classref.py` | The class reference catalog behind `/classes`: every class type built from the generator's own tables (`program_constants`, `physical_constants`, `stellarEvolution.YERKES_CLASS_NAMES`, `nebulaData`, `cometData`), never hand-copied, once per process (`catalog()`, called by `web.init_app` at startup) so the pages always show the values the generator uses. Also `class_url_parts` (is this a real class?) and `star_type_classes` ("G2V Yellow Main Sequence Star" -> `("G", "V")`). Not web-accessible. |
| `../src/planetgen/web/lib/systempage.py` | The system page's Python-built HTML (moved from the old `system.py`): the expandable body list (`system_list_html`) and the Stars/Planets & Moons/Asteroid Belts/Comets tables (`stars_html`, `bodies_html`). Escapes every database value itself; `web/system_pages.py` passes the result through `trusted_html`. Not web-accessible. |
| `../src/planetgen/web/lib/tabledisplay.py` | Computes the same "Star Data"/"Planet Data" display strings once baked into the database's now-removed `table_*`/`binary_table_*` columns, but on demand from the raw numeric columns the system page already has -- reuses `planetgen.util.format`'s formatters directly. Not web-accessible. |
| `../src/planetgen/web/maps/starmap.py` | Builds the Sector Map's scene data (`map_scene_data(link_url, ...)`, where `link_url` is `web.helpers.page_url`; every entry carries a plain `href`), which the sector page answers at `/sector/<id>/scene`: computes every position/size/color/label the map needs (a wedge or fallback-cube outline, one entry per star -- a small core (`_star_dot_radius`: log-scaled from `radius_km` but kept between 1.5 and 6 scene units, so even a supergiant is a point) inside a glow shell whose size, strength and fade come from `luminosity_w` (`_star_glow`: a supergiant gets a big soft halo, a white dwarf almost none, matching how the Galaxy Map draws bright stars; `sectorscene.js` clicks a star through an unseen sphere at least 4 units across), color from `star_type`'s spectral letter (`SPECTRAL_CLASS_COLORS`) shaded by `luminosity_w` and nudged by where `temperature_k` falls in that spectral class's range, so "White Giant" reads white and "Blue Giant" reads blue regardless of temperature -- and one per nearby standalone phenomenon, sized by `radius_ly` (always 0 for the point-like types) and positioned directly in the galaxy frame, no rotation needed unlike a star system's sector-local position -- see `queryDb.phenomena_near_sector`) and returns it as a dict; `static/sectorscene.js` is what actually renders it, this module builds no HTML of its own. Not web-accessible. |
| `../src/planetgen/web/maps/navmap.py` | Builds the NAV page's "NAV Map" panel (`render_nav_map_panel(link_url, ...)`; each point is an SVG `<a href>`): a flat, static, top-down SVG plot of the galactic X-Y plane -- origin and destination as labeled points, a dashed line for the direct course, and (when one was found) a solid polyline through the optimal route's intermediate hops. Auto-scaled to whatever points it's given (no fixed sector size to normalize against), with one uniform light-years-per-pixel ratio on both axes so bearings aren't visually distorted, plus a compass arrow along the origin's bearing 000 (toward the frame's center) and a scale-bar legend. Deliberately blind to altitude/z, same as the flat SVG phenomenon Diagram panel (`planetgen/web/maps/phenomenonmap.py`) -- the course's mark already covers that axis. Not web-accessible. |
| `../src/planetgen/web/maps/phenomenonmap.py` | Builds the phenomenon page's flat, zoomable SVG diagram of a nebula's or supernova remnant's real extent, drawn in astronomical units against an AU-scale yardstick (`static/phenomenonmap.js` adds the zoom and pan). Not web-accessible. |
| `../src/planetgen/web/maps/phenomenonrender.py` | Builds the phenomenon page's "View" panel for a neutron star, black hole, quasar, rogue planet or interstellar comet: the numbers `static/phenomenonrender.js` draws with three.js, plus a static SVG used when WebGL can't start. Not web-accessible. |
| `../src/planetgen/web/warmup.py` | The site's caches built ahead of the first visit (MAP.134, PERF.104): `opening_request` names the tiles the first frame of `/galaxy` asks for (the page uses it too), `warm_opening_view` fetches them and the galaxy stage into the tile cache, and `warm_family` / `run_family_job` visit the pages of one family (`FAMILIES`: galaxy-map, lists, sectors, systems, phenomena, search) through the site's own app so the answers land in the shared page cache; `run_family_job` is the RQ job that shows on the admin runs page. `python -m planetgen.cli.warm_map` runs the families; `update.sh` starts it in the background after every update (at most 30 minutes on Linux, never waited for), and it is the command to run after clearing a cache by hand. Not web-accessible. |
| `../src/planetgen/web/lib/datatable.py`, `../src/planetgen/web/tables.py` | The site's shared data table (UX.41): a page module describes its table once as a `Table`, `register`s it and renders it; the first page is rendered on the server (sort links, a filter form, the shared pager), and `/table/<name>` answers the next pages as JSON for `static/datatable.js`. |
| `../src/planetgen/web/lib/tilepack.py` | The packed form of a Galaxy Map tile (MAP.159) the tile cache stores and `/galaxy/tiles` sends: stars in columns (about 19 bytes each, against about 150 as JSON), the clouds, points and `incomplete` list as JSON. `static/tilepack.js` reads it. |
| `../src/planetgen/web/lib/singleflight.py` | One build of each Galaxy Map tile or stage at a time (PERF.38): a request claims the tiles it is about to fetch with a Redis `SET NX EX`, waits for ones another process claimed, and builds only those whose claimant never delivered. With no Redis it fails open and everything is built. |
| `../src/planetgen/web/lib/estimatecache.py` | The Generate page's size-and-time estimate (`planetgen galaxy ... --estimate-only`, up to two minutes) kept in Redis for five minutes under the job kind, a digest of its command line and the database (PERF.109). It is served only while the generator epoch, the galaxy's data version (`db/dataversion.py`) and the speed statistics (`generation/stats.py` `fingerprint`) are as they were, so a fill, an edit, a new epoch or newly recorded speeds make the next estimate a new one. Two admins asking for the same estimate at once run it once (`singleflight`); with no Redis it just runs. |
| `../src/planetgen/web/lib/habitability_ref.py`, `../src/planetgen/web/lib/facility_display.py` | The text of the habitability levels page (`/classes/habitability`, UX.90), and how a facility's type, affiliation, ratings and field values read (POP.10: rating chips carry a shape, the dimension's name and its word, so colour is never the only signal). |
| `../src/planetgen/web/nearby_page.py`, `population_pages.py`, `system_species.py`, `facility_type_pages.py`, `queue_page.py`, `settings_page.py`, `backups_page.py`, `images_page.py`, `seo.py`, `seo_pages.py`, `assets.py`, `maintenance_page.py`, `maintenance_gate.py`, `job_runner.py` | The What's nearby page (NAV.44); the species and polity pages (POP.1 to POP.4); the system page's "Species that could live here" panel (POP.12); the admin Facility types page (POP.10); the admin work queue page (ADM.10); the admin settings page (ADM.84); the backup archives screen and the settings export and import (ADM.85); the Images screen, the page `<head>`, the indexing matrix, robots.txt, the sitemaps, the small site files and the icon and picture store (ADM.44); the `/maintenance` page and the `before_request` hook that sends everyone but an admin to it while maintenance mode is on (UX.112); and the RQ job that runs one Generate page job's steps in order, writing `output.log` and `state.json` in the job's directory (PERF.24). |
| `../src/planetgen/web/lib/privatedir.py` | The private fallback directory (mode 0700, in the system temp directory) the tile cache and the Generate jobs use when their configured directory can't be created. A directory found there that isn't a real directory owned by this user, or that others can write to, is refused rather than reused. Not web-accessible. |
| `../src/html/static/style.css` | Shared stylesheet (CSS custom properties, light/dark via `prefers-color-scheme` or an explicit `data-theme` on `<html>`, card-style panels, phone layout under 40rem), served directly by the web server. |
| `../src/html/static/theme.js` | Loaded on every page, blocking, before `style.css` (`web/templates/base.html`): applies the saved light/dark/system theme before the first paint and drives the header's theme button. See "The page shell" below. |
| `../src/html/static/favicon.svg` | The site icon (a small ringed planet), linked from every page's `<head>`. |
| `../src/html/static/vendor/` | Vendored third-party JS -- currently just `three.module.min.js` (three.js, bundled+minified from the `three` npm package), used by `sectorscene.js`/`systemmap.js`/`galaxymap3d.js`. Vendored rather than loaded from a CDN so the pages' `Content-Security-Policy: default-src 'self'` needs no exception; see this directory's own `THIRD_PARTY_NOTICES.txt` for the license and how to rebuild it from a newer release. Served directly, same as `style.css`. |
| `../src/html/static/sectorscene.js` | One sector's scene (MAP.66): `buildSectorScene(data, options)` turns `starmap.py`'s JSON into a `THREE.Group` (GPU-billboarded points of light for every star and light-giving phenomenon, textured spheres for the rest, translucent volumes for nebulae, neighbor markers, the cell outline and compass) with the picker layers (`mappick.js`), info-panel specs and tooltip text for what is in it. `options.origin`, `unit` and `flipY` place it in a map's world (the Galaxy Map's parsec frame); with none it is the Sector Map's own scene units. The sector page's map and the Galaxy Map's opened sector both draw it. Each kind of object in it (stars, nebulae, supernova remnants, asteroid fields, black holes, neutron stars, quasars, rogue planets, interstellar comets, neighboring sectors) can be hidden (`setKindHidden`, MAP.79): its points, bodies, rings and labels go, and a hidden kind can't be hovered or picked. Served directly, same as `style.css`. |
| `../src/html/static/galaxysector.js` | The Galaxy Map's last drill-down stage (MAP.66): `createSectorStage(host).open(id, bounds)` fetches a sector's scene (`GET /sector/<id>/scene`, `web/sector_page.py`: the same JSON the sector page embeds, plus `centerPc` and `halfEdgePc`), builds it with `sectorscene.js` at its place in the galaxy, adds its picker layers (with tooltip, hover ring and info-panel select) to the map's picker, and has the galaxy's own tile stars and phenomena left out of its cell; `close()` takes it all away. `galaxystageview.js` decides when a sector is open, flies the camera to it and keeps `?sector=<designation>&open=1` and the map's Back and Forward. Served directly, same as `style.css`. |
| `../src/html/static/generatebuttons.js` | The admin Generate buttons for a sector that isn't generated yet (it alone, its neighborhood, its column, or its whole shell after a confirm), shared by the sector page's map and the Galaxy Map (`galaxymap3d.js`, a sector cell or a neighbor). Each is a plain POST form to the Generate page with the CSRF token from the server's `generate` target (`web.helpers.generate_target`), which the pages only send a logged-in admin. Served directly, same as `style.css`. |
| `../src/html/static/starfade.js` | The Galaxy Map's fade of stars and other objects (MAP.153, MAP.155, MAP.152), free of the page. A star shows at the share of the tile that lists it: a tile gives way to its eight halves as the camera comes within about two edges of it, and a tile's share is what its ancestors leave it, so a tile and its halves always add up to one (`tileShare`); the finest tiles round the target fade by distance alone (`detailShare`). The page works the shares out as the camera moves and hands them to the star shader; black holes, neutron stars and quasars are shared the same way; a nebula sprite's opacity ramps between 1 and 4 px of radius (`cloudSizeShare`); a star arriving is held as a point of light until its zoom share rises (the shader, with `static/starpoint.js`). |
| `../src/html/static/starpoint.js` | MAP.208: how big a star, black hole, neutron star or quasar is drawn: a point of light 1.5 px wide until the camera is about a sector away (`NEAR_SECTOR_PC`, one constant the Galaxy Map and the sector scene share), growing smoothly over the last two zoom steps (`pointGrow`, `pointSizePx`, twinned in GLSL as `STAR_POINT_GLSL`), the far brightness (`pointBrightness`) and the cap that keeps a phenomenon no brighter than the brightest star (`phenomenonLuminositySol`). |
| `../src/html/static/galaxytilecut.js` | The tiles the Galaxy Map fetches (MAP.152): the distance cut, near tiles small and far ones big, inside the cone the camera looks along, and the set to fetch (the cuts either side of it too, for the fade). The same arithmetic as `planetgen.galaxy.viewport.distance_cut_tiles`; `tests/test_js_unit.py` holds the two equal. |
| `../src/html/static/starmagnitude.js` | The Galaxy Map's star visibility law (MAP.148), free of the page so the shader, the page's histogram and the tests share one definition: a star's opacity follows its apparent magnitude from the camera (`M + 5 log10(d / 10 pc)`) against a limiting magnitude with a 1.5 magnitude ramp, its brightness follows tone-mapped flux, and the limit is the larger of the HISTOGRAM limit (where the stars on screen, counted by the MAP.153 zoom opacity, reach 20,000, or 8,000 on a phone) and the CALIBRATION limit (where the dimmest stars the tiles carry still show whole at the focus distance), so it only dims what stands much farther than the thing looked at and never thins the same stars twice with the rank birth radius. `galaxymap3d.js` recounts at most every 150 ms when the view changes and the shader eases toward the new limit over about 0.3 s. Black holes, neutron stars and quasars are never dimmed away. The tile floors (`generated_star_floor_sol`) stay as the server's superset until tiles are chosen by distance. |
| `../src/html/static/nebulalook.js` | How a nebula is painted on the Galaxy Map, per theme (MAP.175): [colour, core opacity, edge opacity] for diffuse, emission, reflection, planetary, dark and unknown nebulae, one set for the light theme and one for the dark, each reaching 3:1 (`CONTRAST_TARGET`) against the map's background; the dark family is a dust colour on the dark theme. Also the WCAG contrast arithmetic the test uses, against the real colours in `style.css`. MAP.181 settles the final dust colour. |
| `../src/html/static/nearfield.js` | The Galaxy Map's near field (MAP.149), free of the page so the shaders, the picker and the tests share one definition: a depth fade (things nearer the camera than 0.15 of the focus distance are gone, from 0.5 on whole), a see-through tube (down to 0.08) round the line from the camera to the focus thinning what stands in front of it, the faint-region rules (`contextOpacity` 0.08-0.3, `prominence`, `limitingMagnitude`) and the polar `blockContains` test for the block the camera is inside. Everything is worked out in view space, and `NEAR_FIELD_GLSL` is the shaders' twin built from the same constants. `galaxymap3d.js` draws blocks and stars with it, and `galaxystageview.js` picks with it: a block whose share is under 0.35 at the hit point, or a back face, is skipped so the pick falls through to what is behind. |
| `../src/html/static/galaxymap3d.js` | Renders `planetgen/web/maps/galaxymap3d.py`'s `#galaxymap3d-data` JSON as a real WebGL scene (three.js, vendored at `static/vendor/`) -- this camera's own orbit target moves freely through the galaxy, so most of what it draws is fetched live from `/galaxy/tiles` (plain GET, no database in the URL) (debounced, on every camera move), one fixed cube of space at a time, rather than server-rendered once. Tiles are kept in memory and in the browser's IndexedDB (MAP.158: database `planetgen-tiles`, store `tiles`, keyed by database, generation and tile key, at most about 96 MB with the least recently used dropped; only the small stamp record stays in `localStorage` under `planetgen:tile:<db>:@stamp`), so panning back over seen space or reloading the page doesn't refetch them. Once the view's tiles are in and the camera has been still for 1.5 s, the tiles one zoom step in are fetched at low priority (not the way out, not with the browser's data saver on or a connection of 3G or slower, at most about 1.5 MB a page visit). The map is driven by the drill-down; the whole galaxy is a 3D disk with no grid lines drawn on it, and in every view drag turns it, right-drag or Shift-drag moves it, the wheel zooms toward the point under the cursor (that point stays under it, MAP.150) and a pinch zooms about the middle, both within limits (Re-center, in the Menu, brings it back; the how-to text is behind Map help in the Menu, UX.50) (`static/galaxystageview.js`, rules in `static/galaxystages.js`): the whole galaxy is a plain disk of blocks of 243 sectors (hovering one outlines it); a click on a block opens a window of nine blocks (the "size" `g`, 3 to 11 in the URL as `g=`) centred on the block clicked, one size down (243, 81, 27, 9, 3, then single sectors, MAP.151b), the sectors in a round region of the block lattice with the cut-off drawn as the window's own edge; Up recentres on the parent block and ends at the galaxy; there are no arc, slab or segment picks and no slab buttons (MAP.150, MAP.151d) (design doc `docs/design/galaxy-drilldown-navigation.md`, sections 4-5). Only the blocks in view are drawn (`static/galaxyprisms.js`, built by `static/galaxyblocks.js`), colored by predicted density; blocks holding generated sectors are amber and grow more solid with their filled share. Hovering a choice dims the others and outlines it along its blocks' own sides ; Back, Forward, Up and Whole galaxy buttons and the round Steps menu (every step so far, and Jump to newest) step through the stages, each of which has its own URL, beside Bookmarks and a Menu holding Re-center, Charted only, Territories and Map help. A block's panel gives its ring, layer and slot ranges, exact sector count and generated count; an unfilled sector's gives its address and designation, plus, for a logged-in admin, the Sector Map's Generate buttons (`static/generatebuttons.js`: this sector, its neighborhood, its column, or its whole shell after a confirm), which post to the Generate page. Nebulae and supernova remnants are drawn as soft translucent spheres their real size (each tile lists the ones reaching into it, `queryDb.galaxy_clouds_in_box`), hidden when only a pixel or two across and faded out as the camera nears them; clicking one while it is small enough to aim at shows its type, class and radius and links to its phenomenon page. Moves between stages fly along a smooth zoom-and-pan path (they jump with `prefers-reduced-motion`), and while idle the page fetches the tiles one zoom step out. A one-line scale readout gives a bar's length in sectors, pc and ly. Served directly, same as `style.css`. |
| `../src/html/static/systemmap.js` | Toggles which System Map `<svg>` scene is visible (the whole-system view, or one per planet's own moon system) and fills the info side panel -- including "Atmosphere"/"Surface composition"/"Surface temperature"/"Life Chemistry" fields -- from a clicked marker's `data-*` attributes. Every visible scene's star/planet/moon markers also get their own live-rendered 3D sphere on `#sysmap-spheres-canvas` (the same vendored three.js build `sectorscene.js` uses): one shared WebGL context, redrawn each frame via a scissored sub-viewport per marker (never one `<canvas>`/context per body -- browsers cap concurrent WebGL contexts), each sized and positioned to exactly cover that marker's own `<circle>` and colored by its `data-color`, banded with a tilted ring for a gas giant (`data-bodytype`), and wrapped in a fresnel-glow atmosphere shell (tinted by `data-surfacetemp`) when `data-hasatmosphere` is set -- the one genuinely 3D layer on this otherwise flat-SVG page, drawn behind the SVG so each marker's own stroke/label/life-badge still shows on top. An appearance layer only (no position data, and no cross-page link of its own to navigate) -- a marker whose sphere renders keeps its flat circle's fill transparent (`sysmap-sphere-active`) but never changes its actual plotted position. Served directly, same as `style.css`. |
| `../src/html/static/copycode.js` | The system page's Copy button: copies the generated Wikitext/Markdown out of the code box named by the button's `data-copy-target` (clipboard API, falling back to a selection copy on a plain-HTTP deployment). A separate file because the Content-Security-Policy allows no inline script; the box stays selectable by hand without it. Served directly, same as `style.css`. |
| `../src/html/static/bodyRendering.js` | Shared three.js building blocks (star granulation texture, fresnel glow shader, planet bands) used by `systemmap.js` and `sectorscene.js` to draw a body as a lit sphere. Served directly, same as `style.css`. |
| `../src/html/static/galaxyprisms.js`, `galaxyblocks.js` | The Galaxy Map's sector grid, block ladder and density shading on shells, bands and wedges (`galaxyprisms.js`, the page's copy of `galaxy/spherical.py`, `galaxy/drill.py` and `galaxy/window.py`), and the block scene for one view packed as GPU-ready arrays (`galaxyblocks.js`, run in a Web Worker). Neither imports three.js, so both run under plain node for the tests. Served directly, same as `style.css`. |
| `../src/html/static/galaxysystem.js` | The Galaxy Map's system stage (MAP.125): `createSystemStage(host).open(href, name, center, mode)` fetches a system's scene (`/system/<id>/scene`), builds it with `systemview3d.js` at the star's place in the sector (compressed or true scale), adds its bodies to the map's picker with tooltip, ring and info panel, and keeps them moving with the orbit clock; `close()` removes it. `galaxystageview.js` decides when it is open (`&system=<id>`). Served directly. |
| `../src/html/static/galaxystages.js`, `galaxystageview.js` | The Galaxy Map's drill-down from the galaxy to one sector (`design/galaxy-drilldown-navigation.md`): the stage rules (`galaxystages.js`) and the scene, camera moves, breadcrumb and stage URLs that draw them (`galaxystageview.js`), backed by `/galaxy/stage`; a generated sector is opened in place as the last stage (`galaxysector.js`, MAP.66). Served directly, same as `style.css`. |
| `../src/html/static/phenomenonmap.js`, `mapzoom.js` | Zoom and pan for the phenomenon page's AU-scale SVG diagram (`mapzoom.js` is the shared viewBox zoom/pan it uses). Served directly, same as `style.css`. |
| `../src/html/static/phenomenonrender.js` | Draws `planetgen/web/maps/phenomenonrender.py`'s "View" panel with three.js. Served directly, same as `style.css`. |
| `../src/html/static/generatefolds.js` | The Generate page's folding sections (remembered per browser) and "Around a sector"'s finder: name search through `/galaxy/locate` and the paged list from `/admin/generate/sectors`. Served directly, same as `style.css`. |
| `../src/html/static/generatejobs.js` | Keeps the Generate page's "Current job" panel live by polling `/admin/generate/status`. Served directly, same as `style.css`. |
| `../src/planetgen/web/maps/systemscene.py` | The 3D system view's data (MAP.69): `build_scene` turns `system_detail` into one JSON scene (every star, planet, moon, belt and comet with its ref, radius, colour, orbit elements, position at the epoch). Served by `GET /api/systems/<id>/scene` and relayed to the page at `/system/<id>/scene`. |
| `../src/html/static/systempage3d.js` | The system page's Diagram / 3D switch (MAP.74): loads the 3D view and the scene the first time 3D is chosen, remembers the last view, keeps `?view=3d&object=<ref>` in the address (MAP.67), wires the time control, scale select and the screen-reader body list, and says so when WebGL is missing. |
| `../src/html/static/systemview3d.js` | The 3D system view (MAP.72, MAP.73): the three.js scene (glowing stars, lit planets and moons, orbit lines and Kepler ellipses, belt particles, heliopause), the camera (drag to turn, pan, wheel and pinch zoom, W/A/S/D and arrow keys, fly to, follow, reset), picking, the info panel and decluttered labels. |
| `../src/html/static/systemscale.js` | True and compressed scale (MAP.71): where and how large each body is drawn, and the note saying which is shown. |
| `../src/html/static/orbitpositions.js`, `../src/html/static/orbitclock.js` | Positions of every body at any time and the time control (MAP.70); `orbitpositions.js` is the twin of `physics/body_positions.py`. |
| `../src/html/static/mapcore.js` | The helpers the Galaxy Map, the Sector Map and the System Map share (MAP.63): reading the panel's scene JSON, theme colors, info-panel fields, a sector's address, the highlight ring's texture, the scale bar's nice numbers and screen span, fitting the renderer to its canvas, and picking a point of light on screen. Served directly, same as `style.css`. |
| `../src/html/static/mappick.js` | The picking, hover and info-panel layer the Galaxy Map and the Sector Map share (MAP.65): a picker that holds each map's layers in priority order (points of light picked on screen, meshes raycast, or a map's own test) and gives back what is under the pointer, a solid layer hiding what is behind it; the hover tooltip (`.map-tooltip`, kept inside the map); the selection and hover rings, a size on screen for a point of light or round a body in the scene; and the info panel both maps fill the same way: a heading, rows, the NAV buttons ("Start Here" and "End Here", or the one for the end being picked, `navpick.js`), the ☆ Bookmark button (`bookmarks.js`), page links, buttons and Generate buttons. The Sector Map gains a hover tooltip and ring and a ☆ for a system, phenomenon or generated neighbor; the Galaxy Map's black holes, neutron stars, quasars and clouds gain a tooltip, a hover ring, the NAV links and a ☆. Served directly, same as `style.css`. |
| `../src/html/static/navpick.js` | The NAV course pick the maps share (NAV.29, NAV.33): which end is being chosen and the end already chosen, the `?pick=` query every map URL carries, the info panel's "Start Here" and "End Here" buttons (a link to `/nav` with both ends once one is chosen, else a button that keeps the user on the map and asks for the other end), the info panel's "What's within 10 pc" link to the What's nearby page with the object, a sector or a cell's middle filled in (NAV.45; none while a course is being picked), and the banner's words. Served directly. |
| `../src/html/static/mapcontrol.js` | The camera and input controller the Galaxy Map's drill-down and the Sector Map share (MAP.64): turning and moving an orbit camera, zoom policies (free, a short range, or locked), the wheel, a two-finger pinch, the arrow keys, and telling a drag from a click. Served directly, same as `style.css`. |
| `../src/html/static/starlight.js` | The brightness boost for faint stars (MAP.87), the Galaxy Map's twin of `planetgen/web/maps/starmap.py`'s `star_light_boost`, which sizes the Sector Map's points: up to four times the halo light at the dim end, none from 1000 L☉ up. Served directly, same as `style.css`. |
| `../src/html/static/unitladder.js` | The maps' and panels' unit formatter, the browser mirror of `planetgen.util.format.format_quantity`: a distance in km, AU, mpc, ly, pc and up, a speed in km/h to c, a time in µs to Gy, each in the largest unit it is at least 1 of. It reads `unitladders.js`, generated from `planetgen/util/ladders.py` by `scripts/build_unit_ladders.py` (a test fails when the file is out of date). Served directly, same as `style.css`. |
| `../src/html/static/numberformat.js` | The browser mirror of `planetgen.util.format.format_number`: a number shows as scientific notation ("1.23 × 10⁶") from 7 whole digits when it has no decimals, and from 5 whole digits when it has decimals (UX.36). Templates use the same rule through the `num()` global. |
| `../src/html/static/units.js`, `unitpreference.js` | The reader's unit system (UX.78): the server writes each measurement as `<span class="qty" data-q=... data-si=...>` with the Automatic text (`web/lib/fmt.py`'s `quantity_html`); `units.js` re-renders it from the SI value with the shared unit ladders and runs the settings menu's "Units" button, and `unitpreference.js` keeps the choice (Automatic, Metric only or Customary) in `localStorage`, guarded so the page works without it. |
| `../src/html/static/datatable.js`, `datatablestate.js` | The site's data table (UX.41): `datatable.js` takes over the table a list page rendered (`partials/datatable.html`, `web/lib/datatable.py`) with TanStack Table (sort and filter state) and TanStack Virtual (row scrolling), fetching further pages from `/table/<name>`; `datatablestate.js` is its logic without the page (the query a state turns into, which pages the visible rows need, the store of fetched pages) so plain node can test it. |
| `../src/html/static/breadcrumb.js`, `pagecrumbs.js` | The one breadcrumb every page and map uses (NAV.14): a single line at any width, with a "…" menu for the steps it hides. `pagecrumbs.js` folds a server-rendered page's breadcrumb list onto that line; without script the list wraps. |
| `../src/html/static/picker.js`, `objectref.js`, `systemnav.js` | The one selection every map and page shares (NAV.13: `select`, `up`, `into`), the browser mirror of `planetgen/galaxy/objectref.py`'s public reference form (`<kind>:<id>`, NAV.7; ids stay strings), and the NAV buttons and the "What's within 10 pc" link of a body's info panel on the System Map (NAV.50, NAV.45). |
| `../src/html/static/mapnavpad.js`, `gridopacity.js` | The Galaxy Map's on-map controls (MAP.200: a zoom slider, a direction star and a readout of the center), and the cap that keeps the drill-down's stacked translucent blocks from ever being more than 50% opaque along a line of sight (MAP.188). |
| `../src/html/static/nebulamesh.js`, `nebulaview.js` | A nebula drawn from its shape instead of a flat sprite (GEN.75, MAP.103): the triangle mesh `/galaxy/nebula/<id>/shape` serves, shared by the maps; and the nebula page's 3D view (MAP.105) with the brightest stars round it from `/galaxy/nebula/<id>/surroundings`. |
| `../src/html/static/tilepack.js` | The reader of the packed Galaxy Map tiles `/galaxy/tiles` sends and the tile store keeps (MAP.159); `web/lib/tilepack.py` writes them and describes the layout, and `tests/test_tilepack.py` with `tests/js/tilepack.test.mjs` hold the two together. |
| `../src/html/static/menuplace.js`, `navtimes.js`, `connectionnote.js` | Button menus that open where they can be seen (UX.85: above the button when there is more room, height capped); the NAV result's Warp/Fold travel-times switch (UX.70); and when to tell an admin a live job stream was lost (ADM.40; the stream ends on purpose every `STREAM_SECONDS`). |
| `../src/html/static/joblog.js`, `generatecustomize.js`, `generateranges.js`, `generatestarmix.js` | The Generate page's job output box (ADM.22: an Xterm.js terminal, falling back to a `<pre>`), its "Customize" window (ADM.28: fieldsets become tabs), the range sliders' read-outs (GEN.183) and the star mix's running total, which must reach 100% (ADM.45). |
| `../src/html/static/components.js`, `shoelaceicons.js`, `shoelace-theme.css` | The shared Shoelace component set (UX.40, vendored under `vendor/shoelace/`; `components.js` is loaded as a module on every page and registers the elements), the icon libraries pointed at files under `vendor/shoelace/assets/` (Shoelace's own data-URI icons would break the Content-Security-Policy), and the theme that maps Shoelace's tokens onto `style.css`'s. `icons.svg` is the site's icon sprite. |
| `../src/html/static/localtime.js` | Rewrites every server-rendered UTC `<time data-local-time>` in the viewer's own time zone. Without script the times stay readable, labelled UTC. Served directly, same as `style.css`. |

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
if a query were buggy -- see `planetgen.db.query`'s module docstring for the same
convention, and [`api.md`](api.md) for where that account is configured.
Database and system names pulled from generated data are HTML-escaped
before being placed in a page; a requested `?db=` schema name is
validated API-side against the actual, prefix-filtered schema listing
(exact match only, `planetgen.db.store.resolve_database`), which is what
prevents it from being used to select a schema this deployment never
meant to expose. `mdrender` escapes every block in full before emitting
any markup, then narrowly re-enables only the one legitimate raw-HTML
pattern generated content ever contains (`<sup>...</sup>`) -- so a
mischievous `--name`/`--star-type` value can't inject live HTML into a
rendered page.

## The page shell

Every page extends `web/templates/base.html`, and every HTML response
carries the security headers (`web.SECURITY_HEADERS`, added in
`web/app.py`: a `Content-Security-Policy` of `default-src 'self';
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
- The favicon, `static/favicon.svg` (an uploaded icon replaces it), and `static/og-default.png`, the built-in 1200 by 630 link-preview picture.
- `static/theme.js`, a small blocking script placed before the
  stylesheet: it applies the visitor's saved theme before the first
  paint (see below).
- `static/style.css`.

Every `static/` URL, here and in the pages' own `<script>` tags, comes
from `planetgen/web/lib/fmt.py`'s `static_url(name)`, which appends `?v=<package
version>` (read from `src/planetgen/_version.py`, which the release
Action stamps). A release therefore changes every static URL, and the
web server can let browsers cache them for a year (see
[`deployment/README.md`](deployment/README.md#what-every-setup-has-in-common);
each platform guide in [`deployment/`](deployment/README.md) sets the same
rule).
The ES modules that import siblings (`galaxymap3d.js`, `systemmap.js`,
`galaxymap3d.js`, `galaxystageview.js` -> `mapcore.js`, `mapcontrol.js`, `mappick.js`,
`starlight.js`, `bodyRendering.js`,
`galaxyprisms.js`, `galaxyblocks.js`, `vendor/three.module.min.js`) use `await import(...)`
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
(`../src/planetgen/web/`, registered by `web/app.py`'s `create_app`). The
pages used to be one CGI script each; the "Replaces" column names the
old script, whose URL now answers with a 301 to the new page (see "Old
URLs" below).

The site has no help or about page yet. DOC.6 in the TODO adds a static help section (a page template, an index,
a help link on each page and a test that every page has one); until then the pages explain themselves.

| URL | Replaces | Shows |
|---|---|---|
| `/` | `index.py`, `browse.py` | The front door (UX.55): the counts of what is here (sectors, systems, phenomena ...), each a link to its list, and links to the Galaxy Map, NAV and Search. |
| `/sectors` | `browse.py#sectors` | Every sector, nearest the galactic core first, as a data table (UX.41). |
| `/systems` | `browse.py#standalone-systems` | Every system with its sector and octant as a data table; the standalone ones (made outside any sector) are a filter on the list (UX.55). `?uncharted=1` (UX.87; a badge on the page counts them) swaps it for every star the brightness scatter placed that has no system yet, brightest first: luminosity, star, temperature, age, sector address as shell, band and wedge (linking to the Galaxy Map there), galaxy and in-sector coordinates; a Members column ("3 members, 2.4 solar masses" for a star with company, a dash for a single star) and each star's name links to its page, `/systems/uncharted/<id>` (UX.114: the star's data, its members and the orbits that bind them; a built star sends the visitor to its system); rows from `/table/uncharted_systems`. An admin gets a Generate button per row (`POST /systems/uncharted/<id>/generate`) that builds just that system and opens it; a note recommends generating the whole sector. |
| `/galaxy` | `galaxy.py` | The 3D Galaxy Map plus the Quadrant summary; `?quadrant=I\|II\|III\|IV` lists that Quadrant's placed sectors in a data table (UX.41), nearest the core first; sort, `zone` filter and `page` in the address, rows from `/table/galaxy-quadrant`. The info panel's "View sector" is a plain link. `?pick=from&to=...` or `?pick=to&from=...` (NAV's "Pick on Galaxy Map") shows a banner with Cancel back to `/nav`, keeps "Charted only" on, and keeps the pick on the sector it opens in place (`?sector=<designation>&open=1`) and on a sector's page (`/sector/<id>?pick=...`). |
| `/galaxy/stage?at=...` | `galaxy_views.py` | JSON for the map's drill-down: one stage's generated counts, `at` a block key `m.shell.band.wedge` or none for the galaxy (see `/api/galaxy/stage`). Through `planetgen/web/lib/tilecache.py`'s disk cache (`fetch_stage`); 400 on a malformed key, 502 on an API failure, both `{"error": ...}`. Shares the `galaxy_tiles` rate limit. `Cache-Control: no-store`. |
| `/galaxy/window?at=...&g=...` | `galaxy_views.py` | JSON for the map's data window (MAP.151): `g` cells across (3, 5, 7, 9 or 11; 9 by default) around the cell `at`, a block key `m.shell.band.wedge` or a sector `1.shell.band.wedge` (see `/api/galaxy/window`). Built from the stages of the parent blocks it touches (`planetgen/web/lib/tilecache.py`'s `fetch_window`); 400 on a malformed key or width, 502 on an API failure, both `{"error": ...}`. Shares the `galaxy_tiles` rate limit. `Cache-Control: no-store`. |
| `/galaxy/locate?q=...` | `galaxy_views.py` | JSON for the map's address bar: sectors and star systems named like `q`, each with its sector address (see `/api/galaxy/locate`). A blank `q` asks the API nothing; 502 on an API failure, as `{"error": ...}`. Shares the `search` rate limit. |
| `/galaxy/territories` | `galaxy_views.py` | JSON for the map's Territories overlay: `/api/territories`' owned systems and capitals, with each polity's name, color and system count folded in from `/api/polities`. 502 on an API failure, as `{"error": ...}`. Shares the `galaxy_tiles` rate limit. `Cache-Control: no-store`. The map shows its Territories button only once `/api/polities` counts at least one polity. |
| `/galaxy/uncharted/<shell>/<band>/<wedge>/scene` | (new) | JSON scene of a sector cell nothing was generated in (MAP.162): the waiting bright stars and unbuilt scattered phenomena the scatters left there, for the map to open in place with an "Uncharted" mark. A generated cell redirects to its sector's scene; a cell outside the grid is a 404 and an API failure a 502. Shares the `galaxy_tiles` rate limit. |
| `/galaxy/nebula/<id>/shape`, `/galaxy/nebula/<id>/surroundings` | (new) | JSON for the maps' nebula meshes (MAP.103, `?lod=low\|full`; `GET /api/nebulae/<id>/shape`) and for a nebula page's 3D view (MAP.105: the brightest stars round it). 404 for an unknown nebula, 502 on an API failure. Shares the `galaxy_tiles` rate limit. |
| `/galaxy/tiles?tiles=...` | `galaxy_tiles.py` | JSON for the map's script: the requested cube tiles (`tiles=level/ix/iy/iz,...`) and, given the browser cache's `stamp`, the changed tiles since. Through `planetgen/web/lib/tilecache.py`'s disk cache, which stores each tile as trimmed JSON bytes (MAP.157: only `clouds`, `stars`, `generated` and `points`; a star keeps `id`, `x`, `y`, `z`, `luminosity_sol`, `temperature_k`, `radius_sol`, the class letter in `star_type` and `system_id`). The body is those bytes joined, not parsed and serialised again (Apache compresses it, Brotli when `mod_brotli` is on). 400 on a malformed request, 502 on an API failure, both `{"error": ...}`. `Cache-Control: no-store`. |
| `/admin/generate/sectors` | (new) | Admins only: JSON, 50 filled sectors a page (`?page=N`), for "Around a sector"'s list (`static/generatefolds.js`). |
| `/admin/generate/status` | (new) | Admins only: JSON for the Generate page's live job panel (`static/generatejobs.js`). |
| `/admin/generate/jobs/<id>` | (new) | Admins only: one past job with its full output. |
| `/admin/generate/jobs/<id>/stream`, `/admin/generate/jobs/<id>/log` | (new) | Admins only: a job's output as Server-Sent Events (`log` events with new output and the byte offset reached as the event id, so a reconnect resumes; `state` events; `done` once all output is sent; the response ends after `STREAM_SECONDS` and the browser reconnects), and its full output as a plain-text download. |
| `/admin/generate/system/download` | (new) | Admins only, POST: returns the one-off system page's text as a `.md` or `.wiki` file. |
| `/system/<id>` | `system.py` | One star system: badges, "Navigate from/to here" (systems in a sector), nearest-neighbour location links, the System Map, the expandable body list, `?code=wikitext\|markdown` code views with a Copy button, the Stars/Planets/Belts/Comets tables, the Facilities panel (with each facility's type, affiliation and rating chips, worst first), the "Species that could live here" panel (POP.12: for each rocky world, the species that could live there unaided and those that could with equipment, each with its six domain chips; hidden when there are none), and for an admin the "Upload to Wiki" form and the facility forms: place one (with a type and its fields, an affiliation and ratings) and edit a facility's affiliation and ratings or have the unset ratings suggested (see below). |
| `/phenomena` | `phenomena.py` | Every exotic phenomenon, paged with `?page=N`. An unbuilt scattered row with company shows its members and total mass in a Members column, and its page (`/uncharted-phenomenon/<id>`) lists the members and the orbits that bind them (UX.114). |
| `/object/<ref>` | `system_pages.py` | NAV.8: redirects any object reference (`planet:12`, `moon:7`, `belt:2`, `comet:9`, `star:3`, `system:9`, `sector:4`, `nebula:3`; a bare number is a system) to the page that shows it. Only systems, sectors and phenomena (and a star, as its system) have pages; a planet, moon, belt or comet opens `/system/<id>#planet-12` (and `#moon-`, `#belt-`, `#comet-`, `#star-`), where `static/bodyanchor.js` opens the row, highlights it, scrolls to it and selects the body on the System Map (a moon enters its planet's scene). A single-star system has no star row, so a star link there just opens the page. Search results link their stars, planets, moons and belts this way. |
| `/phenomenon/<type>/<id>` | `phenomenon.py` | One phenomenon's data table and its view: a three.js render for a neutron star, black hole, quasar, rogue planet or comet (`planetgen/web/maps/phenomenonrender.py`, `static/phenomenonrender.js`, an SVG still without JavaScript), the AU-scale diagram for a nebula or remnant, and none for an asteroid field, with "Navigate from/to here". `<type>` is one of `nebula`, `asteroid_field`, `black_hole`, `neutron_star`, `supernova_remnant`, `rogue_planet`, `interstellar_comet`, `quasar`; anything else is a 404. |
| `/classes` | (new) | The class reference: every class type (star spectral and luminosity classes, planets, nebulae, supernova remnants, asteroid fields, black holes, rogue planets, comets) with how many classes it has. |
| `/classes/<type>` | (new) | One type's classes, each linking to its page, plus notes (an asteroid field's size digit, for one). `<type>` is `star-spectral`, `star-luminosity`, `planet`, `nebula`, `supernova-remnant`, `asteroid-field`, `black-hole`, `rogue-planet` or `comet`; anything else is a 404. |
| `/classes/<type>/<code>` | (new) | One class's facts, e.g. `/classes/planet/M`, `/classes/star-luminosity/IA+`, `/classes/comet/halley_type`. An unknown code is a 404. The system page links a star's type and a planet's or comet's class here, and the phenomenon page its Class (an asteroid field's `C3` by its letter) and a rogue planet's Mass Class. |
| `/classes/habitability` | (new) | What the equipment levels and PHI-6 colours on the planet pages mean (UX.90): the five equipment levels, the six domains with their Blue, Green, Yellow and Red thresholds, and the chips a planet row shows. |
| `/species` | (new) | Every species, paged with `?species_page=N`; `?spacefaring=1` or `0` filters. Like every population page it is a 404, and the header's Species section is hidden, until a population pass has made species (`GET /api/population`). |
| `/species/<id>` | (new) | One species: homeworld, body plan, era, civilization age and polity, and its environment limits (POP.12): each limit of the habitability table with the species' value, the human value and whether they differ. |
| `/polities` | (new) | Every polity (`?polities_page=N`): species, government, capital, systems held, reach. A 404 until there are polities. |
| `/polities/<id>` | (new) | One polity and the systems it holds, nearest its capital first (`?systems_page=N`). |
| `/sector/<id>` | `sector.py` | One sector: badges, the 3D Sector Map, its Profile panel (UX.113: what the galaxy's shape says the cell holds, from `GET /api/galaxy/profile`: region, expected systems, the stellar populations with age, [Fe/H] and velocity dispersion, gas, superbubble and dark matter), and its Contents table (systems, nearby phenomena and the facilities outside its systems, nearest the center first, as a data table (UX.41): `contents_sort`, `contents_order`, `contents_type`, `contents_octant`, `contents_page`; rows from `/table/sector-contents`); for an admin, an Admin menu (wiki upload, generate neighborhood) and an Edit menu, each pick opening its form in a dialog. |
| `/sector/at/<shell>/<band>/<wedge>` | (new) | Any cell of the grid by address (UX.113): a charted cell redirects to its sector page; an uncharted one shows its place, the bright stars and phenomena the scatters left waiting in it, its Profile panel and, for a heavy cell (one that expects more systems than the cell star ceiling, or whose fill is expected to take two minutes or more on one worker), a warning with the size and time before anything is run. A cell that is not in the grid is a 404. |
| `/sector/<id>/galaxy`, `/system/<id>/galaxy` | (new) | Redirect to `/galaxy?sector=<designation>` (the map opened on that sector, selected), or to the plain map for a sector with no galaxy address or a standalone system. Search results link here; the sector page links straight to the map. |
| `/sector/<id>/scene`, `/system/<id>/scene` | (new) | The scene JSON the Galaxy Map opens a sector or a system from (MAP.66, MAP.125): the Sector Map's block (`?pick=...` keeps NAV's pick mode in the links; `?binaries=real` draws a close pair at its real gap, MAP.186), and the 3D system view's scene (`GET /api/systems/<id>/scene`, MAP.69). Both share the `galaxy_tiles` rate limit. |
| `/nav` | `nav.py` | The NAV route planner; see "The NAV page's URLs" below. |
| `/nav/chart` | (new) | Admins only (NAV.48): the uncharted sectors that block a course. Shows what charting would take (the cells of the route's hops through unknown space, whether a route through charted space exists, a one-cell border to tick) and a form that sends the job to the Generate page, which shows the size and time and asks to confirm before it starts. A visitor who is not logged in is sent to log in first. |
| `/nearby` | (new) | What's nearby (NAV.44): pick a place (a sector, then a system, or a point in the galaxy frame), type a distance in parsecs and list what is there, from `GET /api/near`; `kinds` filters by object kind and `page` pages. Every step is a bookmarkable GET. |
| `/table/<name>` | (new) | One page of a data table as JSON, for `static/datatable.js` (UX.41): `?sort=&order=&<filters>&offset=&limit=` (and `after=` for a keyset table: Systems, Uncharted stars and Phenomena, PERF.75), `facets=1` for the filter menus' counts. Answers `{"rows", "total", "capped", "next", "facets"}`. An unknown table is a 404, a read that hit the statement time limit a 504 with `Retry-After`, any other API failure a 502. |
| `/maintenance` | (new) | What a visitor sees while maintenance mode is on (UX.112); with it off the page sends the visitor home. An admin sees the same page with a link to the admin page's switch. |
| `/galaxy-status` | `galaxy_status_page.py` | Public JSON for the banner every page shows while a background run changes the galaxy (UX.3): `{"active": false}`, or `active`, `hours` (whole hours, rounded up by `queue/hourdisplay.py`; `null` while there is no estimate), `percent`, `stage` and `text`. Read from the run record in Redis through a 5 s per-process cache (`web/lib/galaxystatus.py`); `static/galaxybanner.js` polls it every 30 s. |
| `/search` | `search.py` | Faceted search (see below). |
| `/login` | `login.py` | The admin login form (`?next=<local path>` to return to afterwards). An admin with two-factor sign-in on is then asked for the authenticator or recovery code, which `POST`s to `/login/code`. |
| `/logout` | `logout.py` | `GET` asks to confirm and changes nothing; the button `POST`s to end the session. |
| `/account` | `changecreds.py` | Change the admin username and password, and turn two-factor sign-in on or off (QR code, then recovery codes shown once; forms `POST` to `/account/two-factor`). |
| `/admin` | `admin.py` | API keys (list, create, revoke; `?keys_page=N`) and a sector's manual wiki link. |
| `/admin/stats` | `adminstats.py` | Server health (including Redis: answering, append-only file, memory policy, memory used, last save and live workers) and database stats (including about how many bright stars the plan pre-placed, from the `bright_stars` table's row estimate), every name made unique (`?names_page=N`), the galaxy's naming key with Change and Draw buttons (POST `/admin/stats/naming-key`, GEN.70), current login lockouts with Lift buttons (POST `/admin/stats/lockouts`), the newest failed sign-ins, and this server's generation speed per density bucket with the galaxy's size per star system (`GET /api/admin/generation-stats`). |
| `/admin/stats/generation-reset` | (new) | Admins only: POST, the Stats page's "Reset stats" button (deletes every recorded generation rate, PERF.32). |
| `/admin/queue` | (new) | Admins only: the work queue (ADM.10). Workers active and the server's load as "x / x / x" (1, 5 and 15 minutes), who holds the workers' lease and how long since it was refreshed, Pause the queue / Resume the queue, the bulk worker count (Set, or Ludicrous Speed: 95% of the cores; ADM.15), Clear the stale lease, then every job tree newest first (paged) with its state, start, end, duration, progress, task counts and ETA. `/admin/queue/<id>` shows one tree as nested expandable nodes down to single tasks, each with Pause, Resume, Cancel, Retry (a failed or cancelled job, or one failed sector) and Delete (a finished job). Every control first shows `/admin/queue/confirm/<action>?node=<id>`, then posts to `/admin/queue/action`. Retry starts a Generate page job. |
| `/admin/settings` | (new) | Admins only: every option of `config.json` (ADM.84), grouped by category, each with its help text, default, limits, unit, whether it needs a restart, where the value in use comes from (default, `config.json` or an environment variable) and a box fitted to its type (check box, number, text, one-per-line list, choice, write-only secret box). An option the admin may not change (command line only, or the Owner's) is listed read-only with the reason. "Use the default" removes a key from the file; "Clear it" empties a secret. Review the changes (POST) shows the old and new value of each changed option, which need a restart and which an environment variable overrides, or sends the form back with the message at each bad option and nothing saved; Save these changes (POST) then writes the file (a backup archive first, an audit row per option) with the hash of the file the form was drawn from, so a file changed by someone else in between is refused. The page works without script. While a saved option still needs a restart of the web server, every admin page shows a banner naming it. |
| `/admin/settings/backups`, `/admin/settings/backups/<archive>`, `/admin/settings/backups/<archive>/restore` | (new) | Admins only: the backup archives screen (ADM.85). The list shows each archive's time, who made it, why, which parts it holds and its size, with Back up now (POST) and, for the Owner only (from 8.2), Download. An archive's page says whether it passes its hash check, compares its settings part with today's `config.json` (old and new value of each option a restore would change, secrets as `(changed)`; a backup that no longer fits the model is explained and refused), and restores the ticked parts (POST): a safety archive of what is replaced first, an audit row per changed option and per restore, only the options an admin may change from the web. A settings file changed by someone else since the page was drawn is refused. |
| `/admin/settings/images` | (new) | Admins only (ADM.44): the site icon (SVG), the PNG or ICO icon, the home-screen icon and the link-preview picture in use, each with an upload form (checked when uploaded; the reason is shown when a picture is refused) and "Use the built-in one" (POST). The pictures are stored by content name and served at `/site-assets/<name>`. |
| `/robots.txt`, `/sitemap.xml`, `/sitemap-pages.xml`, `/sitemap-systems-<n>.xml`, `/sitemap-sectors-<n>.xml`, `/.well-known/security.txt`, `/humans.txt`, `/site.webmanifest`, `/site-assets/<name>`, `/og/<system|sector>/<id>.png` | (new) | Public (ADM.44): the files crawlers, browsers and link previews ask for, generated from the `seo`, `icons`, `og` and `site` settings at request time (docs/config.md, "Web, link previews and search engines"). The sitemaps and the manifest are 404 when switched off or while `base_url` is the localhost default. |
| `/admin/settings/export`, `/admin/settings/import`, `/admin/settings/import/save` | (new) | Admins only: export downloads the options set in `config.json` that an admin may change as `planetgen-settings.json`, without secrets or command-line-only options; import takes such a file (up to 1 MiB), checks it like an edit, shows the review diff, and `import/save` saves it as `config.import` with a backup first. A bad value, an unknown key, an option an admin may not set or a file from a newer release changes nothing and is named. |
| `/admin/facility-types` | (new) | Admins only (POP.10): every facility type with its fields, and forms to add, change, take out of use and delete types and add, change and delete their fields, or to add the example type; each form posts back to the page, the API holds the rules and a refusal shows as the page's error. |
| `/admin/generate` | (new) | Admins only: generate, plan or reset the galaxy from the browser (see below). |
| `/admin/generate/system` | (new) | Admins only: one star system with every `planetgen system` option, shown as Markdown or wikitext and never saved (see below). |

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

**Delete and Regenerate buttons (ADM.8).** For an admin whose
credentials are current, the sector, system and phenomenon pages each
have one Admin menu (`web/edit_actions.py`,
`templates/partials/admin_menu.html`; ADM.34) listing only that page's
actions, hidden from visitors; no admin form sits inline on the page. The
system page's menu holds the system's own actions, a submenu per planet,
moon and asteroid belt (Regenerate, Delete, Change class), Place a
facility and a Remove a facility submenu; the sector and phenomenon pages
have one pair for the page's own object (not for a system's own black hole
or neutron star). Each item opens a dialog (`static/dialogmenus.js`) with
a short note of what it does, an "also delete facilities" checkbox where
facilities could be lost, and a "Yes, ..." submit button. The form POSTs `edit_action` (`regenerate`
or `delete`) and `edit_target` (`<kind>:<id>`, only targets the page
shows are accepted) to the page itself, which calls the API (see
`docs/api.md`, "Deleting and regenerating"), flashes the outcome, the
bodies the re-validation moved and any warnings, and redirects (303):
back to the page, or after a delete to the system's sector (or the
Systems list), the Sectors list or the Phenomena list, and after a
sector regenerate to the new sector. Every page shows those flashed
lines under its heading (`base.html`), read only when the request
carries a Flask session cookie.

**Change class and Change star (ADM.6, ADM.7).** In the same menu on
the system page, every planet and moon also has a Change class item:
its dialog lists the recommended classes first (from `GET
/api/systems/<id>/class-options`, ones that fit without moving
anything) and every other class under "Needs confirming"; the form posts
`edit_action=class` with `planet_class` (`M`) and, when the "Set it
anyway" box is ticked, `confirm_class=1`. An unticked box with a class
that breaks something flashes what it would break and changes nothing
(ADM.51). A
single-star system's own row has Change star, a spectral-type field
(`edit_action=star`, `star_type`) with the same "also delete
facilities" checkbox, since bodies the new star can't hold are removed.
Both flash the outcome as above, including what moved and what was
removed (see `docs/api.md`, "Changing a class or a star").

**System page facilities.** The system page lists the system's
facilities (starbases, colonies, outposts; `GET /api/systems/<id>/
facilities`) in a Facilities panel (name, kind, host, placement, and an
orbital one's distance, period and speed as stored), and each one again
in its host's row of the body list. The System Map draws each as a small
diamond at its host (`planetgen/web/maps/systemmap.py`'s `_facilities_svg`): a star's on
a dashed orbit on the map's own scale, a planet's or moon's just outside
its marker (and a planet's again around the center of its moon view), a
belt's on the ring. A colony makes its world "Inhabited"
(`queryDb._with_life_fields`). For an admin the panel has a form
(`web/system_facilities.py`): a name, the placement (in orbit, on the
surface, in the belt), the host (only the stars, planets, moons or belts
in the system that take that placement), the kind, and an optional
description. For "in orbit" a logarithmic slider sets the distance, from
just above the host's surface to the edge of its sphere of influence (a
planet's or moon's Hill sphere, a star's heliosphere); a surface or belt
facility has no distance control, and a belt one gets a random spot in
the belt when saved. It POSTs to `/system/<id>` with
`{{ csrf_field() }}` and a `facility_action`: `preview` checks the
placement rules (`planetgen.population.facilities.check_facility`) and shows the
orbit `GET /api/facilities/orbit` works out from the host's mass, saving
nothing; `save` calls `POST /api/facilities` and answers `303` back to
`/system/<id>?facility=added#facilities`; `remove` (each row's Remove
button, behind a `<details>` confirm step) calls `DELETE
/api/facilities/<id>` for one of this system's own facilities. Errors,
including the API's reason for a refusal, show next to the form.
`static/facilityform.js` hides the hosts that don't fit the placement,
shows the slider only for "in orbit" and reads out its distance, period
and speed as it moves; without it the form still works and the server
refuses a host that doesn't fit.

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

The Galaxy Map's disk tile cache (`planetgen/web/lib/tilecache.py`) is written by the
WSGI daemon. A tile cache location or size set only with `SetEnv
PLANETGEN_TILE_CACHE_DIR`/`_MAX_MB` in the vhost does not reach it: put it
in `config.json`'s `tile_cache` (or the app server's own environment) instead.

### The sector page

`/sector/<id>` (`web/sector_page.py`, `templates/sector.html`) shows the
sector's badges (edge, counts, a link to its Galaxy Map quadrant),
the interactive Sector Map and one Contents table of its systems, the
phenomena near it and its facilities outside any system (stand-alone ones
parked in open space and those on its asteroid fields, `GET
/api/sectors/<id>/facilities`; a facility has no page, so its name is not
a link), nearest the center first, 50 per page (`?contents_page=N`).

The Sector Map is the Galaxy Map's own engine locked to the sector
(MAP.68): `render_galaxy_map3d_panel(..., pinned=...)` draws the panel with
no breadcrumb, steps or address bar, and
`static/galaxystageview.js` opens the sector in place (`static/
galaxysector.js`, the scene from `GET /sector/<id>/scene`) and keeps to it:
Up, Home, Escape and the drill-down's picks do nothing, there is no history
entry or query of its own, and a click on a neighboring sector opens that
sector's page. The panel keeps zoom buttons, Re-center and "Mark rogue
planets" (in its Menu), a screen-reader list of the sector's entries, and
the Contents table's "Show on map" buttons. Every entry carries a plain
`href`, so the info panel's "View system/phenomenon/sector" buttons are
ordinary links. A sector with no place in the galaxy has no map.

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
| `/nav?from=system:12&to=system:40&stay=30` | NAV.11: the route's panel also shows the travel time of each hop and of the whole route at every warp and fold factor, with a stay of 30 minutes at each stop between the ends (default 0; the "Recalculate" field on the page sets it). |
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
| `type`, `spectral`, `luminosity`, `class`, `body`, `life`, `equipment`, `moon_class`, `moon_body`, `moon_life`, `moon_equipment`, `density`, `phenomenon`, `phenomenon_class` | Tag facets, repeated once per selected value (`spectral=G&spectral=K`). `equipment` and `moon_equipment` are the equipment level a human needs (0 to 4, linked to the habitability levels page); `phenomenon` and `phenomenon_class` pick exotic phenomena by type and class. |
| `sectors_page`, `systems_page`, `stars_page`, `planets_page`, `moons_page`, `belts_page`, `phenomena_page` | Each result panel's page; each pager keeps the search and the other panels' pages. |

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
  (Galaxy, Sectors, Systems, Phenomena, Nav, Nearby, Classes, and Species once
  a population pass has made species) with `aria-current="page"`
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
| New galaxy | `python3 -m planetgen.cli.reset --yes`, then `planetgen plan --no-bright-stars`, then the bright-star scatter (`planetgen plan --bright-stars-only`), then `planetgen galaxy` around a random start. |
| Generate sectors | `planetgen galaxy` in any of its modes: around a random start; around a sector (`--center-sector --radius-pc` for a filled sector found by name through `/galaxy/locate`, picked from a paged list of filled sectors (`/admin/generate/sectors`) or typed by ID; or `--shell --band --wedge --radius-pc` for a sector address or a galaxy-frame x, y, z in pc, turned into the address of the cell holding it with the plan's sector edge); one address (`--shell --band --wedge`, with an optional neighborhood radius); a span of shells, bands and wedges (`--shells`, `--bands`, `--wedges`, each a number or `FIRST:LAST`, blank for all; ADM.29), with `--limit`, or `--yes` for a very large one; a whole shell or one band of it (`--shell`, `--shell --band`, with `--limit` or `--yes`); or the sectors a course crosses (`--course`, NAV.48, from the NAV chart page). The single-sector neighborhood radius can be given in light-years (13 ly up to the parsec limit). Sent with `Accept: application/json`, the form answers `202 {"job", "url", "status_url"}` (or `{"error"}`) so the Galaxy Map can start a job without leaving the map. The Sector Map's Generate buttons on an unfilled neighbor post straight to this form. Before the job starts the page shows its size, time and the free space on the drive holding the database's data directory (named, or "unknown" if it can't be measured) (`planetgen galaxy --estimate-only`) with a Generate button to confirm, or the refusal when the disk can't hold it; a JSON caller gets `409 {"error", "estimate", "confirm_field"}` and re-sends with `estimate_ok=1`. |
| Plan the galaxy | `planetgen plan --no-bright-stars` with the galaxy shape fields, then the bright-star scatter as its own step (`planetgen plan --bright-stars-only`), so the job shows the scatter's progress bar and the count it placed. On New galaxy and Plan, "Skip the bright-star scatter" leaves that step out. |
| Redo scatters | `planetgen plan --redo-scatters` (GEN.196): one box where the admin ticks the scatters to redo (the massive stars, the bright stars by luminosity floor, the phenomena) and sets each one's new limit (the stellar mass limit, the luminosity floor, the neutron star and black hole limit with its "Use the star mass setting" box). They run as one job with every stage listed and the unticked ones shown as skipped ("not chosen to be redone"). A scatter clears and rewrites only its own rows; the central black hole or quasar is always rebuilt with the phenomena. A new mass limit redoes the luminosity pass too. |
| Add a dimmer layer of bright stars | `planetgen plan --bright-stars-down-to N`: keeps the bright stars already placed and adds only those from N up to the current star-fill level (shown on the panel), leaving filled sectors out. Disabled until a scatter has run. |
| Reset | `python3 -m planetgen.cli.reset --yes`. |

Every section of the page folds (a `<details>` whose summary is the
section's heading, ADM.4). Current job starts open; the others open or
close as the browser last left them (`static/generatefolds.js`, in
`localStorage`), except the section of a form shown again with its error
or estimate, which the server keeps open.

The number fields have upper bounds, the same ones `planetgen` checks
(`src/planetgen/generation/limits.py`): a radius of at most 200 pc,
shells up to the cap the density skeleton walks out to (`MAX_GENERATE_SHELL`), and at most 500 orbital slots on the one-off system
page.

New galaxy and Reset delete every generated row, so both need the
database name typed back. Every job writes the database this site shows,
passed to the child as `PLANETGEN_MYSQL_*` environment variables (so the
MySQL account needs the generator's grants, including `DROP` for
`TRUNCATE`). One job runs at a time; the page shows its step, a progress
bar (from `planetgen`'s `PLANETGEN_PROGRESS_FILE`, see
`planetgen/queue/progress_file.py`; the bright-star scatter's bar shows a
share done, with a second line for slow layers' stars, PERF.4 and PERF.9),
elapsed time and live output
(`static/generatejobs.js` polls `/admin/generate/status`), with a Cancel
button. The last jobs are listed with their full output at
`/admin/generate/jobs/<id>`.

**The one-off system page** (`/admin/generate/system`,
`web/system_page.py` + `templates/generate_system.html`, linked from the
Generate page) offers every `planetgen system` option: the ten
force/forbid choices, name, star type, age, orbital slots, the flavor
overrides, a pasted `--system-file` JSON, Markdown or wikitext, and the
`--debug` narration. It runs `planetgen system --output FILE` in a
temporary directory and waits for it (a system takes about a second), so
nothing touches the database. The result shows in a code box with Copy
and Download buttons (Download posts the text back to
`/admin/generate/system/download`, which returns it as a `.md`/`.wiki`
file), plus a rendered preview for Markdown.

A job is an RQ job on Redis (`redis.url`): the page queues
`planetgen.web.job_runner.run(<job dir>)` on a queue of the job's own and
starts one burst worker for it (`python3 -m planetgen.cli.worker`) in its
own session, so it outlives the request and a graceful reload (Apache's, or
gunicorn's). A full stop or restart of the service under systemd (which
stops everything in the service's cgroup), does stop it; the page then shows it as
interrupted. Cancel writes a `cancel` file into the job's directory; the
runner sees it within a quarter second and stops the running step's
whole process tree (`os.killpg`).
Without a Redis server no job starts and the page says so.
The worker is named after the job, so liveness comes from `/proc` on Linux and `os.kill(pid, 0)` on other POSIX
systems.
So closing the browser never stops a job (ADM.11): the Sector Map's
"Generate the neighborhood" button also starts a Generate page job
rather than running inside the request. To restart the web server
without stopping a running job, use `systemctl reload apache2` (what
`update.sh` suggests); after a full restart, retry the interrupted job
from Admin, Queue. Only one job runs at a time: the `active` lock in
the jobs directory is written with the job id already in it (a hard
link of a finished temporary file, or an exclusive create where there
are no hard links), and a stale lock is cleared under an
`active.clearing` mark, so two admins starting a job at the same moment
start exactly one (TEST.40). A lock naming no readable job is treated
as a job still starting for `STARTING_GRACE_SECONDS`, then as stale.
Pruning old jobs never removes the one the lock names.
Jobs live under `jobs.dir` (`docs/config.md`). A visitor who isn't a
logged-in admin is sent to the login page, and POSTs and status requests
without an admin session get a 403.

**Error pages.** Any failure of a page is an HTML page (`web/errors.py`) with the right status: 404 for a bad id or an unknown URL, 405, 413 and 429 for what they say, 502 "Data unavailable" when the API failed, and 500 with a fixed sentence (the traceback goes to the log, never the page). Two more reload themselves after `Retry-After` seconds, so a page opened while a long job holds the database fills in by itself (PERF.64, ADM.53): **503 "Busy"** when a request waited too long for a row a long job holds locked, and **504 "Took too long"** when a database read hit the statement time limit (`mysql.statement_timeout_seconds`). Under `/api` the same two answers are JSON, `{"error": "BUSY: ..."}` and `{"error": "QUERY_TIMEOUT: ..."}` (see [`api.md`](api.md#errors)), and `/table/<name>` answers a timeout as JSON too.

### How a Flask page is built

```
src/planetgen/web/
  __init__.py         blueprint `web`, template globals, init_app()
  views.py            /, /sectors, /systems, /search
  system_pages.py     /system/<id>, /phenomena, /phenomenon/<type>/<id>
  system_facilities.py  the system page's facility form
  edit_actions.py     the admin Delete, Regenerate, Change class and Change star forms (ADM.6-8)
  sector_page.py      /sector/<id>
  nav_page.py         /nav
  galaxy_views.py     /galaxy and its JSON routes
  class_pages.py      /classes and below
  admin_pages.py      /login, /logout, /account, /admin, /admin/stats
  generate_page.py    /admin/generate (jobs.py runs the jobs)
  system_page.py      /admin/generate/system
  searchpage.py       the search page's data shaping
  helpers.py          db_name, page_url, crumb, render_page, trusted_html,
                      pager, current_admin, generate_target
  transport.py        in-process transport for planetgen/web/lib/apiclient.py
  csrf.py             CSRF tokens for POST forms
  errors.py           HTML 404/502/500 pages, and the 503/504 pages
  old_urls.py         301s from the old /<name>.py URLs
  templates/          base.html + one template per page (+ partials/)
```

A page is a route plus a template (each page module is imported at the
end of `web/__init__.py`, which registers its routes on `bp`):

```python
# web/system_pages.py (simplified)
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
  {{ pager }}              {# Markup from planetgen/web/lib/pagination.py #}
</section>
{% endblock %}
```

The helpers (all in `web/helpers.py`):

- `render_page(template, title=, section=None, breadcrumbs=(),
  description=None, status=200, **context)`: renders a template that
  extends `base.html`. `section` is one of `SECTIONS` (`galaxy`,
  `sectors`, `systems`, `phenomena`, `nav`, `classes`) and gets
  `aria-current`.
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

**Data in-process.** Views call `planetgen/web/lib/apiclient.py` functions
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
  directly with the error and the username typed (with status 401 on
  `/login`, so the web server's access log shows the failed login; every
  failure is also a line in the [activity log](config.md#the-activity-log)).
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
"The page shell" above), set in `web/app.py`; JSON keeps `default-src
'none'`.

**Browser checks.** `src/tests/test_web_a11y.py` loads every GET route
of the `web` blueprint (read from the app's `url_map`, so a new page is
covered automatically) in headless Chromium at 390px and 1280px,
in the light and dark color schemes, against a small database generated
by `planetgen`. Each page must have no serious or critical axe-core
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
| `search` | `/search` and `/galaxy/locate` | 30 per minute |
| `galaxy` | `/galaxy` (the Galaxy Map page) | 60 per minute |
| `galaxy_tiles` | The map's data and scene routes, fetched by its script: `/galaxy/tiles`, `/galaxy/stage`, `/galaxy/window`, `/galaxy/territories`, `/galaxy/uncharted/.../scene`, `/galaxy/nebula/...` and `/sector/<id>/scene` | 600 per minute |
| `health` | `/api/health` | 60 per minute |
| `other` | every other page, all counted together | 300 per minute |

An empty value turns that limit off. These replace the API's default
limits (`ratelimit.default`) for the pages, and the API calls a page makes
in-process are not counted against either. Over a limit, a page answers
`429` with the site's HTML error page ("Too many requests"), while
the `/galaxy/...` JSON routes (read by the map's script) and everything
under `/api/` answer JSON; both carry `Retry-After`. The counts are on
Redis (`ratelimit.storage_uri` empty, the server in `redis.url`), so
every WSGI process shares them (see [`api.md`](api.md#rate-limiting)).

## Locating the database (and the API)

The pages reach the API in-process, so nothing needs pointing at it for
a normal deployment. `PLANETGEN_API_BASE_URL` (or `config.json`'s
`api_base_url`; default `http://127.0.0.1/api`) only matters to
`planetgen/web/lib/apiclient.py` used outside the app (a script). The *database*
server/account is the app's own configuration
(`PLANETGEN_MYSQL_HOST`/`_PORT`/`_USER`/`_PASSWORD`/`_DATABASE`, or
`config.json`'s `mysql` section, via `planetgen.db.store.MySQLConfig`) --
see [`api.md`](api.md#running-locally) for those.

Separately from all of the above, a `config.json`
file at the repo root (next to `src/`, not a file inside `../src/html/`)
holds every deployment-level setting in one place -- MySQL
connection details, the site's own `site_name`/`base_url`, and more --
edited once per deployment rather than passed through the vhost config --
see [`config.md`](config.md) for the full field list
and how it relates to the `PLANETGEN_*` environment variables.

## Deploying

[`../INSTALL.md`](../INSTALL.md) walks through a first install, and
[`deployment/README.md`](deployment/README.md) compares every supported
platform (Apache, nginx or Caddy on Linux; macOS; a VPS) and links a guide for each. The reference setup is
Apache2 with mod_wsgi on Debian or Ubuntu
([`deployment/apache.md`](deployment/apache.md)), where `sudo
./install.sh` does everything except the virtual host.

### Updating an existing deployment

Run `sudo ./update.sh` instead of pulling
manually. `git pull` on its own isn't enough: pulling a changed file
rewrites it with whatever mode is tracked in the repo, silently undoing
any executable bit `install.sh` previously fixed. `update.sh` fetches and
does a `git reset --hard` to origin's branch tip, so the checkout always
matches the branch. It warns about, then overwrites, uncommitted changes
to tracked files, and it never fails on diverged history. It never runs
`git clean`, so untracked files survive, `config.json` (gitignored)
among them. It then checks everything the site needs without
reinstalling anything that's already there: the executable bits, each
Python library (installing only one that's missing, too old or broken,
see [`deployment/apache.md`](deployment/apache.md#planetgens-own-python)), the
NLTK corpus, the schema migration (with the optional population pass
offered afterwards, a y/N question that defaults to no), Apache's
modules, ownership and permissions, the cache/jobs directories and the
debug log, and finally that the web app imports as Apache's user.
`install.sh` remains safe to run directly any time you want a full
reinstall without pulling first.

## Local testing without a web server

`python3 src/html/wsgi.py` serves the pages (and `/static/`) at
`http://127.0.0.1:5000/` alongside the API. Set `admin_cookie_insecure`
in `config.json` to log in over plain HTTP there.
{% endraw %}
