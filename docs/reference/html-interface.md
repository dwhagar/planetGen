{% raw %}
# planetGen Web Interface

A small web interface (in [`../src/html/`](../../src/html/)) for
browsing the MySQL databases described in
[`database-schema.md`](database-schema.md) -- pick a database (schema),
drill into its sectors and star systems, and view (or copy) the
description of each system, written from its database rows. Sectors and systems
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

Each rocky planet's and moon's description also carries the PHI-6 rating (UX.91, GEN.218, `physics/protection.py`): the score and the equipment a visitor needs, then one short sentence per category with its colour, score and the value that sets it.

Planet and moon rows carry habitability chips (GEN.89, design in
`docs/design/habitability/habitability-index.md`): equipment chips (the survival equipment line, GEN.655, UX.154: one chip per item a visitor needs, such as "Pressure suit" and "Dust-tight seals", solid when required and outlined when recommended, each in the colour of the PHI-6 domain it comes from; none for a world that needs nothing), a colour chip
that comes from the six PHI-6 domains (pressure, temperature, air, water,
radiation, suitability; Blue, Green, Yellow, Red) under the gates of GEN.207, and the Habitable,
Habitable moon and Inhabited chips. Planets around pulsars, neutron stars
and black holes follow a host rule. The Equipment search has a facet for
the labels. The explanation page is /help/habitability under the Help menu (UX.90, UX.180); every chip, and the symbols and life scores on the help page itself, link to their own section of it (UX.181); on a planet's or moon's page, in a list and in a tree a body's symbols and life scores link to the Habitability and Life groups of that body's own page (UX.197), whose Habitability group carries a link on to the help page, and a class letter or name links to its class page (UX.182). UX.100: a system page's planet and moon rows show
the overall PHI-6 symbol and the six domain symbols (UX.101) beside the equipment chip, instead of the notes (the
tooltip, the per-domain colour words and the three stored scores); each row names its class (and a "Habitable class"
chip for a habitable one, with no Terrestrial or Gas Giant label) and counts its moons. The Search page's planet and
moon panels have a Habitability column of the same symbols (a `symbols` cell, drawn by the page and by
`datatable.js`). The description below the row says in words where the planet falls in each category (UX.91).

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
| `../src/planetgen/web/` | The pages: routes, templates and helpers (see "Flask pages" below). |
| `../src/planetgen/web/lib/pagination.py` | The site's one pager, used under every paged table (Browse's two tables, Phenomena, a sector's Contents table, a Galaxy Map Quadrant's sector list, each Search result panel, the admin API key list and the admin stats page's duplicate-names list): a "Showing X-Y of Z" summary, then First/Prev, numbered pages and Next/Last, 50 rows a page. Each table has its own page parameter (e.g. `sectors_page`), a plain GET link, and changing a Search filter starts its results back at page 1. |
| `../src/planetgen/web/lib/apiclient.py` | The API client every page calls instead of querying MySQL directly -- one typed wrapper function per read endpoint, plus `auth_*` wrappers (the admin pages) supporting POST/DELETE, a request body, and `Cookie`/`Set-Cookie` relay, and `NotFoundError`/`ApiError` (`web/errors.py` turns these into a 404/502 page; `ApiError.status_code` lets `auth_me`/the admin pages branch on a 401 without string-matching). Inside the Flask app it runs in-process (`web/transport.py`); elsewhere it uses HTTP to `PLANETGEN_API_BASE_URL` (default `http://127.0.0.1/api`). Not web-accessible. |
| `../src/planetgen/web/lib/fmt.py` | HTML-escaping and small formatting helpers (`esc`, `linkify_location`/`nearest_neighbors_location` -- which take a `system_url(system_id)` hook for the neighbour links -- `format_density`, and `static_url`, the versioned `static/` URL every page uses) with nothing to do with fetching data. Not web-accessible. |
| `../src/planetgen/web/maps/galaxymap.py` | Quadrant/Zone classification (`sector_quadrant`/`sector_zone`/`zone_bounds_ly`) shared by the `/galaxy` page's data tables and the sector/sectors pages' "Quadrant N" links back into it -- what's left after this module's former flat-SVG map rendering was superseded by a real 3D scene (`planetgen/web/maps/galaxymap3d.py`). Not web-accessible. |
| `../src/planetgen/web/lib/modulegraph.py` | The modules a map script imports, directly or not (`module_graph('galaxymap3d.js')`), which `galaxy.html` and `sector.html` turn into `<link rel="modulepreload">` hints (MAP.161) so the browser fetches the tree in parallel. |
| `../src/planetgen/web/lib/tilecache.py` | The web layer's on-disk cache (written by the WSGI process serving `/galaxy` and `/galaxy/nodes`) of the 3D Galaxy Map's sector-block nodes (the module keeps the name from when the map was cut into cube tiles): `fetch_nodes` serves each requested node from `<cache dir>/<db>/<generation>/` when present and asks `GET /api/galaxy/nodes` only for the rest. Every 60 s at most it asks `GET /api/galaxy/changes` which nodes changed (from the rows' `modified_at`) and deletes only those, or starts a new generation when the answer is `full`. The changed keys are passed on to the browser's cache too. Size-capped (`cache.max_mb`, default 1024 MB: over it, the files a visit asked for are dropped, the least recently visited first, and the nodes built ahead by the warm-up stay, which stops; PERF.133); location from `PLANETGEN_TILE_CACHE_DIR`/`tile_cache.dir` (see `config.md`). Fails open: an unwritable directory or bad file just means an API call. Not web-accessible. |
| `../src/planetgen/web/lib/pagecache.py` | The pages' cache of the API's public GET answers (`apiclient._request` asks it first; `web.init_app` gives each app its own copy in memory and the answers are shared through Redis, so the warm-up can fill them for every worker). An answer is stored with the data tags it read and served only while those tags' counters (`db/dataversion.py`) are unchanged; an API write, a command-line run and the galaxy's content stamp (checked every `page_cache.stamp_seconds`) bump the tags they touch; capped by age and size (see `config.md`, [cache-warm-up.md](../design/operations/cache-warm-up.md)). Not web-accessible. |
| `../src/planetgen/web/lib/slowpages.py` | The cost rule for slow page requests and what a page sees of them (PERF.105): a route that expects a slow answer runs inside `slowpages.allow_pending()`; a cache miss the rule calls slow raises `slowpages.Pending` carrying the job (`queue/pagejobs.py`) that builds the answer into the page cache, which the page polls at `GET /api/page-jobs/<id>`. Redis down builds it in the request and `take_notice()` says so (see [slow-page-answers.md](../design/operations/slow-page-answers.md)). Not web-accessible. |
| `../src/planetgen/web/lib/slowroute.py` | `@slow_route(title, section, breadcrumbs, description)` on a view whose answers may be slow (UX.118): the view runs inside `slowpages.allow_pending`, and when a job is building the answer the visitor gets the page's frame with a progress bar (`templates/slow_loading.html`, status `202`, never cached or indexed) and `static/slowload.js` swaps the finished page in. `/phenomena` and `/search` use it. Not web-accessible. |
| `../src/planetgen/web/lib/pagetags.py` | What each cached answer reads (`read_tags`) and what each API write changes (`change_for_write`), the two tables behind the page cache's invalidation. Not web-accessible. |
| `../src/planetgen/web/maps/galaxymap3d.py` | Builds the `/galaxy` page's "Galaxy Map" panel: the canvas/controls/info-panel markup, plus the one starting JSON payload (zoom-range numbers, node settings, and the zoomed-all-the-way-out view's nodes, which `initial_node_request` picks the same way the client does for every later view) `static/galaxymap3d.js` reads on first paint -- every later payload, as the camera moves, is fetched by that script directly and never passes through this module. `view_radius_bounds` (the zoom floor/ceiling: a couple of sector-widths up to this galaxy's own real outer edge) is a pure function the `/galaxy` view also calls directly, before this module's own panel-rendering function runs. The panel's JSON names the node URL (`fetchPath`, `/galaxy/nodes`) and a sector-page URL template (`sectorUrl`, from `page_url`) for the info panel's real "View sector" `<a href>`; it carries no database name except `storageKey`, which only namespaces the browser's `localStorage`. Not web-accessible. |
| `../src/planetgen/web/maps/systemmap.py` | Builds the system page's "System Map" panel: a fixed-size, square true-position plot (real angle from `position_x/y_km`, a shared log-radial scale from anchor -- star, or barycenter for a merged/'close' binary pair -- to body) with one colored/sized marker per body, plus a full ring (not a directional band) for each asteroid belt -- distinct from `starmap.py`'s draggable 3D cube, since this only ever needs a flat top-down projection. A once-only pairwise-repulsion pass (`_relax_markers`) nudges apart any two markers real placement happened to put too close together, real position first, decluttering only where needed. Clicking a planet with moons swaps to a "zoom into this planet's moons" scene (`static/systemmap.js` toggles which `<svg>` scene is visible); clicking anything else fills the info side panel, same click-for-info pattern as `starmap.py`/`sectormap.js`. Adds a small green badge to any body with a `life_stage` set. A system's stars beyond the pair (`role` `member`, MAP.290) are drawn as companion markers in the whole-system scene, in the direction of their stored offset at a fixed representative distance and stepped in when two would overlap, each opening a `star-<id>` scene of its own (like a wide pair's secondary) with its own planets and belts. Also hands over each planet/moon marker's resolved class color (`class_color`) plus its `atmosphere`/`composition`/`surface_temperature_k` (whether it has a real atmosphere at all, not just the description text) as `data-*`, plus a star's own spectral-type color (`star_color`, same `data-color` attribute), for `static/systemmap.js`'s own per-marker live sphere rendering. Not web-accessible. |
| `../src/planetgen/web/maps/moonorbits.py` | A moon's orbit about its planet from its stored state vector (MAP.308): `moon_orbit(planet, moon)` returns the osculating ellipse (`physics.state_vectors.elements_from_state`, the two masses for the gravitational parameter) sampled in the map's plane, with its periapsis, apoapsis and eccentricity; a moon with no usable state falls back to the circle of its `distance_km`. `systemmap.py` draws each planet's moon system on one linear scale (`data-pxperkm`, the outermost apoapsis 280 px out) with the moons on those ellipses at their stored positions, at an enlarged body size (`planet-<id>`: one factor for the planet and every moon, the planet 14 to 40 px and never past 40 percent of the way to the nearest moon, the factor said in `data-scale-note`) and at true size (`planet-<id>-true`, the System Map's "True scale" button swaps them; a planet already large against its moons' orbits has the one true-size scene). The planet's object page links `/system/<id>?scene=planet-<id>` when it has moons, and `systemmap.js` opens the scene named in `scene` and keeps it in the URL. Tests: `tests/test_moon_orbits.py`, `tests/test_web_browser_maps.py`. |
| `../src/planetgen/web/lib/classref.py` | The class reference catalog behind `/classes`: every class type built from the generator's own tables (`program_constants`, `physical_constants`, `stellarEvolution.YERKES_CLASS_NAMES`, `nebulaData`, `cometData`), never hand-copied, once per process (`catalog()`, called by `web.init_app` at startup) so the pages always show the values the generator uses. Also `class_url_parts` (is this a real class?), `star_type_classes` ("G2V Yellow Main Sequence Star" -> `("G", "V")`) and `star_type_split` (-> `("G2", "V", " Yellow Main Sequence Star")`, the parts a row links, UX.182). Not web-accessible. |
| `../src/planetgen/web/lib/systempage.py` | The facility pieces the system page and the object pages share (UX.173): `facility_row` (one facility's display values for the Facilities panel) and `facilities_html` (the facilities a body hosts, for its own page). The stars and bodies of a system are the summary tree (`objecttree.py`), not rows built here. Escapes every database value itself; the pages pass the result through `trusted_html`. Not web-accessible. |
| `../src/planetgen/web/lib/tabledisplay.py` | Computes the same "Star Data"/"Planet Data" display strings once baked into the database's now-removed `table_*`/`binary_table_*` columns, but on demand from the raw numeric columns the system page already has -- reuses `planetgen.util.format`'s formatters directly. Not web-accessible. |
| `../src/planetgen/web/maps/starmap.py` | Builds the Sector Map's scene data (`map_scene_data(link_url, ...)`, where `link_url` is `web.helpers.page_url`; every entry carries a plain `href`), which the sector page answers at `/sector/<id>/scene`: computes every position/size/color/label the map needs (a wedge or fallback-cube outline, one entry per star -- a small core (`_star_dot_radius`: log-scaled from `radius_km` but kept between 1.5 and 6 scene units, so even a supergiant is a point) inside a glow shell whose size, strength and fade come from `luminosity_w` (`_star_glow`: a supergiant gets a big soft halo, a white dwarf almost none, matching how the Galaxy Map draws bright stars; `sectorscene.js` clicks a star through an unseen sphere at least 4 units across), color from `star_type`'s spectral letter (`SPECTRAL_CLASS_COLORS`) shaded by `luminosity_w` and nudged by where `temperature_k` falls in that spectral class's range, so "White Giant" reads white and "Blue Giant" reads blue regardless of temperature -- and one per nearby standalone phenomenon, sized by `radius_ly` (always 0 for the point-like types) and positioned directly in the galaxy frame, no rotation needed unlike a star system's sector-local position -- see `queryDb.phenomena_near_sector`) and returns it as a dict; `static/sectorscene.js` is what actually renders it, this module builds no HTML of its own. Not web-accessible. |
| `../src/planetgen/web/maps/navmap.py` | Builds the NAV page's "NAV Map" panel (`render_nav_map_panel(link_url, ...)`; each point is an SVG `<a href>`): a flat, static, top-down SVG plot of the galactic X-Y plane -- origin and destination as labeled points, a dashed line for the direct course, and (when one was found) a solid polyline through the optimal route's intermediate hops. Auto-scaled to whatever points it's given (no fixed sector size to normalize against), with one uniform light-years-per-pixel ratio on both axes so bearings aren't visually distorted, plus a compass arrow along the origin's bearing 000 (toward the frame's center) and a scale-bar legend. Deliberately blind to altitude/z, same as the flat SVG phenomenon Diagram panel (`planetgen/web/maps/phenomenonmap.py`) -- the course's mark already covers that axis. Not web-accessible. |
| `../src/planetgen/web/maps/phenomenonmap.py` | Builds the phenomenon page's flat, zoomable SVG diagram of a nebula's or supernova remnant's real extent, drawn in astronomical units against an AU-scale yardstick (`static/phenomenonmap.js` adds the zoom and pan). Not web-accessible. |
| `../src/planetgen/web/maps/phenomenonrender.py` | Builds the phenomenon page's "View" panel for a neutron star, black hole, quasar, rogue planet or interstellar comet: the numbers `static/phenomenonrender.js` draws with three.js, plus a static SVG used when WebGL can't start. Not web-accessible. |
| `../src/planetgen/web/warmup.py` | The site's caches built ahead of the first visit (MAP.134, PERF.104): `opening_request` names the nodes the first frame of `/galaxy` asks for (the page uses it too), `warm_opening_view` fetches them and the galaxy stage into the node cache, and `warm_family` visits the pages of one family (`FAMILIES`: galaxy-map, lists, sectors, systems, phenomena, search) through the site's own app so the answers land in the shared page cache; each family is a command (`warm-<family>`) whose unit is `warm_<family>_unit` here (PERF.127). `planetgen warm` builds the list of those commands and the wrapper runs them; `update.sh` starts it in the background after every update (at most 30 minutes on Linux, never waited for), and it is the command to run after clearing a cache by hand. Not web-accessible. |
| `../src/planetgen/web/lib/datatable.py`, `../src/planetgen/web/tables.py` | The site's shared data table (UX.41): a page module describes its table once as a `Table`, `register`s it and renders it; the first page is rendered on the server (sort links, a filter form, the shared pager), and `/table/<name>` answers the next pages as JSON for `static/datatable.js`. |
| `../src/planetgen/web/lib/tilepack.py` | The packed form of a Galaxy Map node (MAP.159, MAP.260) the node cache stores and `/galaxy/nodes` sends: stars in columns (about 19 bytes each, against about 150 as JSON), the clouds, points and `incomplete` list as JSON, and the names of the stars that have one as `names` (by id; MAP.298, `TILE_FORMAT` 4). `static/tilepack.js` reads it. |
| `../src/planetgen/web/lib/singleflight.py` | One build of each Galaxy Map node or stage at a time (PERF.38): a request claims the nodes it is about to fetch with a Redis `SET NX EX`, waits for ones another process claimed, and builds only those whose claimant never delivered. With no Redis it fails open and everything is built. |
| `../src/planetgen/web/lib/estimatecache.py` | The Generate page's size-and-time estimate (`planetgen galaxy ... --estimate-only`, up to two minutes) kept in Redis for five minutes under the job kind, a digest of its command line and the database (PERF.109). It is served only while the generator epoch, the galaxy's data version (`db/dataversion.py`) and the speed statistics (`generation/stats.py` `fingerprint`) are as they were, so a fill, an edit, a new epoch or newly recorded speeds make the next estimate a new one. Two admins asking for the same estimate at once run it once (`singleflight`); with no Redis it just runs. |
| `../src/planetgen/web/lib/habitability_ref.py`, `../src/planetgen/web/lib/habitability_symbols.py`, `../src/planetgen/web/lib/facility_display.py` | The text of the habitability page (`/help/habitability`, UX.90, UX.180; the ids its symbols, life scores and chips link to are UX.181), the PHI-6 symbols every page draws through `habitability_symbol(domain, tier)` (UX.101: seven solid shapes in `static/icons.svg` with their insides cut out, no border, in the colour of the `tier-blue` to `tier-red` class at 0.8 opacity so the page shades it in both themes, with the domain and colour in their words; the outline of every other icon is the `.icon` rule of `style.css`, as the sprite is served under a policy that refuses an inline style in it), and how a facility's type, affiliation, ratings and field values read (POP.10: rating chips carry a shape, the dimension's name and its word, so colour is never the only signal). |
| `../src/planetgen/web/lib/life_help.py`, `../src/planetgen/web/templates/life_help.html` | The text of the life help pages (`/help/life`, DOC.84): the kinds of life, the four scales, the stages and each life class, written from `physics/life_readings.py`, `physics/life_classes.py`, `physics/life_class_text.py`, the ledger `physics/data/life_ledger.json` and `tuning.py`, so a page cannot disagree with the score it explains. A planet page links to them: the class names, a "?" beside the Kind, How much, Trend and Security headings, "What these words mean" under the At a glance paragraph, and the stages, suitability and classes in the Life section. |
| `../src/planetgen/web/lib/poi_icons.py`, `../src/html/static/poi/*.svg`, `../src/html/static/poiicons.js`, `../src/planetgen/web/templates/partials/poi_key.html` | The point-of-interest symbols (MAP.320): seven solid silhouettes, one SVG file each under `static/poi/` (habitable world, uninhabitable world, inhabited world, gas giant, asteroid field, sector, region), each a bare `<path>` set on a 24 by 24 viewBox with its name in its `<title>` and no paint of its own (a policy that refuses inline styles would drop it, as with `icons.svg`). A page draws one through `poi_icon(id, tier)` as `<svg class="poi-icon poi-<id>" role="img" aria-label=...><use href="poi/<id>.svg#poi">`; `style.css` paints it (`.poi-icon` fill and edge, the `--poi-brown`, `--poi-grey` and `--poi-gas` tokens, the `poi-tier-0` to `poi-tier-3` classes of the habitability colours). The habitable and the inhabited world take the habitability colour of their score; the uninhabitable world is brown and never red, the asteroid field grey, and brown and grey differ in lightness from all four tiers (`tests/test_poi_icons.py`, `tests/test_poi_icons_browser.py`). `static/poiicons.js` lists the same symbols for the map (a marker kind with `icon` instead of `glyph`; `galaxymarkerlayer.js` reads the file into a texture through `shapesOf` and `Path2D`, as the policy takes no image from a data: or blob: URL; `maplayers.js` draws a legend line with `icon` and `token`), and `partials/poi_key.html` is the key, on the Galaxy Map under the map and on the habitability page, to include on any page with `with context`. The data that says which symbol a system or body gets (the point-of-interest registry, MAP.305) is not built yet; only the habitable-world marker uses its symbol so far. |
| `../src/html/static/galaxymarkers.js`, `../src/html/static/galaxymarkerlayer.js` | The Galaxy Map's marker layer (MAP.219): the registry of marker kinds (black holes; one kind per interest reason, `INTEREST_REASONS` (the planet classes and, from MAP.351, the 15 life classes `LIFE_CLASS_IDS`, each its own Menu group), the symbols of MAP.320), `selectMarkers` (the kinds that are on, the items in view, within the budget of 240, a system once under its first reason that is on), `mergeMarks` (MAP.305: star-level marks closer than `MERGE_PX` = 22 CSS pixels on the screen merge into one mark with the best reason's symbol, a count badge and a card that lists every system and reason; a black hole, a body's symbol and a region's token never merge; zooming in splits the crowd), a kind pulled forward (MAP.306: first call on the budget and on a system's symbol, drawn over the other marks, its names first) and the sprite layer that draws them. No DOM in `galaxymarkers.js`; `tests/js/galaxymarkers.test.mjs` holds it. Served directly. Black holes are chosen for the view (MAP.339, `nearFieldRule`): one within 3 view distances of the point the camera looks at, or beyond it only if its mass is at least 5 solar masses times (distance over that reach) squared, ranked by mass over distance squared, at most one per 120 by 120 pixels of the screen (the `map.blackHoleDensity` preference), with no exemption for a heavy hole. Crowded black holes merge like crowded systems (MAP.340): holes within 22 px of each other and about the same distance from the camera (`MERGE_DEPTH_RATIO`, 1.5) are one mark with a count badge (the heaviest first), whose card lists them, each a button that shows that hole and flies to it, and has a Go to them button; they never merge with a system's mark. With a sector opened in place every mark outside it is drawn at 0.3 of its strength, as the stars are (MAP.341, `setOutsideShare`; the open system's own bodies keep all of it). |
| `../src/html/static/mappreferences.js` | The map's per-browser settings (MAP.339, MAP.336): `readNumber` and `blackHoleDensityPx` read `localStorage` `planetgen.pref.<key>` (the number's text) and fall back to the default when the key is missing, unreadable or out of range (`map.blackHoleDensity`: default 120, the side in pixels of the square one black hole mark may stand on, 30 to 1000). `labelBufferPx` (`map.labelBuffer`, the label engine's margin, MAP.299), `readChoice` and `regionLanding` read a setting that is one of a few words (`map.regionLanding`: `frame`, the default, or `densest`, where Go to a region lands, MAP.336). The preferences screen that sets them is UX.236. No DOM and no three.js. UX.236 adds the other half: `PREFERENCE_TEXT` (each setting's group, name and plain-words description, and the words of a choice), `preferenceList()` (the settings in the order the page lists them), `readSetting`, `writeSetting` (refuses a word or number the reader would not accept), `clearSetting` and `resetSettings`; `preferencespage.js` draws `/preferences` from them. |
| `../src/html/static/regionlanding.js` | Where Go to a region lands (MAP.336), pure: `densestCluster(sectors, region, edgePc)` takes the sectors with their total mass (`stats.mass_kg`, the ranking MAP.322 names labels by), the sphere the region's extent bounds and a sector's edge, and returns the sphere `{center, radius, count, massKg}` of the heaviest group of the region's nine heaviest sectors (within 2.5 edges of one of them, about its mass-weighted middle), or null when nothing is ranked (the landing then stays framed on the region). Tests: `tests/js/regionlanding.test.mjs`. |
| `../src/html/static/findplace.js`, `../src/html/static/findplacebox.js` | The Galaxy Map's find-a-place box (MAP.336). `findplace.js` is the choosing half, pure: what a typed line can be, as groups of results each holding a target for the map's one `goTo` (coordinates `x y z` in pc, a sector designation or shell/band/wedge, an object ID with or without its kind (its first group is the sector, the first two with `-000` the system, a body's `<kind>:<id>` picks the body), charted regions by name, the server's name matches, the places visited before, and the way to the whole galaxy); `createRecents(storage)` is the visit log (`localStorage` `planetgen.recentPlaces`, the last eight places gone to by name, each once, newest first). `findplacebox.js` fills the form: results as you type (the server asked once per line, after a short pause), grouped under headings in a listbox, Up, Down, Enter and Escape, a click, and `/` anywhere outside a field reaching the box; an empty box lists the places visited and the way home. Tests: `tests/js/findplace.test.mjs`, `tests/test_web_browser_go_to.py`. |
| `../src/html/static/maplogruler.js` | The pure half of the Galaxy Map's log ruler (MAP.336): the named scales (Galaxy, Arm and halo 1 kpc, Region 40 pc, Sector 4 pc, Neighbourhood 1 pc, Planetary system 1e-3 pc, Planet 1e-8 pc, Moon 1e-10 pc), `stopsBetween(farPc, floorPc)` with each one's share along the ruler (top 0, the whole galaxy; bottom 1, a moon), `distanceAt` and `shareOf` (the ruler's log scale) and `readoutText` ("Distance from focus: 2.4 pc (7.8 ly); nearest: a star, 0.9 pc", in the reader's units). Tests: `tests/js/maplogruler.test.mjs`. |
| `../src/html/static/poiprefs.js` | The Galaxy Map's point-of-interest choices a browser remembers (MAP.306): each Menu row's state (shown, pulled forward or hidden) is kept in `localStorage` and the next map opens that way, unless its address names a layer (`hide`, `show` or `front`), which wins. The layer registry (`maplayers.js`) holds the third state: a layer with `forward: true` gets a second button, `front=` in the address lists the layers pulled forward, and a group with `bulkGroup` gets a row that sets all its layers at once. No DOM; `tests/js/poiprefs.test.mjs` holds it. Served directly. |
| `../src/html/static/maplayers.js` | The layer registry the Galaxy, Sector and System Maps share (MAP.217): a toggled layer is declared once and its Menu button, URL key, legend line and picker layer come from it. MAP.234: the on-map **Key** is the registry's legend (`renderLegend(container, {collapsible: true})`, a fold the reader can close, kept closed in `localStorage`) plus key lines (`registerKey`: text, swatch, `drawn()`) that show only while their thing is drawn and follow the view (`refreshKey()`); a layer with no `legend`, no key line pointing at it (`layer`) and no `keyless` reason is named by `unkeyed()`, which `tests/test_web_browser_map_key.py` holds empty, so a new layer cannot be added without a line. The help dialog's Key section (`GALAXY_MAP_HELP`) says the same in words. Served directly, same as `style.css`. |
| `../src/html/static/bodymarkers.js` | The symbols of the bodies themselves on the Galaxy Map's system level (MAP.321): which interest reason a body of the system scene is (`reasonOfBody`: inhabited, gas giant, habitable, uninhabitable; a belt is an asteroid field), the marker items for the open system's bodies by the same Menu layers as the stars' marks (`bodyMarkers`, each riding its body through `item.follow`), and the hand-off between the two levels (`createHandoff`: the symbols pass from the star to the bodies within 6 radii of the system and back beyond 10, so a camera resting between does not flicker). No DOM and no three.js; `tests/js/bodymarkers.test.mjs` holds it. |
| `../src/planetgen/web/nearby_page.py`, `population_pages.py`, `system_species.py`, `facility_type_pages.py`, `queue_page.py`, `settings_page.py`, `backups_page.py`, `images_page.py`, `seo.py`, `seo_pages.py`, `assets.py`, `maintenance_page.py`, `maintenance_gate.py`, `job_runner.py` | The What's nearby page (NAV.44); the species and polity pages (POP.1 to POP.4); the system page's "Species that could live here" panel (POP.12); the admin Facility types page (POP.10); the admin work queue page (ADM.10); the admin settings page (ADM.84); the backup archives screen and the settings export and import (ADM.85); the Images screen, the page `<head>`, the indexing matrix, robots.txt, the sitemaps, the small site files and the icon and picture store (ADM.44); the `/maintenance` page and the `before_request` hook that sends everyone but an admin to it while maintenance mode is on (UX.112); and the RQ job that runs one Generate page job's steps in order, writing `output.log` and `state.json` in the job's directory (PERF.24). |
| `../src/planetgen/web/lib/privatedir.py` | The private fallback directory (mode 0700, in the system temp directory) the node cache and the Generate jobs use when their configured directory can't be created. A directory found there that isn't a real directory owned by this user, or that others can write to, is refused rather than reused. Not web-accessible. |
| `../src/html/static/style.css` | Shared stylesheet (CSS custom properties, light/dark via `prefers-color-scheme` or an explicit `data-theme` on `<html>`, card-style panels, phone layout under 40rem), served directly by the web server. |
| `../src/html/static/theme.js` | Loaded on every page, blocking, before `style.css` (`web/templates/base.html`): applies the saved light/dark/system theme before the first paint and drives the header's theme button. See "The page shell" below. |
| `../src/html/static/favicon.svg` | The site icon (a small ringed planet), linked from every page's `<head>`. |
| `../src/html/static/vendor/` | Vendored third-party JS -- currently just `three.module.min.js` (three.js, bundled+minified from the `three` npm package), used by `sectorscene.js`/`systemmap.js`/`galaxymap3d.js`. Vendored rather than loaded from a CDN so the pages' `Content-Security-Policy: default-src 'self'` needs no exception; see this directory's own `THIRD_PARTY_NOTICES.txt` for the license and how to rebuild it from a newer release. Served directly, same as `style.css`. |
| `../src/html/static/sectorscene.js` | One sector's scene (MAP.66): `buildSectorScene(data, options)` turns `starmap.py`'s JSON into a `THREE.Group` (GPU-billboarded points of light for every star and light-giving phenomenon, textured spheres for the rest, translucent volumes for nebulae, neighbor markers, the cell outline and compass) with the picker layers (`mappick.js`), info-panel specs and tooltip text for what is in it. `options.origin`, `unit` and `flipY` place it in a map's world (the Galaxy Map's parsec frame); with none it is the Sector Map's own scene units. The sector page's map and the Galaxy Map's opened sector both draw it. Each kind of object in it (stars, nebulae, supernova remnants, asteroid fields, black holes, neutron stars, quasars, rogue planets, interstellar comets, neighboring sectors) can be hidden (`setKindHidden`, MAP.79): its points, bodies, rings and labels go, and a hidden kind can't be hovered or picked. A supernova remnant's shell is painted from its age (`remnantShell`, MAP.286: a young remnant a thin bright hot-coloured rim, an old one a thick faint cooler shell; the card names the age). Served directly, same as `style.css`. |
| `../src/html/static/galaxysector.js` | The Galaxy Map's last drill-down stage (MAP.66): `createSectorStage(host).open(id, bounds)` fetches a sector's scene (`GET /sector/<id>/scene`, `web/sector_page.py`: the same JSON the sector page embeds, plus `centerPc` and `halfEdgePc`), builds it with `sectorscene.js` at its place in the galaxy, adds its picker layers (with tooltip, hover ring and info-panel select) to the map's picker, and has the galaxy's own tile stars and phenomena left out of its cell; `close()` takes it all away. `galaxystageview.js` decides when a sector is open, flies the camera to it and keeps `?sector=<designation>&open=1` and the map's Back and Forward. Served directly, same as `style.css`. |
| `../src/html/static/generatebuttons.js` | The admin Generate buttons for a sector that isn't generated yet (it alone, its neighborhood, its column, or its whole shell after a confirm), shared by the sector page's map and the Galaxy Map (`galaxymap3d.js`, a sector cell or a neighbor). Each is a plain POST form to the Generate page with the CSRF token from the server's `generate` target (`web.helpers.generate_target`), which the pages only send a logged-in admin. Served directly, same as `style.css`. |
| `../src/html/static/starfade.js` | The Galaxy Map's fade of stars and other objects (MAP.153, MAP.155, MAP.152), free of the page. A star shows at the share of the node that lists it: a node's weight is what the hand-over of its ancestors leaves it (`galaxynodecut.js`'s `weightCalculator`, which the page asks once a frame), so a node and its children always add up to one and a star listed in both stays whole (`starShare`). The page works the shares out as the camera moves and hands them to the star shader; black holes, neutron stars and quasars are shared the same way; a nebula sprite's opacity ramps between 1 and 4 px of radius (`cloudSizeShare`); a star arriving is held as a point of light until its zoom share rises (the shader, with `static/starpoint.js`). |
| `../src/html/static/starpoint.js` | MAP.208: how big a star, black hole, neutron star or quasar is drawn: a point of light 1.5 px wide until the camera is about a sector away (`NEAR_SECTOR_PC`, one constant the Galaxy Map and the sector scene share), growing smoothly over the last two zoom steps (`pointGrow`, `pointSizePx`, twinned in GLSL as `STAR_POINT_GLSL`), the far brightness (`pointBrightness`) and the cap that keeps a phenomenon no brighter than the brightest star (`phenomenonLuminositySol`). |
| `../src/html/static/galaxynodecut.js` | The nodes the Galaxy Map fetches and draws (MAP.261): the cut by projected size, a block opening into its children as it grows on screen, inside the cone the camera looks along and a budget of nodes, nearest first; `weightCalculator` gives each node its weight (a block gives way to its children between 70 and 140 px of projected size) and `parseNodeKey`, `formatNodeKey` the key `m.shell.band.wedge`. The same arithmetic as `planetgen.galaxy.nodecut.node_cut`; `tests/test_js_unit.py` (`tests/js/galaxynodecut.test.mjs`) and `tests/test_nodecut.py` hold the two equal. |
| `../src/html/static/starbuffer.js` | The Galaxy Map's one star buffer (MAP.267), free of the page: relative-to-eye pairs (`splitPair`: each coordinate as a high and a low float32, the shader subtracting the camera's own pair, so a star keeps under a quarter pixel down to about 0.2 AU from the camera at any distance from the origin; `relativeOffset` is that sum as the GPU rounds it) and `SlotStore`, the dense slots the node-by-node fill uses (a node's stars are a chunk at the end; a node leaving is cut out by moving the tail into its place). `tests/js/starbuffer.test.mjs` measures the precision and holds the store. |
| `../src/html/static/galaxyglow.js` | The Galaxy Map's unresolved glow, the choosing half (MAP.266): a node's blob is the light its aggregate stores (`aggregate.lights`, `aggregate.stats.luminosity_sol`) less the stars the node lists, at the centre of that light; `glowLook` gives its strength from its apparent magnitude (judged `GLOW_MAG` deeper than a point) and its radius from the block's size on screen. No DOM and no three.js; `tests/js/galaxyglow.test.mjs` holds it. The glow fades as its block grows on the screen (`GLOW_RESOLVE_PX`: from 30 px radius, on a squared curve, gone at 130 px, where the halo would cover much of the screen over the block's own stars, MAP.330) and is capped at `GLOW_ALPHA_MAX` 0.4 (MAP.324). |
| `../src/html/static/galaxyglowlayer.js` | The glow drawn (MAP.266): one screen-facing quad per blob, sized in pixels in the vertex shader, in the color of the light's temperature, under the stars; each blob shows for the share the cut gives its node, so a node and its children trade places without a pop. A layer of the Menu (`glow`, "Unresolved glow", on by default, `hide=glow` in the address). The Sector Map has none. Overlapping glows are drawn with a maximum blend, not an additive one, so they never add up to white (MAP.324); `strongest()` reports the strongest glow that reaches the screen, for the browser tests. The same layer draws the charted regions (MAP.316) with its `look` option (`regionLook`) and a blob's own `color`, one instance under the glow's: a soft sphere at each region's extent, a layer of the Menu (`regions`, "Charted regions", on by default, `hide=regions` in the address), left out while a sector is opened. A click on a region's sphere shows its card (how much is charted, a Fly to region button) and a double click flies to it, after any star, marker or object under the pointer. |
| `../src/html/static/galaxyregions.js` | The charted regions of the Galaxy Map, the choosing half (MAP.316): `regionBlobs` turns the API's `charted-regions` list into blobs (the extent sphere, the share's color from cool blue for a lone sector through teal to gold and amber, the name), `regionLook` gives each sphere's strength (`REGION_ALPHA` 0.2 to 0.5 by the share's square root) and radius on screen (`REGION_RADIUS_PX` 5 to 280 px), faded as its stars resolve (`REGION_RESOLVE_PX`: from 160 px radius, gone at 600 px), `sortRegions` orders the list by share charted or by name, `regionMarkers` gives each region with an inhabited world one marker item (the inhabited symbol at its middle, in its best world's tier) while its sphere is drawn (`REGION_TOKEN_MIN_ALPHA`) and no system's own mark stands for it (MAP.321). The names on the map are `regionlabels.js` (MAP.314). No DOM and no three.js; `tests/js/galaxyregions.test.mjs` holds it. Served directly. |
| `../src/html/static/galaxynebulae.js` | The nebula regions of the Galaxy Map, the choosing half (MAP.176; `galaxyglowlayer.js` draws them): the blobs of the groups each drawn node lists in `nebula_regions` (MAP.173), a soft sprite in the colour of the group's dominant family (`nebulalook.js`, a remnant orange) with a thin ring when it holds more than one nebula, drawn from 2 px across and gone as it fills the screen, at most 600 at once (the strongest on the screen, the last tenth fading), each at the share the cut gives its node. Six layers of the Menu under "Nebulae", one for each family (`nebulae-<family>`, on by default, `hide=nebulae-dark` in the address), and a legend line that says what a ring is. No DOM and no three.js. |
| `../src/html/static/galaxyregionlist.js` | The Menu's charted-regions list (MAP.316): a button a region ("Fardoo Beez Region, 14 percent charted") that flies to it, sorted by share charted or by name; hidden while no region is charted. Served directly. |
| `../src/html/static/sectorframe.js` | The frame of the sector opened in place on the Galaxy Map (MAP.327): a wireframe round its cell that stays where the cell is (it does not follow the pointer; in the galactic view the same wireframe stands round the charted region the view is in, with no icons and nothing dimmed), and at the centre of each wall that faces another sector a labelled screen-sized icon (arrow, direction, the neighbour's name) whose click goes to that sector (`stageview.openSector`: the camera goes to its window and it opens in place, the address `?sector=...&open=1` naming it). `wallCentres` is the pure geometry (a wall on the galaxy's own boundary, at the centre, past a pole or all the way round in longitude, has no icon); the neighbour is the sector holding the point a quarter of an edge beyond the wall. While a sector is open the rest of the map is dimmed to 30 percent (`outsideDim` in the star shader), its glow, names, comets, rogue planets and black holes are not drawn. Unit tests `tests/js/sectorframe.test.mjs`. Served directly. The region frame the galactic view draws has a legend line, "Outline: <region name> (9 sectors a side)", on the Charted regions layer while it shows (MAP.331). |
| `../src/html/static/systemstep.js` | Stepping from the open system to a neighbouring one on the Galaxy Map (MAP.121): `directionAt` is the unit vector of north, south, east, west, out or in in the galaxy's frame at a point (null on the polar axis and at the centre for the sideways ones), and `nearestInDirection` is the nearest candidate within a 60 degree cone about it (`CONE_HALF_ANGLE_DEG`, a display tweak with no source) or null. The step buttons and keys use it through `stageview.stepToNeighbour` on the open sector's systems and, across the sector's edge, on the systems of the sectors round it (`crossstars.js`, MAP.329). Unit tests `tests/js/systemstep.test.mjs`. Served directly. |
| `../src/html/static/crossstars.js` | The systems of the sectors round the open one (MAP.329), pure: `neighbourRecords` chooses which generated sectors (from `/galaxy/near`) to read, the nearest `CROSS_SECTORS` (8) whose middle is within `CROSS_EDGES` (1.5) sector edges of the system stood at, never the open one; `starsOfScene` reads a sector's scene JSON (`/sector/<id>/scene`) for each system's place (the scene's own units turned to parsecs about the sector's middle, as `sectorscene.js` places them), name, `system:` endpoint and page, a binary's companion folded into its primary. `stageview` reads them once for each system opened and the step buttons and keys then cross the edge: the sector beside opens in place and its system with it (each its own history entry, the scale kept). Nothing of the neighbours is drawn; the faint context of MAP.141 is not built. Unit tests `tests/js/crossstars.test.mjs`, browser test `test_a_step_across_the_sector_edge_opens_the_system_beside_in_its_own_sector`. Served directly. |
| `../src/html/static/starmagnitude.js` | The Galaxy Map's star visibility law (MAP.148), free of the page so the shader, the page's histogram and the tests share one definition: a star's opacity follows its apparent magnitude from the camera (`M + 5 log10(d / 10 pc)`) against a limiting magnitude with a 1.5 magnitude ramp, its brightness follows tone-mapped flux, and the limit is the larger of the HISTOGRAM limit (where the stars on screen, counted by the MAP.153 zoom opacity, reach 20,000, or 8,000 on a phone) and the CALIBRATION limit (where the dimmest stars the tiles carry still show whole at the focus distance), so it only dims what stands much farther than the thing looked at and never thins the same stars twice with the rank birth radius. `galaxymap3d.js` recounts at most every 150 ms when the view changes and the shader eases toward the new limit over about 0.3 s. Black holes, neutron stars and quasars are never dimmed away. The tile floors (`generated_star_floor_sol`) stay as the server's superset until tiles are chosen by distance. |
| `../src/html/static/nebulalook.js` | How a nebula is painted on the Galaxy Map, per theme (MAP.175): [colour, core opacity, edge opacity] for diffuse, emission, reflection, planetary, dark and unknown nebulae, one set for the light theme and one for the dark, each reaching 3:1 (`CONTRAST_TARGET`) against the map's background; the dark family is a dust colour on the dark theme. Also the WCAG contrast arithmetic the test uses, against the real colours in `style.css`. MAP.181 settles the final dust colour. |
| `../src/html/static/nearfield.js` | The Galaxy Map's near field (MAP.149), free of the page so the shaders, the picker and the tests share one definition: a depth fade (things nearer the camera than 0.15 of the focus distance are gone, from 0.5 on whole), a see-through tube (down to 0.08) round the line from the camera to the focus thinning what stands in front of it, the faint-region rules (`contextOpacity` 0.08-0.3, `prominence`, `limitingMagnitude`) and the polar `blockContains` test for the block the camera is inside. Everything is worked out in view space, and `NEAR_FIELD_GLSL` is the shaders' twin built from the same constants. `galaxymap3d.js` draws blocks and stars with it, and `galaxystageview.js` picks with it: a block whose share is under 0.35 at the hit point, or a back face, is skipped so the pick falls through to what is behind. |
| `../src/html/static/galaxymap3d.js` | Renders `planetgen/web/maps/galaxymap3d.py`'s `#galaxymap3d-data` JSON as a real WebGL scene (three.js, vendored at `static/vendor/`) -- this camera's own orbit target moves freely through the galaxy, so most of what it draws is fetched live from `/galaxy/nodes` (plain GET, no database in the URL) (debounced, on every camera move), one sector block (node) at a time, rather than server-rendered once. Nodes are kept in memory and in the browser's IndexedDB (MAP.158: database `planetgen-nodes`, store `nodes`, keyed by database, generation and node key, at most about 96 MB with the least recently used dropped; only the small stamp record stays in `localStorage` under `planetgen:node:<db>:@stamp`), so panning back over seen space or reloading the page doesn't refetch them. Once the view's nodes are in and the camera has been still for 1.5 s, the nodes one zoom step in are fetched at low priority (marked `X-Map-Prefetch`) (not the way out, not with the browser's data saver on or a connection of 3G or slower, at most about 1.5 MB a page visit). The map is driven by the drill-down; the whole galaxy is a 3D disk with no grid lines drawn on it, and in every view drag turns it, right-drag or Shift-drag moves it, the wheel zooms toward the point under the cursor (that point stays under it, MAP.150) and a pinch zooms about the middle, both within one range in every stage (MAP.334: from the whole galaxy in view down to 0.001 pc outside an opened system and 1e-10 pc in one; a stage's own ends are only where it hands over to the next, and a step covers a share of the distance to the nearest thing along the pointer's ray; MAP.335: a sector's scene, and the system's out of it, come up because the eye is near them, with the generated sectors near the eye asked of `/galaxy/near` in any stage, not read from a window of sectors) (Re-center, in the Menu, brings it back; the how-to text is behind Map help in the Menu, UX.50) (`static/galaxystageview.js`, rules in `static/galaxystages.js`): the whole galaxy is a plain disk of blocks of 243 sectors (hovering one outlines it); a click on a block opens a window of nine blocks (the "size" `g`, 3 to 11 in the URL as `g=`) centred on the block clicked, one size down (243, 81, 27, 9, 3, then single sectors, MAP.151b), the sectors in a round region of the block lattice with the cut-off drawn as the window's own edge; Up recentres on the parent block and ends at the galaxy; there are no arc, slab or segment picks and no slab buttons (MAP.150, MAP.151d) (design doc `docs/design/interface/galaxy-drilldown-navigation.md`, sections 4-5). The blocks (`static/galaxyprisms.js`, built by `static/galaxyblocks.js`) are the picking geometry and are not drawn by default (MAP.265: no translucent prisms): they show only through the Menu's Cell grid, a Color by fill (colored by predicted density, generated blocks amber and more solid with their filled share) or Charted only. A click on a block also selects the sector holding the middle of what the click went through, found in closed form (`cellAt`) and so selectable from any distance; once the window arrives that sector is outlined and described in the panel (its record asked for as the smallest window round it). Hovering a choice dims the others and outlines it along its blocks' own sides ; Back, Forward, Up and Whole galaxy buttons and the round Steps menu (every step so far, and Jump to newest) step through the stages, each of which has its own URL, beside Bookmarks and a Menu holding Re-center, Cell grid, Charted only, Territories and Map help. A block's panel gives its ring, layer and slot ranges, exact sector count and generated count; an unfilled sector's gives its address and designation, plus, for a logged-in admin, the Sector Map's Generate buttons (`static/generatebuttons.js`: this sector, its neighborhood, its column, or its whole shell after a confirm), which post to the Generate page. Nebulae and supernova remnants are drawn as soft translucent spheres their real size (each node lists the ones reaching into it, `queryDb.galaxy_clouds_in_box`), hidden when only a pixel or two across and faded out as the camera nears them; clicking one while it is small enough to aim at shows its type, class and radius and links to its phenomenon page. Moves between stages fly along a smooth zoom-and-pan path (they jump with `prefers-reduced-motion`), and while idle the page fetches the nodes one zoom step in. A one-line scale readout gives a bar's length in sectors, pc and ly. Served directly, same as `style.css`. MAP.349: every line of a card that names a container or a place has its own **Go to** (`goToButton`, the map's one `goTo`): the "Region:" and "In:" lines (`regionNode`, on the card of anything inside a charted region and on a region's own card), each system of the nearest-systems list (by its printed ID, as the find box does) and each member of a crowd. |
| `../src/html/static/maplabels.js` | The map label engine (MAP.229, built under MAP.268): `placeLabels` projects each candidate (text, 3D anchor, kind, priority, optional on-screen size range), drops those outside the view or outside their size range, and places the rest by priority on a screen grid so none overlap, up to a budget (60), with hysteresis (a label already shown keeps its place against an equal one and its size range widens by 15%); `createLabelLayer` draws them as pooled `<span class="map-label">` elements in an overlay with a halo that reads on both themes, a short fade (none under reduced motion), `pointer-events: none` (a label never takes a pick or hover itself: `at(x, y)` answers which shown label lies under a screen point, so the map's own pointer code can make a name clickable, MAP.298, and `setHot` marks the one the pointer is over) and `textContent` only (names are database content). `placeBoxes` places labels that already have a preferred box (the System Map diagram's names, MAP.231): a box slides back inside the bounds, steps through its vertical shifts when it overlaps a kept box or a marker circle, and is hidden when none clears. From MAP.299 a name keeps a clear buffer (4 px by default, the `map.labelBuffer` setting of `mappreferences.js`, 0 to 40) from every other, measured on what is drawn: the screen grid reaches each box's buffer so two names in neighbouring cells are still compared, a box is on whole pixels and as wide as its text, and the text's width is read in the font of its own kind (`createKindMeasure`, from a hidden probe in the overlay: a sector's name is semibold, a nebula's italic). Where names of one kind compete for room the heavier keeps the better place (`weightedPriorities`: the kind's priorities dealt out again in weight order, a star's by brightness, an object's or body's by mass, a nebula's by volume), so a cut falls on the lower brightness or mass first. Unit tests `tests/js/maplabels.test.mjs`. Served directly. |
| `../src/html/static/galaxylabels.js` | The Galaxy Map's names (MAP.268): sectors (a generated sector in a window close enough to read, by name), systems (the named stars of the sector open or being come up to, ranked by apparent brightness), bodies (the planets and moons of an open system once a few pixels across, so a planet is named from the sector's far view to the planet itself) nebulae (name, else type, where the nebula is large enough on screen), the named stars the map draws (`label-star`) and the other named objects (`label-object`: comets, rogue planets, black holes, neutron stars), the stars and objects (MAP.298), the inhabited and habitable world markers that are drawn (`label-world`, MAP.317, MAP.321), a planet or moon whose symbol is drawn taking its place in the inhabited slice too, the names of charted regions (`label-region`, MAP.314), all offered through `labelbudget.js`, through `maplabels.js`. A click on a name selects its item and flies to it, as a click and a double click on the item do. Each kind is a layer of the shared registry in the Menu's Labels group (`label-sector`, `label-system`, `label-body`, `label-nebula`, `label-star`, `label-object`, `label-region`; the address's `hide` list keeps them). The Sector Map's and System Map's own labels are not on the engine yet. Served directly. |
| `../src/html/static/regionlabels.js` | The names of charted regions on the Galaxy Map by zoom (MAP.314), without DOM or three.js. From the nodes the view is drawn with (a block that has handed over to its children is not asked): an 81- or 243-block names its family with Reach, once, the family holding the most generated sectors, at the centroid of its charted regions' spheres weighted by their sectors; a 27-block names each charted region it lists and a 9-block its own; blocks of 3 sectors and below name nothing (the info card of MAP.315 says it there). A block with no charted region offers nothing, so an unfilled galaxy shows no name. Every name goes first through `labelbudget.js` (domain `region`, ranked by the share charted) and a click flies to its block (framing a region by its extent, MAP.319). The server puts the lists on the nodes (`regions` on 27-, 81- and 243-block nodes and `charted` on a 9-block's, tile format 7; `galaxy-map-wire-format.md`). Unit tests `tests/js/regionlabels.test.mjs`, `tests/test_region_node_lists.py`, `tests/test_map_region_labels_browser.py`. Served directly. |
| `../src/html/static/labelbudget.js` | The map's label budget (MAP.317): which of the things in view are offered a name. The budget is the screen area over the area one name takes (`labelBudget`, 120 by 24 px); only things the map draws are offered (a marker that is not drawn gets no name); five domains share it by weight (`SLICE_WEIGHTS`: inhabited 30, stars 30, uninhabited 25, nebulae 10, regions 5, a domain with less than its share lending the rest): inhabited worlds per type by population, then by habitability score (the inhabited and habitable world markers until population is stored, MAP.318; an item with `inhabited: true` goes before one that is not and its type is served first, MAP.321); the things that go first (`first: true`: the sectors, a region at the farthest zoom) are taken in the order they came, except that a sector goes by its total mass (`massKg`) heaviest first, so a short budget names the heaviest sectors and a sparse view names all (MAP.322); uninhabited bodies an equal share per category by mass (a belt or field by the volume it occupies); stars a third each by size, brightness and mass, once each; nebulae by volume; regions by share charted. `galaxylabels.js` turns the choice into the label engine's priority, so the first-ranked is placed first and none overlap. Unit tests `tests/js/labelbudget.test.mjs`. Served directly. |
| `../src/html/static/systemmap.js` | Toggles which System Map `<svg>` scene is visible (the whole-system view, or one per planet's own moon system) and fills the info side panel -- including "Atmosphere"/"Surface composition"/"Surface temperature"/"Life Chemistry" fields -- from a clicked marker's `data-*` attributes. Every visible scene's star/planet/moon markers also get their own live-rendered 3D sphere on `#sysmap-spheres-canvas` (the same vendored three.js build `sectorscene.js` uses): one shared WebGL context, redrawn each frame via a scissored sub-viewport per marker (never one `<canvas>`/context per body -- browsers cap concurrent WebGL contexts), each sized and positioned to exactly cover that marker's own `<circle>` and colored by its `data-color`, banded with a tilted ring for a gas giant (`data-bodytype`), and wrapped in a fresnel-glow atmosphere shell tinted by the atmosphere colour of the body's own air and star (the look record's `glow`, MAP.244; [how it is worked out](atmosphere-colour.md)) when the body has air -- the one genuinely 3D layer on this otherwise flat-SVG page, drawn behind the SVG so each marker's own stroke/label/life-badge still shows on top. An appearance layer only (no position data, and no cross-page link of its own to navigate) -- a marker whose sphere renders keeps its flat circle's fill transparent (`sysmap-sphere-active`) but never changes its actual plotted position. The diagram's body names are placed by the label engine's `placeBoxes` (`maplabels.js`, MAP.231: slid back inside the diagram, shifted or hidden when they would cover another name or a marker), and a "Body names" switch in the page's Menu (the `label-body` layer, address `hide=label-body`) shows or hides them. Served directly, same as `style.css`. |
| `../src/html/static/copycode.js` | The system page's Copy button: copies the generated Markdown out of the code box named by the button's `data-copy-target` (clipboard API, falling back to a selection copy on a plain-HTTP deployment). A separate file because the Content-Security-Policy allows no inline script; the box stays selectable by hand without it. Served directly, same as `style.css`. |
| `../src/html/static/bodyRendering.js` | Shared three.js building blocks (star granulation texture, fresnel glow shader, planet bands) used by `systemmap.js`, `systemview3d.js` and `sectorscene.js` to draw a body as a lit sphere (MAP.287 moved the gas giant's band texture and the painted-texture step here from `systemmap.js`). Served directly, same as `style.css`. |
| `../src/html/static/spin.js` | The slow turn of a rendered body on its axis (MAP.328), one function for every view that draws a body as a sphere: the rate comes from the body's stored rotation period (`spinH` of its look record) through one law for every body (one turn in 20 s at the fastest to one in 10 minutes at the slowest, the order of speeds kept), the axis from its stored spin axis (the 3D system view) or its tilt (a flat picture), a body with no stored period left still; the Spin switch (the System Map's Menu, kept in localStorage) and the reader's reduced-motion setting hold every body still; `tick(now)` is the spin time, which runs only while bodies turn; `turn(THREE, mesh, look, seconds)` sets a sphere's orientation. Used by `systemmap.js` (the spheres under the diagram) and `systemview3d.js`; no galactic view uses it (`tests/test_spin_scope.py`). Served directly, same as `style.css`. |
| `../src/html/static/galaxyprisms.js`, `galaxyblocks.js` | The Galaxy Map's sector grid, block ladder and density shading on shells, bands and wedges (`galaxyprisms.js`, the page's copy of `galaxy/spherical.py`, `galaxy/drill.py` and `galaxy/window.py`), and the block scene for one view packed as GPU-ready arrays (`galaxyblocks.js`, run in a Web Worker). Neither imports three.js, so both run under plain node for the tests. Served directly, same as `style.css`. |
| `../src/html/static/galaxysystem.js` | The Galaxy Map's system stage (MAP.125): `createSystemStage(host).open(href, name, center, mode)` fetches a system's scene (`/system/<id>/scene`), builds it with `systemview3d.js` at the star's place in the sector (compressed or true scale), adds its bodies to the map's picker with tooltip, ring and info panel, and keeps them moving with the orbit clock; `close()` removes it. `galaxystageview.js` decides when it is open (`&system=<id>`). Served directly. The bodies of the open system turn on their axes by `spin.js` (MAP.328): `update` advances spin time once a frame and `built.spin` turns each sphere about its stored axis, a body with no stored period (every star so far) left still; the Menu's Spin layer sets the standing preference (`spinSwitch`). The galactic views and the Sector Map's points of light turn nothing. |
| `../src/html/static/galaxystages.js`, `galaxystageview.js` | The Galaxy Map's drill-down from the galaxy to one sector (`design/interface/galaxy-drilldown-navigation.md`): the stage rules (`galaxystages.js`) and the scene, camera moves, breadcrumb and stage URLs that draw them (`galaxystageview.js`), backed by `/galaxy/stage`; a generated sector is opened in place as the last stage (`galaxysector.js`, MAP.66). The sector blocks are not picked by the pointer on the Galaxy Map (MAP.297; only a sector page's map and the choosing of a NAV end still take them): a click selects the star or object under it and shows its panel, a double click flies to it, and a click on nothing puts the panel away (MAP.304; the pick layer "stars" in `galaxymap3d.js`). MAP.336: every move to a named place is one `goTo(target)` (a sphere, a point and distance, an object, a body, a sector with its system and body, a region with its landing, a block, a stage, up one level, home); the crumbs, the card buttons ("Go to"), the names, the double click, the find-a-place box and the ruler's stops all call it, and a target with a `label` goes into the visit log. Served directly, same as `style.css`. |
| `../src/html/static/phenomenonmap.js`, `mapzoom.js` | Zoom and pan for the phenomenon page's AU-scale SVG diagram (`mapzoom.js` is the shared viewBox zoom/pan it uses). Served directly, same as `style.css`. |
| `../src/html/static/phenomenonrender.js` | Draws `planetgen/web/maps/phenomenonrender.py`'s "View" panel with three.js. Served directly, same as `style.css`. |
| `../src/html/static/generatefinder.js` | The "A named sector" tab's finder on the Add sectors page: name search through `/galaxy/locate` and the paged list from `/admin/generate/sectors`. Served directly, same as `style.css`. |
| `../src/html/static/generateplace.js` | The Add sectors form (ADM.108): shows the one panel of the choice made (sectors around a place or by shell, band and wedge; the place tab; whole shell or range) and hides the rest; without script every panel is shown and the radios still decide what is posted. Served directly. |
| `../src/html/static/generatejobs.js` | Keeps the Generate page's "Current job" panel live by polling `/admin/generate/status`. Served directly, same as `style.css`. |
| `../src/html/static/slowload.js` | The loading state of a slow page (UX.118): polls `GET /api/page-jobs/<id>` every second, draws the job's estimated progress with `numberformat.js` and `unitladder.js`, and when the job has succeeded fetches the page again and swaps its content in without a reload (or carries on with the next job the page waits for); Stop waiting, Check again and Try again (`X-Slow-Retry`). Served directly. |
| `../src/planetgen/web/maps/systemscene.py` | The 3D system view's data (MAP.69): `build_scene` turns `system_detail` into one JSON scene (every star, planet, moon, belt and comet with its ref, radius, colour, orbit elements, position at the epoch). Served by `GET /api/systems/<id>/scene` and relayed to the page at `/system/<id>/scene`. |
| `../src/html/static/systempage3d.js` | The system page's Diagram / 3D switch (MAP.74): loads the 3D view and the scene the first time 3D is chosen, remembers the last view, keeps `?view=3d&object=<ref>` in the address (MAP.67), wires the time control, scale select and the screen-reader body list, and says so when WebGL is missing. |
| `../src/html/static/systemview3d.js` | The 3D system view (MAP.72, MAP.73): the three.js scene (glowing stars, lit planets and moons, orbit lines and Kepler ellipses, belt particles, heliopause), the camera (drag to turn, pan, wheel and pinch zoom, W/A/S/D and arrow keys, fly to, follow, reset), picking, the info panel and decluttered labels. The names of the bodies are the label engine's (`maplabels.js`, MAP.231), one DOM name per body, none over another, never picking. From MAP.287 every star, planet and moon is coloured and painted from the `look` record the scene carries (`bodylook.py`, the same record the System Map's `data-look` holds): `paintStar` and `paintSolid` textures at the tier the body's size on screen calls for, painted one per frame, a gas giant's bands and a planet's air shell, with no class-letter colour. |
| `../src/html/static/systemscale.js` | True and compressed scale (MAP.71): where and how large each body is drawn, and the note saying which is shown. |
| `../src/html/static/orbitpositions.js`, `../src/html/static/orbitclock.js` | Positions of every body at any time and the time control (MAP.70); `orbitpositions.js` is the twin of `physics/body_positions.py`. |
| `../src/html/static/mapcore.js` | The helpers the Galaxy Map, the Sector Map and the System Map share (MAP.63): reading the panel's scene JSON, theme colors, info-panel fields, a sector's address, the highlight ring's texture, the scale bar's nice numbers and screen span, fitting the renderer to its canvas, and picking a point of light on screen. Served directly, same as `style.css`. |
| `../src/html/static/mapfit.js` | The size of a map's window on its page (MAP.291): `fitHeight` (pure: the window's height less where the viewport's top lies down the page, the toolbar under it and a gap, never under 320 px) and `watchWindowFit`, which sets a viewport's height from it now, on window resize and when the page above changes height. The viewport is as wide as its column (`style.css`), so its shape follows the browser window's, and `mapcore.js` `watchResize` makes the renderer follow it. `galaxymap3d.js` uses it for the Galaxy Map and sector page; UX.142 reuses it for the other viewers. Served directly, same as `style.css`. The Galaxy Map (not a sector page's map) fills its window instead (MAP.350): `fillHeight` (95% (`FILL_SHARE`) of the browser window's height less what a header held at the top of the window covers, never under 320 px; `fillBuffer` is the half of the rest above and below), `alignScroll` (where to scroll the page so the map has that margin above it under the header) and `watchWindowFit` with `{fill: true}`, which sets that height and scrolls the page to rest with the map there once it has loaded and as what is above it settles, until the reader scrolls the page themselves. Tests: `tests/js/mapfit.test.mjs`, `tests/test_web_browser_map_window.py`. |
| `../src/html/static/flytofit.js` | Where the camera sits for one target to fill the view (MAP.319): `fitDistance(radius, vfovDeg, aspect, {margin, minRadius})` (pure: the radius over the sine of half the smaller field of view, plus an 8 percent margin, a sphere smaller than `minRadius` taken as that size; the Galaxy Map frames nebulae and regions at least `MIN_VIEW_RADIUS_PC` = 2 pc), `reachFrom` and `reachOfPositions` (how far a set of spheres or an orbit line reaches). `galaxystageview.js` uses it for the fly-to of a star with a system and of a cloud, and takes a sphere as `goTo({center, radius, arrive})`; `systemview3d.js` uses it for the whole-system start distance and the fly-to of a planet with moons (`moonSystemReach`, `systemReach`). Tests: `tests/js/flytofit.test.mjs`, the fixture-map and system-page browser tests. |
| `../src/html/static/mappick.js` | The picking, hover and info-panel layer the Galaxy Map and the Sector Map share (MAP.65): a picker that holds each map's layers in priority order (points of light picked on screen, meshes raycast, or a map's own test) and gives back what is under the pointer, a solid layer hiding what is behind it; the hover tooltip (`.map-tooltip`, kept inside the map); the selection and hover rings, a size on screen for a point of light or round a body in the scene; and the info panel both maps fill the same way: a heading, rows, the NAV buttons ("Start Here" and "End Here", or the one for the end being picked, `navpick.js`), the ☆ Bookmark button (`bookmarks.js`), page links, buttons and Generate buttons. The panel's `view` spec is the View button to the item's own page (UX.143, MAP.254), first among the actions on every map and kind; and the panel is a popup inside the map's window at every width (`data-overlay`, MAP.256, MAP.288), never a column beside the map or a second box over it: the whole spec with a Close button, drawn in the corner farthest from where the reader last pressed and not drawn until something is picked (the System Map's panel moves between its diagram's and its 3D view's windows, and `reveal()` draws a map's own text, its measurement, the same way). No map builds a panel of its own (`tests/test_info_box.py`). The Sector Map gains a hover tooltip and ring and a ☆ for a system, phenomenon or generated neighbor; the Galaxy Map's black holes, neutron stars, quasars and clouds gain a tooltip, a hover ring, the NAV links and a ☆. The popup follows the map's window (MAP.295): `style.css` sizes it as a share of the window's width and no taller than it. MAP.347: its first view is the name, the kind (the spec's `kind`, else its "Type" row), the first three rows (`FIRST_FIELDS`) and the buttons, and everything else (the other rows, **More**, the `after` nodes, Generate, a hint) sits in `.map-info-rest` behind an **Expand** button (`fit()`, also `refit()`; the choice is kept from card to card), a region that scrolls on its own and gives the card the pointer (`data-expanded`), never the whole card. The card is set at a reading size (0.95 rem) on a darker see-through background and keeps clear of the zoom controls (MAP.325); a `more` list in the spec is drawn as rows under a collapsed **More**, for the facts a reader rarely needs; an expanded card takes the map's whole height beside the zoom controls (`data-beside`) in a window wide enough for it (MAP.332). The Galaxy Map's layer legend and scale line are one stack in its lower left, the legend above the scale (MAP.323). MAP.348: `createHoverPreview(host)` is the small card a hovered star, system, region, cloud or mark shows after the pointer rests on it: `firstViewOf(spec)` (also what the card draws first) gives its name, kind and rows, it shows the first three, and it takes the first corner of the map, nearest the top and farthest from the points to avoid, that is clear of them, the open card, the legend and scale stack, the zoom controls and the open Menu; it takes no pointer event and goes on a press, the wheel, a key or the pointer leaving. A layer or the sector stage offers its facts through `hover` (the Galaxy Map's `hoverWith`, `previewHover`; the sector stage's `host.preview`). Served directly, same as `style.css`. |
| `../src/html/static/navpick.js` | The NAV course pick the maps share (NAV.29, NAV.33): which end is being chosen and the end already chosen, the `?pick=` query every map URL carries, the info panel's "Start Here" and "End Here" buttons (a link to `/nav` with both ends once one is chosen, else a button that keeps the user on the map and asks for the other end), the info panel's "What's within 10 pc" link to the What's nearby page with the object, a sector or a cell's middle filled in (NAV.45; none while a course is being picked), and the banner's words. Served directly. |
| `../src/html/static/mapcontrol.js` | The camera and input controller the Galaxy Map's drill-down and the Sector Map share (MAP.64): turning and moving an orbit camera, zoom policies (free, a short range, or locked), the wheel, a two-finger pinch, the arrow keys, and telling a drag from a click. Served directly, same as `style.css`. |
| `../src/html/static/starlight.js` | The brightness boost for faint stars (MAP.87), the Galaxy Map's twin of `planetgen/web/maps/starmap.py`'s `star_light_boost`, which sizes the Sector Map's points: up to four times the halo light at the dim end, none from 1000 L☉ up. Served directly, same as `style.css`. |
| `../src/html/static/unitladder.js` | The maps' and panels' unit formatter, the browser mirror of `planetgen.util.format.format_quantity`: a distance in km, AU, mpc, ly, pc and up, a speed in km/h to c, a time in µs to Gy, each in the largest unit it is at least 1 of. It reads `unitladders.js`, generated from `planetgen/util/ladders.py` by `scripts/build_unit_ladders.py` (a test fails when the file is out of date). `formatMeasure(quantity)` shows a Quantity on the ladder of its kind, and the `format*` names (`formatDistanceKm`, `formatPeriodYears`, ...) are thin wrappers that read the registry unit they are for (UX.128); the file declares no unit factor. Served directly, same as `style.css`. |
| `../src/html/static/unitregistry.js` | The browser copy of the unit registry (UX.140), generated from `planetgen/physics/unit_registry.py` by `scripts/build_unit_registry.py` (a test fails when it is out of date): every unit's kind, symbol, name, factor and offset, every kind's base unit and bounds, and `convert`, `toBase`, `fromBase`, `scale` and `check` with the same arithmetic as Python (a test compares every conversion with Python through Node). Scripts read their unit factors from it and declare none. Served directly. |
| `../src/html/static/quantity.js` | The browser's `Quantity` (UX.128), the mirror of `planetgen.physics.quantity.Quantity`: a number of one kind in one registry unit, answered in every other unit of its kind with `as(unit)`, `base()`, `set`, bounds checks that throw as Python's raise; `held(unit, value)` builds one from a unit alone. Served directly. |
| `../src/html/static/numberformat.js` | The browser mirror of `planetgen.util.format.format_number`: a number shows as scientific notation ("1.23 × 10⁶") from 7 whole digits when it has no decimals, and from 5 whole digits when it has decimals (UX.36). Templates use the same rule through the `num()` global. |
| `../src/html/static/units.js`, `unitpreference.js` | The reader's unit system (UX.78): the server writes each measurement as `<span class="qty" data-q=... data-si=...>` with the Automatic text (`web/lib/fmt.py`'s `quantity_html` for a `Quantity`, `measured_html` for a number in a registry unit); `units.js` re-renders it from the SI value with the shared unit ladders and runs the settings menu's "Units" button, and `unitpreference.js` keeps the choice (Automatic, Metric only or Customary) in `localStorage`, guarded so the page works without it. |
| `../src/html/static/datatable.js`, `datatablestate.js` | The site's data table (UX.41): `datatable.js` takes over the table a list page rendered (`partials/datatable.html`, `web/lib/datatable.py`) with TanStack Table (sort and filter state) and TanStack Virtual (row scrolling), fetching further pages from `/table/<name>`; `datatablestate.js` is its logic without the page (the query a state turns into, which pages the visible rows need, the store of fetched pages) so plain node can test it. |
| `../src/html/static/breadcrumb.js`, `pagecrumbs.js` | The one breadcrumb every page and map uses (NAV.14): a single line at any width, with a "…" menu for the steps it hides. `pagecrumbs.js` folds a server-rendered page's breadcrumb list onto that line; without script the list wraps. |
| `../src/html/static/picker.js`, `objectref.js`, `systemnav.js` | The one selection every map and page shares (NAV.13: `select`, `up`, `into`), the browser mirror of `planetgen/galaxy/objectref.py`'s public reference form (`<kind>:<id>`, NAV.7; ids stay strings), and the NAV buttons, the "What's within 10 pc" link and the View button of a body's info panel on the System Map (NAV.50, NAV.45, MAP.254). |
| `../src/html/static/mapnavpad.js`, `gridopacity.js` | The Galaxy Map's on-map controls: a direction star and a readout of the center (MAP.200), and (MAP.336) the log ruler up the map's right side: the zoom slider, top the whole galaxy and bottom a moon, is the camera's distance on one log scale (0 at the whole galaxy to 1 at 1e-10 pc, the same wherever the view is, and a drag, the plus and minus buttons and a pinch take that distance as a goal the stage view carries across levels as the wheel does, MAP.333), with a named stop beside it for each scale (`maplogruler.js`), each a Go to for that scale about the middle of the map (a name too close under the one above keeps its tick, and a map under 40 rem wide keeps only the ticks), and the distance readout in the map's lower-left stack; a card on the right stands to the left of the ruler. Also `gridopacity.js`, the cap that keeps the drill-down's stacked translucent blocks from ever being more than 50% opaque along a line of sight (MAP.188). |
| `../src/html/static/nebulamesh.js`, `nebulaview.js` | A nebula drawn from its shape instead of a flat sprite (GEN.75, MAP.103): the triangle mesh `/galaxy/nebula/<id>/shape` serves, shared by the maps; and the nebula page's 3D view (MAP.105) with the brightest stars round it from `/galaxy/nebula/<id>/surroundings`. |
| `../src/html/static/tilepack.js` | The reader of the packed Galaxy Map nodes `/galaxy/nodes` sends and the node store keeps (MAP.159); `web/lib/tilepack.py` writes them and describes the layout, and `tests/test_tilepack.py` with `tests/js/tilepack.test.mjs` hold the two together. |
| `../src/html/static/menuplace.js`, `navtimes.js`, `connectionnote.js` | Button menus that open where they can be seen (UX.85: above the button when there is more room, height capped); the NAV result's Warp/Fold travel-times switch (UX.70); and when to tell an admin a live job stream was lost (ADM.40; the stream ends on purpose every `STREAM_SECONDS`). |
| `../src/html/static/joblog.js`, `generateranges.js`, `generatestarmix.js` | The Generate page's job output box (ADM.22: an Xterm.js terminal, falling back to a `<pre>`), the range sliders' read-outs (GEN.183) and the star mix's running total, which must reach 100% (ADM.45). |
| `../src/html/static/components.js`, `shoelaceicons.js`, `shoelace-theme.css` | The shared Shoelace component set (UX.40, vendored under `vendor/shoelace/`; `components.js` is loaded as a module on every page and registers the elements), the icon libraries pointed at files under `vendor/shoelace/assets/` (Shoelace's own data-URI icons would break the Content-Security-Policy), and the theme that maps Shoelace's tokens onto `style.css`'s. `icons.svg` is the site's icon sprite. Form fields (UX.124): the set includes `sl-input`, `sl-select`, `sl-option` and `sl-checkbox`, whose edge is the `--line` token (3:1 or more against the panel in each theme; `shoelace-theme.css` sets `--sl-input-border-color` from it, and the native fields that stay native use it too). |
| `../src/planetgen/web/templates/partials/ui_fields.html` | The one form field (UX.124; the Generate task pages, the one-off system page and the recipe name boxes use it, UX.125): macros `ui_input`, `ui_select` and `ui_checkbox` write a labelled `sl-input`, `sl-select` or `sl-checkbox` with its name, required flag and a `help-text` slot for the server's validation error (`error=`, shown as `.field-error` with `role="alert"`, and `data-invalid` on the element), and the same field as a plain control inside `<noscript>`, so a form still posts with script off (the element is hidden until defined). Used by the search, nearby and Plot Course chart forms and the data table's "Starting at" box and filter menus. Hidden, submit, password and one-time-code fields stay native, and so does a field with a `datalist` (a datalist cannot reach into the shadow root). |
| `../src/html/static/fieldsready.js` | `fieldsReady()` (UX.125): a promise that resolves when `sl-input`, `sl-select` and `sl-checkbox` are defined (components.js loads them after the page's own scripts first run) or after 3 s, so a script that reads or listens to the form fields (`generatelevers.js`, `generatestarmix.js`, `generatesize.js`) waits for the fields to have a `value` and events. Those scripts also listen for `sl-change`, which a Shoelace field sends where a plain one sends `change`. Served directly, same as `style.css`. |
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
meant to expose. The Document's HTML renderer
(`describe/html.py`) escapes every string it writes and re-enables nothing, so a
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
[`deployment/README.md`](../deployment/README.md#what-every-setup-has-in-common);
each platform guide in [`deployment/`](../deployment/README.md) sets the same
rule).
The ES modules that import siblings (`galaxymap3d.js`, `systemmap.js`,
`galaxymap3d.js`, `galaxystageview.js` -> `mapcore.js`, `mapfit.js`, `mapcontrol.js`, `mappick.js`,
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
| `/systems` | `browse.py#standalone-systems` | Every system with its sector and octant as a data table; the standalone ones (made outside any sector) are a filter on the list (UX.55). `?uncharted=1` (UX.87; a badge on the page counts them) swaps it for every star the brightness scatter placed that has no system yet, brightest first: luminosity, star, temperature, age, sector address as shell, band and wedge (linking to the Galaxy Map there), galaxy and in-sector coordinates; a Members column ("3 members, 2.4 solar masses" for a star with company, a dash for a single star) and each star's name links to its page, `/systems/uncharted/<id>` (UX.114: the star's data, its members and the orbits that bind them, and a Compact binaries panel when it holds an X-ray binary, a microquasar or a spider pulsar, UX.186; a built star sends the visitor to its system); rows from `/table/uncharted_systems`. An admin gets a Generate button per row (`POST /systems/uncharted/<id>/generate`) that builds just that system and opens it; a note recommends generating the whole sector. |
| `/galaxy` | `galaxy.py` | The 3D Galaxy Map plus the Quadrant summary; `?quadrant=I\|II\|III\|IV` lists that Quadrant's placed sectors in a data table (UX.41), nearest the core first; sort, `zone` filter and `page` in the address, rows from `/table/galaxy-quadrant`. The info panel's "View sector" is a plain link. `?pick=from&to=...` or `?pick=to&from=...` (NAV's "Pick on Galaxy Map") shows a banner with Cancel back to `/nav`, keeps "Charted only" on, and keeps the pick on the sector it opens in place (`?sector=<designation>&open=1`) and on a sector's page (`/sector/<id>?pick=...`). |
| `/galaxy/stage?at=...` | `galaxy_views.py` | JSON for the map's drill-down: one stage's generated counts, `at` a block key `m.shell.band.wedge` or none for the galaxy (see `/api/galaxy/stage`). Through `planetgen/web/lib/tilecache.py`'s disk cache (`fetch_stage`), plus `nebula_cover` (MAP.174: `[shell, band, wedge, cover]` for every child cell with a share above zero inside a nebula, filled or not, and `stats.nebula_cover` on each generated child; read from `/api/galaxy/nebula-cover` on every request, never cached, `null` if that read fails; empty for the cells of a sector); 400 on a malformed key, 502 on an API failure, both `{"error": ...}`. Shares the `galaxy_tiles` rate limit. `Cache-Control: no-store`. |
| `/galaxy/window?at=...&g=...` | `galaxy_views.py` | JSON for the map's data window (MAP.151): `g` cells across (3, 5, 7, 9 or 11; 9 by default) around the cell `at`, a block key `m.shell.band.wedge` or a sector `1.shell.band.wedge` (see `/api/galaxy/window`). Built from the stages of the parent blocks it touches (`planetgen/web/lib/tilecache.py`'s `fetch_window`), plus `nebula_cover` as on `/galaxy/stage`; 400 on a malformed key or width, 502 on an API failure, both `{"error": ...}`. Shares the `galaxy_tiles` rate limit. `Cache-Control: no-store`. |
| `/galaxy/near?x=...&y=...&z=...&r=...` | `galaxy_views.py` | JSON for the map's camera coming near sectors (MAP.335): the generated sectors within `r` parsecs (at most 4.5 sector edges) of the point `x, y, z` (the eye, galaxy frame, parsecs), as `{at, edge_pc, r, sectors, stamp}`: `at` the key of the sector holding the point (`null` outside the grid) and `sectors` the generated ones of the cells round it as a stage lists them (`shell`, `band`, `wedge`, `id`, `name`, `system_count`, `stats`), so a sector's scene can be asked for by its id without a stage window. Cut from the cached stages like a window of sectors (`planetgen/web/lib/tilecache.py`'s `fetch_near`); 400 on a number that is missing or not finite, 502 on an API failure, both `{"error": ...}`. Shares the `galaxy_tiles` rate limit. `Cache-Control: no-store`. |
| `/galaxy/locate?q=...` | `galaxy_views.py` | JSON for the map's find-a-place box: sectors and star systems named like `q`, each with its sector address (see `/api/galaxy/locate`). A blank `q` asks the API nothing; 502 on an API failure, as `{"error": ...}`. Shares the `search` rate limit. |
| `/galaxy/territories` | `galaxy_views.py` | JSON for the map's Territories overlay: `/api/territories`' owned systems and capitals, with each polity's name, color and system count folded in from `/api/polities`. 502 on an API failure, as `{"error": ...}`. Shares the `galaxy_tiles` rate limit. `Cache-Control: no-store`. The map shows its Territories button only once `/api/polities` counts at least one polity. |
| `/galaxy/uncharted/<shell>/<band>/<wedge>/scene` | (new) | JSON scene of a sector cell nothing was generated in (MAP.162): the waiting bright stars and unbuilt scattered phenomena the scatters left there, for the map to open in place with an "Uncharted" mark. A generated cell redirects to its sector's scene; a cell outside the grid is a 404 and an API failure a 502. Shares the `galaxy_tiles` rate limit. |
| `/galaxy/system/<id>/nearest` | (new) | JSON for the Galaxy Map's star system card (MAP.296): `{"nearest": [{id, name, distance_ly, url}]}`, the same stored nearest neighbours the system page's "Nearest:" line names (`nearest_neighbors` of `GET /api/systems/<id>`), nearest first, each with its page's URL. An unknown system is a 404 and an API failure a 502, both `{"error": ...}`; cached privately for five minutes. |
| `/galaxy/nebula/<id>/shape`, `/galaxy/nebula/<id>/surroundings` | (new) | JSON for the maps' nebula meshes (MAP.103, `?lod=low\|full`; `GET /api/nebulae/<id>/shape`) and for a nebula page's 3D view (MAP.105: the brightest stars round it). 404 for an unknown nebula, 502 on an API failure. Shares the `galaxy_tiles` rate limit. |
| `/galaxy/nodes?keys=...` | `galaxy_views.py` | JSON-free packed data for the map's script (MAP.260): the requested sector-block nodes (`keys=m.shell.band.wedge,...`, a sector's `1.shell.band.wedge`; at most 128 a request) and, given the browser cache's `stamp`, the changed nodes since (`history` entries carry `nodes`). The map also sends the node keys it holds nothing of only, and a request marked `X-Map-Prefetch` is a prefetch (MAP.158). Through `planetgen/web/lib/tilecache.py`'s disk cache (`fetch_nodes_wire`, files `n<m>_<shell>_<band>_<wedge>.bin`), packed by `tilepack` (the body is those bytes joined, not parsed and serialised again; Apache compresses it, Brotli when `mod_brotli` is on); builds the blocks around the request ahead (`web/lib/tileprefetch.py`). 400 on a malformed request, 502 on an API failure, both `{"error": ...}`. `Cache-Control: no-store`. |
| `/admin/generate/task/<key>/size` | (new) | Admins only (403 otherwise): JSON for the line under a task's form (`static/generatesize.js`, ADM.126). The form's fields are the query string; the answer is `state` `ready` (`line`, `refused`, `refusal` and `estimate`: the figures of the Review step, from `planetgen galaxy --estimate-only` and its cache), `none` (the job has no estimate) or `invalid` (`message`: why the values cannot be sized yet). Starts and writes nothing. A typed palette command names its own address in `size_url`. |
| `/admin/generate/shape-report` | (new) | Admins only: JSON for the line under the Generate page's shape dials (`static/generatelevers.js`): `?mass=...&flatness=...` (the Plan fields) answer with `line` (the disc's crest over the one-star density, from the built shape and the Census dial, the share of the galaxy's systems in sectors that expect one, integrated from its density, the outline's shells and candidate sectors), `warnings`, `refusal` (set when the outline would be empty) and `sample` (the Auto dials were drawn from a sample seed because no galaxy is planned); 400 `{"error"}` for a value out of range. The same numbers `planetgen plan` prints. |
| `/admin/generate/sectors` | (new) | Admins only: JSON, 50 filled sectors a page (`?page=N`), for the named-sector tab's list (`static/generatefinder.js`). |
| `/admin/generate/status` | (new) | Admins only: JSON for the Generate page's live job panel (`static/generatejobs.js`). |
| `/admin/generate/jobs/<id>` | (new) | Admins only: one past job with its full output. |
| `/admin/generate/jobs/<id>/stream`, `/admin/generate/jobs/<id>/log` | (new) | Admins only: a job's output as Server-Sent Events (`log` events with new output and the byte offset reached as the event id, so a reconnect resumes; `state` events; `done` once all output is sent; the response ends after `STREAM_SECONDS` and the browser reconnects), and its full output as a plain-text download (`/log`; `/log?form=records` is the same job as records, one log line per line with its time in UTC, severity and component). The job log's lines, stage dividers and the records file are described in [`logging.md`](logging.md). |
| `/admin/generate/system/download` | (new) | Admins only, POST: returns the one-off system page's text as a `.md` or `.txt` file. |
| `/system/<id>` | `system.py` | One star system: badges, the "At a glance" panel under the map (UX.190, `web/glance.py`: stars, planets, moons, belts, comets, habitable worlds, the best world with its symbols, facilities, system perimeter; a row that would say none or 0 is left out except planets and habitable worlds), "Navigate from/to here" (systems in a sector), nearest-neighbour location links, the System Map, the system's overview text, the "Stars and bodies" summary tree (UX.173: the system's stars, each already open onto its planets, belts and comets; every row a link to the object's own page and every deeper level loaded on open from `/object/<ref>/tree`), `?code=markdown` code view with a Copy button, the Facilities panel (with each facility's type, affiliation and rating chips, worst first), the "Species that could live here" panel (POP.12: for each rocky world, the species that could live there unaided and those that could with equipment, each with its six domain chips; hidden when there are none), and for an admin the facility forms: place one (with a type and its fields, an affiliation and ratings) and edit a facility's affiliation and ratings or have the unset ratings suggested (see below). |
| `/lists`, `/lists/<kind>` | `object_lists.py` | Every list of objects on the site, and one list per kind (stars, planets, moons, asteroid belts, comets and each phenomenon type) as a data table (UX.41) that reads the object index through `GET /api/object-lists/<kind>`: sort by name, class, mass or radius, jump to a name with the "Starting at" box, and scroll or page by key, so a deep page costs what the first does (UX.141). The stars list shows each star's temperature, luminosity and age, and the planets and moons lists their orbit size, eccentricity and surface temperature (UX.187; stored values, not sortable). The planets and moons lists have a Habitability column of the PHI-6 symbols (UX.100). `/lists` also links `/systems`, `/sectors` and `/phenomena`. |
| `/preferences` | `preferences_page.py` | The Preferences page (UX.236): each per-browser setting of `static/mappreferences.js` listed with its plain-words description and its default, a control that stores a change as it is made, and a Reset to defaults button. The list is drawn in the browser by `static/preferencespage.js`; nothing is kept on the server (the same keys move into the account when accounts land, USR.7). Linked from the gear menu and the Account page; private to search engines. |
| `/phenomena` | `phenomena.py` | Every exotic phenomenon, paged with `?page=N`. An unbuilt scattered row with company shows its members and total mass in a Members column, and its page (`/uncharted-phenomenon/<id>`) lists the members and the orbits that bind them (UX.114). |
| `/object/<ref>` | `item_pages.py` | NAV.8, UX.143, UX.170: any object reference (`planet:12`, `moon:7`, `belt:2`, `comet:9`, `star:3`, `system:9`, `sector:4`, `nebula:3`; a bare number is a system). A sector, a system and each phenomenon type go (302) to the page the site already has for them. A star, planet, moon, belt, comet, facility, globular cluster or superbubble (UX.183: a facility has no example, tree or Navigate; its "About" names the body it is built on, its type, affiliation, description and rating chips. UX.184: a cluster or superbubble has none either; its "About" says what it is, where its centre lies, a cluster's mass, age and metallicity, and lists the star systems, nebulae and other objects inside it, each linked) gets a page of its own, in this order: the Location block (UX.144, `partials/location.html`: sector, system, the star it orbits when its name differs from the system's and, for a moon, its planet; UX.200: no Galaxy link, no line for a cloud it is not in, the block left out when empty), links to the previous and next object at the same level, a rendered example for a star, planet or moon (UX.171, `static/bodyview.js`), the "At a glance" panel (UX.188, `describe/sheet.py` `GLANCE`: what a captain wants first; a planet's or moon's PHI-6 symbols and life scores are its first row on a line of their own, then class, temperature, pressure, gravity, radiation, hours outside, day length, distance and year), the whole property sheet (UX.168; it follows "About", Boss 10-24), every stored value grouped; UX.188: a row that holds none or 0 is left out unless its absence matters, `SHOW_WHEN_ABSENT`, and a line says how many are not listed; UX.203: a planet or moon with no life, under `VACUUM_LIMIT_PA` of surface pressure or effectively no water has that section left out, the atmosphere shown as "None" and one sentence for each in "About", `describe.sheet.missing`; UX.204: the property groups start collapsed and a link to one opens it, `static/opensection.js`), "About" (the description's paragraph, the dominant species of a life world, the facilities the body hosts, each linked to its own page) and the summary tree of what belongs to it (UX.172). A planet's or moon's page also shows its PHI-6 symbols and its habitability chips, and an admin gets that body's own Admin menu (Regenerate, Change class, Delete; a POST here). `item_pages.page_kinds()` is the registry of every kind a reference can name; `tests/test_item_pages.py` holds it equal to `objectref.KINDS`. A reference that names nothing is a 404. A nebula's and a remnant's page lists what is inside it (`partials/contents.html`). A physics row of the property sheet (eccentricity, tide, interior heat, spin, field, rings, the life timeline) carries a one-line source under its value (`describe/sheet_sources.py`, UX.187). The design is [object-pages.md](../design/interface/object-pages.md). |
| `/object/<ref>/tree` | `item_pages.py` | UX.172: one page of an object's children as the rows of a tree level, `{"html", "next"}` (`?after=` is the previous page's token; `Cache-Control: no-store`), for `static/objecttree.js`. Not a page (a 404 for a reference that names nothing, a 400 for a bad `after`). |
| `/phenomenon/<type>/<id>` | `phenomenon.py` | One phenomenon's data table and its view: a three.js render for a neutron star, black hole, quasar, rogue planet or comet (`planetgen/web/maps/phenomenonrender.py`, `static/phenomenonrender.js`, an SVG still without JavaScript), the AU-scale diagram for a nebula or remnant, and none for an asteroid field, with "Navigate from/to here", and an "All stored properties" panel of the same property sheet a body's page shows (UX.175). `<type>` is one of `nebula`, `asteroid_field`, `black_hole`, `neutron_star`, `supernova_remnant`, `rogue_planet`, `interstellar_comet`, `quasar`; anything else is a 404. |
| `/classes` | (new) | The class reference: every class type (star spectral and luminosity classes, planets, nebulae, supernova remnants, asteroid fields, black holes, rogue planets, comets) with how many classes it has. |
| `/classes/<type>` | (new) | One type's classes, each linking to its page, plus notes (an asteroid field's size digit, for one). `<type>` is `star-spectral`, `star-luminosity`, `planet`, `nebula`, `supernova-remnant`, `asteroid-field`, `black-hole`, `rogue-planet` or `comet`; anything else is a 404. |
| `/classes/<type>/<code>` | (new) | One class's facts, e.g. `/classes/planet/M`, `/classes/star-luminosity/IA+`, `/classes/comet/halley_type`. An unknown code is a 404. The system page links a star's type and a planet's or comet's class here, and the phenomenon page its Class (an asteroid field's `C3` by its letter) and a rogue planet's Mass Class. |
| `/help/habitability` | (new) | In the Help menu (UX.180): what the equipment levels and PHI-6 colours on the planet pages mean (UX.90): the survival equipment items and when each is needed, the PHI-6 symbols in each rating colour (UX.101), the six domains with their Blue, Green, Yellow and Red thresholds, and the chips a planet row shows. |
| `/help/properties` | (new) | In the Help menu (UX.187): every value a body's page shows, by kind (planet, moon, star), with its unit and where it comes from (the step that fixes it: the formation, the orbit, the star, the climate). The item page links to it from "Where these values come from". |
| `/help/life` | (new) | In the Help menu (DOC.84): the four kinds of life (microbes, flora, fauna, sentient life), what each word means and what its 0 to 9 score is judged from, with links to the pages below. |
| `/help/life/<topic>` | (new) | DOC.84: one page of the life help. `abundance` (the words absent to dominant and their scores), `trend` (growing, stable, declining and how they are judged), `security` (secure, marginal, precarious, dying out and the time each means), `suitability` (how a life class is scored 0 to 9, and best-suited against living there), `stages` (the stages of life, the three animal grades and the pace of evolution) or `classes` (the nine life classes in one table, with the speculative chemistry). Any other topic is a `404`. Every number is read from the code that uses it. |
| `/help/life/classes/<number>` | (new) | DOC.84: one life class family (1 to 9): each variant with what it drinks, eats, breathes and gives off, the liquid and energy it needs, its temperature, pressure and supply needs, how far it can go, how sure we are it can exist, the words of its stages and whether its air chemistry is speculative. A planet's page links its class names here (`#class-1a` is the variant). |
| `/species` | (new) | Every species, paged with `?species_page=N`; `?spacefaring=1` or `0` filters. Until a population pass has made species (`GET /api/population`) the page answers 200 with a line that says the pass has not run, and the header's Species section is hidden; a species page for an id that does not exist is a 404. |
| `/species/<id>` | (new) | One species: homeworld, body plan, era, civilization age and polity, and its environment limits (POP.12): each limit of the habitability table with the species' value, the human value and whether they differ. |
| `/polities` | (new) | Every polity (`?polities_page=N`): species, government, capital, systems held, reach. Until there are polities the page answers 200 with a line that says so (the pass has not run, or no species founded a polity). |
| `/polities/<id>` | (new) | One polity and the systems it holds, nearest its capital first (`?systems_page=N`). |
| `/sector/<id>` | `sector.py` | One sector: badges, the 3D Sector Map, an "At a glance" panel (UX.191, `web/glance.py`: systems and stars, each kind of phenomenon, the compact-object hazards and the nearest, facilities, a superbubble around it; rows that would say none are left out), its Profile panel (UX.113: what the galaxy's shape says the cell holds, from `GET /api/galaxy/profile`: region, expected systems, the stellar populations with age, [Fe/H] and velocity dispersion, gas, superbubble and dark matter), and its Contents table (systems, nearby phenomena and the facilities outside its systems, nearest the center first, as a data table (UX.41); the arrow before a system's name opens its stars and bodies under it as the summary tree (UX.174, `static/datatable.js` over `static/objecttree.js`), and an opened row stays opened as the table scrolls: `contents_sort`, `contents_order`, `contents_type`, `contents_octant`, `contents_page`; rows from `/table/sector-contents`); for an admin, an Admin menu (generate neighborhood) and an Edit menu, each pick opening its form in a dialog. |
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
| `/account` | `changecreds.py` | Change the admin username and password, and turn two-factor sign-in on or off (QR code, then recovery codes shown once; forms `POST` to `/account/two-factor`), and set or remove the account's picture (PNG, JPEG, ICO or GIF, animated too; `POST` to `/account/avatar`; the caps are under `POST /api/auth/avatar` in docs/reference/api.md). |
| `/admin` | `admin.py` | maintenance mode, a "Check the math" button (POST `action=check_math`: runs the math check and shows its summary, PERF.170), API keys (list, create, revoke; `?keys_page=N`). |
| `/admin/stats` | `adminstats.py` | Server health (including Redis: answering, append-only file, memory policy, memory used, last save and live workers) and database stats (including about how many bright stars the plan pre-placed, from the `bright_stars` table's row estimate), every name made unique (`?names_page=N`), the galaxy's naming key with Change and Draw buttons (POST `/admin/stats/naming-key`, GEN.70), current login lockouts with Lift buttons (POST `/admin/stats/lockouts`), the newest failed sign-ins, and this server's generation speed per density bucket with the galaxy's size per star system (`GET /api/admin/generation-stats`). |
| `/admin/stats/generation-reset` | (new) | Admins only: POST, the Stats page's "Reset stats" button (deletes every recorded generation rate, PERF.32). |
| `/admin/queue` | (new) | Admins only: the work queue (ADM.10). Workers active and the server's load as "x / x / x" (1, 5 and 15 minutes), which run holds the machine and how long since it answered (a stale or ownerless one has a Clear button; Redis down is said in words), Pause the queue / Resume the queue, the bulk worker count (Set, or Ludicrous Speed: 95% of the cores; ADM.15), Clear finished runs, a Start a fill panel (ADM.129) with Regenerate habitability (`planetgen rescore-habitability`: re-scores every planet and moon from its stored conditions, re-draws nothing) and Fill map luminosity (`planetgen fill-map-luminosity`), each behind a confirmation page and started as a Generate page job with a generation run's rules (refused while another run is live, the one-run lock, the activity log entry `job.start`, Stop everything, Retry), then every run in Redis newest first (paged) with its title, state, who started it, start, duration, stage, progress and time left (PERF.99). `/admin/queue/<run id>` shows one run: its facts, its stages with their units and times, the tail of its log and a link to download the log (the last 20 runs keep theirs; `/admin/queue/<run id>/log`), with Cancel (a live run: the waiting units come off the queue and the started ones are stopped at once), Clear (a live run that shows no sign of life), Retry (a failed or cancelled run: starts a Generate page job) and Delete (a finished run). Every control first shows `/admin/queue/confirm/<action>?run=<id>`, then posts to `/admin/queue/action`. The Generate page's job page lists the runs the job started, each linked to its run page and log |
| `/admin/settings` | (new) | Admins only: every option of `config.json` (ADM.84), grouped by category, each with its help text, default, limits, unit, whether it needs a restart, where the value in use comes from (default, `config.json` or an environment variable) and a box fitted to its type (check box, number, text, one-per-line list, choice, write-only secret box). An option the admin may not change (command line only, or the Owner's) is listed read-only with the reason. "Use the default" removes a key from the file; "Clear it" empties a secret. Review the changes (POST) shows the old and new value of each changed option, which need a restart and which an environment variable overrides, or sends the form back with the message at each bad option and nothing saved; Save these changes (POST) then writes the file (a backup archive first, an audit row per option) with the hash of the file the form was drawn from, so a file changed by someone else in between is refused. The page works without script. While a saved option still needs a restart of the web server, every admin page shows a banner naming it. |
| `/admin/settings/backups`, `/admin/settings/backups/<archive>`, `/admin/settings/backups/<archive>/restore` | (new) | Admins only: the backup archives screen (ADM.85). The list shows each archive's time, who made it, why, which parts it holds and its size, with Back up now (POST) and, for the Owner only (from 8.2), Download. An archive's page says whether it passes its hash check, compares its settings part with today's `config.json` (old and new value of each option a restore would change, secrets as `(changed)`; a backup that no longer fits the model is explained and refused), and restores the ticked parts (POST): a safety archive of what is replaced first, an audit row per changed option and per restore, only the options an admin may change from the web. A settings file changed by someone else since the page was drawn is refused. |
| `/admin/settings/images` | (new) | Admins only (ADM.44): the site icon (SVG), the PNG or ICO icon, the home-screen icon and the link-preview picture in use, each with an upload form (checked when uploaded; the reason is shown when a picture is refused) and "Use the built-in one" (POST). The pictures are stored by content name and served at `/site-assets/<name>`. |
| `/robots.txt`, `/sitemap.xml`, `/sitemap-pages.xml`, `/sitemap-systems-<n>.xml`, `/sitemap-sectors-<n>.xml`, `/.well-known/security.txt`, `/humans.txt`, `/site.webmanifest`, `/site-assets/<name>`, `/og/<system|sector>/<id>.png` | (new) | Public (ADM.44): the files crawlers, browsers and link previews ask for, generated from the `seo`, `icons`, `og` and `site` settings at request time (docs/reference/config.md, "Web, link previews and search engines"). The sitemaps and the manifest are 404 when switched off or while `base_url` is the localhost default. |
| `/admin/settings/export`, `/admin/settings/import`, `/admin/settings/import/save` | (new) | Admins only: export downloads the options set in `config.json` that an admin may change as `planetgen-settings.json`, without secrets or command-line-only options; import takes such a file (up to 1 MiB), checks it like an edit, shows the review diff, and `import/save` saves it as `config.import` with a backup first. A bad value, an unknown key, an option an admin may not set or a file from a newer release changes nothing and is named. |
| `/admin/facility-types` | (new) | Admins only (POP.10): every facility type with its fields, and forms to add, change, take out of use and delete types and add, change and delete their fields, or to add the example type; each form posts back to the page, the API holds the rules and a refusal shows as the page's error. |
| `/admin/generate` | (new) | Admins only: generate, plan or reset the galaxy from the browser (see below). |
| `/admin/generate/system` | (new) | Admins only: one star system with every `planetgen system` option, shown as the system's description with Markdown and text downloads, and never saved (see below). |

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
`docs/reference/api.md`, "Deleting and regenerating"), flashes the outcome, the
bodies the re-validation moved and any warnings, and redirects (303):
back to the page, or after a delete to the system's sector (or the
Systems list), the Sectors list or the Phenomena list, and after a
sector regenerate to the new sector. Every page shows those flashed
lines under its heading (`base.html`), read only when the request
carries a Flask session cookie.

**Change star (ADM.7).** In the same menu on
the system page, a
single-star system's own row has Change star, a spectral-type field
(`edit_action=star`, `star_type`) with the same "also delete
facilities" checkbox, since bodies the new star can't hold are removed.
Both flash the outcome as above, including what moved and what was
removed (see `docs/reference/api.md`, "Changing a star").

**System page facilities.** The system page lists the system's
facilities (starbases, colonies, outposts; `GET /api/systems/<id>/
facilities`) in a Facilities panel (name, kind, host, placement, and an
orbital one's distance, period and speed as stored), and each one again
on its host's own page. The System Map draws each as a small
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

The Galaxy Map's disk node cache (`planetgen/web/lib/tilecache.py`) is written by the
WSGI daemon. A node cache location or size set only with `SetEnv
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
and an `action` field: `generate_neighborhood` (an admin with current
credentials, on a galaxy-placed sector). A successful action flashes its message and
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
| `planet_min_gravity`, `planet_max_gravity`, `planet_gravity_unit` (UX.227; the same for each key of `planetgen.db.search_values` with `star_`, `planet_` or `moon_` in front: `surface_temperature`, `surface_pressure`, `surface_dose`, `orbit_distance`, `mass`, `moon_count`, `luminosity`, `gas_<formula>`, ...) | A range on a stored physical value, under "Search by physical values", in the unit beside it (the unit menu offers the common units of its kind; any unit of the kind the API takes is accepted in the URL). A bound that isn't a number, an inverted range and a unit of another kind are ignored. A unit that is the key's own, or has no bound beside it, is left out of the address. Each range shows as a chip with its unit, and the planet, moon and star tables carry a Values column with the values the search asked for. |
| `planet_other_gas`, `planet_other_gas_min`, `planet_other_gas_max`, `planet_other_gas_unit` (and `moon_`) | The form's row for a gas outside its short list: it turns into the key of that gas (`planet_min_gas_C6H6`) and the address is redirected to that form. |
| `planet_value_sort`, `moon_value_sort`, `star_value_sort` | Orders that kind of object by a physical value: `gravity` smallest first, `-gravity` largest first; a body with no value comes last. A table header's own sort wins. An order alone is not a search. |
| `type`, `spectral`, `luminosity`, `class`, `body`, `life`, `equipment`, `moon_class`, `moon_body`, `moon_life`, `moon_equipment`, `density`, `phenomenon`, `phenomenon_class` | Tag facets, repeated once per selected value (`spectral=G&spectral=K`). `equipment` and `moon_equipment` are the survival equipment items a human needs (one tick per item, linked to the habitability levels page); `phenomenon` and `phenomenon_class` pick exotic phenomena by type and class. |
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
  (Galaxy, Sectors, Systems, Phenomena, Nav, Nearby, and Species once
  a population pass has made species; Classes and Habitability are under the Help menu beside the search button, UX.180) with `aria-current="page"`
  on the current one, a search box, and Login or Admin/Stats/Generate/Logout
  plus the theme button. Below 56rem (92rem for a logged-in admin, whose
  header has more links) the sections, search and account links fold into
  a native `<details>` "Menu" (works without JavaScript). A
  "Skip to content" link comes first, and every page but the home page
  has a breadcrumb trail.
- **Fast**: no Python process start per page, and no HTTP call from the
  page back to the API (see below). Visitors without a session cookie
  cost no login lookup at all.

**The Generate screens** (`/admin/generate`, `web/generate_flow.py` over the
logic in `web/generate_page.py`; ADM.107, ADM.108) run the galaxy tools an
admin used to run in a terminal on the server, as background jobs. The Start
screen (`templates/generate.html`) has a card per task: Add sectors, Start a
new galaxy, Change what is scattered, Bring the galaxy to the running generator
epoch (GEN.764: the task page reads the galaxy and says what `planetgen
bring-to-epoch` would redo before Review, and refuses when nothing can be done),
Make one system, Check the database and (last, apart) Wipe the galaxy. Each task is a page of its own
(`/admin/generate/task/<key>`, `templates/generate_task.html`) that posts to
itself and answers with the Review step (`templates/generate_review.html`):
what will run, each step's command line, the size and time of a galaxy run
(and the disk refusal with "Generate anyway"), then Run. Run posts the same
fields to `POST /admin/generate`, which builds the job and starts it, or
queues it behind the running one (ADM.123); Back returns to the form filled
in. The frame (`templates/generate_base.html`) shows the status strip, the
Start and History tabs and the current job on every screen, and the forms stay
usable while a job runs. The jobs below are what the tasks run:

| Form | Runs |
|---|---|
| New galaxy | `python3 -m planetgen.cli.reset --yes`, then `planetgen plan`, then `planetgen sectors` around a random start, then `link`, the scatters (`scatter --mass`, `--luminosity`, `--phenomena`), `backfill` and `settle`, each a step of the job. |
| Generate sectors | `planetgen galaxy` in any of its modes: around a random start; around a sector (`--center-sector --radius-pc` for a filled sector found by name through `/galaxy/locate`, picked from a paged list of filled sectors (`/admin/generate/sectors`) or typed by ID; or `--shell --band --wedge --radius-pc` for a sector address or a galaxy-frame x, y, z in pc, turned into the address of the cell holding it with the plan's sector edge); one address (`--shell --band --wedge`, with an optional neighborhood radius); a span of shells, bands and wedges (`--shells`, `--bands`, `--wedges`, each a number or `FIRST:LAST`, blank for all; ADM.29), with `--limit`, or `--yes` for a very large one; a whole shell or one band of it (`--shell`, `--shell --band`, with `--limit` or `--yes`); or the sectors a course crosses (`--course`, NAV.48, from the NAV chart page). The single-sector neighborhood radius can be given in light-years (13 ly up to the parsec limit). Sent with `Accept: application/json`, the form answers `202 {"job", "url", "status_url"}` (or `{"error"}`) so the Galaxy Map can start a job without leaving the map; behind a running job it answers `202 {"queued": true, "ticket", "position", "waiting_for", "status_url"}` (the job waits in the one waiting list, PERF.154) and `409 {"error"}` naming the list when it is full. The Sector Map's Generate buttons on an unfilled neighbor post straight to this form. Before the job starts the page shows its size, time and the free space on the drive holding the database's data directory (named, or "unknown" if it can't be measured) (`planetgen galaxy --estimate-only`) with a Generate button to confirm, or the refusal when the disk can't hold it; a JSON caller gets `409 {"error", "estimate", "confirm_field"}` and re-sends with `estimate_ok=1`. |
| Plan the galaxy | `planetgen plan` with the galaxy shape fields, then each scatter as a step of its own (`planetgen scatter --mass`, `--luminosity` and `--phenomena`, PERF.121), so the job shows each scatter's progress bar and the count it placed. On New galaxy and Plan, "Skip the bright-star scatter" leaves that step out. |
| Redo scatters | `planetgen scatter` (GEN.196, PERF.121): one box where the admin ticks the scatters to redo (the massive stars, the bright stars by luminosity floor, the phenomena) and sets each one's new limit (the stellar mass limit, the luminosity floor, the neutron star and black hole limit with its "Use the star mass setting" box). They run as one job: `planetgen clear-scatter` for each ticked scatter first (the bright stars before the massive stars), then `planetgen scatter` for each in the order the scatters have to run; an unticked scatter is left as it was. The central black hole or quasar is always rebuilt with the phenomena. The bright stars are split from the massive stars at the mass limit, so redoing the massive stars redoes the bright stars with them. |
| Add a dimmer layer of bright stars | `planetgen scatter --band --down-to N`: keeps the bright stars already placed and adds only those from N up to the current star-fill level (shown on the panel), leaving filled sectors out. Disabled until a scatter has run. |
| Reset | `python3 -m planetgen.cli.reset --yes`. |

Rarely changed settings (prevalence, star mix and the sector minimums) sit
in a folded Advanced block on their task's page, open when a draft that holds
one is shown again.

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
overrides, a pasted `--system-file` JSON and the
`--debug` narration. It runs `planetgen system --format document --output FILE` in a
temporary directory and waits for it (a system takes about a second), so
nothing touches the database. The result shows as the system's description
(the Document, `describe/html.py`) with a Download button for each format
(Markdown and text; a button posts the Document back to
`/admin/generate/system/download`, which returns it as a `.md`/`.txt`
file).

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
from Admin, Queue. Only one job runs at a time: the job lock, a hash in Redis
(`planetgen:job:active`), is written with the job id and the time in
one step and a stale lock is cleared only if it still names the same job,
so two admins starting a job at the same moment start exactly one
(TEST.40). A job that finds another running waits in the one waiting list
(`queue.max_waiting` entries, 8 by default; PERF.154) and starts by itself. A lock naming no readable job is treated
as a job still starting for `STARTING_GRACE_SECONDS`, then as stale.
Pruning old jobs never removes the one the lock names.
Jobs live under `jobs.dir` (`docs/reference/config.md`). A visitor who isn't a
logged-in admin is sent to the login page, and POSTs and status requests
without an admin session get a 403.

**Error pages.** Any failure of a page is an HTML page (`web/errors.py`) with the right status: 404 for a bad id or an unknown URL, 405, 413 and 429 for what they say, 502 "Data unavailable" when the API failed, and 500 with a fixed sentence (the traceback goes to the log, never the page). Two more reload themselves after `Retry-After` seconds, so a page opened while a long job holds the database fills in by itself (PERF.64, ADM.53): **503 "Busy"** when a request waited too long for a row a long job holds locked, and **504 "Took too long"** when a database read hit the statement time limit (`mysql.statement_timeout_seconds`). Under `/api` the same two answers are JSON, `{"error": "BUSY: ..."}` and `{"error": "QUERY_TIMEOUT: ..."}` (see [`api.md`](api.md#errors)), and `/table/<name>` answers a timeout as JSON too.

### How a Flask page is built

```
src/planetgen/web/
  __init__.py         blueprint `web`, template globals, init_app()
  views.py            /, /sectors, /systems, /search
  system_pages.py     /system/<id>, /phenomena, /phenomenon/<type>/<id>
  object_lists.py     /lists, /lists/<kind>
  preferences_page.py /preferences (UX.236)
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
  credentials, create or revoke an API key) is
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
- `/admin/stats` uses the site's configured
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
| `galaxy_tiles` | The map's data and scene routes, fetched by its script: `/galaxy/nodes`, `/galaxy/stage`, `/galaxy/window`, `/galaxy/near`, `/galaxy/territories`, `/galaxy/uncharted/.../scene`, `/galaxy/nebula/...` and `/sector/<id>/scene` | 600 per minute |
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

[`../INSTALL.md`](../../INSTALL.md) walks through a first install, and
[`deployment/README.md`](../deployment/README.md) compares every supported
platform (Apache, nginx or Caddy on Linux; macOS; a VPS) and links a guide for each. The reference setup is
Apache2 with mod_wsgi on Debian or Ubuntu
([`deployment/apache.md`](../deployment/apache.md)), where `sudo
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
see [`deployment/apache.md`](../deployment/apache.md#planetgens-own-python)), the
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
