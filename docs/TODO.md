# planetGen Roadmap

## How to use this document

This file tracks **open work only**. Finished work is not listed here:
`CHANGELOG.md` and git history are the record. For how things work today,
read the code and the reference docs (`src/stellarObjects/schema.sql` and
`docs/database-schema.md` for the schema, `docs/api.md` and
`docs/html-interface.md` for the web interface,
`docs/design/galaxy-coordinate-system.md` for galaxy geometry).

Each open item says what's wrong (the observed symptom), where to start
looking, and what "done" means, so it can be picked up cold. Items with
open questions for Boss state the default taken; the work can start on
that default.

### IDs

- Every item has an ID made of its category and a number, `CAT.N`, for
  example `MAP.16`. The number is a plain count within the category,
  like the database schema version: the first item in a category is
  `.1`, and each new item, subitem or bug takes the next number. There
  are no deeper levels (no `MAP.2.1`); an item's place in the tree is
  shown by indenting it under its parent, not by its ID. The categories
  are in the table below.
- IDs are permanent. Nothing is renumbered, and an ID is never reused.
  When an item ships, delete it here (the changelog is the record) and
  describe the change in the PR's `changes/` note; its ID stays taken. A
  finished item that still has open subitems stays as a checked line so
  its subitems keep their place.
- A new ID is one more than the highest ever issued in its category,
  finished items included.
  [design/todo-number-map.md](design/todo-number-map.md) lists every ID
  issued so far and the next free number in each category. It also maps
  the old running numbers (1 to 111, used until 2026-10-01 in the
  changelog, commits, PR titles and older docs) and the dotted tree IDs
  used briefly on 2026-10-01 (such as `MAP.2.1.1`, in PRs #189 and
  #190) to these IDs.

### Where new work goes

- A new feature is a new item, or a subitem (indented under the item it
  extends, with its own next number).
- A bug is a subitem of the item whose feature it breaks, marked
  "(bug)". A bug in something with no item here is a top-level item in
  its category, marked "(bug)" (for example UX.15).

### Links

- This file is the only place that links an item to its design document,
  on its own line under the item:
  `Design: [docs/design/x.md](design/x.md), section N`. Design documents
  don't list item IDs.
- Code tags cite the ID: `TODO(MAP.16): what to do here`. Grep
  `"TODO(MAP."` for an area, or `"TODO(MAP.16)"` for one item.
- `changes/` notes, PR titles and commit messages cite the ID, for
  example "Galaxy Map: stages (MAP.16)".

### Categories

| ID | Covers |
|---|---|
| UX | Web pages, the site shell, menus, page content |
| MAP | Galaxy Map, Sector Map, System Map |
| NAV | Navigation and travel (NAV page, courses, speeds) |
| GEN | Generation and physics (stars, systems, phenomena, galaxy model) |
| PERF | Speed, caching, bulk generation and parallel work |
| DB | Schema and migrations |
| API | The JSON API |
| ADM | Admin tools: editing, overrides, Generate page |
| SEC | Security |
| TEST | The test suite itself: new tests, test infrastructure, CI test jobs |
| USR | User accounts |
| OPS | Installers, hosting, CI, releases |
| DOC | Documentation and project process |
| VIEW | The view from a planet |
| POP | Population and politics |

NAV, DB, API, SEC, OPS, DOC and POP have no open items today.

## Background

The long-term goal is a fully populated galaxy (every sector, system,
planet, moon and belt) stored in MySQL and browsed through a web
interface served from Apache. The original roadmap phases are all done:
object-graph serialization, the relational schema and migrations, CLI
tools writing to the database, lazy galaxy-scale generation from a
density skeleton, and the Flask API (`src/html/api/`) with server-rendered
pages served by the same app (`src/html/web/`): Galaxy, Sector and System
maps, search, NAV, admin auth and wiki publishing.

## Plan: what to do now

The bug round (MAP.17 to MAP.19, MAP.26, MAP.37, MAP.43 to MAP.51,
UX.15, UX.16, UX.19, UX.20, ADM.9), login security (SEC.1, SEC.20 to
SEC.28), the database-call work (PERF.6, PERF.12 to PERF.17) and
parallel generation (PERF.7, PERF.8) are done. Boss (2026-10-01
14:41Z) set the next round, run as parallel threads:

1. **First:** PERF.5 (scatter bright stars in stages).
2. **Galaxy Map and units:** UX.13, UX.14, MAP.30 (now a layer slider to
   the right of the Galaxy Map), MAP.15, GEN.23 and what is left of
   MAP.2.
3. **Generation estimates and progress:** PERF.3 and PERF.4, then PERF.9
   and PERF.10.
4. **Admin editing:** ADM.1 and its subitems, on the validate module
   (ADM.5, done in PR #235).

Waiting behind those: PERF.11, UX.2, UX.3, ADM.4, GEN.8, GEN.9,
bookmarks (MAP.23, which finishes MAP.22) and user accounts (USR.1,
starting with roles, USR.2). View from a planet (VIEW.1) waits on a
research session with Boss, except the constellation names (VIEW.4).

## UX: Web pages

Boss's UX reference for this section is "Responsive Web Design
Standards" (Boss's notes of 2026-10-01; a copy is in the project's
shared files at `ux-standards/`). The rules from it that apply here:
size classes by width (compact under 600 px, medium 600-839, expanded
840-1199, desktop 1200+), not by device; viewport media queries only for
the top-level frame and container queries for everything inside it; on
medium and wider screens a menu must not stretch across the screen like a
phone's; text columns capped at 45-75 characters, but a map or canvas may
use the full width; touch targets at least 44-48 px on coarse pointers
(`pointer: coarse`), smaller is fine for a mouse; spacing and type sized
with `clamp()`.

- [ ] **UX.2 Menus sized to what they hold**
  Boss (2026-10-01): "I want the
  menus to be proportional to the size needed, I noticed on tablet
  screens that menu acts like a phone screen spanning absurdly across
  the screen." The Menu and gear drop-downs (`.site-menu-panel`,
  `.site-gear-panel` in `static/style.css`) have `min-width: 14rem` and
  `max-width: calc(100vw - 1rem)`, and the header folds the section
  buttons into the Menu below 43rem, so on a tablet the panel can grow
  to nearly the full screen. Done: each panel is as wide as its longest
  entry plus padding (capped, for example `width: max-content` with a
  sensible `max-width`), full width only on compact (phone) screens;
  checked at 390, 600, 768, 820, 1024 and 1280 px in both themes and
  both orientations, with touch targets still at least 44 px on touch
  screens.

- [ ] **UX.3 Warn every visitor while a background job changes the galaxy**
  Boss (2026-10-01): "a warning to all users on the UI when a
  task is running in the background which is modifying the starmap is
  going on with an ETA until it will be finished to the nearest hour
  rounding up." Today only the admin Generate page shows a running job
  (`web/jobs.py`: the `active` lock in the jobs directory, each job's
  `state.json`, and `progress.json` written by `generate.py` through
  `stellarObjects.progressFile`). Done: while a job that writes to the
  galaxy is running (plan, bright-star scatter, sector fill, block and
  neighborhood generation from the map (MAP.20), reset, new galaxy,
  and later ADM.8's regenerate), every page shows a banner to every
  visitor, signed in or not, saying the galaxy is being changed and
  when it should finish, as an ETA rounded up to the next whole hour
  (from `progress.json`, and from PERF.3's measured stars-per-second
  rate). Open questions: a banner at the top of every page, or only on
  the map and list pages? Does it also cover `generate.py` runs started
  from the command line, which don't take the web jobs lock today
  (they would need to write the same lock and progress file)? What does
  it say when there's no ETA yet?

- [ ] **UX.13 One meaningful-unit ladder for speeds**
  Boss (2026-10-01):
  "standardization with speeds, similar to what we do with distances
  to make them always meaningful, speeds should always be meaningful.
  Going from km/h on the low speed end to mm/s on the high speed end."
  Distances already work this way: every page passes them through
  `stellarObjects.utils.format_distance_m` (and its `_km`/`_au`/`_ly`
  wrappers, `html/lib/fmt.py`), mirrored by `static/distance.js` for
  the maps, which picks the largest unit the value is at least 1 of.
  Speeds have no such function; each page formats its own (km/s on
  system and facility pages, multiples of c for warp and fold in
  `stellarObjects/navigation.py`'s `warp_speed_c` and travel tables).
  Done: one shared speed formatter in Python with a JavaScript mirror,
  with a fixed ladder of units, used by every page, API text field
  and map that shows a speed, and the existing call sites converted.
  Open questions: the ladder reads reversed as written (mm/s is slower
  than km/h), so what is the intended order from slowest to fastest?
  For example mm/s, m/s, km/h, km/s, then fractions and multiples of
  c, with warp and fold factors shown alongside rather than replacing
  them. Where it switches units (at 1 of the next unit, as distances
  do, or another rule), and whether it adds a parenthetical in a
  second unit the way distances add ly or AU. Where it lives
  (`stellarObjects/utils.py` next to `format_distance_m`, and a
  `static/speed.js` or a section of `distance.js`).

- [ ] **UX.14 One meaningful-unit ladder for time periods**
  Boss
  (2026-10-01): "Same for orbital periods, galactic, lunar, planetary,
  we should tie all those into a function to do the same. For slowest
  (measured in Gy) to fastest (measured in microseconds). Those are 2
  seperate TODO items." Today orbital periods go through
  `stellarObjects.utils.years_to_time_string` ("x years y days z hours
  m minutes", via `html/lib/tabledisplay.format_period`), which gets
  long for galactic orbits and loses anything under a minute; star
  ages and lifespans are shown in Gy elsewhere, and the admin pages
  have their own `format_duration`. Done: one shared period formatter
  in Python with a JavaScript mirror, picking a meaningful unit from
  Gy at the slow end down to microseconds at the fast end, used for
  planetary, lunar and galactic orbital periods (and rotation periods,
  ages and other durations where it fits), with the existing call
  sites converted. Open questions: the ladder (for example Gy, My, ky,
  years, days, hours, minutes, seconds, ms, µs) and where it switches;
  one unit with decimals ("1.88 years") or a mixed form ("1 year 321
  days") for everyday periods; rounding and significant figures; which
  year length it uses (Julian 365.25 days, as `years_to_time_string`
  does); and whether elapsed-time and ETA displays for jobs (PERF.3,
  UX.3 and PERF.4) and the admin pages use the same function.

- [ ] **UX.21 Clean up the web interface: overlapping buttons and dead controls**
  Boss (2026-10-01 14:58Z): "clean up the web interface, still have
  buttons overlapping, we have +/- buttons that don't do anything
  anymore, etc. Don't start it yet, but it needs to be done." UX.16
  (PR #195) added space between buttons, but some still overlap, and
  some controls survived the map rewrites without anything wired to
  them. Done: a pass over every page (public, admin, and the Galaxy,
  Sector, System and phenomenon maps) at each size class (compact,
  medium, expanded, desktop) in light and dark: no buttons or labels
  overlap or run off their panel; every button, toggle and link does
  something, and dead ones (such as the +/- zoom buttons Boss saw) are
  wired up or removed, along with any hint text that names them;
  controls follow the section's rules above (touch targets, spacing).
  Before/after screenshots of each fixed spot go with the PR. Known
  dead control (found while planning the TEST items): on the nebula and
  supernova remnant diagrams, anything about half a light-year across
  or larger opens already at the 1 ly zoom-out limit
  (`lib/phenomenonmap.py` lines 53 and 127, the clamp in
  `static/mapzoom.js`), so "-" has nowhere to go. The Galaxy Map has no
  +/- buttons, and the Sector Map's buttons work. TEST.55 (every map button
  changes the view) and TEST.56 (no overlapping controls at 390 to
  1280 px) pin this item.

- [ ] **UX.22 Meaningful units for every measurement**
  Boss (2026-10-01 15:10Z): "standardize ALL measurements into trees
  like we have so that we always have meaningful units. From mass, to
  distance, to time, to speed, just everything that can have units.
  Atmospheric pressure and surface conditions should show customary
  units as well as a secondary to help contextualize the metric values
  given." Today only distance has a ladder (`format_distance_m` and
  friends in `stellarObjects/utils.py` and `html/lib/fmt.py`, mirrored
  by `static/distance.js`), with speed (UX.13) and time periods (UX.14)
  being built. Done: one ladder per quantity, in Python with a JavaScript
  mirror, picking a meaningful unit the same way, and every page, map
  panel and form converted to it: mass (kg, Earth, Jupiter and solar
  masses), distance, time, speed, temperature, pressure, gravity,
  density, luminosity, power and any other quantity shown with a unit.
  Surface conditions show temperature in K, °C and °F; atmospheric
  pressure and the other surface conditions show a customary unit
  (such as atm, psi or g) beside the metric value. The UX thread was
  asked (2026-10-01 15:10Z) to add the K/°C/°F temperature display now;
  this item covers the rest. Open questions: the ladder and switch
  points for each quantity; which customary unit goes with each surface
  condition; whether the secondary unit shows in tables or only in
  detail panels.

## MAP: Galaxy Map, Sector Map, System Map

- [ ] **MAP.2 Drill-down navigation**

  Design: [docs/design/galaxy-drilldown-navigation.md](design/galaxy-drilldown-navigation.md)

  Boss's design of 2026-10-01: the Galaxy Map becomes a drill-down. In
  3D, pick a slab; it is pulled out and shown from above; pick a block;
  its contents fill the view as blocks 1/9 the size; repeat until single
  sectors, where a click opens the sector. The ladder is 243 -> 27 -> 3
  -> 1 sectors a side ("the bigger targets"), eight clicks from the
  galaxy to a sector. Admins can generate a sector, a layer or a
  neighborhood (radius asked in light-years) at the sector level, and
  the NAV page can pick its start and destination on the map or in a
  sector. Everything below is specified, with the math, in the design
  doc; each subitem names its section. Shipped: the nested ladder
  (MAP.28, `stellarObjects/galaxyDrill.py`), the stage contents API
  (MAP.29, `GET /api/galaxy/stage`), the stages (MAP.16,
  `static/galaxystages.js`, `static/galaxystageview.js`, stage URLs
  `/galaxy?at=&p=` and `?sector=<designation>`), the Sector Map pick
  mode (MAP.21), the address bar (MAP.24, `/galaxy/locate`), the course
  overlay (MAP.27, `/galaxy?course=<from>,<to>`), generating from the
  map (MAP.20), the "Show on Galaxy Map" links (MAP.25) and the top-down
  drill-down (quarter, layer, arc, ..., sector; MAP.17, MAP.18, MAP.19,
  MAP.44, MAP.26). Since Boss's change of 2026-10-01 09:16Z the camera
  is locked top-down only at the full galaxy and its quarters, and turns
  and moves freely from an arc down (PR #219).

  - [ ] **MAP.22 NAV page picks on the map**

    Design: [docs/design/galaxy-drilldown-navigation.md](design/galaxy-drilldown-navigation.md), section 9

    Boss: "from the nav
    menu select start and destination using either the text dropdowns as
    we have now or the galactic map interface to select. If it's within
    sector then it'll just use the sector interface." Done: beside each
    dropdown, "Pick on Galaxy Map" (`/galaxy?pick=...`, generated-only
    forced on, ending in MAP.21's Sector Map pick mode), "Pick in this
    sector" once the other end is known, and a Bookmarks select. The
    Galaxy Map's side is in: `?pick=` shows the banner with Cancel back
    to NAV, keeps "Generated only" on, and a sector click opens that
    sector in pick mode. The NAV page's side is in too: "Pick on Galaxy
    Map" and "Pick in this sector" at each step. Only the Bookmarks
    select is left, and it waits on MAP.23.

  - [ ] **MAP.23 Bookmarks**

    Design: [docs/design/galaxy-drilldown-navigation.md](design/galaxy-drilldown-navigation.md), section 8.2

    Done: a ☆ on the breadcrumb and
    info panels saves a stage, sector, system or phenomenon in
    `static/bookmarks.js` (per browser, up to 100, storage failures
    tolerated), with a map menu, Ctrl+1-9, rename and delete, and the
    entries offered by the NAV pickers. Shared bookmarks need Boss's
    decision 4 and a migration.

- [x] **MAP.3 A bigger Galaxy Map with controls underneath**
  Done in 7.55.0 (PR #178); kept as the parent of MAP.30.
  - [ ] **MAP.30 Slab list to the left of the map, and a 3:4 map**
    Boss (2026-10-01): "slab selection goes to the left of the galactic
    map if there is room, given the galactic map shoul dhave a 3:4 aspect
    ratio to its window or 1:1 if necessary, like in mobile view
    perhaps." Today the Galaxy Map's viewport (`.galaxymap3d-panel
    .starmap-viewport` in `static/style.css`) is the full width with a
    height of `min(100svh - 9rem, max(20rem, 75vw))`, and the stage
    view's slab list (`.galaxy-slab-row` rows, `static/galaxystageview.js`)
    sits with the other controls below the map. Done: the map keeps a
    3:4 aspect ratio within its window, falling back to 1:1 where the
    window can't fit 3:4 (as on phones); when there's room beside the
    map, the slab list moves to its left; when there isn't, it stays
    under the map with the other controls; no layout jump as stages
    change. The list is now the layer strip of the top-down drill-down
    (MAP.17): it offers layers when a layer is to be picked
    and follows the Responsive Web Design Standards notes in the UX
    section (size classes, container queries). Open questions: is 3:4
    width to height (taller than wide) or height to width (4:3, wider
    than tall, close to today's 75vw height)? At what width does the
    list move to the left (the standards' expanded class, 840 px and up,
    or whenever a readable list column fits)? Does the rest of MAP.3's
    controls row stay under the map, or join the list on the left?

- [ ] **MAP.15 Stars and glowing phenomena as points of light on the Sector Map**
  Boss (2026-10-01): "make the stars in a sector more realistic
  sizes with bright auras, I prefer the tiny point of light in the map
  for stars, same for any stellar phenomena which has a glow / emits
  light, it should be a point of glowing light. Still make brightness
  and size relevant just more like a point of light." Today
  `static/sectormap.js` draws each star (and black hole, neutron star,
  quasar, rogue planet, interstellar comet) as a textured sphere with a
  fresnel glow shell (`bodyRendering.js`), sized in scene units so it
  looks like a ball; the Galaxy Map draws its bright stars as a tiny
  core with a soft halo at a fixed pixel size (`galaxymap3d.js`,
  "Bright stars", `STAR_MIN_PX` to `STAR_MAX_PX`). Done: on the Sector
  Map, stars and every light-emitting phenomenon draw as a tiny bright
  point with a glowing aura, like the Galaxy Map's bright stars, with
  brightness and halo size still scaling with luminosity (and color
  with spectral type), so they read as points of light rather than
  balls; clicking one still selects it. Open questions: does the size
  stay fixed in pixels at every zoom (as on the Galaxy Map) or grow a
  little as the camera gets close? Is the textured sphere kept for a
  close-up (the System Map still shows bodies as spheres)? How is a
  binary's pair kept distinguishable? Which phenomena count as
  emitting light (quasars, neutron stars and accreting black holes
  yes; quiescent black holes and rogue planets, which MAP.46 must
  keep findable, probably not)?

## GEN: Generation and physics

- [ ] **GEN.8 Give rogue planets a planet class, with a rogue flag in the class constants**
  Boss (2026-10-01): "Rogue plants should get a
  planet class, add TODO item to TODO.md that we should add to the
  zone data for planet class constants a flag for if a planet is
  acceptable to be rogue or not (zone r for the purposes of the
  constants)." Today each class in `program_constants.PLANET_CLASSES`
  carries zone flags `"h"`, `"e"` and `"c"` (hot, ecosphere and cold
  zones), and a rogue planet (`stellarObjects/roguePlanetData.py`,
  `RoguePlanet`) has no class: just a `planet_type` of `'t'` or `'g'`
  picked by mass from the rogue mass bins. Done: every class in
  `PLANET_CLASSES` gets an `"r"` flag saying whether it can be a rogue
  planet; a rogue planet is given a class drawn only from the classes
  with `"r": True` that fit its mass and type; the class is stored and
  shown on the rogue planet's page and the Sector Map like any other
  planet's; and the class override and validation items (ADM.5 to
  ADM.7) treat `"r"` as the rogue planet's zone. Open questions: which
  classes are allowed to be rogue (frozen, gas giant and barren classes
  are the obvious ones; does a class with life ever qualify)? Are the
  probabilities `PLANET_CLASS_PROBABILITIES` reweighted for rogues, or
  a separate rogue table? What happens to rogue planets already
  generated: a migration that assigns classes from their stored mass
  and type, or a regenerate? Does a rogue class change its rendering
  (it has no star to light it)?

- [ ] **GEN.9 Plan for more than one galaxy in the database**
  Boss (2026-10-01): "Lay the groundwork for different galaxies within
  the same DB, this is a planning item. Develop a plan to have multiple
  galaxies first as objects in the neighborhood (i.e. we can see
  Andromeda from Earth kind of thing) but also I might need to have
  another galaxy at some point. So just plan and save as a planning
  document." Today the database holds one galaxy: one galaxy skeleton
  (`galaxy_shape`, `galaxy_layer`), one sector grid, and coordinates in
  that galaxy's own frame. This item is a plan only, no code. Done: a
  planning document in `docs/design/` (linked here on a "Design:" line
  once it exists) covering two stages: first, other galaxies as
  objects seen from this one (direction, distance, size, brightness and
  type, so a planet's sky (VIEW) or a map can show them, like Andromeda
  from Earth); later, a second fully generated galaxy with its own
  skeleton, sectors and systems. It says what each stage needs in the
  schema (a galaxy id on which tables), the frames and coordinates
  between galaxies, the URLs and pages, generation and the Galaxy Map,
  and how an existing single-galaxy database migrates. Open questions
  for the plan: are neighboring galaxies real (the Local Group's
  catalogued galaxies) or generated? Does every row get a galaxy id, or
  only the top-level ones (sectors, the skeleton)? Is a second galaxy
  in the same database or a second database chosen at login?

- [ ] **GEN.23 Generate a smaller sphere, then backfill bright stars around it per sector block**
  Boss (2026-10-01): "right now the way it's set up, the default is
  about 100 light-years around a new sector. When the galaxy is
  generated I want to reduce that and make the default sector
  generation sphere about 10 pc, rounded up. Then use the
  100-light-year default to rerun the large star generation on the
  100-light-year area around the generated sector, skipping over
  sectors that already have stars of the required brightness or lower.
  When this happens we're going to generate all the stars that are 100
  luminosities or larger, up to what has already been generated. Again
  it's going to be bound by the lower and upper limit. This means we're
  going to have to extend the star generation storage for the lowest
  level of stars generated for a region. For the sector blocks we'll
  group to make the calculations easier and such. We won't do this per
  sector; we'll do this per sector block, each containing a 3x3 cube.
  For that block we will store the star minimum and maximum generated.
  We'll store it, say, if all stars 100 solar luminosities or higher
  have been generated. I guess we only need to store the data
  pertaining to the dimmest luminosity already generated. It will skip
  over any sectors that have already been filled. This is only for the
  logic for unfilled sectors surrounding the filled sector space."
  Today `program_constants.RANDOM_START_NEIGHBORHOOD_RADIUS_LY` (100 ly,
  about 30.7 pc) is the sphere filled by `generate.py galaxy`'s random
  start, and the default of `POST /api/sectors/<id>/generate-neighborhood`
  (the Generate page and the map's neighborhood generation, MAP.20);
  `--center-sector` has no default and requires `--radius-pc`. The
  plan's bright-star scatter (`brightStars.scatter`) draws every star at
  or above one galaxy-wide threshold
  (`BRIGHT_STAR_MIN_LUMINOSITY_SOL`, default 500 L_sun, stored in
  `galaxy_shape.bright_star_min_luminosity_sol`) into `bright_stars`.
  Done: the default generation sphere is about 10 pc, rounded up; after
  that sphere is filled, a bright-star backfill runs over the 100 ly
  sphere around the same center, adding only the stars from 100 L_sun up
  to (not including) the luminosity already scattered there, so no
  brighter star is drawn twice; the backfill works per sector block (a
  3x3x3 cube of sectors) and stores, for each block, the dimmest
  luminosity already generated in it, so a block already down to 100
  L_sun is skipped; filled sectors are skipped, since this is only for
  the unfilled sectors around the filled space; and the new stars show
  on the Galaxy Map like the existing bright stars. Open questions:
  - "10 pc, rounded up": to the next whole sector (sectors are 4 pc, so
    12 pc), to whole sector blocks, or to a round number of light-years?
    Does the smaller default apply to every generate-around entry point
    (random start, the API and the Generate page, `--center-sector`
    gaining a default) or only to the random start?
  - Where the per-block dimmest luminosity is stored: a new table keyed
    by block (ring/layer/slot of the block's corner or center, and a
    schema version bump), a column on an existing block table if there
    is one, or derived from `MIN(luminosity)` of `bright_stars` in the
    block (which can't tell "no star that dim landed here" from "not
    generated yet")? Boss asked for the minimum and maximum, then only
    the dimmest; is the maximum needed at all, since the upper bound is
    the galaxy-wide scatter level?
  - How this fits PERF.5 (staged bright-star layers, a galaxy-wide
    star-fill level): is the per-block level PERF.5's level stored per
    block instead of per galaxy, so a galaxy-wide "go down to N" and
    this local backfill share one mechanism? A block's level could then
    be lower than the galaxy's, and the galaxy-wide pass skips it.
  - "Up to what has already been generated": the upper limit is the
    galaxy's current scatter level (500 by default), or the block's own
    level when it has one? And "the lower and upper limit" means the
    100 L_sun floor and that level?
  - Is 100 L_sun fixed, a constant in `program_constants`, or a config
    or Generate page option? Is the 100 ly backfill radius the existing
    `RANDOM_START_NEIGHBORHOOD_RADIUS_LY` reused, or its own constant?
  - Must the backfilled stars come from the same random stream as a
    galaxy-wide scatter (PERF.5's question), so a block backfilled now
    matches the stars a later galaxy-wide pass to 100 L_sun would draw?
  - A block that is partly filled: backfill only its unfilled sectors
    and record the block's level anyway, or leave its level unset so a
    later pass knows the filled sectors never got the dimmer stars?
  - Does the backfill run inside the same job and progress bars as the
    sphere's generation (PERF.4, the Generate page), and do PERF.3's
    size and time estimates include it?

- [ ] **GEN.24 Generate the galactic core on layer 0**
  Boss (2026-10-01): "Add a new TODO item to TODO.md (don't start
  implementing yet) so that we have an option to create the galactic
  core sectors (everything in the central core of the galaxy layer 0
  (i.e. middle)." Today layer 0 is the 4 pc slab centered on the
  galactic plane, and ring 0, layer 0, slot 0 holds the galaxy's
  nucleus (the quasar or quiescent supermassive black hole,
  `add_galactic_nucleus`). `generate.py galaxy` can fill one ring on
  one layer (`--ring I --layer J`), a ring through every layer
  (`--shell`), a block, or a sphere around a sector, but nothing fills
  "the core" in one go. Done: one option on the CLI and on the admin
  Generate page generates every not-yet-generated sector of the core on
  layer 0, nucleus sector included, with the same size warning,
  confirmation, `--limit` and progress as the other bulk modes (PERF.3
  estimates). Boss's answers (2026-10-01 14:55Z):
  - Size: the admin chooses. The Generate page offers "core" as one
    more mode alongside the others, where Boss types the size he
    wants, as a radius or a number of rings (the CLI takes the same).
    Suggested default: the bulge scale radius, 200 pc, about 50 rings
    and roughly 7,900 sectors on layer 0.
  - Layer 0 only, not every layer the bulge reaches.
  - GEN.23's smaller generate-around sphere and bright-star backfill
    run around the core too, as they do around any generated sector.
  - Fill from ring 0 outward, so a run that stops early (or hits
    `--limit`) still leaves a solid disc around the nucleus.

- [ ] **GEN.25 A moon reclassified after its planet moves can be too large for its planet (bug)**
  Found by the ADM.1 thread with ADM.5's validator (PR #235):
  `stellarObjects/validation.check_star_system` reports "moon too large
  for its planet" on about 3 of 1,000 generated systems with moons.
  Start in `validation.reconcile_moved_planet` and
  `planetPhysics.reconcile_zone_and_class`, which re-roll a moon's
  class without checking `max_moon_radius_km` (planet radius /
  10^(1/3)) or mass <= planet mass / 10. Done: 1,000 generated systems
  pass `check_star_system` with no moon-size problems.

- [ ] **GEN.26 Rogue planet surface conditions**
  Boss (2026-10-01): "I also want to calculate surface conditions,
  knowing they will be extremely cold with no star to warm the
  surface", with his pasted research: an energy balance with internal
  heat flux plus the cosmic microwave background, radiogenic and
  primordial heat, and three outcomes (frozen atmosphere, hydrogen
  envelope, ocean under an ice lid), with adiabats for gas giants.
  Being built by the GEN.8 thread as schema v48. `has_internal_heat`
  stops being a 40% roll (`ROGUE_PLANET_INTERNAL_HEAT_CHANCE` goes): it
  is computed, true when heat flow is at least 0.04 W/m2 and always for
  giants. Its design document (`docs/design/rogue-planet-surface.md`)
  arrives with its PR.

- [ ] **GEN.27 Class P (glaciated world) only in the habitable zone, and fitting there**
  Boss (2026-10-01 15:26Z): "make sure our frozen world, Class P, only
  appears in the habitable zone and adjust so that it fits there." P is
  already ecosphere-only (`h` False, `e` True, `c` False in
  `program_constants.PLANET_CLASSES`) and placed at `zone_position_mode`
  0.90; sampled P worlds run 198-211 K. Done: no path (generation,
  `reconcile_zone_and_class`, moon regeneration, admin overrides) can
  put P outside the ecosphere, and its placement, albedo, greenhouse
  and atmosphere are tuned so a glaciated world with life is
  consistent at the outer edge of the habitable zone. Evidence: the
  planet class gap report
  (https://claude.ai/code/artifact/549f0ba8-ca6f-4d35-be2d-0c7591b93256).

- [ ] **GEN.28 Seven new planet classes in the letter gaps (R, S, U, W, X, Y, Z)**
  Boss (2026-10-01 15:26Z): "add the other 6 classes filling in the
  letter class gaps sequentially. For the subsurface ocean moon, split
  this so that we have the one the size of a class D (moon / pseudo
  planet) and one similar to a terrestrial world (a modification of
  Class C), for lifeless temperate world, don't we have a class for
  that already? If not, I approve adding one." There is none: every
  ecosphere rocky class with air carries life, and C (the only lifeless
  rocky one) is airless. Proposed mapping, in letter order (the build
  can adjust):
  - R Sub-Neptune: rock and ice core under a hydrogen-helium envelope,
    1.8-4 Earth radii, 3-20 Earth masses, hot, ecosphere and cold. The
    most common real planet type, missing today; consider whether T
    (gas dwarf, 0.05% of planets) merges into it.
  - S Rocky super-Earth: barren, 1.2-1.8 Earth radii, 2-10 Earth
    masses, hot, ecosphere and cold (a hot one can be a lava world). V
    stays the life-bearing super-Earth.
  - U Icy world (ice dwarf or large icy moon): water ice over rock,
    500-3,000 km, cold (Ganymede, Callisto, Triton, Pluto, Eris).
  - W Small subsurface ocean body, Class D sized (moon or pseudo-planet,
    about 50-500 km): Enceladus analog.
  - X Subsurface ocean world, terrestrial sized (a modified Class C,
    about 500-10,000 km): Europa analog and larger.
  - Y Titan-like world: thick nitrogen atmosphere, methane rain,
    hydrocarbon lakes, 1,500-4,000 km, cold.
  - Z Lifeless temperate world: rocky, with an atmosphere, ecosphere,
    no life.
  Done: each class has zone flags, weights, radius and mass ranges,
  atmosphere, moon eligibility and a GEN.8 rogue `"r"` flag decision,
  and shows on the class reference pages. R, W, X and Y belonged to
  classes removed in early September, so no old rows or tests may
  still expect those letters. Evidence: the planet class gap report
  (link in GEN.27).

- [ ] **GEN.29 Sweep every planet class for sense once the new ones are in**
  Boss (2026-10-01 15:26Z): "do a full sweep of planet classes to make
  sure they all make sense logically once the new classes are in
  place." After GEN.28. Done: every class's description, composition,
  zones, sizes, weights, temperatures and atmosphere agree with each
  other. Known oddities to settle: D allowed in the hot zone (icy
  bodies at 265-490 K); C a catch-all for 63% of cold planets; Q never
  generated (weight 0.0001, and orbits are circular); V's composition
  "iron, iridium, tungsten"; L with vegetation at a median 0.02 bar; E
  at 376-414 K, above water's boiling point at 0.6 bar.

## PERF: Speed, caching, bulk generation and parallel work

- [ ] **PERF.1 Generation at scale**
  Boss's notes of 2026-10-01 on bulk generation: estimates before it
  starts, progress while it runs, speed records, and parallel work. The
  code is mostly `generate.py`, `web/generate_page.py`, `web/jobs.py` and
  `stellarObjects/brightStars.py`.

  - [ ] **PERF.3 Estimate size and time before bulk generation, and refuse what won't fit**
    Boss (2026-10-01): "any directive to generate
    sectors in bulk should have a size estimate calculated +10% and make
    sure that it warns the user approximate size of the generated content
    before it generates it. Same for time, which means we'll have to
    store in the control database somewhere how fast stars are generated
    in exact terms as we can, which we'll get from the generate phase,
    not from the plan, but from when we start generating sectors, have
    the program keep track of it's average stars per second and then
    we'll use the expected stellar density of all the sectors being asked
    to generate to determine the amount of time it is estimated to take.
    Also it should refuse to generate anything that would take more than
    1/4 of the total disk space OR would leave less than 5 GB of space
    estimated to be left." Today nothing estimates size or time up front:
    `generate.py` only shows elapsed time and a running ETA once
    generation has started, and `stellarObjects/generationLimits.py`
    caps input sizes (radius, ring counts, orbits), not output size.
    Done, for every bulk path (`generate.py galaxy`/`sector` over many
    sectors, the admin Generate page's jobs, the Galaxy Map's block,
    layer and neighborhood generation (MAP.20), and ADM.8's sector
    regenerate):
    - Before starting, compute the expected star count from the expected
      stellar density of the requested sectors, then an estimated size
      (bytes per star from measured data, plus 10%) and an estimated
      time (stars / measured stars per second), and show both to the
      user, who confirms before anything is written.
    - Measure speed during sector fill, not during the plan: the program
      tracks its average stars per second while generating sectors and
      stores it in the control database (a new table or row next to
      `admin_users` in `control_schema.sql`), updated after each run.
    - Refuse when the estimate would use more than 1/4 of the total disk
      or leave less than 5 GB free, saying why and how much space it
      needs.

    Open questions: which disk counts (the MySQL data directory's volume,
    which may be another machine, or the planetGen host); how bytes per
    star are calibrated (measured from the database's own table sizes
    after each run, or fixed from the 2026-09-30 galaxy-size study in the
    project's `galaxy-studies/`); what rate to use before any run has
    been measured (a conservative default?); whether the rate is kept per
    server or per kind of sector (bright-star and bulge sectors cost
    more per star); and whether an admin can override the refusal.

  - [ ] **PERF.4 A second progress bar for slow layers in the plan**
    Boss
    (2026-10-01): "For building the layers, when the rate drops below 1
    layer per 30 seconds, which is calculated on every star system
    generated, then it should double up the progress bars intelligently
    (they should still not cause any flicker and should stay at the
    bottom while the text above scrolls) that has the ETA for the layer
    being generated based on the estimated number of stars remaining.
    This means the system should estimate the stars remaining every time
    a new star is generated during the plan." Today the plan's bright-star
    scatter (`generate.py`, the "Bright stars (layers)" task) shows one
    bar that moves once per finished layer (`brightStars.scatter`'s
    `on_layer` callback), so a slow layer looks stalled. The display is
    `_generation_progress()` (rich `Progress`, with every log line routed
    through `progress.console` so the bars stay pinned at the bottom
    without flicker). Done: after every star the plan generates, it
    recomputes the layer rate and the current layer's estimated stars
    remaining; while the rate is below 1 layer per 30 seconds, a second
    bar appears under the layers bar showing the current layer's stars
    done of its estimate with that layer's ETA, and it goes away when the
    rate recovers; no flicker, bars stay at the bottom, log lines keep
    scrolling above, and the web job's `progress.json` carries the same
    numbers. Open questions: what "estimated stars remaining" is based on
    (the layer's expected count from the density skeleton and the bright
    fraction, then updated as it goes?); hysteresis so the second bar
    doesn't flash on and off near the 30-second line (the old per-sector
    bar was removed for exactly that); and whether the same rule applies
    to sector fill (`Sectors (ring … layer …)`), where PERF.3's rate is
    measured.

  - [ ] **PERF.5 Scatter bright stars in stages, one luminosity band at a time**
    Boss (2026-10-01): "let's do a default of 100 solar
    luminosities for the star map, but then add to the TODO.md to allow
    us to add another layer down (i.e. so when I generate I do say 500
    solar luminosities because I want to be quick and do testing but then
    after I want to generate down to 100 solar luminosities, so we have
    to make sure when I do that, it only generates between the limits
    (i.e. doesn't generate more brighter stars). Probably add a value for
    the star-fill level." Boss then kept the default at 500 (2026-10-01):
    "OMG, no, so let's make the default 500 then, sorry, I am now down
    with adding 35 gigs to the database." So the default threshold
    (`program_constants.BRIGHT_STAR_MIN_LUMINOSITY_SOL`) is 500, and going
    down to 100 later is the kind of extra layer this item adds. Today
    the plan's scatter (`generate.py`, `--bright-star-min-luminosity`)
    clears `bright_stars` and redraws everything at or above the
    threshold, and refuses when any sector is already filled unless
    `--force` leaves those sectors out. Done: the galaxy stores its
    current star-fill level (the lowest luminosity already scattered); a
    scatter to a lower threshold keeps the existing bright stars and adds
    only stars from the new threshold up to (not including) the stored
    level, then lowers the stored level; asking for a level at or above
    the stored one does nothing (or says so); the Generate page and the
    CLI show the current level and offer "go down to N"; and PERF.3's
    size and time estimates cover just the new band. Open questions:
    where the level is stored (the control database, a `galaxy` row next
    to the skeleton, or derived from `MIN(luminosity_w)` in
    `bright_stars`)? What happens to sectors already filled when new,
    dimmer bright stars land in them: add the stars and build their
    systems in place, skip those sectors, or mark them for ADM.8's
    regenerate? Must the new band draw from the same random stream, so a
    500-then-100 galaxy matches a straight-to-100 one?

  - [ ] **PERF.9 Weight the bright-star ETA by the shape of the galaxy**
    Boss
    (2026-10-01): "so that the bright stars ETA takes into account the
    shape of what's being generated (ie that at layer 0 and 635 take very
    little time)". Today the plan's scatter (`brightStars.scatter`) walks
    the layers in `extents` order and calls `on_layer(done, total)` once
    per layer, so the "Bright stars (layers)" bar in `generate.py` counts
    every layer the same. The thin layers at the edges finish almost at
    once and the dense middle layers take most of the time, so the ETA
    swings badly. Done: the bar and its ETA count expected work, not
    layers. Before the scatter starts, each layer gets an expected star
    count from the same density model the scatter uses
    (`_ring_bins` × `expected_at_density_1` × the bright fraction for the
    chosen threshold), and progress and the ETA are measured in expected
    stars done out of the expected total, so a run through the sparse
    edge layers no longer makes the rest look quick or slow. Ties in with
    PERF.3 (the up-front time estimate uses the same per-layer
    weights), UX.3 and PERF.4 (the banner's ETA, and PERF.4's
    stars-remaining estimate for the current layer), PERF.5 (a staged
    scatter weights only the new luminosity band), and PERF.7 and
    PERF.8 (the decaying-average rate and the parallel tasks). Open
    questions: is the weight the expected star count alone, or does it
    also count rings and slots walked (an empty edge layer still costs
    some loop time)? How does the weighting combine with PERF.7's
    decaying average: the average measured in expected stars per second,
    or in layers per second and then scaled? Is the per-layer expected
    count worked out in a quick pre-pass at the start of every plan, or
    stored with the galaxy skeleton? Once PERF.8 runs layers in
    parallel and out of order, does the ETA add up the expected work
    still queued rather than following the layer order?

  - [ ] **PERF.10 Record generation speed across a log scale of densities**
    Boss (2026-10-01): "record generation stats such as time per star,
    time per sector for a log scale of densities from 0.01 to the max
    expected density / actual density found. ... Both of these will be
    continued to be refined and calculated as long as the galaxy is in
    existence but as a decaying average." Today nothing records how long
    generation takes per star or per sector, and PERF.3's planned
    stars-per-second figure is a single number for the whole server.
    Done: density is split into log-scale buckets from 0.01 up to the
    highest density expected or found; every sector fill adds its time
    per star and time per sector to its density's bucket as a decaying
    average; the buckets keep updating for as long as the galaxy exists;
    and they can be read back by the tools below. These stats feed
    PERF.3 (time estimates before bulk generation, per bucket instead
    of one rate), UX.3 (the banner's ETA), PERF.4 (the slow-layer
    bar's stars-remaining ETA), PERF.5 (the time for a new luminosity
    band), PERF.7 (the decaying-average rate) and PERF.9 (weighting
    the bright-star ETA by expected work per layer). Open questions: how
    many buckets and where their edges sit (per decade, half-decade?);
    whether the top edge is fixed from the density model's expected
    maximum or moves up when a denser sector is found; the decay
    constant (how fast old runs fade); whether the plan's bright-star
    scatter gets its own buckets (its cost per star differs from sector
    fill); whether it lives in the control database (survives a new
    galaxy, as PERF.3 suggests for stars per second) or the galaxy
    database (resets with it), and so what a regenerate or reset does to
    it; and whether PERF.8's parallel workers count wall time or CPU
    time per task.

  - [ ] **PERF.11 Store each sector's expected and actual density**
    Boss
    (2026-10-01): "add stats for each sector's density expected and
    actual in the database in a way that can be easily accessed. Both of
    these will be continued to be refined and calculated as long as the
    galaxy is in existence but as a decaying average." Today `sectors`
    stores only a sector's address and center; its expected density is
    worked out on demand from the galaxy skeleton (`relative_density`
    times `galaxy_layer`'s `expected_system_count_at_density_1`), and
    its actual density means counting its systems. Done: every sector
    has its expected density and its actual density (systems, and stars,
    found when filled) stored where a query can read them directly, as
    columns on `sectors` or a sector-stats table, readable by
    `queryDb`, `adminStats` and the API; the galaxy-wide comparison of
    expected against actual is kept as a decaying average and updated
    after every fill. PERF.10 places each sector in its density bucket
    with these numbers, and PERF.3, PERF.5 and PERF.9 use the
    expected-versus-actual ratio to correct their estimates. Open
    questions: columns on `sectors` (a migration in the Database
    workstream) or a separate table; whether "actual" counts systems,
    stars, or both; what the decaying average is taken over (the ratio
    per density bucket, so it ties in with PERF.10, or one galaxy-wide
    figure); whether existing sectors are backfilled by a migration;
    and what happens to the stats when a sector is regenerated
    (ADM.8) or the galaxy is reset.

## ADM: Admin tools

- [ ] **ADM.1 Admin editing: overrides, delete and regenerate**
  Boss asked for these on 2026-10-01 (quoted where it matters). None is
  designed yet; the open questions are listed in each subitem.

  - [ ] **ADM.6 Admin override of a planet's or moon's class**
    Boss: "it should
    have the option to do 'recommended' which are other classes that fit
    within the given space or I can 'force' which means it sets it to what
    I want no matter what. During validation orbital paths in the way of a
    class change get recalculated and moved around until the system is
    stable. This should be recursive, so that say a moon is changed, then
    that lunar system is changed, then it goes out from there to recheck
    all the planets, etc. It does this until the system can validate."
    "The system will always try to have a stable system and will warn the
    user if that isn't possible." Done: an admin control on a planet or
    moon offering a "recommended" list (classes that fit its current
    space) and a "force" choice; after a change, revalidate outward
    (moon, its lunar system, then every planet of the star system) with
    ADM.5's functions, re-spacing orbits until it validates, and warn
    the admin when no stable layout exists. Open questions: does a forced
    change that can't be made stable still save (with the warning), or is
    it rolled back? May revalidation remove other bodies, or only move
    them?

  - [ ] **ADM.7 Admin override of a star**
    Boss: "that will change the entire
    system but it will change the system to have as many objects as the
    original system had just their orbital positions will change, caveat
    there is if there are too many objects for the star (say a large star
    with a lot of objects changes to a small star that does not have
    orbital space, then it'll be truncated." Done: an admin control to
    change a system's star; the system keeps its planets, moons and belts
    (same count), orbits are re-spaced for the new star with ADM.5's
    validation, and outer objects are dropped when the new star lacks the
    room, telling the admin what was removed. Open questions: do the
    planets keep their classes, or are classes re-checked against the new
    star's zones (which could chain into ADM.6's revalidation)? Does
    this cover companion stars in multiple systems too?

  - [ ] **ADM.8 Delete and regenerate buttons on everything, sector down**
    Boss: "I also want a delete function across the board, so I can
    manually remove a system. With that a regen button ... regenerate a
    system, phenomena, planet, asteroid belt, sector, basically anything
    from a sector to anything in a sector should have a delete and regen
    buttons when admin is logged in." `DELETE /api/sectors/<id>` and
    `DELETE /api/systems/<id>` exist (`src/html/api/routes.py`); there is
    nothing for a single phenomenon, planet, moon or belt, and no
    regenerate for any of them. Done: admin-only Delete and Regenerate
    buttons on the sector, system, phenomenon, planet and asteroid belt
    pages (Boss's list; moons are an open question), with a confirm step,
    an audit-log entry, and ADM.5's validation after a single body is
    removed or regenerated. Open questions: does regenerating a sector
    keep manual overrides (ADM.6 and ADM.7), renamed objects and
    placed facilities (DB.1), or replace everything? Does regenerating
    keep the object's name? Does deleting a sector leave its slot
    unfilled (so it can be filled again) or mark it empty?

- [ ] **ADM.4 Collapsible Generate page sections; pick the center sector**
  Boss (2026-10-01): "In generation screen each section should be
  collapsible and generate around a sector should have the option to
  locate an existing filled sector or put in the coordinates." Today
  the admin Generate page (`web/templates/generate.html`) shows every
  section (Current job, One-off system, New galaxy, Generate sectors,
  Plan the galaxy, Rebuild the bright stars, Reset) open, one after
  another, and "around a sector" (`mode == "center"`) asks for a
  numeric sector ID and a radius. Done: each section can be collapsed
  and expanded (a `<details>` or a heading button, keyboard and screen
  reader friendly); "around a sector" lets the admin either find an
  existing filled sector (search by name or designation, or pick it on
  the Galaxy Map or from a list) or type coordinates (a ring, layer and
  slot address, or galaxy-frame x, y, z). Open questions: which
  sections start open (only Current job, or the last one used,
  remembered per browser)? Which coordinates: a sector address, a
  position in pc or ly, or both? Does "locate" reuse the Sector Map pick
  mode (MAP.21) or the address bar's `/galaxy/locate` (MAP.24)?

## SEC: Security

No open items. The login protection of 2026-10-01 (SEC.1, SEC.20 to
SEC.28: the always-on activity log, per-address and per-username
lockouts, trusted devices, the password blocklist and hashing cost,
two-factor sign-in and the fail2ban example) shipped in PRs #217, #220
and #221.

Design: [docs/design/login-brute-force-protection.md](design/login-brute-force-protection.md)

## TEST: The test suite

From the test suite plan of 2026-10-01 (Boss: "build me a list of TEST
items (TEST.x) to build our testing suite to cover all edge cases, and
prepare to do another solid full on bug hunt. Don't implement just plan
for it"). Each item ends with its areas in brackets. "Suspected bug"
means read from the code, not reproduced yet; the bug hunt confirms or
clears each one.

### Infrastructure and CI

- [ ] **TEST.1 Test category and suite markers**
  Add TEST to `test_todo_tags.py`; register `db`, `slow`, `browser`
  markers in `pytest.ini` so fast/no-DB, DB and browser runs can be
  picked separately. [infra]

- [ ] **TEST.2 Parallel test runs**
  Add pytest-xdist, give each worker its own control database name
  (today `configured_control_database()` defaults to one fixed name),
  and a session-scoped template schema so each test doesn't rebuild
  `schema.sql` from scratch; goal: the 17-22 minute suite well under 10.
  [infra]

- [ ] **TEST.3 MariaDB in CI**
  Add MariaDB 10.11 and 11.x legs (and MySQL 8.4) to `ci.yml`; today CI
  is MySQL 8.0 only, so the engine-specific paths in `_db.py` (statement
  timeout via `max_statement_time`, the ALGORITHM fallbacks, the
  CTE-in-UPDATE workaround) only run on Boss's server. [infra, DB]

- [ ] **TEST.4 Revive and widen the known-bug tests**
  The "Real bugs (strict xfail)" block in
  `test_fuzz_system_generation.py` (about line 562) now passes on 5-200
  seeds; widen the seed counts and fix the stale docstring; promote the
  three Tier-2 "tracked-not-fixed" reports (NaN radius in the system
  map, `SERIALIZABLE_FIELDS` drift, galaxy-boundary report) to hard
  asserts or strict xfails. [infra, GEN]

- [ ] **TEST.5 Real 4 pc in boundary tests**
  `test_bughunt_galaxy_boundaries.py` uses `EDGE_PC = 10.0`; run it at
  the real 4 pc sector edge too. [infra, GEN]

### Database and migrations

- [ ] **TEST.6 SQL portability lint**
  Parse `schema.sql`, `control_schema.sql` and the SQL strings in
  `_db.py`/`queryDb.py` against the reserved-word lists of MySQL 8.0,
  8.4 and MariaDB 10.x/11.x; flag deprecated `VALUES()` in ON DUPLICATE
  KEY. [DB]

- [ ] **TEST.7 Strict sql_mode on both engines**
  Run the DB tests with `ONLY_FULL_GROUP_BY` and `STRICT_TRANS_TABLES`
  forced on (MySQL 8's default, not MariaDB's), so a GROUP BY or
  truncation that MariaDB forgives fails locally too. [DB]

- [ ] **TEST.8 Migrate from real old schemas**
  Checked-in historic schema dumps (for example v8, v20, v33, v44) each
  migrated to the latest and compared with a fresh `schema.sql` database
  (tables, columns, types, indexes, FKs). Steps with no test today:
  14-15, 15-16, 18-19, 22-23, 23-24, 24-25, 25-26, 30-31, 43-44, 44-45.
  [DB]

- [ ] **TEST.9 Migration crash and re-run**
  Fail at step N (DDL has already auto-committed), re-run, and get a
  correct database; every step idempotent when applied twice; empty or
  gapped `schema_migrations`. [DB]

- [ ] **TEST.10 Database newer than the code**
  A galaxy or control database with a version above the code's is
  refused with a clear message (suspected bug: `migrateDb.py` says "no
  migration path available yet" and the app uses it anyway). [DB]

- [ ] **TEST.11 Every column round-trips**
  Walk `INFORMATION_SCHEMA.COLUMNS` and prove each column is written by
  an insert and read back by a loader, so a new column left NULL or
  never loaded fails. [DB]

- [ ] **TEST.12 Boundary values round-trip**
  Float extremes, NaN/inf refused, DECIMAL precision, VARCHAR length
  limits under strict mode, 4-byte UTF-8 names, NULL tristates. [DB]

- [ ] **TEST.13 Collation collisions**
  Names equal under `utf8mb4_unicode_ci` but not byte-equal (case,
  accents) against the unique name registries; a database created with
  the server's default collation joined to pinned tables (error 1267).
  [DB]

- [ ] **TEST.14 CHECK constraints enforced**
  Each CHECK in the schema rejects a bad row on both engines. [DB]

- [ ] **TEST.15 Sector save fails halfway**
  Inject an error after the systems are written and before phenomena or
  neighbours; nothing persists, no orphan rows, no stale name
  reservations or `GET_LOCK`; retry exhaustion at
  `SECTOR_SAVE_ATTEMPTS`; lock-wait timeout (1205). [DB]

- [ ] **TEST.16 Id blocks after reset and rollback**
  Id block cache across `resetDb`, a manual insert, two processes
  exhausting blocks; no duplicate primary key. [DB]

- [ ] **TEST.17 Batched writes at the limits**
  Multi-row inserts near `max_allowed_packet`; FK ordering for
  self-referencing tables. [DB]

- [ ] **TEST.18 Full-text search edge cases**
  Names with `+ - " * '`, words shorter than `innodb_ft_min_token_size`,
  stopwords that differ between engines, `%`, `_` and backslash in
  `/search` fields; results checked, not just "no 500". [DB, UX]

### Generation and the work queue

- [ ] **TEST.19 Same galaxy at any worker count**
  One seed generates identical sectors with `--workers 1`, 2 and N
  (suspected bug: the one-worker path in `workQueue.submit` never calls
  `random.seed(task_seed(...))`; the parallel path does, and the only
  test compares 2 with 3). [GEN, PERF]

- [ ] **TEST.20 Work queue failure paths**
  A worker dies (BrokenProcessPool), `on_done` raises, a payload won't
  pickle, a result isn't JSON, two tasks share a key, heartbeat fails,
  the control database drops mid-run, lease expiry under clock skew.
  [PERF]

- [ ] **TEST.21 Cancelling a run**
  SIGTERM, Ctrl+C and SystemExit during a parallel run end the job as
  cancelled, free the lease and leave no half-written sector. [PERF]

- [ ] **TEST.22 Every bulk mode in parallel**
  `--shell`, `--block`, `--column`, `--center-sector`, random start and
  `sector --num-sectors N` with `--workers 2`, checking run counts and
  that no sector is filled twice. [PERF, GEN]

- [ ] **TEST.23 Resume after an interrupted fill**
  Stop a ring, shell or block run partway, run it again, and get the
  same result as one uninterrupted run with no duplicates. [GEN]

- [ ] **TEST.24 Bright-star scatter edge cases**
  `_place_one` running out of redraws, a zero-weight bin picked by float
  rounding, a layer where nothing qualifies, `outer_ring=0`, empty
  extents, a threshold below every white dwarf. [GEN, PERF]

- [ ] **TEST.25 Interrupted bright-star scatter**
  A worker fails mid-scatter after some 10,000-row commits (the seed is
  written only at the end); a re-plan or later fill handles the partial
  table. [GEN, PERF]

- [ ] **TEST.26 `--force` scatter then fill**
  Sectors skipped by a forced scatter fill correctly afterwards. [GEN]

- [ ] **TEST.27 Progress and ETA under bad clocks**
  `DecayingRate` with time going backwards, NaN or infinite amounts,
  many adds at one instant, a tiny rate; workers never write the
  progress file. [PERF]

- [ ] **TEST.28 CLI errors by message**
  Every `parser.error` in `generate.py` (block, column, shell,
  center-sector, limit, plan shape, workers, population, phenomenon,
  mysql-port) asserted by its text, and each limit tested at exactly its
  maximum (`MAX_GENERATE_RING`, `MAX_GENERATE_LIMIT`, first and last
  `--block-layer`). [GEN, OPS]

- [ ] **TEST.29 Limits stay consistent**
  `MAX_GENERATE_LIMIT` matches `ring_sector_count(MAX_GENERATE_RING)`
  whatever `DEFAULT_MAX_RING` is. [GEN]

- [ ] **TEST.30 Grid seams and the nucleus**
  Points at θ just under 2π and at -0.0 on the 4 pc grid; the outermost
  planned ring and layer against `galaxy_bounds`; the nucleus sector
  (ring 0, slot 0) across layers 0 and -1. [GEN]

- [ ] **TEST.31 Sector placement exhaustion**
  "could not place a new object", the Poisson cap, explicit positions on
  the cell boundary, `nearest_neighbors` with a bad count. [GEN]

- [ ] **TEST.32 System builder internals**
  Direct tests for `generate_slot_object`,
  `calculate_distance_for_class`, `_forced_habitable_distance`,
  `_trim_to_orbit_ceiling`, `_reconcile_moved_planet`,
  `_clear_circumbinary_floor` and the `from_dict` error; these are what
  ADM.5's validate module (`stellarObjects/validation.py`) wraps. [GEN]

- [ ] **TEST.33 Moon stability helpers**
  `moon_orbit_bounds_km`, `drop_unstable_moons`, `update_hill_sphere`
  and the "no valid planet class" errors, tested directly. [GEN]

- [ ] **TEST.34 Kepler solver extremes**
  Eccentricity above 0.99, negative mean anomaly and above 2π,
  non-convergence detected rather than silently returned,
  `_real_cube_root` at 0 and negative. [GEN]

- [ ] **TEST.35 Star and evolution helpers**
  `calculate_heliosphere`, the population-model star path, age windows,
  radius and temperature helpers, white dwarf radius, Yerkes class, and
  their raise messages. [GEN]

- [ ] **TEST.36 Phenomenon class helpers**
  Direct tests for nebula, remnant, black hole, rogue planet, comet and
  asteroid-field class and designation helpers, and every
  `get_table_properties`. [GEN]

- [ ] **TEST.37 Names under parallel saves**
  Two workers saving systems and sectors with the same base name at
  once; the diminutive tier filling under parallel saves; the species
  name race and "could not find a free species name". [GEN, DB]

- [ ] **TEST.38 Population incremental rescans**
  `scan_life_worlds`, `refresh_civilizations`, `refresh_territories` and
  the watermark path; a rescan after new fills adds only the new worlds.
  [POP]

- [ ] **TEST.39 Navigation graph**
  K-d tree neighbours checked against brute force, duplicate
  coordinates, k = 0 and k >= n, NaN positions, travel time table. [NAV]

### Web, API and jobs

- [ ] **TEST.40 Two admins start a job at once**
  Exactly one job runs (suspected bug: `_take_lock` creates an empty
  lock before writing the job id, so a second caller can read it as
  stale, delete it and start its own job). [ADM]

- [ ] **TEST.41 Job files damaged**
  Corrupt or truncated `job.json`, `state.json`, `progress.json`; a lock
  holding garbage; job id collision; unwritable jobs directory; prune
  never removes the running job; cancel with unknown, malformed or
  finished job ids. [ADM]

- [ ] **TEST.42 Pages fresh after a CLI write**
  After `generate.py` writes straight to the database, `/galaxy`,
  `/sector`, tiles and lists show the new data (only API writes are
  tested today). [UX, PERF]

- [ ] **TEST.43 Auth sweep over every route**
  Generated from `app.url_map`: every API write route gives 401 to
  anonymous, garbage Bearer and revoked keys; every admin route gives
  403 to an admin who must still change credentials. [SEC, API]

- [ ] **TEST.44 What an API key may do**
  Whether a Bearer key can change credentials, set up or turn off TOTP,
  make keys, or log out, pinned to the intended answer. [SEC, API]

- [ ] **TEST.45 More than one admin**
  Admin B can't revoke admin A's key, lifting another admin's lockout is
  audited, two admins editing the same system. [SEC, ADM]

- [ ] **TEST.46 Trusted device and TOTP edge cases**
  Expired, tampered and other-user device cookies; turning TOTP off
  voids trust; a code reused across the API and `/login/code`; a pending
  login that expires. [SEC]

- [ ] **TEST.47 Oversized requests**
  Multi-megabyte JSON and form bodies to `/api/systems`, `/login` and
  the facility form get 413 (there is no `MAX_CONTENT_LENGTH` set
  today). [SEC]

- [ ] **TEST.48 Security headers everywhere**
  CSP and the other headers on JSON responses, 404/405/500 pages and
  redirects, not only pages. [SEC]

- [ ] **TEST.49 Thin API routes**
  Unknown ids, empty galaxy, paging limits and wrong-system ids for
  `/api/galaxy/sectors`, `/shape`, `/phenomena`, `/bright-stars`,
  star/planet/moon PATCH, facilities POST/PATCH/DELETE,
  `/api/admin/login-failures`, `/api/population`, `/api/species/<id>`,
  `/api/systems/<id>/owner`; deleting a sector that has facilities or
  wiki pages. [API]

- [ ] **TEST.50 Galaxy URLs combined**
  `?at=` with `?p=` and `?sector=` together; `?course=` to deleted
  objects; `/galaxy/locate` with unicode, very long input, NaN/inf and
  out-of-range coordinates, ambiguous names. [MAP]

- [ ] **TEST.51 Page-number sweep gaps**
  `/species?species_page=` and `/polities?polities_page=` join the
  page-clamping sweep. [UX]

- [ ] **TEST.52 Old URLs and error codes**
  Unknown `/<name>.py`, case variants, redirect chains; 400 for
  malformed form encoding; HEAD and OPTIONS on pages. [UX]

- [ ] **TEST.53 Formatters with bad numbers**
  Every `fmt` and `tabledisplay` formatter with NaN, inf, negative, zero
  and None; empty tables; huge values. [UX]

- [ ] **TEST.54 Caches under threads**
  Page cache fill and clear from real threads; two writers to the same
  tile file. [PERF]

### Browser and JavaScript

- [ ] **TEST.55 Map buttons do something**
  Playwright clicks every map control (Sector Map, System Map,
  phenomenon diagram, Galaxy Map) and asserts the view changes. Found
  while planning: on the phenomenon diagram, any nebula or remnant about
  half a light-year across or larger opens already at the 1 ly zoom-out
  limit (`phenomenonmap.py` lines 53 and 127), so "-" does nothing; this
  test pins the UX cleanup item. [MAP, UX]

- [ ] **TEST.56 No overlapping controls**
  Playwright compares the bounding boxes of every button and control on
  every page at 390, 600, 820 and 1280 px, both themes; no two
  intersect, none off-screen. [UX]

- [ ] **TEST.57 Galaxy Map JavaScript logic**
  Node tests for `galaxystageview.js` (zoom, pan and tilt clamps) and
  `galaxymap3d.js` (`sectorDesignation` BigInt packing, address form,
  history, control handlers); neither has any test today. [MAP]

- [ ] **TEST.58 Other map JavaScript**
  Node tests for `mapzoom.js` (`zoomedBox` and its clamp),
  `sectormap.js` and `systemmap.js` zoom and selection logic,
  `generatejobs.js` polling, `facilityform.js`. [MAP]

- [ ] **TEST.59 Galaxy Map drill-down in a browser**
  Playwright walks quarter, layer, arc, block and sector by clicks,
  checks the URL and breadcrumb at each step, Back/Forward, and the free
  camera from an arc down. [MAP]

### Scripts and ops

- [ ] **TEST.60 Admin script command lines**
  `main()` tests for `queryDb`, `adminStats`, `checkRenderParity`,
  `dedupeNames`; bad port, unknown database, empty password vs
  environment for every script; `resetDb --yes --dry-run`;
  `updateOrbits` with the clock moved back; `loginLockouts` with bad
  IPv6. [OPS]

- [ ] **TEST.61 SQLite import script**
  `migrateSqliteToMysql.py` has no tests and (suspected bug) can't run:
  it requires the SQLite file to be at today's schema version (46),
  which no SQLite database ever was. Test it from a v12 file, or retire
  the script. [OPS, DB]

- [ ] **TEST.62 update.sh against a real database**
  A CI job runs `update.sh` (and `install.sh`'s database step) on Linux
  against a live database: up to date, needs migrating, newer than the
  code, unreachable, failed migration. Today they're only
  syntax-checked. [OPS]

- [ ] **TEST.63 Math check that runs first**
  Boss (2026-10-01 15:29Z): "I want a specific way that runs first
  before other tests that basically validates the math works, before
  batch generation, we need to verify the actual math works." One
  module, `src/stellarObjects/mathCheck.py`: pure functions, no
  database, network or files, fixed seeds, under 5 seconds. Each check
  has a name, the function it calls, the expected value, a tolerance and
  the source of the expected value (a textbook figure, a paper's table,
  or an exact identity). The same module is used three ways: pytest runs
  it first, `generate.py` and the Generate page run it before bulk
  generation, and `update.sh` runs it after an update. The coverage
  tests in TEST.4 and TEST.32 to TEST.36 stay as they are; this item is
  the reference-value gate in front of them, and TEST.34's Kepler
  reference values move here.

  Default taken: no skip switch for the bulk gate, since it costs under
  5 seconds. Open questions for Boss: should there be an emergency skip
  flag anyway? Should the web app also run it at startup and show admins
  a warning if it fails? Should the one-off system generator run it too,
  or only bulk paths?

  - [ ] **TEST.64 Reference values**
    Known answers from real astronomy, each within a stated tolerance.
    For example: the Sun (1 M_sun gives 1 L_sun, about 10 Gy on the main
    sequence, about 5,772 K from L and R through Stefan-Boltzmann);
    Earth's orbit (1 AU around 1 M_sun is 1 year at 29.78 km/s by
    vis-viva; Jupiter 11.86 years); Earth's Hill sphere about 1.5
    million km; habitable zone and snow line at 1 L_sun; a 0.6 M_sun
    white dwarf about Earth-sized; the Sun's Schwarzschild radius 2.95
    km; the Sun's galactic orbit (about 8 kpc, about 220-230 km/s, about
    230 My); Holman-Wiegert critical radii from the paper's table; the
    Kepler and Barker equations against known solutions.

  - [ ] **TEST.65 Identities and invariants**
    Things that must be exactly or nearly true for any input. Every unit
    conversion round-trips (pc, ly, AU, km, mpc) and the constants agree
    with each other (found while planning: `SPEED_OF_LIGHT_M_S` is
    2.998e8 while `LIGHTYEAR_M` uses the exact 299,792,458 m/s, a 0.003%
    mismatch); luminosity rises and lifetime falls with mass; orbital
    energy is conserved around a Kepler orbit; the sector grid's cell
    volumes add up to each ring's annulus,
    `sector_address_at(sector_position_pc(...))` returns the same
    address, and ring sector counts match `ring_sector_count`; density
    is normalised to 1 where the code says it is; no NaN or infinity
    over a fixed sweep of inputs.

  - [ ] **TEST.66 Distributions match their targets**
    With fixed seeds, a few thousand draws of the IMF, star ages, the
    Poisson sector counts, the bounded bell and the planet class table
    land on their intended shares within a statistical tolerance (for
    example a chi-square test), so a broken sampler fails even when
    every single value looks fine.

  - [ ] **TEST.67 Runs first in the suite and in CI**
    A `mathcheck` marker (TEST.1 adds the markers); `conftest.py` moves
    those tests to the front and stops the run if any fails, saying the
    math is broken and the rest would be noise; CI runs it as its own
    quick first job that the other jobs wait on.

  - [ ] **TEST.68 Gate before bulk generation**
    `generate.py check-math` runs it by hand; every bulk path (`galaxy`,
    `sector` over many sectors, `plan`, `population`, the Generate
    page's jobs and the map's block and neighbourhood generation) runs
    it first and refuses to start if a check fails, naming the failed
    check and writing nothing; the Generate page shows the result as the
    job's first step; `update.sh` runs it after updating and warns on
    failure.

## USR: User accounts

- [ ] **USR.1 User accounts**
  Boss (2026-10-01): "a full user level interface to allow users to
  bookmark this will, of course, require an email loop for password
  setting / resetting, invite only, so only an admin can invite a user
  which is done by unique link". Today there are only admin accounts:
  `admin_users`, `admin_sessions`, `admin_api_keys` and
  `admin_audit_log` in `stellarObjects/control_schema.sql`, managed by
  `stellarObjects/adminAuth.py`. Install seeds one admin with a random
  first password and `must_change_credentials`
  (`bootstrap_control_schema`), and nothing in the web interface or the
  CLI adds another account. There is no email support and no
  saved-bookmark feature (pages only have bookmarkable URLs). Order:
  USR.2, then USR.3 and USR.4, then USR.5, USR.6 and USR.7.
  Login protection already exists per username
  (`src/html/api/loginbackoff.py`, PR #131) and is planned per IP
  address (SEC.1); both must cover user logins, password resets and
  invite links too, as must SEC.20's logging and SEC.26's second
  factor.

  - [ ] **USR.2 Accounts with roles: user, admin and Owner**
    Boss: "Admin can
    then make a user admin or take away admin rights on everything but the
    primary 1st admin account generated at install, that'll have a
    designation of 'Owner' and no other account can override that." Done:
    one accounts table in the control database (with an email address)
    holding a role per account (user, admin, Owner); exactly one Owner,
    which is the account install creates; an admin page where admins
    promote a user to admin or demote an admin to user, refused for the
    Owner; sessions and API keys that work for every role while the
    admin-only pages and API routes stay admin-only; each role change
    written to the audit log. Open questions: does `admin_users` become
    this table (renamed, with a role column) or do users get their own
    table? Which existing account becomes Owner on a server that already
    has several admins (the lowest id?)? Can an admin demote themselves,
    or another admin? Can an admin delete or disable a user account? Do
    users get API keys?

  - [ ] **USR.3 SMTP settings in the admin config**
    Boss: "We'll use SMTP for
    email which means admin config needs SMTP settings." Done: SMTP host,
    port, security (STARTTLS or TLS), username, password and From address,
    set from an admin page (and the config file/installer), with a "send
    test email" button; one small mail module that every flow in USR.4
    to USR.6 uses, which logs failures and never shows the SMTP
    password. Open questions: is the SMTP password kept in the config
    file (like the database password) or in the control database, and is
    it encrypted there? Who can change SMTP settings: any admin, or only
    the Owner? What do invites and resets do when SMTP isn't configured
    (show the link to the admin to pass on by hand?)?

  - [ ] **USR.4 Invite-only sign-up by unique link**
    Boss: "only an admin can
    invite a user which is done by unique link, admin can select how many
    uses the link has or if it expires in 1 hour, 4, 6, 12, 24, 3 days, 7
    days, 1 month, 1 year, or never. Likewise invite # of uses the link
    gets can be infinite but infinite and never expire in combo ask for
    confirmation as this is not recommended." Done: an admin page that
    creates invite links (a random token stored hashed, like session
    tokens) with a use count (a number, or infinite) and one of the
    expiry choices above; infinite uses together with never expires asks
    the admin to confirm and says it isn't recommended; a list of invites
    with uses left, expiry and who made them, and a way to revoke one;
    opening a valid link lets someone register (username, email), then
    USR.5's email loop sets their password; the account is a user, not
    an admin. Open questions: does an admin optionally type the invitee's
    email so the link is sent for them, or only copy the link? Is the
    invite page rate-limited, and is there a cap on open invites? Does
    a multi-use link record who used it?

  - [ ] **USR.5 Email loop for setting and resetting passwords**
    Boss: "an
    email loop for password setting / resetting". Done: a new account
    sets its first password from an emailed link; "forgot password" on
    the login page emails a reset link; links are single-use, short-lived
    (stored hashed) and end the account's other sessions once used; the
    page never says whether an email address has an account; resets are
    rate-limited per address and per IP alongside the existing login
    backoff and SEC.1. Open questions: how long a reset link lasts
    (30 minutes? 1 hour?); does changing the email address also need an
    email confirmation to the old and new addresses; does the Owner's
    reset need anything extra?

  - [ ] **USR.6 Owner transfer**
    Boss: "The owner CAN (with specific approval
    and double password confirmation and email loop confirmation) assign
    owner to someone else who then has to accept via email loop." Done:
    only the Owner can start a transfer, from a page that asks them to
    confirm the choice explicitly, enter their password twice, and then
    confirm from an email link; the chosen account then gets an email and
    must accept from its own link; only when both are done does Owner
    move; every step goes to the audit log. Open questions: what the old
    Owner becomes (admin?); does the new Owner have to be an admin
    already; how long the pending transfer lasts and whether the Owner
    can cancel it; what happens if the Owner loses their email or
    password (recovery from the server's command line?).

  - [ ] **USR.7 A user-level interface with bookmarks**
    Boss: "a full user
    level interface to allow users to bookmark". Done: signed-in users
    (any role) get an account page and can bookmark sectors, systems,
    planets, phenomena and NAV courses, see them in a list, name them and
    remove them; bookmarks are stored per account in the control
    database. Open questions: which objects can be bookmarked, and can
    users also add notes? What else a user can do that an anonymous
    visitor can't (is the site still public to read, or sign-in only?)?
    Do bookmarks survive a galaxy regenerate (object ids change), and if
    not, what does a broken bookmark show?

## OPS: Installers, hosting, CI, releases

No open items; OPS.1 shipped with the version scheme in `changes/README.md`.

## VIEW: The view from a planet

- [ ] **VIEW.1 View from a planet**
  **Research first.** Boss: "view-from-planet will have to do
  calculations on colors and A LOT Of stuff, so make special note of
  that, it will need a full research pass." Before any code for VIEW.2
  and VIEW.3, Boss wants a research session with him "into exactly how
  one would do that". VIEW.2 and VIEW.3 are blocked on it; VIEW.4
  is not.

  - [ ] **VIEW.2 A starmap seen from a planet. RESEARCH WITH BOSS FIRST**
    Boss:
    "Build a function to select a planet and generate an effective starmap
    from that planet based on all visible stars, this will have to include
    a lot, A LOT, of math so remind me to do research when we get there
    into exactly how one would do that, but it would have to account for
    where each star would have been at that light years back in time?"
    Done: pick a planet (or moon), and get every star visible from it
    with its direction and brightness as seen there, placed where it was
    when the light now arriving left it (light-travel time back along
    its galactic orbit; the correlative update now moves everything along
    galactic orbits, GEN.6, PR #157). The research pass covers at least:
    - which stars are visible (apparent magnitude from luminosity and
      distance, a magnitude cut, interstellar extinction and reddening by
      dust, and whether stars beyond the generated sectors are included,
      for example from the density skeleton or the `bright_stars` table
      from PR #159);
    - star colors as seen from the planet (colour from temperature,
      reddening, the planet's atmosphere and its own star's glare; Boss
      singled out colors as needing real work);
    - light-time positions (the star's position at "now minus distance /
      c" along its galactic orbit), and where the planet is in its own
      orbit and its sky orientation (axial tilt, rotation, latitude);
    - nebulae, the galactic band, the companion stars of the planet's own
      system, and performance (millions of stars per view).

  - [ ] **VIEW.3 Render the view as a PNG, with constellations**
    Boss: "when
    it does that it will generate a PNG and it will generate
    constellations." Done: VIEW.2's view is drawn to a PNG (a sky
    projection, star size and colour by apparent brightness), and the
    brighter stars are grouped into constellations with lines and names
    from VIEW.4, stored so a planet keeps the same constellations each
    time. Blocked on VIEW.2's research pass. Open questions: whole-sky
    or a horizon view from a point on the surface? How are constellations
    chosen (bright-star patterns, by clustering, a set number per sky)?
    Are the PNGs cached on disk and served by the web interface, or made
    on request?

  - [ ] **VIEW.4 Constellation names in the name generator**
    Boss: "add to our
    name generator constellation name support based on constellation
    names throughout all known languages and then slice it up like we do
    for all our naming". Done: a constellation name list in
    `stellarObjects/names.py` gathered from constellation and star-group
    names across the world's languages and sky cultures (not just the 88
    IAU ones), and a constellation name generator that slices and
    recombines them into new names the same way stars, planets and
    sectors are named (`split_into_syllables` in `utils.py`, the
    prefix/suffix lists, the `offensive_words.txt` filter). Used by
    VIEW.3. Open questions: what counts as a source list (licensing of
    sky culture data such as Stellarium's), transliteration of non-Latin
    scripts, and whether names are unique per planet or galaxy-wide.
