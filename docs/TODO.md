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
4. **Admin editing:** ADM.1 and its subitems, starting with the validate
   module (ADM.5).

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

  - [ ] **ADM.5 Central validate module in `stellarObjects`**
    Boss: "we should
    get a whole set of validate functions in their own file within the
    stellarObjects class (if we don't already) so that we can have a
    central place to validate a system (lunar system, star system) which
    includes validating planets." There isn't one today; validation is
    spread out: `SystemData.validate_system`,
    `_validate_cross_star_clearance` and `_trim_to_orbit_ceiling`
    (`systemData.py`) for orbits, the `_validate_*` checks in
    `planetPhysics.py` for a planet's class/radius/mass, and
    `moon_orbit_bounds_km`/`drop_unstable_moons` (`planetPhysics.py`)
    for moons. Done: one module (for example
    `stellarObjects/validation.py`) that validates a planet, a lunar
    system and a star system, which generation and ADM.6 and ADM.7
    both call, with existing behavior unchanged. Prerequisite for
    ADM.6 and ADM.7.

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
