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
- The phase plans in [plan/](plan/) are the other place that lists
  item IDs: each phase file has a table of its items (see "Plan:
  phases" below). File a new item into one phase's table in the same
  PR, and delete its row when it ships.
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

SEC, DOC and POP have no open items today.

## Background

The long-term goal is a fully populated galaxy (every sector, system,
planet, moon and belt) stored in MySQL and browsed through a web
interface served from Apache. The original roadmap phases are all done:
object-graph serialization, the relational schema and migrations, CLI
tools writing to the database, lazy galaxy-scale generation from a
density skeleton, and the Flask API (`src/html/api/`) with server-rendered
pages served by the same app (`src/html/web/`): Galaxy, Sector and System
maps, search, NAV, admin auth and wiki publishing.

## Plan: phases

Boss (2026-10-01 23:53Z) asked for the open work to be planned in
phases, one plan file per phase with everything that phase needs, and
this file as the master and index. On 2026-10-02 the phases were
rebuilt as phases 0 to 3+ from the dependency report (each item's
prerequisites, the files it shares, and Boss's decisions of that
night). Every open item below is in exactly one phase, and every
item's prerequisites are in its own phase or an earlier one. Each phase
file gives the goal, the build threads with their order and
prerequisites, and the open questions;
[plan/notes.md](plan/notes.md) keeps the research notes, the files
several items share and the judgment calls. This file keeps each
item's full text. A new item goes into a phase's table in the same PR
that files it.

| Phase | Plan | Goal | Items |
|---|---|---|---|
| 0 | [phase-0-roots.md](plan/phase-0-roots.md) | Fix the bugs nothing else depends on and lay the groundwork everything later builds on: the parallel path first (Boss: top priority), the galaxy seed, the physics and database bugs, the map groundwork through the arc pick (MAP.85, which Boss wants in phase 0), object references and the small page and ops fixes. | TEST.76, TEST.74, PERF.22, PERF.21, TEST.73, PERF.23, GEN.32, GEN.39, DB.3, DB.2, GEN.34, GEN.35, GEN.36, GEN.25, GEN.37, GEN.45, GEN.53, GEN.54, GEN.49, GEN.50, GEN.31, MAP.90, GEN.46, TEST.70, MAP.87, MAP.83, MAP.82, MAP.84, NAV.30, MAP.81, MAP.63, MAP.64, MAP.57, MAP.88, NAV.7, OPS.6, OPS.7, OPS.8, TEST.71, TEST.72, UX.34, OPS.9, UX.28, UX.23, UX.2, UX.24, UX.29, ADM.14, MAP.60, MAP.55, MAP.85, MAP.52, API.15, TEST.79, NAV.34, DB.6, OPS.10, TEST.78, NAV.38 |
| 1 | [phase-1-built-on-roots.md](plan/phase-1-built-on-roots.md) | Work that needs phase 0 in place: galaxy generation on the parallel path with the per-sector stats table (GEN.44 and PERF.11), prevalence controls, planet classes, reproducible galaxies up to the golden-seed test, sector colors, routing, the first picker pieces and the queue. | GEN.33, GEN.28, GEN.38, GEN.27, GEN.51, GEN.52, TEST.75, ADM.16, GEN.48, GEN.24, GEN.44, PERF.11, PERF.1, GEN.41, GEN.47, MAP.89, MAP.86, MAP.80, NAV.8, NAV.9, NAV.13, NAV.14, NAV.10, NAV.12, UX.35, NAV.11, UX.22, UX.25, UX.26, UX.31, UX.27, UX.3, PERF.19, ADM.15, API.4, API.7, API.9, OPS.11, GEN.56, GEN.57, DB.7, GEN.58, TEST.77 |
| 2 | [phase-2-maps-picker-backfill.md](plan/phase-2-maps-picker-backfill.md) | The Galaxy Map built out around the arc pick, the shared picker and courses, the parallel backfill and density pass, and the API pieces remote generation needs first. | MAP.56, MAP.53, MAP.58, MAP.78, MAP.76, MAP.54, MAP.75, MAP.59, MAP.77, MAP.65, MAP.79, NAV.15, NAV.29, NAV.33, NAV.31, NAV.16, NAV.20, NAV.21, NAV.17, NAV.18, NAV.4, NAV.24, GEN.29, UX.32, UX.30, UX.33, GEN.42, GEN.43, PERF.18, GEN.40, PERF.20, API.5, API.10, API.11, API.12, MAP.69, MAP.70, API.16, ADM.17, NAV.36, NAV.39 |
| 3 | [phase-3-engine-3d-remote.md](plan/phase-3-engine-3d-remote.md) | The three maps on one engine, the 3D system view, courses that bend around gravity wells, and remote generation through the API, reproducing what the server would make. | MAP.66, MAP.67, MAP.68, MAP.61, MAP.71, MAP.72, MAP.73, MAP.74, MAP.62, NAV.32, NAV.3, NAV.22, NAV.23, NAV.5, NAV.25, NAV.26, NAV.27, NAV.28, NAV.6, API.13, API.14, API.8, ADM.13, API.3, UX.21, API.17, GEN.59 |
| 3+ | [phase-3plus-accounts-sky-galaxies.md](plan/phase-3plus-accounts-sky-galaxies.md) | The open-ended tail: user accounts (with API.6 keys and saved courses), the view of the sky from a planet, the plan for more galaxies, and the end state of reproducible galaxies (`generate.py reproduce`). | API.6, USR.2, USR.3, USR.4, USR.5, USR.6, USR.7, USR.1, NAV.19, VIEW.1, VIEW.4, VIEW.2, VIEW.3, GEN.9, GEN.55, OPS.12 |

Phases overlap: a phase's later threads can start while the next
phase's first ones run, as long as the order inside each phase holds.

Boss's list of 2026-10-01 23:53Z (`new todos.txt`, with research notes;
the files are in the project's shared files under `todo-tasks/research/`)
became these items:

| # | Boss's item | ID |
|---|---|---|
| 1 | Store every sector's backfill level (-1, lowest L_sun, 0) | GEN.44 |
| 2 | Rogue planets after systems and phenomena; expanded rows span the table | UX.24 |
| 3 | Rogue planet gas giant vs terrestrial probability | GEN.45 |
| 4 | Rogue planet octant and a map symbol link | UX.25 |
| 5 | Edit and admin buttons as a button menu | UX.26 |
| 6 | Rogue planets dim, barely noticeable | MAP.82 |
| 7 | "Mark Rogue Planets" keeps its highlight while on | MAP.83 |
| 8 | Rogue planets smaller unmarked, bigger and clickable marked | MAP.84 |
| 9 | System and nav buttons in one row, nav folding into a menu | UX.27 |
| 10 | Icons instead of words on buttons (investigate) | UX.28 |
| 11 | Star names at most two words | GEN.46 |
| 12 | Every comet's type shown with a link | UX.29 |
| 13 | Planet information without the Markdown render | UX.30 |
| 14 | System edit as a quick menu, not a long panel | UX.31 |
| 15 | Planet rows show class only | UX.32 |
| 16 | Filter phenomena by type and class | UX.33 |
| 17 | No nebulae being created | GEN.47 |
| 18 | 3D galaxy, arc pick, no sector lines (other galaxy-map items edited to match) | MAP.85 |
| 19 | Filled sectors translucent, colored by their stars; blocks averaged | MAP.86 |

Done before this plan: the bug round (MAP.17 to MAP.19, MAP.26, MAP.37,
MAP.43 to MAP.51, UX.15, UX.16, UX.19, UX.20, ADM.9), login security
(SEC.1, SEC.20 to SEC.28), the database-call work (PERF.6, PERF.12 to
PERF.17), parallel generation (PERF.7, PERF.8), the round Boss set at
14:41Z (PERF.5, UX.13, UX.14, MAP.15, MAP.30, MAP.2, GEN.23, PERF.3,
PERF.4, PERF.9, PERF.10, ADM.5 to ADM.8) and the OPS, ADM and TEST
build of the evening.

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

- [ ] **UX.2 Menus sized to what they hold (bug)**
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

- [ ] **UX.21 Clean up the web interface: overlapping buttons and dead controls (bug)**
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
  1280 px) shipped in PR #307: the nebula "-" no-op is pinned by a
  strict xfail in `test_web_browser_maps.py` (fix it and drop the
  xfail), and the overlap check found no overlapping controls on the
  pages as they were then.
  Order (pre-planning thread): run it after MAP.55 and MAP.60, which
  already remove some dead controls on the Galaxy Map.

- [ ] **UX.22 Meaningful units for every measurement**
  Boss (2026-10-01 15:10Z): "standardize ALL measurements into trees
  like we have so that we always have meaningful units. From mass, to
  distance, to time, to speed, just everything that can have units.
  Atmospheric pressure and surface conditions should show customary
  units as well as a secondary to help contextualize the metric values
  given." Today only distance has a ladder (`format_distance_m` and
  friends in `stellarObjects/utils.py` and `html/lib/fmt.py`, mirrored
  by `static/distance.js`), speed (UX.13: `format_speed_kms`,
  `static/speed.js`) and durations (UX.14: `format_duration_seconds`,
  `format_period_years`, `static/period.js`). Done: one ladder per quantity, in Python with a JavaScript
  mirror, picking a meaningful unit the same way, and every page, map
  panel and form converted to it: mass (kg, Earth, Jupiter and solar
  masses), distance, time, speed, temperature, pressure, gravity,
  density, luminosity, power and any other quantity shown with a unit.
  Surface conditions show temperature in K, °C and °F; atmospheric
  pressure and the other surface conditions show a customary unit
  (such as atm, psi or g) beside the metric value. Surface
  temperature (K, °C, °F) and pressure (kPa, atm, psi) shipped in PR
  #234; this item covers the rest. Open questions: the ladder and switch
  points for each quantity; which customary unit goes with each surface
  condition; whether the secondary unit shows in tables or only in
  detail panels.

  - [ ] **UX.23 A shared unit-ladder module**
    Done: one Python ladder module (in `stellarObjects/utils.py` or a new
    `units.py`) and its JavaScript twin, generalising the existing
    distance, speed and duration ladders, with a test that the Python and
    JavaScript twins agree for a table of values. Then one sub-item per
    quantity family: mass; temperature; pressure and gravity; density,
    luminosity and power.

- [ ] **UX.24 Sector contents: rogue planets after systems and phenomena, and expanded rows the full table width (bug)**
  Boss (2026-10-01 23:53Z): "rogue planets should come after star
  systems and other phenomena, and each line item should span the table
  when expanded." Today `_contents` in `web/sector_page.py` sorts every
  row (systems, phenomena, facilities, the folded "N rogue planets" row)
  by distance from the sector's center only, and the only row that
  expands is the rogue group, a `<details>` inside its Name cell
  (`templates/sector.html`), so its members squeeze into one column.
  Done: the table lists star systems first, then the other phenomena,
  then rogue planets (each group still nearest first); an expanded row's
  content spans every column of the table (a full-width detail row
  under it, `colspan` of the whole table), on phones too.

- [ ] **UX.25 Rogue planets: octant and a small map symbol beside each name (bug)**
  Boss (2026-10-01 23:53Z): "Rogue planets each one has a location and
  octant, show on map should be a little map symbol next to the name as
  a link." Today a lone rogue planet's row shows its octant and a
  "Show on map" text button, but the members of the folded rogue group
  show no octant or location, only a "Show on map" button
  (`sector_page.py` `_rogue_group_row`, `sector.html`). Done: every
  rogue planet, alone or in the group, shows its octant and location,
  and "Show on map" is a small map icon link beside its name (with a
  text label for screen readers and a tooltip) that selects it on the
  Sector Map, as the button does today. Applies to the same link on
  other rows once UX.28 settles the icon set.

- [ ] **UX.26 Edit and admin actions as a button that opens a menu (bug)**
  Boss (2026-10-01 23:53Z): "Edit buttons and admin buttons in any view
  (such as sectors or star systems) should appear as a button menu
  (click the button the menu appears)." Today the sector page has an
  Admin panel (generate-neighborhood forms) and an Edit panel, and the
  `edit_buttons` macro (`templates/partials/edit_controls.html`) shows
  Regenerate and Delete as inline `<details>` confirm forms. Done: in
  every view, the admin and edit actions sit behind one button that
  opens a menu (keyboard and screen-reader friendly, closes on Escape
  or a click outside, doesn't shift the page); picking an action opens
  its confirm step or form on top (a popover or dialog), not inline
  below the page. UX.31 is the system page's case.

  - [ ] **UX.31 Editing a star system: an edit button with a quick menu, not a long panel (bug)**
    Boss (2026-10-01 23:53Z): "Edit on a star system shouldn't be a big
    long menu under the interface but instead an edit icon or button
    that when clicks opens a little quick menu showing the edit
    options." Today the system page has a separate Edit panel below the
    System panel (`system.html`, `_edit_rows` in `system_pages.py`): a
    table of every body with Regenerate, Delete, "Change star" and
    "Change class" `<details>` forms. Done: the system, and each body
    in the System panel, has an edit icon button that opens a small
    menu of its actions (Regenerate, Delete, Change star or Change
    class); choosing one opens just that form; the long Edit panel is
    gone.

- [ ] **UX.27 System page: the system and navigation buttons on one row that doesn't overlap (bug)**
  Boss (2026-10-01 23:53Z): "When viewing a star from sector view system
  and nav buttons should be in a row and should not overlap. If room is
  needed the nav buttons can collapse into a nav button that opens a
  from here or to here menu." Today the system page's subhead
  (`system.html`, `.page-subhead` and `.page-actions` in `style.css`)
  wraps the badges, "Navigate from here", "Navigate to here", "Show on
  Galaxy Map" and the bookmark button. Done: these buttons sit on one
  row with no overlap at every size class; when there isn't room, the
  two navigate buttons fold into one "Navigate" button whose menu holds
  "From here" and "To here" (container query, not a device check).
  Ties in with NAV.29 (the Start Here / End Here wording) and UX.21.

- [ ] **UX.28 Investigate icons instead of words on buttons**
  Boss (2026-10-01 23:53Z): "Investigate using symbols for buttons
  instead of words where appropriate. Adhering to modern web interface
  standards." Done: a short survey of every button and link on the site
  saying which should become an icon (for example edit, delete, show on
  map, filter, navigate, bookmark, menu), the icon set to use (inline
  SVG, one sprite, no icon font or external CDN), and the rules
  (visible tooltip, `aria-label`, 44 px touch target on coarse
  pointers, a text label kept where an icon alone is unclear); Boss
  approves the list, then the changes are filed as their own items or
  folded into UX.25, UX.26, UX.27 and MAP.55.

- [ ] **UX.29 Every comet in a system shows its type as a link (bug)**
  Boss (2026-10-01 23:53Z): "Comets in a star system some show the type
  with a link and some don't. Type / link should be there for all."
  Today `_comet_row_html` (`lib/systempage.py`) shows "Elliptical comet"
  or "Parabolic comet" as plain text and links only the period class
  (`class_url("comet", period_class)`); parabolic comets have no period
  class (`cometData.py`), so they get no link, and an unknown code gets
  none either. Done: every comet row shows its type as a link to its
  class reference page, parabolic comets included (a class entry for
  them, or a link to the comet class overview), and a test checks every
  comet kind has a working link.

- [ ] **UX.30 Planet information without the Markdown render**
  Boss (2026-10-01 23:53Z): "Rework planet information displays to pull
  away from the markdown render and instead fits with the modern web UI
  as we've seen it so far." Today each body's expanded row on the
  system page is the generator's Markdown (`systemRender.render_system_sections`)
  turned into HTML by `lib/mdconvert.py` (`systempage._row_html`), and
  so is the system overview; only the stat chips are built as HTML.
  Done: the system page builds each body's details from its stored
  values as structured HTML (property grids, atmosphere composition
  bars, orbit figures, moons, life), styled like the rest of the site in
  both themes and every size class, using the unit ladders (UX.22);
  `mdconvert` is no longer used for them. The wikitext and Markdown
  views (the Wikitext/Markdown toggle and wiki upload) keep using the
  text render. Open question: does the overview text stay as prose?

- [ ] **UX.32 Planet rows show the class only, without the type and moon labels**
  Boss (2026-10-01 23:53Z): "Remove type of planet (Terrestrial / gas)
  and moon from the row display just giving class instead of class and
  type." Today each planet row's chips (`systempage.py`, around line
  215) are "Class X", a type chip ("Terrestrial" or "Gas Giant" from
  `_type_chip`, or "Habitable"), optional "Habitable moon" and
  "Inhabited", then distance, period and gravity, and moons sit in a
  "N moons of X" group under the row. Done: the row shows the class
  (linked to its reference page) and no Terrestrial/Gas Giant chip or
  moon label; the moons stay reachable from a small count on the row.
  Open question: does "Habitable" stay as a chip?

- [ ] **UX.33 Filter phenomena by their classes and types (bug)**
  Boss (2026-10-01 23:53Z): "Searching for phenomena should have the
  ability to filter rogue plants by type, comets by types, really and
  class related to the phenomena should be searchable." Today the
  `/phenomena` list (`system_pages.py`, `phenomena.html`) has no filters,
  and the Search page's "Phenomenon Class" facet covers only nebulae,
  supernova remnants and asteroid fields (`_SEARCH_PHENOMENON_CLASS_COLUMNS`
  in `queryDb.py`), not rogue planets (`rogue_planets.planet_class`) or
  interstellar comets. Done: phenomena can be filtered by kind and by
  every class or type that kind has (rogue planet class and terrestrial
  or gas giant, interstellar comet type, nebula class, remnant type,
  black hole and neutron star kinds), on both the Search page and the
  `/phenomena` list, through the API, with the counts per option.

- [ ] **UX.34 The sector summary calls white dwarfs "B-type" and "A-type" systems (bug)**
  Low. Found by the debug-mode bug hunt (2026-10-02; report and evidence in the project's shared files under `bug-hunt/`). The "Systems: 4 B-type, ..." line in `generate.py`
  (around line 1449) groups by spectral letter only; in a 40-sector run,
  all 4 "B-type" and 6 of 9 "A-type" systems were white dwarfs. Done:
  white dwarfs (and giants) are counted under their own label.

- [ ] **UX.35 The NAV page route shown horizontally, wrapping onto several lines on narrow screens**
  Boss (2026-10-02 01:53Z): "display the path horizontally and find a
  way to split it between multiple lines for mobile or limited
  displays." Today the route is a vertical list (`<ol class="nav-route">`
  in `nav.html`). Done: the stops run left to right, each a link, with
  the hop distance between them; on phones and narrow panels it wraps
  onto several lines (a container query, not a device check), never
  splitting a stop across lines; screen readers still get an ordered
  list. Runs alongside NAV.12; NAV.36 styles its unknown-space hops.


## MAP: Galaxy Map, Sector Map, System Map

MAP.2 with MAP.22 and MAP.23, MAP.15 and MAP.30 shipped in PR #234.

- [ ] **MAP.52 Galaxy Map highlights the wrong area; pick a 40-degree wedge around the cursor (bug)**
  Boss (2026-10-01 19:35Z): "galaxy map still doesn't highlight
  correctly. highlights should be tween valid beginning and end points,
  adjust the quadrant philosophy to selecting a wedge of the galaxy in
  40 degree arcs and not set, the cursor point will be the center of
  the arc, so it will always be +/- 20 degrees from the cursor's
  merdidian snapped to the available meridians from the center that go
  from the center to the edge." Today the first pick on the Galaxy Map
  is a fixed quarter of the disk (`galaxystages.js`, `kind:
  "quadrant"`, four arcs, MAP.19), later picks are fixed regions, and
  the hover highlight (`galaxystageview.js`) can start or end where no
  wedge line is. Done: the first pick is a 40-degree wedge centered on
  the cursor's angle from the galaxy's center, running from the center
  to the edge, its two sides snapped to the nearest meridians (the
  wedge lines that run from the center to the edge, MAP.42 to MAP.44);
  the highlight always starts and ends on such valid lines and follows
  the cursor as it moves; clicking zooms to that wedge (the wedge zoom
  of PR #243) with no gaps between blocks (PR #224); the URL and the
  breadcrumb label name the wedge by its angles rather than "Quarter
  n". Boss (19:37Z): "Doesn't have to be +/- 20 so long as it fits
  into the wedge from center (ring 1) to edge." So 40 degrees is the
  target, not an exact width: the wedge snaps to lines that run all the
  way from ring 1 to the edge, and may come out a little wider or
  narrower. Open questions: how snapping works where meridians stop
  short of ring 1 (the inner rings have fewer slots); do the later
  picks (regions inside the wedge) follow the same cursor-centered
  rule?
  Order (pre-planning thread): These nine all change
  `galaxystageview.js` and `galaxystages.js`, so they suit one build
  thread in this order: MAP.60 (scale line) and MAP.55 (buttons into a
  menu) first (small, and they free space); MAP.52 (the 40-degree
  wedge); MAP.56 (drop the 3x3 pick); MAP.53 (rotate and fit); MAP.58
  (zoom limits); MAP.54 (slab buttons and leader lines); MAP.59 (ghost,
  mini map, header). MAP.57 (System Map NaN) is independent. MAP.61's
  first two sub-items (shared helpers, one camera controller) should
  come before or with MAP.53 and MAP.58, which add camera rules.
  Arc pick (MAP.85, Boss 2026-10-01 23:53Z): the first pick is an arc,
  not a whole wedge from the center to the edge: the 40-degree width and
  meridian snapping here still apply to its bearing, and it is also
  bounded in radius; the highlight shows the arc and its neighbors'
  boundaries only.

- [ ] **MAP.53 Rotate a zoomed-in wedge, and zoom it to fit the window (bug)**
  Boss (2026-10-01 19:40Z): "allow the user to rotate the galaxy wedge
  around once it's zoomed into the wedge and zoom in more based on
  window size." Today the wedge zoom of PR #243 fits the real wedge at
  a fixed orientation and scale. Done: once zoomed into a wedge, the
  user can rotate the view around it (drag, keys, and a touch gesture,
  like the free camera below quarter level), and the zoom fits the
  wedge to the map's actual size, so a bigger window shows it larger;
  it refits when the window is resized or rotated. Ties in with MAP.52
  (the 40-degree wedge pick). Boss (19:44Z): "Slabs will rotate
  around their immediate center", so the view turns about the middle
  of what is shown (the wedge, or the slab), not the galaxy's center.
  Open question: is the rotation kept in the URL and bookmarks?
  Arc pick (MAP.85, Boss 2026-10-01 23:53Z): the view rotates about the
  picked arc (then the slab), fitted to the window.

- [ ] **MAP.54 Slab leader lines instead of the slab slider (bug)**
  Boss (2026-10-01 19:40Z): "don't use a slider for the slab, instead
  have a line going from each slab on the map (dynamically rendered to
  always point where it needs to) from the button for that slab to the
  slab itself on the map." Today slabs are picked with the slab slider
  beside the map (MAP.30, shipped as a slider in PR #234,
  `galaxystageview.js`). Done: the slider is replaced by one button per
  slab, and each button has a line drawn from it to its slab on the
  map; the lines are redrawn whenever the view rotates, zooms, pans or
  the window resizes, so they always point at the slab; hovering or
  focusing a button highlights its line and slab, and clicking picks
  the slab as the slider does today. MAP.30 stays done; this item
  replaces its slider. Open questions: how the lines stay readable with
  many slabs (thin lines, only the hovered one drawn bright, or
  grouping); how crossing or overlapping lines are kept apart; and
  where the buttons sit on a phone-width screen. Default taken: a slab
  should always be on screen after the zoom-fit, but if one ever falls
  outside the frame (the isometric tilt, a very tall stack of slabs, a
  small window), its line ends at the window's edge with an arrow
  pointing toward it and its button still works; this may never happen
  in practice. A slab hidden behind another keeps its line, drawn to
  the visible part.
  Arc pick (MAP.85, Boss 2026-10-01 23:53Z): the slabs are height bands
  of the picked arc.

  - [ ] **MAP.76 Leader-line layout**
    An SVG overlay above the canvas, recomputed on every camera change
    from the slabs' projected centres; lines kept from crossing by
    ordering the buttons by the slabs' screen height; a phone layout
    with the buttons in one column below the map. Picks MAP.54's
    defaults for its open questions.

- [ ] **MAP.55 Galaxy Map buttons: a menu, with only back, forward, up, reset and bookmark showing**
  Boss (2026-10-01 19:44Z): "the button row in the galaxy view should
  be under a menu except for back, forward, up, and reset. Reset and
  whole galaxy do the same thing. Remove the wedges button entirely,
  and bookmarks should be next to the back, forward, up, reset, and
  bookmarks. Sector sell should go next to the galaxy map if there's
  room right next to the slab buttons." Today the Galaxy Map's controls
  (`galaxystageview.js`, `galaxymap3d.js`) are one row of buttons.
  Done: only back, forward, up, reset and the bookmark button stay in
  view, in that row; every other control moves into one menu button
  beside them (keyboard and screen-reader friendly); "Whole galaxy" is
  removed, since reset does the same; the Wedges button is removed
  entirely; the "Sector cell" info panel (`#galaxymap3d-info`, showing
  a picked sector's address and designation) sits beside the map next
  to the slab buttons (MAP.54) when there is room, and below the map
  when there isn't. Ties in with UX.21 (overlapping buttons). Open
  question: what is in the menu and in what order?
  Arc pick (MAP.85, Boss 2026-10-01 23:53Z): the wedge lines go entirely
  (no sector lines on the galaxy), as well as the Wedges button.

- [ ] **MAP.56 Drop the 3x3 block pick: select a slab, zoom in, select a segment (bug)**
  Boss (2026-10-01 19:47Z): "I was wrong when before I said a 3x3 cube
  of blocks should be selectable by the user, it isn't, so lets take
  that back to the select-slab, zoom in, select segment of slab."
  Today the drill-down (`galaxystages.js`) alternates layer and region
  picks (MAP.19, PR #208: "quadrant, layer, region, layer, region, ...,
  layer, sector"), and a region pick offers up to 3 x 3 options (a
  third of the rings across, a third of the arc along, `PICK_SPLIT`).
  Done: the region (3 x 3) pick is removed; the ladder becomes the
  wedge pick (MAP.52), then select a slab (the slab buttons and lines
  of MAP.54), then the view zooms to that slab (fitted to the window
  and rotatable, MAP.53), then select a segment of the slab, repeating
  slab and segment inside each smaller block down to a sector. Default
  taken: a segment is one drill block of the next level inside the
  slab (27 or 3 sectors a side), picked directly on the zoomed slab
  with the same hover highlight as today's blocks; no existing item
  defines it further. The URL and breadcrumb forms of a region pick
  ("r4") go away; old links with one open at the nearest valid stage.
  MAP.19's big targets still apply. Open question: should a segment be
  one block, or a run of blocks along the arc when a block is too small
  to click on a small screen?
  Arc pick (MAP.85, Boss 2026-10-01 23:53Z): the ladder is arc, slab,
  segment, then slab and segment again down to a sector.

- [ ] **MAP.57 The System Map writes NaN or infinite positions into its SVG (bug)**
  Found by the generation tests (2026-10-01): a body whose computed
  position is NaN or infinite is written straight into the System
  Map's SVG. Done: such a body is left out or drawn at a safe place
  with a note, the SVG never holds NaN or inf, and the strict xfail
  test for it passes.

- [ ] **MAP.88 Parts of a star system run off the edge of the System Map (bug)**
  Boss (2026-10-02 00:49Z): "sometimes the star systems go off the edge
  for planets and such, the window should be scaled down so that we
  don't run any part of the star system past the outer edge of the
  viewing space." Today the System Map (`lib/systemmap.py`) draws into a
  fixed 700 px square viewBox: each body's distance from its anchor is
  log-scaled out to `_MIN_RADIUS_PX + _RADIUS_SPREAD_PX` (335 px from the
  center), then moons are placed around their planets, belts get a band,
  and `_relax_markers` pushes crowded markers apart, all of which can
  carry an outer planet, its moons or its marker past the frame; only
  the labels are slid back inside it (MAP.50). Done: after everything is
  placed, the map works out the drawn extent of every star, planet,
  moon, belt, facility and marker (with its radius) and scales the
  whole scene down to fit the frame with a small margin, so nothing
  ever runs past the edge, with a test over many generated systems that
  every drawn element sits inside the viewBox. Separate from the orbit
  spacing study (the "Orbit spacing options" thread), which may change
  how the same map spaces orbits; whichever lands second keeps this
  fit.

- [ ] **MAP.89 System Map: space orbits with a fitted scale and a minimum ring gap instead of plain log**
  Boss (2026-10-02 00:45Z): "investigate different ways to space orbits
  visually, right now we use logarithmic spacing because of the vast
  distances involved, but investigate other methods of spacing to make
  better use of the visual space and make more visual sense to someone
  looking at it. Then pick the best option and add the plan to the
  TODO.md file." Today `lib/systemmap.py` places every orbit on one log
  scale per scene (`_radial_scale_bounds`, `_radial_px`: `55 + frac *
  280` px in a 700 px square), used for planets, belt edges
  (`_belt_band`), close-binary stars, facilities and moon scenes, and
  `static/systemmap.js` `kmToPx()` repeats the formula for the Measure
  distance route. In a study of 240 generated systems, 13 methods were
  scored (the "Orbit spacing options" thread; plan and comparison
  renders in the project's shared files under `orbit-spacing/`). Plain
  log leaves the tightest ring gap at a median 2.4 px, with 27% of gaps
  under 6 px. Square root, cube root and linear knot the inner planets
  together. Even (rank) spacing throws distance away. Frost-line zones
  leave half the map empty. Chosen: a fitted scale with a gap floor,
  which gives a median tightest gap of 11.3 px with 15% under 6 px, and
  keeps big real gaps looking big better than log does.
  1. Fit the scale per scene: power curves from linear (p = 1) down to
     log (p = 0) between today's `lo`/`hi` bounds; pick the most linear
     one that needs the gap floor on no more than a quarter of its ring
     gaps (compact systems go linear, huge-span systems stay log).
  2. Gap floor: any gap between neighboring rings (orbits, belt edges,
     close-binary star rings) under 12 px widens to 12 px, and the other
     gaps shrink in proportion so everything fits; with too many orbits
     for 12 px each, the floor drops to 60% of an even share.
  3. One mapping for everything: a list of knots (km, px), interpolated
     in log km between knots, written onto the scene (for example
     `data-knots`) so `kmToPx()` reads it instead of its own log formula.
  Order, real angles and "farther out is drawn farther out" don't
  change. Done:
  - `_radial_scale_bounds`/`_radial_px` are replaced by the fitted scale
    with the gap floor, in the system scene, both wide-binary scenes and
    moon scenes.
  - Belts, facilities and close-binary star placement use the same
    mapping; no orbit is drawn inside a belt ring (MAP.49's test still
    passes).
  - `systemmap.js` `kmToPx()` reads the scene's knots; the measure route
    still meets marker centers and bends around bodies as today.
  - The panel hint stops saying "log-scaled distance" (for example
    "distance scaled to fit, order and angle true").
  - Tests: the scale only increases with distance; no neighboring-ring
    gap under 12 px when the orbits allow it; a compact system gets
    p = 1; a huge-span system falls back to log; JavaScript and Python
    give the same pixel radius for the same knots.
  - Before and after screenshots of the six system types in the study.
  Separate from MAP.88 (fit the whole drawn system inside the frame);
  both change `lib/systemmap.py`, and whichever lands second keeps the
  other working. Decided (Boss, 2026-10-02 01:01Z: "I agree, we'll go
  with fitted for the orbital spacing in MAP.89"): fitted scale with the
  12 px minimum ring gap, and no "fitted / even" spacing toggle.

- [ ] **MAP.58 Galaxy Map zoom limits: a short manual range on the galaxy wedge, locked below it**
  Boss (2026-10-01 20:45Z): "It may be necessary for users to zoom in
  and out manually. This should only be within a short range ... they
  can zoom in to about, say, twice as close as it starts out and they
  can zoom out back to the full galaxy, but no further. When we get into
  smaller chunks like blocks and wedges that aren't as large as the full
  galactic wedge, we're going to lock the zoom ... The system will still
  be able to zoom in stages as we've discussed but the user won't be
  able to arbitrarily zoom in and out." Today every view below the whole
  galaxy and its quarters can be zoomed freely with the wheel, a pinch
  or the zoom keys (`isFree` in `galaxystageview.js`), from
  `MIN_ZOOM` = 1/8 of the stage's fitted camera distance (8 times
  closer) to `MAX_ZOOM` = 2.5 times it, and the galaxy and its quarters
  can't be zoomed at all. Done:
  - On the full galaxy wedge (the 40-degree wedge of MAP.52, fitted to
    the window by MAP.53), the user can zoom in to about twice as close
    as the fitted view, and out only until the whole galaxy fits, with
    the wheel, pinch, keys and any zoom buttons.
  - On every smaller view (slabs, segments, blocks, the sector cube),
    user zoom is locked: the wheel scrolls the page, and pinch and the
    zoom keys do nothing. The staged zoom of each pick (MAP.53, MAP.56)
    still animates to its fitted view.
  - Rotation (MAP.53) and panning are not affected.
  - Reset (MAP.55) returns to the fitted zoom.
  Ties in with MAP.53, MAP.55 and MAP.56. Open questions: does the
  whole-galaxy view itself zoom (today it doesn't)? Should a locked view
  keep panning, or only rotate?
  Arc pick (MAP.85, Boss 2026-10-01 23:53Z): the short manual zoom range
  applies on the whole galaxy and the picked arc; it locks below the
  arc.

- [ ] **MAP.59 Make it plain that a zoomed-in slab is a slab, not a wedge**
  Boss (2026-10-01 20:45Z): "We need to make it clearer, when we've
  zoomed into a specific slab, that we're viewing a specific slab and
  not a wedge. I'm not sure how to do that so do some research on that
  and then add the to-do items to make it happen."
  Arc pick (MAP.85, Boss 2026-10-01 23:53Z): the ghost is the rest of
  the picked arc.

  Why it looks like a wedge today: once a slab is picked, only that
  slab's blocks are drawn (`galaxystageview.js`). A slab is a thin
  layer of the wedge (a sector layer inside a level-3 block, otherwise
  9 or 27 sector layers, `slabLayers` in `galaxystages.js`), so seen
  from the isometric tilt it looks like a flat wedge. Nothing on screen
  shows the rest of the stack, and the slab's height range appears only
  as "Slab 3" in the breadcrumb.

  Options (from how 3D map, CAD and volume viewers show a selected
  slice):
  1. **Ghost of the parent wedge.** The other slabs of the wedge stay in
     view as a faint, see-through outline (wireframe edges only), and
     the picked slab is the one solid layer inside it. This is the cut-away
     or "section view" of CAD tools and of floor pickers in building
     maps. It shows at a glance that this is one layer of a taller stack.
  2. **Slab thickness edges.** Draw the slab's top and bottom faces and
     its vertical side edges in a distinct line color, so it reads as a
     slice with depth rather than a surface.
  3. **A labelled header over the map**, such as "Slab 3 of 9 · 120 to
     160 pc above the plane (layers 25 to 33)", replacing the bare
     "Slab 3". The breadcrumb crumb says "Slab 3 of 9" too.
  4. **Side-view inset.** A small fixed diagram in a corner of the map
     shows the wedge edge-on as a stack of bars, with the picked slab
     highlighted and the galactic plane marked. This is the slice
     indicator of medical and volume viewers. It doubles as a slab
     picker if clicks on it are allowed.
  5. **Tint.** The picked slab gets a color band that differs from a
     whole wedge, kept the same at every depth.

  Boss's answers (2026-10-01 20:50Z): "Ghost of the other wedges should
  be just wire lines and faint and yes I want to implement a mini map
  that shows the segment of the whole galaxy. When we zoom in to a slab
  off to the side, have a locked view in isometric form of the block,
  highlighting which slab we're in. Then add navigation tools so that
  if the user goes and clicks on another slab in the isometric view, it
  switches to the slab in the main view", and for the height: "both".
  So the build is options 1, 2, 3 and 4.

  Done:
  - **Ghost:** with a slab picked, the rest of the wedge or block is
    drawn as faint wire lines only, with no fill, around the solid,
    edge-lined picked slab. The ghost is not clickable and doesn't block
    clicks on the slab or its segments.
  - **Mini map:** beside the main view sits a small, locked isometric
    view of the whole block (or wedge) the slab belongs to, showing
    every slab of it with the current one highlighted, and where that
    block sits in the whole galaxy. It doesn't rotate or zoom. Clicking
    (or tapping, or picking with the keyboard) another slab in it
    switches the main view to that slab, the same way picking a slab
    does today, with the URL and breadcrumb following.
  - **Header and breadcrumb:** "Slab N of M" with its height both ways:
    the distance above or below the galactic plane in the map's chosen
    units (pc or ly, the units setting of PR #234) and its sector layer
    numbers, e.g. "Slab 3 of 9 · 120 to 160 pc above the plane · layers
    25 to 33". The same applies one level down ("Layer N of M" inside a
    block).
  - Everything stays aligned when the main view rotates or zooms.

  Ties in with MAP.53 (rotation keeps the ghost aligned), MAP.54 (the
  picked slab's button and line stay highlighted; the mini map is a
  second way to pick a slab), MAP.55 (where the mini map sits next to
  the slab buttons and the Sector cell panel, and below the map on a
  phone), MAP.56 (the segment pick happens on the solid slab) and
  MAP.58 (the mini map is never zoomable).

  - [ ] **MAP.75 The mini map as a second engine view**
    The mini map is a second, locked camera on the same scene data;
    built on MAP.61's controller with a "locked" policy rather than its
    own renderer (one WebGL context, scissored like the System Map's
    sphere overlay).

- [ ] **MAP.60 Galaxy Map scale readout: one scale line**
  Boss (2026-10-01 20:50Z): "I want to trim the scale information from
  the galactic map so that it just has one scale line." Today the
  readout under the Galaxy Map (`updateScaleBar` in `galaxymap3d.js`,
  `#galaxymap3d-scale`) stacks three lines: "1 px ≈" (what one screen
  pixel spans), "1 block =" (the size of one drawn block) and a scale
  bar of about 70 px with its length. Done: only the scale bar and its
  length remain, on one line, in the map's chosen units (sectors and pc
  or ly, as now).

- [ ] **MAP.61 One map engine and control set for the Galaxy Map and the Sector Map**
  Boss (2026-10-01 20:55Z): "unify the sector view with the galactic
  view so it's all the same code and control set, because right now
  they're different." Today the Galaxy Map (`galaxymap3d.js`,
  `galaxystageview.js`, `galaxystages.js`) and the Sector Map
  (`sectormap.js`) are separate three.js pages, with their own camera,
  controls, picking, tooltips and scale readouts. Done: one shared
  engine (scene, camera, controls, picking, hover, info panel, scale
  line, bookmarks, keys and touch) draws both. The sector is the
  deepest stage of the galaxy drill-down, with the same buttons and
  gestures, and the only differences are the data each level shows.
  Ties in with NAV.3 (the shared picker), MAP.53 to MAP.60 (the
  Galaxy Map controls being reworked now) and MAP.62 (the system view
  joins the same engine).
  Order: the first two sub-items change no behaviour and should land
  before (or as the first PR of) the MAP.52 to MAP.60 work, because
  those items rewrite the same files (`galaxymap3d.js`,
  `galaxystageview.js`); the rest follow MAP.52 to MAP.60.

  - [ ] **MAP.63 Shared map helpers in one module**
    The same helpers are copied between `galaxymap3d.js`,
    `sectormap.js` and `systemmap.js`: `readSceneData`, `cssVar`,
    `isLightBackground`, `addField`, `formatAddress`,
    `makeRingTexture`, `niceScaleValue`/`updateScaleBar`, `resize`, and
    the screen-space point pick (`starAtClientPoint` and
    `pointAtClientPoint`). Done: one `static/mapcore.js` exports them
    and all three maps import it; no visible change.

  - [ ] **MAP.64 One camera and input controller**
    The Galaxy Map's stage view (`galaxystageview.js`: drag, pan,
    wheel, pinch, two-tap select, keys) and the Sector Map
    (`sectormap.js`: its own `THREE.Spherical` orbit, pointer and arrow
    key handlers, zoom buttons) each have their own. Done: one
    controller module (orbit, pan, zoom, pinch, keys, drag-or-click
    threshold) with a zoom policy each view sets (free, a short range,
    or locked, which is MAP.58's rule), used by both maps.

  - [ ] **MAP.65 One picking, hover and info-panel layer**
    Done: one module for raycast and screen-space picking, the hover
    highlight and tooltip (the Sector Map has none today) and the info
    panel (fields, Nav from/to, Use as destination, Generate buttons,
    bookmark ☆), fed by each view's objects. The Sector Map's info
    panel gains the ☆ the drill-down design left for later.

  - [ ] **MAP.66 The sector as the drill-down's last stage, on the same page**
    Today clicking a generated sector leaves `/galaxy` for
    `/sector/<id>` (a full page load), and the Sector Map's data is
    baked into its page. Done: a sector-contents endpoint (systems,
    phenomena, clouds, neighbours, in the Sector Map's existing JSON
    shape) that the Galaxy Map fetches, so the sector opens as one more
    stage with the same buttons, breadcrumb, Back/Forward, URL
    (`/galaxy?sector=<designation>`) and bookmarks; neighbouring
    sectors are a sideways step. `/sector/<id>` stays as the sector's
    page (tables, text, edit tools) with the same engine embedded.

  - [ ] **MAP.67 One URL and history scheme for every level**
    Galaxy stages, a sector, and (with MAP.62) a system and a body in
    one URL form (for example `?at=` for stages, `?sector=`, `?object=<ref>`),
    so Back, Forward, reload and bookmarks work the same at every level.

  - [ ] **MAP.68 Remove the old Sector Map code**
    Once the sector stage matches it (tests from TEST.70), `sectormap.js` is deleted and the docs
    (`docs/html-interface.md`, the drill-down design) describe the one
    engine.

- [ ] **MAP.62 A full 3D star system view with a free camera**
  Boss (2026-10-01 20:55Z): "rendering a star system as a full 3D
  movable free-camera motion view." Today the System Map
  (`systemmap.js`, `lib/systemmap.py`) is a flat SVG diagram. Done: a
  star system is drawn in 3D (its stars, planets, moons, belts and
  comets on their orbits, sizes and distances shown legibly, with a
  scale option), and the camera can be moved freely (orbit, pan, zoom,
  fly to a body), using the shared engine of MAP.61 and the shared
  picker of NAV.3. Clicking a body opens or selects it. The flat
  diagram stays available.

  - [ ] **MAP.69 A system scene endpoint with 3D orbits**
    Today the System Map is drawn in Python as SVG (`lib/systemmap.py`):
    orbits are circles, z is dropped, distances are log-scaled per
    scene, and there is no JSON for a system. The database already has
    what 3D needs: each planet's and moon's `orbital_inclination_deg`,
    `orbital_ascending_node_deg`, `orbital_phase_deg`, `position_*_km`
    and period; the binary pair's mutual orbit; comets' Kepler elements.
    Done: `GET /api/systems/<id>/scene` returning every star, planet,
    moon, belt and comet with its reference, radius, colour, orbit
    elements, current position and the epoch they are valid for
    (`orbit_simulation_state.last_updated_at`).

  - [ ] **MAP.70 Positions at any time**
    Planets and moons move on circular orbits (phase plus 360 times
    elapsed time over period), comets on Kepler orbits
    (`keplerMotion.py`), binaries on their mutual orbit. Done: one
    JavaScript module (with a Python twin used by NAV.6) that gives any
    body's position at a time, so the view can animate and NAV.6 can
    plan around where bodies will be. A time control (now, play,
    faster, pause) in the view.

  - [ ] **MAP.71 Scale modes that keep everything visible**
    True scale makes planets invisible dots. Done: true scale, plus a
    compressed distance scale (the System Map's log scale) and enlarged
    bodies, switchable, with a note saying which is shown; moons stay
    outside their planet's drawn size.

  - [ ] **MAP.72 Rendering at system scale**
    Distances run from kilometres to light-days, beyond float
    precision near the camera. Done: a camera-relative (floating
    origin) scene; orbit lines as 3D ellipses; stars with the glow of
    `bodyRendering.js`; belts as particle rings; labels that declutter;
    the heliopause shown as a faint sphere.

  - [ ] **MAP.73 Free camera on the shared engine**
    Done: orbit, pan, zoom and a fly mode (keys and touch) on MAP.61's
    controller, plus "fly to" a body (the drill-down's smooth flight),
    following a moving body, and a reset. Clicking a body selects it
    (NAV.3's picker); double-click flies to it.

  - [ ] **MAP.74 The 3D view on the system page, the flat diagram kept**
    Done: the system page offers 3D and Diagram (the current SVG,
    default the last one used); a screen-reader list of bodies stands in
    for the canvas; a fallback to the diagram where WebGL is missing.
    Open question: should 3D be the default?

- [ ] **MAP.77 Galaxy Map draws block divisions inside a picked slab before zooming to it (bug)**
  Boss (2026-10-01 21:11Z): "when zooming in and navigating from the
  galactic map, when selecting a slab, don't show the divisions between
  interior blocks, only show divisions between the slabs. Then on a
  slab show the divisions between the blocks." Done: while slabs are
  being picked (the wedge view, MAP.52), the map draws only the
  boundaries between slabs, with no lines between the blocks inside
  each slab; once the view is on one slab (MAP.53, MAP.56), it draws
  the divisions between that slab's blocks, which are the segments the
  user picks next. The faint wire ghost of the other slabs (MAP.59)
  shows their outlines only, never their blocks. This repeats at every
  level of the slab and segment ladder of MAP.56.
  Arc pick (MAP.85, Boss 2026-10-01 23:53Z): the arc view draws only the
  slab boundaries; nothing is outlined on the whole galaxy except the
  hovered arc and its neighbors.

- [ ] **MAP.78 Zooming into a wedge must show the whole wedge at every drill-down level (bug)**
  Boss (2026-10-01 21:11Z): "when the system zooms into a wedge, make
  sure it is the entire wedge as you drill down." Done: whenever the
  view zooms to a wedge (MAP.52) or to a slab or segment inside it
  (MAP.56), the zoom frames all of what was picked, with no part cropped
  by the map's edges or by the controls over it; the fit uses the map's
  actual size (MAP.53) and holds while the view rotates and when the
  window is resized. Under MAP.58's locked zoom, the locked level is
  this whole-wedge fit, not a closer one.
  Arc pick (MAP.85, Boss 2026-10-01 23:53Z): read "wedge" as the picked
  arc: the zoom frames the whole arc, then the whole slab or segment.

- [ ] **MAP.79 Rogue planets clog the Sector Map: dim them, and a show/hide button per kind of object (bug)**
  Boss (2026-10-01 21:15Z): "rogue plants are just, everyhere and clog up the screen,
  make each dim, visible but the points for stars, comets, and other
  objects should shine through. Or let's say provide a button that
  turns each phenomena on and off in the sector map." Today every rogue
  planet gets a bright glow and a fixed-size ring on the Sector Map
  (`sectormap.js`, MAP.46), so in a busy sector they cover the stars.
  Done (default taken, Boss's second wording): the Sector Map has one
  toggle button per kind of object it draws (stars, rogue planets,
  interstellar comets, black holes, neutron stars, nebulae and so on),
  all on by default; turning one off hides those points, their rings
  and their labels, and hidden kinds can't be hovered or picked. The
  choice is kept in the URL so a bookmark keeps it. Ties in with
  MAP.61 (one control set for both maps). Boss (21:17Z): "Toggle and
  dim", so rogue planets are also drawn dim by default (a faint point,
  no bright glow or ring) while they are on, and stars, comets and
  other objects show through them.

  - [ ] **MAP.82 Unmarked rogue planets barely visible (bug)**
    Boss (2026-10-01 23:53Z): "Rogue planet detail is for them to be dim
    barely noticeable." Today every rogue planet has a bright violet
    core and glow (`sectormap.js`, `{power 2.0, strength 1.6}`) and at
    least a 4 px radius (`lib/starmap.py`). Done: unmarked rogue
    planets are a dim, small point with no glow, barely noticeable
    against the background, while stars and other objects show through.

  - [ ] **MAP.83 The "Mark rogue planets" button shows when it is on (bug)**
    Boss (2026-10-01 23:53Z): ""Mark Rogue Planets" button should
    retain a highlight if it is "on" and loose the highlight when it is
    "off" (default)". Today the button (`lib/starmap.py`) starts with
    `aria-pressed="true"`, and `sectormap.js` flips `aria-pressed` but no
    style follows it (`.starmap-btn-active` exists in `style.css` but is
    never applied), so on and off look the same. Done: off is the
    default; while on, the button keeps a clear highlight in both themes
    (styled from `aria-pressed`), and it loses it when turned off.

  - [ ] **MAP.84 Marked rogue planets grow and become clickable; unmarked ones stay small (bug)**
    Boss (2026-10-01 23:53Z): "Rogue planets should not only be dim and
    hard to see when not "marked" but also should be physically
    smaller. When "marked" they get bigger and more prominent and
    clickable." Today the toggle only shows or hides a ring sprite; the
    planet itself never changes. Done: unmarked rogue planets are drawn
    smaller than stars (with MAP.82's dimness); marking them makes them
    bigger, brighter and ringed, and only then easy to hover and pick
    (a larger hit area). Proposed values from Boss's research notes
    (tune on screen): unmarked about 1.5 px at 0.2 opacity, marked about
    5 px at full opacity with a glow ring.

- [ ] **MAP.80 Sector-level zoom on the Galaxy Map should show almost every star in the sector (bug)**
  Boss (2026-10-01 21:15Z): "as zooming into the sector level, when a sector is shown on
  the galactic arc it is close enough to see almost all stars in the
  sector including pulsars, quasars, and black holes." Today the Galaxy
  Map's level of detail (MAP.14, MAP.51) thins out the points drawn in
  a filled sector, so at the last drill-down stage a sector shows only
  some of its stars and its phenomena may be missing. Done: once the
  view is zoomed to sector level, the sector is drawn with nearly all
  of its stars and every pulsar, quasar and black hole in it, as close
  as the Sector Map shows them; the thinning only applies farther out.
  Ties in with MAP.66 (the sector as the drill-down's last stage).
  Boss (21:17Z) confirmed it is a fix: "fix it".

- [ ] **MAP.81 Ctrl+1 to Ctrl+9 bookmark keys clash with the browser's tab switching (bug)**
  Found by the bug audit (2026-10-01, `bug-audit.md`), first reported in PR #234's thread and left for Boss to decide:
  the map bookmark shortcuts (`static/bookmarks.js`, lines 23 and 325,
  MAP.23) use Ctrl+1 to Ctrl+9, which Chrome and Firefox on Windows and
  Linux take for switching tabs, so the shortcuts don't work there.
  Done: the bookmark keys use a combination no major browser reserves
  (for example Alt+Shift+1 to 9, or plain 1 to 9 while the map has
  focus), the help text says which, and a test pins it. Decided (Boss,
  2026-10-02): plain 1 to 9 while the map has focus.

- [ ] **MAP.87 Stars on the Sector Map and Galaxy Map need to be brighter, most of all the dim ones (bug)**
  Boss (2026-10-02 00:42Z): "ALL stars need to become about 4 times as
  bright in the sector maps but it's bright enough in the large galactic
  map so scale so that the dimmest red dwarf stars are 4 times as bright
  as they are now and when we approach the 1000+ sol lum mark it evens
  out to be not any brighter. That's just in how it's displayed."
  Today each star's point of light on the Sector Map comes from
  `_star_light` in `lib/starmap.py` (MAP.15): luminosity mapped on a log
  scale from 1e-4 to 1e6 L_sun onto the halo's size (13 to 40 px), its
  strength (`_LIGHT_GLOW`, 0.55 to 0.85) and the core's opacity
  (`_LIGHT_BRIGHT`, 0.9 to 1.0), drawn by `sectormap.js`; the Galaxy
  Map's stars use the same log range (`STAR_LOG_LUMINOSITY` in
  `galaxymap3d.js`). Boss (00:43Z): "Adjust TODO above to also use the
  same logic in the galactic map, on 2nd though". Done: display only,
  on both the Sector Map and the Galaxy Map, through one shared
  brightness curve (in Python with its JavaScript twin, or computed
  once and sent with the star data): the dimmest red dwarfs look about
  four times as bright as now, the boost shrinks smoothly with
  luminosity, and stars of about 1000 L_sun and up look as they do now;
  nothing stored changes. Before and after screenshots of a busy sector
  and of a zoomed Galaxy Map view, in both themes, go with the PR.

- [ ] **MAP.85 The galaxy pick is an arc, on a 3D galaxy with no sector lines**
  Boss (2026-10-01 23:53Z): "Redo the galactic selection, so that the
  galaxy map is 3D, we can manipulate it. The user doesn't select an
  entire wedge, just a large arc, then zoom in to select slab, and go
  from there. (edit all other TODO's about galactic interface to match
  this). Don't show any sector lines at so that the star map really
  comes through and the spiral pattern. When the user mouses over an
  arc they can select (each arc goes from top to bottom so slab of that
  we'll do in the next part) then they can see it's boundaries and the
  boundaries of the other segments." Today the first pick is a quarter
  of the disk (MAP.52 would make it a 40-degree wedge from the center to
  the edge), the map draws the wedge lines of the sector grid with
  bearing labels (`galaxymap3d.js`, the "Wedges" toggle), and every
  block is a shaded prism. Done:
  - The whole-galaxy view is a 3D galaxy the user can turn and tilt
    (rotation as MAP.53, zoom as MAP.58), drawn as its stars and spiral
    structure with no sector, block or wedge lines.
  - The first pick is an arc: a large piece of the disk bounded by
    bearing and by distance from the center, running the full height of
    the disk from top to bottom. Hovering shows the arc under the cursor
    with its boundary, and the boundaries of the neighboring arcs
    faintly; nothing else is outlined.
  - Clicking zooms to the arc (fitted to the window, MAP.78), where the
    user picks a slab (a height band of the arc, MAP.54 and MAP.59),
    then a segment of the slab (MAP.56), and on down to a sector.
  Default taken: an arc spans about 40 degrees of bearing (MAP.52's
  width, snapped to the grid's meridians) and a third of the disk's
  radius (inner, middle or outer), so the disk has about 27 arcs.
  Decided (Boss, 2026-10-02 01:53Z): this default. The other
  galaxy-map items carry an "Arc pick (MAP.85)" note saying how this
  changes them.

- [ ] **MAP.86 Sector and block colors from what is in them: filled sectors translucent (bug)**
  Boss (2026-10-01 23:53Z): "Filled in sectors should be translucent,
  just a hair more solid than the unfilled sectors, since they are a
  different color. Also make the color based on density averaged out
  with average star color and brightness of the stars within the
  sector. Then, that color will be averaged with the other sectors in a
  block (or mega block) to come up with that region's color. The scale
  is unfilled (color and opacity of an unfilled sector), filled but
  empty (more opaque and a shade more saturation), then the scale goes
  from empty to full (max possible density) in saturation and average
  color of the stars within the sector for hue. Average Luminosity
  compared to the sun to set the luminosity of the color of the sector.
  Blocks / Mega Blocks are then set by averaging the color and opacity
  of every sector in the block." Today (`galaxyblocks.js`,
  `galaxymap3d.js`) unfilled blocks use a density ramp at opacity 0.1 to
  0.3; a block with filled sectors is lifted to at least 0.6 opacity,
  fully opaque when all of it is filled; a one-sector block is colored
  bronze to gold by density; nothing uses the stars' colors or
  luminosity. Done:
  - Unfilled sector: today's unfilled color and opacity.
  - Filled but empty sector: a little more opaque and a shade more
    saturated than unfilled, still translucent.
  - Filled sector with stars: saturation from its star density (empty to
    the highest possible density), hue from the average color of its
    stars (by temperature), lightness from their average luminosity
    compared to the Sun; still translucent, only a little more solid
    than unfilled.
  - A block or mega block's color and opacity are the average of its
    sectors' (unfilled sectors counted as unfilled).
  - The per-sector color, saturation and lightness are stored or served
    with the tiles (a new API field, computed when a sector is saved), so
    the map doesn't need each star.
  Proposed numbers from Boss's research notes (tune on screen):
  unfilled opacity about 0.03, filled-empty 0.15, densest filled up to
  about 0.45. Ties in with MAP.85 (no lines, so color carries the
  structure) and MAP.59 (the ghost).

- [ ] **MAP.90 The tile-level helper crashes on a subnormal view radius (bug)**
  Low. Found by the debug-mode bug hunt (2026-10-02; report and evidence in the project's shared files under `bug-hunt/`) (deep fuzz). `galaxyViewport.tile_level_for_view_radius(2.2e-311)`,
  and its copy in `lib/galaxymap3d.py`, raise `OverflowError` because
  `log2` is infinite. Not reachable from a request as far as the hunt
  found (the radius comes from the galaxy shape). Done: it returns the
  finest level for any tiny positive radius.

## NAV: Navigation and courses

- [ ] **NAV.3 One shared picker for the Galaxy, Sector and System displays**
  Boss (2026-10-01 20:55Z): "completely functionalize all functions
  for the Galactic Picker and join it up with functionalizing the
  Sector Display and Star System Display so that the user, upon
  clicking on things in the navigation segment, can: go from one
  specific item (say, a moon in a star system) and go back out;
  navigate visually via the UI; select another sector, another star
  system, stellar phenomena, or anything like that, and vice versa to
  go from one to the other, so that we don't have to repeat code for
  the visual interfaces." Today the NAV page (`web/nav_page.py`) picks
  endpoints from dropdowns (a sector, then a system), and the Galaxy
  Map's drill-down, the Sector Map and the System Map each have their
  own picking code. Done: one picker module, shared by every visual
  display and the NAV page, that can:
  - pick any object at any level: a sector, a star system, a
    phenomenon, a star, a planet or a moon;
  - step out from any item to its parents (moon to planet to system to
    sector to the galaxy) and back in, visually and through a
    breadcrumb;
  - move sideways from one item to another of any kind.
  The NAV page uses it for both endpoints, so a course's ends are
  picked on the maps. Ties in with MAP.61 and MAP.62, and NAV.4 to
  NAV.6.
  Needs NAV.7. Built on MAP.61's engine; the parts that don't need the
  engine (the picker module, the breadcrumb, pick mode) can start first.

  - [ ] **NAV.13 A picker module: select, step out, step in, step sideways**
    The stage view already does this inside the galaxy (Esc and
    Backspace go up, arrows move to a sibling). Done: one
    `static/picker.js` holding the current selection as an object
    reference, with `select(ref)`, `up()`, `into(ref)` and
    `sideways(dir)` (siblings from the resolver), firing one change
    event every display listens to; keys, clicks and the breadcrumb all
    go through it.

  - [ ] **NAV.14 One breadcrumb for every level**
    Galaxy, stages, sector, system, star or planet, moon: one
    breadcrumb component built from the reference's parent chain, the
    same on the Galaxy Map, the sector, the system page and the NAV
    page.

  - [ ] **NAV.15 Pick mode everywhere**
    Today pick mode (choose a NAV start or destination) exists on the
    Galaxy Map and the Sector Map only, and stops at systems and
    phenomena. Done: the same "Choosing a destination · Cancel" mode on
    every display, the System Map included, picking any object down to
    a moon, and returning to the NAV page with it.

  - [ ] **NAV.16 NAV endpoints can be any object**
    Today an endpoint is a whole system or a phenomenon, and a course
    always leaves the heliopause. Done: the NAV page and `/api/nav`
    take any object reference; a course between bodies has legs: from
    the body out of its system, between systems, and in to the body
    (System Local Frame inside a system, as navigation-frames.md
    describes; the frame exists in `navigation.py` but is never chosen
    today); a course inside one system is one in-system leg. Open
    question: does an in-system leg use the same warp and fold speeds,
    or sublight speeds (impulse)? Default: the same tables, with a note.

- [ ] **NAV.4 Save a course**
  Boss (2026-10-01 20:55Z): "add a to-do item where I can save a course
  as a user ... Actually just have it save both so the user has either,
  no matter what they wanted in the first place." Today a course
  (`/nav?from=...&to=...`, `queryDb.nav_between`) can only be
  bookmarked as a URL. Done: a signed-in user can save a course under a
  name. A saved course always keeps both forms:
  - the direct, point-to-point line (its bearing and mark, NAV.1);
  - the system-to-system route (the chain of systems it hops through).
  The user can list, open, rename and delete their saved courses. The
  saved course shows whichever form the user views, and they can switch
  between the two. Needs user accounts (USR.1) and their storage. Open
  question: until user accounts exist, should saving be admin-only, or
  per browser like bookmarks?
  Default taken for the open question (pre-planning thread): Default for
  the open question: per browser now (the same storage and menu as
  bookmarks), moved into the account when USR.7 lands. NAV.4 then does
  not wait for user accounts.

  - [ ] **NAV.17 A saved course record with both forms**
    Done: a saved course stores its two end references, a name, when it
    was saved, both forms (the direct line's bearing, mark and distance;
    the route's list of stop references and distance) and, with NAV.6,
    the adjusted path. Opening it recomputes both from the current
    galaxy and says if anything changed since it was saved (the
    correlative update moves systems along their galactic orbits, and
    a regenerated system can vanish).

  - [ ] **NAV.18 Save, list, open, rename and delete, per browser**
    Done: a Save course button on the NAV result; a Courses list
    (beside Bookmarks, from `bookmarks.js` or a sibling module) with
    open, rename and delete; the viewer switches between Direct and
    Route on a saved course.

  - [ ] **NAV.19 Saved courses in the account (after USR.7)**
    Done: a `user_courses` table in the control database, owned by an
    account; courses saved per browser can be imported into the account
    once; the API lists and edits a user's own courses only.

- [ ] **NAV.5 Show a course on the Galaxy Map**
  Boss (2026-10-01 20:55Z): "In the navigation screen I want to be able
  to view the route in the context of the galactic map, zoomed in as far
  as it can be zoomed in and still show the entire path. The course
  path, direct and system-to-system, should then be specially
  highlighted as a course." Today the NAV page draws its own flat map
  of the route (`lib/navmap.py`), and the Galaxy Map only takes
  `?course=` for an end point. Done: from the NAV page (and a saved
  course, NAV.4), the course opens on the Galaxy Map, zoomed in as far
  as it can be while showing the whole path. Both the direct line and
  the system-to-system route are drawn in a distinct course style,
  told apart from each other, with their end points and hops marked
  and clickable. The same works inside a sector. Ties in with MAP.58
  (zoom limits: a course view may need a fitted zoom outside the user
  range) and MAP.61.
  Much of this exists: `/galaxy?course=<from>,<to>` (MAP.27, 7.52.0)
  draws one line through the route's stops with the ends named. Missing:
  the direct line drawn apart from the route, a fitted zoom, and any
  course inside a sector.

  - [ ] **NAV.20 Draw the direct line and the route apart**
    Done: the direct line and the system-to-system route as two styles
    (for example solid and dashed, in the course colour), a legend, the
    hops ringed and clickable (each opens its system, through NAV.3's
    picker), and the course readout (distance, bearing and mark, times)
    beside the map.

  - [ ] **NAV.21 Fit the view to the whole course**
    Done: the course opens zoomed in as far as possible with every
    point of both paths on screen, fitted to the window (the same fit as
    MAP.53), refitted on resize. This view is exempt from MAP.58's zoom
    lock; the user can still rotate it.

  - [ ] **NAV.22 Courses inside a sector and a system**
    Today a course that stays in one sector opens the sector with no
    line. Done: the same course drawing at the sector stage (MAP.61)
    and in the 3D system view (MAP.62), and a course that spans levels
    shows its in-system legs when zoomed in.

  - [ ] **NAV.23 Open a saved course on the map**
    Done: a saved course (NAV.4) opens in this view, with Direct, Route
    and (NAV.6) Adjusted switchable.

- [ ] **NAV.6 Courses that steer clear of gravity wells**
  Boss (2026-10-01 20:55Z): "we need to factor gravitational bodies into
  the course. A ship piloting would adjust the course to avoid falling
  into the gravitational field of objects it knows about. We'll need to
  have the system automatically adjust the course to avoid objects,
  attempting to stay out of the Hill sphere of each object. We're also
  going to use this within the sector and within the star system."
  Today a course is a straight line between its ends (NAV.1), and the
  route is a chain of systems. Done: course planning finds the bodies
  near the path that it knows about, and bends the path to stay outside
  each one's Hill sphere (or a safe radius where a Hill sphere doesn't
  apply, such as a star in the galaxy or a black hole). This works:
  - between systems in the galaxy (stars, black holes, neutron stars,
    nebulae and other phenomena);
  - inside a sector;
  - inside a star system (planets and moons, around the star).
  The adjusted path, its extra length and its time at each speed (NAV.2)
  are shown with the course, and the straight line stays available for
  comparison. Open questions: a Hill sphere needs an orbit around a
  heavier body, so what radius applies to a star or a lone object? And
  does the system level need the bodies' positions at a given time?
  Pre-planning (default taken for the first open question): The data is
  mostly there: planets and moons store `hill_radius_km`; each star
  stores `system_perimeter_km`, its Hill radius against the galaxy's
  tide (`spaceSector.hill_radius_ly`, already used to keep systems apart
  when they are placed); black holes, neutron stars and quasars have
  masses, so the same formula gives theirs. That answers the item's
  first open question: a star or lone object uses its galactic Hill
  radius.

  - [ ] **NAV.24 A keep-out radius for every kind of object**
    Done: one function giving each object's keep-out radius: a planet
    or moon its Hill radius; a star or system its stored
    `system_perimeter_km` (the wider of the pair for a binary); black
    holes, neutron stars and quasars the galactic Hill radius from their
    mass (`utils.calculate_hill_sphere`); a rogue planet the same.
    Nebulae, remnants and asteroid fields have no mass stored, only
    `radius_ly`. Open question for Boss: should courses avoid them too
    (as hazards, not gravity), or pass through? Default: pass through,
    with a note on the course.

  - [ ] **NAV.25 Find the obstacles along a path**
    Done: from the corridor query of NAV.10 (and NAV.38's
    line-to-sectors helper), every
    object whose keep-out sphere comes within reach of the path, at
    each scale: between systems, inside a sector, inside a system.

  - [ ] **NAV.26 Bend the path around keep-out spheres**
    Done: a path planner that leaves the straight line only where it
    crosses a keep-out sphere, going around it on the shortest detour
    (a tangent arc, or waypoints just outside the sphere), checked
    again against the other spheres; the end bodies' own spheres are
    exempt (a ship has to enter them to arrive). Unit tests with
    hand-placed spheres.

  - [ ] **NAV.27 Moving bodies inside a system**
    Inside a system the planets move. Done: the in-system planner uses
    positions at the time of travel (MAP.62's positions-at-time module,
    Python twin), from the departure time and the leg's speed. That
    answers the item's second open question: yes. Default departure:
    now.

  - [ ] **NAV.28 Show and save the adjusted course**
    Done: the adjusted path, its extra length and its times beside the
    straight line on the NAV page and on the map (NAV.5); saved courses
    (NAV.4) keep it as a third form.

- [ ] **NAV.7 One reference for every object, with its parents**
  Today only systems and phenomena can be named in a URL or a NAV
  endpoint (`nav_page.endpoint`, `<kind>:<id>`); stars, planets, moons,
  belts and comets have no page, no URL and no reference of their own
  (search links them to their system page), and nothing returns an
  object's chain of parents. Done: one reference form for every object
  kind (`sector:<id>`, `system:<id>`, `star:<id>`, `planet:<id>`,
  `moon:<id>`, `belt:<id>`, `comet:<id>` and the phenomenon types, the
  same strings NAV uses today, so old links keep working); one resolver
  (`queryDb` plus `GET /api/objects/<ref>`) that returns the object's
  kind, name, parent chain up to the galaxy (moon, planet, system,
  sector, galaxy), its siblings' references, and its position in each
  frame that applies (galaxy pc, sector-local ly, system-local km); and
  Python and JavaScript helpers that parse and print references. The
  picker (NAV.3), saved courses (NAV.4), the course on the map (NAV.5),
  gravity-aware courses (NAV.6) and the 3D system view (MAP.62) all use
  it.

  - [ ] **NAV.8 Pages and anchors for stars, planets, moons and belts**
    A star, planet, moon, belt or comet reference opens something:
    by default the system page scrolled to and highlighting that body
    (`/system/<id>#planet-<id>`), with its own System Map scene
    selected, rather than a new page per body. Search results, the
    locate box and bookmarks link this way. Open question: should
    planets and moons get pages of their own later?

  - [ ] **NAV.9 Search and locate return references for every kind**
    `/api/search` and `/galaxy/locate` (`queryDb.galaxy_locate`) return
    each hit's reference and parent chain, so any picker can jump to a
    star, planet or moon by name.

- [ ] **NAV.10 Routing that scales past a few thousand systems**
  Today `queryDb.nav_between` rebuilds the whole k-nearest-neighbour
  graph (`navGraph.build_knn_adjacency`, k = 6, an in-memory k-d tree)
  from every placed system on every galaxy-scope request, then runs
  Dijkstra; no position column has an index, and nothing finds the
  systems or bodies near a line. Done: the route search loads only the
  systems in a corridor around the direct line (a box query on indexed
  sector centers, widened if no route is found), or reads the stored
  `nearest_systems` table instead of rebuilding the graph; A* with the
  straight-line distance as its heuristic; a query that returns every
  system, star and phenomenon within a given distance of a line segment
  (used by NAV.6); and a measured time on a 100,000-sector database.
  Needs a galaxy schema migration for the position indexes. The
  hop-length study measured the full graph rebuild at 16 s for 200,000
  systems. NAV.34 (joining the graph's islands) comes first, and NAV.12
  (a route always exists, no hop limit) is built with it.

  - [ ] **NAV.11 Travel times for the system-to-system route too**
    Today warp and fold times are shown only for the direct distance;
    the route shows only its length. Done: the route gets the same warp
    and fold tables, per hop and in total.

- [ ] **NAV.29 Replace "Nav from here" and "Nav to here" with "Start Here" and "End Here" while picking (bug)**
  Boss (2026-10-01 21:15Z): "when navigating the "nav from and have to" buttons take you
  to different pages, so should be replaced by "Star Here" Or "End
  Here" and then the user goes to the next stage navigating back out
  from where their start or end is." Today a system or phenomenon's
  info panel on the maps offers "Nav from here" and "Nav to here"
  (`appendNavActions` in `sectormap.js`), which jump to the NAV page's
  own pickers and leave the map. Done: while picking a course, the panel
  offers "Start Here" (or "End Here" once the start is set); choosing
  it keeps the user on the map, and they pick the other end by stepping
  back out from where the first end is (NAV.13's step out and step in)
  and in again, without a page change. Outside pick mode the panel
  offers the same two buttons, which start the course from that object.
  Ties in with NAV.3, NAV.13 and NAV.15.

- [ ] **NAV.30 Hide "View phenomenon" and "View system" links while picking a course (bug)**
  Boss (2026-10-01 21:15Z): "don't show the view phenomena when navigating as it'll take
  you out of the page." Today the info panel shows "View phenomenon →"
  (and "View system →") in pick mode too (`sectormap.js`,
  `galaxymap3d.js`), and following it drops the course being built.
  Done: in pick mode the panel shows only the pick buttons (NAV.29),
  no link that leaves the picking flow. Ties in with NAV.15.

- [ ] **NAV.31 Galaxy wedges don't highlight on the navigation screens (bug)**
  Boss (2026-10-01 21:15Z): "in the navigation screen the wedges of the galaxy do not
  highlight at all and they should." Done: when picking a course on
  the Galaxy Map, hovering highlights the wedge under the cursor the
  same way the Galaxy Map does outside pick mode (MAP.52), and every
  later stage's hover highlight works too. Ties in with NAV.32.
  Arc pick (MAP.85, Boss 2026-10-01 23:53Z): the hover highlight on the
  navigation screens is the arc highlight of MAP.85.

- [ ] **NAV.32 Every Galaxy and Sector Map control works on the navigation screens (bug)**
  Boss (2026-10-01 21:15Z): "All the same UX from the galaxy screen and sector screens
  should be functional in the nav screens." Done: picking a course
  uses the same maps with the same controls as browsing them: hover
  highlight, wedge, slab and segment picks, zoom, rotate, the slab
  buttons, toggles (MAP.79), breadcrumb and bookmarkable URLs; pick
  mode only adds the Start Here and End Here buttons (NAV.29) and hides
  the links that leave the page (NAV.30). Best done by building pick
  mode on the one engine (MAP.61) and the shared picker (NAV.3, NAV.15)
  rather than as a separate copy; until then each map fix must be
  checked in pick mode too.

- [ ] **NAV.33 After picking one end of a course, stay at that zoom level (bug)**
  Boss (2026-10-01 21:17Z): "when nevigating via the picker, when we
  pick a start or destination first, it should keep us at that zoom
  level and let the user zoom out to find their destination via the
  picker." Today picking one end of a course moves the user to the NAV
  page's own pickers, away from the map view where they picked it.
  Done: after the first end is picked (with NAV.29's Start Here or End
  Here, or the existing "Use as start" and "Use as destination"
  buttons), the view stays where it was, at the same zoom level, with
  that end marked; the user zooms or steps out from there (NAV.13) to
  find the other end with the same picker, and the course is shown
  once both ends are set. Ties in with NAV.3, NAV.29 and NAV.32.

- [ ] **NAV.12 No maximum hop length: a route always reaches the nearest star it can, across any number of sectors**
  Boss (2026-10-02 01:53Z): "routs will always find the nearest star
  they can even if it crosses sector boundaries even across multiple
  sectors. This is for a game mechanic I need in place and I also want
  a jump through unknown space is marked in red and glows to draw
  attention to it. Also UX change here to display the path
  horizontally and find a way to split it between multiple lines for
  mobile or limited displays." This replaces the hop-length study's
  "optional ship range" (the study has reported; report in the
  project's shared files under `nav-hop-length/`). Today two things
  stop it: a same-sector route uses only that sector's own systems
  (`queryDb.nav_between`, sector scope), and the k = 6
  nearest-neighbour graph splits into pieces (NAV.34), so the NAV page
  says "No route via adjacent systems could be found". Done: same-sector
  routes may leave the sector; the route graph always joins its pieces
  with the shortest link between them; the longest hop is shown; each
  hop is flagged as a jump through unknown space when its line crosses
  one or more unfilled (ungenerated) sectors (the default reading of
  "unknown space"; NAV.38 finds the sectors), and `/api/nav` returns
  the flag per hop. Built with NAV.10, which already rebuilds the
  routing. Phase 1 is its anchor: NAV.34, NAV.38 and TEST.79 come
  before it (phase 0), UX.35 runs alongside it (phase 1), and NAV.36
  and NAV.39 need it first (phase 2).

- [ ] **NAV.34 Courses between separately generated areas find no route: the route graph splits into islands (bug)**
  Found by the hop-length study (2026-10-02; report in the project's
  shared files under `nav-hop-length/`). `navGraph.build_knn_adjacency`
  links each system to its 6 nearest generated systems, so any
  separately generated area with 7 or more systems becomes an island
  with no links out: 2,000 generated sectors over the disk split into
  714 islands and a cross-galaxy course found no route; two generated
  neighbourhoods 1 kpc apart found none either; and where a route did
  cross a gap it could hide one huge hop (a 6,504 ly last hop). Done:
  the islands of the route graph are joined (each island linked to its
  nearest few islands; in the study, 6 nearest joined all 714 in 0.2 s
  and gave a cross-galaxy route 1.29 times the direct distance), so two
  placed endpoints always have a route; a test with separated generated
  areas finds one. Lands with or before NAV.10; NAV.12 builds on it.

- [ ] **NAV.36 Unknown-space jumps drawn red and glowing**
  Boss (2026-10-02 01:53Z): "a jump through unknown space is marked in
  red and glows to draw attention to it." Done: a hop NAV.12 flags as
  an unknown-space jump is drawn red with a glow in UX.35's horizontal
  route strip, on the NAV page's own map (`lib/navmap.py`) and on the
  Galaxy Map course (with NAV.20), with a legend entry; with
  `prefers-reduced-motion` it stays red without the pulse; it reads in
  both themes, and the route list also labels it in text so it isn't
  shown by color alone. Prerequisites: NAV.12, UX.35.

- [ ] **NAV.38 Every sector a straight line passes through**
  A line-to-sectors helper, needed before NAV.12: given a straight
  segment between two points in the galaxy, every sector address it
  passes through, in `galaxyGeometry.py` with a JavaScript twin
  (`galaxyprisms.js`) and tests that the two agree. NAV.12 uses it for
  the unknown-space flag, and NAV.25 later for obstacles. Done: the
  helper, exact at sector faces and edges (after GEN.31 fixes the
  layer boundary), with tests along an axis, diagonally, through the
  core and out into the halo. Prerequisite: GEN.31.

- [ ] **NAV.39 Saved courses remember their unknown-space jumps and check them again**
  Done: a saved course (NAV.17) keeps which hops were unknown-space
  jumps when it was saved, and opening it checks again, since sectors
  may have been filled in the meantime; a hop that is now known shows
  as ordinary, and the course says what changed. Prerequisites: NAV.12,
  NAV.17.

## GEN: Generation and physics

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
  - The generate-around sphere and bright-star backfill (GEN.23, PR
    #226) run around the core too, as they do around any generated
    sector.
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
  GEN.27 and GEN.29 follow GEN.28.

  - [ ] **GEN.33 One class per PR, each with its tests**
    Suggested split: R and S (the commonest missing types) first; then
    U, W, X, Y and Z. Each adds the class to `PLANET_CLASSES`, the
    reference pages, the rogue flag, moon eligibility, and a
    distribution test over 1,000 generated systems.

- [ ] **GEN.29 Sweep every planet class for sense once the new ones are in (bug)**
  Boss (2026-10-01 15:26Z): "do a full sweep of planet classes to make
  sure they all make sense logically once the new classes are in
  place." After GEN.28. Done: every class's description, composition,
  zones, sizes, weights, temperatures and atmosphere agree with each
  other. Known oddities to settle: D allowed in the hot zone (icy
  bodies at 265-490 K); C a catch-all for 63% of cold planets; Q never
  generated (weight 0.0001, and orbits are circular); V's composition
  "iron, iridium, tungsten"; L with vegetation at a median 0.02 bar; E
  at 376-414 K, above water's boiling point at 0.6 bar.

- [ ] **GEN.31 A point just under layer 0's top face lands in layer 1 (bug)**
  Found by the generation tests (TEST.4-36 work, 2026-10-01): a point
  one float step below layer 0's top face is put in layer 1, both in
  the Python grid code (`galaxyGeometry`, `sector_address_at`) and in
  the map's `galaxyprisms.js`. Done: a point inside a layer's own
  height range always maps to that layer, in Python and JavaScript
  alike, and the strict xfail test for it passes. [MAP]

- [ ] **GEN.32 Re-running an interrupted bright-star band draws it twice (bug)**
  Found by the generation tests (2026-10-01): if `generate.py plan
  --bright-stars-down-to N` stops part way and is run again, the layers
  it already finished get the band a second time. Done: a re-run adds
  only the layers the interrupted run didn't finish (or starts the band
  over cleanly), never the same stars twice, and the strict xfail test
  for it passes.

- [ ] **GEN.34 Gas and ice giants come out too light, so there are no super-Jupiters (bug)**
  Found by the bug audit (2026-10-01, `bug-audit.md`), from the planet class gap report
  (2026-10-01); the thread held its physics fixes back because Boss
  hadn't approved them. Not re-measured on current main. The median bulk density of generated giants is about
  0.25 g/cm³, against 0.69 to 1.64 for real gas and ice giants, so
  massive giants (super-Jupiters) never appear. Done: giant masses and
  radii give real densities, super-Jupiters occur, and a test checks the
  density range over many seeds. Ties in with GEN.28 and GEN.29.

- [ ] **GEN.35 Rocky planets only ever get Class D moons (bug)**
  Found by the bug audit (2026-10-01, `bug-audit.md`), from the planet class gap report
  (2026-10-01); the thread held its physics fixes back because Boss
  hadn't approved them. Not re-measured on current main. `generate_moons` applies its moon size rule across the whole
  size range, so a rocky planet's moons all come out Class D. Done:
  rocky planets get the moon classes their size and zone allow, with a
  test over many seeds.

- [ ] **GEN.36 Moon regeneration can produce gas-giant or blacklisted moon classes (bug)**
  Found by the bug audit (2026-10-01, `bug-audit.md`), from the planet class gap report
  (2026-10-01); the thread held its physics fixes back because Boss
  hadn't approved them. Not re-measured on current main. `reconcile_zone_and_class` can regenerate a moon as a gas
  giant or as a class moons are never meant to have. Related to GEN.25
  (a reclassified moon too large for its planet), but a different
  failure. Done: a regenerated moon only ever gets a moon-eligible
  class, with a test.

- [ ] **GEN.37 97% of planets land in the cold zone (bug)**
  Found by the bug audit (2026-10-01, `bug-audit.md`), from the planet class gap report
  (2026-10-01); the thread held its physics fixes back because Boss
  hadn't approved them. Not re-measured on current main. Nearly every generated planet is in the cold zone, so hot
  and temperate planets are rare. Done: the zone mix is measured on
  current main, the orbit or zone placement is fixed to give a
  plausible spread, and a test checks the share over many seeds.

- [ ] **GEN.38 Rocky rogue planets over 10,000 km are still classed C (bug)**
  Found by the bug audit (2026-10-01, `bug-audit.md`), from the GEN.8/GEN.26 thread report (PR #263): Class C's size range
  tops out at 10,000 km and no other rogue-eligible rocky class exists,
  so bigger rocky rogues are classed C anyway. Done: they get a class
  that fits their size. May be solved by GEN.28's S class (rocky
  super-Earth) if S is rogue-eligible.

- [ ] **GEN.39 The same seed can't reproduce the same galaxy (bug)**
  Found by the bug audit (2026-10-01, `bug-audit.md`), from the parallel, population and navigation tests thread: star
  positions and star draws use the operating system's random source by
  design, so even at the same worker count one seed gives a different
  galaxy each run, and TEST.19 (retired with PR #321) couldn't check
  its goal of the same sectors at any worker count. Decided (Boss,
  2026-10-02 01:34Z): "one seed reproduces the same galaxy." Design:
  each sector gets its own seed, derived from the galaxy seed and the
  sector's address, and every draw for that sector (star positions,
  star draws, systems, phenomena, backfill) comes from it, so the order
  sectors run in and the worker count don't matter. Reproducible on the
  same version only: a release that changes generation may change what
  a seed makes. Storage and the log line are DB.6 and OPS.10 below; the
  end goal, a version and a seed rebuilding the same galaxy, is GEN.55.
  Seed shape (Boss, 2026-10-02 01:36Z: "I want to use a 128bit seed"): one 128-bit galaxy seed stored with the galaxy (a
  `BINARY(16)` column or a 32-character hex string, since `BIGINT
  UNSIGNED` holds only 64 bits) and shown as 32 hex digits, a new
  `--seed` option on `plan` and new galaxies that accepts it, and each
  sector, layer and backfill block seeded from the galaxy seed plus its
  address. Each unit's seed is SHA-256(128-bit galaxy seed || "kind:"
  address), for example `sector:12/3/0`; the version is not mixed into
  the hash but stored alongside (DB.6), and the bright-star scatter's own 63-bit seed
  is derived from the 128-bit galaxy seed too, so no step throws bits
  away. Builds after PERF.21, in the parallel path thread. Done: no
  generation draw uses the operating system's random source; a test
  checks that one galaxy seed gives the same sectors at 1, 2 and N
  workers (with TEST.74); and a test checks that two seeds differing
  only in their high 64 bits give different output. Decided (Boss,
  2026-10-02 01:46Z: "I'm going to nuke the galaxy anyway so let's say
  GEN.39 requires the galaxy to be nuked and a fresh start"): GEN.39
  starts from a wiped galaxy, which Boss does himself; no support for an
  unseeded galaxy and no migration of old galaxies is needed.

  - [ ] **DB.6 Store the galaxy's 128-bit seed, the version that made it, and every generation run**
    The storage half of GEN.39. Boss (2026-10-02 01:40Z): "Ok use a 128 bit value and store the seed in
    the database, and put it in the log at the top of any generation, also
    populate the TODO upward from here to eventually build a system that a
    version number and a seed value would reproduce the same galaxy by the
    end of the phases." Today the
    only stored seed is the bright-star scatter's 63-bit one
    (`_db.record_bright_star_scatter`); the run seed in `generate.py`
    main (`secrets.randbits(128)`) and the work queue's `run_seed` are
    never stored. Done: a galaxy schema migration adds the galaxy's seed
    to `galaxy_shape` (`BINARY(16)`, written once, when the galaxy is
    first planned, set with `generate.py plan --seed` as 32 hex digits or
    drawn at random) and the
    PlanetGen version, Python version and platform that made it, plus a
    `generation_runs` table with one row per run that changes the
    galaxy: the command and its options, the version, start and end
    times, and the outcome. Simplest default: the run history is what
    makes "replay" possible (OPS.12), since a galaxy is built by a series
    of commands (plan, sectors, scatter, backfill), not by the seed
    alone. A test checks the seed reads back bit for bit. No migration of
    existing galaxies: GEN.39 starts from a wiped galaxy.
    The version: Boss (2026-10-02 01:46Z): "I would say the version
    number all added together into a number stored in hex", so the
    stored version value combines MAJOR, REVISION and BUILD into one
    number in hex. A plain sum collides (7.127.352 and 7.128.351 both
    sum to 486), so the default, the collision-free form of "all added
    together", packs the parts: MAJOR<<32 | REVISION<<16 | BUILD, shown
    as 12 hex digits (7.127.352 is `0007007F0160`), stored with the
    galaxy (here) and with each sector (DB.7), with the full
    MAJOR.REVISION.BUILD string kept next to it. Decided (Boss,
    2026-10-02 01:53Z, "you are correct on all points"): the packed form.

  - [ ] **OPS.10 The galaxy seed and version at the top of every generation log**
    The log half of GEN.39 (Boss: "put it in the log at the top of any
    generation"). Today `generate.py` main logs "Seeded the random number
    generator with ... (no --seed option exists to reproduce this run)"
    at debug level only. Done: every `generate.py` subcommand, every job
    the web site or API starts, and every work queue run writes one line
    first, at normal level, to the console, the job log and the debug
    log: the galaxy seed as 32 hex digits, the PlanetGen version, and the
    run's command (for example `Galaxy seed 3f2a...c901, PlanetGen
    7.127.352, run: sector 12 3 0`); a test checks the line is first.

- [ ] **GEN.40 Weed out sectors by star density before the bright-star backfill**
  Boss (2026-10-01 22:22Z): "see if we can cut down the number of
  sectors that need to be backfilled by taking the star density
  calculation and using it to weed out sectors from the star
  scattering. Might manage probability of getting a region through
  eliminating x number of sectors from the region so that the region
  probability matches what we're looking for. One pass would eliminate
  sectors based on probability, then the back fill would literally have
  less work to do. We would need to make sure that we don't over-filter
  though because bright stars in odd places is really neat and
  realistic." Today the GEN.30 backfill (`backfill_bright_stars_around`
  in `generate.py`, `brightStars.backfill_cells`) works out the
  densities and makes a Poisson draw for every unfilled sector of every
  block in range, most of which get no star. Investigate first; the
  subitems split the work.

  - [ ] **GEN.41 Investigate: how much backfill work a density pre-pass would save**
    Measure how many sectors the backfill visits per run and how many
    end up with a star, by tier and by density, and work out how many a
    pass that drops sectors by probability would skip. Done: numbers
    written up for Boss, with a go or no-go for GEN.42.

  - [ ] **GEN.42 A pass that drops sectors from a region by probability**
    One pass over a region (a sector block, or a tier's shell) uses the
    star density calculation (`galaxyDensity.population_densities`) to
    eliminate some of its sectors, chosen so that the region's chance of
    getting each bright star is still what the density says; the
    backfill then draws only in the sectors left. Done: the pass runs
    before the backfill, the backfill visits fewer sectors, and the same
    seed still gives the same stars wherever it lands (see GEN.39).

  - [ ] **GEN.43 Don't over-filter: keep bright stars in odd places**
    The pass in GEN.42 must not strip the rare bright star from
    low-density places (between the arms, the halo, far out in the
    disk), which Boss calls "really neat and realistic". Done: no
    sector's chance drops to zero just for being low density, and a
    test over many seeds checks that the number and spread of bright
    stars, including in low-density regions, match the backfill without
    the pass.

  - [ ] **GEN.44 Store each sector's backfill level so finished sectors drop out of any backfill**
    Boss (2026-10-01 23:53Z): "stores every sector it's backfill level
    from -1 (infinite, no backfill), and if it has been back-filled but
    is not generated entirely it will store the lowest solar lum value
    it was backfilled to. If it is totally generated it stores a 0.,
    this way we can eliminate sectors from any backfill dynamically".
    Today the backfill's depth is kept per 3x3x3 block, not per sector
    (`bright_star_blocks.min_luminosity_sol`, schema v49, read and set
    by `_backfill_block` in `generate.py`), and a sector counts as
    generated only because its `sectors` row exists. Done: every sector
    address has a backfill level: -1 never backfilled, a positive
    number for the dimmest luminosity (L_sun) it was backfilled down to,
    0 once it is fully generated. A backfill at a given floor skips
    every sector at 0 or already at or below that floor, and sets the
    level of each sector it draws. Since unfilled sectors have no
    `sectors` row, the level lives in its own table keyed by
    (ring, layer, slot), or replaces `bright_star_blocks` (a schema
    migration, its row in `database-schema.md`). Decided (Boss,
    2026-10-02 01:46Z): "gen.44/perf.11 is one shared per sector stats
    table", so the level lives in the same per-sector stats table as
    PERF.11's densities, keyed by sector address (ring, layer, slot), so
    unfilled sectors have rows for their backfill level, built in one
    migration. Default: MAP.86's per-sector color, saturation and
    lightness go in the same table. No migration of existing galaxies (GEN.39 starts fresh).
    GEN.40 to GEN.43 and PERF.18 use it to skip work.

- [ ] **GEN.45 Check the rogue planet mix of terrestrial and gas giants (bug)**
  Boss (2026-10-01 23:53Z): "investigate probability for a gas giant
  rogue planet vs terrestrial, verify we are actually seeing the
  results we expect." Today `roguePlanetData.py` picks a mass bin by
  weight (`ROGUE_PLANET_MASS_BINS`: terrestrial 0.1-2 Earth masses
  weight 5, sub-Neptune 2-20 weight 1, Saturn and Jupiter bins 0.25
  each), draws a log-uniform mass in it, and calls anything at or over
  0.05 Jupiter masses (about 16 Earth masses) a gas giant: about 91%
  terrestrial and 9% gas giants by arithmetic. The only test
  (`test_rogue_planets_are_mostly_terrestrial`) checks the terrestrial
  bin's share, not the terrestrial/gas split. Done: the split is
  measured over a large sample and compared with the expected one
  (microlensing surveys: free-floating planets are mostly Earth-mass to
  Neptune-mass, Jupiter-mass ones rarer); if they differ the weights
  are fixed; a test pins the split. Boss's research notes propose a
  power-law mass function, dN/dM proportional to M^-0.65 from 0.01
  Earth masses to 13 Jupiter masses, as the reference to check against.

- [ ] **GEN.46 Star system names of at most two words (bug)**
  Boss (2026-10-01 23:53Z): "Name generation should not produce star
  names that are more than 2 words long. This keeps planet names from
  getting too long to be reasonable." Today a base name is one word,
  or two after `split_long_word` (`utils.py`), but the uniqueness
  decorations add words (`nameUniqueness.py`, applied in `_db.py`): a
  Greek prefix (3 words), "Alpha <base> <Roman>" (4), a diminutive
  stacked outside ("Little Alpha Xy Zz IV", 5), and diminutives can
  stack again. Done: no star system name is longer than two words: a
  collision is resolved within two words (for example a Greek letter
  with a one-word base, or a new base name), and a test over many
  generated and decorated names checks it. Planet and moon names keep
  their numeral and letter (`bodyNames.py`). Decided (Boss,
  2026-10-02): only new names follow the rule; existing names stay as
  they are, with no renaming migration.

- [ ] **GEN.47 Nebulae almost never appear (bug)**
  Boss (2026-10-01 23:53Z): "No nebulae are being created at all."
  Checked in the code: nebulae are generated (`generate_sector_phenomena`
  in `generate.py` runs for galaxy-placed sectors too), but at rates
  that make them vanishingly rare in 4 pc sectors: about 3e-4 molecular
  clouds and 2e-6 planetary nebulae per sector (densities 5e-6 and 3e-8
  per pc^3 in `program_constants.py`), star-hosted nebulae only around
  O stars and some B and A stars (`NEBULA_HOST_RULES`), and diffuse gas
  (classes A and B) deliberately not generated. Nebulae are also many
  parsecs across, far bigger than a sector, so per-sector rolls don't
  fit them. Done: nebulae appear at realistic numbers across a
  generated region: large clouds placed once per region at galaxy
  scale (for example a density field with seeded centers, more in the
  arms), spanning the sectors they cover, plus the star-hosted ones; a
  test over a generated neighborhood finds them. Boss's research notes
  suggest a 3D noise density field with Poisson-seeded centers and
  per-class thresholds and radii (emission 15-45 pc near O/B stars,
  reflection 10-30 pc, dark 5-25 pc, planetary 0.1-2 pc, remnants
  5-20 pc). The field comes from the galaxy seed through GEN.39's
  per-region seeds, so every worker and every later run agrees where a
  nebula is (GEN.55). Design: [docs/design/nebula-and-asteroid-field-classes.md](design/nebula-and-asteroid-field-classes.md)

- [ ] **GEN.48 Forcing options are impractical for whole sectors; replace them with prevalence controls (bug)**
  Found by the debug-mode bug hunt (2026-10-02; report and evidence in the project's shared files under `bug-hunt/`). Boss (2026-10-02): "the forcing options are impractical
  for an entire sector, so those options only apply to generating a star
  system, a single star system. Replace it with a control for the
  prevalence of each as a probability adjustment ... as a % deviation
  from the normal probability. Split that into several TODO's. Log as a
  bug." The `+x`/`-x` options (habitable_world, asteroid_belt, comets,
  large_star, moons, max_planets, intelligent_life, binary_system,
  wide_binary, planets) force every system in a sector or galaxy run,
  which is impractical and on some stars impossible (GEN.49).

  - [ ] **GEN.49 `+habitable_world` silently fails on hot stars (bug)**
    After 8 tries (`MAX_SYSTEM_GENERATION_ATTEMPTS`), `StarSystem` keeps a
    system with no habitable world, exits 0 and saves it. With 40 systems
    per type it failed for O5V 17/40, B3V 8/40, M2IA 6/40 and A0V 2/40;
    the CLI `system --star-type O5V +habitable_world` failed 3 times in 6.
    The debug log shows `Planet generation attempt 8/8: retrying
    (habitable world required and found: False ...)` and then the save.
    Done: for a single system, a forced option the star can't meet is
    refused up front (as `-large_star +habitable_world +asteroid_belt`
    already is) or fails with a clear error; it never saves silently.

  - [ ] **GEN.50 `-planets +asteroid_belt` still makes an asteroid belt (bug)**
    `-planets` is documented as "skips the planet generation process
    entirely" (`config.py`), but `+asteroid_belt` still places a belt (3
    of 3 CLI runs; deep fuzz `test_random_valid_configs_generate_sane_systems`).
    Done: contradictory forcing options are rejected for a single system.

  - [ ] **GEN.51 Forcing options only for single-system generation**
    Remove the `+x`/`-x` options from `generate.py sector` and `galaxy`,
    keep them on `system` (and the one-off system page). A saved system
    config or file that still carries them gets a clear message.

  - [ ] **GEN.52 Prevalence controls for sector and galaxy runs**
    Each former forcing option gets a prevalence setting: a percentage
    deviation from its normal probability (+50% means 1.5 times the usual
    chance that a system has comets; -100% means never). It applies to
    every system the run generates, through the CLI (`--prevalence
    comets=+50`, or one option per feature) and the stored system config.

  - [ ] **ADM.16 Prevalence controls on the Generate page**
    The Generate page's sector, galaxy and "around a sector" forms swap
    their forcing checkboxes for the prevalence fields (GEN.52),
    defaulting to 0%. The one-off system page keeps forcing.

  - [ ] **TEST.75 Tests for forcing and prevalence**
    Single systems: every forced option holds or is refused. Sector runs:
    the measured shares across many seeds move by the requested
    percentage within a tolerance. [generation]

- [ ] **GEN.53 The two stars of a binary don't share one age (bug)**
  Found by the debug-mode bug hunt (2026-10-02; report and evidence in the project's shared files under `bug-hunt/`). (a) With `--star-type`, the secondary is a new
  `Star(mass_override=...)` with its own random age: 194 of 200 pairs
  differed (a G2V of 5.26 Gy next to 9.39 Gy). (b) Population-model pairs
  start with one age, but `adjust_age_for_planets` then ages each star
  of a wide pair separately (12 of 60 wide pairs, for example M2V 13.80
  Gy and M8V 0.45 Gy), and a close pair's adjusted proxy age goes to the
  primary only (1 of 60). Done: both stars always share one age, checked
  over many seeds.

- [ ] **GEN.54 A `--star-type` secondary gets a mass that doesn't fit its type (bug)**
  Found by the debug-mode bug hunt (2026-10-02; report and evidence in the project's shared files under `bug-hunt/`). The secondary's mass is `primary.mass * uniform(0.1,
  0.8)`, clamped only to the whole Yerkes-class range, while its
  temperature and luminosity are still drawn for the requested type: for
  example a "G2V" of 0.17 Msun and 0.64 Lsun, and 191 of 200 secondaries
  more than 10% off the mass-luminosity relation. The comment in
  `systemData.py` says the mass is "clamped into its own class's range",
  but it isn't. Done: the secondary's type comes from its mass, or its
  mass is held to its type's range, and its luminosity follows from the
  mass.

- [ ] **GEN.55 A version number and a seed reproduce the same galaxy (end goal)**
  Boss (2026-10-02 01:40Z): "Ok use a 128 bit value and store the seed in
  the database, and put it in the log at the top of any generation, also
  populate the TODO upward from here to eventually build a system that a
  version number and a seed value would reproduce the same galaxy by the
  end of the phases." This item is the
  chain that gets there; it is done when its last sub-item is. Order:
  - Phase 0: GEN.39 (per-unit seeds, after PERF.21) with DB.6 and
    OPS.10.
  - Phase 1: OPS.11 (what "the same galaxy" means), GEN.56 (every
    draw seeded), GEN.57 (a sector's contents depend only on the seed,
    the version and its address), DB.7 (the version kept with each
    sector), GEN.58 (a fingerprint) and TEST.77 (the golden-galaxy
    test). GEN.47's nebula field uses the derived seeds.
  - Phase 2: PERF.18 and GEN.42 give the one-process stars for one
    seed (already in their Done text); API.16 and ADM.17 show the seed
    and version.
  - Phase 3: API.17 makes remote generation reproduce the server's;
    GEN.59 records edits and time evolution as layers on top of the
    seed.
  - Phase 3+ (the end state): OPS.12, `generate.py reproduce`.
  Anything that draws new randomness later (GEN.47, GEN.42, PERF.18,
  API.12, API.13) uses the derived seeds and keeps TEST.77 green.

  - [ ] **OPS.11 Define "the same galaxy" and which versions stay reproducible**
    A short design note (`docs/design/reproducible-galaxies.md`) that
    the other items build to. Simplest defaults, unless Boss changes
    them: "the same galaxy" means every generated object has the same
    address, position, properties and name; database ids, timestamps,
    population data rebuilt later and admin edits (the edit log keeps
    those) are not compared. A seed reproduces a galaxy only on the
    exact PlanetGen version that made it (recorded by DB.6); an older
    galaxy is reproduced by checking out its version. Python's version
    and platform are recorded too, and TEST.77 runs on every CI Python
    leg so any difference between them shows up. Known risk: the math
    library (libm) can differ in the last digit between machines, which
    can change a result near a threshold; the golden-galaxy test
    (TEST.77) watches for it. A release whose generation output changes
    says so in its `changes/` note. The math behind "seed + version
    gives the same galaxy" is in the seed math report (dependency tree
    thread, 2026-10-02).

  - [ ] **GEN.56 Every random draw in generation comes from the derived seeds**
    The sweep GEN.39's "no draw uses the operating system's random
    source" needs; found by reading the code on 2026-10-02.
    `secrets.choice` and `secrets.randbelow` draw generation choices in
    `planetPhysics.py` (planet class, 3 places), `planetLife.py` (5),
    `planetData.py` (moons) and `utils.py` (name syllables, prefixes and
    suffixes, 5); `spaceSector._rng` is a `secrets.SystemRandom` for
    every star position; `utils.reseed_rng()` reseeds the global
    `random` from `secrets` and is called from 13 generator modules
    (stars, comets, asteroids, nebulae, remnants, quasars, rogue planets
    and others); the bright-star scatter and band seeds are
    `random.SystemRandom().getrandbits(63)` (`generate.py`); and the run
    seed is `secrets.randbits(128)`. Done: each unit (sector fill,
    scatter layer, backfill block, phenomenon and population pass) gets
    seed = SHA-256(galaxy seed + unit address), the full digest seeds a
    `random.Random` that is passed down explicitly, and no generator
    module calls the module-level `random` or `secrets`; and a test scans the
    generation modules and fails on `secrets`, `SystemRandom`,
    `os.urandom`, `uuid4`, `time` or an unseeded `random.seed()` used
    for a draw, with an allowlist for login, CSRF, API keys and job ids.
    Python only promises that `random.random()` gives the same numbers
    across Python versions, not `choice`, `uniform`, `gauss` or
    `shuffle`, so generators draw through small in-house helpers built
    on `random()` alone (choice, uniform, normal, shuffle), and the test
    above also fails on direct calls to the others. Prerequisite: GEN.39.

  - [ ] **GEN.57 A sector's contents depend only on the seed, the version and its address**
    Seeding every draw (GEN.56) isn't enough when the result depends on
    what ran first, at what worker count or how fast. Known order
    dependencies: name-collision decorations are resolved by how many
    names already exist (`resolve_greek_roman_collision(base_name,
    existing_count)` in `nameUniqueness.py` and the reservation in
    `_db.py`), so the sector that reserves first gets the plain name
    (settle with GEN.46); population seeds key on database ids
    (`random.Random(planet_id)` and `random.Random(species_id * 7919)` in
    `population.py`), and ids come from per-worker id blocks; the
    backfill depends on which sectors are already filled and where the
    run started (with GEN.44); and nearest-system links depend on which
    neighbours exist yet. Done: a sector's contents (as OPS.11 defines
    them; ids and timestamps excluded) come out the same whichever
    sectors were generated before it, at any worker count; names, seeds
    and skips key on addresses and the galaxy seed, not on ids or arrival
    order: in a name collision the sector with the lower address keeps
    the name and the other is renamed from its own seeded stream (with
    GEN.46);
    nearest-system links are rebuilt from content, so they are left out
    of the comparison; and a test generates the same sectors at 1 and 4
    workers and in two orders and compares fingerprints (GEN.58).
    Prerequisites: PERF.21, GEN.39, GEN.56.

  - [ ] **DB.7 The version that generated each sector, and a warning for mixed-version galaxies**
    Done: each sector row records the PlanetGen release that generated
    it, as DB.6's packed hex version and the full string (DB.6 keeps the
    galaxy's first version); extending a galaxy with
    a different release warns before it starts (CLI and Generate page),
    because a mixed-version galaxy reproduces only sector by sector, each
    on its own version. Prerequisite: DB.6.

  - [ ] **GEN.58 A fingerprint of a galaxy's generated content**
    A way to tell whether two builds are the same in the sense OPS.11
    defines. Done: `generate.py fingerprint` (for the galaxy or a region)
    prints a canonical SHA-256 digest per sector and one for the region over the compared content in
    a fixed order (address order, canonical number formatting), skipping
    ids, timestamps and edits; the same function backs GEN.57's test,
    TEST.77 and OPS.12. Prerequisite: OPS.11.

  - [ ] **TEST.77 A golden-seed regression test**
    Done: a fixed 128-bit seed builds a small galaxy (plan, a few
    sectors, a scatter and a backfill) at 1 and 4 workers, and its
    fingerprint (GEN.58) must match the one pinned in
    `src/tests/golden/`, pinned per release with the version it was made
    on. When generation
    output changes, the test fails until the PR updates the pinned
    fingerprint and its `changes/` note says generation output changed
    (`bump_version.py --check` checks the two go together). It runs on
    every CI Python leg. Prerequisites: GEN.57, GEN.58. [generation, infra]

  - [ ] **OPS.12 `generate.py reproduce`: a version and a seed rebuild a galaxy and check it**
    The end state Boss asked for (phase 3+). Done: `generate.py reproduce
    --seed X --version Y` rebuilds a galaxy, or a region, into a fresh
    database from the seed and the stored run history (DB.6), and checks
    it against the fingerprint (GEN.58) of the live galaxy or a given
    one, listing any sector that differs; it refuses, naming the version
    to check out, when the running release isn't Y. Without the edit and
    epoch layers (GEN.59) it rebuilds the galaxy as first generated; with
    them, as it is now. A test runs it on a small galaxy, and on one with
    a deliberately changed sector. Simplest default: no automatic
    migration of old galaxies to a new release's output. Prerequisites:
    DB.6, DB.7, GEN.57, GEN.58, TEST.77, GEN.59.

  - [ ] **ADM.17 The Generate page shows the galaxy's seed and version**
    Done: the Generate page shows the galaxy seed (32 hex digits, with a
    copy button), the version that made the galaxy and the run history
    (DB.6), and the new-galaxy form takes an optional seed (blank means a
    random one); the admin's System page shows the same for one system.
    Prerequisite: DB.6.

  - [ ] **API.16 The API reports the galaxy's seed, version and run history**
    Done: an API route returns the galaxy seed, the version that made it
    and the run history (DB.6), documented with the API; API.12's
    download uses the same fields. Prerequisites: DB.6, API.5.

  - [ ] **API.17 Remote generation reproduces what the server would make**
    Done: a remote run (API.12's download, API.13's generation without a
    database) with the same seed and release produces exactly what the
    server would for those sectors, checked by fingerprint (GEN.58); and
    API.8 can verify an upload by re-running a sample of its sectors on
    the server and comparing. Prerequisites: API.12, API.13, GEN.57,
    GEN.58.

  - [ ] **GEN.59 Edits and time evolution recorded as layers on top of the seed**
    What sits on top of generation is recorded separately so it can be
    replayed: admin edits, overrides and regenerations (`editStore.py`,
    the edit log) and time evolution (`updateOrbits.py`, the correlative
    update) are logged as ordered layers with the release that applied
    them. Done: "seed + version" rebuilds a galaxy as first generated,
    and "seed + version + edit log + epoch" rebuilds it as it is now;
    a regeneration draws from SHA-256(sector seed || edit number), not
    `random.seed()` (today in `api/edits.py`), so replaying the edit log
    replays it. Prerequisites: GEN.56,
    GEN.58.

## PERF: Speed, caching, bulk generation and parallel work

- [ ] **PERF.1 Generation at scale**
  Boss's notes of 2026-10-01 on bulk generation: estimates before it
  starts, progress while it runs, speed records, and parallel work. The
  code is mostly `generate.py`, `web/generate_page.py`, `web/jobs.py` and
  `stellarObjects/brightStars.py`.

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
    expected-versus-actual ratio to correct their estimates. Decided
    (Boss, 2026-10-02 01:46Z): one per-sector stats table shared with
    GEN.44's backfill level, built with it; no backfill of
    existing sectors (GEN.39 starts fresh). Open questions: whether
    "actual" counts systems, stars, or both; what the decaying average
    is taken over (the ratio per density bucket, so it ties in with
    PERF.10, or one galaxy-wide figure); and what happens to the stats when a sector is regenerated
    (ADM.8) or the galaxy is reset.

- [ ] **PERF.18 Run the GEN.30 bright-star backfill in parallel on the work queue**
  Boss (2026-10-01 22:13Z): "make the GEN.30 backfill parallelized and
  use the work queue." Today `backfill_bright_stars_around`
  (`generate.py`) runs in one process, one block at a time: after a
  sector run, in the block top-up loop in `add_bright_star_band`, and
  inline inside web requests (the map visit and the API neighborhood
  route). Each block takes a row lock and commits on its own, so blocks
  can be separate work-queue tasks. Done: the backfill's blocks run as
  tasks on the work queue (`stellarObjects/workQueue.py`, PERF.8) with
  the run's worker count, giving the same stars as the one-process
  backfill for the same seed, and the in-request backfills queue their
  work instead of doing it inside the request. The bright-star scatter
  already runs per layer in parallel (PERF.7) but stays serial within a
  layer.

- [ ] **PERF.19 Everything the API or web site starts runs on the work queue (investigate)**
  Boss (2026-10-01 22:13Z): "I want _everything_ that communicates
  with the API to use the work queue whenever possible. Add a TODO item
  about that and about we need to investigate that." Investigate first:
  list every generation or database-write path the API and web site
  trigger (API routes, remote generate, uploads, map-visit and
  neighborhood backfills, admin edits and regenerates, Generate page
  jobs), whether each runs inside the request or on the work queue
  today, and which should move. Done: the audit, written up with a plan
  that Boss approves; the moves themselves are filed as their own items
  from that plan.

- [ ] **PERF.20 Short-term caching through the work queue and API (needs planning)**
  Boss (2026-10-01 22:17Z): "Everything should go through the work
  queue so that we can implement short-term caching where possible, add
  a TODO item to investigate that so that the work queue and API work
  together to optimize API calls back and forth. This is an item that
  needs planning." Builds on PERF.19's audit of what the API and web
  site start and on the page cache (PERF.2). Plan first, with Boss:
  which API calls and queued results can be cached for a short time, and
  where (the queue's task results, the API's responses, or both); how
  long entries live and what clears them (a write, a regenerate, a
  galaxy reset); how a request finds a queued or just-finished task that
  already answers it instead of starting the same work again; and how
  the API and the queue hand results back and forth with fewer round
  trips. Done: an approved plan, with the build work filed as its own
  items from it.

- [ ] **PERF.21 Generation works with any worker count: the parallel path is built, used and tested (bug)**
  Top priority, and the first thread of phase 0 (Boss, 2026-10-02). Found by
  the debug-mode bug hunt (2026-10-02; report and evidence in the project's shared files under `bug-hunt/`). Boss: "We need a working parallel path, we need to make
  sure the code properly builds and uses the work queue and the worker
  pool to run no matter how many workers it is given. This is the top
  priority as it impacts generation performance significantly." Today
  `src/tests/conftest.py` pins `PLANETGEN_WORKERS=1`. With
  `PLANETGEN_WORKERS=2` or `4` the same 14 generation tests fail every
  time: 6 in `test_gen_resume.py`, 4 in `test_gen_bright_scatter_edges.py`,
  3 in `test_galaxy_gen.py`, and
  `test_fuzz_cli.py::test_galaxy_absurd_radius_is_clean[200.0-...]`
  ("Can't pickle local object ...<lambda>"). The tests' monkeypatches
  never reach the spawned workers, and some tests depend on the seeded
  single-process stream (GEN.39), so nothing checks the path the server
  really uses. Done: every galaxy, sector, scatter and backfill mode runs
  through the work queue and the pool at 1, 2 and N workers; the
  generation tests run at several worker counts (TEST.74); those 14
  tests pass at every count, or are rewritten so their fault injection
  works across processes. Related: TEST.19 (same galaxy at any worker
  count, done), TEST.73, GEN.39.

  - [ ] **PERF.22 On Python 3.12 a run hangs forever when a worker process dies (bug)**
    High. After a worker dies (killed or out of memory),
    `workQueue.WorkQueue._dispatch` calls `ProcessPoolExecutor.submit`. On
    Python 3.12.3 (Ubuntu 24.04's stock Python) that call sometimes blocks
    for good on the executor's `_shutdown_lock`
    (`concurrent/futures/process.py:811`): the run never fails, the lease
    is never freed, and nothing is logged. Repro: `python3.12 -m pytest -p
    no:timeout src/tests/test_work_queue_failures.py::test_a_dead_worker_fails_the_run_and_frees_the_lease`
    hung 4 times in 15 runs; Python 3.11 hung 0 times in 15 (py-spy stack
    in the bug hunt report). Done: a dead worker fails the run cleanly on
    Python 3.11 to 3.13, and the test runs on each.

  - [ ] **TEST.74 Generation tests at more than one worker**
    The test half of PERF.21: a CI leg or a parametrized fixture that runs
    the generation tests at 1, 2 and 4 workers instead of only the
    `PLANETGEN_WORKERS=1` that `conftest.py` pins today, with fault
    injection that reaches the spawned workers. [infra, PERF]

- [ ] **PERF.23 The bright-star progress bar can end at 101% (bug)**
  Found by the debug-mode bug hunt (2026-10-02; report and evidence in the project's shared files under `bug-hunt/`); Boss asked to file it only if it was a real bug and not
  just more stars than expected, and it is real. The bar's total is the
  expected stars per layer and each layer's credit is capped at its own
  share, so extra stars can't push it past 100%. But
  `_LayerTracker.layer_progress` in `generate.py` accepts a report for a
  layer that `layer_done` already credited (reports come through the
  channel-draining thread), so that layer is counted twice: in 3 of 9
  runs with 4 workers under load, a late report arrived and the run ended
  above 100% (for example 101.3% at layer 9); runs with no late report
  ended at exactly 100%. Repro: `generate.py plan --disk-scale-length-pc
  150 --disk-scale-height-pc 20 --bulge-scale-radius-pc 15 --max-ring 40
  --workers 4 --force`, three at once. Done: reports for finished layers
  are ignored, the displayed and `progress.json` percentage is capped at
  100 as Boss asked, and a test feeds a late report.

## DB: Database and schema

DB.1 shipped in 7.35.0 (PR #152).

- [ ] **DB.2 Asteroid field and comet composition rows are written but never read (bug)**
  Found by the Database tests thread (TEST.11, PR #315, 2026-10-01):
  `asteroid_field_composition` and `interstellar_comet_composition` rows
  are saved, but nothing reads them back, so pages show the parent's
  `composition_summary` instead. Decided (Boss, 2026-10-02): keep the
  rows and show them. Done: the asteroid field and interstellar comet
  pages and their API responses read and show the composition rows, and
  a test checks a saved composition comes back on both.

- [ ] **DB.3 resetDb while another process holds id blocks can duplicate primary keys (bug)**
  Found by the Database tests thread (PR #315, 2026-10-01): running
  `resetDb` from one process while another long-lived process (the web
  app or a generation worker) still holds cached id blocks lets the old
  process hand out ids the reset database gives out again, so inserts
  fail on duplicate primary keys. Done: a reset can't lead to reused
  ids, probably by not restarting ids at 1 after a reset (or by making
  holders drop their blocks), with a test that runs both processes.

## API: The JSON API

- [ ] **API.3 Remote generate: generate on a local machine, upload through the API**
  Boss (2026-10-01 19:22Z): "using an API key add a remote-generate
  option where I can use generate.py with the right command lines to
  have my local system (much faster than the server-side) generate
  whatever I want and use the API to add it into the system once it's
  generated. Do not implement yet." Today `generate.py` writes straight
  to a MySQL/MariaDB database it can reach, and the API's write routes
  (`api/edits.py`, API.1's create-a-system route from PR #170) take one
  system or body at a time with an admin's API key. Done: `generate.py`
  gets a remote mode (for example `--remote URL --api-key KEY`) for the
  same commands and options as a local run (`sector`, `galaxy` modes,
  `plan`, the bright-star scatter and backfill, `population`); it asks
  the server for what it needs first (the galaxy's seed and skeleton,
  which sectors are already filled, the id and name state), generates
  on the local machine with all its workers, then uploads the results
  in batches to new API routes that check the key, validate each batch
  (ADM.5's validator), refuse sectors that were filled meanwhile, save
  them the same way a server-side run does (names, ids, bright-star
  levels, caches and tiles invalidated), and report what was added. A
  remote run appears in the job tree and job page (ADM.10, PR #294; ADM.12, PR #285) like
  a server-side one. Boss's answers (2026-10-01 19:28Z):
  - Generate and keep in memory: the local machine needs no database of
    its own.
  - The name indexes and id state are downloaded from the API once, for
    local use, so the run doesn't keep calling the API for them (the
    server reserves id blocks and claims the sectors for the run).
  - Uploads and downloads are compressed (gzip, bz2 or similar): some
    compression time is worth faster transfers.
  - Fully resumable: the local machine caches on disk every API call it
    means to send until the server confirms it; the server keeps a
    buffer of everything it receives, stages it, and writes to the
    database in a controlled way only complete units (a complete star
    system, a complete sector). This probably needs database changes
    so staged data can be stored and flagged incomplete until it is.
  - The local version must match the server's API version (API.4 and
    API.5).
  Boss's answers (2026-10-01 19:32Z):
  - Only an admin's API key can upload; user-level keys can read but
    never upload (API.6).
  - Upload limits are investigated and planned separately (API.7).
  - Reserved ids and claimed sectors stay reserved until the upload
    finishes or an admin clears it; nothing else may use them in the
    meantime. Admins see and clear incomplete uploads on their own page
    (ADM.13).
  - The server checks every upload and has the final say on what is
    written to the database (API.8).
  Order: API.4, API.5 and API.6 first (version and key scopes), then
  API.7's plan, then API.3's sub-items below with API.8 and ADM.13.

  - [ ] **API.9 Key scopes**
    `admin_api_keys` has no scope column today (every key is an admin
    key). Done: a scope column (read, admin, upload), a control schema
    migration; API.6's user-level keys and the upload right use it.

  - [ ] **API.10 Reservations: claimed sectors and id blocks per run**
    `id_blocks` exists (one next id per table, used by parallel
    generation) but nothing reserves sectors. Done: a run record
    holding its claimed sectors and id ranges until it finishes or an
    admin clears it (ADM.13); server-side runs and other uploads skip
    claimed sectors.

  - [ ] **API.11 Staging tables**
    Done: staging storage (a galaxy schema migration) for received
    batches, flagged incomplete until a whole unit (a system, a sector)
    has arrived and passed API.8's checks, then copied into the real
    tables in one transaction.

  - [ ] **API.12 The download: seed, skeleton and name state**
    Done: a route that returns, compressed, what a local run needs
    (galaxy seed, skeleton, filled sectors, name registries, reserved
    id ranges), so `generate.py` needs no database.

  - [ ] **API.13 Generation without a database**
    Today `generate.py` writes through `_db` as it goes. Done: a mode
    where the generators write to an in-memory or on-disk outbox in
    the upload format instead; a client-side cache so an interrupted
    run resumes and resends only what the server hasn't confirmed.

  - [ ] **API.14 Upload routes, compressed, in batches**
    Done: upload routes for gzip batches, with their own body limit
    (separate from the 2 MB `MAX_CONTENT_LENGTH` of PR #293, from
    API.7's plan), answering what was received and verified.

- [ ] **API.4 API compatibility data in the docs**
  Boss (2026-10-01 19:28Z): "let's make API compatibility data and put
  that in the docs as another TODO item to add but not do yet." Done:
  the API has a version of its own, separate from the release number,
  and the docs carry a compatibility table: each API version, the
  release that introduced it, the routes and payload formats it
  changed, and which client versions (`generate.py` remote mode,
  API.3, and the web pages' `lib/apiclient.py`) can talk to which
  server versions. The table is updated in the same PR as any API
  change. Open question: is the API version bumped by hand, or by a
  release note kind like `changes/<name>.api.md`?

- [ ] **API.5 API version and compatibility checking**
  Boss (2026-10-01 19:28Z): "Another TODO item will be API version /
  compatibility checking so that we make sure the server knows how to
  take data from different client versions." Done: every API request
  carries the client's API version (a header); the server answers with
  its own version and the range it accepts (`/api/health`), refuses a
  client outside that range with a clear message saying which version
  to install, and, inside the range, reads older clients' payloads
  through converters so it knows how to take data from each supported
  version. `generate.py`'s remote mode (API.3) checks before it starts
  generating, not after. Built on API.4's compatibility data. Open
  question: how many older versions the server keeps accepting.

- [ ] **API.6 User-level API keys, owned by the account that created them, that can read but not upload**
  Boss (2026-10-01 19:32Z): "Only admin can upload, and Admin can
  create user-level API keys that can access but not upload." Today
  API keys (`admin_api_keys`) belong to admins and carry every admin
  right. Done: an admin can create, name, list and revoke user-level
  API keys; a user-level key can use the read routes (and anything
  later granted to users) but every write route, the remote upload
  routes of API.3 above all, answers 403 for it; an admin key keeps
  every right. Ties in with USR.1 and USR.2 (user accounts and roles)
  and PR #293 (TEST.44), which already answers 403 to
  any API key that makes keys, changes credentials or 2FA, or logs
  out. Decided (Boss, 2026-10-02 01:31Z): "all keys are attached to
  the account that created them". Every API key, admin or user-level,
  records the account that created it, and requests made with it act as
  that account. Today that is an admin account; once user accounts exist
  (USR.1, USR.2) a user's keys belong to that user, so this touches the
  user account work (USR): the key-to-account link for user accounts
  builds after USR.2. Logging every call is its own item, API.15,
  which doesn't wait for user accounts.

- [ ] **API.15 Log every API call with its user, how it came in, and its HTTP response code**
  Boss (2026-10-02 01:31Z): "all api call logs should log the user that
  called for it and if it was called through an api key or someone
  logged into the web site or from the console (console user will just
  be called god). Standard will have all API calls logged and what
  response they generated back (use HTTP responses codes and note in
  documentation that you're using HTTP response codes for parity with
  web logs)." Done: by default every API call is logged with the time,
  the route, the account that made it (the key's owner for an API key,
  the signed-in user for the web site, and "god" for the console), how
  it came in (API key, web session or console), and the HTTP response
  code it got back; the API docs say that HTTP response codes are used
  for parity with the web server's logs; a test makes a call each way
  and checks the log rows. Needs no user accounts: it logs today's admin
  accounts (and "god"), so it is placed in phase 0, ahead of API.6; once
  user accounts exist (USR.1) the logged user can be any account.

- [ ] **API.7 Investigate and plan upload limits**
  Boss (2026-10-01 19:32Z): "upload limits add that as a TODO.md item to
  investigate and plan." Done: a plan, agreed with Boss, for the limits
  on API.3's uploads: the largest request and batch (compressed and
  uncompressed), how many uploads may run at once per key and in all,
  the rate per key, the server's staging-space cap, and what the
  client does when it hits one (wait and retry, shrink the batch).
  Measured against real sectors and the web server's own limits
  (Apache `LimitRequestBody`, IIS `maxAllowedContentLength`, Flask
  `MAX_CONTENT_LENGTH`, set to 2 MB by PR #293) and the API limiter
  (`api/limiter.py`). The plan only; building the limits is a later
  item.

- [ ] **API.8 Verify uploaded data before it is finalized**
  Boss (2026-10-01 19:32Z): "the API will have to have a reliable method
  of syncing and making sure uploaded data is verified good before it
  is finalized in the DB, the server has the final say in how items are
  added to the database." Done: the client and server agree, by batch
  and by unit (a star system, a sector), on what has been sent and
  received, with checksums, so nothing is lost or written twice; every
  staged unit is checked on the server before it is finalized
  (complete, well formed, inside its reserved sector and id block,
  names unique, and the physics checks of ADM.5's validator); the
  server finalizes a unit in one transaction or rejects it with the
  reasons, and may correct what it can (renaming a clashing name,
  re-running derived values) rather than trusting the client's copy.
  Open question: which corrections the server makes on its own, and
  which reject the unit so the client regenerates it.


## ADM: Admin tools

- [ ] **ADM.13 Incomplete uploads page**
  Boss (2026-10-01 19:32Z): "Admin will have to have a page where they
  can see incomplete uploads and clear them but reserved sectors by ID
  cannot be used if an upload isn't finished." Done: an admin-only page
  lists the remote uploads of API.3 that have not finished: who started
  them, when, the sectors and id blocks they reserved, how much is
  staged and verified (API.8), and when the client last sent anything.
  Clearing one, confirmed and written to the admin activity log, throws
  away its staged data and releases its sectors and id blocks; until
  then nothing else (a server-side run, another upload) may use them.
  Part of, or linked from, the job page (ADM.10, PR #294). Open question: should
  an upload with no contact for a long time be flagged as stale on the
  page?

- [ ] **ADM.14 Line up the Generate page's text boxes, not their headings (bug)**
  Boss (2026-10-01 20:20Z): "on the generate screen, line up the text
  boxes not the headings. Text boxes should all be even with each
  other". Today each field on the admin Generate page
  (`web/templates/generate.html`, the `field` macro inside
  `search-fields`) puts its label above its input, and the fields flow
  side by side, so inputs start at different heights and widths
  wherever a label wraps or is longer. Done: down every form on the page
  (New galaxy, Generate sectors and its modes, Plan, Rebuild the bright
  stars, Add a dimmer layer, One-off system), the text boxes share one
  left edge and width and sit level with each other, however long their
  labels are, at desktop and phone widths, in both themes. This includes
  ADM.4's sections (PR #279) and GEN.30's "Bright stars from" field (PR
  #295).

- [ ] **ADM.15 Change the worker count from the Queue page, with a "Ludicrous Speed" mode**
  Boss (2026-10-01 22:11Z): "have admin in the queue menu able to
  change the worker count including a "Ludicrous Speed" that will
  basically change the system to run max CPU power on the system up to
  95% of available CPU power but periodized so that the web interface
  still works (even though they will be slow) and DB calls still work
  (even though they will be slow). This speed mode will not care if
  mysql is running and will have the goal of saturating the host CPU
  safely." Today the worker count is fixed when a run starts
  (`--workers`, `PLANETGEN_WORKERS`, or `worker_count()` in
  `stellarObjects/workQueue.py`: 80% of the cores, one fewer when MySQL
  is on this machine), with workers at lowered priority, and the admin
  Queue page (`admin_queue.html`, ADM.10/11) can't change it. Done: an
  admin can set the worker count on the Queue page (for running and
  future work) and pick "Ludicrous Speed", which aims to saturate the
  host's CPU up to 95% of its capacity without keeping a core back for a
  local MySQL, while still giving the web site and database calls enough
  time that they keep working, if slowly. The change is written to the
  admin activity log.

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

- [ ] **TEST.71 Intermittent failure in the admin planet-regenerate test (bug)**
  `test_admin_edits.py::test_system_page_regenerates_a_planet` fails
  about 1 run in 12 on main (seen by the TODO thread while testing
  GEN.30, PRs #295 and #296, 2026-10-01). Done: the failing case is
  found (loop the test over seeds or runs), the cause is fixed in the
  test or in the code it found, and the test passes on every run tried.
  [infra, ADM]

- [ ] **TEST.72 Intermittent failure in the two-step (2FA) sign-in test (bug)**
  Found by the bug audit (2026-10-01, `bug-audit.md`): the two-step sign-in test failed once in a full run for PR #263
  and passed 3 of 3 times alone. Done: the failing case is found (loop
  it, including near a time-step boundary of the one-time code), the
  cause is fixed in the test or in the code it found, and the test
  passes on every run tried. [infra, SEC]

- [ ] **TEST.73 Intermittent failure in the parallel galaxy-run interrupt test (bug)**
  `test_work_queue_failures.py::test_interrupting_a_parallel_galaxy_run_leaves_no_half_written_sector[2-True]`
  failed once in a full suite run under load and 1 time in 12 targeted
  runs, with a job not in state `cancelled` (reported by the parallel,
  population and navigation tests thread, PR #321, 2026-10-01). It's a
  timing problem in the parallel (two-worker) path; #321 only changed
  the one-worker path. Done: the race is found, the cause is fixed in
  the test or in the work queue's interrupt handling, and the test
  passes on every run tried. [infra, PERF]

- [ ] **TEST.78 A resume test's sector query fails under ONLY_FULL_GROUP_BY on MariaDB 10.11 (bug)**
  Reported by the Database thread (PR #342, 2026-10-02):
  `_sector_rows` in `src/tests/test_gen_resume.py` selects
  `s.ring_index`, `layer_index` and `ring_slot_index` with `GROUP BY
  s.id`, which ONLY_FULL_GROUP_BY rejects on MariaDB 10.11 (error 1055),
  so `test_sectors_a_forced_scatter_skipped_fill_correctly_afterwards`
  fails on main there. Done: the query is valid under ONLY_FULL_GROUP_BY
  on both engines (group by every selected column, or no GROUP BY), and
  the test passes on MariaDB 10.11 and SQLite. [infra, DB]

- [ ] **TEST.79 Route edge cases, written before NAV.12**
  Tests that pin the cases NAV.12 must handle, taken from the
  hop-length study (`nav-hop-length/report.md` in the project's shared
  files): an isolated system; empty and unfilled sectors between the
  endpoints; the galaxy edge and the halo (a lone system 2 kpc above
  the disk); both endpoints in one sector when the best route leaves
  it; and a route graph in separate pieces (the study's 714 islands
  from 2,000 sectors). Done: the tests exist, marked expected-to-fail
  where today's code fails them, so NAV.34 and NAV.12 turn them green.
  [web, NAV]

- [ ] **TEST.76 A bright-star test breaks on Python 3.9 and 3.10 (bug)**
  Found by the debug-mode bug hunt (2026-10-02; report and evidence in the project's shared files under `bug-hunt/`). `_ScriptedRandom(random.Random)` in
  `test_gen_bright_scatter_edges.py` passes a list as its first
  argument; before 3.11 `Random.__new__` seeds with it, so
  `test_a_zero_weight_bin_is_never_picked_by_float_rounding` fails with
  `TypeError: unhashable type: 'list'`. `python_requires` is `>= 3.9` and
  CI has a py3.9 leg. Done: the helper works on Python 3.9 to 3.13.
  [infra, generation]

### Database and migrations

### Generation and the work queue

### Web, API and jobs

### Scripts and ops

- [ ] **TEST.70 Tests for the map JavaScript**
  Today only the pure modules (`galaxystages.js`, `galaxyprisms.js`,
  number and distance formatting) have node tests, run from pytest; the
  stage view, the Sector Map, the System Map and `bookmarks.js` have
  none, and only the accessibility check drives a real browser. Done: a
  Playwright harness (Chromium is already installed for the a11y test)
  that loads each map against fixture data with no database, and tests
  for what MAP.61 will move: picking, hover, keys, Back/Forward and URL
  state, bookmarks and the scale line on the Galaxy Map and the Sector
  Map, written before the refactor so it can't change behaviour
  unnoticed. Builds on TEST.55 to TEST.59 (browser and JavaScript tests,
  PR #307): reuse their harness.

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
    not, what does a broken bookmark show? Per-browser bookmarks
    (`static/bookmarks.js`, MAP.23) shipped in PR #234; storing them in
    the database still needs decision 4 of the drill-down design and a
    migration.

## OPS: Installers, hosting, CI, releases

OPS.1 shipped with the version scheme in `changes/README.md`.

- [ ] **OPS.6 Admin scripts accept impossible `--mysql-port` values (bug)**
  Found while building TEST.60 (PR #288): the admin scripts take
  `--mysql-port 0`, `-1` or `70000` and only fail later with a
  connection error, because `_db.add_mysql_connection_args` doesn't
  range-check the port. Done: every script that takes `--mysql-port`
  rejects anything outside 1 to 65535 with a clear argument error
  before connecting, with a test in `test_admin_script_cli.py`.

- [ ] **OPS.7 Update asks to fill a wiped database with population data (bug)**
  Boss (2026-10-01): "if in the update the user selects to wipe the DB,
  don't then ask to fill it with population data, in fact, remove that
  question entirely from the update." Today step 5 of `update.sh` runs
  `migrate_or_reset_db` and then `offer_population_pass`, so a user who
  just chose to delete the database is asked to run the population pass
  over an empty galaxy; `update.ps1` does the same through
  `Invoke-OptionalPopulation`. Done: neither `update.sh` nor `update.ps1`
  asks about or runs the population pass (the prompt, the `POPULATION=1`
  variable and the `-Population` switch are gone from the update, with
  their usage and header comments), and the closing message says to run
  `generate.py population` by hand when wanted. The installers keep their
  own prompt unless Boss says otherwise.

- [ ] **OPS.8 Update reloads Apache itself when run as root**
  Boss (2026-10-01): "it should just automatically reload apache2 if
  it's running as root." Today `update.sh` ends by printing
  `sudo systemctl reload apache2` (or `restart` after enabling a module)
  for the user to run. Done: when the update runs as root on Linux and
  apache2 is running, it reloads Apache itself (restarts it when a module
  was just enabled) and says so; when not root, or Apache isn't running,
  it prints the command as today. macOS (gunicorn) and Windows stay as
  they are unless Boss asks.

- [ ] **OPS.9 Multi-line messages lose their prefix in the debug log (bug)**
  Low. Found by the debug-mode bug hunt (2026-10-02; report and evidence in the project's shared files under `bug-hunt/`). The sector summary is one log record with embedded
  newlines, so its "Systems:", "Phenomena:" and "Star density:" lines
  land in the debug log with no timestamp, process or source. Done: each
  line is prefixed (or each is logged separately).

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
