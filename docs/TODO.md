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

DOC has no open items today.

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
rebuilt as phases 0 to 3+ from the dependency report, and on
2026-10-07 rebuilt again from Boss's lists of 2026-10-03 and
2026-10-07: every bug in phase 0, which has two lanes (bugfixes and
groundwork). Every open item below is in exactly one phase, and every
item's prerequisites are in its own phase or an earlier one. Each phase
file gives the goal, the build threads with their order and
prerequisites, and the open questions;
[plan/notes.md](plan/notes.md) keeps the research notes, the files
several items share and the judgment calls. This file keeps each
item's full text. A new item goes into a phase's table in the same PR
that files it.

| Phase | Plan | Goal | Items |
|---|---|---|---|
| 0 | [phase-0-roots.md](plan/phase-0-roots.md) | Two lanes (Boss 2026-10-03 05:38Z: "the most fundamental and needed changes first in phase 0 along with the bug fixes in 2 lanes, bugfixes and groundwork"). Bugfixes: the CI failures first, then generation, console and progress, maps and pages, prevalence (GEN.48), and ops and test flakes. Groundwork: the package layout and the move to third-party libraries with Redis (Boss's explicit directive), the data model (SQLAlchemy and Alembic, values in columns not JSON, one point-in-space object, Pydantic, scipy and astropy), names from IDs, the RQ queue with streamed logs and progress, Shoelace and TanStack components, the shared map engine, nebula shapes, and the UX sweep. Bugs that the groundwork fixes are folded into it and listed under it. GEN.65 is done (PR #476); GEN.117, GEN.116 keeps watch for Boss's error text. | GEN.116, TEST.92, TEST.93, TEST.94, SEC.29, TEST.72, SEC.30, TEST.83, UX.39, OPS.20, DB.11, DB.13, ADM.21, GEN.66, GEN.74, GEN.68, GEN.69, GEN.70, GEN.71, GEN.72, GEN.73, GEN.67, PERF.24, OPS.19, ADM.22, ADM.24, ADM.25, ADM.26, UX.40, UX.26, UX.31, UX.27, UX.2, ADM.14, UX.41, UX.33, ADM.34, NAV.7, MAP.95, NAV.13, NAV.14, MAP.65, MAP.79, NAV.15, NAV.29, NAV.33, MAP.66, MAP.67, MAP.68, MAP.61, NAV.32, MAP.102, MAP.107, MAP.106, MAP.108, MAP.109, MAP.110, MAP.111, MAP.112, MAP.113, MAP.116, NAV.46, GEN.75, MAP.103, MAP.104, MAP.105, UX.38, UX.37, UX.21 |
| 1 | [phase-1-built-on-roots.md](plan/phase-1-built-on-roots.md) | Object references and the picker finished, the database check and repair, resumable runs, the reworked Generate page (spans, radial fills, random neighborhoods, directives, show on the map), galaxy generation changes (phenomena placed galaxy-wide first, fill order, backfill from the run's edge, nebula volume backfill, star-type tiers), the habitability index, tech levels and facility types, spin and orbital-update thresholds with their limits, routing with unknown-space stops and charting along a course, the nearby search, Select mode and view filters on the maps, bookmarks, admin control from every screen, in-universe wording, the API log and scopes, and reproducible galaxies up to the golden-seed test. | NAV.8, NAV.9, NAV.16, NAV.3, DB.8, DB.9, PERF.29, PERF.30, ADM.29, ADM.30, ADM.31, GEN.96, GEN.97, ADM.28, GEN.24, GEN.41, GEN.98, GEN.100, GEN.99, GEN.101, GEN.102, GEN.103, GEN.84, GEN.85, GEN.86, GEN.87, GEN.88, GEN.89, GEN.83, POP.8, POP.9, POP.7, POP.10, GEN.94, GEN.104, GEN.106, GEN.107, GEN.108, NAV.10, NAV.12, UX.35, NAV.11, NAV.42, NAV.47, NAV.48, NAV.43, NAV.44, MAP.119, MAP.122, MAP.123, MAP.124, MAP.89, UX.46, UX.47, ADM.32, MAP.120, ADM.35, UX.23, UX.22, UX.3, UX.42, ADM.15, API.15, API.4, API.7, API.9, GEN.56, GEN.57, DB.7, GEN.58, TEST.77, OPS.8, OPS.13, OPS.14, ADM.18, GEN.59 |
| 2 | [phase-2-maps-picker-backfill.md](plan/phase-2-maps-picker-backfill.md) | The planet class refactor around the habitability index (with GEN.33, GEN.28, GEN.27 and GEN.29), nebula conditions on planets, n-body orbital updates and rogue collisions, editable trajectories, courses and waypoints, generate-by-recipe in the API, the backfill density pass, daily maintenance, the pilot-style visual design, light-travel positions, the asteroid-field and anomaly plans, and the API pieces remote generation needs first. | GEN.33, GEN.28, GEN.27, GEN.91, GEN.92, GEN.29, GEN.90, MAP.58, MAP.75, MAP.59, NAV.20, NAV.21, NAV.17, NAV.18, NAV.4, NAV.24, NAV.36, NAV.39, NAV.49, UX.32, UX.30, UX.43, GEN.42, GEN.43, PERF.18, GEN.40, PERF.20, API.5, API.10, API.11, API.12, API.16, MAP.69, MAP.70, ADM.17, OPS.15, GEN.61, OPS.18, OPS.16, OPS.17, ADM.19, NAV.45, UX.48, UX.45, GEN.115, GEN.109, ADM.36, GEN.110, GEN.105, GEN.95, GEN.93, GEN.112, GEN.113, MAP.121, API.18, VIEW.5 |
| 3 | [phase-3-engine-3d-remote.md](plan/phase-3-engine-3d-remote.md) | The 3D system view and infinite zoom from galaxy to moon on one interface, orbital trajectories in their own frame, courses that bend around gravity wells, remote generation through the API (with galaxy-scale recipes) reproducing what the server would make, repair that reads the newest settings JSON, and the anomalies chosen in phase 2. | MAP.71, MAP.72, MAP.73, MAP.74, MAP.62, MAP.126, NAV.22, NAV.23, NAV.5, NAV.25, NAV.26, NAV.27, NAV.28, NAV.6, API.13, API.14, API.8, ADM.13, API.3, API.17, DB.10, ADM.20, GEN.114, MAP.125, API.19 |
| 3+ | [phase-3plus-accounts-sky-galaxies.md](plan/phase-3plus-accounts-sky-galaxies.md) | The open-ended tail: user accounts (with API.6 keys, saved courses and Hill-radius emails), the view of the sky from a planet, the plan for more galaxies, and the end state of reproducible galaxies (`generate.py reproduce`). | USR.2, API.6, USR.3, USR.4, USR.5, USR.6, USR.7, USR.8, USR.1, NAV.19, GEN.111, VIEW.1, VIEW.4, VIEW.2, VIEW.3, GEN.9, OPS.12, GEN.55 |

Phases overlap: a phase's later threads can start while the next
phase's first ones run, as long as the order inside each phase holds.

Boss's list of 2026-10-01 23:53Z (`new todos.txt`, with research notes;
the files are in the project's shared files under `todo-tasks/research/`)
became these items:

| # | Boss's item | ID |
|---|---|---|
| 1 | Store every sector's backfill level (-1, lowest L_sun, 0) | GEN.44 (done, PR #425) |
| 2 | Rogue planets after systems and phenomena; expanded rows span the table | UX.24 (done, PR #487) |
| 3 | Rogue planet gas giant vs terrestrial probability | GEN.45 |
| 4 | Rogue planet octant and a map symbol link | UX.25 (done, PR #490) |
| 5 | Edit and admin buttons as a button menu | UX.26 |
| 6 | Rogue planets dim, barely noticeable | MAP.82 |
| 7 | "Mark Rogue Planets" keeps its highlight while on | MAP.83 |
| 8 | Rogue planets smaller unmarked, bigger and clickable marked | MAP.84 |
| 9 | System and nav buttons in one row, nav folding into a menu | UX.27 |
| 10 | Icons instead of words on buttons (investigate) | UX.28 (done, PR #490) |
| 11 | Star names at most two words | GEN.46 (done, PR #370) |
| 12 | Every comet's type shown with a link | UX.29 (done, PR #487) |
| 13 | Planet information without the Markdown render | UX.30 |
| 14 | System edit as a quick menu, not a long panel | UX.31 |
| 15 | Planet rows show class only | UX.32 |
| 16 | Filter phenomena by type and class | UX.33 |
| 17 | No nebulae being created | GEN.47 (done, PR #419) |
| 18 | 3D galaxy, arc pick, no sector lines (other galaxy-map items edited to match) | MAP.85 |
| 19 | Filled sectors translucent, colored by their stars; blocks averaged | MAP.86 (done, PR #432) |

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
  Prerequisite: UX.40.
  Plan (2026-10-07): Folded into UX.40: the menus become Shoelace
  dropdowns.

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
  Prerequisite: ADM.22.
  Plan (2026-10-07): The banner reads the RQ job's published progress
  (ADM.22), not progress.json.

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
  The nebula "-" no-op is split out as UX.38 (phase 0).
  Order (pre-planning thread): run it after MAP.55 and MAP.60, which
  already remove some dead controls on the Galaxy Map.
  It also runs after UX.37 (Boss's sweep for redundant and duplicate
  controls), so this layout pass checks the controls that are kept.
  Prerequisites: MAP.68, UX.26, UX.27, UX.31, NAV.32, UX.37.
  Plan (2026-10-07): Moved into phase 0 (all bugs in phase 0), as the
  last item of the groundwork lane.

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
    Prerequisite: GEN.66.
    Plan (2026-10-07): The ladder is built on astropy.units (GEN.66).

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
  Prerequisite: UX.40.
  Plan (2026-10-07): Folded into UX.40, and the first screen of ADM.34
  (one admin menu per screen).

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
    Prerequisite: UX.26.
    Plan (2026-10-07): Folded into UX.40; the system page's case of
    ADM.34.

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
  Prerequisite: UX.40.
  Plan (2026-10-07): Folded into UX.40.

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
  The filter options are built from each kind's class list, so classes
  GEN.28 adds later appear on their own; GEN.47 (nebulae that exist)
  is done (PR #419), so nothing blocks it.
  Prerequisite: UX.41.
  Plan (2026-10-07): Folded into UX.41: the phenomena list gets faceted
  filters there.

- [ ] **UX.35 The NAV page route shown horizontally, wrapping onto several lines on narrow screens**
  Boss (2026-10-02 01:53Z): "display the path horizontally and find a
  way to split it between multiple lines for mobile or limited
  displays." Today the route is a vertical list (`<ol class="nav-route">`
  in `nav.html`). Done: the stops run left to right, each a link, with
  the hop distance between them; on phones and narrow panels it wraps
  onto several lines (a container query, not a device check), never
  splitting a stop across lines; screen readers still get an ordered
  list. Runs alongside NAV.12; NAV.36 styles its unknown-space hops.
  Design: [docs/design/course-routing.md](design/course-routing.md)

- [ ] **UX.37 A UX sweep: remove redundant and duplicate controls so the interface gets out of the way**
  Boss (2026-10-02 04:42Z): "Add a UX sweep TODO item, we want to sweep
  the user interface and get rid of redundancies and duplicate buttons
  and the like, the main goal is for the Interface to get out of the way
  but be there when we need it." The principle: the interface gets out
  of the way but is there when it is needed. Done, in two steps. First,
  an audit of every public and admin page and every map (Galaxy, Sector,
  System, phenomenon, NAV) at phone and desktop widths, listing each
  duplicate button or control (the same action offered twice on one
  screen, such as a map button that repeats a Menu or breadcrumb entry),
  each redundant link or panel (the same information or route shown
  twice), and each control that could move into a menu or appear only
  when it applies (selection, pick mode, admin). The audit is a list
  with a proposal per entry (remove, merge, move into a menu, show only
  when relevant), filed in the project's shared files; Boss reviews it
  and nothing is removed until he has. Second, the removals and merges
  he approves, with the hint text and docs that name them, and a browser
  test per page that the kept controls still work. Related items, linked
  rather than repeated: UX.21 (overlapping buttons and dead controls) is
  the final layout pass and runs after this sweep; UX.28 (icons) and
  UX.2 (menu sizes) set how the kept controls look; MAP.55 (done, PR
  #369) already folded the Galaxy Map's buttons into one Menu. Why phase
  2: the sweep audits controls that are still changing, so it waits for
  the pages and maps they live on to settle: UX.28's icons, the Galaxy
  Map breadcrumb and history buttons (MAP.93 and MAP.94 done in PR #399, MAP.95), the NAV
  page layout (NAV.41), the slab and segment pick (MAP.56), and the page
  action menus (UX.26, UX.27, UX.31). Prerequisites: MAP.95, UX.26, UX.27, UX.31, UX.40, ADM.34.
  Plan (2026-10-07): Moved into phase 0 (groundwork lane, last) because
  UX.21 is a bug and needs it.

- [ ] **UX.38 The nebula and remnant diagrams' "-" button does nothing at the 1 ly limit (bug)**
  Split out of UX.21 (replan, 2026-10-02), its one known dead control:
  on the nebula and supernova remnant diagrams, anything about half a
  light-year across or larger opens already at the 1 ly zoom-out limit
  (`lib/phenomenonmap.py` lines 53 and 127, the clamp in
  `static/mapzoom.js`), so "-" has nowhere to go. TEST.55 pins it with a
  strict xfail in `test_web_browser_maps.py`. Done: the diagram opens
  with room to zoom out (or the "-" button is disabled at the limit,
  with its hint text), the xfail is dropped and the test passes.
  Prerequisite: MAP.105.
  Plan (2026-10-07): Folded into MAP.105: the nebula view becomes the 3D
  shape, which replaces this diagram.

- [ ] **UX.39 Markdown rendered by the markdown library**
  Today `mdconvert.py` converts the system and wiki Markdown by hand.
  Done: the `markdown` package renders it with the same output for the
  pages and wiki export (a test compares a sample), and `mdconvert.py`
  is deleted.
  Design: [docs/design/library-migration.md](design/library-migration.md)

- [ ] **UX.40 Buttons, menus and dialogs from Shoelace web components**
  From "Web UX Development Notes.md": buttons, menus, dropdowns, dialogs
  and form fields use Shoelace components, vendored as ES modules (no
  CDN at runtime). Done: a shared component set in the base template;
  the menus size to their contents (UX.2), admin and edit actions open
  from a button menu (UX.26, UX.31), the system page buttons sit on one
  row (UX.27), and the Generate page text boxes line up (ADM.14), each
  closed in this item's PRs; both themes and keyboard use work.
  Design: [docs/design/library-migration.md](design/library-migration.md)

- [ ] **UX.41 Tables on TanStack Table and TanStack Virtual**
  Today every list page is a server-rendered table paged 50 rows at a
  time (`tabledisplay.py`). Done: list pages use TanStack Table with
  virtual scrolling, sorting and faceted filters fed by the API, keeping
  50-row pages for the server; the phenomena page filters by class and
  type (UX.33) through it.
  Design: [docs/design/library-migration.md](design/library-migration.md)

- [ ] **UX.42 In-universe wording across the interface**
  Boss (2026-10-03 05:38Z): "Begin replacing wording to make the
  interface in-universe (i.e. the interface and wording should not be
  designed so they sound like they know it's for a game)." Done: a word
  list (for example "Generate" becomes "Chart", "Generated only" became
  "Charted only" in MAP item MAP.111) agreed with Boss, then pages,
  buttons and messages changed page by page; admin-only tools may keep
  plain wording.
  Prerequisite: UX.37.

- [ ] **UX.43 A visual design built like a pilot's starmap and navigation console**
  Boss (2026-10-07 11:47Z): "All design elements should be built from
  the ground up to look like, as much as possible, a starmap and
  navigational system that a space pilot might use." Done: a design
  brief (colours, type, panel shapes, iconography) approved by Boss,
  then the base template and shared components restyled to it in both
  themes, with contrast and reduced motion kept.
  Prerequisites: UX.42, UX.37.

- [ ] **UX.45 Bookmark management**
  Boss (2026-10-07 11:47Z): "Bookmark management system." Today
  bookmarks are kept per browser by `bookmarks.js` (NAV.40). Done when
  its subitems are.
  Prerequisites: UX.46, UX.47, UX.48.

  - [ ] **UX.46 Wiping the galaxy wipes the bookmarks**
    Boss (2026-10-07 11:47Z): "When a galaxy is wiped all bookmarks are
    wiped too." Done: bookmarks are stamped with the galaxy seed, and a
    page load that finds bookmarks for another seed (after a reset or
    new galaxy) drops them; a test covers it.

  - [ ] **UX.47 A bookmark manager: list, rename, sort, group and delete**
    Done: a Bookmarks page listing every bookmark with its kind and
    place, where each can be renamed, grouped, sorted, opened on its map
    and deleted, and the list can be exported and imported as a file.
    Prerequisite: UX.46.

  - [ ] **UX.48 Charted regions as bookmarks that frame and outline the region**
    Boss (2026-10-07 11:47Z): "Clusters of contiguous generated content
    gets added to a system bookmark which will zoom in as close as it
    can to see all sectors in the group and highlight the shape of the
    region specifically." Done: every contiguous group of charted
    sectors gets a system bookmark (kept up to date as sectors are
    charted) that opens the Galaxy Map fitted to the group with its
    outline highlighted.
    Prerequisites: UX.47, MAP.65.

## MAP: Galaxy Map, Sector Map, System Map

MAP.2 with MAP.22 and MAP.23, MAP.15 and MAP.30 shipped in PR #234.

- [ ] **MAP.95 A "Forward to current" button next to the map's Back and Forward**
  Boss (2026-10-02 04:29Z, with MAP.93, done in PR #399): "we'll keep track of back and
  forth so we can always undo our last zoom, and we'll use that for the
  forward if we just went back we get to go back forward again and that
  also need a 'forward to current' button." Checked on main: the map's
  Back and Forward (MAP.26) already walk its own history (`mapIndex`
  and `maxIndex` on browser history entries in
  `static/galaxystageview.js`), so undoing a zoom and redoing it work
  today. Done: a "Forward to current" button beside Forward jumps
  straight to the newest view in that history (`maxIndex`), disabled
  when already there, on desktop and in MAP.94's phone layout. A
  browser test goes in three levels, back two, and forward to current.

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
  Builds on MAP.88 (fit the whole drawn system inside the frame; done,
  PR #405) in the same `lib/systemmap.py`, and keeps that fit working. Decided (Boss, 2026-10-02 01:01Z: "I agree, we'll go
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
  Design: [docs/design/galaxy-drilldown-navigation.md](design/galaxy-drilldown-navigation.md), section 15

- [ ] **MAP.59 Make it plain that a zoomed-in slab is a slab, not a wedge**
  Boss (2026-10-01 20:45Z): "We need to make it clearer, when we've
  zoomed into a specific slab, that we're viewing a specific slab and
  not a wedge. I'm not sure how to do that so do some research on that
  and then add the to-do items to make it happen."
  Arc pick (MAP.85, Boss 2026-10-01 23:53Z): the ghost is the rest of
  the picked arc.
  Its ghost keeps MAP.77's rule (slab outlines only; MAP.77 done, PR
  #410). Prerequisite: MAP.75.

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
  Design: [docs/design/galaxy-drilldown-navigation.md](design/galaxy-drilldown-navigation.md), section 15

  - [ ] **MAP.75 The mini map as a second engine view**
    The mini map is a second, locked camera on the same scene data;
    built on MAP.61's controller with a "locked" policy rather than its
    own renderer (one WebGL context, scissored like the System Map's
    sphere overlay).

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
  Prerequisites: MAP.65, MAP.66, MAP.67, MAP.68.
  Plan (2026-10-07): Moved into phase 0 with its subitems. MAP.125
  (infinite zoom) is its end state in phase 3.

  - [ ] **MAP.65 One picking, hover and info-panel layer**
    Done: one module for raycast and screen-space picking, the hover
    highlight and tooltip (the Sector Map has none today) and the info
    panel (fields, Nav from/to, Use as destination, Generate buttons,
    bookmark ☆), fed by each view's objects. The Sector Map's info
    panel gains the ☆ the drill-down design left for later.
    Plan (2026-10-07): Moved into phase 0: the shared picking layer is
    what fixes MAP.108, MAP.107, MAP.112 and NAV.46.

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
    Prerequisite: MAP.69.
    Plan (2026-10-07): Positions go through the point-in-space object
    (GEN.74).

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
  Prerequisite: MAP.65.
  Plan (2026-10-07): Moved into phase 0 with the engine. Its per-kind
  toggles include nebulae, which closes half of MAP.113; MAP.123 extends
  them to star types and the Galaxy Map.

- [ ] **MAP.102 Galaxy Map streaming with a BVH and 3D tiles, and camera-relative rendering**
  From "Web UX Development Notes.md": the Galaxy Map streams its stars
  as 3D tiles with level of detail, picks through a BVH (three-mesh-bvh)
  instead of scanning points, and draws camera-relative so far-out
  coordinates don't jitter (`galaxymap3d.js`, `galaxyprisms.js`). Done:
  zooming loads only the tiles in view at the right detail, picking
  stays fast in dense sectors, and the existing map tests pass.
  Design: [docs/design/library-migration.md](design/library-migration.md)

- [ ] **MAP.103 Nebulae don't show on the Galaxy Map or any other map (bug)**
  Boss (2026-10-03 05:38Z): "I cannot see nebula on the galaxy map or on
  any other map." GEN.47 (PR #419) made nebulae exist. Done: nebulae are
  drawn on the Galaxy Map, the Sector Map and the System Map from their
  shape mesh, and a browser test finds one.
  Prerequisite: GEN.75.

- [ ] **MAP.104 Nebula shading is missing on unfilled sectors (bug)**
  Boss (2026-10-03 05:38Z): "Nebula shading should show up on unknown /
  unfilled sectors as in the galactic map as well as for filled
  sectors." Done: a nebula placed galaxy-wide shades the Galaxy Map over
  unfilled sectors as well as filled ones.
  Prerequisite: MAP.103.

- [ ] **MAP.105 The nebula view should show its whole shape with the dimmed galaxy around it (bug)**
  Boss (2026-10-03 05:38Z): "Nebula phenomena view should show the total
  shape of the nebula's area as a oblong region with the larger galaxy
  map around it but dimmed (just there for reference with the brightest
  stars showing)." Done: the nebula page shows its mesh in the 3D view
  with the surrounding galaxy dimmed and only the brightest stars drawn.
  Prerequisite: GEN.75.

- [ ] **MAP.106 The breadcrumb trail falls out of sync with the map (bug)**
  Boss (2026-10-03 05:38Z): "Breadcrumb trail should always be in sync
  with the actual map interface and view." Done: the breadcrumb is drawn
  from the map's current state on every change (zoom, pick, back and
  forward, URL load), and a browser test walks the levels and checks
  them.
  Prerequisites: NAV.14, MAP.67, MAP.107.

  - [ ] **MAP.107 Selecting an empty slab near the core says "There is no layer x here" (bug)**
    Boss (2026-10-03 05:38Z): "When selecting an empty slab near the
    galactic core results in an error where it says "There is no layer x
    here" (x is whatever layer it is).  This happens because it gets out
    of sync, you select the slap, or try to, and nothing happens except
    it moves the breadcrumb trail along." Done: picking any slab, empty
    or not, moves the map and breadcrumb together with no error.
    Prerequisite: MAP.65.

- [ ] **MAP.108 Empty slabs and wedges near the core can't be selected, and the side buttons block clicks (bug)**
  Boss (2026-10-03 05:38Z): "Galacit view, some slabs are unselectable
  near the core when zooming in.  The buttons to the side block being
  clicked on and the map interface does not respond to those clicks in
  that region. The system seems unable to allow me to select any slab or
  wedge in which nothing exists. I should be able to navigate to any
  space freely because I might be selecting a sector or space to fill."
  Done: every slab, wedge and block can be picked whether or not
  anything is in it, and the slab buttons never cover the map's pick
  area.
  Prerequisite: MAP.65.

- [ ] **MAP.109 Zooming in and out loads slowly (bug)**
  Boss (2026-10-03 05:38Z): "Loading between zooms in and out gets very,
  very sluggish." Boss (2026-10-07 16:26Z), the first of three
  problems: "the galaxy map doesn't update very fast because of all the
  data it's receiving from the server.  We're dealing with sometimes
  1000's or 10's of 1000's or 100's or 1000's of stars." Done: the
  server sends only the tiles in view at the detail the zoom needs
  (MAP.102), the payload per view has a stated cap, and a measured zoom
  in and out across levels stays under a stated time on a full test
  galaxy.
  Prerequisite: MAP.102.

- [ ] **MAP.110 Slab button lines come out of numerical order (bug)**
  Boss (2026-10-03 05:38Z): "Slab button lines don't always appear in
  the button order, sometimes appearing out of numerical order." Done:
  slab buttons and their leader lines are always in slab-number order
  (`galaxystageview.js`), with a test.
  Progress (2026-10-07): PR #491 fixed the crossing lines only; the
  buttons still follow on-screen order, which can differ from slab
  number order after the view turns.

- [ ] **MAP.111 "Generated only" should be "Charted only" and dim the stars too (bug)**
  Boss (2026-10-03 05:38Z): "In the full galaxy view when I select
  "generated only" (change the name to charted only), it should dim the
  stars as well and make it blatantly obvious where all the generated
  sectors are within each wedge." Today the toggle (`galaxymap3d.py`
  line 581, `galaxyblocks.js`) dims blocks only. Done: renamed "Charted
  only", it dims stars outside charted sectors as well and outlines the
  charted sectors in each wedge.

- [ ] **MAP.112 Nothing can be selected while "Generated only" is on (bug)**
  Boss (2026-10-03 05:38Z): "Cannot select things from the galacitc map
  while the only generated filter is on." Done: picking works the same
  with the filter on or off.
  Prerequisites: MAP.65, MAP.111.

- [ ] **MAP.113 A nebula covering the whole sector can't be unselected, and nebulae need a show/hide toggle (bug)**
  Boss (2026-10-03 05:38Z): "When I select a nebula either on purpose of
  on accident in the sector map, I cannot unselect it if it covers the
  entire sector space, as if often does.  Nebula should also have a
  toggle button to turn on or off so if they get in the way we can turn
  them off from the display." Done: Escape or a click on the selected
  nebula clears it, and the Sector Map's nebula toggle hides them.
  Prerequisite: MAP.79.

- [ ] **MAP.116 Dense generated sectors crowd the Galaxy Map when zoomed out (bug)**
  Boss (2026-10-03 05:38Z): "When sectors are filled in and we're in the
  galaxy view, we have to be more intentional and strict with how many
  stars to show from those sectors as we zoom out. Recent test run
  indicates that for dense sectors it becomes unreadable at about 2-3
  zoom levels outside of the sector zoom in the galactic view. Fix with
  an algorithm that is more choosey about what to show even if we have
  to allow a fudge factor for when it's a very dense sector vs. when
  it's a very sparse sector.  Sparse sectors have more room to show more
  and be usable." Boss (2026-10-07 16:26Z), the second of three
  problems: "When viewing and zooming in, the fine tuning between what
  star brightnesses should be displayed at a particular zoom level to be
  useful.  Comets, neutron stars, etc, do get displayed on the galaxy
  map if sufficiently zoomed in." Folds MAP.115 (Boss 2026-10-03
  05:38Z: "Comets don't need to show up in a sector map until we get to
  the full sector, don't show comets, rogue planets, or asteroid fields
  from the galactic map."). Done: each zoom level has a rule for which
  stars and objects it shows: a screen-density budget that takes the
  brightest stars first, with a per-sector allowance that gives sparse
  sectors more room; a brightness floor per level; remnants such as
  neutron stars only where the budget reaches them; and comets, rogue
  planets and asteroid fields never on the Galaxy Map, only on the
  Sector Map. The rules are one table in the code, and a test checks
  each level's counts on a dense and a sparse sector. It must also fix
  GEN.117: brightest-first picking per tile keeps only the young blue
  stars on the plane, so the budget must keep stars spread across
  height (a share reserved by layer, or picks by luminosity band).
  Prerequisite: MAP.102.

- [ ] **MAP.119 Expected star density editable by admins on the Galaxy Map**
  Boss (2026-10-03 05:38Z): "In the galaxy view, the expected star
  density should be an editable field for admin and a static field for
  users/guests.  Once edited it changes it for the current view (i.e. an
  entire slab or entire block)." Done: the info panel shows expected
  density; admins can edit it for the slab or block in view, stored as
  an override that later fills use and logged.
  Prerequisite: MAP.65.

- [ ] **MAP.120 Bright-star backfill from the Galaxy Map's block, slab and wedge menus**
  Boss (2026-10-03 05:38Z): "Add bright star backfill to block menus in
  galactic view." Done: the admin menu on a block, slab or wedge runs
  the bright-star backfill for it as a job.
  Prerequisite: MAP.65.

- [ ] **MAP.121 Every map shows and steps to its neighbouring regions, on one map engine**
  Boss (2026-10-03 05:38Z): "In the sector level map, we should show
  adjacent sectors dimmed around it, showing all surrounding sector
  interiors dimmed and minus comets and rogue planets for anything but
  the sector in question.  Implement so that the sectors between the
  camera eye and the sector we are viewing become invisible when they
  would block the view in a fade, the more they obstruct the view the
  more invisible they become." Boss (2026-10-07 16:26Z), on the Galaxy
  Map: "you cannot see what is in adjacent blocks or slabs nor can you
  move to adjacent regions without zooming back and forth, which is
  combersome." and (16:27Z): "edit MAP.121 to include the new item as
  we want to have 1 mapping engine for sectors, galaxy, and systems
  eventually". Folds MAP.127. The goal is one map engine shared by the
  Galaxy, Sector and System Maps. Done: on every map and at every
  Galaxy Map drill-down stage, the neighbouring regions (sectors,
  blocks, and the slabs above and below) are drawn dimmed around the
  current one, without comets or rogue planets and within MAP.116's
  budget; regions between the camera and the one in view fade the
  more they block it; arrow buttons and keys step to a neighbour at
  the same level without zooming out, updating the URL and breadcrumb;
  the code for this is one shared engine piece, not one per map; and
  a browser test steps sideways and up and back on each map.
  Prerequisites: MAP.66, MAP.102.

- [ ] **MAP.122 A Select mode on every galaxy view: Galaxy (blocks and sectors) or Star**
  Boss (2026-10-07 11:47Z): "Button in Galaxy display to allow selecting
  a star / star system or other object that is normally blocked from
  being selected in the galaxy interface across all interfaces." And:
  "At all galactic views we have a 'Select Mode' it defaults to Galaxy
  (wedges, blocks, slabs, etc like we have right now), but one can have
  a Star Select mode where the stars become selectable, when a star is
  selected you can bookmark it or set a waypoint." Done: a Select mode
  control on every galaxy view (default Galaxy); in Star mode stars and
  other objects pick, hover and open their info panel with Bookmark and
  Waypoint actions.
  Prerequisite: MAP.65.

- [ ] **MAP.123 Show or hide star types and phenomena, and set the luminosity floor, on the Galaxy and Sector Maps**
  Boss (2026-10-07 11:47Z): "Sector display can hide systems by star
  type or phenomena, same in galaxy views, users can select to hide or
  show classes of stars, adjust the luminosity at which stars appear on
  the map, and highlight phenomena or hide them in the galactic view."
  Done: both maps have a filter panel for star classes and phenomenon
  kinds (show, hide, highlight) and a luminosity slider, kept in the
  URL.
  Prerequisites: MAP.65, MAP.79.

- [ ] **MAP.124 The Galaxy Map opens zoomed to fit all charted space**
  Boss (2026-10-07 11:47Z): "Galactic map should automatically zoom in
  to the nearest zoom possible to see as much of the charted space as
  possible." Done: with no view in the URL, the map opens at the closest
  zoom that shows every charted sector.
  Prerequisite: MAP.67.

- [ ] **MAP.125 Infinite zoom: one 3D interface from the galaxy down to a moon**
  Boss (2026-10-07 11:47Z): "Infinite Zoom!  Able to use the same
  interface to go from galaxy all the way down to a moon by zooming in
  and out.  Make everything the same 3D interface." Done: one view zooms
  continuously from the galaxy through block, sector and system to a
  planet and its moons, with the frames switched by the position object
  (GEN.74).
  Prerequisites: MAP.61, MAP.62.

- [ ] **MAP.126 Show the orbital trajectories of selected objects in their frame of reference**
  Boss (2026-10-07 11:47Z): "Show orbital trajectories for selected
  objects in their frame of reference." Done: selecting a body draws its
  path (orbit, or galactic trajectory) in the frame of what it moves
  around.
  Prerequisites: MAP.70, MAP.73.

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
  Prerequisites: NAV.13, NAV.14, NAV.15, NAV.16, NAV.32.
  Plan (2026-10-07): Moved from phase 3 to phase 1: its subitems are now
  in phases 0 and 1.

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
    It keeps MAP.93's one-line collapse and MAP.94's phone layout (PR #399) on
    every page (Boss 04:29Z: "the breadcrumbs should never be more
    than a single line").
    Prerequisite: NAV.13.
    Plan (2026-10-07): Moved into phase 0; the one breadcrumb fixes
    MAP.106.

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
    Prerequisite: NAV.7.
    Plan (2026-10-07): Moved from phase 2 to phase 1 so the picker
    parent (NAV.3) closes there.

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
  Prerequisites: NAV.20, NAV.21, NAV.22, NAV.23.
  Plan (2026-10-07): Boss (2026-10-07): "Plotted courses should appear
  on the galactic map and stay until cleared." See NAV.49.

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
  Plan (2026-10-07): Moved into phase 0 (groundwork lane): the map
  engine items that fix the 2026-10-03 map bugs need it.

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
  systems. NAV.34 (joining the graph's islands) is done (PR #427), and NAV.12
  (a route always exists, no hop limit) is built with it.
  Design: [docs/design/course-routing.md](design/course-routing.md)
  Prerequisite: DB.11.
  Plan (2026-10-07): The position indexes are an Alembic migration
  (DB.11).

  - [ ] **NAV.11 Travel times for the system-to-system route too**
    Today warp and fold times are shown only for the direct distance;
    the route shows only its length. Boss (2026-10-02 04:19Z, with
    NAV.41): "it should calculate the travel time using that route
    assuming each planet gets stopped at." Done: the route gets the
    same warp and fold tables, per hop and in total, the total being
    the sum of the hops with a stop at every system on the route (each
    hop timed from rest to rest at the chosen warp or fold factor, the
    constant speeds `warp_speed_c` and `fold_speed_c` of
    `navigation-frames.md`; no acceleration model exists, and no ship
    range since NAV.37 was dropped). Default taken: no time spent at a
    stop. Open question: should each stop add a fixed stay, and does
    "each planet" mean only the systems on the route (the default) or
    a visit to every planet inside each system?
    Design: [docs/design/course-routing.md](design/course-routing.md)

  - [ ] **NAV.42 Each route stop shows the course and distance to the next stop**
    Boss (2026-10-02 04:19Z, with NAV.41): "each stop has the course and
    distance to the next stop in xxx mark yyy zzz distance". Done: every
    stop in the route (UX.35's strip, and the course map's tooltip)
    shows the course to the next stop in the existing notation from
    `navigation.format_course`, three-digit bearing, "mark", three-digit
    mark, then the hop's distance on the site's distance ladder, for
    example `045 mark 012, 3.2 ly`; each hop's course is worked out with
    `course_between` in the frame NAV uses for that pair (the Sector
    Local Frame when both stops share a sector, the Galactic Frame
    otherwise, as `navigation-frames.md` sets out), and the last stop
    shows none. Unknown-space jumps (NAV.36) show theirs the same way.
    A test checks the bearing, mark and distance of a known hop.
    Prerequisite: UX.35.
    Design: [docs/design/navigation-frames.md](design/navigation-frames.md)

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
  Prerequisites: MAP.68, NAV.15, NAV.29, NAV.33, MAP.79.
  Plan (2026-10-07): Moved into phase 0 with the engine (all bugs in
  phase 0).

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
  nearest-neighbour graph split into pieces (fixed by NAV.34, PR #427:
  `navGraph.join_islands` joins them, so a route always exists). Done:
  same-sector routes may leave the sector; the longest hop is shown; each
  hop is flagged as a jump through unknown space when its line crosses
  one or more unfilled (ungenerated) sectors (the default reading of
  "unknown space"; NAV.38's `galaxyGeometry.sectors_along_segment`,
  done in PR #357, finds the sectors), and `/api/nav` returns the flag
  per hop. Built with NAV.10, which already rebuilds the routing.
  Three strict-xfail tests in `src/tests/test_route_edge_cases.py` pin
  it (TEST.79, PR #427): `test_nav_between_reports_the_longest_hop`
  (`route["longest_hop_ly"]`), `test_nav_between_flags_a_hop_through_unfilled_sectors`
  (`route["hops"][i]["unknown_space"]`) and
  `test_nav_between_same_sector_route_uses_nearer_stars_next_door` (a
  same-sector route that leaves the sector); NAV.12 turns them green and
  removes the xfail marks (it may rename the keys).
  Phase 1 is its anchor: NAV.34 and TEST.79 are done (PR #427), UX.35 runs alongside it (phase 1), and NAV.36
  and NAV.39 need it first (phase 2).
  Design: [docs/design/course-routing.md](design/course-routing.md)

- [ ] **NAV.36 Unknown-space jumps drawn red and glowing**
  Boss (2026-10-02 01:53Z): "a jump through unknown space is marked in
  red and glows to draw attention to it." Done: a hop NAV.12 flags as
  an unknown-space jump is drawn red with a glow in UX.35's horizontal
  route strip, on the NAV page's own map (`lib/navmap.py`) and on the
  Galaxy Map course (with NAV.20), with a legend entry; with
  `prefers-reduced-motion` it stays red without the pulse; it reads in
  both themes, and the route list also labels it in text so it isn't
  shown by color alone. Prerequisites: NAV.12, UX.35, NAV.20.
  Design: [docs/design/course-routing.md](design/course-routing.md)

- [ ] **NAV.39 Saved courses remember their unknown-space jumps and check them again**
  Done: a saved course (NAV.17) keeps which hops were unknown-space
  jumps when it was saved, and opening it checks again, since sectors
  may have been filled in the meantime; a hop that is now known shows
  as ordinary, and the course says what changed. Prerequisites: NAV.12,
  NAV.17.
  Design: [docs/design/course-routing.md](design/course-routing.md)

- [ ] **NAV.43 Find everything within a distance of a place: the query and the API**
  Boss (2026-10-02 04:39Z): "Add TODO items for I want the ability to
  select a location and ask the sytem what stuff (of any kind) is within
  some distance in parsecs." Checked on main: the only proximity search
  is `queryDb.systems_within_radius` (the `queryDb.py` CLI and `GET
  /api/systems/<id>/near?radius=`), which takes a system, a radius in
  light-years, and looks only at systems in that system's own sector.
  Done: one query, `queryDb.objects_within` with a matching `GET
  /api/near` (and a `queryDb.py near` command), takes a place (an object
  reference from NAV.7, such as `system:<id>`, `planet:<id>` or a
  phenomenon, or a bare galaxy-frame point in parsecs) and a distance in
  parsecs, and returns everything of any kind within that distance:
  systems with their stars, planets, moons, belts and comets (placed at
  their system's position), facilities, phenomena of every type, quasars
  and rogue objects, each with its NAV.7 reference, kind, name, parent
  and distance, nearest first, paged at 50 rows, with an optional kind
  filter. Approach: `galaxyGeometry.enumerate_sectors_within_radius`
  (the cube geometry used by the neighborhood generator) names the
  sectors the sphere can reach, the query reads only those sectors' rows
  by their indexed `sector_id` (systems, rogue objects) or bounding box
  (placed phenomena, as `_placed_phenomenon_rows` does, which already
  accounts for an extended phenomenon's own radius), then measures exact
  distances in the galaxy frame; no new index is needed, and NAV.10's
  position indexes speed it up later if they land first. The old
  `systems/<id>/near` endpoint keeps working, answered by the new query.
  Tests: a known neighborhood returns exactly the objects inside the
  sphere, sorted, across a sector boundary. Open questions for Boss,
  with the defaults taken: the largest distance allowed (default 50 pc, about
  8,000 sectors at 4 pc a side); whether ungenerated sectors inside the
  sphere are reported (default: listed once as a count of sectors "not
  generated yet", never generated by the search); and whether only
  generated objects count (default yes). Prerequisite: NAV.7.

- [ ] **NAV.44 A "What's nearby" page: pick a place, enter a distance in parsecs, list what is there**
  Boss (2026-10-02 04:39Z, with NAV.43): "Add TODO items for I want the
  ability to select a location and ask the sytem what stuff (of any
  kind) is within some distance in parsecs." Done: a "What's nearby"
  page (linked from the menu and from every system, phenomenon and
  sector page) where the place is chosen with the same start and
  destination controls as the NAV page, or arrives from a link, and the
  distance is typed in parsecs; the result table lists every object
  NAV.43 returns with its kind, name (a link to its page), parent and
  distance on the site's distance ladder, filters by kind, pages at 50
  rows like every list, and says how many sectors in range are not
  generated yet. Prerequisite: NAV.43.

- [ ] **NAV.45 "What's within N pc" from the Galaxy Map and Sector Map**
  Boss (2026-10-02 04:39Z, with NAV.43): "Add TODO items for I want the
  ability to select a location and ask the sytem what stuff (of any
  kind) is within some distance in parsecs." Done: on the Galaxy Map and
  the Sector Map, a selected sector, system, phenomenon or picked point
  gets a "What's within N pc" action on the shared control panel that
  opens NAV.44 with that place filled in. Optional (default left out
  unless Boss asks): the map draws the search sphere and highlights the
  objects inside it. It uses MAP.65's shared control panel and NAV.15's picking of
  any object on the maps. Prerequisites: NAV.44, MAP.65, NAV.15.

- [ ] **NAV.46 The NAV picker can't click galaxy wedges to zoom in (bug)**
  Boss (2026-10-03 05:38Z): "When trying to select a destination from
  the picker in the galaxy nav screen, I cannot click on any galactic
  wedges to begin zooming in." Done: in course-pick mode the wedges,
  slabs and blocks pick and zoom as on the Galaxy Map.
  Prerequisite: NAV.15.

- [ ] **NAV.47 Unknown-space jumps stop at scattered stars, black holes, neutron stars and quasars**
  Boss (2026-10-03 05:38Z): "navigational jump path (system to system)
  should use scattered stars in unfilled sectors as potential stops
  along the path.  Same for black holes, neutron stars, quasars, but the
  by far preference is system-to-system, only if a jump crosses unknown
  space should it start looking for other things." Done as he says; such
  stops are marked as their kind in the route.
  Prerequisites: NAV.12, GEN.100.
  Design: [docs/design/course-routing.md](design/course-routing.md)

- [ ] **NAV.48 Offer to generate the uncharted sectors that block a course**
  Boss (2026-10-07 11:47Z): "From navigation menu if there are
  unexplored / uncharted / unfilled sectors in the way that cannot be
  bypassed, it should have the option to generate all sectors between
  the two points." Done: when a route has to cross uncharted sectors,
  the NAV page offers (to admins) a job that charts the sectors along
  the line, with the estimate first, and re-plots after.
  Prerequisite: NAV.12.
  Design: [docs/design/course-routing.md](design/course-routing.md)

- [ ] **NAV.49 Waypoints: pick objects in Star select mode and plot a course through them, kept on the map until cleared**
  Boss (2026-10-07 11:47Z): "Then you can go find another location
  (star, planet, whatever) and select it and you can add a 2nd waypoint,
  waypoints are highlighted at all map levels with an icon it will plot
  a course." And: "Plotted courses should appear on the galactic map and
  stay until cleared." Done: waypoints set from Star select mode are
  marked with an icon at every map level, two or more plot a course, and
  the course stays drawn on every map until cleared.
  Prerequisites: MAP.122, NAV.16, NAV.17.
  Design: [docs/design/course-routing.md](design/course-routing.md)

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
  Prerequisite: GEN.33.
  Plan (2026-10-07): Moved to phase 2 under GEN.90.

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
  Class S landed first, with GEN.38, as this item's first one-class PR
  (PR #415): barren rocky super-Earth, 1.2-1.8 Earth radii, 2-10 Earth
  masses, not allowed as a moon, weight 0.03. The other six classes
  (R, U, W, X, Y, Z) follow in phase 1.
  Prerequisite: GEN.33.
  Plan (2026-10-07): Moved to phase 2 under GEN.90. Class Z covers
  Boss's 2026-10-03 ask: "We need a planet that is like earth sized but
  never had any life so it has an atmosphere that was never matabolized
  by life."

  - [ ] **GEN.33 One class per PR, each with its tests**
    Suggested split: R and S (the commonest missing types) first; then
    U, W, X, Y and Z. Each adds the class to `PLANET_CLASSES`, the
    reference pages, the rogue flag, moon eligibility, and a
    distribution test over 1,000 generated systems.
    Prerequisite: GEN.89.
    Plan (2026-10-07): Moved to phase 2 under GEN.90 (the class
    refactor), after the habitability score GEN.89.

- [ ] **GEN.29 Sweep every planet class for sense once the new ones are in (bug)**
  Boss (2026-10-01 15:26Z): "do a full sweep of planet classes to make
  sure they all make sense logically once the new classes are in
  place." After GEN.28. Done: every class's description, composition,
  zones, sizes, weights, temperatures and atmosphere agree with each
  other. Known oddities to settle: D allowed in the hot zone (icy
  bodies at 265-490 K); C a catch-all for 63% of cold planets; Q never
  generated (weight 0.0001, and orbits are circular); V's composition
  "iron, iridium, tungsten"; L with vegetation at a median 0.02 bar; E
  at 376-414 K, above water's boiling point at 0.6 bar. Rocky rogues of
  10-16 Earth masses come out up to 17,600 km, past class S's
  11,468 km top, and get S as the nearest fit (PR #415); Boss to decide
  whether those sizes make sense or need their own range.
  Prerequisites: GEN.28, GEN.27, GEN.91, GEN.92.
  Plan (2026-10-07): Under GEN.90: the refactor is the sweep, so it
  stays with it in phase 2.

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

- [ ] **GEN.55 A version number and a seed reproduce the same galaxy (end goal)**
  Boss (2026-10-02 01:40Z): "Ok use a 128 bit value and store the seed in
  the database, and put it in the log at the top of any generation, also
  populate the TODO upward from here to eventually build a system that a
  version number and a seed value would reproduce the same galaxy by the
  end of the phases." This item is the
  chain that gets there; it is done when its last sub-item is. Order:
  - Phase 0: GEN.39 (per-unit seeds, after PERF.21) with DB.6 and
    OPS.10.
  - Phase 1: GEN.56 (every
    draw seeded), GEN.57 (a sector's contents depend only on the seed,
    the version and its address), DB.7 (the version kept with each
    sector), GEN.58 (a fingerprint), TEST.77 (the golden-galaxy
    test), ADM.18 (the creation settings as a JSON file), GEN.59 (admin
    changes as a net difference in that file), OPS.13 (each
    update records the key, keeping the last 10) and OPS.14 (a warning when the running key differs from the
    galaxy's). GEN.47's nebula field (done, PR #419) uses the derived seeds.
  - Phase 2: PERF.18 and GEN.42 give the one-process stars for one
    seed (already in their Done text); API.16 and ADM.17 show the seed
    and version; OPS.15 has each update say whether it changes
    generated output.
  - Phase 2: the daily maintenance run (OPS.16, with OPS.17's
    schedule) merges pending admin changes into a new JSON file
    (GEN.61) and keeps 18 backups (OPS.18, listed by ADM.19).
  - Phase 3: API.17 makes remote generation reproduce the server's;
    DB.10 repairs from the newest JSON plus pending deltas; ADM.20 adds
    a "merge now" button (low priority).
  - Phase 3+ (the end state): OPS.12, `generate.py reproduce`.
  Anything that draws new randomness later (GEN.42, PERF.18,
  API.12, API.13) uses the derived seeds and keeps TEST.77 green.
  Design: [docs/design/reproducible-galaxies.md](design/reproducible-galaxies.md)

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
    above also fails on direct calls to the others. Generators also must
    not depend on the iteration order of sets or dicts keyed by strings
    (which follows `PYTHONHASHSEED`) or on locale-dependent sorting; the
    test runs a small generation under two `PYTHONHASHSEED` values and
    two locales and compares. GEN.39 is done (PR #381).
    Design: [docs/design/reproducible-galaxies.md](design/reproducible-galaxies.md)

  - [ ] **GEN.57 A sector's contents depend only on the seed, the version and its address**
    Seeding every draw (GEN.56) isn't enough when the result depends on
    what ran first, at what worker count or how fast. Known order
    dependencies: name-collision decorations are resolved by how many
    names already exist (`resolve_greek_roman_collision(base_name,
    existing_count)` in `nameUniqueness.py` and the reservation in
    `_db.py`), so the sector that reserves first gets the plain name
    (GEN.46, done in PR #370); population seeds key on database ids
    (`random.Random(planet_id)` and `random.Random(species_id * 7919)` in
    `population.py`), and ids come from per-worker id blocks; the
    backfill depends on which sectors are already filled and where the
    run started (with GEN.44); and nearest-system links depend on which
    neighbours exist yet. Done: a sector's contents (as docs/design/reproducible-galaxies.md defines
    them; ids and timestamps excluded) come out the same whichever
    sectors were generated before it, at any worker count; names, seeds
    and skips key on addresses and the galaxy seed, not on ids or arrival
    order: in a name collision the sector with the lower address keeps
    the name and the other is renamed from its own seeded stream;
    nearest-system links are rebuilt from content, so they are left out
    of the comparison; and a test generates the same sectors at 1 and 4
    workers and in two orders and compares fingerprints (GEN.58).
    GEN.44's per-sector luminosity levels (Boss, 2026-10-02 03:25Z) let
    a sector be filled in steps (down to 1000 L_sun, later down to 500);
    the stars drawn between two floors must not depend on which steps
    ran, so draws key on the sector and the luminosity range, not on the
    order of runs.
    Since GEN.64 (done, PR #406), interstellar objects and bright-sweep
    systems are named by a position ID and drop out of the name-collision
    rule, which now covers named star systems. The ID's collision counter
    (`_db._claim_object_ids`) checks rows already stored, so when one ID
    cell holds objects from two sectors, the counter depends on which
    sector saved first; this item keys it on the lower address too.
    GEN.44 (done, PR #425) made per-sector backfill draws exact: a sector
    taken down to 500 L_sun in two steps gets the same stars as in one.
    The galaxy-wide layer scatter is not split by band, because that was
    5x slower (8.8 s against 1.7 s a layer), so a galaxy-wide
    `--bright-stars-down-to` band is statistically equal to one deep
    scatter but not star for star. Open question for Boss: is that good
    enough, or must the layer scatter be star-identical too (at that
    cost)? Default: statistically equal.
    Prerequisite: GEN.56.
    Design: [docs/design/reproducible-galaxies.md](design/reproducible-galaxies.md)
    Plan (2026-10-07): Name collisions disappear with GEN.67 (names come
    from IDs), so the name part of this item is dropped, and GEN.63 with
    it.

  - [ ] **DB.7 The version that generated each sector, and a warning for mixed-version galaxies**
    Done: each sector row records the PlanetGen release that generated
    it, as DB.6's packed hex version and the full string (DB.6, done in
    PR #387, keeps the galaxy's first version in `galaxy_shape.version_key`
    with `planetgen_version`, `python_version` and `platform`, galaxy
    schema v52; each sector gets the same four columns); extending a
    galaxy with
    a different release warns before it starts (CLI and Generate page),
    because a mixed-version galaxy reproduces only sector by sector, each
    on its own version. Prerequisite: none (DB.6 done, PR #387).
    Design: [docs/design/reproducible-galaxies.md](design/reproducible-galaxies.md)
    Prerequisite: DB.11.
    Plan (2026-10-07): Written as an Alembic migration (DB.11).

  - [ ] **OPS.13 Every update records the version key, keeping the last 10**
    Boss (2026-10-02 02:08Z): "I approve the plan, every time the script
    updates it will have to make a sweet to recalculate the seed value
    based on the current code version and it'll need to keep the
    seed/version combo for say the last 10 updates, that won't take up
    much space in the control database." Done: `update.sh` and
    `update.ps1`, after updating the code, compute the current key
    (DB.6's helper, `versionKey.version_key`, PR #387) and add one row
    per galaxy to a new control-database
    history table: the galaxy seed, the key, the SHA-256 of the nltk
    `words` corpus files, `offensive_words.txt` and any name lists, the
    SHA-256 of `requirements.lock`, and the date. Only the last 10 rows
    per galaxy are kept, and `generate.py` can list them. The galaxy
    seed never changes on update. A test runs the step 11 times and
    checks the oldest row goes. Lands after OPS.7 and OPS.8, which
    change the same two scripts. Open question: Boss wrote "recalculate
    the seed value". The default reads that as recording the current
    key next to the unchanged seed, because a changed seed makes a
    different galaxy. Prerequisite: OPS.8.
    Design: [docs/design/reproducible-galaxies.md](design/reproducible-galaxies.md)
    Plan (2026-10-07): The nltk corpus, offensive_words.txt and
    name-list hashes are dropped once GEN.71 lands; the lock-file hashes
    stay.

  - [ ] **OPS.14 A warning when the running version key differs from the galaxy's**
    Done: one check compares the running key (DB.6) and the corpus and
    lock hashes (OPS.13) with the ones stored for the galaxy, and for
    each sector (DB.7) when it has them, and names each field that
    differs (for example "Python 3.12.3 now, 3.11.9 when generated").
    `generate.py` and the Generate page warn before a run that extends
    the galaxy; GEN.58's fingerprint output and OPS.12's reproduce
    report print the same comparison, so a mismatched fingerprint says
    whether the platform or the corpus changed too. A test checks each
    field is named. Prerequisites: DB.7, OPS.13.
    Design: [docs/design/reproducible-galaxies.md](design/reproducible-galaxies.md)

  - [ ] **OPS.15 Each update says whether it changes generated output**
    The second half of Boss's update "sweep". Done: after OPS.13's row
    is written, the update script fingerprints a small fixed region
    (GEN.58) from the galaxy seed under the new code, compares it with
    the fingerprint stored in the previous history row, stores the new
    one, and reports "generated output unchanged" or names the sectors
    that differ. The live galaxy is not touched. A test runs it across a
    change that alters a sector. Prerequisites: OPS.13, GEN.58.
    Design: [docs/design/reproducible-galaxies.md](design/reproducible-galaxies.md)

  - [ ] **OPS.16 A daily maintenance script for Linux, macOS and Windows**
    Boss (2026-10-02 02:28Z): "a JSON file is only changed with the
    deltas at the end of the day. We're going to have to build a
    maintenance script for powershell and bash that will run the
    positional update script, then kick off this delta script that will
    update the JSON so it's only updated once every 24 hours, old JSON
    files are kept in the following order, 1 year ago, 6 months ago, 4
    weeks ago, 7 days ago. A total of 18 backup slots so that we have a
    good span of the different deltas."
    Done: `scripts/maintenance.sh` (Linux and macOS) and
    `scripts/maintenance.ps1` (Windows) run once a day: first the
    positional update (`updateOrbits.py`), then the delta merge
    (GEN.61), then the backup rotation (OPS.18). A lock keeps two runs
    from overlapping, every step logs to the normal log, and the script
    exits non-zero on any failure. Optionally it also runs OPS.15's
    fingerprint check, so changed output is noticed daily. A test runs
    it on a small galaxy and checks a second run started during the
    first exits at once. Prerequisites: GEN.61, OPS.18.
    Design: [docs/design/reproducible-galaxies.md](design/reproducible-galaxies.md)

  - [ ] **OPS.17 Install and update set up the daily maintenance schedule**
    Done: `install.sh`, `install.ps1`, `update.sh` and `update.ps1`
    (with `deploy-common.*`) set up OPS.16's daily run: a systemd timer
    or cron entry on Linux, launchd on macOS, Task Scheduler on
    Windows. Update leaves an existing schedule as it is and adds a
    missing one. The deployment docs say how to change the time or turn
    it off. Same scripts as OPS.7, OPS.8 and OPS.13, so it lands after
    them. Prerequisites: OPS.16, OPS.8, OPS.13.
    Design: [docs/design/reproducible-galaxies.md](design/reproducible-galaxies.md)

  - [ ] **OPS.18 Settings JSON backups kept in 18 slots: 7 daily, 4 weekly, 6 monthly, 1 yearly**
    Boss (2026-10-02 02:28Z): "old JSON files are kept in the following
    order, 1 year ago, 6 months ago, 4 weeks ago, 7 days ago. A total of
    18 backup slots". Done: after each merge, the JSON files are kept by
    grandfather-father-son rotation: the newest 7 daily files, then 4
    weekly, 6 monthly and 1 yearly, 18 in all, and older files are
    deleted. The current file is always kept. A unit test runs the
    rotation over a simulated year of dates and checks which files
    survive. Prerequisite: GEN.61.
    Design: [docs/design/reproducible-galaxies.md](design/reproducible-galaxies.md)

  - [ ] **GEN.58 A fingerprint of a galaxy's generated content**
    A way to tell whether two builds are the same in the sense
    `docs/design/reproducible-galaxies.md`
    defines. Done: `generate.py fingerprint` (for the galaxy or a region)
    prints a canonical SHA-256 digest per sector and one for the region over the compared content in
    a fixed order (address order, canonical number formatting), skipping
    ids and timestamps, either as first generated or with the JSON
    file's edits and regenerations applied (GEN.59); the same function backs GEN.57's test,
    TEST.77 and OPS.12.
    Design: [docs/design/reproducible-galaxies.md](design/reproducible-galaxies.md)

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
    Design: [docs/design/reproducible-galaxies.md](design/reproducible-galaxies.md)

  - [ ] **OPS.12 `generate.py reproduce`: a version and a seed rebuild a galaxy and check it**
    The end state Boss asked for (phase 3+). Done: `generate.py reproduce
    --seed X --version Y` rebuilds a galaxy, or a region, into a fresh
    database from the seed and the stored run history (DB.6), and checks
    it against the fingerprint (GEN.58) of the live galaxy or a given
    one, listing any sector that differs; it refuses, naming the version
    to check out, when the running release isn't Y. Without the edit and
    epoch layers (GEN.59) it rebuilds the galaxy as first generated; with
    the JSON file's net changes, regeneration seeds and epoch, as it is
    now. It reads the newest JSON file plus any pending deltas still in
    the control database (GEN.59), or rebuilds from any of the 18
    backups kept by OPS.18 (`--as-of DATE` picks one). A test runs it on a small galaxy, and on one with
    a deliberately changed sector. Simplest default: no automatic
    migration of old galaxies to a new release's output. It prints
    OPS.14's comparison of the stored and running key and hashes, and
    warns when they differ. It reads the galaxy's settings from ADM.18's
    JSON file. Prerequisites: DB.7, GEN.57, GEN.58, TEST.77,
    GEN.59, OPS.14, ADM.18, GEN.61, OPS.18.
    Design: [docs/design/reproducible-galaxies.md](design/reproducible-galaxies.md)

  - [ ] **ADM.17 The Generate page shows the galaxy's seed and version**
    Done: the Generate page shows the galaxy seed (32 hex digits, with a
    copy button), the version that made the galaxy and the run history
    (DB.6), and the new-galaxy form takes an optional seed (blank means a
    random one); the admin's System page shows the same for one system.
    Prerequisite: none (DB.6 done, PR #387).
    Design: [docs/design/reproducible-galaxies.md](design/reproducible-galaxies.md)

  - [ ] **ADM.18 The galaxy's creation settings saved as a JSON file, downloadable from the Admin dashboard**
    Boss (2026-10-02 02:13Z): "Inject into Phase 1 that the system
    should also take all the variables that were set at galaxy creation
    time and when the galaxy is created generate a backup file with all
    the settings in JSON format with the seed and version data. The file
    will be downloadable as a json file from the Admin dashboard at any
    time and a backup is made if the file ever has to be changed." Then
    (02:20Z): "The wordlist will also be in the JSON file", and (02:21Z):
    "JSON file name will be seed-version-date-time". Done:
    - When a galaxy is created, a JSON file is written holding only what
      reproduction needs: every creation setting (all `plan` and
      new-galaxy options, defaults included: disk scale length and
      height, bulge radius, max ring, bright-star floor, prevalence
      settings and the rest), the 128-bit seed, DB.6's 22-digit version
      key with its parts spelled out, the corpus and lock hashes
      (OPS.13), the word list, and GEN.59's regenerate seeds and net
      diff.
    - The word list is the filtered list the name generator actually
      uses (from nltk's `words` corpus, about 2.5 MB raw, after
      `offensive_words.txt` and any name lists are applied), stored
      gzip-compressed and base64-encoded with its SHA-256. A rebuild
      reads names from this list. This supersedes keeping only a hash of
      the corpus; the hashes stay for OPS.14's check.
    - The file name is `<32-hex seed>-<22-hex version key>-<YYYYMMDD>-<HHMMSS>Z.json`
      in UTC, with no colons so it is valid on Windows, for example
      `9F3A07C2E81B44D5A1C06E7B3D2F9081-0007007F000160030C0300-20261002-022133Z.json`.
      It lives in the site's data directory (path in `config.json`).
    - The file is written once when the galaxy is created. After that
      it changes only through the daily merge (GEN.61), which writes a
      new file under the new date-time and leaves the old one as a dated
      backup, so the newest file is the current one. Setting changes
      and admin changes wait in the pending-delta table until then.
      Backups are kept by OPS.18's 18-slot rotation.
    - The Admin dashboard offers the current file as a `.json` download
      at any time, with the control database's key history (the last 10
      updates, OPS.13) alongside it. ADM.19 later lists every backup.
    OPS.12 rebuilds from this file. A test creates a small galaxy,
    downloads the file and checks the name, settings, seed, key and word
    list.
    Prerequisites: DB.7, OPS.13, GEN.70.
    Design: [docs/design/reproducible-galaxies.md](design/reproducible-galaxies.md)
    Plan (2026-10-07): The JSON stores the galaxy's naming key (GEN.70)
    instead of the word list; the word list and corpus hashes go with
    GEN.71.

  - [ ] **ADM.19 The Admin dashboard lists the 18 settings backups for download**
    Done: the Admin dashboard lists every kept JSON file (OPS.18's 18
    slots) with its date, slot (daily, weekly, monthly, yearly) and
    version key, each downloadable as a `.json` file, next to the
    current one (ADM.18). Admin only. A test checks the list matches
    the files on disk and a download returns the file. Prerequisites:
    ADM.18, OPS.18.
    Design: [docs/design/reproducible-galaxies.md](design/reproducible-galaxies.md)

  - [ ] **ADM.20 A "merge now" button on the Admin dashboard (low priority)**
    Boss (2026-10-02 02:31Z), on an on-demand merge button: "Let's put
    that part of Phase 3, low priority." Done: an admin-only button on
    the Admin dashboard runs the delta merge (GEN.61) now, under the
    same lock and rules as the daily run (OPS.16), and writes a new
    seed-key-date-time JSON file. That file counts toward the day's
    daily slot in OPS.18's rotation. If the daily run holds the lock,
    the button says so and does nothing. A test presses it with pending
    deltas and checks the new file and the cleared rows. Prerequisites:
    OPS.16, GEN.61, OPS.18.
    Design: [docs/design/reproducible-galaxies.md](design/reproducible-galaxies.md)

  - [ ] **API.16 The API reports the galaxy's seed, version and run history**
    Done: an API route returns the galaxy seed, the version that made it
    and the run history (DB.6), documented with the API; API.12's
    download uses the same fields. Prerequisite: API.5.
    Design: [docs/design/reproducible-galaxies.md](design/reproducible-galaxies.md)

  - [ ] **API.17 Remote generation reproduces what the server would make**
    Done: a remote run (API.12's download, API.13's generation without a
    database) with the same seed and release produces exactly what the
    server would for those sectors, checked by fingerprint (GEN.58); and
    API.8 can verify an upload by re-running a sample of its sectors on
    the server and comparing. Prerequisites: API.12, API.13, GEN.57,
    GEN.58.
    Design: [docs/design/reproducible-galaxies.md](design/reproducible-galaxies.md)

  - [ ] **GEN.59 Admin changes stored as a net difference from the generated galaxy**
    Boss (2026-10-02 02:20Z): "Now track changes from original to new
    (skipping everything inbetween) made through the admin system,
    regenerate will generate a new seed for that specific whatever it is
    beingr regenerated and store that in the JSON storing only as mcuh
    as is required in the JSON to reproduce the identical data in the
    databse." Admin changes are stored in ADM.18's JSON file as a net
    difference from what seed + key would produce, not as a history.
    Done:
    - For each changed object the JSON keeps only its final state:
      changed fields with their current values, and deleted objects as
      tombstones.
    - A regenerate draws a fresh random 128-bit seed for that object,
      which replaces its derived seed, and stores it; edits made after
      the regenerate are recorded on top. A regenerate clears earlier
      entries for that object and its children. This replaces the
      SHA-256(sector seed || edit number) idea and today's
      `random.seed()` in `api/edits.py`.
    - An edit that puts a value back to the original drops its entry.
    - Objects are named by a stable address path (sector ring, layer
      and slot, then system, body and moon by generated index), never by
      database id, through a small stable-path helper (it may share
      NAV.7's reference work but must not use row ids).
    - Admin changes (edits, deletes and regenerate seeds, by stable
      path) go into a pending-delta table in the control database as
      they happen. The JSON file is not rewritten per change: Boss
      (2026-10-02 02:28Z): "a JSON file is only changed with the deltas
      at the end of the day". The daily merge (GEN.61) folds the
      pending deltas into a new file with the rules above.
    - The positional-update epoch (when `updateOrbits.py` last moved
      systems) is recorded with the deltas and in the JSON, and
      `updateOrbits.py` is checked to give the same positions for the
      same epoch.
    Rebuild = seed + key + JSON + epoch: "seed + key" gives the galaxy
    as first generated, and the JSON's diff and epoch give it as it is
    now. Changes made since the last daily merge live only in the
    database until the next one; DB.9's parity file protects them in
    between. A test
    edits, regenerates and deletes objects, merges the deltas, rebuilds
    a fresh database from seed + key + JSON + epoch, and gets the same
    content by fingerprint
    (GEN.58). Moved to phase 1, beside ADM.18. Uses the admin edit code
    (`adminEdits.py`, `editStore.py`). Prerequisites: GEN.56, GEN.58,
    ADM.18.
    Design: [docs/design/reproducible-galaxies.md](design/reproducible-galaxies.md)

  - [ ] **GEN.61 The daily merge folds pending admin changes into a new JSON file**
    Boss (2026-10-02 02:28Z): "a JSON file is only changed with the
    deltas at the end of the day. We're going to have to build a
    maintenance script for powershell and bash that will run the
    positional update script, then kick off this delta script that will
    update the JSON so it's only updated once every 24 hours, old JSON
    files are kept in the following order, 1 year ago, 6 months ago, 4
    weeks ago, 7 days ago. A total of 18 backup slots so that we have a
    good span of the different deltas."
    Done: a delta-merge step reads the pending-delta table (GEN.59) and
    the newest JSON file, applies the net-diff rules (latest value per
    field, tombstones, a regenerate clearing earlier entries for that
    object and its children, values back at the original dropped),
    records the positional-update epoch, and writes a new JSON file
    named by ADM.18's seed-key-date-time rule. Pending rows are cleared
    only after the new file is written and read back. With nothing
    pending and no epoch change it writes nothing. A test merges a set
    of deltas, rebuilds from the new file, and matches the live galaxy
    by fingerprint (GEN.58); another kills the merge before the check
    and finds the pending rows still there. Prerequisites: GEN.59,
    ADM.18.
    Design: [docs/design/reproducible-galaxies.md](design/reproducible-galaxies.md)

- [ ] **GEN.117 The Galaxy Map shows bright stars only in a thin band on the galactic plane (bug)**
  Boss (2026-10-07 18:43Z), after PR #461 (GEN.79) merged: "Here is the
  latest run of the bright star sweep, screenshots taken.  As you can
  see from these, when we look at a cross-section aligning with the
  galactic plane, all the bright stars are strictly kept in this narrow
  band." Screenshots: /mnt/project-files/bugs/gen-117/. Cause (Bugfixes
  lane 1, 2026-10-07 18:47Z, read from the code and scatter numbers;
  Boss's fresh run confirms): the generator does place bright stars
  above and below the plane, but the Galaxy Map draws only the 400 most
  luminous bright stars per tile (`GALAXY_TILE_MAX_BRIGHT_STARS`,
  `planetgen/db/query.py` around line 2424). Young blue stars on the
  plane outshine the old giants above it (about 1,000 to 2,500 Lsun),
  so the giants never make the cut. The fix belongs with MAP.116's
  per-zoom density and brightness rules: for example, reserve part of
  each tile's budget for stars spread across height, or pick by
  luminosity band instead of brightest first. Done: a tile's drawn
  bright stars span the layers its stars occupy, a test checks that
  spread on a tile with a crowded plane, and Boss's sweep shows stars
  above and below the plane. A generator test of the z spread lands
  separately with Bugfixes lane 1's flake PR.
  Prerequisite: MAP.116.

- [ ] **GEN.116 Error when generating a neighbourhood centred near the galaxy edge, awaiting Boss's error text (bug)**
  Back burner. Boss (2026-10-07 17:11Z), on GEN.65: "I do not have the
  error message, keep an eye out for it, but put it on the back burner
  for something to watch out for ... I think this error occurred when I
  was attempting to generate  a neighborhood when the center was close
  to the edge of the galaxy." GEN.65's tests (PR #476: 18 neighbourhood
  runs at the rim, just inside and outside it, the top and bottom
  layers, the core, mid-disk and around a filled sector, through both
  the Generate page job and the CLI) all pass, and nothing web-only was
  reproduced. Done: the error is reproduced from Boss's error text and
  fixed with a test, or Boss confirms it no longer happens.

- [ ] **GEN.66 Physics on scipy, and astropy constants and units**
  Today `keplerMotion.py` solves Kepler and Barker by hand and
  `physical_constants.py` holds 450 lines of constants. Done: Kepler
  solving uses `scipy.optimize` (vectorised where many bodies move at
  once), constants and unit conversions come from `astropy.constants`
  and `astropy.units`, results agree with today's within stated
  tolerances (tests), and the hand-rolled solvers are deleted.
  Design: [docs/design/library-migration.md](design/library-migration.md)

- [ ] **GEN.67 Names from IDs: replace word-salad name generation**
  Boss (2026-10-03 05:38Z): "Do away with name generation for any
  object, do away with our word-salad code entirely add a value in the
  control database that will be generated by random generation at the
  galaxy generation time (add the ability to change it via admin consol)
  so when the galaxy is generated we'll generate a random numbe." Boss's
  `gatedPhonemeCodec.py` (repo root) turns an ID and a key into
  pronounceable words and back. Done: every name is the codec's output
  for the object's ID under the galaxy's naming key, word-salad naming
  is gone, and the subitems are done.
  Prerequisites: GEN.68, GEN.69, GEN.70, GEN.71, GEN.72, GEN.73.
  Design: [docs/design/object-ids.md](design/object-ids.md)

  - [ ] **GEN.68 Research: the cheapest unique IDs for every object, unfilled sectors included**
    Boss (2026-10-07 11:47Z): "Research different ways to generate an
    object ID (star, sector, planet, etc) everything in the galaxy needs
    a unique ID, but I also need to create unique ID's for sectors that
    haven't been filled yet, etc.  Research the most computationally
    efficient way to do this." Done: a short comparison (GEN.64's packed
    position ID, address-derived IDs for sectors, hash-derived IDs,
    database sequences) for cost, collision risk, stability when content
    moves, and whether an unfilled sector can have its ID before any row
    exists; Boss picks one.

  - [ ] **GEN.69 A unique ID for every object, star systems and unfilled sectors included**
    Today only interstellar objects and bright-sweep systems carry
    GEN.64's 76-bit position ID; star systems keep generated names and
    sectors have only an address. Done: every sector (filled or not),
    system, star, planet, moon, belt, comet and phenomenon has an ID by
    GEN.68's method, stored and indexed, and stable across regeneration
    in place.
    Prerequisite: GEN.68.
    Design: [docs/design/object-ids.md](design/object-ids.md)

  - [ ] **GEN.70 A naming key in the control database, made at galaxy creation and changeable by admin**
    Done: a new galaxy draws a random naming key (from the galaxy seed's
    stream) into the control database; the admin console can change it,
    which renames everything at once with no rows rewritten (names are
    computed from ID and key); ADM.18's settings JSON stores it instead
    of the word list.
    Design: [docs/design/object-ids.md](design/object-ids.md)

  - [ ] **GEN.71 Remove the word-salad name code, nltk and the name registries**
    Done: names come from the codec everywhere (systems, stars, planets
    as "<system> I", moons, belts, sectors, phenomena), and `names.py`,
    `bodyNames.py`, `nameUniqueness.py`, `dedupeNames.py`, the nltk
    corpus and the name registry tables go, with an Alembic migration.
    Wide binaries: Boss (2026-10-07 17:11Z): "No, we should never have
    A I or such for planet names.  Adjust the algorithm to produce 2
    words from the name.  A says word 1 I, word 1 II, etc...  B planets
    say word 2 I, word 2 II, etc..." So a wide pair's codec name has two
    words; star A's planets are "<word 1> I", "<word 1> II" and so on,
    and star B's are "<word 2> I", "<word 2> II". No planet name carries
    "A" or "B", and a test checks it.
    Prerequisites: GEN.69, GEN.70.
    Design: [docs/design/object-ids.md](design/object-ids.md)

  - [ ] **GEN.72 A backfilled bright star should get a name only when its sector is generated (bug)**
    Boss (2026-10-03 05:38Z): "When the system generates a sector with
    an existing bright star that was created during a bright star back
    fill, the start should then, and only then, get a unique name rather
    than only it's ID, we'll still keep the ID of course but we need to
    change the star's name when we generate the star system's content."
    Done: a bright-sweep star shows its position ID until its sector is
    generated; at that point it gets its codec name and keeps the ID.
    Prerequisite: GEN.69.

  - [ ] **GEN.73 Nebulae don't get unique names (bug)**
    Boss (2026-10-03 05:38Z): "Nebulae should get unique names." Done:
    every nebula and remnant has a unique codec name from its ID, shown
    on the maps, lists and its page.
    Prerequisite: GEN.69.

- [ ] **GEN.74 One point-in-space object that keeps every coordinate system in step, used by every object**
  Boss (2026-10-07 11:47Z): "Introduce a centralized point-in-space
  object.  When any of it's coordinates are changed in any coordinate
  system (cartesian, spherical, cylindrical, etc...) it automatically
  updates itself to the new position, and integrate that into all of our
  objects, and we'll store with that the point-mass information with
  that moving forward." Boss's prototype is `spacial-position.py` (repo
  root, `SpatialPosition3D`: galactic, sector and system frames,
  Cartesian, cylindrical and spherical forms, velocity,
  observable-movement thresholds). Done: the class moves into the
  physics package with tests, every positioned object (sector, system,
  star, planet, moon, belt, comet, phenomenon) holds one, with its mass
  and `mu` beside it, and the stored columns map to it.
  Design: [docs/design/orbital-updates.md](design/orbital-updates.md)

- [ ] **GEN.75 A nebula shape from metaballs and warped noise, as a mesh**
  Boss (2026-10-03 05:38Z): "Generate a realistic shape mathematically
  for a nebula so on the map it will appear as a shape, not a sphere,
  not an ellipse, but a kind of bulbous region with an irregular shape."
  His recommended pipeline: 4 to 8 centres in an anisotropic ellipsoid
  as polynomial metaballs, coordinates warped by 3D domain-warped
  simplex noise, then marching cubes (or dual contouring) at isovalue T
  for a mesh, or raymarching with density max(0, F(x) - T) for shading.
  Done: each nebula has a seeded shape (centres, radii, warp settings)
  stored with it, a mesh built from it (marching cubes) with a low-poly
  level for the Galaxy Map, and a containment test that tells whether a
  point is inside; nebula sectors and `surrounding_cloud` use that test.
  Design: [docs/design/nebula-and-asteroid-field-classes.md](design/nebula-and-asteroid-field-classes.md)

- [ ] **GEN.83 A planetary habitability index (PHI)**
  Boss (2026-10-03 05:38Z): "Create a habitability index based on the
  pressure, temperature, composition, etc...  This will use several
  documents in the design docs folder: Atmospheric Toxicity.md, Chemical
  Habitability.md, Naturally Occurring Ionizing Radiation.md, Planetary
  Habitability and Speculative Xenobiology.md, Planetary Habitability
  Index.md, Speculative Xenobiology Extremes.md, Speculative Xenobiology
  Examples.md, Mathematical and Algorithmic Implementation of the
  Planetary Habitability Index.md" Done when the subitems are.
  Prerequisites: GEN.84, GEN.85, GEN.86, GEN.87, GEN.88, GEN.89.
  Design: [docs/design/habitability-index.md](design/habitability-index.md)

  - [ ] **GEN.84 Habitability design: one score structure and reconciled thresholds**
    The docs give two structures (PHI-4's Pressure, Temperature,
    Chemistry, Radiation with colour tiers, and the Xenobiology doc's
    PHI_bio, PHI_cpx and Phi_tech), conflicting CO, CO2 and
    mask-pressure thresholds, undefined constants, and a PHI_bio formula
    that doesn't reproduce its own table (habitability-index.md lists
    them). Done: one structure chosen by Boss, every constant and
    threshold fixed with its source, worked examples that reproduce,
    written into habitability-index.md.
    Design: [docs/design/habitability-index.md](design/habitability-index.md)

  - [ ] **GEN.85 Atmosphere species, partial pressures and mantle redox for every planet**
    Today `atmosphere` and `composition` are free text. Done: each
    planet and moon stores its mantle redox state and the partial
    pressures of O2, CO2, CO, N2, Ar, H2, H2O, CH4, H2S and SO2
    (columns, DB.13), drawn by class and redox, with the free text
    generated from them.
    Prerequisite: GEN.84.
    Design: [docs/design/habitability-index.md](design/habitability-index.md)

  - [ ] **GEN.86 Stellar activity (XUV, flares) and planetary magnetic fields**
    Done: stars store activity (saturation phase, L_XUV/L_bol, flare and
    particle-event rates by mass and age), and planets a magnetic field
    strength from mass, rotation and tidal locking (the spin from GEN
    item GEN.104 when it lands).
    Prerequisite: GEN.84.
    Design: [docs/design/habitability-index.md](design/habitability-index.md)

  - [ ] **GEN.87 Surface radiation dose**
    Done: each planet stores its surface dose from column mass (P0/g),
    magnetic field, cosmic rays (more inside a compressed heliosphere,
    using the nebula containment of GEN.75) and flares, plus an
    ozone-loss flag, following the radiation doc's tables.
    Prerequisites: GEN.85, GEN.86.
    Design: [docs/design/habitability-index.md](design/habitability-index.md)

  - [ ] **GEN.88 Hydrosphere and ocean chemistry**
    Done: each planet stores water fraction, ocean and land fraction,
    ocean depth or ice shell, and an ocean class (acid sulfate, neutral,
    soda, chloride brine, ice-sealed) with pH, water activity and
    phosphorus flags.
    Prerequisite: GEN.85.
    Design: [docs/design/habitability-index.md](design/habitability-index.md)

  - [ ] **GEN.89 The habitability score for every planet and moon**
    Done: every planet and moon gets the scores GEN.84 defines, with a
    colour tier and a human equipment profile (shirtsleeve, mask, mask
    and scrubber, pressure suit, full life support), shown on its page
    and searchable.
    Prerequisites: GEN.85, GEN.86, GEN.87, GEN.88.
    Design: [docs/design/habitability-index.md](design/habitability-index.md)

- [ ] **GEN.90 Refactor the planet classes around the habitability index**
  Boss (2026-10-03 05:38Z): "Refactor all planetary classes to better
  align with our habitability index system and the research for the
  different kinds of life, etc." Known conflicts with the research: N
  (Venus analog, about 737 K) carries life with no stage cap, E sits
  above the 122 C ceiling, K and L are below the Armstrong limit while L
  has land animals, P is uncapped, and no class is an Earth-size world
  with air that life never touched (GEN.28's Z). Done when its subitems
  are.
  Prerequisites: GEN.33, GEN.28, GEN.27, GEN.29, GEN.91, GEN.92.
  Design: [docs/design/habitability-index.md](design/habitability-index.md)

  - [ ] **GEN.91 Classes like S and V in the hot and cold zones**
    Boss (2026-10-03 05:38Z): "We need a planet class for like Class S
    and V but in the other zones, so that will require research on when
    an atmosphere is possible and when it isn't." Candidates from the
    research: a desiccated post-runaway super-Earth (hot), a
    hydrogen-rich Hycean super-Earth (cold or outer). Done: atmosphere
    retention rules (escape velocity, XUV, temperature) written down,
    then the classes added one per PR under GEN.33's rule.
    Prerequisites: GEN.33, GEN.85.
    Design: [docs/design/habitability-index.md](design/habitability-index.md)

  - [ ] **GEN.92 Life and its highest stage follow the habitability score**
    Done: whether a world has life, and how far it gets, comes from its
    scores instead of its class alone; life chemistry follows the
    atmosphere, radiation and solvent rather than the star's letter
    only; existing galaxies are unchanged until regenerated.
    Prerequisites: GEN.89, GEN.28.
    Design: [docs/design/habitability-index.md](design/habitability-index.md)

- [ ] **GEN.93 Nebula conditions in planet generation**
  Boss (2026-10-03 05:38Z): "Account for nebula temperature and
  conditions when generating planets and calculating surface conditions.
  You'll have to do a fesability study for our nebula classes in if a
  planet could even form there and what would happen if one did, how it
  would change the planets development and envirionment according to the
  science." Done when the subitems are.
  Prerequisites: GEN.94, GEN.95.
  Design: [docs/design/nebula-and-asteroid-field-classes.md](design/nebula-and-asteroid-field-classes.md)

  - [ ] **GEN.94 Feasibility study: can planets form in each nebula class, and what changes**
    Done: for each nebula and remnant class (A to W), whether disks
    survive (photoevaporation near O and B stars), what enrichment does
    (26Al and 60Fe heating), how extinction and cosmic rays change
    surface conditions, and a rule table for generation, written into
    nebula-and-asteroid-field-classes.md and reviewed by Boss.
    Design: [docs/design/nebula-and-asteroid-field-classes.md](design/nebula-and-asteroid-field-classes.md)

  - [ ] **GEN.95 Nebula conditions applied when planets and surfaces are generated**
    Done: generation knows a system's surrounding cloud while it runs
    (today `surrounding_cloud` is set only when a system is loaded) and
    applies GEN.94's rules to planet formation and surface conditions.
    Prerequisites: GEN.94, GEN.89, GEN.75.
    Design: [docs/design/nebula-and-asteroid-field-classes.md](design/nebula-and-asteroid-field-classes.md)

- [ ] **GEN.96 Generation directives for a sector (an override button)**
  Boss (2026-10-03 05:38Z): "When generating a sector should have the
  ability to click an override button and give directive (must have a
  certain density, must have at least x types of x stars, needs to have
  at least x habitable worlds within the sector, etc etc)." Done: an
  Override button on sector generation opens directives (density, at
  least N stars of a type, at least N habitable worlds and the like);
  the run draws until they hold or reports which it couldn't meet;
  `generate.py` takes the same as `--directive`.

- [ ] **GEN.97 Generate N random neighborhoods**
  Boss (2026-10-03 05:38Z): "Add the ability to tell the system to
  produce x number of random neighborhoods in the generation process."
  Done: the Generate page and `generate.py` take a count and a
  neighborhood radius and fill that many neighborhoods around random
  qualifying centres.

- [ ] **GEN.98 Bright-star backfill from the farthest generated boundary outward**
  Boss (2026-10-03 05:38Z): "Make star backfill after a sector or group
  of sectors is generated go from the edge of the farthest sector
  boundary in all directions by the specified distance in the
  algorithm." Done: after a run, the backfill region is the run's outer
  boundary pushed out by the backfill distance in every direction, not a
  radius around the starting sector.

- [ ] **GEN.99 Nebula volume backfill with the star types the nebula needs**
  Boss (2026-10-03 05:38Z): "Some Nebulae must have certain stars in
  them by their definition, so backfill nebula volume sectors to seed
  with stars of required types.  Any nebulae that requires backfill
  after it's created will trigger a back fill for the nebula volume down
  to 750 Solar Luminisities with the probability skewed to ensure
  dominance of the correct types of stars.  Existing backfilled stars
  should also be regenerated in-place (star location doesn't change)
  with the skewed probabilities favoring the correct class and similar
  luminocity as the original if possible (it won't always be compatible
  between what is required and what brightness to use, if it is NOT
  possible, then don't regenerate the star." Done as he says, using the
  nebula's shape for its volume.
  Prerequisites: GEN.75, GEN.100.
  Design: [docs/design/nebula-and-asteroid-field-classes.md](design/nebula-and-asteroid-field-classes.md)

- [ ] **GEN.100 Phenomena placed galaxy-wide first and kept when sectors fill**
  Boss (2026-10-03 05:38Z): "Move black hole (all types), neutron star,
  and quasar scattering to just before the bright star scatter during
  galaxy creation, use the same table as the star scatter for ease just
  add type data if needed." Boss (2026-10-07 11:47Z): "Have phenomena
  like neutron stars, nebula, quasars, and black holes generated through
  the entire galaxy first.  Thus do not regenerate those when a sector
  is filled in.  Account for all contents in a sector to stay where they
  have already been put if more content is being added." Done: black
  holes, neutron stars, quasars and nebulae are scattered galaxy-wide
  during creation, just before the bright-star scatter, in the
  bright-star table with a type column; a sector fill keeps everything
  already placed in it and only adds.

- [ ] **GEN.101 Fill order: nearest sectors first along a pruned Hilbert octree curve**
  Boss (2026-10-07 11:47Z): "Change fill algorithm to fill sectors
  nearest to the original sector first, circling outward using a
  geometric 3D space-filling curve via a pruned Hilbert octree over the
  $\mathbb{Z}^3$ integer lattice. Recursively traverse the octree from
  the origin sector outward to depth $D = \lceil\log_2(2R)\rceil$,
  restricting traversal to cells satisfying $\Vert{}v - v_0\Vert{}_2 \le
  R$ to achieve deterministic $\mathcal{O}(N)$ generation without the
  NP-complete overhead of general Hamiltonian path searches. The
  generated sequence must maintain face-adjacent unit steps
  ($\Vert{}v_{i+1} - v_i\Vert{}_1 = 1$) visiting every valid voxel
  centroid exactly once with zero backtracking, preserving radial
  locality as it scales to the boundary." Note: a pruned ball is not
  always walkable with unit steps and no backtracking; the item reports
  where that fails and what it does then. Done: neighborhood and radial
  runs fill in this order, deterministically, with a test of the
  unit-step and visit-once rules.

- [ ] **GEN.102 Investigate filling all near-zero-density void space at once**
  Boss (2026-10-07 11:47Z): "Investigate filling all void space (space
  that plots to have a 0 or nearly 0 star density) at once, since those
  sectors have the least in them, we can afford to fill them." Done: a
  count of such sectors, the time and storage a fill would take, and a
  recommendation for Boss.

- [ ] **GEN.103 Research where each star type and phenomenon belongs in the galaxy's structure**
  Boss (2026-10-07 11:47Z): "Research and add tiered structure
  information about where stars of different types would be in the
  galactic structure, what kinds of phenomena have similar restrictions,
  etc." Done: a design note (thin and thick disk, bulge, bar, halo,
  arms; ages and metallicity; where O and B stars, white dwarfs,
  remnants, clusters and nebulae sit) and the rules added to the density
  and type draws.

- [ ] **GEN.104 A spin vector and a realistic axial tilt for every rotating object**
  Boss (2026-10-07 11:47Z): "Add to the database every object gets a
  spin value as a vector containing direction of spin, velocity of spin,
  in 3 dimensional space.  Each rotating object will be generated a
  realistic axis tilt. To be clear when I say spin I am talking
  literally about the movement of a body around it's center of gravity
  (i.e. the rotation of the Earth)." Done: stars, planets, moons, small
  bodies and black holes store a spin axis and rate, drawn by the doc's
  cascade (tidal locking, gyrochronology for cool stars, log-normal
  speeds under the breakup limit for hot ones, the 2.2-hour spin barrier
  for small bodies, black hole spin distributions).
  Prerequisite: GEN.74.
  Design: [docs/design/orbital-updates.md](design/orbital-updates.md)

- [ ] **GEN.105 Orbital updates**
  Boss (2026-10-03 and 2026-10-07) asked for an orbital update that
  moves only what has visibly moved, counts what changed, and lets
  nearby masses bend paths, with "reasonable limitations for when the
  math breaks down at the edge cases". `updateOrbits.py` is today's
  positional update. Time step, Boss (2026-10-07 17:11Z): "No 1 year
  per orbital update or turn, rather, once set up and configured we
  follow orbital paths in real time.  The update script should have an
  option to update for more time in 1 go if specified.  Default is 1
  day = 1 day." So each run advances the galaxy by the real time since
  its last update (one day of real time is one day of orbit), and an
  option advances it by a stated extra span in one go. Done when its
  subitems are.
  Prerequisites: GEN.106, GEN.107, GEN.108, GEN.115, GEN.109, GEN.110.
  Design: [docs/design/orbital-updates.md](design/orbital-updates.md)

  - [ ] **GEN.106 Movement thresholds and a next-update-due column**
    Boss (2026-10-03 05:38Z): "orbital update script should only count
    items that it had to change the position for.  It should only update
    the position of items that changed location in a noticable way.  On
    a galactic scale that's going to be only if it has had sufficient
    time to move by at least 0.01 mpc (milliparsec), that includes
    stars, black holes, neutron stars, brown dwarves, basically anything
    star-like that isn't a planet.  For system scale (planets, binary
    stars, etc) it will need to have moved at least 0.01 AU to get an
    update.  For planetary scale (moons and other satellites) the
    position will only be updated if it has moved at least 100,000 km
    since the lsat positional update. Update the database schema so that
    we can as fast as possible check for items due for an update.  When
    we update the position on any item for any reason, we will update
    the table with when the next update should be due for that object."
    Done as he says: an indexed `next_update_due` per object, set from
    its speed and its threshold whenever it moves.
    Prerequisites: GEN.74, DB.11.
    Design: [docs/design/orbital-updates.md](design/orbital-updates.md)

  - [ ] **GEN.107 The update reports how many objects moved, changed sector, or entered or left a nebula**
    Boss (2026-10-03 05:38Z): "orbital update script should also say how
    many objects changed sector or entered / exited a nebulae." Done:
    the update's summary counts objects moved (only those past their
    threshold), sector changes and nebula entries and exits (by the
    nebula shape test).
    Prerequisites: GEN.106, GEN.75.
    Design: [docs/design/orbital-updates.md](design/orbital-updates.md)

  - [ ] **GEN.108 Orbital math limits: where each method breaks down and what happens there**
    Boss (2026-10-07 11:47Z): "Orbital path implementation, we need to
    ensure that reasonable limitations for when the math breaks down at
    the edge cases." Done: orbital-updates.md lists each edge case
    (near-parabolic orbits, e close to 1, time steps longer than an
    orbit, bodies inside a Roche limit or a Hill sphere, the galactic
    centre's singular potential, float precision at galactic distances)
    with the guard the code applies, and tests hit each guard. The
    guards follow "Computational Astrodynamics.md": the universal
    variable for near-parabolic orbits (with the corrected Kepler
    equation in orbital-updates.md 10.6), modified equinoctial elements
    for circular and equatorial orbits, phase wrapping for steps longer
    than an orbit, and Plummer softening at the centre.
    Prerequisite: GEN.106.
    Design: [docs/design/orbital-updates.md](design/orbital-updates.md)

  - [ ] **GEN.115 The galaxy's own gravity: a smooth disk, bulge and halo potential**
    Boss (2026-10-07 12:25Z): "We'll have to add a galactic gravitational
    gradient but we need to make sure that it's consistent with actual
    science." His research ("Computational Astrodynamics.md") gives the
    model: a Hernquist bulge, a Miyamoto-Nagai disk and an NFW halo (a
    flattened log halo as a setting), with the nearby point masses as
    perturbations on top and Plummer softening of 1 pc. Done: every
    galactic step adds the potential, its masses and scales are settings
    with the document's values at the default galaxy shape (lengths scale
    with the disk scale length otherwise), and tests check the rotation
    curve (229.3 km/s at 8.128 kpc, 213 to 231 km/s from 4 to 20 kpc),
    the vertical pull near the disk and a star staying at its radius over
    many runs.
    Prerequisites: GEN.106, GEN.108.
    Design: [docs/design/orbital-updates.md](design/orbital-updates.md)

  - [ ] **GEN.109 N-body influence from the nearest 10 bodies of equal or larger mass, with a Hill-radius warning**
    Boss (2026-10-03 05:38Z): "Add robust vector geometry updates for
    orbital path influenced by the nearest 10 gravitational bodies that
    are as large or larger than the object in question so that object
    effect each other's orbital path.  Again, we only change
    trajectories when the movement is updated.  Should fire off a
    warning if any object is foudn to be within the hill radius of
    another object." Following "Orbital Update Full Algorithm.md":
    Velocity Verlet over the top influencers from per-sector point-mass
    tables, swept-sphere collision checks, then Hill-sphere crossings
    handled per system. "Computational Astrodynamics.md" adds the
    routines: flybys deflect analytically, and a pair inside a mutual
    Hill sphere is sub-stepped by a 4th-order integrator with the Roche
    limit and contact checked each sub-step. Boss (2026-10-07 17:11Z)
    approved keeping the ring, layer and slot sectors with forces summed
    in galactic coordinates, "but verify the algorithm will work with
    our sector geometry." Rings hold different slot counts, so a
    sector's neighbours across a ring or layer don't line up one to one.
    Done: trajectories change only at updates; warnings go to the debug
    and activity logs; neighbours are found by distance from the
    sector's bounds, not by index; and a test confirms that the
    influencer search finds the same nearest bodies as a brute-force
    search for sectors at the core, at mid radius, at the rim, at the
    top and bottom layers and across the slot-0 seam.
    Prerequisites: GEN.106, GEN.108, GEN.115.
    Design: [docs/design/orbital-updates.md](design/orbital-updates.md)

  - [ ] **GEN.110 Rogue planet collisions: asteroid fields, merged giants and new stars**
    Boss (2026-10-03 05:38Z): "If rogue planets ever hit and at least
    one is terrestrial both are destroyed and they turn into an asteroid
    field.  If two rogue gas giant collide we assume they become a
    bigger gas giant losing only a small % of their mass, then we have
    to find out if they would have sufficient mass for nuclear fusion to
    occur and then add a star to the map if that's yes.  If a star does
    get added to the map that gets reported to admin as well and the
    star should get no planets and instead should be a planetary nebula.
    Anything like that" Done as he says, with momentum kept (the
    collision doc's inelastic merger) and every event in the admin
    report.
    Prerequisite: GEN.109.
    Design: [docs/design/orbital-updates.md](design/orbital-updates.md)

- [ ] **GEN.111 Email the admin when two objects are inside each other's Hill radius**
  Boss (2026-10-03 05:38Z): "If email is configured and that happens,
  call the API to send an email to admin the exact ID and location
  information for where the objects are located." Done: with SMTP
  configured, each Hill-radius warning from GEN.109 emails the admin the
  two IDs and their positions.
  Prerequisites: GEN.109, USR.3.
  Design: [docs/design/orbital-updates.md](design/orbital-updates.md)

- [ ] **GEN.112 Plan asteroid fields and belts as object systems for rendering**
  Boss (2026-10-03 05:38Z): "Begin planning for actually generating
  asteroid fields and asteroid belts as a system of objects but we're
  still going to treat them like a single object for positional updating
  and orbital speed, only for rendering will be generate a system of
  objects." Done: a plan (how many bodies, seeded from the field's ID so
  they are the same every time, size and spacing laws, level of detail)
  reviewed by Boss.
  Design: [docs/design/nebula-and-asteroid-field-classes.md](design/nebula-and-asteroid-field-classes.md)

- [ ] **GEN.113 Analyze the anomaly docs: which anomalies to add and how**
  Boss (2026-10-07 11:47Z): "Two additional files have been added, but
  these are part of adding more anomalies to the starmap so it'll be
  probably Phase 2 and 3 just to analyze those." The files are "Star
  Trek Anomalies and Science.md" and "Star Trek Anomaly
  Probabilities.md" (anomaly classes, per-sector rates with a lore
  multiplier, sensor-path encounters). Done: a list of anomalies worth
  adding, each with its real-science basis, rate, placement and how it
  is drawn, reviewed by Boss; note the docs' 20 ly sectors against our 4
  pc (about 13 ly) sectors.

- [ ] **GEN.114 Add the chosen anomalies to the starmap**
  Done: the anomalies GEN.113 chose are generated, stored, shown on the
  maps and searchable, one kind per PR.
  Prerequisite: GEN.113.

## PERF: Speed, caching, bulk generation and parallel work

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
  Prerequisite: PERF.24.
  Plan (2026-10-07): The backfill blocks become RQ jobs (PERF.24).

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
  Prerequisite: PERF.24.
  Plan (2026-10-07): Plans on RQ job results and the
  cachetools caches (PERF.24, PERF.25).

- [ ] **PERF.24 The work queue and web jobs on Redis with RQ**
  Today `stellarObjects/workQueue.py` (1,600 lines) runs its own process
  pool with database leases (`work_lease`), heartbeats and stale-run
  reclaim, and `jobRunner.py` runs the web jobs as process trees with
  state in the jobs folder. Boss chose Redis (2026-10-03), which
  overrides the job guide's broker-less `ProcessPoolExecutor`. Done:
  generation units (sectors, layer scatters, backfill blocks) and web
  jobs run as RQ jobs on Redis workers; the worker count, retries,
  timeouts and status come from RQ; a dead worker's job is retried, not
  lost (PERF.22's case); progress and log lines are published per job
  for ADM.22; results match the old queue for the same seed at any
  worker count; and `workQueue.py`, `jobRunner.py` and the `work_lease`
  table go. PERF.19's audit is its first step.
  Design: [docs/design/library-migration.md](design/library-migration.md)

- [ ] **PERF.29 Record which runs a partly filled sector still needs**
  Boss (2026-10-07 11:47Z): "Add a way to record and resume which
  sectors might be only partially filled by identifying what runs
  haven't been made yet." Done: each sector records which steps of its
  fill have finished (systems, phenomena, bright-star levels,
  population), so a check lists every partly filled sector and the steps
  it still needs; the level rule of GEN.44 (-1 untouched, a positive
  L_sun, 0 generated) stays the summary.

- [ ] **PERF.30 Finish an interrupted block or sector run on the next start**
  Boss (2026-10-07 11:47Z): "Add a mechanism to finish generating a
  block or sector on next start." Done: when the site or `generate.py`
  starts, unfinished runs found by PERF.29 are offered (web) or listed
  with a `--resume` flag (CLI) and finish from the step they reached,
  with the same result as an uninterrupted run.
  Prerequisite: PERF.29.

## DB: Database and schema

DB.1 shipped in 7.35.0 (PR #152). DB.2 to DB.5 done (PR #342, PR #347).

- [ ] **DB.8 Check a galaxy database and say whether it is damaged**
  Boss (2026-10-02 02:13Z): "add to TODO and inject into Phase 0 build DB
  consistency checking so that we can check a DB and tell if it's been
  damaged, then Phase 1 inject a DB repair using parity data stored in a
  file." Done: `generate.py check-db`, and a button on the Admin
  dashboard that runs it as a job, check the galaxy and control
  databases without changing anything. The checks:
  - Schema: tables, columns and indexes match the recorded migration
    level (sharing DB.4's shape detection).
  - Orphans: no planet without its system, moon without its planet,
    system without its sector, composition row without its field or
    comet, and so on.
  - Ids: every row id sits below `id_blocks`' next ids (the DB.3 case),
    and no blocks overlap.
  - Names: the name registries match the names in use.
  - Values: every sector address is inside the galaxy's bounds, no
    value is NaN or infinite, and every system passes
    `validation.check_star_system`.
  - Counts: the per-sector stats table agrees with the rows, once GEN.44
    and PERF.11 exist; the version keys are present and well formed,
    once DB.6 and DB.7 exist.
  The report lists each problem with the rows involved and ends with a
  pass or fail line per check; the exit code is non-zero when damage is
  found. `--sector` and `--region` limit it to part of the galaxy. A test
  damages a copy of a small galaxy in each of these ways and checks the
  matching problem is reported, and an undamaged one passes.
  Design: [docs/design/reproducible-galaxies.md](design/reproducible-galaxies.md)
  Prerequisite: DB.11.
  Plan (2026-10-07): The schema check reads Alembic's revision (DB.11);
  the name-registry check is dropped, since GEN.71 removes the
  registries.

- [ ] **DB.9 Repair a damaged galaxy database from a parity file**
  Boss (same message): "then Phase 1 inject a DB repair using parity
  data stored in a file." Done:
  - Every sector gets a content checksum, the canonical hash of GEN.58's
    fingerprint.
  - A parity file kept outside the database (path in `config.json`)
    holds Reed-Solomon parity over groups of sector exports, so any one
    damaged sector per group can be rebuilt; the group size sets the
    overhead. It is updated whenever sectors are saved or edited, and
    carries its own checksum so damage to it is caught too.
  - `generate.py repair-db` uses DB.8's check to find damaged sectors
    and rebuilds each from the parity file. Where parity can't, and the
    galaxy's seed and version key match the running code (OPS.14), it
    regenerates the sector from its seed and replays its edits from the
    edit log. It re-runs the check and lists anything it could not
    repair.
  - A test damages rows in a copy of a small galaxy, repairs them, and
    gets a passing check.
  Prerequisites: DB.8, GEN.57, GEN.58, OPS.14.
  Design: [docs/design/reproducible-galaxies.md](design/reproducible-galaxies.md)

- [ ] **DB.10 Repair reads the newest settings JSON and the pending deltas**
  From Boss's 02:28Z daily-merge rule (GEN.61). Done: DB.9's repair,
  where it regenerates a sector from its seed, applies the newest JSON
  file's diff and epoch plus any pending deltas still in the control
  database, so admin changes since the last daily run survive a repair.
  If the newest JSON file is damaged it falls back to the next backup
  (OPS.18) and replays the pending deltas on top, saying so. A test
  repairs a sector with both merged and pending changes. Prerequisites:
  DB.9, GEN.61, OPS.18.
  Design: [docs/design/reproducible-galaxies.md](design/reproducible-galaxies.md)

- [ ] **DB.11 The database layer and migrations on SQLAlchemy and Alembic**
  Today `_db.py` (9,800 lines) builds SQL strings over PyMySQL with its
  own pooling, and `migrateDb.py` replays `schema_vNN.sql.gz` fixtures.
  Done: tables are declared as SQLAlchemy models (MySQL 8.4 and MariaDB
  11.4 both supported), connections and transactions use SQLAlchemy's
  pool, Alembic carries migrations from the current schema (v53) forward
  with a baseline that recognises existing databases, and the id-block
  and batching behaviour is kept (TEST.81 and TEST.87 stay fixed).
  Design: [docs/design/library-migration.md](design/library-migration.md)

- [ ] **DB.13 Every stored value in its own column, not in JSON blocks, and indexed for search**
  Boss (2026-10-07 11:47Z): "Ensure to have all database information to
  be stored directly in the database and not in JSON blocks.  Optimize
  to be able to use any bit of information to search the database
  quickly with the right query.  This is to be considered a Phaser 0
  priority." Done: an audit lists every JSON or serialized column in the
  galaxy and control schemas and what reads it; each becomes real
  columns or child tables (as Alembic migrations), the searches in
  `queryDb.py` and the API use indexed columns, and a test checks no
  galaxy table keeps a JSON blob.
  Prerequisite: DB.11.

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
    Prerequisite: API.5.
    Plan (2026-10-07): Downloads the naming key (GEN.70) instead of name
    registries.

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

- [ ] **API.18 Generate by recipe: JSON for sectors, systems, planets, moons and phenomena**
  Boss (2026-10-07 11:47Z): "API should contain a 'generate recipe'
  function, the API would accept JSON for configurations.  Basically you
  need to make an adaptable JSON system for the API to accept to
  generate exactly what the user wants with all possible items
  configurable and varied degrees of randomness.  There should even be
  an API JSON input to create a galaxy, sectors, planets, moons.  More
  than that, phenomena of all types should be creatable by recipe.
  Validation errors should return error objects in JSON from the API but
  also have details, log output for the command, why it failed, HTTP
  code should be 400 or 422 depending.  A 400 error would be for a
  gibberish request, something that the system cannot handle because it
  makes no sense.  A 422 error should be returned for validation
  failures." Done: a documented recipe schema where any field can be
  fixed, ranged or left random; recipes for a sector, a system, a
  planet, a moon and each phenomenon run as jobs; errors come back as
  JSON with details and log output, 400 or 422 as he says.
  Prerequisites: ADM.21, GEN.96, API.9.

- [ ] **API.19 Galaxy-scale recipes: build a whole galaxy, piece by piece, from JSON**
  Boss (2026-10-07 11:47Z): "The idea is that one could go so far as to
  create an approximation of the known galaxy if they were willing to
  build the JSON file for it and the server had enough data to store it.
  In reality, I want the framework and theoretical capability, but I
  know the system cannot really handle that.  I DO however want the
  ability to perhaps (one bit at a time) build a simulated Milky Way but
  that's like so far down the road it doesn't even matter." Done: a
  galaxy recipe (shape, settings, seed) creates a galaxy, and recipes
  can be applied region by region to build one up over time.
  Prerequisite: API.18.

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
  Prerequisite: UX.40.
  Plan (2026-10-07): Folded into UX.40: the Generate page fields become
  Shoelace inputs.

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
  Prerequisite: PERF.24.
  Plan (2026-10-07): The worker count now means RQ worker processes
  (PERF.24).

- [ ] **ADM.21 Input validation on Pydantic models**
  Today `validation.py` checks systems and admin input by hand (ADM.5).
  Done: request bodies, admin forms and generation settings are Pydantic
  models with the same limits, validation errors list every field at
  once, and the old checks are deleted.
  Design: [docs/design/library-migration.md](design/library-migration.md)

- [ ] **ADM.22 Job logs streamed over SSE into Xterm.js, with native progress bars**
  From "Web UX and Job Management Guide.md": the Generate page and Queue
  page show a running job's log in an Xterm.js terminal fed by
  Server-Sent Events from the RQ job, and its progress in native
  `<progress>` bars (`generatejobs.js`, `generatefolds.js`, `jobs.py`,
  `progressRate.py`). SSE needs threaded or gevent workers under Apache
  (mod_wsgi settings documented). Done: logs and progress stream live,
  reconnect after a drop, and keep the full log for download.
  Prerequisite: PERF.24.
  Design: [docs/design/library-migration.md](design/library-migration.md)

  - [ ] **ADM.24 A failed action's log closes before it can be read (bug)**
    Boss (2026-10-03 05:38Z): "If action fails that has a log screen,
    pause on the log output waiting for user input to continue so they
    can read the log." Done: when a job with a log fails, its log stays
    open with the error in view until the user clicks Continue (web) or
    presses Enter (interactive console; non-interactive runs exit as
    now).
    Prerequisite: ADM.22.

  - [ ] **ADM.25 Error tracebacks don't reach the console and the web log window (bug)**
    Boss (2026-10-03 05:38Z): "Error traces should be output to the
    console and to the web-log-window when an action is running so they
    can be readily copy-pasted into documents for tracing." Done: an
    exception in a running action prints its full traceback to the
    console and to the job's web log, with a Copy button, as well as to
    the debug log.
    Prerequisite: ADM.22.

  - [ ] **ADM.26 The bright-star backfill shows no progress bar on the web (bug)**
    Boss (2026-10-07 11:47Z): "During star backfill progress bars do not
    appear on the web at all." Done: the backfill publishes progress
    like the sector run does, and the Generate page shows its bar and
    ETA.
    Prerequisite: ADM.22.

- [ ] **ADM.28 A simpler Generate page: layer specs, a Customize window and plain controls**
  Boss (2026-10-07 11:47Z): "Actual specs on layers on the generation
  screen.  Let's also have a customize button that brings up a special
  window with all the settings, the generate screen is getting a bit
  complex, we need to do a full rework to make it easier to use." Done:
  the page shows the galaxy's layer specs (count, height, extent, how
  many charted), the common actions stay on the page and every other
  setting moves into a Customize dialog, and its subitems are done.
  Prerequisites: UX.40, ADM.29, ADM.30, ADM.31, GEN.97.

  - [ ] **ADM.29 Fill a span of layers, rings or columns**
    Boss (2026-10-07 11:47Z): "Ability to specify a span of layers to
    fill in or a span of rings to fill in.  Or a span of columns to fill
    in." Done: the Generate page and `generate.py` take a layer range, a
    ring range or a column range (with the estimate shown first).

  - [ ] **ADM.30 Radial generation: a cylinder of N sectors around a point**
    Boss (2026-10-07 11:47Z): "From generate menu specify a radial
    generation from a point in a direct cylinder x sectors radius (1 =
    minimum for contiguous orthogonal connection between each sector and
    it's adjacent sectors." Done: a point (sector or coordinates) and a
    radius in sectors fill a cylinder around it; radius 1 fills the
    point and its face neighbours.

  - [ ] **ADM.31 Every generate action offers to show what it made on the Galaxy Map**
    Boss (2026-10-03 05:38Z): "Add a button from the generate screen to
    show generated sectors in galaxy view." Boss (2026-10-07 11:47Z):
    "Generating any space gives you the option to see that space in the
    galaxy viewer." Done: the Generate page has a "Show charted sectors"
    button, and every finished generate action (page, map menus, API
    job) offers "Show on Galaxy Map", which opens the map fitted to what
    it made with those sectors highlighted.

- [ ] **ADM.32 Add a star system to a sector: at the emptiest spot, at given coordinates, or at random outside every Hill sphere**
  Boss (2026-10-07 11:47Z): "Need a way to add a single star system to a
  sector, the computer will place in the area of lowest density if it
  can without crossing the hill sphere of another object, if the object
  cannot be placed then it will gracefully error out." And: "Add the
  ability to insert a sector or move a stand-alone sector into some
  location in a sector of admin's choosing.  They can specify
  coordinates (in 3 different ways) of the star system relative to the
  center of the sector or center of the galaxy (either needs to be an
  option).  The admin should also be able to select random placement
  within the sector so long as it is not in a hill sphere, the system
  will not move existing stars without admin approval. Admin can also
  choose to have the computer put the new star system in the least
  populated area, basically the spot that is farthest from everything
  that you can get within the sector." Done: from the sector page, an
  admin adds a new or existing stand-alone system by (1) the emptiest
  point (farthest from every object), (2) coordinates in Cartesian,
  cylindrical or spherical form relative to the sector centre or the
  galaxy centre (GEN.74), or (3) random; no placement may fall inside
  another object's Hill sphere, nothing existing moves without the admin
  confirming, and a placement that can't be made fails with a message
  and changes nothing.
  Prerequisite: GEN.74.

- [ ] **ADM.34 One admin menu per screen, holding only that screen's actions**
  Boss (2026-10-07 11:47Z): "All admin items are hidden under a simple
  menu and each menu list is customized to the actual screen we're in."
  Done: every page and map view has one Admin menu button (hidden for
  visitors) listing only the actions that apply there; inline admin
  panels are gone.
  Prerequisites: UX.40, UX.26, UX.31.

- [ ] **ADM.35 Full control from every screen: edit anything, regenerate with every input, backfill or erase what is in view**
  Boss (2026-10-07 11:47Z): "Everything can be edited, regeneration asks
  me for all possible inputs to regen any item in the system from the
  item detail screen, from the galactic screen (regen entire sectors),
  backfill bright stars again to any currently viewed wedge or slap,
  erase all contents if I ask, etc, full control from every screen."
  Done: every object's detail screen can edit each stored field and
  regenerate it with a form holding every generation input (defaults
  filled in); the Galaxy Map can regenerate or erase whole sectors and
  blocks and re-run the bright-star backfill for the wedge or slab in
  view; erasing always asks first.
  Prerequisites: ADM.34, MAP.120.

- [ ] **ADM.36 Change an object's trajectory vector**
  Boss (2026-10-07 11:47Z): "Give the user the ability to change an
  object's trajectory vector." Done: an admin sets a body's velocity
  vector (in any frame the position object supports) from its detail
  screen; the change is validated (escape, collisions), logged, and used
  by the next orbital update.
  Prerequisites: GEN.74, GEN.109.

## SEC: Security

The login protection of 2026-10-01 (SEC.1, SEC.20 to
SEC.28: the always-on activity log, per-address and per-username
lockouts, trusted devices, the password blocklist and hashing cost,
two-factor sign-in and the fail2ban example) shipped in PRs #217, #220
and #221.

Design: [docs/design/login-brute-force-protection.md](design/login-brute-force-protection.md)

- [ ] **SEC.29 Two-step sign-in on pyotp, QR codes on segno**
  Today `totp.py` (100 lines) and `qrcodegen.py` (900 lines, a copy of
  Nayuki's generator) do this by hand. Done: TOTP codes come from
  `pyotp` with the same secrets, step and window, so enrolled users keep
  working; QR codes come from `segno` as SVG; both old modules are
  deleted; and the two-step tests pass every run, including near a step
  boundary (TEST.72).
  Design: [docs/design/library-migration.md](design/library-migration.md)

- [ ] **SEC.30 Login and request rate limits on Flask-Limiter with Redis storage**
  Today `loginThrottle.py`, `api/limiter.py` and `api/loginguard.py`
  keep their own counters. Done: request limits use Flask-Limiter with
  Redis storage; the per-address and per-username lockouts keep their
  rules (docs/design/login-brute-force-protection.md) on the same
  storage; the activity log and fail2ban lines are unchanged; and the
  rate-limit tests pass alone and under `-n auto` (TEST.83).
  Design: [docs/design/library-migration.md](design/library-migration.md)

## TEST: The test suite

From the test suite plan of 2026-10-01 (Boss: "build me a list of TEST
items (TEST.x) to build our testing suite to cover all edge cases, and
prepare to do another solid full on bug hunt. Don't implement just plan
for it"). Each item ends with its areas in brackets. "Suspected bug"
means read from the code, not reproduced yet; the bug hunt confirms or
clears each one.

### Infrastructure and CI

- [ ] **TEST.72 Intermittent failure in the two-step (2FA) sign-in test (bug)**
  Found by the bug audit (2026-10-01, `bug-audit.md`): the two-step sign-in test failed once in a full run for PR #263
  and passed 3 of 3 times alone. Done: the failing case is found (loop
  it, including near a time-step boundary of the one-time code), the
  cause is fixed in the test or in the code it found, and the test
  passes on every run tried. [infra, SEC]
  Prerequisite: SEC.29.
  Plan (2026-10-07): Folded into SEC.29: the TOTP code moves to pyotp,
  and the test is fixed there.

- [ ] **TEST.83 Rate-limit tests fail under parallel load (bug)**
  `test_web_admin.py::test_real_login_keeps_rate_limit` passes alone
  but fails under `pytest -n auto` load, which points at timing (seen
  by the Galaxy map picker and arc thread, PR #369, 2026-10-02). The
  Routing groundwork thread (PR #427, 2026-10-02) saw 6 more rate-limit
  tests in `test_web_pages.py`, `test_web_request_limits.py` and
  `test_web_security_limits.py` fail under `-n auto` on MariaDB 10.11
  and all pass when rerun alone. Done: the failing cases are found
  (loop them under load), the tests or the rate limits stop depending
  on wall-clock speed, and they pass on every run tried, alone and
  under `-n auto`. [infra, SEC]
  Prerequisite: SEC.30.
  Plan (2026-10-07): Folded into SEC.30: the limits move to
  Flask-Limiter on Redis, and the tests are fixed there.
  Seen again (2026-10-07): Foundations' full local run after PRs #501
  and #504 (13,587 passed) had rate-limit tests among 10 load-only
  failures; the timing and dropped-connection part is TEST.93.

- [ ] **TEST.93 Timing tests fail and MariaDB drops connections under full-suite load (bug)**
  Foundations' last full local run (2026-10-07, after PRs #501 and
  #504; 13,587 passed) had 10 rate-limit and timing tests fail under
  the 13-minute load, each passing alone, and MariaDB dropped
  connections during the run; earlier runs saw the same kind of
  load-only failures. The rate-limit cases are TEST.83 (folded into
  SEC.30); this item covers the timing tests and the dropped
  connections. Done: the failing tests are named (loop the suite under
  load), the dropped connections are explained (pool size, timeouts or
  server limits on the test database) and fixed, timing tests stop
  depending on wall-clock speed, and the full suite passes on repeated
  runs. [infra]
  Cause found (Bugfixes lane 1, 2026-10-08): GET /api/databases opens a
  connection pool per listed schema, including every other worker's
  test database, and never closes them.

- [ ] **TEST.94 test_old_jobs_are_pruned raises JobBusy again: the job lock outlives the finished job (bug)**
  `src/tests/test_web_generate.py::test_old_jobs_are_pruned` failed once
  in Bugfixes lane 1's full run (2026-10-08) with JobBusy "job ... is
  still running" at `start_job`: `_wait_finished` saw the previous job
  finished before the runner released the job lock. It needs Redis
  (`PLANETGEN_TEST_REDIS_URL`) and passes 3 of 3 alone. TEST.90 fixed
  the same symptom before the move to RQ (done, PR #484); this is a new
  cause in the web jobs code (Foundations' area, PERF.24). Done: a job
  counts as finished only once its lock is released (or the test waits
  for the lock), and the test passes repeatedly under `pytest -n auto`.

- [ ] **TEST.92 The web job runner's first-failure test reports the job as interrupted under full-suite load (bug)**
  Seen once by Foundations (2026-10-07) in a full local run of PR #496
  (PERF.24 steps 1-2): `test_runner_stops_at_the_first_failure`
  (`src/tests/test_web_generate.py`) reported the job as interrupted
  instead of stopped at the first failure. It passes 3 of 3 runs alone.
  Since PR #496, generation with more than one worker needs Redis and
  falls back to one worker without it; tests that monkeypatch
  generation must pin `PLANETGEN_WORKERS=1` or use
  `tests/worker_patches`. Done: the cause is found and fixed, and the
  test passes repeatedly under `pytest -n auto`.

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

  - [ ] **USR.8 Every signed-in user can generate a one-off system**
    Boss (2026-10-02 05:16Z): "TODO Item, all users can generate a
    one-off system, but you have to be logged in." Checked on main: the
    one-off system page (`/admin/generate/system` and its `/download`,
    `html/web/system_page.py`) is admins only (`_admin_or_redirect` /
    `_admin_or_403` from `generate_page.py`), and the only accounts
    today are admin accounts, so anonymous visitors already can't use it
    and nothing needs changing until user accounts exist. Done: once
    USR.2 adds the user role, any signed-in account (user, admin or
    Owner) can open the one-off page from a link outside the admin menu
    (for example `/generate/system`, with the old admin URL
    redirecting), with the same options and output; it still never saves
    to the database; anonymous visitors are sent to sign in; the rest of
    `/admin/generate` stays admin-only; a test covers anonymous, user
    and admin. Open question for Boss, with the default taken: a
    per-user limit, since each system runs the generator in a separate
    process for about a second (default: 30 one-off systems per hour per
    user account, admins and the Owner unlimited, a clear message when
    the limit is reached, the count kept in the control database).
    Prerequisite: USR.2.

## OPS: Installers, hosting, CI, releases

OPS.1 shipped with the version scheme in `changes/README.md`.

- [ ] **OPS.8 Update reloads Apache itself when run as root**
  Boss (2026-10-01): "it should just automatically reload apache2 if
  it's running as root." Today `update.sh` ends by printing
  `sudo systemctl reload apache2` (or `restart` after enabling a module)
  for the user to run. Done: when the update runs as root on Linux and
  apache2 is running, it reloads Apache itself (restarts it when a module
  was just enabled) and says so; when not root, or Apache isn't running,
  it prints the command as today. macOS (gunicorn) and Windows stay as
  they are unless Boss asks.

- [ ] **OPS.19 The Generate jobs folder is /var/lib/planetgen while the checkout is /var/lib/planetGen (bug)**
  Boss (2026-10-02 03:40Z): "the jobs folder goes under
  /var/lib/planetgen while the main folder is /var/lib/planetGen so it
  ends up in 2 folders." Linux file systems are case-sensitive, so a
  default install has two folders whose names differ only by case.
  The lowercase jobs default is set in two places: `DEFAULT_JOBS_DIR`
  in `src/html/web/jobs.py` (line 59) and in
  `examples/apache/deploy-paths.py` (line 25, which `install.sh` and
  `update.sh` use through `scripts/deploy-common.sh` to create the
  folder). It is repeated in `examples/apache/create-cache-dir.sh`
  (comment, line 79), `src/tests/test_deploy_scripts.py` (lines 80, 94
  and 123), `INSTALL.md` (line 147), `docs/deployment/README.md`
  (line 56), `docs/deployment/macos.md` (lines 18 and 110),
  `docs/deployment/windows.md` (line 86), `docs/server-checklist.md`
  (line 86) and the `jobs.dir` row of `docs/config.md`. Everything else
  (the checkout in the Apache, nginx, Caddy, systemd and macOS
  examples, `ci.yml`, `deploy-common.sh`'s path rewrite) already says
  `/var/lib/planetGen`. Not part of this: the lowercase
  `/var/cache/planetgen/tiles`, `/usr/local/planetgen/venv` and
  `/var/log/planetgen-*` are separate folders, not the checkout's
  parent. Done: the jobs default is `/var/lib/planetGen/jobs` in code,
  deploy scripts, tests and docs (the folder sits inside the checkout,
  so `jobs/` is added to `.gitignore`); `update.sh` moves an existing
  `/var/lib/planetgen/jobs` into the new folder (keeping job history,
  skipping a running job) and removes the empty old folder, saying
  what it did; a `jobs.dir` set in `config.json` is left alone; and a
  test checks the default and the move.
  Prerequisite: PERF.24.
  Plan (2026-10-07): Folded into PERF.24: the web jobs move to RQ, and
  the jobs folder default is fixed in the same PR.

- [ ] **OPS.20 Move the code base from zero dependencies to third-party open-source libraries**
  Boss (2026-10-03 05:38Z): "Transition the code to open source 3rd
  party libraries and simplify the code base, also in this step we'll
  use a redis server for use in managing and executing the work queue.
  Read document "Library Migration Workflow.md" for python libraries
  that we'll be switching to and we'll also begin using libraries for
  our web interface as well from files "Web UX Development Notes.md" and
  "Web UX and Job Management Guide.md".  To be clear, we have an
  explicit directive to move the code base from a 0-dependency model
  into using 3rd party open source libraries to simplify our own code
  deployment." This ends the zero-dependency rule (design decision 11).
  Each library swap is its own subitem and PR, after the package layout
  (OPS.23, done in PR #435, then OPS.24) so code moves once. Done when every subitem is done
  and the hand-rolled modules they replace are deleted.
  Prerequisites: SEC.29, UX.39, PERF.24, SEC.30, DB.11,
  ADM.21, GEN.66, UX.40, UX.41, ADM.22, MAP.102.
  Design: [docs/design/library-migration.md](design/library-migration.md)

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
    `planetgen/names/wordlists.py` gathered from constellation and star-group
    names across the world's languages and sky cultures (not just the 88
    IAU ones), and a constellation name generator that slices and
    recombines them into new names the same way stars, planets and
    sectors are named (`split_into_syllables` in `utils.py`, the
    prefix/suffix lists, the `offensive_words.txt` filter). Used by
    VIEW.3. Open questions: what counts as a source list (licensing of
    sky culture data such as Stellarium's), transliteration of non-Latin
    scripts, and whether names are unique per planet or galaxy-wide.
    Plan (2026-10-07): Constellation names come from the codec under the
    naming key (GEN.67), not from sliced word lists.

- [ ] **VIEW.5 Light-travel positions: where an object appears to a distant observer**
  Boss (2026-10-03 05:38Z): "Begin laying groundwork to produce a
  snapshot of what any point in the galaxy would see.  The groundwork
  that needs to be laid is calculations for calculating the previous
  position of an object say 1 ly away from where it would be to the
  viewer from 1 year ago (i.e. we never see where a distant star really
  is at this moment in time, only where it was)." Done: a function
  giving an object's apparent position from an observer point (position
  at now minus distance over c, solved iteratively for moving bodies),
  tested against a hand case.
  Prerequisites: GEN.74, MAP.70.

## POP: Population and politics

Population, species and polities (POP.1 to POP.6) shipped in PRs #169 to
#184.

- [ ] **POP.7 Tech levels for technological species**
  Boss (2026-10-03 05:38Z): "Integrate a technolog leveling system for
  planets that have a technological species present using the document:
  "Technological Assessment.md"  Use the data to develop a method of
  generating a tech level for every technological species." Done when
  the subitems are.
  Prerequisites: POP.8, POP.9.
  Design: [docs/design/population-and-politics.md](design/population-and-politics.md)

  - [ ] **POP.8 A tech-level design from the six domain indices**
    The doc scores six domains 0 to 7 (Energy, Materials, Information,
    Medical, Propulsion, Defense) and weights them: TL = 0.25 EI + 0.20
    MI + 0.20 II + 0.15 MeI + 0.10 PI + 0.10 DI. Done: how each index is
    drawn from a species' age, era and world (population.py's
    civilization age and era), written into population-and-politics.md.
    Design: [docs/design/population-and-politics.md](design/population-and-politics.md)

  - [ ] **POP.9 A tech level generated for every technological species**
    Done: every technological species stores its six indices and its
    tech level (columns), shown on the species page and searchable;
    existing species get theirs on the next population pass.
    Prerequisite: POP.8.
    Design: [docs/design/population-and-politics.md](design/population-and-politics.md)

- [ ] **POP.10 Facility types: programmable, picked from a dropdown, with affiliation and Green/Yellow/Red ratings**
  Boss (2026-10-03 05:38Z): "Add facility types that are programable and
  can be selected from a dropdown when adding facilities.  Facility data
  should include affiliation and a place to store data like a
  Green/Yellow/Red system for saying the facilities crime, housing,
  resources, maintenance, and health." Today facilities
  (`facilities.py`, ADM.9) have no type list. Done: admins define
  facility types (name, icon, fields) and pick one from a dropdown when
  placing a facility; each facility stores an affiliation and
  Green/Yellow/Red ratings for crime, housing, resources, maintenance
  and health.
  Design: [docs/design/population-and-politics.md](design/population-and-politics.md)
