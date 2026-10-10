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
rebuilt from the dependency report, and on 2026-10-07 rebuilt again
from Boss's lists of 2026-10-03 and 2026-10-07. Phase 0 (every bug and
the groundwork) was finished on 2026-10-09 and its plan file deleted;
the open bugs now sit in phase 1. Every open item below is in exactly one phase, and every
item's prerequisites are in its own phase or an earlier one. Each phase
file gives the goal, the build threads with their order and
prerequisites, and the open questions;
[plan/notes.md](plan/notes.md) keeps the research notes, the files
several items share and the judgment calls. This file keeps each
item's full text. A new item goes into a phase's table in the same PR
that files it.

| Phase | Plan | Goal | Items |
|---|---|---|---|
| 1 | [phase-1-built-on-roots.md](plan/phase-1-built-on-roots.md) | Object references and the picker finished, the database check and repair, resumable runs, the reworked Generate page (spans, radial fills, random neighborhoods, directives, show on the map), galaxy generation changes (phenomena placed galaxy-wide first, fill order, backfill from the run's edge, nebula volume backfill, star-type tiers), the habitability index, tech levels and facility types, spin and orbital-update thresholds with their limits, routing with unknown-space stops and charting along a course, the nearby search, Select mode and view filters on the maps, bookmarks, admin control from every screen, in-universe wording, the API log and scopes, and reproducible galaxies up to the golden-seed test. | NAV.8, NAV.9, DB.9, PERF.29, PERF.30, ADM.28, GEN.24, GEN.99, GEN.101, GEN.102, GEN.103, GEN.87, GEN.89, GEN.83, POP.8, POP.9, POP.7, POP.10, NAV.11, NAV.42, NAV.47, NAV.48, MAP.119, MAP.122, UX.46, UX.47, ADM.32, MAP.120, ADM.35, UX.23, UX.22, UX.3, UX.42, ADM.15, API.15, API.4, API.7, API.9, GEN.57, TEST.77, UX.49, PERF.31, GEN.128 |
| 2 | [phase-2-maps-picker-backfill.md](plan/phase-2-maps-picker-backfill.md) | The planet class refactor around the habitability index (with GEN.33, GEN.28, GEN.27 and GEN.29), nebula conditions on planets, n-body orbital updates and rogue collisions, editable trajectories, courses and waypoints, generate-by-recipe in the API, the backfill density pass, daily maintenance, the pilot-style visual design, light-travel positions, the asteroid-field and anomaly plans, and the API pieces remote generation needs first. | GEN.33, GEN.28, GEN.27, GEN.91, GEN.92, GEN.29, GEN.90, MAP.75, MAP.59, NAV.21, NAV.17, NAV.18, NAV.4, NAV.36, NAV.39, NAV.49, UX.32, UX.30, UX.43, GEN.42, PERF.18, GEN.40, PERF.20, API.5, API.10, API.11, API.12, OPS.15, OPS.16, OPS.17, NAV.45, UX.48, UX.45, GEN.115, GEN.109, ADM.36, GEN.110, GEN.105, GEN.95, GEN.93, GEN.112, GEN.113, MAP.121, MAP.132, API.18, PERF.33, MAP.139, MAP.140, MAP.141, MAP.142, MAP.143, GEN.129, GEN.130, ADM.43, ADM.44 |
| 3 | [phase-3-engine-3d-remote.md](plan/phase-3-engine-3d-remote.md) | The 3D system view and infinite zoom from galaxy to moon on one interface, orbital trajectories in their own frame, courses that bend around gravity wells, remote generation through the API (with galaxy-scale recipes) reproducing what the server would make, and the anomalies chosen in phase 2. | NAV.22, NAV.23, NAV.5, NAV.25, NAV.26, NAV.27, NAV.28, NAV.6, API.13, API.14, API.8, ADM.13, API.3, API.17, GEN.114, API.19 |
| 3+ | [phase-3plus-accounts-sky-galaxies.md](plan/phase-3plus-accounts-sky-galaxies.md) | The open-ended tail: user accounts (with API.6 keys, saved courses and Hill-radius emails), the view of the sky from a planet, the plan for more galaxies. | USR.2, API.6, USR.3, USR.4, USR.5, USR.6, USR.7, USR.8, USR.1, NAV.19, GEN.111, VIEW.1, VIEW.4, VIEW.2, VIEW.3, GEN.9, GEN.55 |

Phases overlap: a phase's later threads can start while the next
phase's first ones run, as long as the order inside each phase holds.

Version: the next revision is 8.1 (Boss, 2026-10-09 18:12Z: no 8.1 until
Phase 1 is complete). Boss (23:07Z) put the fly-through Galaxy Map
(MAP.146 and MAP.147 to MAP.155) into phase 1 as the main feature that
earns the 8.1 bump. Phase 1 is complete only when those are done too;
until then every release stays on 8.0 (`REVISION_HOLD` in
`scripts/bump_version.py`), and the TODO thread flips the hold when Boss
declares phase 1 complete.

Boss's list of 2026-10-01 23:53Z (`new todos.txt`, with research notes;
the files are in the project's shared files under `todo-tasks/research/`)
became these items:

| # | Boss's item | ID |
|---|---|---|
| 1 | Store every sector's backfill level (-1, lowest L_sun, 0) | GEN.44 (done, PR #425) |
| 2 | Rogue planets after systems and phenomena; expanded rows span the table | UX.24 (done, PR #487) |
| 3 | Rogue planet gas giant vs terrestrial probability | GEN.45 |
| 4 | Rogue planet octant and a map symbol link | UX.25 (done, PR #490) |
| 5 | Edit and admin buttons as a button menu | UX.26 (done, PR #558) |
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
  Plan (2026-10-07): The banner reads the RQ job's published progress
  (ADM.22), not progress.json.
  Research (2026-10-09, performance-eta-queue-and-caching.md): answer
  its open questions: the banner shows on every page, polled every 30 s
  from a public endpoint that reads a Redis progress record (5 s
  per-process cache); command-line runs publish the same record and
  refresh `planetgen:active:<db>` (60 s expiry) so they appear (open
  question for Boss, default yes); "no ETA yet" shows the PERF.3
  estimate labelled "from earlier runs" or "estimating". The hour rule:
  show ceil(1.15 x eta), fall freely, rise only when 0.88 x eta exceeds
  it (upward jumps of 9 to 18 per 10 h run fall to about 0; the
  under-promise share of 15 to 25% falls to under 1%, except 9% on
  two-cost-class runs).

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
  Research (2026-10-09, units-and-number-formatting.md): answer the open
  questions with the defaults in section 1 of the doc (ladder and switch
  points per quantity; atm, psi and bar for pressure; g and m/s² for
  gravity, ft/s² under Customary; °F for temperature; lb and mi only
  under the Customary preset). Secondary units show in detail panels,
  property grids and generator text; table rows and chips show the
  primary unit with the secondary in a tooltip. The rogue-planet page
  (`web/system_pages.py` `_temperature_text`, `_pressure_text`) still
  bypasses the shipped K/°C/°F and pressure formatters. Open questions
  for Boss (defaults taken): lunar-mass rung wanted (0.01 lunar mass to
  0.1 Earth mass; Jupiter masses from 0.1 M_J; solar masses from 0.075
  M☉); keep "AU" and accept "au" on input; no fourth in-universe unit
  set.

  - [ ] **UX.23 A shared unit-ladder module**
    Done: one Python ladder module (in `stellarObjects/utils.py` or a new
    `units.py`) and its JavaScript twin, generalising the existing
    distance, speed and duration ladders, with a test that the Python and
    JavaScript twins agree for a table of values. Then one sub-item per
    quantity family: mass; temperature; pressure and gravity; density,
    luminosity and power.
    Plan (2026-10-07): The ladder is built on astropy.units (GEN.66).
    Research (2026-10-09, units-and-number-formatting.md): re-scope to a
    data-driven ladder table in Python, one generated JS data file and
    one algorithm per language, built on astropy constants. Step 0: the
    module, a generated `unitladders.js`, wrappers keeping the old
    names; behaviour-preserving, proven by a 600,000-value equivalence
    run against the current formatters, plus the rounding fix as the one
    deliberate output change. Replace the sub-item list with the
    migration order in section 5: distance, time,
    temperature/pressure/gravity, mass, luminosity/power/flux,
    density/number density/magnetic field/angles, JS stragglers,
    preference menu. Add the seven Hypothesis/parity tests of section 4
    to "Done". The text names `stellarObjects/utils.py` and
    `html/lib/fmt.py`; the files are now `src/planetgen/util/format.py`
    and `src/planetgen/web/lib/fmt.py`.

- [ ] **UX.30 Planet information without the Markdown render**
  Boss (2026-10-01 23:53Z): "Rework planet information displays to pull
  away from the markdown render and instead fits with the modern web UI
  as we've seen it so far." Today each body's expanded row on the
  system page is the generator's Markdown (`systemRender.render_system_sections`)
  turned into HTML by `lib/mdrender.py` (`systempage._row_html`), and
  so is the system overview; only the stat chips are built as HTML.
  Done: the system page builds each body's details from its stored
  values as structured HTML (property grids, atmosphere composition
  bars, orbit figures, moons, life), styled like the rest of the site in
  both themes and every size class, using the unit ladders (UX.22);
  `mdconvert` is no longer used for them. The wikitext and Markdown
  views (the Wikitext/Markdown toggle and wiki upload) keep using the
  text render. Open question: does the overview text stay as prose?
  Research (2026-10-09, units-and-number-formatting.md): prerequisite
  UX.23 step 0, so the new property grids emit `<span class="qty" data-q
  data-si>` from the start.
  Prerequisite: UX.23.

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

- [ ] **UX.49 Form fields as Shoelace components (sl-input, sl-select, sl-checkbox) across the site**
  UX.40 (done, PRs #528, #533, #539, #544) moved the buttons, menus and
  dialogs onto Shoelace and lined up the Generate page's text boxes with
  plain CSS; the form fields themselves are still native `<input>`,
  `<select>` and checkbox elements. Done: every text box, number box,
  drop-down and checkbox in the templates (Generate and one-off system
  pages, the edit dialogs, the search and filter fields, sign-in and
  account forms, admin pages) is an `sl-input`, `sl-select` or
  `sl-checkbox` from the vendored set, submitting and validating
  exactly as before (labels, required fields, keyboard use, values
  posted by name), both themes and phone widths look right, and a
  browser test fills and submits one form of each kind. Not a bug.
  Research (2026-10-09, map-ui-and-frontend-libraries.md): re-scope to
  "form fields use one console-token look; sl-input / sl-select /
  sl-checkbox only where a component feature helps; hidden, submit,
  password and one-time-code fields stay native". Done criteria:
  `--sl-input-border-color` set from the new `--line` token (3:1),
  server validation errors shown on the right field (ADM.21 errors), a
  no-JS check that each form still posts. Open question for Boss
  (default): Shoelace for select, checkbox, switch and the Generate and
  admin forms; native for sign-in, account, 2FA code, hidden and submit
  fields.
  Prerequisites: none.
  Design: [docs/design/library-migration.md](design/library-migration.md)

- [ ] **UX.42 In-universe wording across the interface**
  Boss (2026-10-03 05:38Z): "Begin replacing wording to make the
  interface in-universe (i.e. the interface and wording should not be
  designed so they sound like they know it's for a game)." Done: a word
  list (for example "Generate" becomes "Chart", "Generated only" became
  "Charted only" in MAP item MAP.111) agreed with Boss, then pages,
  buttons and messages changed page by page; admin-only tools may keep
  plain wording.
  Research (2026-10-09, map-ui-and-frontend-libraries.md): schedule with
  UX.43 so each page is restyled once. Shoelace 2.20.1 is its final
  release ("sunset"); open question for Boss (default: stay, and do the
  `sl-` to `wa-` rename inside UX.43's restyle if it rewrites
  `shoelace-theme.css` anyway).

- [ ] **UX.43 A visual design built like a pilot's starmap and navigation console**
  Boss (2026-10-07 11:47Z): "All design elements should be built from
  the ground up to look like, as much as possible, a starmap and
  navigational system that a space pilot might use." Done: a design
  brief (colours, type, panel shapes, iconography) approved by Boss,
  then the base template and shared components restyled to it in both
  themes, with contrast and reduced motion kept.
  Research (2026-10-09, map-ui-and-frontend-libraries.md): reuse the
  rating chip and palette for the habitability tiers so the site has one
  status component. Add the measured token table
  (map-ui-and-frontend-libraries.md section 3.3) as the starting brief,
  the rule that pairing every colour with a symbol, word or shape is
  part of done, and "B612 Mono or the system monospace stack for
  readouts, tabular figures" (open question for Boss, default: self-host
  B612 Mono for numbers and designations only); `--line` is a new token
  beside `--border`. Selection and caution are both amber, by aviation
  custom and never in the same place (default: accept).
  Prerequisite: UX.42.

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
    Research (2026-10-09, map-ui-and-frontend-libraries.md): change
    "stamped with the galaxy seed" to "seed plus a per-wipe instance id"
    (system ids restart after a wipe, seeds can repeat); `galaxy_shape`
    has no `created_at` or instance column, so the item needs a small
    migration.

  - [ ] **UX.47 A bookmark manager: list, rename, sort, group and delete**
    Done: a Bookmarks page listing every bookmark with its kind and
    place, where each can be renamed, grouped, sorted, opened on its map
    and deleted, and the list can be exported and imported as a file.
    Research (2026-10-09, map-ui-and-frontend-libraries.md): add import
    validation (size, count, kind whitelist, same-origin relative `url`
    only; `bookmarks.js` `urlOf` returns `entry.url` unvalidated into an
    `href`) and a move up / move down path (WCAG 2.5.7); drag and drop
    optional.
    Prerequisite: UX.46.

  - [ ] **UX.48 Charted regions as bookmarks that frame and outline the region**
    Boss (2026-10-07 11:47Z): "Clusters of contiguous generated content
    gets added to a system bookmark which will zoom in as close as it
    can to see all sectors in the group and highlight the shape of the
    region specifically." Done: every contiguous group of charted
    sectors gets a system bookmark (kept up to date as sectors are
    charted) that opens the Galaxy Map fitted to the group with its
    outline highlighted.
    Research (2026-10-09, map-ui-and-frontend-libraries.md): regions are
    server-derived (`GET /api/galaxy/regions`), not localStorage;
    contiguity is face adjacency (open question for Boss, default:
    sectors touching only along an edge or corner are not contiguous);
    the outline is the coplanar-merged edge set.
    Prerequisite: UX.47.

- [ ] **UX.78 Unit preference: Automatic, Metric only or Customary**
  Stored in `localStorage` key `planetgen.units` and applied by
  re-rendering `.qty` spans from `data-si`; a header menu on
  `sl-dropdown`. UX.42 is wording only; do not fold units into it. Open
  question for Boss (default: the three presets, per browser).
  Prerequisite: UX.23.
  Design: [docs/design/units-and-number-formatting.md](design/units-and-number-formatting.md)

- [ ] **UX.81 Time symbols Gyr, Myr, kyr in place of Gy, My, ky; AU from 1,000,000 km; scientific text below mantissa 1e-3**
  Rename before GEN.87 introduces the gray (touches `PERIOD_LADDER`,
  `format_age_string`, `period.js`, tests). Start the AU rung at
  1,000,000 km (0.0067 AU) so inner-system distances never print as
  "5.79 × 10⁷ km", and use scientific text below mantissa 1e-3 for every
  rung (pressure prints "0.00000000000025 Pa" today). Open questions for
  Boss (defaults yes): the rename and the AU start. Also pass
  `aria-valuetext` the unit name, not the symbol, in `facilityform.js`,
  and route direct `toFixed`/`toPrecision` calls in `systemview3d.js`,
  `galaxysystem.js`, `galaxymap3d.js`, `galaxystages.js` and
  `galaxystageview.js` through the ladders.
  Prerequisites: none.
  Design: [docs/design/units-and-number-formatting.md](design/units-and-number-formatting.md)

- [ ] **UX.82 Theme checks after PR #800: SVG currentColor, two Shoelace contrast failures, alpha in --bg-subtle**
  Carry into UX.43 or a small follow-up: a check that every file under
  `vendor/shoelace/assets` and `icons.svg` uses `currentColor`; recheck
  the two light-theme Shoelace failures measured before #800 (primary
  hover 4.23, placeholder 4.41, both below 4.5:1); `--bg-subtle` is
  still a colour with alpha in the dark theme. UX.76's cause was
  `a.btn:visited` out-specifying the secondary-button rules.
  Prerequisites: none.
  Design: [docs/design/map-ui-and-frontend-libraries.md](design/map-ui-and-frontend-libraries.md)

- [ ] **UX.87 The system list shows uncharted systems: every scattered star, with its location and a way to generate it**
  Boss (GitHub issue
  [#929](https://github.com/dwhagar/planetGen/issues/929), 2026-10-10
  01:51Z): "Every star that is scattered throughout in the brightness
  scatter needs to be also listed or able to be listed as 'uncharted' in
  the star system list. Information about that star and its location is
  displayed, its coordinates, sector coordinates, and other information
  including layer, shell, and slot that it occupies. This interface
  should also allow the user to generate that star system by itself,
  though the system will recommend generating the entire sector." Done:
  the system list has an 'uncharted' filter (off by default) that lists
  scattered stars with their coordinates, sector coordinates, layer,
  shell and slot and the star's own data; each row has a Generate button
  for that one system, with a note recommending the whole sector. Decided (Boss, 2026-10-10 04:06Z, "your defaults are confirmed"): scattered stars from the mass-limit and
  luminosity passes both count; generating one system fills only that
  system and leaves the sector's other contents ungenerated.
  Prerequisites: none. Related: MAP.162, ADM.32, ADM.35, NAV.48, DOC.9.

## MAP: Galaxy Map, Sector Map, System Map

MAP.2 with MAP.22 and MAP.23, MAP.15 and MAP.30 shipped in PR #234.

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
  the rule that the mini map is never zoomable (MAP.58, superseded by MAP.125).
  Design: [docs/design/galaxy-drilldown-navigation.md](design/galaxy-drilldown-navigation.md), section 15

  - [ ] **MAP.75 The mini map as a second engine view**
    The mini map is a second, locked camera on the same scene data;
    built on MAP.61's controller with a "locked" policy rather than its
    own renderer (one WebGL context, scissored like the System Map's
    sphere overlay).
    Research (2026-10-09, map-ui-and-frontend-libraries.md): no change
    to the decision (scissored second view); record the 2D fallback, the
    dirty-flag and `camera.layers` rules.

- [ ] **MAP.119 Expected star density editable by admins on the Galaxy Map**
  Boss (2026-10-03 05:38Z): "In the galaxy view, the expected star
  density should be an editable field for admin and a static field for
  users/guests.  Once edited it changes it for the current view (i.e. an
  entire slab or entire block)." Done: the info panel shows expected
  density; admins can edit it for the slab or block in view, stored as
  an override that later fills use and logged.
  Research (2026-10-09, map-ui-and-frontend-libraries.md): adopt design
  doc section 8: a multiplier, not E; an append-only `density_overrides`
  table; `density_multiplier_used` on each filled sector; range 0.1 to
  10; an activity-log entry; live preview, confirm and undo.

- [ ] **MAP.120 Bright-star backfill from the Galaxy Map's block, slab and wedge menus**
  Boss (2026-10-03 05:38Z): "Add bright star backfill to block menus in
  galactic view." Done: the admin menu on a block, slab or wedge runs
  the bright-star backfill for it as a job.

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
  Research (2026-10-09, course-avoidance.md): draw the course on a layer
  that is not faded with the sectors in front of the camera.
  Folded (2026-10-09, fly-through-view-distance.md): the blocker fade
  and the faint context around the focus are folded into MAP.149 (the
  near field), which is built as part of the fly-through (MAP.146).

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
  Research (2026-10-09, map-ui-and-frontend-libraries.md): specify a
  radiogroup control, the URL parameter `select=stars`, 10 px mouse and
  22 px touch reach with a disambiguation list, and a boundary filter
  (MAP.141). Picking cost measured: 5.4 ms at 70,000 stars with a
  typed-array loop; the BVH is not needed (remove three-mesh-bvh
  wording, see the stale-text item).

- [ ] **MAP.132 Overlay markers for black holes, nebulae and habitable worlds**
  Boss (2026-10-08 01:59Z): "Density determines transparency with more solid meaning more dense. Stellar age measures hue. Brightness measures luminosity of the sector or block. Blocks will never be fully opaque for coloring that is the color of the translucent fill."
  Later. The sector and block fills carry only age, density and
  luminosity, so the notable things inside get a small marker: black
  holes, nebulae and habitable worlds (by the habitability index, GEN.84).
  Done: markers show at the zooms where MAP.116's budget allows them,
  each kind can be turned off (MAP.123), and a browser test finds one of
  each on a seeded sector.
  Research (2026-10-09, anomalies.md): its list includes the new kinds
  (magnetar glyph, Wolf-Rayet bubble); keep the per-kind toggle
  (MAP.123) and the MAP.116 budget; nebula markers cover the far zoom,
  sprites under 2 px stay hidden.

- [ ] **MAP.139 The Galaxy View uses its spare space: an info box with a Details link, and menu items**
  Boss (GitHub issues [#758](https://github.com/dwhagar/planetGen/issues/758) and [#715](https://github.com/dwhagar/planetGen/issues/715)): "there's a lot of wasted space
  that I want to fill up with something useful, maybe menu items? Maybe
  information about what is selected?" and "A simple information box
  should appear overlaid on the map when there is room that displays the
  basics. Then a details link will jump you to the part of the page the
  details are listed on." Done: on a wide Galaxy View the space beside
  the map holds the selection's basics and the map's Menu items; where
  there is room on the map itself a small info box shows the basics with
  a Details link that scrolls to the full details; it builds on MAP.65's
  shared panel.

- [ ] **MAP.140 Double-click on a selected object goes there and opens its information**
  Boss (GitHub issue [#714](https://github.com/dwhagar/planetGen/issues/714), 2026-10-09 00:26Z): "Double clicking on a
  selected object goes there and pulls up its information." Done: on
  every map a double-click on the selected object moves the view to it
  (the next stage down, or the object's page) and opens its information;
  the single click keeps selecting; a browser test covers it on the
  Galaxy, Sector and System views.
  Research (2026-10-09, map-ui-and-frontend-libraries.md): specify the
  pointerup-based recogniser in `createPointerControl` (350 ms, 24 px),
  the 350 ms deferral of the MAP.113 deselect (open question for Boss,
  default yes), and a non-gesture "Go to" alternative (MAP.139's button
  beside the Details link).
  Extended (2026-10-09, fly-through-view-distance.md): the go-to flight
  is extended by MAP.150 (the free camera): double-click flies to the
  thing clicked, with the wheel zooming to the cursor.

- [ ] **MAP.141 Context around the selection: faint neighbours, and the sectors above and below**
  Boss (GitHub issues [#718](https://github.com/dwhagar/planetGen/issues/718) and [#716](https://github.com/dwhagar/planetGen/issues/716)): "When selecting a slab we should
  still be able to see what is behind it and around it faintly just like
  in the other views" and "When viewing a sector you should also be able
  to see what is above and below it in the same style as what is around
  it, and the ability to select the sector below it. Also remember not
  to put sector clicking on the absolute boundaries of the galaxy."
  Done: while a slab, block or sector has focus, the objects around it
  are drawn faintly in every direction, including the sectors above and
  below, which can be selected; no click lands on the galaxy's absolute
  boundary. Issue #716 was labelled bug but is a new view, so it sits
  here.
  Research (2026-10-09, map-ui-and-frontend-libraries.md): use one
  shared `contextOpacity(region, focus, camera)` with MAP.121; the
  blocker fade formula is in design doc section 5.4.
  Folded (2026-10-09, fly-through-view-distance.md): the faint context
  around the selection is folded into MAP.149 (the near field), built as
  part of the fly-through (MAP.146).

- [ ] **MAP.142 Nebulae have fuzzy, fading boundaries**
  Boss (GitHub issue [#713](https://github.com/dwhagar/planetGen/issues/713), 2026-10-09 00:25Z): "Can we make the nebula
  boundaries fuzzy and kind of fade a bit so it looks more gaseous?"
  Done: nebulae drawn from their shape (GEN.75) fade out toward the edge
  with a soft falloff instead of a hard outline, in the Galaxy, Sector
  and nebula views, and the falloff keeps the nebula's extent readable.
  Research (2026-10-09, nebula-and-asteroid-field-classes.md): use the
  per-family edge widths and shell profiles in "Appearance by class":
  nested iso-shells with edge fade now, a bounded ray-march of a served
  3D field for the focused nebula. Done criteria add: the extent is the
  50% contour (open question for Boss, default hard extent at today's
  outline); dark nebulae get an outline reaching 3:1 (they composite at
  1.17:1 today); classes are distinguished by a non-colour cue.

- [ ] **MAP.143 Color sectors by their number of habitable locations**
  Boss (GitHub issue [#717](https://github.com/dwhagar/planetGen/issues/717), 2026-10-09 00:44Z): "Coloring option
  available for sectors based on the number of habitable locations."
  Done: the Galaxy Map's Color by switch (MAP.131) gains a Habitable
  worlds choice that colours each sector and block by how many planets
  and moons score as habitable, with a legend; it needs the habitability
  score of GEN.89 and the sector stats to hold the count.
  Research (2026-10-09, map-ui-and-frontend-libraries.md): a single-hue
  sequential ramp that survives greyscale, 4 to 5 bins with a legend,
  zero habitable worlds drawn with no fill; blocks never fully opaque
  (Boss).

- [ ] **MAP.145 A sky and Galaxy Map drawing rule for neighbour galaxies**
  An extended sprite with computed magnitude and size and a
  surface-brightness cut at about 23 mag/arcsec^2; the Magellanic Clouds
  are recomputed per observer, M31 and farther need only one stored
  direction.
  Prerequisite: GEN.157.
  Design: [docs/design/multiple-galaxies.md](design/multiple-galaxies.md)

- [ ] **MAP.146 Fly through the galaxy: scroll-zoom, double-click flight, distance-based visibility and a see-through near field**
  Boss (2026-10-09 22:24Z): "I want the zoom drill down to be less
  specific, instead of set wedges, blocks, and slabs pre-determined, have
  them based on the center of where the cursor is clicked. So that we
  always are drilling down exactly where the user wants. We'll have to
  convert block/slab measures to ranges of layers/shells/slots for
  filling on demand." Boss (22:49Z and 22:56Z): the end goal is to
  scroll-zoom and fly smoothly through the interface to a location; the
  system does not do view distance well and renders and keeps clickable
  what is right in front of the camera and the object the camera is in;
  stars should be shown by a smooth distance-and-brightness gradient so
  only what the observer could see is drawn; inside a block or sector the
  observer should see its contents without looking around what is in the
  way, the nearer things becoming more transparent as they approach;
  double-click zooms or flies to a clickable thing, the wheel zooms in and
  out, the user is always in the 3D galaxy even when viewing a sector;
  the sectors looked at are rendered well and the surroundings stay
  visible but quiet; smooth from the galaxy down to a star system.
  Done: the Galaxy Map has one free camera from the whole galaxy to a
  star system, built by its five sub-items: MAP.148 (the star visibility
  law), MAP.149 (the near field), MAP.150 (the free camera), MAP.151
  (the region data layer) and MAP.152 (scale hand-offs). The wheel zooms
  toward the point under the cursor, double-click flies to the thing
  clicked, stars fade by apparent magnitude against an on-screen limit,
  things near the camera dissolve, and the camera position names the
  container. The arc, slab and segment picks stop being the way to move;
  old stage URLs stop working (no backward compatibility). The open
  questions for Boss, each with a default, are on the sub-items.
  Phase 1, and the headline of the 8.1 release (Boss, 2026-10-09
  23:07Z): "Let's inject into phase 1 to full build ... This will be our
  major version bump to 8.1 later when we finish Phase 1, this is the
  main feature to move that." So MAP.146 and every item under it
  (MAP.147 to MAP.155) are part of what "Phase 1 complete" means: the
  8.1 stamp waits until all of them are done, and the 8.0 hold stays
  until then. The open questions on the sub-items stand, each with its
  default.
  Order (Foundations lane 2, the map engine lane): MAP.157 and MAP.158
  (the small wire format steps, no prerequisites), then MAP.153 and
  MAP.149 (client side only, no prerequisites), then MAP.148 (needs
  MAP.153), MAP.150 (needs MAP.149) and MAP.155 (needs MAP.153), then
  MAP.154 (needs MAP.153), MAP.159 (packed binary tiles; needs MAP.154
  and MAP.158; one cache stamp bump with MAP.154 and MAP.151's tile
  keys), MAP.151 (needs ADM.29 from Foundations lane 1) and MAP.152 (needs MAP.148,
  MAP.150 and MAP.154); this umbrella closes last.
  Prerequisites: MAP.148, MAP.149, MAP.150, MAP.151, MAP.152, MAP.153,
  MAP.154, MAP.155, MAP.157, MAP.158, MAP.159. Related:
  MAP.120, MAP.121, MAP.141, MAP.140, MAP.59, MAP.116, MAP.122, MAP.125,
  MAP.131, MAP.134, MAP.147, ADM.29, ADM.30, GEN.101, GEN.126, NAV.13,
  NAV.14.
  Overlap (2026-10-09, zoom-star-visibility.md): Star visibility while
  zooming is also covered by MAP.153 to MAP.155
  (docs/design/zoom-star-visibility.md): MAP.153 is the first client
  stage of MAP.148's law.
  Design: [docs/design/fly-through-view-distance.md](design/fly-through-view-distance.md)
  Design: [docs/design/drilldown-region-sizes.md](design/drilldown-region-sizes.md)
  Design: [docs/design/galaxy-drilldown-navigation.md](design/galaxy-drilldown-navigation.md)

- [ ] **MAP.147 The Galaxy Map wire format: the investigation (done) and the record of what was built from it**
  Boss (2026-10-09 22:41Z): "File section 5.6 as an item please, and
  start an investigation thread to measure current payload and compare
  options please, use research lane 3 after it's current research is
  done." Investigation finished by Research Lane 3 (2026-10-09, report
  docs/design/galaxy-map-wire-format.md). Measured facts: on the wire a
  star is about 270 bytes of JSON (45 gzipped), so a full 70,000-star
  view (MAP.109) is about 19 MB raw and 3.2 MB gzipped, not the 3.6 MB
  section 5.6 assumed; the 52 bytes a star in section 5.6 is the GPU
  buffer, so quantising it would not change the download; the opening
  view pre-fetches 28 tiles (533 KB gzipped), about 18 times what it
  shows; `placed`, `planned` and `filled` are 9 to 35% of every tile and
  are never read; localStorage caps at about 5 MB, so revisits mostly
  miss; the first star waits on the 173 static files, not on a tile.
  Recommendation, now built as separate items: (1) MAP.157 trims the JSON
  (-48% gzipped, no client change) and serves prebuilt, brotli
  precompressed bytes; MAP.158 makes the prefetch gentler and moves the
  tile cache to IndexedDB; (2) MAP.159 packs the tiles as binary (12.9 to
  16.7 bytes a star, -78% gzipped) together with the nested lists of
  MAP.154 and MAP.151's tile keys, on one cache stamp bump; (3) MAP.160
  defers quantising the GPU buffers. This item stays open as the record of
  the investigation until MAP.157 to MAP.159 are done.
  Open question for Boss (default each sub-item goes ahead as written;
  MAP.160 is deferred): other?
  Prerequisites: MAP.157, MAP.158, MAP.159.
  Linked (2026-10-09, fly-through-view-distance.md): the tile and stage
  cache keys follow MAP.151 (the region data layer), so decide the wire
  format and those keys together; the star visibility law (MAP.148) sets
  how many stars a view needs.
  Design: [docs/design/map-ui-and-frontend-libraries.md](design/map-ui-and-frontend-libraries.md)

- [ ] **MAP.148 The star visibility law: apparent-magnitude opacity, flux-based brightness and an on-screen limit from a histogram**
  Source: docs/design/fly-through-view-distance.md (section 7 and 8),
  written at Boss's request of 2026-10-09 22:56Z. Open questions are for
  Boss; defaults stand until he answers. Nothing is built until he
  decides.
  Done: a star's opacity follows its apparent magnitude from the camera,
  against a limit chosen so about 20,000 stars are on screen, found from
  a histogram of the stars in view. Brightness is by flux. This replaces
  the step floors in `generated_star_floor_sol` as the rule and folds
  MAP.116's budget table into it. First on the tiles already fetched
  (client side, no schema change), then tiles chosen by distance (with
  MAP.152).
  Open question for Boss (default 20,000 stars on screen, 8,000 on a
  phone, with a 1.5 magnitude ramp, tuned after a first build): other
  numbers?
  Overlap (2026-10-09, zoom-star-visibility.md): MAP.153 is the first
  client stage of this same visibility rule
  (docs/design/zoom-star-visibility.md, Research Lane 1): a rank birth
  radius on the tiles as fetched today, with the apparent-magnitude law
  here as the end state; the two are one rule in two stages, not
  competitors. MAP.154 and MAP.155 carry the server list nesting and the
  other objects.
  Dependency (2026-10-09, fly-through-view-distance.md): The law
  multiplies MAP.153's rank birth radius in one shader: a = a_rank(R) *
  a_mag(d) * a_near, built after MAP.153. Calibrate m_lim so the
  magnitude factor is about 1 for a star at the target distance at any
  camera radius; it only dims stars much farther than the target, so the
  two rules never thin the same stars twice. The distance-cut tiles it
  leads to also need MAP.154 (nested lists).
  Prerequisite: MAP.153. Related: MAP.116, MAP.146, MAP.147.
  Design: [docs/design/fly-through-view-distance.md](design/fly-through-view-distance.md)

- [ ] **MAP.149 The near field: depth fade, a see-through focus tube, drawing from inside a container, and picking that matches what is drawn**
  Source: docs/design/fly-through-view-distance.md (section 7 and 8),
  written at Boss's request of 2026-10-09 22:56Z. Open questions are for
  Boss; defaults stand until he answers. Nothing is built until he
  decides.
  Done: things nearer the camera than a fraction of the focus distance
  dissolve; a soft see-through tube thins what stands between the camera
  and the focus; the container the camera is in is drawn from the
  inside; only what is visible enough can be picked, so picking agrees
  with drawing. One shared function does the fade for drawing and
  picking. The region looked at and the container are drawn at full
  strength and the rest faintly. Folds in MAP.121's blocker fade and
  MAP.141's faint context. First client-side, no schema change; it
  improves today's map.
  Open question for Boss (default: context regions at opacity 0.08 to
  0.3 and 2.5 magnitudes shallower than the focus): other strengths?
  Prerequisites: none. Related: MAP.121, MAP.141, MAP.146.
  Design: [docs/design/fly-through-view-distance.md](design/fly-through-view-distance.md)

- [ ] **MAP.150 The free camera: wheel zoom to the cursor, double-click flight, and the observer inside, with the container named from position**
  Source: docs/design/fly-through-view-distance.md (section 7 and 8),
  written at Boss's request of 2026-10-09 22:56Z. Open questions are for
  Boss; defaults stand until he answers. Nothing is built until he
  decides.
  Done: one free camera from the whole galaxy to a star system. The
  wheel zooms toward the point under the cursor with clearance from what
  is in front; double-click flies to the thing clicked (MAP.140's go-to
  extended); the observer is always in the 3D galaxy, also when looking
  at a sector; the container is named from the camera position and the
  breadcrumb is derived from it; URLs and bookmarks hold the camera. The
  arc, slab and segment picks (MAP.85, MAP.56 and MAP.17 flow) stop
  being the way to move, and old stage URLs stop working (no backward
  compatibility).
  Open question for Boss (default yes): retire the arc, slab and segment
  picks as the navigation flow?
  Open question for Boss (default: keep the slab strip as an optional
  section plane, not a stage): or drop it?
  Prerequisite: MAP.149. Related: MAP.140, MAP.85, MAP.59, NAV.13,
  NAV.14, MAP.146.
  Dependency (2026-10-09, fly-through-view-distance.md): If the camera
  radius R used by MAP.153's rank rule is redefined for a free camera
  (for example distance to the nearest sector instead of to the target),
  keep it continuous in the camera position: a jump in R is a pop for
  every star at once.
  Design: [docs/design/fly-through-view-distance.md](design/fly-through-view-distance.md)

- [ ] **MAP.151 The region data layer: exact-centred frame, aligned cells, per-level aggregates and slot-wrap ranges**
  Source: docs/design/fly-through-view-distance.md (section 7 and 8),
  written at Boss's request of 2026-10-09 22:56Z. Open questions are for
  Boss; defaults stand until he answers. Nothing is built until he
  decides.
  Done: region data is read in aligned cells with an exact-centred
  frame, with per-level aggregates, and a region is described as layer,
  ring and slot-arc ranges (with slot wrap) that stats, ADM.29 and
  ADM.30 fills, MAP.120 backfills and fill on demand all use. Tile and
  stage cache keys follow the cells and bump the cache stamp (clear
  /var/cache/planetgen/tiles on update); decide them together with
  MAP.147's wire format. This is the data layer of the original MAP.146
  text (see docs/design/drilldown-region-sizes.md for region sizes,
  shapes and the cell pyramid). It touches the per-sector stats and
  colouring (MAP.131), the Galaxy Map opening view (MAP.134) and the
  settle step (GEN.126), and `galaxy/drill.py` and `galaxyprisms.js` are
  replaced, not kept beside it.
  Open question for Boss (default: a 3-ary pyramid accepting cells 0.84
  to 1.25 of an edge across): or keep 9-ary and accept a 2.45x gap in
  sizes?
  Prerequisite: ADM.29. Related: MAP.120, MAP.147, ADM.30, GEN.101,
  GEN.126, MAP.131, MAP.134, MAP.146.
  Dependency (2026-10-09, fly-through-view-distance.md): The region data
  layer and its aggregates (3-ary pyramid, per-star id, per-level
  aggregates) wait on MAP.147 (wire format, Research Lane 3).
  Wire format (2026-10-09, galaxy-map-wire-format.md): The wire format
  decision is made (MAP.147's report): the tile and stage cache keys are
  decided here together with MAP.159 (packed binary tiles) and MAP.154,
  on one cache stamp bump.
  Design: [docs/design/drilldown-region-sizes.md](design/drilldown-region-sizes.md)
  ADM.29 is merged (PR #926): `planetgen/galaxy/span.py` (`Span`,
  `parse_range`: layer, ring and slot-arc ranges with wrap and an
  instant count from prefix sums). This item is unblocked on the ADM.29
  side and reuses `span.py`; the wire-format steps it still waits on
  are in its other notes (MAP.157 to MAP.159).
  Note (2026-10-09): From ADM.31 (PR #942, 2026-10-10): the 'made' rings
  on the Galaxy Map cover only the first 2,000 sectors of a run; the
  region data layer should replace that cap with ranges.

- [ ] **MAP.152 Scale hand-offs: galaxy, sector and system cross-fade with hysteresis, and per-tile camera-relative origins**
  Source: docs/design/fly-through-view-distance.md (section 7 and 8),
  written at Boss's request of 2026-10-09 22:56Z. Open questions are for
  Boss; defaults stand until he answers. Nothing is built until he
  decides.
  Done: the galaxy, sector and system scales cross-fade with hysteresis
  so the camera never flickers between them, and each tile has its own
  camera-relative origin below about 100 pc so stars stay precise at
  close range. Includes the distance-cut tile choice that completes the
  visibility law.
  Dependency (2026-10-09, fly-through-view-distance.md): Distance-cut
  tiles need MAP.154 (nested server lists: a child tile must contain the
  parent's stars in its box).
  Prerequisites: MAP.148, MAP.150, MAP.154. Related: MAP.102, MAP.125,
  MAP.146.
  Design: [docs/design/fly-through-view-distance.md](design/fly-through-view-distance.md)

- [ ] **MAP.153 Stars fade in with the zoom: a birth radius from each star's rank in its tile list (first client stage)**
  Source: docs/design/zoom-star-visibility.md (section 3.1 and 5, stage
  1), written at Boss's request of 2026-10-09 22:30Z: "I want a very
  smooth transition where stars and objects are slowly added as one
  zooms in." Today the drawn set is frozen between 11 tile-level radii
  and then up to 8.8 times as many stars appear in one frame. Done: in
  `galaxymap3d.js` every star gets a birth radius from its rank in its
  tile's list (most luminous first, as the server already sorts it), R_b
  = R* 2^W (N0/r)^(1/3), and its opacity is a smoothstep of the camera
  radius, cross-faded from the parent tile's rank across the tile
  level's octave; zooming out removes stars as smoothly as zooming in
  adds them, a late tile changes nothing visible, and a browser-test
  hook returns the opacity sum at a given radius so a test bounds the
  step. Measured on the same tile data the worst single 9% step falls
  from +883% to +59% (dense) and from +775% to +74% (thin). No server or
  schema change. This is the FIRST STAGE of the visibility rule that
  MAP.148 ends with: both decide when a star shows while zooming, so
  they are not two rules. Until MAP.148 lands, the rank rule leads;
  MAP.148 then replaces the rank with apparent magnitude (rank stays as
  the cap and tiebreak).
  Decided (Boss, 2026-10-09 23:11Z, "default options are approved"):
  the rank birth-radius fade is stage 1 and the apparent-magnitude law
  of MAP.148 is the end state.
  Decided (Boss, 2026-10-09 23:11Z, "default options are approved"): each
  star fades over one halving of the camera radius (W = 1); a dense
  sector's stars arrive in rank order between about 35 pc and 8 pc of
  view radius, so they wait for sector zoom; the GEN.30 backfill shells
  stay as they are.
  Prerequisites: none. Related: MAP.148, MAP.146, MAP.149, MAP.147,
  MAP.116.
  Design: [docs/design/zoom-star-visibility.md](design/zoom-star-visibility.md)

- [ ] **MAP.154 Nested bright-star lists on the server, so every parent list is a subset of its child's**
  Source: docs/design/zoom-star-visibility.md (3 scheme J, section 5
  stage 2). Done: the bright lists use one key for every tile level (a
  population weight in the key instead of equal-share picking) and the
  per-level budgets never fall, so a child list always holds the
  parent's stars inside its box and the rank fade of the previous item
  is exact instead of degrading gracefully. This changes the tile cache
  stamp. The key and list shapes are decided together with MAP.147 (wire
  format) and MAP.151 (region data layer), and with MAP.148 if the
  magnitude law changes what a tile lists. Open question for Boss
  (default build it only after MAP.153 is seen working): go ahead?
  Wire format (2026-10-09, galaxy-map-wire-format.md): MAP.147's
  investigation is done (docs/design/galaxy-map-wire-format.md): the
  packed binary tile format ships with these nested lists as MAP.159, on
  one cache stamp bump, so list order is rank and finer tiles omit stars
  a coarser tile sent.
  Prerequisite: MAP.153. Related: MAP.147, MAP.148, MAP.151.
  Design: [docs/design/zoom-star-visibility.md](design/zoom-star-visibility.md)

- [ ] **MAP.155 Other objects fade in too: point objects from level 8, a size ramp for cloud sprites, and stars that grow from a faint dot**
  Source: docs/design/zoom-star-visibility.md (section 5 stage 3,
  schemes I and E). Done: black holes, neutron stars and quasars are
  listed from tile level 8 (160 pc) instead of 10 (40 pc) and fade in by
  the same rank rule, cloud sprites ramp opacity between 1 and 4 px
  instead of switching on at 2 px (`CLOUD_MIN_PX`), and a new star
  starts one pixel and faint and grows to its size with its opacity.
  Rows with no natural rank use a hash tiebreak (a fixed random number
  from the object id). Sector blocks and fills keep their own
  level-of-detail question (the mega-block plan). Decided (Boss, 2026-10-09 23:11Z, "default options are
  approved"): point objects are listed from level 8.
  Prerequisite: MAP.153. Related: MAP.153, MAP.148, MAP.149.
  Design: [docs/design/zoom-star-visibility.md](design/zoom-star-visibility.md)

- [ ] **MAP.156 View one layer or a range of layers top-down from the galaxy view, as a secondary option**
  Boss (2026-10-09 23:17Z): "I should be able to from the galaxy view
  select a single layer (or range of layers) to view top-down. It should
  not be a main option." Done: from the Galaxy Map the user can pick one
  layer or a range of layers (layers are the vertical slices of the
  galaxy, the same layers the Generate page fills) and see only those
  looking straight down on the galaxy plane. The choice sits in a
  secondary place (a menu entry or the advanced section of the control
  panel), never as a main control or a step the user must take to move
  around, and it is off by default so the free camera of MAP.146 is what
  opens. Leaving it restores the previous view. Decided (Boss, 2026-10-09 23:18Z, "default is approved"): a two-handle
  layer range in the "Show" menu, top-down camera with the same
  vertical-fade rules as the rest of the map, stars and sectors outside
  the range hidden.
  Prerequisites: none. Related: MAP.146, MAP.150, MAP.141, MAP.122,
  ADM.29.

- [ ] **MAP.157 Trim the Galaxy Map tile JSON and serve it from prebuilt, precompressed bytes**
  Source: docs/design/galaxy-map-wire-format.md (section 1, steps 1 and
  2; MAP.147's recommendation). Measured by Research Lane 3
  (2026-10-09): a star is about 270 bytes of JSON on the wire (45
  gzipped); the `placed`, `planned` and `filled` sections are 9 to 35%
  of every tile and the client never reads them. Done: the tile JSON
  drops the sections and star fields the client never reads, writes
  `star_type` as its letter code only and rounds the floats; the cache
  stores the finished bytes (not parsed dicts that are serialised again
  on every hit) and a brotli copy made at cache-write time, which Apache
  serves. Measured on 400 real tiles the gzipped size falls 48% (3.52 MB
  to 1.83 MB) with no change to what the page draws; a warm hit falls
  from 22-75 ms to about 2 ms. No client change, so it can land at once
  and does not wait for the fly-through work; it changes the tile cache
  stamp. Open question for Boss (default yes, do it now and first in the
  fly-through step): go ahead?
  Prerequisites: none. Related: MAP.147, MAP.158, MAP.159, PERF.38,
  PERF.41, MAP.109.
  Detail (2026-10-09, galaxy-map-wire-format.md): detail from Research
  Lane 3 (Boss said yes, 2026-10-09 23:26Z): drop `placed`, `planned`
  and `filled` and the star fields `ring_index`, `layer_index`,
  `ring_slot_index`, `population`, `yerkes_class` and the generated
  name; `star_type` to its class letter; round x/y/z to 3 decimals,
  luminosity to 4 and radius to 3 significant digits, temperature to 10
  K. Re-grep `static/` for readers before removing any field. Fix the
  `GALAXY_VIEW_MAX_STARS` docstring (270 bytes a star, not about 114).
  Serve stored response bytes in `fetch_tiles` and add `mod_brotli` (or
  serve the `.br` copy) to `examples/apache`. Lands after or with
  PERF.38.
  Design: [docs/design/galaxy-map-wire-format.md](design/galaxy-map-wire-format.md)

- [ ] **MAP.158 A gentler tile prefetch and an IndexedDB tile cache instead of localStorage**
  Source: docs/design/galaxy-map-wire-format.md (sections 2.3, 4.4, 1
  step 2 and 5). Measured: the opening view downloads 28 tiles (533 KB
  gzipped) to show one, about 18 times what is on screen; localStorage
  holds about 5 MB, which one sector link fills (40 tiles), so a revisit
  mostly misses. Done: the page prefetches only the next zoom step in,
  only after the view has been idle, and not when the browser asks to
  save data; tiles are cached in IndexedDB (or by immutable HTTP caching
  per tile) and a revisit finds them; the old localStorage cache is
  removed (no compatibility). Open question for Boss (default yes): go
  ahead?
  Detail (2026-10-09, galaxy-map-wire-format.md): detail from Research
  Lane 3: prefetch in the zoom-in direction only, after idle, skipped on
  `navigator.connection.saveData` or a slow `effectiveType`, with a byte
  cap (`prefetchTiles` in `galaxymap3d.js`). Interim cheap fix for the
  cache before IndexedDB: store only the tiles the next view needs and
  cap by bytes with least-recently-used removal, instead of `storeTile`
  wiping every stored tile when the quota fails.
  Prerequisites: none. Related: MAP.147, MAP.157, MAP.109.
  Design: [docs/design/galaxy-map-wire-format.md](design/galaxy-map-wire-format.md)

- [ ] **MAP.159 Packed binary Galaxy Map tiles (quantised planes) with the nested tile lists, on one cache stamp bump**
  Source: docs/design/galaxy-map-wire-format.md (section 1 step 2, 4.1,
  4.3). Measured: a packed binary tile is 12.9 to 16.7 bytes a star, 78%
  fewer gzipped bytes than today's JSON (3.52 MB to 0.79 MB over 400
  real tiles; the worst 70,000-star view falls from 3.2 MB to about 0.8
  MB). Done: the server writes tiles as packed typed arrays (quantised
  position, colour and scalar planes; `clouds` and `points` stay JSON),
  the client decodes them in `tileStars`, list order is rank so the rank
  fade of MAP.153 reads it directly, and a finer tile omits the stars a
  coarser tile already sent (31 to 43% of bright records at levels 3 to
  6, all at 7 to 12). It ships together with MAP.154 (the nested lists)
  and MAP.151's tile keys so the cache is invalidated once. Open
  question for Boss (default yes, with MAP.154 and MAP.151 on one bump):
  go ahead?
  Prerequisites: MAP.154, MAP.158. Related: MAP.147, MAP.154, MAP.151,
  MAP.153, MAP.157.
  Detail (2026-10-09, galaxy-map-wire-format.md): detail from Research
  Lane 3 and Research Lane 2: x/y/z uint16 inside the tile cube, log
  luminosity uint16, log temperature and radius uint8, class and flags
  uint8, id uint32 (a lean record without ids is 11 bytes); `clouds` and
  `points` stay JSON; the 80-bit object ID stays off the tile wire. Keep
  list order, so rank is implied (use the unsorted 16.7 byte form, or
  add a 2 byte rank plane). The format fits any octree edge from 16 pc
  to 65,536 pc; aggregates cost about 10 bytes a cell; block responses
  are tiny, so changing block keys costs nothing on the wire. Depends on
  MAP.157 and MAP.158 and on Boss's yes on the design.
  Design: [docs/design/galaxy-map-wire-format.md](design/galaxy-map-wire-format.md)

- [ ] **MAP.160 Quantise the Galaxy Map GPU buffers (deferred)**
  Source: docs/design/galaxy-map-wire-format.md (section 5, from section
  5.6 of map-ui-and-frontend-libraries.md). The 52 bytes a star in
  section 5.6 is the GPU buffer, not the download, so quantising it (to
  about 20 to 24 bytes) saves GPU memory and the re-upload on each tile
  arrival, not network bytes; not measured on a real GPU. Done: colour,
  scalars and flags in `setStars` are byte-quantised, and a client that
  appends the new tile's stars instead of rebuilding every star is
  considered with it. Open question for Boss (default defer until a real
  phone or GPU measurement shows the upload matters): build it?
  Detail (2026-10-09, galaxy-map-wire-format.md): also considered here:
  the client appends the new tile's stars instead of rebuilding and
  re-uploading every star on each tile arrival (3,200 stars is a 166 KB
  upload, 11,000 is 572 KB each time), and `tileStars` copies an object
  per generated star. Optional.
  Prerequisites: none. Related: MAP.147, MAP.159, MAP.109.
  Design: [docs/design/galaxy-map-wire-format.md](design/galaxy-map-wire-format.md)

- [ ] **MAP.161 Load the Galaxy Map faster on a first visit: bundle or preload its scripts**
  Source: docs/design/galaxy-map-wire-format.md (sections 2.3 and 1).
  Measured by Research Lane 3: the first star does not wait on a tile
  (the opening tile is inside the 35 KB gzipped page); it waits on 173
  static files (543 KB gzipped, 1.8 MB decoded; the map scripts and
  three.js), which is 7.9 to 8.1 s on a slow 4G link and 1.0 to 1.3 s
  locally. Production serves `/static` immutable, so repeat visits are
  fine. Done: the map scripts load faster on a first visit, by bundling
  them or `modulepreload` hints, and HTTP/2 in the Apache example
  (`examples/apache/planetgen.conf.example`); this keeps the no-bundler
  decision (vendored ES modules) unless Boss says otherwise, so the
  default is `modulepreload` plus HTTP/2. Phase 1 (Boss, 2026-10-09 23:29Z).
  Open question for Boss
  (default `modulepreload` and HTTP/2, no bundler): or bundle?
  Prerequisites: none. Related: MAP.147, MAP.157, MAP.158.
  Design: [docs/design/galaxy-map-wire-format.md](design/galaxy-map-wire-format.md)

- [ ] **MAP.162 A sector holding scattered objects but never generated can still be opened, marked uncharted**
  Boss (GitHub issue
  [#928](https://github.com/dwhagar/planetGen/issues/928), 2026-10-10
  01:44Z): "When a scatter places an object inside a sector that is
  otherwise ungenerated, the user should still be able to select that
  sector and view it, so they can view the placed item(s) within the
  sector. But the sector should have some indicator that the sector is
  uncharted." Done: on the Galaxy Map and the sector view, a sector with
  scattered stars, black holes or other phenomena but no generated
  contents opens, shows those objects, and carries a clear uncharted
  mark (in the title, the info panel and the sector view's frame);
  generating it removes the mark.
  Prerequisites: none. Related: MAP.122, MAP.120, UX.87, NAV.48.

## NAV: Navigation and courses

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
  and clickable. The same works inside a sector. Ties in with MAP.61 (a course view
  may need a fitted zoom outside the user range; MAP.58's zoom limits
  were superseded by MAP.125).
  Much of this exists: `/galaxy?course=<from>,<to>` (MAP.27, 7.52.0)
  draws one line through the route's stops with the ends named. Missing:
  the direct line drawn apart from the route, a fitted zoom, and any
  course inside a sector.
  Prerequisites: NAV.21, NAV.22, NAV.23.
  Plan (2026-10-07): Boss (2026-10-07): "Plotted courses should appear
  on the galactic map and stay until cleared." See NAV.49.

  - [ ] **NAV.21 Fit the view to the whole course**
    Done: the course opens zoomed in as far as possible with every
    point of both paths on screen, fitted to the window (the same fit as
    MAP.53), refitted on resize. This view is exempt from any zoom
    limit (MAP.58 was superseded by MAP.125); the user can still rotate it.
    Research (2026-10-09, course-avoidance.md): fit the view to the
    bounding box of the adjusted `polyline`, not of the waypoints.

  - [ ] **NAV.22 Courses inside a sector and a system**
    Today a course that stays in one sector opens the sector with no
    line. Done: the same course drawing at the sector stage (MAP.61)
    and in the 3D system view (MAP.62), and a course that spans levels
    shows its in-system legs when zoomed in.
    Research (2026-10-09, course-avoidance.md): define the scale
    hand-off: the galactic path ends at the system's star; the in-system
    leg starts where the path crosses the heliopause
    (`heliosphere_radius_km`) and ends at the target body; sector-scale
    courses use the same keep-out table in light-years.

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
  Research (2026-10-09, course-avoidance.md): inside a system a star's
  keep-out is its radiation radius `max(10 R_star, sqrt(L / (4 pi
  F_lim)))` with F_lim = 50 kW/m^2, and a neutron star or black hole
  inside a system uses its tide, field and radiation radius; the first
  open question is answered as already decided (galactic Hill radius for
  a lone star or object). Hops are bent one at a time; the route's stops
  do not change in practice (overhead 0.04% to 0.2% on average).

  - [ ] **NAV.25 Find the obstacles along a path**
    Done: from the corridor query of NAV.10 (and NAV.38's
    line-to-sectors helper), every
    object whose keep-out sphere comes within reach of the path, at
    each scale: between systems, inside a sector, inside a system.
    Research (2026-10-09, course-routing.md): point at the shared
    segment query (course-routing.md section 6); no separate work.
    Research (2026-10-09, course-avoidance.md): specify the query as
    "every object whose keep-out sphere is within R_i of the straight
    segment" in two parts: stars, neutron stars, rogue planets and small
    fields by the sectors the segment crosses plus a one-sector margin
    (their largest radius, 2.1 pc, is under one 4 pc sector edge); large
    kinds (intermediate black holes up to 25 pc, SMBH and quasar 2 to 10
    pc and up, nebulae up to 200 ly, remnants, large fields) by a
    bounding-box query on centre and radius, as
    `_placed_phenomenon_rows` does. Return `ref`, centre, radius, basis
    and note per object. An ungenerated sector has no objects to avoid;
    the answer says so (NAV.12's unknown-space flag already marks it).

  - [ ] **NAV.26 Bend the path around keep-out spheres**
    Done: a path planner that leaves the straight line only where it
    crosses a keep-out sphere, going around it on the shortest detour
    (a tangent arc, or waypoints just outside the sphere), checked
    again against the other spheres; the end bodies' own spheres are
    exempt (a ship has to enter them to arrive). Unit tests with
    hand-placed spheres.
    Asteroid fields count as keep-out areas too (Boss, 2026-10-09
    01:02Z: "Courses should avoid asteroid fields."), which NAV.51
    adds; nebulae and remnants stay pass-through with a note until Boss
    decides.
    Research (2026-10-09, course-avoidance.md): use the "insert wrap at
    first hit, then validated polish" method (section 5 of the design
    doc; the prototype needs only `math` and numpy). Exempt every sphere
    that contains an end, not only the end body's own. Add the overlap
    policy (cluster detection, enclosing sphere, flagged straight-line
    fallback `overlap_fallback`). Hand-placed tests: one sphere on the
    line (closed form), collinear start-centre-end, start inside a
    sphere, end on a surface, two spheres with a gap and with contact, a
    chain of wraps, 40 tiny spheres on a long line, an overlapping pair,
    a scene shifted by 2e4 pc (same answer). Say "short detour", not
    "shortest".

  - [ ] **NAV.27 Moving bodies inside a system**
    Inside a system the planets move. Done: the in-system planner uses
    positions at the time of travel (MAP.62's positions-at-time module,
    Python twin), from the departure time and the leg's speed. That
    answers the item's second open question: yes. Default departure:
    now.
    Research (2026-10-09, course-avoidance.md): warp and fold speeds
    make departure-time positions sufficient, with the target aimed at
    its arrival-time position (`t0 + L/v`, 2 to 3 fixed-point steps);
    add the 3-to-6-pass iteration only if sub-0.01c speeds are
    introduced. `body_positions.positions_at` is the position source; no
    new ephemeris code. Galactic-scale motion is negligible at warp, so
    this stays system-only.

  - [ ] **NAV.28 Show and save the adjusted course**
    Done: the adjusted path, its extra length and its times beside the
    straight line on the NAV page and on the map (NAV.5); saved courses
    (NAV.4) keep it as a third form.
    Research (2026-10-09, course-avoidance.md): use the `adjusted`
    object of design section 7 (exact legs, waypoints, polyline,
    `avoided`, `passed_through`, `ignored_below_resolution`, `ruleset`,
    `warnings`); keep `waypoints`, `legs` and the radii used in the
    saved record (NAV.17).

  - [ ] **NAV.51 Courses route around asteroid fields**
    Boss (2026-10-09 01:02Z): "Courses should avoid asteroid fields."
    Follows NAV.24's open question (it ruled on asteroid fields only).
    Done: course planning treats an asteroid field as a keep-out area
    from its stored `radius_ly`, finds it with NAV.25's obstacle query and
    goes around it with NAV.26's planner, at every scale. Nebulae and
    supernova remnants stay undecided: the course passes through them
    with a note. Parked with the rest of the NAV chain (Boss, 2026-10-09).
    Research (2026-10-09, course-avoidance.md): keep as written.
    `keepout.py` currently contradicts it (no keep-out for fields); the
    whole stored sphere is the radius, no margin (the nebula document
    suggested 1.1 times; the planner already offsets waypoints by 1e-6;
    open question for Boss, default exactly the stored `radius_ly`). Say
    plainly that it is a game rule: a real field would need millions of
    times the belt's density to threaten a ship. Collision risk stays
    keep-out only, with expected hits (about 1e-11 per AU in the main
    belt for bodies over 1 km) as an optional note on the course.
    Prerequisites: NAV.25, NAV.26.

- [ ] **NAV.8 Pages and anchors for stars, planets, moons and belts**
  A star, planet, moon, belt or comet reference opens something:
  by default the system page scrolled to and highlighting that body
  (`/system/<id>#planet-<id>`), with its own System Map scene
  selected, rather than a new page per body. Search results, the
  locate box and bookmarks link this way. Decided (Boss, 2026-10-10 02:48Z, via Foundations lane 1): only stars get pages; planets and moons are anchors on the system page.

- [ ] **NAV.9 Search and locate return references for every kind**
  `/api/search` and `/galaxy/locate` (`queryDb.galaxy_locate`) return
  each hit's reference and parent chain, so any picker can jump to a
  star, planet or moon by name.
  Research (2026-10-09, course-routing.md): build `ref` and the parent
  chain from each panel's own JOIN, not `resolve_object` per row.

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
  stop. Decided (Boss, 2026-10-10 02:48Z, via Foundations lane 1): the stay per stop defaults to 0 minutes but is a user-changeable parameter; stops are the systems on the route, not every planet inside each system.
  Research (2026-10-09, course-routing.md): stops are the systems on the
  route; default stay 0; an optional stay per stop; do not model
  visiting every planet (the research answers the open question this
  way).
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
  Research (2026-10-09, course-routing.md): add per-stop sector id and
  sector-local position to the route data; a stone stop in an unfilled
  sector uses the Galactic frame. For an adjusted course (NAV.28) the
  readout bearing and mark are the first leg's, not the direct line's.
  Built (2026-10-09, course-routing.md): NAV.12 (PR #838) built the
  unbounded route but left three things out, all already in the design
  text of course-routing.md section 2: per-stop sector_id and local
  position (this item needs them), the packed filled-set cache, and the
  adjacent-cell shortcut.
  Design: [docs/design/navigation-frames.md](design/navigation-frames.md)

- [ ] **NAV.36 Unknown-space jumps drawn red and glowing**
  Boss (2026-10-02 01:53Z): "a jump through unknown space is marked in
  red and glows to draw attention to it." Done: a hop NAV.12 flags as
  an unknown-space jump is drawn red with a glow in UX.35's horizontal
  route strip, on the NAV page's own map (`lib/navmap.py`) and on the
  Galaxy Map course (with NAV.20), with a legend entry; with
  `prefers-reduced-motion` it stays red without the pulse; it reads in
  both themes, and the route list also labels it in text so it isn't
  shown by color alone. Prerequisite: UX.35.
  Design: [docs/design/course-routing.md](design/course-routing.md)

- [ ] **NAV.39 Saved courses remember their unknown-space jumps and check them again**
  Done: a saved course (NAV.17) keeps which hops were unknown-space
  jumps when it was saved, and opening it checks again, since sectors
  may have been filled in the meantime; a hop that is now known shows
  as ordinary, and the course says what changed. Prerequisite: NAV.17.
  Design: [docs/design/course-routing.md](design/course-routing.md)

- [ ] **NAV.45 "What's within N pc" from the Galaxy Map and Sector Map**
  Boss (2026-10-02 04:39Z, with NAV.43): "Add TODO items for I want the
  ability to select a location and ask the sytem what stuff (of any
  kind) is within some distance in parsecs." Done: on the Galaxy Map and
  the Sector Map, a selected sector, system, phenomenon or picked point
  gets a "What's within N pc" action on the shared control panel that
  opens the "What's nearby" page (NAV.44, done) with that place filled in. Optional (default left out
  unless Boss asks): the map draws the search sphere and highlights the
  objects inside it. It uses MAP.65's shared control panel and NAV.15's picking of
  any object on the maps. GitHub issue [#677](https://github.com/dwhagar/planetGen/issues/677) (Boss, 2026-10-08 22:24Z): "Add search options so that when I select a star from anywhere I can take it directly to the search tab to search for items / objects within x distance of that point. Should work for stars, phenomena, even planets and moons." So the action is offered for any pickable object, down to planets and moons (NAV.15), and the result lists every kind.
  Research (2026-10-09, course-routing.md): page 1 comes from shells and
  totals are capped ("300+"); show the "not generated yet" counts.
  NAV.43 (built) filters by exact distance after asking for `radius +
  one edge`, with a keyset cursor `(distance, kind, id)`; it should also
  count the pre-placed bright stars, black holes and neutron stars of
  unfilled cells. Open question for Boss (default yes): list those as
  "not generated yet" rows.

- [ ] **NAV.47 Unknown-space jumps stop at scattered stars, black holes, neutron stars and quasars**
  Boss (2026-10-03 05:38Z): "navigational jump path (system to system)
  should use scattered stars in unfilled sectors as potential stops
  along the path.  Same for black holes, neutron stars, quasars, but the
  by far preference is system-to-system, only if a jump crosses unknown
  space should it start looking for other things." Done as he says; such
  stops are marked as their kind in the route.
  GEN.100 now also scatters hypervelocity stars and supernova remnants;
  the list of stops should include hypervelocity stars.
  Research (2026-10-09, course-routing.md): refine only flagged hops;
  stones come from unfilled cells only (filled-cell black holes are not
  stops); reach starts at twice the local stone spacing and doubles;
  define "position at the time" once (`p0 + v (t - t_plan)`) and add
  `galaxy_shape.phenomenon_scatter_at`; add a sub-item for the `kind`
  lookup (a small `phenomenon_scatter_special` side table, preferred
  over a `kind` index on about 5e8 rows). Open question for Boss
  (default: no): do planetary nebulae and supernova remnants count as
  stops? Needs the galactic-motion bug below fixed first.
  Design: [docs/design/course-routing.md](design/course-routing.md)

- [ ] **NAV.48 Offer to generate the uncharted sectors that block a course**
  Boss (2026-10-07 11:47Z): "From navigation menu if there are
  unexplored / uncharted / unfilled sectors in the way that cannot be
  bypassed, it should have the option to generate all sectors between
  the two points." Done: when a route has to cross uncharted sectors,
  the NAV page offers (to admins) a job that charts the sectors along
  the line, with the estimate first, and re-plots after.
  Research (2026-10-09, course-routing.md): add the bypass test (flagged
  edges removed), "chart the cells of the unknown hops only" as the
  default, a one-cell-border checkbox, outside-the-galaxy cells
  excluded, the estimate from `generation.stats.estimate` shown first
  (not a fixed figure; PERF.3's default is 0.2 s per system while the
  research costed charting at about 1 s per system), the confirmation at
  the existing 5,000-sector threshold, and a call to the bright-star
  backfill (GEN.30) per block. Open question for Boss (default: the
  unknown hops only, then re-plot): or every cell on the straight line?
  And how long may a charting job be before the page refuses (default:
  the 5,000-sector confirmation plus the PERF.3 disk refusal)?
  Design: [docs/design/course-routing.md](design/course-routing.md)

- [ ] **NAV.49 Waypoints: pick objects in Star select mode and plot a course through them, kept on the map until cleared**
  Boss (2026-10-07 11:47Z): "Then you can go find another location
  (star, planet, whatever) and select it and you can add a 2nd waypoint,
  waypoints are highlighted at all map levels with an icon it will plot
  a course." And: "Plotted courses should appear on the galactic map and
  stay until cleared." Done: waypoints set from Star select mode are
  marked with an icon at every map level, two or more plot a course, and
  the course stays drawn on every map until cleared.
  Prerequisites: MAP.122, NAV.17.
  Design: [docs/design/course-routing.md](design/course-routing.md)

- [ ] **NAV.52 Port `join_islands` and the k-d tree to cKDTree**
  `galaxy/nav_graph.py`'s pure-Python k-d tree is 44 times slower than
  cKDTree at 10,000 points, and `join_islands` takes 15 to 18 s for
  10,000 points in 4 islands and 6.9 s for 53,000 points in 394
  (vectorised: 0.07 s and 0.23 s, identical edges). Done: cKDTree
  versions with the same edge set, compared on the 106,529-point set. A
  performance cliff, not a correctness bug; NAV.10 shipped on the old
  code. Avoid `np.unique` on very large arrays in new code (9.4 s on 6e6
  int64 against 0.11 s for sort-plus-diff).
  Prerequisites: none.
  Design: [docs/design/course-routing.md](design/course-routing.md)

- [ ] **NAV.54 Keep-out radii for asteroid fields, supermassive holes, moons and nebulae (NAV.24 built)**
  NAV.24 shipped with these gaps: (a) asteroid fields get a hard
  keep-out of `radius_ly` in km, basis "radius" (NAV.51 asks for this;
  `galaxy/keepout.py` and `db/query.py` `keep_out_radius` return none
  today); (b) a supermassive hole or quasar keeps out a sphere of
  influence `G M / sigma^2` (sigma = 100 km/s, 1.7 pc at 4e6 solar
  masses), a quasar the larger of that and the 10 kW/m^2 radiation
  radius, not the event horizon (0.08 AU, which a course would graze);
  (c) a moon keeps out its Hill radius against the planet's mass; (d)
  the pass-through note for nebulae and remnants ("No mass is stored for
  it...") gives the wrong reason: the drawn shape fills 13% of its
  bounding sphere and the gas is thin. A test for each. Open questions
  for Boss (defaults taken): nebulae and remnants pass through with a
  note (chord length, column density, extinction) and no hard keep-out,
  with an optional soft cost later; a star inside its own system keeps
  out its radiation radius `max(10 R_star, sqrt(L / (4 pi F_lim)))`,
  F_lim = 50 kW/m^2 (0.17 AU for the Sun, 165 AU for an O3 star), hard
  and exempt when an end lies inside it; belts inside a system pass
  through with a note; the galactic centre keep-out of about 1.7 pc is
  exempt for any course to a star inside it.
  Prerequisite: NAV.51.
  Design: [docs/design/course-avoidance.md](design/course-avoidance.md)

- [ ] **NAV.55 A tuning block for the keep-out knobs**
  Put the keep-out knobs in `tuning.py`: F_lim = 50 kW/m^2, the star
  floor of 10 R_star, the below-resolution fraction (0.002 of the leg)
  and the polish tolerances.
  Prerequisites: none.
  Design: [docs/design/course-avoidance.md](design/course-avoidance.md)

- [ ] **NAV.56 Census of overlapping keep-out spheres in a generated galaxy**
  Count how many pairs of keep-out spheres overlap in a generated
  galaxy. The answer decides whether NAV.26's cluster fallback ever
  fires.
  Prerequisite: NAV.51.
  Design: [docs/design/course-avoidance.md](design/course-avoidance.md)

- [ ] **NAV.57 The Intergalactic Frame in navigation-frames.md and `navigation.py`**
  A frame type and the rule for choosing it (stage 2 of GEN.9).
  Prerequisite: GEN.157.
  Design: [docs/design/multiple-galaxies.md](design/multiple-galaxies.md)

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
  Research (2026-10-09, multiple-galaxies.md): point its "Design:" line
  at docs/design/multiple-galaxies.md. Answer the three questions as in
  section 8: the real Local Group; no `galaxy_id`; a second database
  chosen by `?db=`, not at login. Stage 1 (`neighbor_galaxies` table,
  the data file, the `planetgen plan` step) is separate from stage 2
  (`linked_database`, the Intergalactic Frame, `/galaxies`), which costs
  a galaxy switch in the UI and a `galaxy_links` control table. Open
  questions for Boss (defaults taken): real neighbours (about 25 rows, a
  hand-built data file) with a per-plan option for generated ones; the
  proper embedding (see the handedness item); database per galaxy with
  the host's frame as the intergalactic origin; no second fully
  generated galaxy soon; the data file is edited by file and a re-plan,
  not from a page.
  Research (2026-10-09, globular-clusters.md): the globular-cluster
  system for generated galaxies is GEN.164.

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
  Research (2026-10-09, fill-order-curves-and-core.md): correct "bulge
  scale radius, 200 pc" (the shipped default `--bulge-scale-radius-pc`
  is 1,580 pc along the bar, 620 and 430 across and up; 200 pc was the
  revision-2 sphere). State the 200 pc cost honestly: 7,857 sectors,
  7.45 million systems, 17 days on one worker, about 403 GB, refused by
  PERF.3 on a disk under about 1.6 TB. Change the suggested default to
  50 pc (rings 0..12, 531 sectors, 0.51 million systems, 28 GB), scaled
  as min(50 pc, 0.25 x bulge scale radius) for other shapes, and keep
  200 pc as a typed value. Name the ring serpentine as the fill order.
  The model has no central cusp (the core sector holds about 1,000
  systems, 5 to 15 times below the real nuclear disc, 1,000 times or
  more below the Sgr A* sector). Open question for Boss (default 50 pc;
  leave the model alone).

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
  Research (2026-10-09, atmospheres-retention-and-classes.md): R becomes
  1.6 to 3.2 Earth radii and 3 to 10 Earth masses (cap 4 Re); S and R
  share mass 2 to 10 and split by bulk density relative to R = M^0.27
  (rocky if R <= 1.12 M^0.27); add the Hycean window (M 1 to 10, R 1.4
  to 2.6, T_eq 150 to 500 K, density 1.5 to 3); X goes to 500 to 11,500
  km; W is moon-only (not a rogue); Y needs shoreline ratio under 0.7; Z
  has M-dwarf and G/K sub-types (design doc 5.1). Open question for Boss
  (default R): rocky rogues of 10 to 16 Earth masses become class R
  rather than a widened S (cap S at 13,500 km); a 1-bar "dune world" is
  a 10 percent variant of N.

  - [ ] **GEN.33 One class per PR, each with its tests**
    Suggested split: R and S (the commonest missing types) first; then
    U, W, X, Y and Z. Each adds the class to `PLANET_CLASSES`, the
    reference pages, the rogue flag, moon eligibility, and a
    distribution test over 1,000 generated systems.
    Prerequisite: GEN.89.
    Plan (2026-10-07): Moved to phase 2 under GEN.90 (the class
    refactor), after the habitability score GEN.89.
    Research (2026-10-09, atmospheres-retention-and-classes.md): confirm
    the PR order: retention module first, then R, the S fix, U, W/X, Y,
    N, Z, the Hycean flag, the GEN.29 sweep.

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
  Research (2026-10-09, atmospheres-retention-and-classes.md): known
  oddities to fold in: classes A and B hold 13 to 17 kPa at S = 4.4 and
  v_esc 3.7 km/s, about 50 times beyond the cosmic shoreline; N sits in
  the ecosphere at S about 1.05 (below runaway); T and I overlap; rogue
  S reaches 17,600 km (recommend class R at 10 to 16 Me, or cap S at
  13,500 km); K, L and P fail the Armstrong and 20 kPa limits (P needs
  the stage cap from GEN.92).

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
  Research (2026-10-09, sampling-backfill-and-resume.md): fix the stale
  code pointers: the backfill is `generation/run_galaxy.py`
  `backfill_bright_stars_around` and `generation/bright_stars.py`
  `backfill_cells` (not `generate.py` and `brightStars`), drawn cell by
  cell in chunks of 200 sectors since GEN.44.

  - [ ] **GEN.42 A pass that drops sectors from a region by probability**
    One pass over a region (a sector block, or a tier's shell) uses the
    star density calculation (`galaxyDensity.population_densities`) to
    eliminate some of its sectors, chosen so that the region's chance of
    getting each bright star is still what the density says; the
    backfill then draws only in the sectors left. Done: the pass runs
    before the backfill, the backfill visits fewer sectors, and the same
    seed still gives the same stars wherever it lands (see GEN.39).
    Research (2026-10-09, sampling-backfill-and-resume.md): re-scope
    from "a pass that drops sectors" to an exact block-first draw for
    large volumes: a Poisson candidate count from the ring-and-layer
    density bound, then thinning, with the per-cell path kept for dense
    blocks (chosen per block: candidate path when expected candidates
    are under about 5% of the block's cells). The literal drop-and-boost
    form is rejected (dispersion 1.136 against 1.0, multiplicity
    chi-square 743,000 on 14 degrees of freedom). Done-test: the
    stratified multiplicity chi-square (half-decade strata of lambda,
    3,000+ replicates) run on the real `backfill_cells`, an assertion
    that every cell's rate is at most its block bound, and the order,
    partition and sub-region identity test. GEN.43's floor test becomes
    the lowest-lambda strata of this test (no cell's probability is
    altered, and `MIN_RELATIVE_DENSITY` keeps every cell above zero), so
    GEN.43 is folded in here. Answer to the GEN.41 question: no-go for
    the literal form; revisit when a run touches 10^5 or more cells. Ask
    PERF.31 to measure the SQL part at scale (2 ms per sector warm,
    about 0.8 ms of it SQL, plus a `sector_stats` row per visited
    sector). Include one block-first backfill in TEST.77's golden
    galaxy; the bound uses the `detmath` helpers.

- [ ] **GEN.55 Same seed, same data: a sector's contents depend only on the seed, the version and its address (internal)**
  Boss (2026-10-02 01:40Z): "Ok use a 128 bit value and store the seed in
  the database, and put it in the log at the top of any generation, also
  populate the TODO upward from here to eventually build a system that a
  version number and a seed value would reproduce the same galaxy by the
  end of the phases."
  Changed by Boss (2026-10-09 20:42Z): "remove the idea of us letting
  the USER regenerate an entire galaxy from the seed and version, we'll
  use that internally, but no need for it to go anywhere else." The seed
  stays as an internal mechanism: parallel workers, fills, backfills, the
  phenomenon scatter and its sector fill, and the settle step all rely on
  the same seed and address giving the same data. There is no user-facing
  rebuild of a galaxy from a seed and a version: no `generate.py
  reproduce`, no seed and version pages, no JSON net-difference file, no
  daily merge or backup slots. This item is done when its last internal
  sub-item is. Order:
  - Phase 0: GEN.39 (per-unit seeds, after PERF.21) with DB.6 and
    OPS.10 (done).
  - Phase 1: GEN.56 (every draw seeded), GEN.57 (a sector's contents
    depend only on the seed, the version and its address), DB.7 (the
    version kept with each sector), GEN.58 (a fingerprint), TEST.77 (the
    golden-seed test), OPS.13 (each update records the key, keeping the
    last 10) and OPS.14 (a warning when the running key differs from the
    galaxy's). GEN.47's nebula field (done, PR #419) uses the derived
    seeds.
  - Phase 2: PERF.18 and GEN.42 give the one-process stars for one seed
    (already in their Done text); OPS.15 has each update say whether it
    changes generated output.
  - Phase 3: API.17 makes a remote run produce what the server would
    (kept for remote generation; open question for Boss).
  Anything that draws new randomness later (GEN.42, PERF.18, API.12,
  API.13) uses the derived seeds and keeps TEST.77 green.
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
    scatter but not star for star. Boss decided that statistically
    equal is good enough; the layer scatter need not be star-identical.
    Design: [docs/design/reproducible-galaxies.md](design/reproducible-galaxies.md)
    Plan (2026-10-07): Name collisions disappear with GEN.67 (names come
    from IDs), so the name part of this item is dropped, and GEN.63 with
    it.
    Note (2026-10-09): order: it must land before TEST.77 (open order
    dependencies: name collisions by save order, the ID-cell counter,
    population seeds keyed on database ids, nearest links). Research
    Lane 1 also notes the admin regenerate paths for planets, moons and
    belts (`admin/edits.py`) draw from the ambient process stream, so an
    admin-edited sector is not rebuildable from the seed alone; DB.17's
    repair therefore replays the edit log rather than the seed.

  - [ ] **OPS.15 Each update says whether it changes generated output**
    The second half of Boss's update "sweep". Done: after OPS.13's row
    is written, the update script fingerprints a small fixed region
    (GEN.58) from the galaxy seed under the new code, compares it with
    the fingerprint stored in the previous history row, stores the new
    one, and reports "generated output unchanged" or names the sectors
    that differ. The live galaxy is not touched. A test runs it across a
    change that alters a sector.
    Research (2026-10-09, reproducible-galaxies.md): the compared
    fingerprint is the battery digest from the epoch item
    (`generator_battery_digest()`), compared with the epoch and the
    previous history row.
    Design: [docs/design/reproducible-galaxies.md](design/reproducible-galaxies.md)

  - [ ] **OPS.16 A daily maintenance script for Linux, macOS and Windows**
    Boss (2026-10-02 02:28Z): "We're going to have to build a
    maintenance script for powershell and bash that will run the
    positional update script, then kick off this delta script ..."
    Scope cut by Boss (2026-10-09 20:42Z): the delta merge and the
    18-slot JSON backups (GEN.61, OPS.18) are dropped with the user-facing
    galaxy rebuild, so this script runs the positional update only (and
    any later daily step). Decided (Boss, 2026-10-09 20:52Z): keep the
    daily positional update.
    Done: `scripts/maintenance.sh` (Linux and macOS) and
    `scripts/maintenance.ps1` (Windows) run once a day: the
    positional update (`updateOrbits.py`). A lock keeps two runs
    from overlapping, every step logs to the normal log, and the script
    exits non-zero on any failure. Optionally it also runs OPS.15's
    fingerprint check, so changed output is noticed daily. A test runs
    it on a small galaxy and checks a second run started during the
    first exits at once. Prerequisites: none.
    Research (2026-10-09, db-check-and-parity-repair.md): the parity
    export and checksum must not cover the columns the orbit update
    rewrites, or the daily run invalidates all of them; the run also
    refiles systems and phenomena into other sectors, so sector
    membership of placed objects changes daily.
    Research (2026-10-09, ops-scheduling-and-rotation.md): restructure:
    implement the logic once in `planetgen.cli.maintenance` (steps,
    lock, `--if-due`, `--dry-run`, `--only STEP`, rotating log file,
    exit codes 0/1/2/75) and make `scripts/maintenance.sh` and `.ps1`
    thin launchers that find the right Python and run it as the web
    user. Take `locks/maintenance.lock` plus a check of
    `jobs.active_job()`; exit 75 when skipped; alert after N skipped
    days. Treat Redis as optional. Iterate every galaxy database. Run
    orbits in-process or as `planetgen.cli.orbits` (this item's text
    says `updateOrbits.py`, which no longer exists). Open question for
    Boss (default: skip the day while a Generate job runs, exit 75,
    alert after 3 skipped days).
    Note (2026-10-09): The Windows parts are dropped by OPS.39 (Boss,
    2026-10-10 04:25Z).
    Design: [docs/design/reproducible-galaxies.md](design/reproducible-galaxies.md)

  - [ ] **OPS.17 Install and update set up the daily maintenance schedule**
    Done: `install.sh`, `install.ps1`, `update.sh` and `update.ps1`
    (with `deploy-common.*`) set up OPS.16's daily run: a systemd timer
    or cron entry on Linux, launchd on macOS, Task Scheduler on
    Windows. Update leaves an existing schedule as it is and adds a
    missing one. The deployment docs say how to change the time or turn
    it off. Same scripts as OPS.7, OPS.8 and OPS.13, so it lands after
    them. Prerequisite: OPS.16.
    Research (2026-10-09, ops-scheduling-and-rotation.md): re-scope to
    extending `examples/maintenance/install-maintenance-timer.sh` and
    `install-maintenance-task.ps1` (they exist) and calling them from
    install and update with the rules create if missing, leave if
    present, refresh the service text by version stamp. Default per OS:
    systemd timer (UTC `OnCalendar`, `Persistent=true`), cron.d
    fallback, launchd with `RunAtLoad` and `--if-due`, Task Scheduler as
    SYSTEM with `StartWhenAvailable`. Retire the monthly
    `planetgen-orbits@` timers, tasks and plists when the daily run is
    installed. Keep the monthly `planetgen-update.timer` (it takes the
    maintenance lock and waits up to 30 minutes; open question for
    Boss). Add `maintenance.schedule: auto|off` to `config.json`. Fix
    and test `org.planetgen.update.plist`, and correct
    `docs/deployment/macos.md` (an asleep Mac runs at wake; only a
    powered-off Mac skips).
    Note (2026-10-09): The Windows parts are dropped by OPS.39 (Boss,
    2026-10-10 04:25Z).
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
    every CI Python leg. Prerequisite: GEN.57. [generation, infra]
    Research (2026-10-09, reproducible-galaxies.md): replace "pinned per
    release with the version it was made on" with "pinned per
    `generator_epoch`, one golden for all legs". Add tier A (the
    database-free battery, 3 to 5 s) to the Windows and macOS CI jobs
    (`ci.yml`'s matrix is Linux only; `macos-latest` is Apple Silicon);
    keep 1 against 4 workers as tier B and the hash-seed and locale
    probe as tier C; add one full-text "Rosetta" sector for readable
    diffs. The epoch item owns the `golden_update` tool and the
    `--check-pr` rule. Needs the deterministic math helpers first. Once
    GEN.42 lands, include one block-first backfill in the golden galaxy.
    Design: [docs/design/reproducible-galaxies.md](design/reproducible-galaxies.md)
    Prerequisites: GEN.135, OPS.28.

  - [ ] **API.17 Remote generation reproduces what the server would make**
    Done: a remote run (API.12's download, API.13's generation without a
    database) with the same seed and release produces exactly what the
    server would for those sectors, checked by fingerprint (GEN.58); and
    API.8 can verify an upload by re-running a sample of its sectors on
    the server and comparing. Prerequisites: API.12, API.13, GEN.57.
    Decided (Boss, 2026-10-09 20:52Z): keep the fingerprint check. This
    is remote generation matching the server, not a user rebuild of a
    galaxy from a seed and version.
    Research (2026-10-09, reproducible-galaxies.md): the remote path
    calls the same pure generation function; compare fingerprints at 9
    digits; per api-design-standards.md the handshake compares epoch and
    battery digest rather than the full version key. A remote machine
    with a different OS, architecture or Python micro version cannot
    match the fingerprint unless the key is loosened or the check is
    advisory.
    Design: [docs/design/reproducible-galaxies.md](design/reproducible-galaxies.md)

- [ ] **GEN.83 A planetary habitability index (PHI)**
  Boss (2026-10-03 05:38Z): "Create a habitability index based on the
  pressure, temperature, composition, etc...  This will use several
  documents in the design docs folder: Atmospheric Toxicity.md, Chemical
  Habitability.md, Naturally Occurring Ionizing Radiation.md, Planetary
  Habitability and Speculative Xenobiology.md, Planetary Habitability
  Index.md, Speculative Xenobiology Extremes.md, Speculative Xenobiology
  Examples.md, Mathematical and Algorithmic Implementation of the
  Planetary Habitability Index.md" Done when the subitems are.
  Prerequisites: GEN.87, GEN.89.
  Design: [docs/design/habitability-index.md](design/habitability-index.md)

  - [ ] **GEN.87 Surface radiation dose**
    Done: each planet stores its surface dose from column mass (P0/g),
    magnetic field, cosmic rays (more inside a compressed heliosphere,
    using the nebula containment of GEN.75) and flares, plus an
    ozone-loss flag, following the radiation doc's tables.
    Research (2026-10-09, activity-magnetism-radiation-hydrosphere.md):
    add UV/ozone and the galactic-hazard flag to the done criteria; the
    SEP term is rate-dependent and uses R0 = 0.2 GV with a polar-cap
    term; limit the heliosphere multiplier to 2.5x and fade it by X =
    300 g/cm2; correct the dose table in the radiation document (cosmic
    0.39 not 2.4 mSv/yr; Mars 0.64 mSv/day dose equivalent). Open
    question for Boss (default on): the supernova and GRB flag lowers
    the score only if the flagged event rate exceeds one lethal event
    per 100 Myr.
    Design: [docs/design/habitability-index.md](design/habitability-index.md)

  - [ ] **GEN.89 The habitability score for every planet and moon**
    Done: every planet and moon gets the scores GEN.84 defines, with a
    colour tier and a human equipment profile (shirtsleeve, mask, mask
    and scrubber, sealed suit, full life support), shown on its page
    and searchable.
    Research (2026-10-09, atmospheres-retention-and-classes.md): add an
    `equipment_tier` column (0 to 4) and the four domain tiers (the five
    names match this text). Stored per star: `log_lx_lbol`,
    `lxuv_erg_s`, `xuv_saturated`, `flare_n33_per_yr`, `flare_alpha`;
    per planet: `xuv_flux_erg_cm2_s`, `xuv_exposure_index`,
    `magnetic_moment_a_m2`, `magnetopause_rp`, `surface_dose_msv_yr`
    (plus the GCR, SEP, crust and UV parts), `ozone_loss_flag`,
    `ocean_class`, `ocean_depth_km`, `hp_ice_km`, `ice_shell_km`,
    `land_fraction`. It returns the lowest tier, with a reason string,
    for any pulsar, black hole or X-ray-binary planet. Open question for
    Boss (default: store the two flare numbers, derive the rest).
    Prerequisite: GEN.87.
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
  Research (2026-10-09, atmospheres-retention-and-classes.md): the
  "Known conflicts" list gains: A and B atmospheres, S in non-hot zones,
  N zone placement, N and Q life data hidden but stored.

  - [ ] **GEN.91 Classes like S and V in the hot and cold zones**
    Boss (2026-10-03 05:38Z): "We need a planet class for like Class S
    and V but in the other zones, so that will require research on when
    an atmosphere is possible and when it isn't." Candidates from the
    research: a desiccated post-runaway super-Earth (hot), a
    hydrogen-rich Hycean super-Earth (cold or outer). Done: atmosphere
    retention rules (escape velocity, XUV, temperature) written down,
    then the classes added one per PR under GEN.33's rule.
    Research (2026-10-09, atmospheres-retention-and-classes.md): (1)
    `planetgen/physics/atmosphere_retention.py` with the shoreline,
    Jeans-species, condensation and energy-limited functions and a
    solar-system plus exoplanet test table (design doc section 3.6); the
    XUV history and exposure index come from GEN.86's activity module,
    not the stand-in. (2) Classes follow in the PR order of section 5.2.
    Record the answer to "when is an atmosphere possible": shoreline
    ratio below 30, species by lambda (30 kept), the cold side by vapour
    pressure, H2 and water by energy-limited XUV. "S and V in the other
    zones" is mostly S itself with air plus a re-scoped N (hot) and an
    extended X or an R Hycean flag (cold). Letters: none are free after
    R, U, W, X, Y, Z; Q must not be reused until the database is checked
    for Q rows. Open question for Boss (default Option A): re-scope N,
    extend X, Hycean as an R flag and retire T, or two-character class
    codes.
    Prerequisite: GEN.33.
    Design: [docs/design/habitability-index.md](design/habitability-index.md)

  - [ ] **GEN.92 Life and its highest stage follow the habitability score**
    Done: whether a world has life, and how far it gets, comes from its
    scores instead of its class alone; life chemistry follows the
    atmosphere, radiation and solvent rather than the star's letter
    only; existing galaxies are unchanged until regenerated.
    Research (2026-10-09, atmospheres-retention-and-classes.md): life
    stage is the latest stage reached on the pace and gated by the score
    (rule table in design doc 7.2); add `t_hab` (needs the
    pre-main-sequence rule and optionally main-sequence brightening,
    which `stellar_evolution.py` lacks); flag the `fast` and `slow`
    scales as fiction (open question for Boss, default: keep them and
    cap `fast` at stage 3); add planetary tidal locking (needs GEN.104).
    Open question for Boss (default yes): dose above 0.1 Sv/yr forces
    equipment tier 3 and above 1 Sv/yr tier 4.
    Prerequisites: GEN.89, GEN.28.
    Design: [docs/design/habitability-index.md](design/habitability-index.md)

- [ ] **GEN.93 Nebula conditions in planet generation**
  Boss (2026-10-03 05:38Z): "Account for nebula temperature and
  conditions when generating planets and calculating surface conditions.
  You'll have to do a fesability study for our nebula classes in if a
  planet could even form there and what would happen if one did, how it
  would change the planets development and envirionment according to the
  science." Done when the subitems are.
  Prerequisite: GEN.95.
  Design: [docs/design/nebula-and-asteroid-field-classes.md](design/nebula-and-asteroid-field-classes.md)

  - [ ] **GEN.95 Nebula conditions applied when planets and surfaces are generated**
    Done: generation knows a system's surrounding cloud while it runs
    (today `surrounding_cloud` is set only when a system is loaded) and
    applies GEN.94's rules to planet formation and surface conditions.
    Research (2026-10-09, nebula-and-asteroid-field-classes.md): use the
    plan in "Applying the rules at generation": `StarSystem.apply_cloud`
    after the position is chosen; extend `surrounding_cloud` with
    extinction, radius, age and hosts (it is set only at load today);
    add `age_myr` and `q_total` to `nebulae` and `age_years` to
    remnants. Depends on GEN.94 and the habitability index (radiation
    tier).
    Research (2026-10-09,
    exotic-environments-planets-and-compact-binaries.md): Add to its
    Done: (1) the photoevaporation cut radius scales with host mass (cut
    x M/Msun; 200/50/10 AU is the solar-mass calibration); (2) planetary
    nebulae H to L allow a flagged circumbinary young disc only when the
    central star is a close binary, never planets; (3) an engulfed giant
    of over 5 Jupiter masses sets a "swallowed giant" flag (a flag only,
    for the red-nova anomaly); (4) surviving planets near the engulfment
    radius get an eccentricity draw (up to about 0.3). Open question for
    Boss (default yes): the swallowed-giant flag.
    Decided (2026-10-09): Boss (2026-10-09 18:37Z): the defaults are
    accepted, including the swallowed-giant flag. The orbit-expansion
    law for adiabatic mass loss is verified: a (M_star + M_planet) =
    constant with the eccentricity unchanged (Jeans 1924; Hadjidemetriou
    1963; Veras et al. 2011), valid when the mass-loss time is much
    longer than the orbital period and the loss is isotropic; impulsive
    loss uses vis-viva and unbinds an orbit when more than half the mass
    goes; engulfment overrides. Worked example: 1 AU round a 1 Msun star
    that drops to 0.6 Msun moves to 1/0.6 = 1.67 AU.
    Prerequisite: GEN.89.
    Design: [docs/design/nebula-and-asteroid-field-classes.md](design/nebula-and-asteroid-field-classes.md)

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
  Research (2026-10-09, nebula-and-asteroid-field-classes.md): the draw
  must use a young, age-capped population (measured 76% K/M giants
  otherwise); Hill-sphere separation needs a spread-out placement
  (default; classes D and E) or a flat 0.5 ly override for class C if a
  second star is wanted; "regenerate in place" must also rebuild the
  system's planets (use the same `apply_cloud` as GEN.95). Specify
  mark-only skewing (keep each cell's count, tilt the type law; a weight
  of about 40 to 50 gives about 80% O/B from the natural 8%) and
  conditional redraw on luminosity keyed by (star seed, nebula id);
  leave the star alone when no compatible type exists. Open questions
  for Boss (defaults taken): backfill down to 750 Lsun only inside the
  nebula; only 750+ stars plus normal field density, not the realistic
  crowd of low-mass stars.
  Design: [docs/design/nebula-and-asteroid-field-classes.md](design/nebula-and-asteroid-field-classes.md)

- [ ] **GEN.101 Fill order: nearest sectors first along a pruned Hilbert octree curve**
  Boss (2026-10-07 17:11Z) approved keeping the Hilbert order, with a
  logged jump where the ball cuts the curve.
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
  Research (2026-10-09, fill-order-curves-and-core.md): re-scope the
  "Done" wording to what can be true. Boss's 2026-10-07 17:11Z decision
  (keep the Hilbert order with a logged jump) stands, but "every step
  face-adjacent, visit once, zero backtracking" cannot hold for a pruned
  ball on any curve (only 16 of 1,335 balls with R squared up to 1,600
  pass the parity test; radius 1 forces at least 4 jumps). Write: the
  order is deterministic and face-adjacent for about 99% of steps, and
  every other step is logged with its length (count and mean at the end
  of the run); the unit-step test asserts the adjacent share (at least
  95% in the 50 pc sphere and 200 pc disc cases) and determinism, not
  100%. Choice of order: ring serpentine for the core and axis-centred
  discs and cylinders (0 jumps), a greedy nearest-first walk for any
  other centre, Hilbert only as a sort key or tie-break inside 4-voxel
  blocks (plain pruned Hilbert starts at the far edge of the cube and
  fails nearest-first; the sector grid is not a voxel lattice). The
  depth formula `D = ceil(log2(2R))` is one level short when R is a
  power of two. Fix `_neighborhood_batch` in `run_galaxy.py` (its
  docstring says nearest first as enumerated; the enumeration is
  ring-major, and off the axis the first 10% of the list holds none of
  the nearest 10%): sort before applying a limit. Open question for Boss
  (default as written): keep Hilbert as approved, using the
  4-voxel-block variant so the first sector is the centre, or take the
  serpentine and greedy walk?

- [ ] **GEN.102 Investigate filling all near-zero-density void space at once**
  Boss (2026-10-07 11:47Z): "Investigate filling all void space (space
  that plots to have a 0 or nearly 0 star density) at once, since those
  sectors have the least in them, we can afford to fill them." Done: a
  count of such sectors, the time and storage a fill would take, and a
  recommendation for Boss.
  Research (2026-10-09, sampling-backfill-and-resume.md): re-scope the
  "count": there are no zero-density sectors inside the outline (minimum
  density 0.069). Define "void" as expected systems below 1: 12.4% of
  sectors (2.97 billion) and 0.77% of systems; the density-below-0.1
  subset is 2.44% of sectors (5.8e8). Cost (design doc 5.2): about 330
  days on one worker for density below 0.1, about 5 years for expected
  systems below 1; the fixed cost per sector is the lever (the timing
  report's three levers could roughly halve the 40 ms). Recommendation:
  do not pre-fill; keep lazy fill. Open question for Boss (default no):
  do you still want the sparse rim filled?

- [ ] **GEN.103 Research where each star type and phenomenon belongs in the galaxy's structure**
  Boss (2026-10-07 11:47Z): "Research and add tiered structure
  information about where stars of different types would be in the
  galactic structure, what kinds of phenomena have similar restrictions,
  etc." Done: a design note (thin and thick disk, bulge, bar, halo,
  arms; ages and metallicity; where O and B stars, white dwarfs,
  remnants, clusters and nebulae sit) and the rules added to the density
  and type draws.
  Research (2026-10-09, anomalies.md): adopt the population table in
  "Placement by population" (young disk, thin disk, bulge and thick
  disk, halo) and the young-star weight for `phenomenon_scatter`; add
  the population-dependent multiplicity rules of
  multistar-and-compact-systems.md 6.6 (metallicity-dependent
  close-binary fraction, massive multiples in young regions, a wide-pair
  survival cap shrinking in dense regions).

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
  Prerequisites: GEN.115, GEN.109, GEN.110.
  Design: [docs/design/orbital-updates.md](design/orbital-updates.md)
  Research (2026-10-09, orbital-solvers-and-integrators.md): the Full
  Algorithm's one-year macro step and "nightly" language are superseded
  by real time: perturbation integration is per due object over its
  whole interval with an adaptive integrator (`DOP853`; IAS15-class for
  encounters), not Velocity Verlet in day steps. REBOUND is not a
  dependency (GPL-3.0-only against CC0). The text "`updateOrbits.py` is
  today's positional update" is stale: the code is `planetgen orbits`
  (`cli/orbits.py`) with `db/store.py` `advance_*` and
  `physics/kepler.py`. Replace the collision-detection routine with the
  closest-approach form in collisions-and-mergers.md section 5.2: fix
  `new_radius` for giants (`giant_radius_km(M_new)`), state the
  pair-window rule and window cap (5.3), compute differences from
  sector-relative values, say that R_eff with focusing is a capture
  cross-section, and that the 2.44 Roche factor applies to fluid bodies
  only (rocky 1.26). Cross-link collisions-and-mergers.md.

  - [ ] **GEN.115 The galaxy's own gravity: a smooth disk, bulge and halo potential**
    Boss (2026-10-07 17:11Z) approved scaling the potential with the
    galaxy's shape.
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
    Research (2026-10-09, galactic-potential.md): (a) replace
    `calculate_galactic_orbit` and its two constants
    (`GALACTIC_ROTATION_FLAT_VELOCITY_KMS`,
    `GALACTIC_ROTATION_CORE_RADIUS_PC`) with the potential's circular
    speed in the same change, and reseed or migrate the stored galactic
    velocities (stars, systems, phenomena, stand-alone facilities); the
    old curve gives 206 km/s at the Sun and would leave every star
    eccentric in the new potential. (b) Lengths follow the shape, masses
    scale by s^p (default p = 1), with b_d from the shape's scale height
    and c_b from its bulge radius; correct the shipped values to a disk
    scale length of 2,600 pc, bulge 1,580 pc and radius about 15.3 kpc
    (not 2,800 / 200 / 15 kpc). (c) The flattened log halo as written
    (v0 = 175, R_c = 2.5) fails the item's own tests (243.8 km/s at
    8.128 kpc): NFW is the default and the only halo the 229.3 test
    covers; the log halo needs v0 about 177 km/s, R_c about 5.1 kpc and
    its own tests, or is dropped. Plummer softening of 1 pc is for the
    central black hole only. (d) Test tolerances: 229.3 +- 0.1 at 8.128
    kpc; the 213 to 231 band on a 2 kpc grid (continuous extremes 213.47
    and 230.76); Lz conserved to 1e-10; energy drift below 1e-3 over 1
    Gyr at the recommended step; a circular star keeps its radius to 1
    pc over 5 Gyr with dt <= 1 Myr; Kz(1.1 kpc) about 75 Msun/pc2 and
    local density about 0.10 Msun/pc3 at 5%; finite acceleration at (0,
    0, 0). (e) Integration is kick-drift-kick with a per-star step (eta
    0.05 of the pericentre crossing time, dt_max 5 Myr, block steps, a
    step cap with a flag); seeding uses Jeans asymmetric drift and a
    dispersion table from the potential; these could become sub-items. A
    facility has no dispersion and uses the same circular speed
    (integrated in the potential once it has a vertical offset). Open
    questions for Boss (defaults taken): constant 229.3 km/s whatever
    the galaxy's size; Sun distance 8.128 kpc as the test radius and 8.2
    kpc for placing the Sun in the density model (change
    `GALACTIC_CENTER_DISTANCE_LY` from 25,800 ly); NFW only for now.
    Design: [docs/design/orbital-updates.md](design/orbital-updates.md)

  - [ ] **GEN.109 N-body influence from every object inside the largest nearby Hill sphere plus the galactic gradient, with a Hill-radius warning**
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
    Boss (2026-10-09 09:01Z, relayed): drop "nearest 10". The influence
    radius is the maximum Hill sphere of the largest nearby object; the
    point masses of every object inside it, plus the galactic
    gravitational gradient, give each visited object's new vector. This
    replaces the nearest-10 rule above (the Hill-radius warning, the
    update-only trajectory changes and the sector-geometry check stand).
    The influencer search test therefore compares the objects found
    inside that radius against a brute-force search.
    Research (2026-10-09, orbital-solvers-and-integrators.md): replace
    the title and the "nearest 10" wording with Boss's 2026-10-09 rule
    (the item's last paragraph); the "top 10" in "Orbital Update Full
    Algorithm.md" and the old orbital-updates.md text is outdated. Add:
    per-sector `max_mass_kg` and `max_mass_id` columns; a Jacobi
    constant `k(R)` per ring; cover lists for objects with a Hill radius
    above 1.5 sector edges (without them 30% of "target inside a heavy
    object's sphere" cases are missed); a search radius widened by the
    cell reach (about 4 pc; 21.6% of neighbours are missed otherwise); a
    `max_exact_members` cap (default 64, to agree with
    collisions-and-mergers.md; the rest as one monopole per sector);
    `geometry.enumerate_sectors_within_radius` in place of the
    `ResolveNeighborSectors` stencil. Plummer softening (1 pc) on the
    central black hole only; ordinary point masses unsoftened (or 1e-3
    to 1e-2 pc); delete the `a_min` filter (it removes every star at
    day-scale steps). Warnings use the gravitational criterion of
    collisions-and-mergers.md section 3.3 (bound or capture candidate,
    deflection of at least 1 degree, or same system), not a plain Hill
    entry (about 200 a day on a full-size galaxy; plain entries go to
    debug); hierarchy members are Kepler-tree nodes exempt from mutual
    warnings and captures, one point mass per system at its barycentre,
    and a flyby inside the outermost orbit's Hill zone is flagged
    "disrupts hierarchy". Collision candidates are a separate mass-blind
    distance query. Extend the test per
    orbital-solvers-and-integrators.md section 6.4 (grid against brute
    force and `cKDTree`, the exact pairwise Hill criterion, a negative
    control without the margin, a cap test, a seed-reproducibility
    test); a lone rogue has an empty set, a rogue within 1 pc of a star
    sees that star. Open questions for Boss (defaults taken): full Hill
    spheres with cover lists rather than clipping to the sector stencil;
    a star cluster counts as its members; no warnings for unbound flybys
    or hierarchical members.
    Prerequisite: GEN.115.
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
    Research (2026-10-09, collisions-and-mergers.md): keep Boss's
    wording as the intent and re-scope to the outcome table in the
    design doc section 5.1. A terrestrial body hitting a giant, brown
    dwarf or star is absorbed (no field); hit-and-run and partial
    erosion leave both bodies alive; the fusion check (0.075 Msun, 78.6
    Mjup, a setting) applies to any substellar pair, since two giants
    alone cannot reach it; drop "planetary nebula" (a merger product is
    not one): Option A, no nebula and the star flagged "merger remnant"
    (default), or Option B, a short-lived dark ejecta shell, or reserved
    class X "Merger remnant nebula" text. Such collisions are about
    1e-11 per object per 10 Gyr, so the handler matters for dense
    nebulae and hand-made scenarios. Event-made fields need `mass_kg`,
    `formed_at`, `origin_event_id` and a velocity (add a nullable
    `mass_kg` to `asteroid_fields`); the radius grows from about 1e4 km
    at 1 to 3 km/s, may fall below the 0.001 ly minimum (store the real
    value, display "<63 AU") and is class U (collisional family). Let
    the star generator accept 0.075 to 0.08 Msun (the IMF floor
    `IMF_BREAKS_SOL[0]` is 0.08) or clamp merger stars to 0.08. Open
    questions for Boss (defaults taken): the physical table, with an
    `ALWAYS_DESTROY_TERRESTRIAL` setting for the literal rule (about a
    quarter of Earth-Earth hits at galactic speeds destroy both);
    hit-and-run keeps both bodies, eroded; the fusion check covers brown
    dwarf and star involvement.
    Prerequisites: GEN.109, GEN.142.
    Design: [docs/design/orbital-updates.md](design/orbital-updates.md)

- [ ] **GEN.111 Email the admin when two objects are inside each other's Hill radius**
  Boss (2026-10-03 05:38Z): "If email is configured and that happens,
  call the API to send an email to admin the exact ID and location
  information for where the objects are located." Done: with SMTP
  configured, each Hill-radius warning from GEN.109 emails the admin the
  two IDs and their positions.
  Research (2026-10-09, orbital-solvers-and-integrators.md):
  event-based, entry only, bound-or-deflecting pairs, no hierarchical
  members, a cooldown and a digest, as in
  orbital-solvers-and-integrators.md section 8; stdlib `smtplib` through
  an RQ job after USR.3. Do not mail on bare Hill-sphere entries; mail
  on collisions, star creation, brown dwarf creation and the
  gravitational-criterion warnings.
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
  Research (2026-10-09, nebula-and-asteroid-field-classes.md): the plan
  is the asteroid section of nebula-and-asteroid-field-classes.md; the
  Done clause is met pending Boss's review; the first implementation PR
  is item 1 of its PR order. The 3D belt seed must change from `belt.id`
  to the `uid` (`systemview3d.js` seeds from the database id, which
  changes when the database is rebuilt), and `phenomenonrender.py`'s
  `_NO_VIEW = {"asteroid_field"}` leaves a field page without a 3D view.
  Open questions for Boss: render exaggerated body sizes with a stated
  factor (default) or true to scale.
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
  Research (2026-10-09, anomalies.md): the Done clause can point at
  `docs/design/anomalies.md`: Tier 1 magnetar subtype, Einstein radius
  display and the pulsar fraction fix; Tier 2 Wolf-Rayet, superbubble,
  jets, symbiotic, radiation score, Fermi bubbles; Tier 3 fiction.
  Replace any instruction to apply the documents' lore multipliers to
  real kinds with the conversion rule (point-like objects times 0.2776,
  extended objects equal the filling factor). Record the document errors
  found (6 to 10 systems per sector, the black hole "0.08%" (it is 80%,
  3.3% per 20 ly sector), G mu 7.4e-10 already excluded, deficit-angle
  shear, ergosphere shear). Open question for Boss (default no): are the
  subspace and other fiction anomalies wanted? Tier 3 only through a
  default-off `lore-anomaly` kind with `K_lore` 0.

- [ ] **GEN.114 Add the chosen anomalies to the starmap**
  Done: the anomalies GEN.113 chose are generated, stored, shown on the
  maps and searchable, one kind per PR.
  Research (2026-10-09, anomalies.md): split into one PR per kind in the
  order of the ranked list; add prerequisites GEN.104 (jets) and GEN.128
  to GEN.130 (symbiotic binaries).
  Prerequisite: GEN.113.

- [ ] **GEN.128 Design: multi-star hierarchies and compact-object primaries**
  Boss (GitHub issues [#777](https://github.com/dwhagar/planetGen/issues/777) and [#778](https://github.com/dwhagar/planetGen/issues/778), 2026-10-09 06:43Z): "scientifically
  accurate star systems with up to 7 stars, this is going to be complex
  but that is the highest number of stars we've seen in orbit around
  each other" and "exotic star systems that have black holes, neutron
  stars, or similar as the central star for systems, binary systems."
  Done: a design note in docs/design covers how a hierarchy of up to
  seven stars is stored (a tree of pairs, each pair's orbit around its
  barycentre), the stability limits it must satisfy, how it fits
  GEN.62's naming and the binary code, how often each shape occurs, what
  a black hole, neutron star or similar primary changes for the planets
  around it, where such systems sit in the galaxy (GEN.103), and a go or
  no-go list for GEN.129 and GEN.130.
  Research (2026-10-09, multistar-and-compact-systems.md): the design
  note is the deliverable and it exists now. The Done text can keep its
  list and add the binary period distribution redraw, eccentric pair
  orbits, the Kozai-Lidov screen and the compact-object slices.
  `compact_remnant.py` cites a `docs/design/exotic-phenomena.md` that
  does not exist: do not create it; point the docstring at
  `docs/design/anomalies.md` and the new note.

- [ ] **GEN.129 Multi-star systems of up to seven stars**
  Boss (GitHub issue [#777](https://github.com/dwhagar/planetGen/issues/777)): "We need to add scientifically accurate star
  systems with up to 7 stars, this is going to be complex but that is
  the highest number of stars we've seen in orbit around each other."
  Done as GEN.128's design says: generation draws hierarchies of three
  to seven stars that satisfy the stability limits, stores them, names
  them (GEN.62), and the maps and the orbital update handle them.
  Research (2026-10-09, multistar-and-compact-systems.md): replace
  "three to seven stars" with "two to seven" so the new model replaces
  the two-star classes. Sub-steps: (a) the period and mass-ratio redraw
  and eccentric pair orbits on the existing two-star code (no schema
  change; ships first and alone, it changes every new binary: a 0.26 to
  50 AU range appears and the close/wide coin flip goes); (b) the
  `orbit_nodes` schema, the N=2 migration and readers and writers; (c)
  the naming extension of GEN.62 (path codes such as Aa, Ab, Ba with
  branch words at wide splits, planet stems like `Pikkita Ba I`; branch
  words come from the seeded stream in a fixed order); (d) the hierarchy
  generator for N=3, then N=4 to 7, with the Mardling-Aarseth,
  Eggleton-Kiseleva and Kozai-Lidov checks (recreate the generator from
  design 5.1 to 5.3; the research scripts live only in a session
  scratchpad); (e) the scene, maps (per-node zoom, log-radial overview)
  and the JavaScript twin; (f) the orbital update
  (`advance_orbital_phases` for nodes) and the GEN.109 hooks; (g) planet
  zones and flux sums (design 5.6, 6.5). Raise the B and A binary
  probabilities in `BINARY_SYSTEM_PROBABILITY_BY_SPECTRAL_CLASS` (0.65,
  0.55) to about 0.75 and 0.60 when the N table replaces the binary
  yes/no; the table's mean is 1.38 stars per system. Open questions for
  Boss (defaults taken): faithful rarity (about 1 in 33,000 systems have
  seven stars) with a directive to force N and a prevalence knob; the
  path-code naming; REBOUND stays out (dev script outside the package,
  leapfrog or `DOP853` regression in the suite).
  Research (2026-10-09,
  exotic-environments-planets-and-compact-binaries.md): Extend the Kozai
  screen with the relativistic-precession quench: reject the hierarchy,
  or mark it `kl_active`, only if the inclination is 39.2 to 140.8
  degrees, the Kozai-Lidov timescale is shorter than the system's age,
  and it is also shorter than the inner orbit's general-relativity
  precession period.
  Verified (2026-10-09,
  exotic-environments-planets-and-compact-binaries.md): Boss's engine
  verified the triple-stability criteria (2026-10-09): Mardling and
  Aarseth 2001 and Eggleton and Kiseleva 1995 stay the implemented
  forms. Acceptance tests: (1) P-type, m1 = m2 = 1, m3 close to 0, outer
  eccentricity 0, inclination 0 gives R_p / a_in = 2.80; (2) S-type, m1
  = 1, m2 close to 0, m3 = 1, outer eccentricity 0.5, inclination 0
  gives R_p / a_in = 2.8 x 2^0.4 x (1.5 / sqrt(0.5))^0.4 = 4.99 and
  a_out / a_in = 9.98.
  Update (2026-10-09, after Boss supplied the paper): the implementable
  Vynatheya et al. 2022 form (MNRAS 516, 4146; arXiv 2207.03151v2) is
  its Equation 4: q_out = m3 / (m1 + m2); e_in,max = sqrt(1 - 5/3 cos^2
  i); e_in,avg = 0.5 e_in,max^2; the effective inner eccentricity is
  max(e_in, e_in,avg), and equals e_in where cos^2 i > 3/5 (our choice;
  the paper is silent); Y = a_out (1 - e_out) / [a_in (1 + e_in_eff)];
  Y_crit = 2.4 [(1 + q_out) / ((1 + e_in_eff) (1 - e_out)^(1/2))]^(2/5)
  x [((1 - 0.2 e_in_eff + e_out) / 8) (cos i - 1) + 1]; stable if Y >
  Y_crit. It was fitted for 0.01 <= q_in <= 1, 0.01 <= q_out <= 100 and
  1e-4 < alpha < 1. Accuracy against N-body runs (paper Table 4): EK95
  0.86, MA01 0.90, Equation 4 0.93, the paper's neural net 0.95 (its
  weights were not reachable, so it is not implementable). Decided by
  Boss (2026-10-09 18:56Z): one function with a switch; Mardling-Aarseth
  for planets (q about 1e-3 is outside Equation 4's fitted range),
  Equation 4 for star-only triples inside the fitted range, EK95 as a
  cross-check; unit tests for both formulas. Extra tests, a_out / a_in needed (MA01
  / Equation 4): P-type coplanar circular m1 = m2 = 1, m3 close to 0:
  2.80 / 2.40; S-type m1 = 1, m2 close to 0, m3 = 1, e_out = 0.5, i = 0:
  9.98 / 7.28; q_out = 0.5, all e = 0, i = 0: 3.29 / 2.82; i = 90
  degrees: 2.80 / 3.19; i = 180 degrees: 2.31 / 2.12; q_out = 0.5, e_in
  = 0.6, i = 0: 3.29 / 3.74. Do not repeat the earlier engine
  descriptions (a constant 2.8 replaced by f(e_in), a piecewise f(i),
  97% accuracy): they were wrong. The Tory-Grishin-Mandel boundary as
  pasted contradicts the test-particle limit. Full text and transcribed
  code: exotic-environments-planets-and-compact-binaries.md section 2.3
  and /mnt/project-files/research/scripts/exotic/vynatheya2022.py.

- [ ] **GEN.130 Exotic star systems: a black hole, neutron star or similar at the center**
  Boss (GitHub issue [#778](https://github.com/dwhagar/planetGen/issues/778)): "We need to add exotic star systems that
  have black holes, neutron stars, or similar as the central star for
  systems, binary systems." Done as GEN.128's design says: generation
  can make a system or binary whose primary is a black hole, neutron
  star or similar compact object, with the planets, belts and radiation
  environment that follow; GEN.100's galaxy-wide scatter supplies the
  compact objects.
  Research (2026-10-09, multistar-and-compact-systems.md): split into
  (a) collapsed leaves in the tree with the survival check, (b) NS/BH/WD
  presets, (c) radiation fields; planets around young pulsars and black
  hole debris disks are out of scope. Compact-object habitability: the
  lowest tier for every pulsar, black hole and X-ray-binary planet, a
  normal score for cool white dwarf planets in the narrow zone (open
  question for Boss, default as stated). Open question for Boss
  (default: only through GEN.100's scatter, which promotes a scattered
  remnant to a system when it draws a bound companion or planets, at the
  design 7.3 fractions): do neutron star and black hole systems appear
  in normal fills at their real rate?
  Research (2026-10-09,
  exotic-environments-planets-and-compact-binaries.md): Millisecond
  pulsars get planets with probability 0.7% (upper bound 1%), in three
  types: disc 25%, captured circumbinary giant 25% (only in globular
  clusters), ablated companion 50%; young pulsars get none; pulsar
  planets are carbon rich with no atmosphere and radiation tier 3.
  Checks to add: reject a compact-compact binary whose
  gravitational-wave inspiral time (Peters) at its drawn separation and
  eccentricity is shorter than the system's age; allow a bound companion
  beyond about 1,000 AU around a black hole or neutron star only for
  direct-collapse black holes. Slices (a) to (c) as in
  multistar-and-compact-systems.md section 8. Open questions for Boss
  (defaults taken): a globular-cluster model does not exist, so the
  captured-giant type is left out until clusters are modelled; the 0.7%
  rate for millisecond pulsars is accepted.
  Decided (2026-10-09): Boss (2026-10-09 18:37Z): the defaults are
  accepted (millisecond-pulsar rate 0.7%, no globular-cluster model so
  the captured-giant type stays out). Where a formula is needed, use
  the verified published forms: Peters for the gravitational-wave
  inspiral time, and the Mardling-Aarseth and Eggleton-Kiseleva
  criteria for triple stability (see GEN.129).
  Research (2026-10-09, globular-clusters.md): the captured-giant type
  (type B) is enabled once clusters exist, as GEN.163; Boss's 18:37Z
  decision (left out until clusters are modelled) holds today.

- [ ] **GEN.134 Tune the star populations to the observed star-formation profile by galactic radius**
  From GEN.133 (Boss's 09:20Z request; the analysis is `docs/design/star-types-by-galactic-radius.md`, PR #807), whose four proposals were left unbuilt. Done: (1) the young and intermediate populations are weighted by the observed star-formation profile (peak at 5 kpc, about -0.28 dex per kpc beyond, a dip inside 3 kpc, the Central Molecular Zone as its own small young region); (2) the young population's B share is lowered, or its weight cut, so local B stars come to about 0.04%; (3) the bulge gets a small young tail (about 10% under 5 Gyr, between the HST and microlensing figures); (4) a metallicity gradient is added only if planet occurrence is later tied to it, otherwise the note records why not. The changes are reproducible (GEN.56), and a test compares the star type shares at the core, mid radius and rim against the note's table. Boss to confirm which of the four proposals he wants before the build starts.
  Research (2026-10-09, globular-clusters.md): item 4 now has its
  reason: planet occurrence is tied to metallicity in clusters (see
  GEN.158).
  Prerequisites: none. Related: GEN.133 (done).

- [ ] **GEN.135 Deterministic math helpers for every stored float, with a lint and frozen constants**
  Python 3.9 gives different sectors from Python 3.10 to 3.13 on all 30
  probe sectors (the research measured it); with `sqrt(x*x+y*y)` in
  place of `math.hypot` all five interpreters agree. Done:
  `planetgen/util/detmath.py` with `hypot(*c)`, `dist(a, b)` and `fsum`;
  `math.hypot`, `math.dist` and float `sum()` replaced at about 25
  stored-value sites in `galaxy/`, `population/`, `physics/` and `db/`
  (`geometry.sector_address_at` already does it); the AST scan in
  `test_reproducible_draws.py` extended to reject the old calls; the
  astropy-derived constants in `physics/constants.py` frozen as literals
  (`HYDROGEN_ATOM_MASS_KG` differs by 1.35e-9 between astropy 6.1.7 and
  8.0.1) with a test that the astropy values equal them; `@` (BLAS)
  removed from `galaxy/nebula_shape.py` (lines 107, 163) and
  `physics/sector_path.py` (line 232). It changes stored values, so it
  ships with a generator epoch bump, and it blocks TEST.77 on the Python
  3.9 CI leg. Related: GEN.66 follow-up on pinning astropy's constant
  sets.
  Prerequisites: none.
  Design: [docs/design/reproducible-galaxies.md](design/reproducible-galaxies.md)

- [ ] **GEN.136 Fingerprint encoding: floats to 9 significant digits, a stored leaf digest per sector and a ring-and-layer tree**
  GEN.58 as built hashes floats exactly; the research recommends
  rounding to 9 significant digits (generation-determinism.md section
  4). Done: the canonical encoder (floats as `format(x, ".8e")`,
  integers exact and above 2^53 as hex, no Unicode normalisation, sorted
  keys, generation order for lists, no ids or timestamps, a `pgfp1`
  format tag), SHA-256 (not RFC 8785, not BLAKE3, which needs Python
  3.11); a per-sector leaf digest stored at save time (DB.9 keeps its
  own storage checksum); a (ring, layer) node table; ring and galaxy
  roots on demand; a region digest that includes the sector count;
  separate structural and float digests over named sections. Replaces
  the `tests/galaxy_fingerprint.py` stand-in. Open question for Boss
  (default: round to 9 digits): keep exact float hashing as built, or
  round?
  Prerequisite: GEN.135.
  Design: [docs/design/generation-determinism.md](design/generation-determinism.md)

- [ ] **GEN.139 Orbit-update thresholds: per-object epoch, path-length rule and what the 0.01 mpc applies to (GEN.106 built)**
  Follow-up to GEN.106: (a) a per-object epoch column beside
  `next_update_due`; (b) define the movement rule on path length, not
  displacement (Phobos and Deimos can never move 100,000 km); (c) state
  whether 0.01 mpc applies to the galaxy-frame velocity (every star due
  every 16 days at 230 km/s, so "1 day = 1 day" runs touch about 1/15 of
  the stars, 6.2e7 writes a day per 1e9 stars) or to peculiar motion
  only; `next_update_due` can use the current speed; (d) use
  orbital-solvers-and-integrators.md section 5 as the justification; (e)
  "visible" could mean half a pixel at the deepest zoom where the object
  type is drawn (needs a map-scale value). Open questions for Boss
  (defaults taken): stars follow their analytic orbit and are rewritten
  at sector exit or when the influence set changes, planets and moons
  are computed from phase on read with thresholds only scheduling
  perturbation re-fits; the deepest star zoom is the 4 pc sector view.
  Prerequisites: none.
  Design: [docs/design/orbital-solvers-and-integrators.md](design/orbital-solvers-and-integrators.md)

- [ ] **GEN.140 Orbital math guards the edge-case table adds (GEN.108 built)**
  Follow-up to GEN.108: adopt the edge-case table of
  orbital-solvers-and-integrators.md section 7 as the checklist. Add the
  guards not in the built item: Stumpff functions by series for abs(z) <
  1, the `fmod(t, P) / P` phase form, DOUBLE-only period and phase
  columns, radial orbits (h = 0), retrograde equinoctial singularities
  (use state vectors), the r = 0 and b = 0 guards and a sub-step cap for
  the galactic centre, and a test that the acceleration at (0, 0, 0) is
  finite. The corrected universal equation is in orbital-updates.md
  10.6.
  Prerequisites: none.
  Design: [docs/design/orbital-solvers-and-integrators.md](design/orbital-solvers-and-integrators.md)

- [ ] **GEN.141 Faster Kepler solver (Mikkola or Markley) with brentq as fallback**
  Swap `solve_eccentric_anomalies` for Mikkola or Markley with `brentq`
  as fallback (20 to 30 times faster for comets). Low priority. When
  planets gain eccentricity (MAP.70 and MAP.62 are built), reuse the
  same non-iterative solver in `body_positions.py` and
  `orbitpositions.js` so the parity test stays simple.
  Prerequisites: none.
  Design: [docs/design/orbital-solvers-and-integrators.md](design/orbital-solvers-and-integrators.md)

- [ ] **GEN.142 Peculiar velocity for rogue planets and asteroid fields**
  `rogue_planets` and `asteroid_fields` have no velocity, so the
  co-rotation model gives neighbouring rogues near-zero relative speed
  and any collision test on current data finds nothing. Done: velocity
  columns on both, drawn per body at generation from the age and
  population dispersion (default thin disc (35, 25, 20) km/s per axis,
  thick disc (60, 40, 35), bulge isotropic about 110 km/s). Rogue
  positions are absolute float64 parsec values good to about 14 km at 8
  kpc, unusable for small-body collision tests; compute differences from
  sector-relative values. Open question for Boss (default yes).
  Prerequisites: none.
  Design: [docs/design/collisions-and-mergers.md](design/collisions-and-mergers.md)

- [ ] **GEN.143 A collision_events table, an admin report and a test that runs the whole collision path**
  `collision_events` table, the admin report in design section 7 and a
  "last collisions" admin view, plus a test that seeds two rogues on a
  collision course through the whole update to the field, giant merger
  and star branches. Cross-link collisions-and-mergers.md from
  nebula-and-asteroid-field-classes.md and interstellar-object-rates.md.
  Prerequisite: GEN.110.
  Design: [docs/design/collisions-and-mergers.md](design/collisions-and-mergers.md)

- [ ] **GEN.144 One distance for the Sun from the galactic centre across the constants, the density model and the design docs**
  `physics/constants.py` has 7.91 kpc (`GALACTIC_CENTER_DISTANCE_LY =
  25800`), the GEN.115 text and design docs 8.128 kpc, and
  `tuning.GALAXY_SOLAR_RADIUS_TO_SCALE_LENGTH` 8.2 kpc. Done: one named
  value for the test radius (8.128 kpc, the Gaia curve) and one for
  placing the Sun in the density model (8.2 kpc), used everywhere. Also
  correct the comment in `physics/constants.py` (lines 95 to 103)
  claiming 206 km/s and 236 Myr: the Sun's angular speed from Sgr A*
  gives about 203 Myr, and the new model's circular period at 8.128 kpc
  is 218 Myr. `MILKY_WAY_MASS = 1.15e12 Msun` against the potential's
  6.8e11 inside 100 kpc is a follow-up for
  `keepout.galactic_hill_radius_km`.
  Prerequisites: none.
  Design: [docs/design/galactic-potential.md](design/galactic-potential.md)

- [ ] **GEN.145 Class S atmosphere rule: S keeps air unless the shoreline ratio is over 30**
  Class S landed with `atmosphere = None` in the ecosphere and cold
  zones, which breaks the shoreline by a factor of 50 or more. S is
  airless only for shoreline ratio over 30 (T_eq above 650 to 1,200 K).
  It changes existing output, so it needs GEN.92's "unchanged until
  regenerated" guard.
  Prerequisite: GEN.91.
  Design: [docs/design/atmospheres-retention-and-classes.md](design/atmospheres-retention-and-classes.md)

- [ ] **GEN.146 Teff-dependent habitable zone from the Kopparapu table**
  Replace the fixed 1.1 and 0.53 flux limits in
  `orbits.calculate_habitable_zone` with the Kopparapu table (design doc
  4.1): the fixed limits put the outer edge 33 percent too near for a
  3,200 K star and at 1.37 AU against 1.68 AU for the Sun. It changes
  the zone of every star, so it needs reproducibility handling (the
  GEN.57 family) and the epoch bump.
  Prerequisite: GEN.91.
  Design: [docs/design/atmospheres-retention-and-classes.md](design/atmospheres-retention-and-classes.md)

- [ ] **GEN.148 Habitability index follow-ups from the research (GEN.84 built)**
  GEN.84 is built; the research adds: record Boss's decision as PHI-4
  domains for display plus three Xenobiology scores behind (2026-10-07
  17:11Z, confirmed 2026-10-09 07:48Z); `PHI_bio` becomes the Liebig
  blend (habitability-index.md 7.1; the settled product in the built
  item contradicts it, a point for Boss to settle), the Master doc's
  example table recomputed once the L terms exist, the display-band
  table (2.3), the equipment-tier table (4), the 16 conflicts and
  resolutions (3) and the constants (2.4). Numeric radiation tiers fill
  the gap between 100 mSv/yr and 10 Gy/yr (design doc 4.6); fix the Gy
  and Sv units in the L_rad item (Mars is 0.077 Gy/yr, 0.236 Sv/yr); use
  D10 5 to 10 kGy for Deinococcus. A summed-insolation input with a
  `flux_variation` term for multi-star systems. Expose a harshness
  number (0 to 1, default 0.5, no effect until it exists) for the tech
  draw and facility rating suggestions. The dose ladder cannot lean on
  astropy on Python 3.9 (6.0.1 has no `Gy` or `Sv`); keep Gy/yr and
  mSv/yr as two quantities. Open question for Boss (default: 20 mSv/yr
  for the Blue tier, 50 as the single-year limit).
  Design: [docs/design/habitability-index.md](design/habitability-index.md)

- [ ] **GEN.149 Planetary-nebula central stars: 0.5 to 0.7 Msun, 1e2 to 1e4 Lsun, up to 2e5 K**
  Central stars are generated at 1.1 to 1.44 Msun and at most 100 Lsun;
  real ones are 0.5 to 0.7 Msun (median 0.59), 1e2 to 1e4 Lsun and up to
  2e5 K. Done: fix the hot-white-dwarf branch of `generation/star.py`
  and `physics/constants.py` (`HOT_WHITE_DWARF_MIN_MASS_SOL`,
  `TEMP_RANGES`, `WD_LUMINOSITY_RANGE_SOL`). Affects the
  planetary-nebula scatter and `PLANETARY_NEBULA_CENTRAL_STAR_TYPES`.
  Prerequisites: none.
  Design: [docs/design/nebula-and-asteroid-field-classes.md](design/nebula-and-asteroid-field-classes.md)

- [ ] **GEN.150 H II region radius and density from the ionizing photon rate, and IMF-based nebula hosts**
  Derive the H II radius and density from the ionizing photon rate Q
  (Stromgren radius plus Spitzer expansion), add `age_myr`, `q_total`
  and the host list to the nebula record, and replace the
  `NEBULA_HOST_RULES` chances with the IMF-based cluster draw; drop B0
  to B2 as hosts of class D (keep for C and G). Today a B0 to B2 star
  has a 50% chance of a class D nebula (3 to 30 pc) though one B2V
  ionizes about 0.8 ly, and class radii are drawn independently of the
  stars. The nebula age window must join the `_bright_table` cache key
  `(threshold, population)`.
  Prerequisites: none.
  Design: [docs/design/nebula-and-asteroid-field-classes.md](design/nebula-and-asteroid-field-classes.md)

- [ ] **GEN.151 Supernova remnant sizes from the density-dependent Sedov-Taylor law (GEN.10 follow-up)**
  Replace `SEDOV_TAYLOR_RADIUS_COEFFICIENT_LY = 0.35` (an ambient
  density of about 220 cm^-3, while the remnant classes list 0.1 to 100)
  with `1.03 ly (E51/n)^0.2 (t/yr)^0.4`, capped where the snowplow phase
  begins.
  Prerequisites: none.
  Design: [docs/design/nebula-and-asteroid-field-classes.md](design/nebula-and-asteroid-field-classes.md)

- [ ] **GEN.152 Nebula cloud field is 10 to 40 times too full; lower it to the observed filling (GEN.47 rate check)**
  GEN.47's cloud field fills 3% to 46% of volume against 0.5% to 1%
  observed. Recommended about 1e-7 pc^-3 for class M with small dark
  clouds as separate rows; `GMC_ARM_FILLING_FACTOR` is still read by
  nothing. Open question for Boss (default: lower it to about 1%).
  Prerequisites: none.
  Design: [docs/design/anomalies.md](design/anomalies.md)

- [ ] **GEN.153 Magnetar subtype of neutron star, and an age-dependent pulsar fraction**
  `NEUTRON_STAR_PULSAR_CHANCE` is 0.7 against about 1e-4 to 1e-3 real.
  Done: a magnetar subtype and a pulsar fraction that follows age, with
  a `PHENOMENON_RATE_SCALE`-style override. Open question for Boss
  (default physical): keep 70% as a gameplay choice?
  Prerequisites: none.
  Design: [docs/design/anomalies.md](design/anomalies.md)

- [ ] **GEN.154 Show the Einstein radius on compact-object pages**
  Tier 1 of GEN.113: show the Einstein radius on the pages of compact
  objects.
  Prerequisites: none.
  Design: [docs/design/anomalies.md](design/anomalies.md)

- [ ] **GEN.155 A nuclear-cluster object for the Sgr A* sector (optional)**
  If Boss wants the real galactic centre: a nuclear-cluster object for
  the Sgr A* sector. Open question for Boss (default: leave the model
  alone and note it).
  Prerequisites: none.
  Design: [docs/design/fill-order-curves-and-core.md](design/fill-order-curves-and-core.md)

- [ ] **GEN.156 Pin astropy to CODATA 2018 and IAU 2015 constants before its first import (GEN.66 follow-up)**
  `requirements.lock` resolves astropy 6.0.1 (Python 3.9), 6.1.7 (3.10)
  and 8.0.1 (3.11+), and 8.0 defaults to CODATA 2022 (the proton and
  electron masses differ by about 1.4e-9). Done:
  `astropy.physical_constants` pinned to CODATA 2018 and
  `astropy.astronomical_constants` to IAU 2015 before the first astropy
  import, with tests that the constants agree to 1e-12 across the
  supported interpreters, and the `HYDROGEN_ATOM_MASS_KG` tolerance in
  `test_astropy_constants.py` tightened from 1e-4 to 1e-12 (the loose
  one hides the difference). May be done by the deterministic math
  helpers' frozen constants.
  Prerequisites: none.
  Design: [docs/design/units-and-number-formatting.md](design/units-and-number-formatting.md)

- [ ] **GEN.157 The `neighbor_galaxies` table and a verified data file of about 25 real galaxies**
  One Alembic revision, no change to existing tables, and the data file
  checked against McConnachie 2012 (VizieR J/AJ/144/4) and NED. Hold
  absolute magnitude and distance, never a hand-keyed apparent magnitude
  (the draft's recalled values were off by 0.7 to 1.9 mag). The layout
  already has Mpc and Gpc units, so a galaxy's centre fits the 76-bit
  position ID; add type code 13 `galaxy` in `names/object_id.py` only if
  the neighbours use it as their `uid`.
  Prerequisites: none.
  Design: [docs/design/multiple-galaxies.md](design/multiple-galaxies.md)

- [ ] **GEN.158 Add a metallicity value to stars**
  Research (2026-10-09, globular-clusters.md, PR #824; handoff in
  /mnt/project-files/research/handoff/globular-clusters.md): one [Fe/H]
  per star, defaulting from position (a mild disc gradient of about
  -0.06 dex per kpc, recalled and to verify; clusters carry their own).
  The first slice feeds only the planet-occurrence factor. Prerequisite
  for the cluster items. It gives GEN.134 item 4 its reason: giant
  planets are multiplied by 10^(2 [Fe/H]) in clusters. The repository's
  lifetime law gives a 0.93 Msun turn-off at 12 Gyr (the note's 0.85
  needs metallicity, which the code lacks).
  Prerequisites: none.
  Design: [docs/design/globular-clusters.md](design/globular-clusters.md)

- [ ] **GEN.159 Globular clusters: cluster table, King tables and the Milky Way catalogue**
  Research (2026-10-09, globular-clusters.md, PR #824; handoff in
  /mnt/project-files/research/handoff/globular-clusters.md): one row per
  cluster (galaxy, name, source, position and velocity in the galaxy
  frame, epoch, mass, half-mass radius, King W0, age, [Fe/H],
  core-collapse flag, tidal radius). Ten King tables for W0 = 3 to 12
  built once, resampled to about 2,000 log-grid points and scaled by
  length. Tests: the solver reproduces concentrations c = 0.67, 1.03,
  1.53, 2.12 and 2.74 for W0 = 3, 5, 7, 9 and 12; the total expected
  count is M / 0.4. Decided by Boss (2026-10-09 19:02Z, default taken): wait and build the synthetic generator first; the options were to supply the Harris catalogue (about 157
  rows: position, distance, [Fe/H], c, r_c, r_h, M_V, sigma_v,
  core-collapse flag) or approve a one-time download script. Decided by Boss (2026-10-09 19:02Z, default taken): derived: a cluster star's sector address
  is derived from the cluster centre's sector path, not stored per star.
  Prerequisite: GEN.158.
  Design: [docs/design/globular-clusters.md](design/globular-clusters.md)

- [ ] **GEN.160 Cluster density in the sector gate, with a "cluster" population**
  Research (2026-10-09, globular-clusters.md, PR #824; handoff in
  /mnt/project-files/research/handoff/globular-clusters.md): add the
  summed cluster density to relative_density and predicted_star_count
  with a spatial index (a sector tests only the clusters whose tidal
  sphere meets it), include the cluster peak in the skeleton's bound
  (_bound_raw_density_at), and add a "cluster" population (age 10 to 13
  Gyr, from the cluster row) to tuning.STELLAR_POPULATION_AGE_RANGES_GY.
  The lifetime redraw removes stars above the turn-off by itself. Tests:
  the sampled radial profile matches the table; no sector above the gate
  is skipped; counts match the note's table (a typical cluster 5.0e5
  systems in 965 sectors, a 47 Tuc-like one 2.0e6 in 5,848).
  Prerequisite: GEN.159.
  Design: [docs/design/globular-clusters.md](design/globular-clusters.md)

- [ ] **GEN.161 Bright-first fill for cluster sectors**
  Research (2026-10-09, globular-clusters.md, PR #824; handoff in
  /mnt/project-files/research/handoff/globular-clusters.md): giants,
  blue stragglers and horizontal-branch stars first, by luminosity level
  (GEN.44); the dim main sequence on demand. The cost reason: every
  system of a 157-cluster Milky Way is 1.8e8 systems, 10 to 34 days of
  single-worker time at 5 to 16.6 ms each.
  Prerequisite: GEN.160.
  Design: [docs/design/globular-clusters.md](design/globular-clusters.md)

- [ ] **GEN.162 Planet cull and blue stragglers in clusters**
  Research (2026-10-09, globular-clusters.md, PR #824; handoff in
  /mnt/project-files/research/handoff/globular-clusters.md):
  giant-planet probability multiplied by 10^(2 [Fe/H]), capped at 1;
  hard cut: drop planets with a > a_h = G M_host / sigma(r)^2 (2 AU at
  15 km/s, 4.4 at 10, 18 at 5, 440 at 1), with sigma(r) from the King
  model. Blue straggler count N_BSS ~ M_core^0.4 (main-sequence stars of
  1 to 1.7 Msun in the core). Decided by Boss (2026-10-09 19:02Z, default taken): a hard cut at a_h rather than a smooth exponential. Open question for Boss (default taken): blue stragglers
  scale as M_core^0.4 while millisecond pulsars and X-ray binaries scale
  with the encounter rate (his text says linearly with Gamma for both).
  Prerequisites: GEN.158, GEN.160.
  Design: [docs/design/globular-clusters.md](design/globular-clusters.md)

- [ ] **GEN.163 Type-B pulsar planets in globular clusters (GEN.130 follow-on)**
  Research (2026-10-09, globular-clusters.md, PR #824; handoff in
  /mnt/project-files/research/handoff/globular-clusters.md): enable
  GEN.130's captured-circumbinary-giant type only for millisecond
  pulsars inside a cluster: pulsar plus white-dwarf companion plus
  circumbinary giant, exempt from the metallicity cull. Rate 0.7% x 25%
  = 0.175% of cluster millisecond pulsars (PSR B1620-26 b: 1 of about
  340 = 0.29%, computed). Pulsar counts per cluster scale with the
  encounter rate Gamma ~ rho_c^2 r_c^3 / sigma; default 10 times the
  field rate per star, by Gamma rank.
  Prerequisites: GEN.130, GEN.158, GEN.159.
  Design: [docs/design/globular-clusters.md](design/globular-clusters.md)

- [ ] **GEN.164 Synthetic globular-cluster systems for generated galaxies**
  Research (2026-10-09, globular-clusters.md, PR #824; handoff in
  /mnt/project-files/research/handoff/globular-clusters.md): for the
  galaxies of GEN.9 and VIEW.2: N_GC = S_N 10^(-0.4 (M_V + 15)) from
  each galaxy's M_V; S_N by morphology (spirals 0.5 to 1.5, ellipticals
  2 to 8, central giants 8 to 15, recalled); a lognormal mass function
  with peak 2e5 Msun and width 0.5 dex (recalled, not in Boss's text),
  truncated to 1e4 to 3e6; [Fe/H] a two-Gaussian mixture (-1.5 and -0.5,
  width 0.3, recalled) truncated near [-2.4, 0.0]; metal-rich positions
  drawn from the bulge plus thick-disc density, metal-poor from a
  spherical r^-3.5 profile with a core of a few kpc. Open question for
  Boss (default: 0.5 dex and the S_N ranges above): the mass-function
  width and S_N ranges.
  Prerequisites: GEN.9, GEN.158, GEN.159.
  Design: [docs/design/globular-clusters.md](design/globular-clusters.md)

- [ ] **GEN.169 Decide the phenomenon scatter rates: regional factors and the 0.1% intermediate-mass black holes**
  Research (2026-10-09, phenomenon-scatter-mass-cut.md, PR #837; handoff
  in /mnt/project-files/research/handoff/phenomenon-mass-cut.md): Boss
  decided at 19:54Z that the GEN.100 scatter keeps only objects above a
  lowest mass and the sector fill draws the rest below it, like the
  bright stars. Open question for Boss (default: keep both as they are):
  the regional factors give neutron stars 0.80 times and black holes
  1.49 times their nominal numbers (his retune said 1e9 neutron stars
  and 1e8 black holes; the scatter gives 9.5e8 and 2.2e8). Renormalise
  the factors so the totals match his figures? And is the 0.1% share of
  intermediate-mass black holes intended? They carry 53% of the
  black-hole mass.
  Note (2026-10-09): Superseded in part (Boss, 2026-10-10 03:01Z rush
  job): the mass cut is now a user preset between 8 and 20 solar masses
  (GEN.183, default 20) and the scatter runs in the five passes of
  GEN.185; the rates above are still open.
  Prerequisites: none.
  Design: [docs/design/phenomenon-scatter-mass-cut.md](design/phenomenon-scatter-mass-cut.md)

- [ ] **GEN.170 Object ID layout: an 80-bit ID of birth sector, serial and body number, with pack, unpack, format and parse functions**
  Source: docs/design/object-id-options.md section 0 (Boss decided
  2026-10-09 22:39Z: birth location plus serial, galaxy-wide; the same
  length for every object; always identifies that one object; up to 128
  bits but shorter preferred; fix the deficits; no backward
  compatibility). Filed from the object-ID research thread. Nothing is
  built until Boss asks.
  Done: pure functions that pack, unpack, format and parse the 80-bit ID
  (20 hex digits, stored as BINARY(10)): a 40-bit birth sector address
  (ring 12 bits, biased layer 12 bits, slot 16 bits), a 28-bit serial in
  that sector (top two bits 00 generated rank, 01 added at run time, 10
  field-drawn) and a 12-bit body number (0 is the top-level object; 1
  and up number the stars, planets, moons, belts and comets of that
  system on one counter). The kind is not in the ID. Field widths are
  chosen at plan time and stored in galaxy_shape; they grow to 96 or 128
  bits only when a galaxy's bounds do not fit. Replaces the hash in
  galaxy/uid.py. Tests for fixed length and round trip. As a generator
  change it bumps the generator version (OPS.37) when it ships.
  Decided (Boss, 2026-10-09 22:44Z): 80 bits, 20 hex digits, with system
  and body fields, not a 64-bit flat per-sector counter.
  Prerequisites: none.
  Design: [docs/design/object-id-options.md](design/object-id-options.md)

- [ ] **GEN.171 The sector fill gives object IDs by generation rank**
  Source: docs/design/object-id-options.md section 0 (Boss decided
  2026-10-09 22:39Z: birth location plus serial, galaxy-wide; the same
  length for every object; always identifies that one object; up to 128
  bits but shorter preferred; fix the deficits; no backward
  compatibility). Filed from the object-ID research thread. Nothing is
  built until Boss asks.
  Done: the fill assigns generated serials by generation rank inside the
  sector (replaces `_UidIssuer` hashing and the `assign_uids` rank
  recount); one worker owns a sector, so no coordination and no database
  read. Re-measure the insert rate on the real fill order and on MySQL
  8.4 and MariaDB 11.4 (research model: 63,000 to 67,000 rows a second
  against 46,000 for the hash, one run on MariaDB 10.11).
  Prerequisite: DB.20.
  Design: [docs/design/object-id-options.md](design/object-id-options.md)

- [ ] **GEN.172 Run-time births get object IDs from the counters**
  Source: docs/design/object-id-options.md section 0 (Boss decided
  2026-10-09 22:39Z: birth location plus serial, galaxy-wide; the same
  length for every object; always identifies that one object; up to 128
  bits but shorter preferred; fix the deficits; no backward
  compatibility). Filed from the object-ID research thread. Nothing is
  built until Boss asks.
  Done: an admin-added body or system, an ejected object, a merger
  remnant, split fragments and a stand-alone facility take their IDs
  from the `id_counters` allocator. An ejected planet keeps its ID and
  becomes a rogue-planet row; a merge keeps the heavier body's ID and
  retires the other; a split gives each fragment a new run-time serial;
  a deleted ID is never reused.
  Open question for Boss (default yes): an ejected planet keeps its ID?
  Prerequisite: DB.20.
  Design: [docs/design/object-id-options.md](design/object-id-options.md)

- [ ] **GEN.173 Deleting a body and then adding one fails with IntegrityError 1062 on uq_planets_uid (bug)**
  Source: docs/design/object-id-options.md section 0 (Boss decided
  2026-10-09 22:39Z: birth location plus serial, galaxy-wide; the same
  length for every object; always identifies that one object; up to 128
  bits but shorter preferred; fix the deficits; no backward
  compatibility). Filed from the object-ID research thread. Nothing is
  built until Boss asks.
  Reproduced on MariaDB with `save_system_edits` plus `assign_uids`: the
  new row is ranked by position among the current rows and takes a uid
  already used. `add_system_to_sector` runs the sector-wide pass (by
  code reading). Fixed by the layout, fill and run-time birth items;
  keep a regression test.
  Prerequisites: GEN.171, GEN.172.
  Design: [docs/design/object-id-options.md](design/object-id-options.md)

- [ ] **GEN.174 Bodies an admin adds are saved with a NULL uid (bug)**
  Source: docs/design/object-id-options.md section 0 (Boss decided
  2026-10-09 22:39Z: birth location plus serial, galaxy-wide; the same
  length for every object; always identifies that one object; up to 128
  bits but shorter preferred; fix the deficits; no backward
  compatibility). Filed from the object-ID research thread. Nothing is
  built until Boss asks.
  `db/edits.py` `save_system_edits` goes through `store.insert_planet`,
  `insert_moon` and `insert_belt`, and nothing assigns a uid. Fixed by
  the run-time birth item.
  Prerequisite: GEN.172.
  Design: [docs/design/object-id-options.md](design/object-id-options.md)

- [ ] **GEN.176 A nebula or remnant is born in the sector holding the centre of the space it occupies**
  Source: docs/design/object-id-options.md section 0 (Boss decided
  2026-10-09 22:39Z: birth location plus serial, galaxy-wide; the same
  length for every object; always identifies that one object; up to 128
  bits but shorter preferred; fix the deficits; no backward
  compatibility). Filed from the object-ID research thread. Nothing is
  built until Boss asks.
  Done: the birth sector of a nebula or supernova remnant is the sector
  holding the geometric centre of the space it occupies: the centroid of
  the interior of its metaball field on the fixed 24-cell grid, in
  integer arithmetic, rounded to 1 mpc before the sector is taken (for a
  remnant, the centre of the shell), not the shape origin. Replaces
  "first sector saved that the cloud reaches" in
  `_insert_field_nebulae`, so the home sector no longer depends on save
  order. The field-drawn serial is the cloud's rank among the clouds of
  its field cell (`nebula_field.cell_clouds`) whose centroid is in that
  sector, so no shared counter is needed and the centre sector need not
  exist yet. Test the sector-face tie case. The same rule holds for any
  object stored by a sector other than its own.
  Open question for Boss (default 1 mpc): the rounding used before the
  birth sector is taken from the centroid?
  Prerequisite: DB.20.
  Design: [docs/design/object-id-options.md](design/object-id-options.md)

- [ ] **GEN.177 Planetary magnetic fields: a stagnant-lid factor**
  Left over from GEN.86 (PR #908, Foundations lane 2): the planetary
  magnetic field model does not yet apply the stagnant-lid factor (a
  planet with a single rigid lid and no plate tectonics cools its core
  differently, which changes whether a dynamo runs). Done: the factor is
  applied where the field strength is computed, with a test on a
  stagnant-lid and a plate-tectonic planet. Decided with GEN.86 (Boss,
  decision card, 2026-10-10): unlocked rocky planets draw their day
  length from 8 to 48 h log-uniform ("Realistic days"), shipped as a
  small follow-up patch without an item.
  Prerequisites: none. Related: GEN.86, GEN.88.

- [ ] **GEN.178 Magnetic fields: the induced field of an ocean moon**
  Left over from GEN.86 (PR #908): the induced magnetic field of a moon
  with a subsurface ocean in its parent's changing field is not
  modelled. Done: the induced field is computed for ocean moons and
  stored beside the intrinsic field, with a test on a known case (a
  Europa-like moon of a gas giant).
  Prerequisites: none (GEN.88 merged, PR #963). Related: GEN.86.

- [ ] **GEN.179 Store each sector's generation directive and attempt record with the sector**
  Left over from GEN.96 (PR #915, Foundations lane 1). Done: the
  directive a sector was generated under and the record of its attempts
  are saved with the sector (a schema migration, numbered after the
  migrations already in flight), so the sector page and the check can
  show what was asked and what was drawn.
  Prerequisites: none. Related: GEN.96, GEN.97.

- [ ] **GEN.180 Directives: a forced fill (met_forced) after K failed draws**
  Left over from GEN.96 (PR #915, Foundations lane 1). Done: when a
  directive is not met after K draws (K set in tuning.py), the generator
  takes the fallback fill marked `met_forced` instead of failing or
  looping, and the result says it was forced; a test covers a directive
  that cannot be met by drawing.
  Prerequisites: none. Related: GEN.96.

- [ ] **GEN.181 Directives: refuse impossible requests up front from the compound-Poisson tables**
  Left over from GEN.96 (PR #915, Foundations lane 1). Done: before any
  drawing, a directive whose requested mix cannot occur in the sector
  (judged from the compound-Poisson tables) is refused with the reason,
  in the Generate page and the CLI.
  Prerequisites: none. Related: GEN.96, GEN.180.

- [ ] **GEN.186 Random neighborhoods: an option to keep away from filled space**
  Left over from GEN.97 (merged, PR #950; Foundations lane 1 report,
  2026-10-10 03:37Z). The research note
  (sampling-backfill-and-resume.md) made it optional that the qualifying
  list excludes centres near filled space; neighborhoods already avoid
  each other always. Done: the Generate page and generate.py take a flag
  that drops random neighborhood centres whose neighborhood touches
  sectors that are already filled, with a cost line that says how many
  centres qualify.
  Prerequisites: none. Related: GEN.97, ADM.28.

- [ ] **GEN.187 Bright-star back scatter by mass: rings of 1, 2, 5 and 8 solar masses around filled space**
  Boss (GitHub issue
  [#952](https://github.com/dwhagar/planetGen/issues/952), 2026-10-10
  03:43Z): "Back-scatter of bright stars needs to be changed to filtered
  by mass. Within 1 sector on all sides (orthogonal only, no diagonals)
  should fill to: 1. 1 Solar Masses (1 Sector around filled region) 2. 2
  Solar Masses (1 Sector around step 1) 3. 5 Solar Masses (1 Sector
  around step 2) 4. 8 Solar Masses (1 Sector around step 3)." Done: the
  back scatter that fills around generated sectors (GEN.30's luminosity
  tiers) is replaced by four rings, each one sector further out by face
  adjacency only (no diagonals): the first ring fills stars down to 1
  solar mass, the second down to 2, the third down to 5 and the fourth
  down to 8; beyond the fourth ring only the scatter's own mass limit
  (GEN.183) applies. Decided (Boss, 2026-10-10 04:06Z, "your defaults are confirmed"): the nearest ring
  takes the lowest mass cut, as above, and each ring is counted from the
  previous ring's outer edge; the GEN.30 luminosity tiers go away.
  Prerequisites: none (GEN.184 merged, PR #965). Related: GEN.30, GEN.40, GEN.99, GEN.183,
  GEN.184, MAP.120, PERF.18.

- [ ] **GEN.188 Mass limit default is 8, and the mass slider and luminosity dropdown sit side by side on Generate and New galaxy**
  Decided (Boss, 2026-10-10 04:38Z, to Bugfixes lane 2): change the
  Generate-galaxy defaults so the mass limit default is 8 solar masses
  (GEN.183 shipped 20); put the mass slider and the luminosity floor
  dropdown (GEN.184, default 3000 L_sun) next to each other on the
  Generate page, and show both in the New galaxy section as well. Done:
  a fresh galaxy plan, the CLI default and the stored default use 8, the
  two controls sit together in both places with the same presets, and
  tests cover the default and the controls. A galaxy that stored 20
  keeps its stored limit.
  Prerequisites: none. Related: GEN.183, GEN.184, GEN.185.

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
  Plan (2026-10-07): The backfill blocks become RQ jobs (PERF.24).
  Research (2026-10-09, sampling-backfill-and-resume.md): pre-warm or
  persist the band-fraction tables (`star_population._bright_table`, an
  `lru_cache`): a forked RQ work horse rebuilds them on every job (5 s
  with a 1000 L_sun ceiling, 40 s without). Make the task unit a chunk
  of at least about 10^4 cells. The CPU gain is small (draws are about
  0.06 ms per cell; SQL and neighbour locks dominate, 1.87 times at 2 to
  4 workers). Discovery enumerates a 60-fold redundant list (531 calls
  and 1.0 million enumerated cells for a 50 pc core, about 70 s
  extrapolated for 200 pc): enumerate once around the run's boundary or
  bounding disc and keep cells within the tier radius of the region.
  Research (2026-10-09, generation-performance-study.md): the remark
  that the CPU gain is small is wrong for the scatter: 95% of its CPU
  was per-job import and band-table rebuild, and pre-warming the worker
  took 322 s to 39.5 s (PERF.42).

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
  Plan (2026-10-07): Plans on RQ job results and the
  cachetools caches (PERF.24, PERF.25).
  Research (2026-10-09, performance-eta-queue-and-caching.md): the
  planning is done by design doc 4.4: do not cache generation jobs;
  cache estimates, tile and stage builds, stats and search aggregates
  keyed by (kind, argument digest, epoch, speed-stats version); a Redis
  epoch counter incremented by writers; single-flight keys; no RQ round
  trip for cache reads (2.4 s start-up against microseconds); RQ does
  not deduplicate by job id, so dedupe with `SET NX EX`. Build items:
  (a) the epoch counter, (b) single-flight for tiles and stages, (c) a
  cached estimate, (d) XFetch only for age-expired costly entries.

- [ ] **PERF.29 Record which runs a partly filled sector still needs**
  Boss (2026-10-07 11:47Z): "Add a way to record and resume which
  sectors might be only partially filled by identifying what runs
  haven't been made yet." Done: each sector records which steps of its
  fill have finished (systems, phenomena, bright-star levels,
  population), so a check lists every partly filled sector and the steps
  it still needs; the level rule of GEN.44 (-1 untouched, a positive
  L_sun, 0 generated) stays the summary.
  Research (2026-10-09, sampling-backfill-and-resume.md): specify a
  `sector_pending(ring, layer, slot, mask, run_id)` table whose row is
  inserted in the sector's own transaction and whose bits are cleared by
  each step in that step's transaction; name the two steps whose "which
  sectors" is in memory only today (paths via `sector_ids_since`, and
  the backfill around the run via `sector_centers_since(started_at)`).
  Derive steps that have a natural artifact (bright level, population
  watermark) instead of duplicating them.

- [ ] **PERF.30 Finish an interrupted block or sector run on the next start**
  Boss (2026-10-07 11:47Z): "Add a mechanism to finish generating a
  block or sector on next start." Done: when the site or `generate.py`
  starts, unfinished runs found by PERF.29 are offered (web) or listed
  with a `--resume` flag (CLI) and finish from the step they reached,
  with the same result as an uninterrupted run.
  Research (2026-10-09, sampling-backfill-and-resume.md): build
  `--resume` on `generation_runs` plus `sector_pending`; no new outcome
  value is needed (`generation_runs` already stores the galaxy and run
  seeds, `generation_run_arguments` the command line, and a killed run
  leaves `finished_at` NULL; outcomes today are `ok`, `failed`,
  `interrupted`). Add failure-mode tests: kill between sector commit and
  settle, kill inside a backfill chunk, kill a forked work horse, a
  duplicate task, a clock change. Check whether the uid claim and
  name-registry confirmation (#657) are inside the sector save
  transaction; if not, a crash leaves an orphan reservation.
  Prerequisite: PERF.29.

- [ ] **PERF.31 Investigate: where generation spends its time, from the plan to a finished galaxy**
  Boss (GitHub issue [#761](https://github.com/dwhagar/planetGen/issues/761), 2026-10-09 04:43Z): "We need to do a timed
  analysis in full debug from the planning of the galaxy to the galaxy
  being ready and generation finished. We need to know where the system
  spends the most time in each phase." and (issue [#750](https://github.com/dwhagar/planetGen/issues/750)) "there should be
  a method to benchmark." Done: a repeatable benchmark command runs a
  small galaxy from the plan to the finished fill in full debug and
  prints the time spent in each phase and sub-phase; a report names the
  largest costs and what to do about each, and files an item per fix.
  PERF.34 and PERF.32 build on it.
  Research (2026-10-09, performance-eta-queue-and-caching.md): specify
  `planetgen benchmark [--sectors N] [--profile]` as in design doc 5.3
  (wall-clock `perf.phase()` timers always on, database counters before
  and after, `pyinstrument` only with `--profile`, runs at 1, 2 and 4
  workers and with a page-request thread). Tools for Python 3.9:
  cProfile, pyinstrument 5.1.3 (cp39 wheels), py-spy 0.4.2. In-memory
  generation is about 5 ms per system, 16.6 ms per system in a filled
  sector, so the database path (neighbour lock, uid pass, registry
  upsert) is the target.
  Research (2026-10-09, generation-performance-study.md): this study
  answers it for the scatter and fill phases (see PERF.42, PERF.43,
  PERF.44, PERF.45, PERF.46, PERF.47 and DB.19). Open for Boss: which
  phase took the 10 hours (the log will say).
  Measurement (2026-10-09, corrected): Bugfixes lane 1 first reported
  that one dense core sector took 84 s with 30 s in
  `reserve_system_names`; that first profile ran while another run was
  saving into the same database. Foundations lane 1 re-measured (PR #870,
  PERF.49) with nothing else writing: reservation itself is 0.12 s a
  sector, but at 4 workers the name registry's row locks, held until the
  sector commit, made a save wait 6 to 25 s. Names are now claimed in a
  short transaction of their own; 8 core sectors on 4 workers take 54 s
  against about 90 s.

- [ ] **PERF.33 Progress bars and ETAs from measured performance**
  Boss (GitHub issue [#661](https://github.com/dwhagar/planetGen/issues/661)): "time remaining on all progress bars should
  be calculated from this performance metric averaged with the actual
  performance at the time" and "if it is expected to take longer than 15
  seconds to complete, it gets a progress bar for that sub-task from the
  main task. This includes generating stars inside a layer, star systems
  in a sector, etc." Done: every job's expected time is its recorded
  rate (PERF.32) averaged with the live rate; a sub-task expected to
  take over 15 seconds gets its own bar under the main one; the Queue
  page ETA (ADM.39) and the banner of UX.3 use the same estimate.
  Research (2026-10-09, performance-eta-queue-and-caching.md): replace
  "recorded rate averaged with the live rate" by the estimator in design
  doc section 2.3: a ratio of two decayed sums weighted by each unit's
  predicted cost (the PERF.3 density sum), blended with the recorded
  rate by weight n/(n+15) on the live side, tau = max(60 s, 20 x mean
  task seconds / workers). Simulated error about 12% against about 50%
  for the present count-based EWMA, worst case 26% against 181%. Done
  also: a multi-step job adds each not-yet-started step's PERF.3
  estimate (`work._roll_up` adds nothing for such steps today and drops
  them from the parent's sum); the web shows a range or "estimating",
  never `-:--:--`; one pure function serves the terminal bar, the Queue
  page, the banner and DB.15.
  Left over (2026-10-09): left over from UX.83 (PR #895, Bugfixes lane
  1): the phenomenon scatter, the neighbour-linking steps and the
  population pass now draw their own bars; a bar inside one sector's
  save (the slowest sub-step in a dense sector) is now PERF.50.
  Prerequisite: PERF.32. Related: PERF.51, UX.84, PERF.50, DB.15.
  Lane (Boss, 2026-10-09 23:38Z): Bugfixes lane 1, order PERF.32, PERF.33, PERF.51, PERF.50, DB.15, then UX.84.
  Built with PERF.32 (PR #905): the control schema is v12; the
  recorded kinds are sector, scatter and phenomena; the `sums` and
  `rate_cost_per_s` columns from the research sketch were not added and
  belong to this estimator; the `bench:` prefix is reserved for PERF.31's
  rows.
  Partly built (PR #910, Bugfixes lane 1): one `blended_rate`,
  `time_constant` and `eta_range` in `queue/progress_rate.py` (the
  recorded rate of PERF.32 blended with the live rate, weight n/(n+15),
  held back for 5 units or 20 s, tau = max(60, 20 x mean task seconds /
  workers), a stall hold); the sector and bright-star bars start from
  the recorded rate; the Generate and Queue pages show a range or
  "estimating", never dashes; the Queue page ETA blends the recorded
  task seconds. Remaining in this item: (1) a multi-step job adds each
  not-yet-started step's estimate, which needs the estimate stored on
  the step node when the job is planned; (2) the UX.3 banner reads the
  same `blended_rate`; (3) the phenomenon scatter bar has no recorded
  rate to start from (its recorded units are placed phenomena, the
  bar's are weights), so it needs a recorded rate in the bar's own
  units, or it keeps the live rate. The migration bar (DB.15) stays its
  own item.

- [ ] **PERF.35 An interval or chunk ledger for untouched sectors once block-first backfill lands**
  Replace the one-`sector_stats`-row-per-visited-cell ledger of
  untouched sectors with an interval or chunk ledger once GEN.42's
  block-first draw lands; at 24 billion cells the ledger is 24 billion
  rows.
  Prerequisite: GEN.42.
  Design: [docs/design/sampling-backfill-and-resume.md](design/sampling-backfill-and-resume.md)

- [ ] **PERF.36 Memory and request guard: never list more than about 50,000 candidate cells, and refuse huge enumerations in a web request**
  Stream the enumeration instead of materialising it; refuse to
  enumerate more than about 2 million cells in a web request. Add
  sampling-based estimation above about 20,000 sectors to the PERF.3
  estimate (200 samples land within 1% of the full 7.455 million-system
  core); the span and core modes need it.
  Prerequisites: none.
  Design: [docs/design/fill-order-curves-and-core.md](design/fill-order-curves-and-core.md)

- [ ] **PERF.38 Cache fixes for the Galaxy Map under a fill: single-flight tile builds, a busy rule for the page cache, a deletion epoch in place of COUNT(*)**
  PERF.34 is built; the research reorders its suspects by evidence
  (design doc 5.4): no single-flight on tile builds; the page cache is
  emptied at every stamp check during a fill; the stamp's linear
  `COUNT(*)` in `db/query.py` `galaxy_content_state` (0.13 s per million
  placed sectors, 2.7 s at 20 million, per web process per database
  every 15 to 60 s); five API threads held 2.4 s or more by queued
  waits. Done: `busy` handling in `pagecache.py` like `tilecache.py`; a
  Redis `SET NX EX` single-flight around tile and stage builds; a
  deletion epoch counter in place of the `COUNT(*)`. The first four can
  be fixed before the benchmark. Re-run PERF.34's page-time test under
  Ludicrous Speed. The `"""int: How long a database's stamp is trusted
  ..."""` docstring in `web/lib/tilecache.py` sits after
  `FAILED_CHECK_RETRY_SECONDS` instead of under `STAMP_TTL_SECONDS`.
  Prerequisites: none.
  Design: [docs/design/performance-eta-queue-and-caching.md](design/performance-eta-queue-and-caching.md)

- [ ] **PERF.39 Every API job costs 2.4 s and 195 MB: import lazily and cap the burst workers**
  `queue/api_jobs.py` `execute` imports `planetgen.api.common` (2.3 s)
  before running anything, and `submit` starts one uncapped burst worker
  per job. Fix: catch refusals via a light `Refused` base class or a
  lazy lookup; import `nltk` and `scipy.stats` lazily; start a worker
  only when fewer than `worker_count()` are alive, on shared queues.
  Required before Boss's batch wiki uploads. Reuse one Redis connection
  in `api_jobs.status` and `wait`.
  Research (2026-10-09, generation-performance-study.md): the same fix
  as PERF.42 (warm the worker before the fork) removes about 2 to 3 s
  per queue job.
  Prerequisites: none.
  Design: [docs/design/performance-eta-queue-and-caching.md](design/performance-eta-queue-and-caching.md)

- [ ] **PERF.40 Two shared queues, a reserved interactive worker and a real "cancel now"**
  `planetgen-interactive` (one reserved worker) and `planetgen-bulk`
  (workers serve `[interactive, bulk]`); "cancel now" for running tasks
  via `send_stop_job_command` (safe: sector and backfill-chunk saves are
  single transactions; `cancel()` alone leaves a running job running);
  keep the interval-less `Retry(max=1)` and add no retry intervals
  without a scheduler (a retry with an interval stayed `scheduled` on a
  burst worker). Open question for Boss (default yes): a reserved
  interactive worker costs one process (about 195 MB) while any bulk run
  is active; start it on demand and exit when its queue is empty.
  Prerequisite: PERF.39.
  Design: [docs/design/performance-eta-queue-and-caching.md](design/performance-eta-queue-and-caching.md)

- [ ] **PERF.41 Stamp the tile and page caches with a TILE_FORMAT constant instead of the version (optional)**
  A `TILE_FORMAT` constant in the stamp base in place of `__version__`,
  with a golden test that fails when the tile or page structure changes
  without a bump. Releases are 26 to 123 a day, so a development server
  never keeps a cache. Open question for Boss (default): decide after
  PERF.31 shows the cold rebuild cost; switch if a full rebuild takes
  more than a minute.
  Prerequisite: PERF.31.
  Design: [docs/design/performance-eta-queue-and-caching.md](design/performance-eta-queue-and-caching.md)

- [ ] **PERF.46 Planets and moons: set the position once per body**
  Research (2026-10-09, generation-performance-study.md, PR #835;
  handoff in
  /mnt/project-files/research/handoff/generation-performance.md; from
  Boss's requests of 19:08Z and 19:21Z, generation being his slowest
  point): `SpatialPosition3D._sync` was called 83,000 times for 4
  sectors. Set a planet's or moon's position once per body. Optional:
  skip `util/checks.finite_domain` in bulk fills (3 to 4%, but it loses
  a safety net). Open question for Boss (default: keep the check):
  accept skipping it in bulk fills?
  Boss (2026-10-09, generation-performance-study.md): Boss (2026-10-09
  20:02Z) approved the position-once saving. The `finite_domain` half
  stays open: it costs 2 to 4% of a fill (0.85 microseconds a call,
  about 535 calls a system) and Research Lane 1 recommended keeping it;
  Boss's answer is pending, so do not skip the check yet.
  Built (PR #863, 2026-10-09): the position half; bodies work out their
  coordinates when first read. Only the `finite_domain` question is
  open, on Foundations lane 2.
  Prerequisites: none.
  Design: [docs/design/generation-performance-study.md](design/generation-performance-study.md)

- [ ] **PERF.47 The PERF.31 benchmark records the buffer pool, table sizes and worker start-up cost**
  Research (2026-10-09, generation-performance-study.md, PR #835;
  handoff in
  /mnt/project-files/research/handoff/generation-performance.md; from
  Boss's requests of 19:08Z and 19:21Z, generation being his slowest
  point): record `innodb_buffer_pool_size`, the table sizes and the
  worker start-up cost with every benchmark run, so a result can be
  compared with the next.
  Prerequisite: PERF.31.
  Design: [docs/design/generation-performance-study.md](design/generation-performance-study.md)

- [ ] **PERF.48 Low priority: a numeric-only INSERT formatter or C driver for bright_stars and phenomenon_scatter**
  Research (2026-10-09, generation-performance-study.md, PR #835;
  handoff in
  /mnt/project-files/research/handoff/generation-performance.md; from
  Boss's requests of 19:08Z and 19:21Z, generation being his slowest
  point): 13.8 down to 10.1 microseconds a row. Low priority.
  Prerequisites: none.
  Design: [docs/design/generation-performance-study.md](design/generation-performance-study.md)

## DB: Database and schema

DB.1 shipped in 7.35.0 (PR #152). DB.2 to DB.5 done (PR #342, PR #347).

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
  Research (2026-10-09, db-check-and-parity-repair.md): (1) replace "the
  canonical hash of GEN.58's fingerprint" with a separate storage
  checksum: SHA-256 over a lossless, id-preserving dump of the stored
  rows, read back after the write. (2) The checksum and the parity
  export omit the columns the daily orbit update rewrites (phases,
  positions, velocities, reflex offsets, binary positions, comet
  anomalies, and a placed object's `sector_id`; the list is in the
  design doc section 5.2). (3) Groups are built from sectors of similar
  export size, far apart in sector id, for an overhead of about 7
  percent. (4) The codec is an in-tree stdlib GF(256) Reed-Solomon
  module (`planetgen/util/rscodec.py`), no new dependency; default G =
  32, m = 2, both in `config.json`. (5) Add the delta update, the record
  format (per-record CRC, per-member checksums) and the stale-record
  decision table. (6) Repair rebuilds the table when a clustered-index
  page is unreadable (dump survivors with `innodb_force_recovery=1`,
  recreate, reload with FK checks off, insert rebuilt sectors with their
  ids). (7) The regenerate-from-seed fallback is split off as the new
  fallback item; this item needs DB.8 only. (8) The test flips bytes in
  an `.ibd` page of a copy while the server is stopped, not `UPDATE`s
  rows. (9) A `dirty` flag on the sector row, set by every writer of
  exported columns (`store.py` direct UPDATEs, `api/edits.py`,
  `db/edits.py`) and cleared after the parity update; listing those
  writers is part of the item. Before DB.19's mass cut `phenomenon_scatter`
  would be 1.17e9 rows and about 161 GB (the earlier 1.6e8 rows and 21 GB
  was an unverified estimate); after the cut it is about 2.7e5 rows. Open questions for Boss (defaults taken): after a repair
  the orbit-updated columns come back as generated and the next orbit
  update carries on (a small daily snapshot of ids, phases and sector
  only if exact positions matter); G = 32, m = 2; repairing a bad
  clustered-index page may mean restarting the database server with
  `innodb_force_recovery=1`.
  Prerequisite: GEN.57.
  Design: [docs/design/reproducible-galaxies.md](design/reproducible-galaxies.md)

- [ ] **DB.16 Store the generator epoch and run id on each sector instead of four version text columns**
  DB.7 stored `version_key`, the version, python and platform as text on
  every sector; at 12 billion sectors those strings repeat. Store
  `generator_epoch` (SMALLINT) and `run_id` (foreign key to the existing
  `generation_runs`) instead, with the key text looked up from the run.
  Open question for Boss (default: epoch plus run id): or keep the full
  key text per sector as DB.7 shipped. ALGORITHM=INSTANT, idempotent
  migration, per db-check-and-parity-repair.md.
  Prerequisite: OPS.28.
  Design: [docs/design/reproducible-galaxies.md](design/reproducible-galaxies.md)

- [ ] **DB.17 Repair by regenerating a damaged sector from its seed when parity cannot rebuild it**
  The second half of DB.9: when the parity file cannot rebuild a sector,
  regenerate it from the galaxy seed, the sector's stored directives and
  the settings file (ADM.18), then replay its edits from the edit log.
  (Boss, 2026-10-09 20:42Z, dropped the pending-delta JSON and daily
  merge this item used to read; the seed is still used internally for
  repair. Boss, 2026-10-09 20:52Z: keep the repair from the seed, replaying
  the edit log.) Compare the result with the stored leaf digest. Open question for Boss (default
  yes): ship the parity half first, without GEN.57, GEN.58 and OPS.14.
  Note (2026-10-09): admin regenerate of a planet, moon or belt
  (`admin/edits.py`) draws from the process stream, so an edited sector
  cannot be rebuilt from its seed alone, which is why this repair
  replays the edit log; the web handlers were not checked (Research Lane
  1).
  Prerequisites: DB.9, GEN.57.
  Design: [docs/design/db-check-and-parity-repair.md](design/db-check-and-parity-repair.md)

- [ ] **DB.18 Migration helpers for slow DDL: online indexes, instant columns and batched updates**
  `planetgen/db/migrations/helpers.py` with `add_index_online`,
  `add_column_instant`, `batched_update` (keyset-paged, commit per
  batch, progress callback), `exists_*` checks and a `lock_wait_timeout`
  context. Every slow DDL is written through `op.execute` with explicit
  `ALGORITHM=` and `LOCK=NONE`; never `op.batch_alter_table` on the
  large tables. Budget 11 to 31 minutes per index per 1.6e8 rows (a
  sizing figure, not the scatter's row count).
  Revision 0063 (DB.13) builds five indexes with plain `CREATE INDEX`
  and a per-row INSERT loop: fine today, not a pattern for the 10^8-row
  tables.
  Design: [docs/design/db-check-and-parity-repair.md](design/db-check-and-parity-repair.md)

- [ ] **DB.20 Object IDs in the schema: uid becomes BINARY(10), unique on its own, plus an id_counters table**
  Source: docs/design/object-id-options.md section 0 (Boss decided
  2026-10-09 22:39Z: birth location plus serial, galaxy-wide; the same
  length for every object; always identifies that one object; up to 128
  bits but shorter preferred; fix the deficits; no backward
  compatibility). Filed from the object-ID research thread. Nothing is
  built until Boss asks.
  Done: one Alembic revision. `uid` becomes BINARY(10) on the object
  tables that belong to a sector (star_systems, stars, planets, moons,
  asteroid_belts, comets, every phenomenon table, facilities) and is
  UNIQUE on `uid` alone (drops UNIQUE (star_system_id, uid)). A new
  `id_counters` table (per sector for run-time serials, per system for
  body numbers) whose counters only grow and survive the deletion of the
  object or the sector, with an upsert-and-LAST_INSERT_ID allocator like
  `_reserve_id_block`. The row `id` stays as the foreign key. Existing
  rows are migrated (rank = row id order within each sector or system,
  which is how the old ranks were counted) unless the combined reseed
  has not run yet, in which case the change rides it; nebulae need their
  centroid recomputed (see the nebula item). Boss resets by hand anyway,
  so a fresh galaxy is acceptable. Takes the next free Alembic revision.
  Prerequisite: GEN.170.
  Design: [docs/design/object-id-options.md](design/object-id-options.md)

- [ ] **DB.21 A deep pass for the database check: validate every star system, with the estimated time shown first**
  Boss (2026-10-09 23:32Z): "Yes, add an option for a deep pass but warn
  the user the estimated time it will take." Foundations lane 1 left
  `validation.check_star_system` for every system out of `planetgen
  check-db` (DB.8, PR #893) because of its cost. Done: `planetgen
  check-db --deep` also runs `validation.check_star_system` on every
  star system (and the slower table checks DB.8's research proposed,
  `CHECK TABLE ... EXTENDED` and `CHECKSUM TABLE`, where the lane judges
  them worth it), and the "Check the database" section of the Generate
  page offers a "Deep check" option. Before it runs, both show the
  estimated time, from the sector count and the recorded generation statistics (`generation_stats`, PERF.32): a stats kind for the per-system validation is recorded by the check itself, so the second deep check on a version predicts from the first; a conservative fallback is used only when no history exists yet, and then the estimate says so, and the Generate page asks for
  confirmation; the CLI prints the estimate and waits for a yes unless
  `--yes` is given. The deep pass draws its own progress bar (UX.83's
  rule: any sub-step over 15 seconds) and its report lists each failing
  system with its rows, like the plain check. Open question for Boss
  (default `--deep` flag and a "Deep check" option with the estimate and
  a confirmation, as written): other?
  Prerequisite: PERF.33. Related: DB.9, PERF.33, PERF.32, PERF.50.
  From UX.84 (PR #930): the deep check's passes register on the shared
  mechanism: use `steps.Step` with a `STEP_KINDS` entry (a test in
  `tests/test_step_registry.py` fails on an unregistered kind).

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
  Research (2026-10-09, api-design-standards.md): the upload unit is
  `StarSystem.to_dict()` objects in a versioned envelope, not table rows
  (rows named after `schema.sql` would force an API bump on nearly every
  migration). Measured 58.6 KB raw and 13.8 KB gzip per system.

  - [ ] **API.9 Key scopes**
    `admin_api_keys` has no scope column today (every key is an admin
    key). Done: a scope column (read, admin, upload), a control schema
    migration; API.6's user-level keys and the upload right use it.
    Research (2026-10-09, api-design-standards.md): replace "a scope
    column (read, admin, upload)" with four scopes (`read`, `generate`,
    `upload`, `admin`; `admin` implies all, `upload` and `generate`
    imply `read`) in a join table `admin_api_key_scopes`, plus
    `expires_at` and `key_prefix`; control migration v13 or later (control v11 is OPS.13's, v12 is PERF.32's; the text says "v8"); existing keys grandfathered as
    `admin`; throttle `last_used_at` writes; per-key rate-limit buckets
    (the 50/hour IP default must not apply to keys). Add the sweep
    checks of design doc 5.4 to the Done text. The scope work is a
    prerequisite for API.15's `key_id` field and for API.18.

  - [ ] **API.10 Reservations: claimed sectors and id blocks per run**
    `id_blocks` exists (one next id per table, used by parallel
    generation) but nothing reserves sectors. Done: a run record
    holding its claimed sectors and id ranges until it finishes or an
    admin clears it (ADM.13); server-side runs and other uploads skip
    claimed sectors.
    Research (2026-10-09, api-design-standards.md): reservation on the
    existing `_allocate_id` hi/lo mechanism, with a new
    `upload_id_ranges` table and a "reserve more" call; sizes from the
    rows per system in design doc 6.1 plus 25%.

  - [ ] **API.11 Staging tables**
    Done: staging storage (a galaxy schema migration) for received
    batches, flagged incomplete until a whole unit (a system, a sector)
    has arrived and passed API.8's checks, then copied into the real
    tables in one transaction.
    Research (2026-10-09, api-design-standards.md): re-scope from
    "staging storage (a galaxy schema migration)" to upload bookkeeping
    tables in the control database plus a disk spool (`uploads.dir`),
    unless Boss specifically wants rows staged in the galaxy schema
    (about 30 content tables, 65 migrations). Open question for Boss
    (default: spool plus bookkeeping). Verify and finalize of a sealed
    unit run as an RQ job (200 if done within the short wait, else 202
    with a job id; default yes).

  - [ ] **API.12 The download: seed, skeleton and name state**
    Done: a route that returns, compressed, what a local run needs
    (galaxy seed, skeleton, filled sectors, name registries, reserved
    id ranges), so `generate.py` needs no database.
    Prerequisite: API.5.
    Plan (2026-10-07): Downloads the naming key (GEN.70) instead of name
    registries.
    Superseded (2026-10-08): the star and sector name registries are
    still needed (those names stay word-salad, GEN.67); it downloads
    them and the naming key.
    Research (2026-10-09, api-design-standards.md): a paginated, gzip'd
    name-registry download (cursor, 50,000 per page). Whether the seed
    route is public is open (default: it requires an `upload` or
    `generate` key).

  - [ ] **API.13 Generation without a database**
    Today `generate.py` writes through `_db` as it goes. Done: a mode
    where the generators write to an in-memory or on-disk outbox in
    the upload format instead; a client-side cache so an interrupted
    run resumes and resends only what the server hasn't confirmed.
    Research (2026-10-09, api-design-standards.md): same unit as API.3:
    `StarSystem.to_dict()` objects in a versioned envelope.

  - [ ] **API.14 Upload routes, compressed, in batches**
    Done: upload routes for gzip batches, with their own body limit
    (separate from the 2 MB `MAX_CONTENT_LENGTH` of PR #293, from
    API.7's plan), answering what was received and verified.
    Research (2026-10-09, api-design-standards.md): Done gains "per-path
    body limits in every deployment example (nginx, Caddy, IIS, Apache,
    macOS)" and "`Content-Length` required, `Content-Encoding: gzip`
    only, decompressed through a bounded `zlib.decompressobj`".
    Sub-item: the decompression guard with tests (bomb, truncated
    stream, trailing garbage, wrong encoding) and an endpoint-aware
    `_reject_oversized_body`.

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
  Research (2026-10-09, api-design-standards.md): answer the open
  question: the version is bumped by hand against a CI rule. A committed
  OpenAPI file for the remote contract (`docs/api/openapi.json`, built
  from the Pydantic models by `planetgen/api/spec.py`), `python -m
  planetgen.api.spec --check` in tests, and `oasdiff breaking` and
  `changelog` between PR base and head (a breaking diff means a higher
  MAJOR, any other diff a higher MINOR, none unchanged). No
  `changes/<name>.api.md` note kind. The compatibility contract covers
  the remote-client routes only. The table is generated from
  `planetgen/api/version.py` history, and the stamp job fills in the
  release number.

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
  Research (2026-10-09, api-design-standards.md): answer "how many older
  versions": the current major, plus the previous major for one release
  cycle after a major bump, for upload payloads only (no read
  converters). Add `GET /api/version` (public, database-free, returns
  `api_version`, `min_client`, `max_client`, `release`, `version_key`,
  `galaxy_schema`) and the `PlanetGen-API-Version` header both ways;
  `generate.py --remote` checks before reading config. A refusal is `400
  api_version_unsupported` (or `410` past a sunset), not `426`; the
  check does not depend on `/api/health`. Add the version-key comparison
  (design doc 4.5). Open questions for Boss (defaults taken): a remote
  client may differ in release within the accepted range, with a
  warning; a version-key mismatch (OS, architecture, Python) warns and
  does not refuse.

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
  Research (2026-10-09, settings-seo-and-accounts.md): effective rights
  are the intersection of the key's scope and the owner's current role,
  evaluated per request (open question for Boss, default: a user-level
  key made by an admin acts as the creating account with read-only
  scope).

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
  Research (2026-10-09, api-design-standards.md): log through the
  activity log first (a new `API` category in `CATEGORIES` and
  `docs/config.md`), with `X-Request-Id`, an `after_request` hook,
  `internal=1` for in-process page calls (about 6 per page view) and an
  `api_log.internal` setting, and the privacy and retention paragraph
  (open question for Boss, default 30 days full IP, then truncated, 13
  months total). An `api_calls` table with a background writer only when
  an admin page needs queries. "Console" in `god` means API-shaped
  actions run by CLI tools without a request; existing `GEN` lines
  unchanged.

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
  Research (2026-10-09, api-design-standards.md): adopt the numbers in
  design doc 6.6 for Boss's agreement; API.7 can then close. Open
  question for Boss (default as written).

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
  Research (2026-10-09, api-design-standards.md): adopt the auto-correct
  versus reject rule (design doc 7.4) as the default answer: correct
  only a clashing name and pure derived columns; reject everything else;
  the server returns a rename map; a rejection does not release the
  claim.

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
  Research (2026-10-09, api-design-standards.md): record the 400 versus
  422 split exactly (the envelope layer gives 400, the kind model and
  physics checks give 422), the RFC 9457 problem body keeping the old
  `error` and `errors` members, the size caps (256 KiB, depth 12, 500
  slots, 100 moons per planet, 5 queued jobs per key) and the `seed` and
  `seed_used` rule. Existing routes keep 400 (open question for Boss,
  default: only recipes and uploads move to 422). The
  `require_json_body` depth fix is a prerequisite.
  Prerequisite: API.9.

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
  Research (2026-10-09, api-design-standards.md): recipes stay one
  region per request.
  Prerequisite: API.18.

- [ ] **API.22 An API version number: one sequential integer, shown in admin and in the status response**
  Boss (2026-10-09 20:59Z): "I want ... API version number (same) by the
  end of phase 1", a plain sequential integer like the DB schema number.
  Done: `API_VERSION = 1` in one module, bumped by a PR that makes a
  breaking change to an endpoint (a removed or renamed route or field, a
  changed meaning or type, a new required parameter); additive changes
  do not bump it. It is returned by the API status response and every
  `/api` response header, shown on the admin status page in place of the
  release string now labelled "API version", and recorded in
  docs/api.md's change list. A test fails when the route table or
  response shapes change without the integer moving (a stored schema
  snapshot). API.4's compatibility data and the remote-run handshake
  (API.17) compare this integer.
  Decided (Boss, 2026-10-09 21:02Z): bump on any breaking change to an
  endpoint; additive changes do not bump.
  Prerequisites: none. Related: API.4, API.17.

- [ ] **API.23 The object ID as the public reference: pages, URLs, the API, wiki links and objectref use it in place of row ids**
  Source: docs/design/object-id-options.md section 0 (Boss decided
  2026-10-09 22:39Z: birth location plus serial, galaxy-wide; the same
  length for every object; always identifies that one object; up to 128
  bits but shorter preferred; fix the deficits; no backward
  compatibility). Filed from the object-ID research thread. Nothing is
  built until Boss asks.
  Done: pages, URLs, API routes and payloads, wiki links and `objectref`
  use the object ID in place of row ids (lookup probes the object tables
  by `uid`; a kind prefix such as planet:ID is only a hint). Row ids
  stay internal. No compatibility shim. This is a breaking API change,
  so it bumps API.22's API version number.
  Decided (Boss, 2026-10-10 02:48Z, via Foundations lane 1): yes, the 80-bit object ID replaces row ids in pages, URLs and the API, and Boss accepts the API break. Cleared to build once API.22, GEN.171 and GEN.172 are in.
  Prerequisites: API.22, GEN.171, GEN.172.
  Design: [docs/design/object-id-options.md](design/object-id-options.md)

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
  Research (2026-10-09, api-design-standards.md): the stale flag
  defaults to 3 days without contact and only flags; the page reads
  `upload_runs`, `upload_units` and `upload_claims`.

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
  Plan (2026-10-07): The worker count now means RQ worker processes
  (PERF.24).
  Research (2026-10-09, performance-eta-queue-and-caching.md): keep
  lowered worker priority, but it does not protect the database: add one
  connection per worker, a cap on total worker connections and the
  reserved interactive worker; re-run PERF.34's page-time test under
  Ludicrous Speed.

- [ ] **ADM.28 A simpler Generate page: layer specs, a Customize window and plain controls**
Boss (2026-10-07 11:47Z): "Actual specs on layers on the generation
  screen.  Let's also have a customize button that brings up a special
  window with all the settings, the generate screen is getting a bit
  complex, we need to do a full rework to make it easier to use." Done:
  the page shows the galaxy's layer specs (count, height, extent, how
  many charted), the common actions stay on the page and every other
  setting moves into a Customize dialog, and its subitems are done.
  GitHub issue [#736](https://github.com/dwhagar/planetGen/issues/736) (Boss, 2026-10-09 02:26Z): "Each set of settings should be a tab for the generate screen so the user only sees the ones relevant to what they are looking at." So the Customize dialog groups its settings into tabs, one per kind of generation.

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
  Prerequisite: MAP.120.

- [ ] **ADM.36 Change an object's trajectory vector**
  Boss (2026-10-07 11:47Z): "Give the user the ability to change an
  object's trajectory vector." Done: an admin sets a body's velocity
  vector (in any frame the position object supports) from its detail
  screen; the change is validated (escape, collisions), logged, and used
  by the next orbital update.
  Research (2026-10-09, collisions-and-mergers.md): after an admin edits
  a velocity vector, recompute the osculating elements,
  `next_update_due` and Hill and cover-list membership at once (the next
  run could be far off); report if the new vector makes the orbit
  unbound. The validation step ("escape, collisions") computes the swept
  path against stored bodies and shows the predicted event; admin-aimed
  hits are the realistic way to see the collision code run.
  Prerequisite: GEN.109.

- [ ] **ADM.43 A full configuration page under Admin**
  Boss (GitHub issue [#515](https://github.com/dwhagar/planetGen/issues/515), 2026-10-08 00:39Z): "A page under admin is
  needed to configure everything EXCEPT for the database information.
  All other settings should be accessed and changeable from there."
  Done: an Admin page generated from ADM.42's model lists every option
  with its help text, edits and validates it, saves config.json, and
  says which changes need a restart.
  Research (2026-10-09, settings-seo-and-accounts.md): re-scope to three
  tiers (admin, Owner, read-only host-level). `jobs.python`, `jobs.dir`,
  log paths, `secret_key`, `proxy_fix`, `admin_cookie_insecure`,
  `redis.url` and `api_base_url` are read-only in the browser: they are
  code-execution or request-forgery paths. Add: a diff before save,
  timestamped backups with restore, audit rows, per-field errors, a
  pending-restart banner from a startup snapshot, import and export, a
  recent-password prompt for the Owner tier. `set-permissions.sh` and
  the installers must create the settings folder for the web user. Open
  question for Boss (default): show everything, edit in the browser only
  for the admin and Owner tiers.
  Read side built with ADM.42 (PR #912): the `settings.json` overlay
  is read, nothing writes it yet; options marked `x-editable` false are
  file-only; `x-min_role` owner marks Owner-tier options;
  `util/appconfig.py` is gone (log paths are in `util/logpaths.py`);
  `docs/config.md` has a generated table (`python -m planetgen.cli.config
  docs`).

- [ ] **ADM.44 Web, Open Graph and SEO settings**
  Boss (GitHub issue [#743](https://github.com/dwhagar/planetGen/issues/743), 2026-10-09 02:55Z): "Add config options to
  set up icons, open graph fields for discord and other link previews,
  SEO fields for description and keywords, and web documents for things
  like robots.txt and such." Done: the settings page gains site icons,
  Open Graph and link-preview fields, a description and keywords, and
  editable robots.txt and similar documents, all served by the site.
  Research (2026-10-09, settings-seo-and-accounts.md): replace the
  "description and keywords" wording with the settings table of A8
  (Google ignores meta keywords; keep the field). `web.base_url`
  validation is a prerequisite; add the indexing matrix, a generated
  robots.txt with extra and override, a sitemap index capped at 50,000
  URLs per file, security.txt and humans.txt, the Open Graph and Twitter
  tag block, upload of icons and the default OG image, and the pluggable
  card renderer (`render_card(kind, object_id) -> bytes`,
  `og.generated_cards` off by default). It does not depend on VIEW.3.
  Open question for Boss (default): generated system and sector pages
  are `noindex,follow`, top pages only in the sitemap, switchable with
  `seo.detail_pages`.

- [ ] **ADM.45 Prevalence fields take the override share directly and must total 100%**
  Boss (2026-10-09 07:48Z): "Prevalence fields instead of being +/- % they will be just type in the override % number, the form should force the user to make sure the whole thing =100% and should make it clear what to do so it isn't confusing. I don't think we need a density dependent share."
  Follows ADM.37 (PR #514), which made the Generate page show each
  feature's real default share. Done: every prevalence field on the
  Generate page is a plain number box holding the share, in percent,
  that the feature should have, starting at its real default; the user
  types the override instead of a plus or minus change. The page shows a
  running total of the shares that belong together, says in words what
  to do when it is not 100% (which fields to raise or lower and by how
  much), and will not start the run until the set adds up to exactly
  100%. A test submits sets that add up and sets that do not, and the
  CLI's `--prevalence` accepts the same shares and rejects a set that is
  not 100%. A share that depends on the local star density is not part
  of this item; the shares stay the same in every sector.

## SEC: Security

The login protection of 2026-10-01 (SEC.1, SEC.20 to
SEC.28: the always-on activity log, per-address and per-username
lockouts, trusted devices, the password blocklist and hashing cost,
two-factor sign-in and the fail2ban example) shipped in PRs #217, #220
and #221.

Design: [docs/design/login-brute-force-protection.md](design/login-brute-force-protection.md)

- [ ] **SEC.32 Argon2id password hashing, a 1,024-character password limit and the `__Host-` cookie prefix**
  Argon2id with the OWASP minimum (19 MiB, t=2, p=1) via `argon2-cffi`;
  a prefix-dispatching `verify_password`; `needs_rehash` true for old
  schemes and weaker Argon2 parameters (upgraded at login); a maximum
  password length of 1,024 (none exists in `admin/auth.py`
  `validate_password_policy`; only the 2 MB body limit bounds it); cost
  parameters in the settings model; the `__Host-` session cookie prefix
  when Secure. Ideally before USR.2 so new accounts start on it. Open
  question for Boss (default yes): move from PBKDF2-600k to Argon2id (a
  new dependency); minimum length 12 with the blocklist, 15 only if
  NIST's password-only guidance is wanted.
  Prerequisites: none.
  Design: [docs/design/settings-seo-and-accounts.md](design/settings-seo-and-accounts.md)

## TEST: The test suite

From the test suite plan of 2026-10-01 (Boss: "build me a list of TEST
items (TEST.x) to build our testing suite to cover all edge cases, and
prepare to do another solid full on bug hunt. Don't implement just plan
for it"). Each item ends with its areas in brackets. "Suspected bug"
means read from the code, not reproduced yet; the bug hunt confirms or
clears each one.

### Infrastructure and CI

- [ ] **TEST.110 Object ID tests: identical IDs on 1 and 4 workers, none reused, none missing**
  Source: docs/design/object-id-options.md section 0 (Boss decided
  2026-10-09 22:39Z: birth location plus serial, galaxy-wide; the same
  length for every object; always identifies that one object; up to 128
  bits but shorter preferred; fix the deficits; no backward
  compatibility). Filed from the object-ID research thread. Nothing is
  built until Boss asks.
  Done: tests that 1 worker and 4 workers give identical IDs; no ID is
  reused after a delete; an ejected planet keeps its ID; every object in
  a saved sector has a 20-digit ID; the golden fill digests (TEST.77)
  use the new IDs.
  Prerequisites: GEN.171, GEN.172, GEN.176.
  Design: [docs/design/object-id-options.md](design/object-id-options.md)

- [ ] **TEST.111 test_ensure_sector_generated_creates_then_reuses_the_same_sector fails in a busy parallel run (bug)**
  Reported by Bugfixes lane 1 (2026-10-09 23:32Z):
  `test_ensure_sector_generated_creates_then_reuses_the_same_sector`
  failed once in a busy parallel full-suite run and passed alone, with
  no change to the code it covers. Done: the cause is found (a shared
  sector address or timing under load is the first thing to check) and
  the test is made robust without skipping or loosening it, or the
  product bug it hides is fixed. Related: TEST.71, TEST.73, OPS.19.
  Note (2026-10-09): Bugfixes lane 1 (PR #955, 2026-10-10): no repro in
  15 loaded runs; the test's asserts now print the results, so the next
  failure names the cause. Stays open until it recurs and is fixed, or
  Boss closes it.
  Prerequisites: none.

- [ ] **TEST.115 Three tests fail on the Windows CI leg in every recent run: Redis in WSL is unreachable (bug)**
  Reported by Bugfixes lane 1 (2026-10-10 03:57Z, PRs #957 and #958):
  test_a_slow_runner_still_alive_is_starting_not_interrupted,
  test_without_redis_no_job_starts and
  test_without_redis_windows_runs_the_job_itself fail in all 8 recent
  main runs on the Windows leg, because Redis in WSL is not reachable
  from the Windows side (127.0.0.1:6379 refused). Done: the Windows leg
  runs these tests with a reachable Redis or skips the ones that need
  none with a stated reason, and the tests that test the missing-Redis
  path set up that state themselves; the Windows leg passes on main. CI
  now runs only by hand (Actions, CI, Run workflow), so the leg is
  checked when someone runs it.
  Superseded (2026-10-09): Superseded by OPS.39 (Boss, 2026-10-10
  04:25Z): Windows support is being removed. Do not start; retire with
  OPS.39.
  Prerequisites: none. Related: TEST.111, OPS.19, PERF.24.

- [ ] **TEST.116 test_bughunt_end_to_end stores M2V or M6V where it expects K2V on the MySQL 8.4 and MariaDB legs once (bug)**
  Reported by Bugfixes lane 1 (2026-10-10 03:57Z): the K2V check in
  test_bughunt_end_to_end failed once on MySQL 8.4 and once on MariaDB
  in run 37882710453 (the stored star was M2V or M6V instead of K2V); it
  did not reproduce on current main in 83 local runs. Done: the cause is
  found (a seed or ordering dependence in the test, or a real bug in
  star-type selection) and the test is made robust without loosening it,
  or the product bug is fixed. Open question for Boss (default: leave
  open until it recurs, then investigate with the failing run's data).
  Prerequisites: none. Related: TEST.111, TEST.71.

- [ ] **TEST.118 test_star_scatter_passes.py fails twice on main since GEN.184 raised the luminosity floor to 2500 or more (bug)**
  Reported by Foundations lane 2 (2026-10-10 04:23Z, GEN.183 merge, PR
  #969): test_sampled_stars_stay_inside_their_mass_range and
  test_the_scatter_runs_the_mass_pass_then_a_luminosity_pass_that_skips_marked_sectors
  fail on main since GEN.184's luminosity floors (2500 L_sun at the
  lowest) met GEN.185's test setup. Done: the tests' setup uses a floor
  the new ladder allows, and the full test file passes on main.
  Prerequisites: none. Related: GEN.184, GEN.185, GEN.183.

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
  Research (2026-10-09, settings-seo-and-accounts.md): order stays; add
  one shared table `account_tokens` used by USR.4, USR.5 and USR.6
  (itsdangerous is not used for these tokens). The text names stale
  paths (`stellarObjects/control_schema.sql`,
  `stellarObjects/adminAuth.py`, `src/html/api/loginbackoff.py`); the
  files are now `src/planetgen/db/control_schema.sql`,
  `src/planetgen/admin/auth.py` and `src/planetgen/admin/throttle.py`
  with `api/loginguard.py`.

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
    Research (2026-10-09, settings-seo-and-accounts.md): answer the open
    questions as in B7: one table with a role column (renamed to
    `accounts` in a follow-up; open question for Boss, default: add
    columns now); the lowest id becomes Owner; a database-level
    single-Owner unique index on a generated column; admins cannot
    demote themselves; admins disable users, only the Owner disables
    admins; delete by self or Owner; admins create API keys only at
    first. Add the console CLI `planetgen.cli.accounts
    --set-owner/--reset-password`. The control schema is the next
    Alembic revision after the current one.

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
    Research (2026-10-09, settings-seo-and-accounts.md): answer as in
    B9: stdlib smtplib, Owner-only (SMTP, `base_url`, From address; open
    question for Boss, default yes), password in the web-owned 0600
    file, no SMTP means invites show the link only; sent from an RQ job.

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
    Research (2026-10-09, settings-seo-and-accounts.md): add the
    `invite_uses` table, atomic redeem, an open-invite cap of 50 and
    rate limits (B5).

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
    Research (2026-10-09, settings-seo-and-accounts.md): lifetimes:
    reset 30 minutes, first-password 24 hours, email change 1 hour; a
    notification email on password and email change; an Owner reset asks
    for the TOTP code when enrolled.

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
    Research (2026-10-09, settings-seo-and-accounts.md): the state
    machine of B10: the old Owner becomes admin; the target must be an
    admin with TOTP; it expires after 48 hours; the Owner can cancel;
    the confirm link lands on a neutral page because of
    `SameSite=Strict`.

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
    Research (2026-10-09, settings-seo-and-accounts.md): nothing new
    beyond B11 (forced sign-in would force `seo.indexing=off`).

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
    Research (2026-10-09, settings-seo-and-accounts.md): change "the
    count kept in the control database" to Redis through Flask-Limiter
    (`moving-window`, key `user:<id>`), add a concurrency cap, and make
    both numbers settings (`limits.one_off_per_hour` 30,
    `limits.one_off_concurrent` 3). The text names
    `html/web/system_page.py`; the file is
    `src/planetgen/web/system_page.py`. Open question for Boss
    (default): user sessions last 30 days while admins stay at 12 hours.
    Prerequisite: USR.2.

- [ ] **USR.9 `seo.privacy_note` text, "download my data" and "delete my account" on the account page**
  A `seo.privacy_note` text setting shown on the account page, plus
  "download my data" and "delete my account" (B11). Not asked for by
  Boss; drop if unwanted.
  Prerequisites: none.
  Design: [docs/design/settings-seo-and-accounts.md](design/settings-seo-and-accounts.md)

## OPS: Installers, hosting, CI, releases

OPS.1 shipped with the version scheme in `changes/README.md`.

- [ ] **OPS.28 Generator epoch and battery digest: say whether two checkouts generate the same galaxy**
  An integer `generator_epoch` bumped only by a PR that changes
  generated output, kept in `src/tests/golden/epochs.json` (epoch, first
  release of that epoch). `bump_version.py --check-pr` enforces all or
  none: golden changed, epoch bumped, and a `generation-output: changed`
  line in the `changes/` note. A `golden_update` tool bumps the epoch
  and refuses under `CI`. `generator_battery_digest()` is a pure
  database-free function over 6 to 8 tiny canned cases (about 3 to 5 s).
  The OPS.13 history row also holds the epoch, the battery digest, the
  lock hash and an environment JSON (libc, numpy, astropy, scipy,
  scikit-image versions), the last 10 per galaxy as already specified;
  the seed itself is not changed. ADM.18's settings file records the
  epoch and `fp_spec`, and is written with LF only. Open questions for
  Boss (defaults taken): (1) is the same galaxy on Linux, macOS and
  Windows a goal, exact for integers, strings and structure and equal to
  9 significant digits for floats? Default yes; if not, drop the Windows
  and macOS CI additions and keep the epoch. (2) Is an integer epoch in
  the version record acceptable? Default yes. (3) Does a mismatch refuse
  a remote or rebuild run? Default: refuse on an epoch mismatch, warn on
  a battery-digest mismatch naming the first differing case.
  Note (2026-10-09): The Windows parts are dropped by OPS.39 (Boss,
  2026-10-10 04:25Z).
  Prerequisite: GEN.135.
  Design: [docs/design/reproducible-galaxies.md](design/reproducible-galaxies.md)

- [ ] **OPS.29 Update reload: cover the gunicorn units and non-Apache hosts, and say what a reload aborts**
  Follow-up to OPS.8 (built in PR #808): check it against the research.
  Skip if Apache is absent or stopped (never start it), `configtest`
  first, `systemctl reload apache2` (or `httpd`), restart only if the
  script itself installed mod_wsgi and no Generate job runs, check
  Apache is up and `/api/health` answers, exit non-zero on failure.
  Cover the gunicorn units (`planetgen-gunicorn` on Linux, `launchctl
  kill SIGHUP` on macOS) in the same function, since nginx and Caddy
  hosts have no Apache. State plainly in the docs that a reload aborts
  in-flight daemon requests after about 4 s while `touch
  src/html/wsgi.py` does not. Update `docs/deployment/apache.md` (lines
  22, 29, 164), the `docs/server-checklist.md` "reload the app" row and
  the closing text of `update.sh` and `install.sh`; `update.sh` prints
  "restart" after `a2enmod` though a graceful reload loads new modules.
  Open question for Boss (default: print the manual command and end
  non-zero on a failed reload).
  Prerequisites: none.
  Design: [docs/design/ops-scheduling-and-rotation.md](design/ops-scheduling-and-rotation.md)

- [ ] **OPS.30 A lock helper for the maintenance run**
  `planetgen/util/locks.py` with `filelock` (pinned `==3.19.1` for
  Python 3.9, current otherwise), a sidecar info file, `timeout=0`
  helpers and `PermissionError` handling, plus a two-process test that
  runs on the Windows CI leg.
  Note (2026-10-09): The Windows parts are dropped by OPS.39 (Boss,
  2026-10-10 04:25Z).
  Prerequisites: none.
  Design: [docs/design/ops-scheduling-and-rotation.md](design/ops-scheduling-and-rotation.md)

- [ ] **OPS.31 Lint every example plist, XML and service file in CI**
  `plistlib` for `examples/**/*.plist`, an XML parser for `*.xml` and
  `systemd-analyze verify` where available, so a malformed example
  cannot ship again.
  Prerequisites: none.
  Design: [docs/design/ops-scheduling-and-rotation.md](design/ops-scheduling-and-rotation.md)

- [ ] **OPS.34 Windows Redis in WSL: fix the keep-alive advice and add a Start-RedisInWsl remedy**
  `docs/deployment/windows.md` (Redis, step 4) says a logon task running
  `wsl -d Ubuntu` keeps Redis alive; it does not (a WSL instance idles
  out after about 15 s and systemd services do not hold it). Document a
  hidden keep-alive process (`wsl -e sleep infinity`) started at logon,
  marked unverified until tested on Windows; add a `Start-RedisInWsl`
  remedy to `Test-Redis` (only when the account can see a distro); state
  that unattended Windows servers cannot rely on WSL. Open question for
  Boss (default: keep telling unattended servers to use a Linux VM).
  Superseded (2026-10-09): Superseded by OPS.39 (Boss, 2026-10-10
  04:25Z): Windows support is being removed. Do not start; retire with
  OPS.39.
  Prerequisites: none.
  Design: [docs/design/ops-scheduling-and-rotation.md](design/ops-scheduling-and-rotation.md)

- [ ] **OPS.35 A vendored-version lock file for the static libraries**
  `static/vendor/VENDORED.json` (version, npm `dist.integrity`, SHA-256
  per shipped file, esbuild version, licence) written by the vendor
  scripts and checked by a pytest; add the esbuild version to
  `THIRD_PARTY_NOTICES.txt` and list the transitive licences bundled in
  Shoelace (Lit, floating-ui, tinycolor).
  Prerequisites: none.
  Design: [docs/design/map-ui-and-frontend-libraries.md](design/map-ui-and-frontend-libraries.md)

- [ ] **OPS.37 A Generator version number: one sequential integer, shown in admin and the API**
  Boss (2026-10-09 20:59Z): "I want a Generator version number (like DB
  number just sequential integers) ... by the end of phase 1." Done: one
  plain integer, 1, 2, 3 ..., named the generator version, that says
  which rules made a galaxy. It is OPS.28's `generator_epoch` under the
  name Boss asked for, so there is one counter, not two: bumped by the
  same PR rule (a change that alters generated output for the same seed
  bumps it), stored with each galaxy's history row (OPS.13), per sector
  (DB.16) and in the settings file (ADM.18), shown on the admin status
  page beside the DB schema number, and returned by the API status
  response. The release version (MAJOR.REVISION.BUILD) stays as it is.
  Decided (Boss, 2026-10-09 21:02Z): the generator version and OPS.28's
  `generator_epoch` are the same number, bumped only when output changes
  for the same seed.
  Prerequisite: OPS.28.

- [ ] **OPS.39 Remove Windows support; keep only a simple docs/WINDOWS.md**
  Decided (Boss, 2026-10-10 04:25Z, in the Foundations lane 1 thread):
  rip out all Windows support. If someone wants it to work on Windows
  they do that work themselves; the most the project provides is a
  simple docs/WINDOWS.md with basic instructions for a typical Windows
  setup. Done: the Windows CI leg, update.ps1, install.ps1, every other
  .ps1 script and scheduled-task helper, Windows branches in the code,
  Windows-only tests and test setup, and Windows mentions in the docs
  are removed; docs/reference/deployment/windows.md (and any
  docs/deployment/windows.md) is replaced by docs/WINDOWS.md; the README
  and install docs name Linux (and macOS where it still applies) only; a
  search for 'windows', 'ps1', 'WSL' and 'win32' finds nothing outside
  docs/WINDOWS.md and history. This supersedes TEST.115 (Windows-leg
  Redis failures) and OPS.34 (Windows Redis in WSL): retire both with
  this item, and drop the Windows halves of OPS.16, OPS.17, OPS.28 and
  OPS.30 (the epoch question about Windows reproducibility is answered
  by this). No backward compatibility.
  Prerequisites: none. Related: TEST.115, OPS.34, OPS.16, OPS.17,
  OPS.28, OPS.30.
  Lane (Boss, 2026-10-10 04:26Z): Foundations lane 3 (not lane 1).

## DOC: Documentation

- [ ] **DOC.4 Correct the stale statements the research found in docs, docstrings and comments**
  The research (PR #812 and its handoff) listed these as out of date;
  none is a behaviour change. Docs: `docs/cli.md` (bright-star default
  is 1,000 L_sun, not 500, and the default shape drew 26.9 million
  stars, not about 60 million); `docs/html-interface.md` line 94
  (`numberformat.js`: scientific from 7 whole digits with no decimals, 5
  with decimals, UX.36); `docs/database-schema.md` (Versioning and the
  `migrate_database` bullet still describe `_migrate_vN_to_vN+1`,
  `_MigrationConnection` and `_migrate_v52_to_v53`; its "Planned
  changes" DB.9 row; old module names on lines 105, 1176, 1182:
  `systemData.py`, `doubleStar.py`, `wideBinary.py`, `starData`,
  `compactRemnant`); `docs/deployment/macos.md` ("Monthly orbit update";
  an asleep Mac runs the job at wake, only a powered-off Mac skips);
  `docs/design/architecture.md` ("Meant for a monthly timer" at lines
  125, 268 and 808; lines 138 to 141 and 472 to 473 name old modules;
  the shared-infrastructure and static-file tables use pre-restructure
  paths; check `architecture.html` too); `docs/api.md` "Not done yet"
  (API.11 is no longer staging tables, API.9 no longer read, admin or
  upload); the BVH wording left over from MAP.102 in
  `docs/design/todo-number-map.md`, `design-decisions.md` and
  `galaxy-drilldown-navigation.md` section 16 (use "screen-space
  typed-array pick"); `docs/design/orbital-updates.md` sections 1, 4, 5,
  10.1, 10.3 to 10.6 and 11 (nearest-10, volume-sum radius, 2.44 for all
  bodies, 2,800 / 200 / 15 kpc, the log halo, the macro step);
  `reproducible-galaxies.md` and `api-design-standards.md` (no numpy:
  scipy and astropy are used); `generation-determinism.md` 4.3 (DB.9's
  checksum is not the GEN.58 leaf); `sky-view.md` 3.2 (the backfill cost
  per cell is now measured); `galaxy-coordinate-system.md` header and
  `geometry.py` docstring (`skeleton.py`, `density.py`, `geometry.py`;
  ring 0 into ring 1 gives three overlapping slots);
  `docs/plan/notes.md` (`bright_star_blocks` was replaced by
  `sector_stats`, schema v53); `bright-star-timing/report.md` (the "26B
  (x966)" column matches nothing); `docs/design/Library Migration
  Workflow.md` section 3 (hand-written revisions and `schema_v61`
  onward). Code text: `db/alembic_runner.py` module docstring (no
  `_migrate_vN_to_vM` steps remain), `cli/orbits.py` header ("once a
  month or so"), the monthly comments in
  `examples/maintenance/planetgen-{update,orbits@}.timer`,
  `generation/phenomena/compact_remnant.py` line 9 (cites a
  `docs/design/exotic-phenomena.md` that does not exist; point it at
  `anomalies.md` and `multistar-and-compact-systems.md`),
  `tuning.CIVILIZATION_CHANCE` (mean spacing is 134 ly, nearest
  neighbour 74 ly, not 150 to 200), `population/model.py` (the pass runs
  only with `--population`), the `population --rescan` help (ages and
  traits come out the same), `web/__init__.py` line 157 and
  `static/galaxystages.js` line 720 (scientific notation thresholds),
  `facilityform.js` `aria-valuetext`. TODO.md text still naming old
  paths: GEN.55 and GEN.40 (`generate.py`), GEN.105
  (`random.seed()` in `api/edits.py`, `updateOrbits.py`), UX.22 and
  UX.23, USR.1, USR.2 and USR.8. Cross-links to add:
  collisions-and-mergers.md from orbital-updates.md section 5,
  nebula-and-asteroid-field-classes.md and interstellar-object-rates.md;
  multistar-and-compact-systems.md from orbital-updates.md (section 4),
  anomalies.md, interstellar-object-rates.md (stars per system 1.38) and
  database-schema.md; api-design-standards.md from `docs/api.md` and
  `reproducible-galaxies.md` section 10; course-avoidance.md from
  course-routing.md; galactic-potential.md from orbital-updates.md 10.1;
  the three map-ui design docs. Open the Boss-facing decisions these
  touch as they are made.
  Prerequisites: none.

- [ ] **DOC.5 Rewrite the object ID docs: object-ids.md, database-schema.md and api.md**
  Source: docs/design/object-id-options.md section 0 (Boss decided
  2026-10-09 22:39Z: birth location plus serial, galaxy-wide; the same
  length for every object; always identifies that one object; up to 128
  bits but shorter preferred; fix the deficits; no backward
  compatibility). Filed from the object-ID research thread. Nothing is
  built until Boss asks.
  Done: docs/design/object-ids.md (the GEN.68 and GEN.69 section) is
  rewritten for the 80-bit ID; GEN.69's hash scheme is marked
  superseded; database-schema.md and api.md describe the new column and
  the ID as the public reference; GEN.72 and GEN.73 and the position-ID
  naming (GEN.64 stays as the name of interstellar objects) are updated
  to match.
  Prerequisite: GEN.170.
  Design: [docs/design/object-id-options.md](design/object-id-options.md)

- [ ] **DOC.6 A static help section in the web interface: page template, index, per-page help links and a coverage test**
  Boss (2026-10-09 23:53Z): "Add to-do items to build static
  documentation pages for all features accessible through the web
  interface." Done: the web interface serves a Help section of static
  pages (no database, no JavaScript needed to read them), with one
  shared template that matches the site shell, an index page that groups
  the pages below by feature, a "?" link in the header of every
  interface page that opens the matching help page, a search box over
  the help text, and a plain page of the current version. The pages are
  written as Markdown files in the repository (docs/help/) and built
  into HTML by the update script, so they can be edited without touching
  code. A test lists every user-facing route and fails when a route has
  no help page or a help page has no route, and a second test checks
  that every link inside the help section resolves. Open question for
  Boss (default: Phase 2, after the features they describe have settled;
  Markdown source in docs/help/, built at update time, served at /help):
  other?
  Prerequisites: none. Related: DOC.7, DOC.8, DOC.9, DOC.10, DOC.11,
  DOC.12, DOC.13, DOC.14, DOC.15, DOC.16.

- [ ] **DOC.7 The Galaxy Map help page: layers, zoom, fly-through, Color by, select modes, bookmarks and the locate box**
  Boss (2026-10-09 23:53Z): "Add to-do items to build static
  documentation pages for all features accessible through the web
  interface." Done: every control and gesture on /galaxy: the opening
  view, the layer and ring menus, zooming and the fly-through camera
  (MAP.146 when built), Color by, Select mode, the block, slab and wedge
  menus, nebula and territory overlays, bookmarks, and how big a tile is
  and why a sector may be empty. Screenshots are captured by script so
  they can be regenerated. The page is a Markdown file in docs/help/
  built into the Help section, linked from each page it describes, and
  updated by any later change to those features (a feature item is not
  done until its help text is).
  Prerequisite: DOC.6. Related: DOC.6.

- [ ] **DOC.8 The sector help pages: the sector list, a sector page, the sector map and the sector scene**
  Boss (2026-10-09 23:53Z): "Add to-do items to build static
  documentation pages for all features accessible through the web
  interface." Done: the sectors list (/sectors), the sector page
  (/sector/<id>) with its star-system table, the sector map and 3D
  scene, sector paths and neighbours, what the sector address means, and
  the admin edit panel that appears for admins. The page is a Markdown
  file in docs/help/ built into the Help section, linked from each page
  it describes, and updated by any later change to those features (a
  feature item is not done until its help text is).
  Prerequisite: DOC.6. Related: DOC.6.

- [ ] **DOC.9 The star system help pages: the system list, a system page, the system map and the planets, moons and belts shown**
  Boss (2026-10-09 23:53Z): "Add to-do items to build static
  documentation pages for all features accessible through the web
  interface." Done: the systems list (/systems), the system page
  (/system/<id>), the system map and scene (orbits, scale, time), the
  tables of stars, planets, moons and belts, the habitability index
  (GEN.83/GEN.89) and what each column means with its units. The page is
  a Markdown file in docs/help/ built into the Help section, linked from
  each page it describes, and updated by any later change to those
  features (a feature item is not done until its help text is).
  Prerequisite: DOC.6. Related: DOC.6.

- [ ] **DOC.10 The search and navigation help pages: search, nearby, the nav page and routes**
  Boss (2026-10-09 23:53Z): "Add to-do items to build static
  documentation pages for all features accessible through the web
  interface." Done: search (/search), nearby search (/nearby), the nav
  page (/nav), routes between systems and across sectors,
  within-N-parsec search, travel times (NAV.11 when built), and how
  unknown space and asteroid fields affect a course. The page is a
  Markdown file in docs/help/ built into the Help section, linked from
  each page it describes, and updated by any later change to those
  features (a feature item is not done until its help text is).
  Prerequisite: DOC.6. Related: DOC.6.

- [ ] **DOC.11 The reference browser help pages: species, polities, object classes and phenomena**
  Boss (2026-10-09 23:53Z): "Add to-do items to build static
  documentation pages for all features accessible through the web
  interface." Done: /species, /polities, /classes and /phenomena with
  their detail pages: what each field means, how classes and codes are
  named, and how a phenomenon differs from a star system. The page is a
  Markdown file in docs/help/ built into the Help section, linked from
  each page it describes, and updated by any later change to those
  features (a feature item is not done until its help text is).
  Prerequisite: DOC.6. Related: DOC.6.

- [ ] **DOC.12 The account help pages: signing in, two-factor, the account page, API keys and bookmarks**
  Boss (2026-10-09 23:53Z): "Add to-do items to build static
  documentation pages for all features accessible through the web
  interface." Done: /login, the one-time code step, /account, two-factor
  setup, API keys and their scopes (API.9 when built), bookmarks (UX.47
  when built) and what is stored about a signed-in user. The page is a
  Markdown file in docs/help/ built into the Help section, linked from
  each page it describes, and updated by any later change to those
  features (a feature item is not done until its help text is).
  Prerequisite: DOC.6. Related: DOC.6.

- [ ] **DOC.13 The Generate page help pages: layer specs, spans, radial fills, directives, one-off systems and jobs**
  Boss (2026-10-09 23:53Z): "Add to-do items to build static
  documentation pages for all features accessible through the web
  interface." Done: /admin/generate and its job pages: the layer specs,
  Customize window, spans and radial fills (ADM.28, ADM.29, ADM.30),
  directives (GEN.96), random neighbourhoods (GEN.97), the one-off
  system generator (/admin/generate/system) and its download, reading a
  job log and the progress bars, and what a regeneration changes. The
  page is a Markdown file in docs/help/ built into the Help section,
  linked from each page it describes, and updated by any later change to
  those features (a feature item is not done until its help text is).
  Prerequisite: DOC.6. Related: DOC.6.

- [ ] **DOC.14 The admin help pages: the queue, the stats page, settings, lockouts and the naming key**
  Boss (2026-10-09 23:53Z): "Add to-do items to build static
  documentation pages for all features accessible through the web
  interface." Done: /admin, /admin/queue and its tree, confirm and
  action pages, /admin/stats (galaxy settings, naming key, lockouts),
  the worker count and the performance statistics (PERF.32), and which
  actions cannot be undone. The page is a Markdown file in docs/help/
  built into the Help section, linked from each page it describes, and
  updated by any later change to those features (a feature item is not
  done until its help text is).
  Prerequisite: DOC.6. Related: DOC.6.

- [ ] **DOC.15 A glossary and units help page: coordinates, scales, sector paths, object IDs, time and the in-universe wording**
  Boss (2026-10-09 23:53Z): "Add to-do items to build static
  documentation pages for all features accessible through the web
  interface." Done: the coordinate frames (galactic, sector, system),
  the unit ladder and when a unit switches (UX.22, UX.23, UX.78, UX.81),
  sector addresses and paths, the object ID (GEN.170), the epoch and
  time scales, and the in-universe terms (UX.42), each with a worked
  example. The page is a Markdown file in docs/help/ built into the Help
  section, linked from each page it describes, and updated by any later
  change to those features (a feature item is not done until its help
  text is).
  Prerequisite: DOC.6. Related: DOC.6.

- [ ] **DOC.16 The API help page for visitors: what the API is, how to get a key, and where the reference lives**
  Boss (2026-10-09 23:53Z): "Add to-do items to build static
  documentation pages for all features accessible through the web
  interface." Done: a short non-developer page on the API (/api): what
  it can do, getting a key, rate and upload limits, versioning (API.22),
  and a link to the full reference in docs/api.md; kept in step with
  that file. The page is a Markdown file in docs/help/ built into the
  Help section, linked from each page it describes, and updated by any
  later change to those features (a feature item is not done until its
  help text is).
  Prerequisite: DOC.6. Related: DOC.6.

## VIEW: The view from a planet

- [ ] **VIEW.1 View from a planet**
  **Research first.** Boss: "view-from-planet will have to do
  calculations on colors and A LOT Of stuff, so make special note of
  that, it will need a full research pass." Before any code for VIEW.2
  and VIEW.3, Boss wants a research session with him "into exactly how
  one would do that". VIEW.2 and VIEW.3 are blocked on it; VIEW.4
  is not.
  Research (2026-10-09, sky-view.md): split VIEW.2 into (a) photometry
  (magnitude, colour, extinction: `physics/photometry.py`), (b) the sky
  catalogue and its backfill tiers (`galaxy/sky.py`), (c) the galactic
  band (`galaxy/sky_band.py`), (d) orientation and horizon; each ships
  alone, in that order. Replace the "millions of stars per view"
  performance worry by "about 9,000 resolved stars at V <= 6.5 plus a
  band".

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
    Research (2026-10-09, sky-view.md): `bright_stars` (the 1,000 L_sun
    floor) plus the 100 ly backfill cover only about 8% to 10% of the
    naked-eye sky: a tier table (`BRIGHT_STAR_SKY_TIERS`, floors from
    visibility) or an angle-bin scatter like `scatter_layer` is needed.
    GEN.98 is built, so measure `backfill_cells` per cell on a real
    database first (the 2026-10-08 timing run measured the sector fill,
    not the backfill). Light-time: the position effect is under a tenth
    of a pixel except for hypervelocity stars; the part worth doing is
    the star's age and luminosity at the retarded time and
    supernova-remnant visibility; add the vectorised first-order
    `sky_vectors` (sky-view.md 2.2) beside VIEW.5's per-body
    `apparent_position`. GEN.104 is a prerequisite for horizon views and
    time of day (a whole-sky chart does not need it). Open questions for
    Boss (defaults taken): exact stars for V <= 6.5 with tiers built by
    an angle-bin scatter and a background glow for the rest, the
    planet's page saying "sky incomplete" until its tiers ran; limiting
    magnitude 6.5 for every planet with horizon extinction and twilight
    only; physical colour from blackbody, desaturated for faint stars,
    with a "boost colour" option; galactic plane dust A_V 1.0 per kpc;
    the planet's north for the horizon is its spin pole.
    Research (2026-10-09, globular-clusters.md): the globular clusters
    of generated galaxies come from GEN.164.

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
    Research (2026-10-09, sky-view.md): VIEW.3 can register a planet-sky
    card with ADM.44's card renderer when it lands. Answer the open
    questions as in sky-view.md: whole-sky first (Hammer-Aitoff),
    horizon dome (stereographic) second; constellations by the
    seeded-partition algorithm over the 500 to 900 brightest stars,
    stored per system; PNGs cached (key: system uid, epoch bucket,
    projection, size, m_lim, version key) and served by the API, written
    with numpy and `zlib`; constellation lines and names as an SVG
    overlay. Neighbour galaxies are extended sprites in a planet's sky
    (sky-view.md 2.10).

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
    naming key (GEN.67), not from sliced word lists. Constellations are
    neither stars, sectors nor bodies in a system, so they use the codec
    (GEN.67's default).
    Research (2026-10-09, sky-view.md): add the sub-items (1)
    `constellation` type code 12 in `names/object_id.py` and its ID
    layout (type bits plus the first 70 bits of SHA-256(galaxy seed ||
    "constellation:" system uid || index)), (2) the Bayer lettering
    rule, (3) record that galaxy-wide uniqueness comes from the codec's
    injectivity. Mark the "languages and sky cultures" source list as
    superseded by the 2026-10-07 plan, and the licence worry (Stellarium
    CC BY-SA 4.0) as moot; the open questions on licensing,
    transliteration and per-planet uniqueness are moot under the codec
    plan.

- [ ] **VIEW.6 Settle the handedness of the generated galaxy before mapping real sky coordinates**
  The generated galaxy rotates counterclockwise about +Z and the real
  one clockwise about the north pole. Settle which frame sky charts and
  the neighbour galaxies use. Open question for Boss (default): keep the
  project's rotation; for neighbour galaxies use the proper embedding of
  multiple-galaxies.md 2.2 (longitude 90 along the rotation, b flipped);
  a sky chart's own l and b are measured in the generated frame and
  labelled as such.
  Prerequisites: none.
  Design: [docs/design/sky-view.md](design/sky-view.md)

- [ ] **VIEW.7 Sky tiers and the backfill cost measurement**
  Build the tiers that complete the naked-eye sky (see VIEW.2) and
  measure `backfill_cells` per cell on a real database first.
  Prerequisite: VIEW.1.
  Design: [docs/design/sky-view.md](design/sky-view.md)

- [ ] **VIEW.8 A dust model for the sky: the galactic extinction law and `nebulae.extinction_av`**
  The galactic law of sky-view.md 2.5 plus the use of
  `nebulae.extinction_av`.
  Prerequisite: VIEW.1.
  Design: [docs/design/sky-view.md](design/sky-view.md)

- [ ] **VIEW.9 Replace the BC_V table with a published relation (Flower 1996 or Torres 2010)**
  Before any star count depends on the bolometric correction.
  Prerequisites: none.
  Design: [docs/design/sky-view.md](design/sky-view.md)

- [ ] **VIEW.10 Declare `pillow` in `setup.py` if it is used, and add a HYG calibration test**
  Pillow is already installed through scikit-image (pillow 11.3.0 below
  Python 3.10, 12.3.0 above) and no repo code imports it; declare it if
  Pillow is used (open question for Boss, default: SVG overlay, Pillow
  only for a self-contained PNG with text). A calibration test that
  reads HYG from a developer machine (not committed).
  Prerequisites: none.
  Design: [docs/design/sky-view.md](design/sky-view.md)

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
  Research (2026-10-09, population-and-politics.md): the weights barely
  matter (equal weights differ by 0.03 TL on average, 0.27 at most, same
  ranking), so no calibration work on them.

  - [ ] **POP.8 A tech-level design from the six domain indices**
    The doc scores six domains 0 to 7 (Energy, Materials, Information,
    Medical, Propulsion, Defense) and weights them: TL = 0.25 EI + 0.20
    MI + 0.20 II + 0.15 MeI + 0.10 PI + 0.10 DI. Done: how each index is
    drawn from a species' age, era and world (population.py's
    civilization age and era), written into population-and-politics.md.
    Research (2026-10-09, population-and-politics.md): add to Done that
    `civilization_age_years` is years since the industrial transition
    (Energy index about 1.5 to 2), that indices are integers 0 to 7, and
    that Earth 2026 is the calibration point (TL about 3.1; E 3.0, M
    3.0, I 3.4, Me 3.2, P 3.0, D 3.2, a judgment for Boss to review).
    The parameter table, band names and limiting-domain decision are in
    the design doc, so the written-into-population-and-politics.md
    clause is met once Boss accepts them. Open questions for Boss
    (defaults taken): age zero is the industrial transition; old
    civilizations are not all near the top (the oldest median is about
    6.1); the weighted sum is the only stored TL, with the limiting
    domain and a lopsided flag shown; label "Tech band" and "Era" on the
    page since both use "Industrial".
    Design: [docs/design/population-and-politics.md](design/population-and-politics.md)

  - [ ] **POP.9 A tech level generated for every technological species**
    Done: every technological species stores its six indices and its
    tech level (columns), shown on the species page and searchable;
    existing species get theirs on the next population pass.
    Research (2026-10-09, population-and-politics.md): six TINYINT
    UNSIGNED indices with CHECK 0..7, `tech_level` DECIMAL(3,2)
    (indexed) and `tech_source` ('generated' or 'admin'). The draw uses
    its own seeded stream (`draw.Stream(f"tech:{planet_id}")`) so stored
    ages do not change; `refresh_civilizations` fills NULL indices,
    which gives "existing species get theirs on the next population
    pass". Searchable by `tech_level` range and by each index. The
    limiting domain (lowest of E, M, I) and the lopsided flag (max minus
    min at least 3) are derived, not stored; show the tech band name
    beside the age era on the species page.
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
  Research (2026-10-09, population-and-politics.md): may split into (a)
  types and typed fields with admin CRUD and API, (b) an affiliation FK
  plus name snapshot, (c) five ratings, the chip component and the icon
  allowlist, (d) optional "suggest ratings". The rating chip shows an
  icon and the word, not colour alone (palette in the design doc). The
  icon allowlist and vendoring (60 to 100 Bootstrap Icons, MIT licence
  file kept) is a small separate task under UX.43 or UX.49. Open
  questions for Boss (defaults taken): facility types refine the five
  existing kinds rather than replace them; ship the mechanism plus one
  example type.
  Design: [docs/design/population-and-politics.md](design/population-and-politics.md)
