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
- The project is in early alpha (Boss, 2026-10-15 5:27 PM EDT): the
  current database schema is not a constraint on any item. Redesign
  tables, columns and layouts as the work needs, with one migration or a
  fresh galaxy.

### Links

- This file is the only place that links an item to its design document,
  on its own line under the item:
  `Design: [docs/design/x.md](design/x.md), section N`. Design documents
  don't list item IDs.
- The version files in [plan/](plan/) are the other place that lists
  item IDs: each version file has a table of its items per thread and
  phase (see "Plan: versions" below). File a new item into one thread's
  table of one phase in the same PR, and delete its row when it ships.
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

## Plan: versions

Boss (2026-10-01 23:53Z) asked for the open work to be planned in
phases, with this file as the master and index. Phase 0 (every bug and
the groundwork) was finished on 2026-10-09 and is now version 8.0. Boss
(2026-10-16, 1:07 PM EDT) rebuilt the tiers as a version map: one version
per release, one plan file per version, starting at 8.1, which is the old
Phase 1 we were building. A version has one or more phases (Boss, 1:09 PM
EDT: one version can have several phases when there is too much, with a
logical breakdown), each with a theme and its build threads. Every open
item below is in exactly one thread of one phase of one version, and a
blocker is done in the phase before: in an earlier position of its own
thread, in an earlier phase of its version or in an earlier version; the
threads of one phase are independent work trees that lanes can build side
by side. Each version file gives the goal, the phases with their themes,
the threads with their order and prerequisites, and any open questions;
[plan/notes.md](plan/notes.md) keeps the research notes, the files
several items share and the judgment calls, and
[plan/release-roadmap.md](plan/release-roadmap.md) holds each version's
headline, what it needs and the ideas not yet filed. The whole-sphere work
(GEN.233 onwards, the phenomena catalogue and the fill) has its own plan,
[plan/intergalactic-space-plan.md](plan/intergalactic-space-plan.md). This
file keeps each item's full text. A new item goes into a version file's
table in the same PR that files it.
When an item ships, its row leaves the version file and is copied into
[plan/progress-ledger.md](plan/progress-ledger.md) in the same PR: the
progress bar on every version of the version plan counts finished items
from that ledger (Boss, 2026-10-16), and the build fails if a retired
item has no row there.

| Version | Plan | Headline |
|---|---|---|
| 8.1 | [version-8-1-zoom-and-scatters.md](plan/version-8-1-zoom-and-scatters.md) | Infinite zoom and the new generate scatters |
| 8.2 | [version-8-2-accounts-and-invitations.md](plan/version-8-2-accounts-and-invitations.md) | User accounts with email invitations |
| 8.3 | [version-8-3-moving-galaxy.md](plan/version-8-3-moving-galaxy.md) | The moving galaxy |
| 8.4 | [version-8-4-new-worlds.md](plan/version-8-4-new-worlds.md) | New worlds |
| 8.5 | [version-8-5-courses-around-gravity-wells.md](plan/version-8-5-courses-around-gravity-wells.md) | Courses that bend around gravity wells |
| 8.6 | [version-8-6-edit-everything.md](plan/version-8-6-edit-everything.md) | Edit everything |
| 8.7 | [version-8-7-console-and-help.md](plan/version-8-7-console-and-help.md) | The pilot's console and built-in help |
| 8.8 | [version-8-8-deep-sky.md](plan/version-8-8-deep-sky.md) | The deep sky |
| 8.9 | [version-8-9-generate-anywhere.md](plan/version-8-9-generate-anywhere.md) | Generate anywhere |
| 9.0 | [version-9-0-view-from-a-planet.md](plan/version-9-0-view-from-a-planet.md) | The view from a planet |
| 9.1 | [version-9-1-gravity-map.md](plan/version-9-1-gravity-map.md) | The gravity map |
| 9.2 | [version-9-2-politics-and-territories.md](plan/version-9-2-politics-and-territories.md) | Politics and territories |
| 9.3 | [version-9-3-ships-and-stations.md](plan/version-9-3-ships-and-stations.md) | Ships and stations |
| 10.0 | [version-10-0-roles-and-permissions.md](plan/version-10-0-roles-and-permissions.md) | Roles and permissions |
| later | [version-later-future-features.md](plan/version-later-future-features.md) | Future features (not scheduled) |

Version 8.0 is the finished Phase 0 and has no file. Version "later" holds
possible future features that belong to no version yet (multiple
galaxies).

Version: the next revision is 8.1 (Boss, 2026-10-09 18:12Z: no 8.1 until
Phase 1 is complete; Phase 1 is now version 8.1). Boss (23:07Z) put the
fly-through Galaxy Map (MAP.146 and MAP.147 to MAP.155) into it as the
main feature that earns the 8.1 bump, and on 2026-10-16 (1:07 PM EDT) added
the new generate scatters, accurate galactic generation and scattering
including the halo (GEN.301), as the second 8.1 headline. 8.1 is complete
only when all its phases are done; until then every release stays on 8.0
(`REVISION_HOLD` in `scripts/bump_version.py`), and the TODO thread flips
the hold when Boss declares it complete.

Headline features (Boss, 2026-10-16, 1:07 PM EDT): 8.1 infinite zoom
(MAP.146) and the new generate scatters (GEN.301), 8.2 user accounts with
sign-up by an email invitation only an admin can send (USR.1), 8.3 the
orbital update system (GEN.105), 8.5 courses that bend around gravity wells
(NAV.59), 8.6 edit everything (ADM.57), 9.1 the gravity map (MAP.167),
10.0 roles and permissions with user edit requests and moderators
(USR.25, Boss 1:43 PM EDT). 9.0 has no chosen headline yet (the view from
a planet is proposed); multiple galaxies is a possible future feature. A
headline feature carries a `Headline: VERSION (name).` line in its item;
`scripts/build_todo_docs.py` reads it and shows, in every plan page, each
headline with all the open items it needs in build order. Every release in
CHANGELOG.md opens with a **Headline:** line (Boss, 2026-10-10); the TODO
thread adds one to each newly stamped release.

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
| 13 | Planet information without the Markdown render | UX.30 (absorbed into UX.105) |
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
  rate). Question (decided below): a banner at the top of every page, or
  only on
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
  Decided (2026-10-10): Boss (19:04-19:12Z cards, recommended option):
  the banner shows on every page, and it also covers command-line runs,
  so `generate.py` and the CLI must write the same lock and publish
  progress.
  Plan change (2026-10-15): the banner and its public endpoint are built on
  the Redis run record and ETA function of PERF.96 (done, PR #1553), never on the MySQL job
  tree that PERF.100 deletes; PERF.99 only repoints the Generate page,
  Queue page and job log. PERF.96 has merged, so it is unblocked (held until Boss releases the lanes).
  Narrowed (2026-10-16 review): the public status endpoint over the run
  record, the banner on every page and the hour-rounding rule;
  command-line runs appear because PERF.98 puts them on the run record.

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
  Question (decided below): does "Habitable" stay as a chip?
  Decided (2026-10-10): Boss (19:12Z): keep the Habitable chip; revisit
  with a habitable-class-only chip later.

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
  no-JS check that each form still posts. Decided (Boss, 2026-10-10
  18:44Z, defaults approved; default): Shoelace for select, checkbox,
  switch and the Generate and
  admin forms; native for sign-in, account, 2FA code, hidden and submit
  fields.
  Split (2026-10-16, Boss 1:07 PM EDT: look at each item for smaller
  chunks and shared work): this item is now the umbrella and is done
  when its subitems are: UX.124, UX.125, UX.126. About 170 form controls
  in 19 templates (generate.html 27, admin_facility_types.html 26,
  system.html 22, account.html 15, generate_system.html 12 and others),
  both themes and phone width, plus three cross-cutting pieces (a common
  look from one token, server validation errors on the right field, a
  no-JS posting check).
  Prerequisites: UX.124, UX.125, UX.126.
  Design: [docs/design/library-migration.md](design/library-migration.md)

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
  Split (2026-10-16, Boss 1:07 PM EDT: look at each item for smaller
  chunks and shared work): this item is now the umbrella and is done
  when its subitems are: UX.130, UX.131, UX.132. UX.43 is a
  Boss-approved design brief followed by a restyle of the base template
  and every shared component in both themes, with contrast and reduced
  motion kept.
  Scope after the review (2026-10-16 review): this item is the umbrella
  of the design brief Boss approves (child a), the tokens, type and base
  template in both themes (child b) and the look of the shared
  components (child c).
  Prerequisites: UX.130, UX.131, UX.132.

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
    Needs the saved-list store, which holds the galaxy instance stamp
    that this item compares (2026-10-16 review).
    Prerequisites: UX.129.

  - [ ] **UX.47 A bookmark manager: list, rename, sort, group and delete**
    Done: a Bookmarks page listing every bookmark with its kind and
    place, where each can be renamed, grouped, sorted, opened on its map
    and deleted, and the list can be exported and imported as a file.
    Research (2026-10-09, map-ui-and-frontend-libraries.md): add import
    validation (size, count, kind whitelist, same-origin relative `url`
    only; `bookmarks.js` `urlOf` returns `entry.url` unvalidated into an
    `href`) and a move up / move down path (WCAG 2.5.7); drag and drop
    optional.
    Scope after the review (2026-10-16 review): the Bookmarks manager
    page over the saved-list store (sort, group, open on its map, export
    and import with validation, move up and down).
    Prerequisites: UX.46, UX.129.

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
    Scope after the review (2026-10-16 review): a system bookmark per
    region kept up to date, framing the region and highlighting its
    outline; the grouping and outline come from the endpoint of MAP.190,
    so this item needs MAP.190.
    No grouping or API of its own: the bookmark per charted group (kept
    up to date as sectors are charted) and the highlighted outline on
    open, from the library of MAP.190 (2026-10-16 review).
    Boss (2026-10-16, 7:24 PM EDT): dislikes the outlined blocks marking
    charted space and wants colour and opacity highlighting instead, so
    the region highlight here is a fill (the same as MAP.238), not an
    outline.
    Prerequisites: UX.47, MAP.190. Related: MAP.238.

- [ ] **UX.82 Theme checks after PR #800: SVG currentColor, two Shoelace contrast failures, alpha in --bg-subtle**
  Carry into UX.43 or a small follow-up: a check that every file under
  `vendor/shoelace/assets` and `icons.svg` uses `currentColor`; recheck
  the two light-theme Shoelace failures measured before #800 (primary
  hover 4.23, placeholder 4.41, both below 4.5:1); `--bg-subtle` is
  still a colour with alpha in the dark theme. UX.76's cause was
  `a.btn:visited` out-specifying the secondary-button rules.
  Prerequisites: none.
  Design: [docs/design/map-ui-and-frontend-libraries.md](design/map-ui-and-frontend-libraries.md)

- [ ] **UX.92 A Bookmark button for planets, moons and belts**
  Foundations lane 1 (2026-10-10, after NAV.8, PR #1095): bookmarks
  cover star systems only today. Add a Bookmark button on each body's
  anchor on the system page, so a planet, moon or belt can be
  bookmarked; the bookmark's link target is /object/<ref>, which
  redirects to the system page and the body's anchor. Done: a body can
  be bookmarked and listed in the bookmark manager (UX.47) with its
  kind; a test covers it. Small. Phase 2 unless Foundations lane 1 says
  otherwise.
  Prerequisite: UX.45. Related: UX.47.

- [ ] **UX.94 Gravity map wording: the legend, the mode names and the Map help text**
  Done: the legend names the unit (metric or Customary through the UX.36
  number formatting), the three modes have plain in-universe names with
  a one-line explanation each, the Map help dialog explains what the
  colours mean and what a saddle point is, and no TODO code appears in
  any of it (UX.93, done).
  Prerequisite: MAP.168. Related: MAP.167.

- [ ] **UX.100 Planet and moon lists show habitability as symbols, not notes**
  GitHub issue [#1217](https://github.com/dwhagar/planetGen/issues/1217)
  (Boss, 2026-10-11 01:34Z): "The habitability notes are too much in the
  planet and moon list. Just use just symbols instead." Overall rating a
  coloured star, pressure a coloured trapezoid, radiation a square,
  temperature a thermometer, chemistry a beaker, all simple vector
  shapes from the icon library. Superseded in part by the PHI-6 symbols
  of GitHub issue
  [#1221](https://github.com/dwhagar/planetGen/issues/1221) (UX.101);
  this item removes the notes from the lists and builds the symbol
  column, UX.101 swaps in the final shapes. Boss, 2026-10-11 01:35Z: the
  description text at the end should be in paragraph form, simply saying
  for each category where the planet falls and why, with no extra
  information (GEN.208). Owner: UX (moved from the
  retired Bugfixes lane 2, Boss 2026-10-15 05:39Z), in its queue order.
  Prerequisites: UX.101. Related: UX.30 (absorbed into UX.105), UX.32,
  UX.101, GEN.208 (done, PR #1433).
  Scope after the review (2026-10-16 review): removes the notes, builds
  the symbol column from UX.101's shapes and sets the description in
  paragraph form (GEN.208); needs UX.101.

- [ ] **UX.101 The PHI-6 symbols as vector shapes in the icon library, both themes, coloured by rating**
  GitHub issue [#1221](https://github.com/dwhagar/planetGen/issues/1221)
  (Boss, 2026-10-11 03:24Z), which supersedes the symbols of #1217.
  Symbols: Radiation, a small central circle with three radial arcs 120
  degrees apart; Temperature, an upright equilateral triangle bisected
  by a vertical line; Pressure, a central square flanked by two opposing
  inward-pointing equilateral triangles; Air Chemistry, an upright
  regular hexagon enclosing three horizontal lines; Water Chemistry, a
  downward-pointing equilateral triangle enclosing a centred circle;
  Suitability, a circle enclosing a nested diamond with a centre dot;
  Overall Habitability, a solid circle inside an equilateral triangle,
  inside an outer ring. Each symbol takes the colour of its rating
  (Blue, Green, Yellow, Red) so the list reads at a glance. Done: the
  six symbols and the overall symbol as vector shapes in the icon
  library, drawn on both themes with enough contrast, used in the planet
  and moon lists, and a test checks every rating colour. Needs the PHI-6
  model.
  Builds the six symbols and the overall symbol as vector shapes, both
  themes, one colour per rating (Blue, Green, Yellow, Red), with a test
  of every rating colour; it comes first and no longer needs UX.100, so
  no interim shapes are built and thrown away (2026-10-16 review).
  Prerequisites: none. Related: UX.100, GEN.207 (done, PR #1429).

- [ ] **UX.102 One logging standard for every output the user or an admin reads: the same timestamps, format, dividers and step markers everywhere**
  (Boss, 2026-10-14 14:16Z): "Add a TODO item to standardize all user
  output logs, timestamps, same formatting, same style, dividers, step
  markers, etc, etc, etc, adopt an industry based standard for file and
  log output." Today the Generate page log, the CLI text, the downloaded
  log, the job log files and the other commands each format their own
  lines (UX.98, done in #1369, was the first visible symptom: stage numbering and section
  separators differ). Done: a design note that lists every place the
  program writes a log or report a person reads (Generate log, CLI
  commands, job log files, the downloadable log, update.sh output, admin
  pages that show logs), picks the standard, and a build that moves all
  of them onto one logging layer, with tests that fail if a new output
  path formats its own lines. Decided (Boss, 2026-10-14 14:17Z, approved as proposed): use
  Python's standard logging module as the one entry point; every record
  carries an RFC 3339 UTC timestamp, a severity from the RFC 5424 syslog
  levels (debug, info, notice, warning, error, critical), a short
  component name and a message, with machine-readable key=value fields
  (logfmt style) for numbers such as counts and elapsed time; files are
  one record per line in that format, with an optional JSON Lines form
  for machines; the readable display in the web log and CLI is rendered
  from the same records, with one divider style, one stage marker form
  ("Stage N of M: name", "Stages N to M of M: name" for a run) and one
  wording for start, finish and failure, so a downloaded log and the
  screen always agree; timestamps are always UTC in files and local time
  only in the page display, labelled). Bugfixes lane 1 built
  UX.98 as the first slice (done, #1369: one stage formatter); the rest is
  a Phase 2 item for the interface builders.
  Split (2026-10-16, Boss 1:07 PM EDT: look at each item for smaller
  chunks and shared work): this item is now the umbrella and is done
  when its subitems are: UX.121, UX.122, UX.123, TEST.183. A design
  note, a library and a migration of every output in the program
  (Generate page log, CLI commands, job log files, the downloadable log,
  update.sh output, admin log pages) with an enforcement test: weeks of
  work across Python, web and shell, none of it one PR.
  Prerequisites: UX.121, UX.122, UX.123, TEST.183. Related: PERF.66.

- [ ] **UX.105 The system page builds every body's detail from the document model as HTML tables**
  Audit: /mnt/project-files/wording/text-output-audit-2026-10-15.md.
  Owner lane: Foundations Wording (2026-10-15), which writes the code
  for this item. Facts table, rating table (with the six symbols UX.101
  draws), life, orbit and moons, styled like the rest of the site in
  both themes and every size class, quantities as qty spans so the unit
  preference applies. Removes web/lib/mdrender.py, GET
  /api/systems/<id>/sections, render_system_sections, the markdown
  dependency from the web lock file, and the Markdown tests. The
  overview stays prose (Boss, 2026-10-10). The one-off system page
  previews the same HTML for a system that is never saved. This is the
  build of UX.30, absorbed here; UX.30's text and Boss's 2026-10-01
  quote are kept below.
    UX.30, Boss (2026-10-01 23:53Z): "Rework planet information displays
  to pull away from the markdown render and instead fits with the modern
  web UI as we've seen it so far." Plan change (Boss, 2026-10-15
  04:21Z): markdown is no longer the web default for any planet or
  system text; the richer tables replace it, and markdown and wikitext
  are generated on the fly only when the user asks.
  Prerequisites: UX.119, UX.120. Related: UX.30 (absorbed into this
  item), UX.32, UX.101.
  Split (2026-10-16, Boss 1:07 PM EDT: look at each item for smaller
  chunks and shared work): this item is now the umbrella and is done
  when its subitems are: UX.119, UX.120. A renderer and a page rewrite
  in one item, plus a removal list that UX.109 repeats: UX.105 'removes
  web/lib/mdrender.py, GET /api/systems/<id>/sections,
  render_system_sections, the markdown dependency and the Markdown
  tests' and UX.109 removes 'the /sections endpoint, mdrender,
  db/render.py render_system_sections' again.
  Deletions moved to UX.109 (2026-10-16 review): this item builds the
  renderer and moves the system page and the one-off page onto it; it
  removes nothing.

- [ ] **UX.107 Units in the description go through the ladders**
  Audit: /mnt/project-files/wording/text-output-audit-2026-10-15.md.
  Owner lane: Foundations Wording (2026-10-15), which writes the code
  for this item. The describer, the map side panel (gravity, star
  temperature, pressure of an airless body as 0 kPa in line with the
  habitability text), the CLI query lines, the class page radius ranges,
  and the wide-binary, comet, nebula and quasar numbers. The ladders
  themselves (mass, temperature, pressure, Gyr rename) are UX.22, UX.23
  and UX.81 (UX, all done); this item is the call sites outside the
  files they edit.
  Prerequisites: none. Related: UX.22 (done, PR #1504), UX.81 (done, PR #1504).

- [ ] **UX.108 Describe phenomena, sectors, galaxy and facilities in the same Document**
  Nebulae, remnants, quasar, rogue planets, asteroid fields, sector
  summaries and facilities are still described by their own per-class
  text. Build them into the same document model as UX.104 so one set of
  renderers (text, HTML tables, Markdown, wikitext) covers everything
  the program describes. Foundations Wording (filed 05:50Z, split out of
  API.25 (done, PR #1527), which was too big for one PR).
  Done also: each kind's text is written in short plain prose as it
  moves into the Document, so the phenomenon prose is not rewritten
  twice (GEN.221 drops it).
  Prerequisites: UX.104 (done, PR #1509).

- [ ] **UX.109 Remove the old per-class description text**
  Remove SystemConfig.MARKDOWN, the per-class to_paragraph_list and
  __str__ text, db/render.py render_system_sections, the /sections
  endpoint and mdrender, once nothing reads them. Foundations Wording
  (filed 05:50Z, split out of API.25, done PR #1527).
  Owns the one deletion of the Markdown machinery
  (`web/lib/mdrender.py`, the `/sections` route,
  `render_system_sections`, the markdown dependency and the Markdown
  tests) that UX.105 also listed; UX.105 now builds only the renderer
  and the page move (2026-10-16 review).
  Prerequisites: UX.120, UX.108.

- [ ] **UX.110 The phenomena count on the overview page is wrong (bug)**
  Boss (2026-10-15 06:41Z): "Add a bug to TODO the phenomena count in
  the overview page is wrong." The front page (UX.55) shows one count
  per kind, each linking to its list; the phenomena figure comes from
  `count_phenomena` in `db/query.py` (cached by `countcache`, PERF.69),
  and the Phenomena page shows its own total. Done: the overview figure
  equals the number of rows the Phenomena list shows, for the whole
  galaxy and after every fill, scatter or edit; a test fails when the
  two differ. Cause not yet found (it may be a stale cache, a missing
  table in the union or double counting). Owner: Database
  (count and cache). Queued behind the Spherical Sector Coordinates stop
  (Boss, 2026-10-15 06:39Z).
  Owner: Core 1, in its bug phase after the global top chain (Boss,
  2026-10-16, 2:43 PM EDT: "Core 1 will start working on bugs after";
  order set at 3:00 PM EDT: non-test bugs, then the database-version
  test failures, then the failures that hit full runs, then the one-off
  load flakes). Fix the cause; no skipping or quarantine.
  Prerequisites: none.

- [ ] **UX.113 The sector page shows the profile of any address in the reach, and heavy cells and far neighbourhoods warn before they run**
  Step 8 of the whole-sphere train. An uncharted address shows its
  profile computed from the address alone: expected systems, the
  population mix, the dark-matter density and mass in the cell, the
  ambient medium (disk ISM, hot halo, intergalactic) with its density
  and gas mass, and the expected counts per kind, like the interstellar
  debris count already shown. An asked, empty sector shows the same with
  "no objects". A heavy cell (a globular cluster core cell is 14 minutes
  to 2.5 hours) is a queued job with an ETA and a warning, never a
  blocking request, and a neighbourhood request out in the halo shows
  its expected system count and cost line before it runs. Units go
  through the ladders and no TODO codes appear in the text. Done: any
  address inside the reach opens a page that says what is expected there
  and why, nothing is stored by opening it, and the warnings show before
  a heavy request starts. Owner: Interface; the numbers come from Core's
  profile (GEN.235) through API.27.
  Built-on note (Core 1, 2026-10-16, 3:17 AM EDT, GEN.235 merged as
  #1784): call galaxy.profile.sector_profile(address, skeleton) for the
  profile and GalaxyBounds.describe_miss for the "outside the outline"
  line.
  Owner: Core 4, "Core 4 (2026-10-16)" (Boss approved the lane at 4:31
  PM EDT, 2026-10-16; queue in the lane report, queue position 3).
  Prerequisite: API.27. Related: GEN.237.
  Design: [docs/design/intergalactic-and-halo-contents.md](design/intergalactic-and-halo-contents.md)

- [ ] **UX.114 Scatter-row pages list the stored members and the map shows a system mass**
  Step 8 (design note 11 item 8). The sector and bright-star pages list
  the stored members of an unbuilt scatter row and the Galaxy Map shows
  the system mass; a `Members` column in the Phenomena table. Done: an
  unbuilt row with companions shows them and its total mass on its page
  and in the table. Owner: Interface (2026-10-15); asks Core for any
  schema change. Boss (2026-10-15 5:51 PM EDT): stellar companions of
  stars and star-like objects, including neutron star and black hole
  pairs, are generated at scatter time because they matter for gravity;
  planets, moons and belts stay at sector fill. Defaults taken
  (reversible): compact pairs promoted to the scatter, black hole
  companion fraction 1.25%, cap 12 stars (7 before the 6:08 PM EDT
  catalogue), hierarchical trees only, a row's members deleted when its
  sector builds it, pair rates from the middle of the evidence.
  Owner: Core 4, "Core 4 (2026-10-16)" (Boss approved the lane at 4:31
  PM EDT, 2026-10-16; queue in the lane report, queue position 6).
  Design: [docs/design/multi-star-systems.md](design/multi-star-systems.md)

- [ ] **UX.115 Interface update: every page and control works with the spherical grid, 96-bit IDs, the new phenomena and the region queries (stage after the sphere and generation stages)**
  Boss (2026-10-16, 12:52 AM EDT): "also we need to run an interface
  update on everything once this is done so all the controls work
  properly, then run the ADM update chain since our admin pages need
  desperately to be streamlined and fixed, they aren't very usable."
  Done: after the sphere and generation stages finish, an audit walks
  every public and map page (Galaxy Map, Sector Map and sector page,
  System Map and system page, planet, moon and belt pages, phenomena and
  nebula pages, search, nearby, navigation and courses, bookmarks,
  account, help, API pages) and, for each control, checks it against the
  spherical address (shell, band, wedge), 96-bit object IDs (inputs,
  URLs, copy and paste, display width, bookmarks, search), the new
  phenomena, scatter rows, satellites, clusters and wanderers beyond the
  disk, and the region queries (NAV.45, NAV.60, NAV.61); the result is a
  page-by-page list of controls with working, broken or missing, and
  every defect is filed as an item or fixed in the same pass when it is
  a one-line fix. The open map, navigation and interface items of the
  Interface queue (the MAP, NAV, UX and API items) finish before the
  audit starts, so it checks the finished pages. Queued, not started:
  Interface and Quality stay idle until the sphere and generation stages
  finish (Boss, 12:27 AM EDT) unless Boss says earlier.
  Prerequisites: GEN.240, NAV.61, API.27, UX.113, NAV.45.
  Related: UX.114, MAP.207, MAP.121, API.26, NAV.60.

- [ ] **UX.116 Sliders for every scatter and fill limit: one limit_slider with a number box, in three new segments under the scatter-limits row of the Generate page**
  Boss (2026-10-16, 1:09 AM EDT, in Core 3's thread): "Each control will
  be a slider and clearly labeled with a text box at the end of the
  slider that shows the number in the slider. Use existing layout and
  segments to design all the slider controls for scatter floors and
  ceilings." Done (spec: section 5.2 of
  docs/design/realistic-galactic-core.md): one limit_slider macro (a
  range with a number box at the end that shows the value and accepts a
  typed one; log scale for luminosities) replacing mass_slider; the
  baseline luminosity select becomes a slider too; three new segments
  under the existing scatter-limits row (Core, Halo and intergalactic
  space, Dense cells); core sliders follow the baseline until touched;
  defaults shown and a Reset per segment; the halo luminosity floor
  slider runs from 50 (default 500; Boss, 2026-10-15 10:06 PM EDT, was
  100). Owner Core 3 (Boss, 1:09 AM EDT).
  Prerequisites: DB.54. Related: GEN.234, ADM.43, UX.49.

- [ ] **UX.118 Slow pages show a progress bar while their data loads: the search page and the phenomena list pages**
  Boss (2026-10-16, 12:21 PM EDT): "have a progress bar in the search
  page so the user knows it is loading, that is one of the main issues".
  Owner Core 2, after PERF.105. When a page's data is being built by a
  queued job (PERF.105), the page shows a loading state with a progress
  bar (done of total, percentage and ETA from the status route, in the
  ladder units and the wording of the description layer), keeps working
  controls, can be cancelled, and swaps in the finished list without a
  full reload; the search page and the phenomena list pages (black
  holes, neutron stars and the rest of the kinds) are the first users,
  and every other slow list takes it by marking its route slow. The bar
  is accessible (role progressbar, a live text status, reduced-motion
  respected) and works at phone width. A failed job shows its error and
  a retry. Done: the progress bar component, the search page and the
  phenomena list pages on it, and browser tests for the loading,
  finished, cancelled and failed states in light and dark.
  Prerequisites: PERF.105. Related: PERF.104, UX.117, UX.3.

- [ ] **UX.119 An HTML renderer for the Document: facts, rating and composition tables, with unit spans, in both themes**
  Split from UX.105 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work). Add the fourth renderer beside text,
  markdown and wikitext in `planetgen/describe/`: the Document's facts,
  ratings (with the six symbols UX.101 draws), composition, life and
  orbit-and-moon blocks as HTML tables, every quantity as a `qty` span
  so the unit preference applies, text escaped, styled like the rest of
  the site in both themes and every size class. No page changes yet.
  Done: the renderer with snapshot and escaping tests for every block
  type, and the CSS, checked in light and dark and at phone width in one
  browser test over a fixture Document.
  Prerequisites: none. Related: UX.105.

- [ ] **UX.120 The system page and the one-off system page use it: every body's detail from the Document**
  Split from UX.105 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work). `web/lib/systempage.py` builds each
  body's detail (facts table, rating table, life, orbit and moons) from
  the Document through child a instead of from Markdown; the overview
  stays prose (Boss, 2026-10-10) but is plain paragraphs, not Markdown;
  the one-off system page previews the same HTML for a system that is
  never saved (the Markdown branch of `web/system_page.py` goes).
  Nothing is deleted here: the deletions are UX.109's. Done: the system
  page and the one-off page render from the Document; browser tests in
  light and dark, desktop and phone; no call to `markdown_to_html`
  remains in the two page modules.
  Prerequisites: UX.119. Related: UX.105.

- [ ] **UX.121 The logging layer and the written standard**
  Split from UX.102 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work). The design note that lists every
  place the program writes a log or report a person reads, and the one
  layer: Python's `logging` as the entry point; every record has an RFC
  3339 UTC timestamp, an RFC 5424 severity (debug, info, notice,
  warning, error, critical), a short component name, a message and
  logfmt `key=value` fields; files are one record per line with an
  optional JSON Lines form; one display renderer (used by the CLI and
  the web log) with one divider style, one stage marker form ("Stage N
  of M: name", "Stages N to M of M: name") and one wording for start,
  finish and failure; timestamps in files always UTC, local time only in
  the page display and labelled. Done: the note and the module with
  tests that the same records give a file line and a display line that
  agree.
  Prerequisites: none. Related: UX.102.

- [ ] **UX.122 The command line, the job log files and the downloadable log use the layer**
  Split from UX.102 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work). Move the CLI commands' output, the
  job log files (including the file the Redis run record writes at the
  end of a run, PERF.99) and the download onto the layer. Done: a run's
  CLI text, its log file and its download agree line for line; tests for
  each command family.
  Prerequisites: UX.121. Related: UX.102.

- [ ] **UX.123 The Generate page log, the admin log pages and update.sh output use the layer**
  Split from UX.102 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work). Move the web displays that show a log
  (the Generate page, the Queue page, the activity log pages) and
  `update.sh`'s own output onto the shared renderer and record format.
  Done: the screen and the downloaded log show the same lines; update.sh
  prints the same timestamp and marker forms; browser and shell tests.
  Prerequisites: UX.121. Related: UX.102.

- [ ] **UX.124 One form-field look and error slot, piloted on the search, nearby and list-filter fields**
  Split from UX.49 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work). Set `--sl-input-border-color` (and
  the select and checkbox equivalents) from the new `--line` token (3:1
  contrast), and add one Jinja macro that renders a labelled `sl-input`,
  `sl-select` or `sl-checkbox` with its name, required flag and an error
  slot for the server's validation error (ADM.21). Hidden, submit,
  password and one-time-code fields stay native (decided). Pilot on
  `search.html`, `nearby.html`, `nav_chart.html` and the datatable
  filters. Done: the macro and token, the pilot forms converted, values
  posted by name as before, a test that a validation error appears on
  the right field, a no-JS check that each pilot form still posts, and a
  browser test that fills and submits one form of each kind in light and
  dark and at phone width.
  Prerequisites: none. Related: UX.49.

- [ ] **UX.125 The Generate page and the one-off system page on the field macro**
  Split from UX.49 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work). Convert the controls of
  `generate.html` and `generate_system.html` to the macro (select,
  checkbox, switch and the numeric boxes), submitting and validating as
  before. Done: converted; the page's browser tests pass in both themes
  and at phone width; keyboard use and labels unchanged.
  Prerequisites: UX.124. Related: UX.49.

- [ ] **UX.126 The admin pages and the edit dialogs on the field macro**
  Split from UX.49 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work). Convert `admin.html`,
  `admin_facility_types.html`, `admin_stats.html`, `admin_queue*.html`,
  the system and sector edit dialogs (`system.html`, `sector.html`) to
  the macro; the sign-in, account, 2FA code, hidden and submit fields
  stay native (decided). Done: converted; admin browser tests pass; a
  form of each kind is submitted in a test.
  Prerequisites: UX.124. Related: UX.49.

- [ ] **UX.127 The Galaxy Map block menu offers Generate again**
  Split from GEN.231 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). Put `blockGenerateButtons` back
  in `generatebuttons.js` and show Generate in the Galaxy Map's block
  menu, starting child a's job and linking to its progress; disabled
  with a reason while another run is active. Done: a browser test over
  fixture data clicks Generate on a block menu and sees the job start,
  and sees it disabled while a run is active.
  Prerequisites: GEN.310. Related: GEN.231.

- [ ] **UX.128 Map side panels and other browser call sites take Quantity**
  Split from GEN.300 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). The browser side: gravity,
  temperature and pressure in the map side panels and every other script
  that keeps a value next to its unit use the JavaScript Quantity built
  from the same ladder table. Done: no browser call site multiplies by a
  unit factor; the Python and JavaScript outputs for the same value
  agree in the existing parity test.
  Prerequisites: GEN.299, UX.107. Related: GEN.300.

- [ ] **UX.129 Saved lists in the browser and in the account: one store with a local and an account backend**
  Found in the split review (2026-10-16, Boss 1:07 PM EDT: look at each
  item for smaller chunks and shared work, 1:33 PM EDT: shared
  foundations). Split `static/bookmarks.js` (447 lines; synchronous
  read() and write() over localStorage, one list per database) into a
  store and its views. The store offers list, find, add, rename, remove
  and onChange on a list kind (`bookmark` now, `course` for NAV.18)
  behind a backend interface that may be asynchronous: the local backend
  (today's per-database localStorage key) and, later, an account backend
  over the API; a hook for the galaxy instance stamp that UX.46 needs.
  The toggle button and menu keep working on top. Done: the existing
  bookmark toggle, menu, rename and delete behave as before through the
  store, a second list kind can be created without new storage code, and
  a fake asynchronous backend passes the same tests.
  Prerequisites: none. Related: UX.46, UX.47, USR.7, NAV.18, NAV.19.

- [ ] **UX.130 Console design brief: colours, type, panel shapes and iconography, approved by Boss**
  Split from UX.43 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work). A design note in docs/design (with
  its README row in the same PR) that fixes the pilot's-console look:
  the measured token table of map-ui-and-frontend-libraries.md section
  3.3 as the starting brief, `--line` as a new token beside `--border`,
  B612 Mono or the system monospace stack for readouts (default:
  self-host B612 Mono for numbers and designations only), selection and
  caution both amber but never in the same place, the rule that every
  colour is paired with a symbol, word or shape, and the habitability
  tiers reusing the rating chip so the site has one status component.
  Done: the note, with a mock page in both themes, approved by Boss.
  Prerequisites: none. Related: UX.43.

- [ ] **UX.131 Console tokens, type and base template in both themes**
  Split from UX.43 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work). Put the approved brief into the base
  template and the token files: the colour and `--line` tokens for light
  and dark, the monospace readout type and its self-hosted font, panel
  and frame shapes in the shell, reduced motion kept. Done: base.html
  and the token files use the brief, a contrast test (the project
  checker) passes on both themes for every token pair, and nothing
  outside the shell changes shape yet.
  Prerequisites: UX.130, UX.82. Related: UX.43.

- [ ] **UX.132 Console look for the shared components: buttons, tables, panels, status chips, icons and the maps' chrome**
  Split from UX.43 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work). Restyle the shared components onto
  the tokens of child b: buttons and form fields, tables, panels and
  dialogs, the one rating and status chip (with its symbol or word), the
  icon sprite and the map controls' chrome, each with a visible
  non-colour cue. Done: every shared component uses the tokens, the
  browser layout and colour tests pass in both themes, reduced motion is
  honoured, and a before-and-after screenshot of the main pages is
  attached.
  Prerequisites: UX.131. Related: UX.43.

- [ ] **UX.133 `GET /polities` and `GET /species` answer 404 before the population pass has run: show an empty page instead (bug)**
  Found by a lane (2026-10-16), filed by the TODO thread. On a galaxy
  that has had no population pass, `/polities` and `/species` return 404
  instead of a page. Done: both routes answer 200 with an empty list and
  a line that says the population pass has not run; the 404 stays for a
  polity or species id that does not exist. A test covers both routes on
  a galaxy with no pass. Fails on main: yes. Frequency: always on a
  galaxy before the first population pass. Reproducer: start a fresh
  galaxy and request both URLs.
  Prerequisites: none.

## MAP: Galaxy Map, Sector Map, System Map

MAP.2 with MAP.22 and MAP.23, MAP.15 and MAP.30 shipped in PR #234.

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
  Split (2026-10-16, Boss 1:07 PM EDT: look at each item for smaller
  chunks and shared work): this item is now the umbrella and is done
  when its subitems are: DB.60, ADM.94, MAP.218. MAP.119 stores an
  override 'that later fills use', but the table for it is DB.38, which
  sits in the Regions thread of the same phase (8.6-1), so MAP.119 has a
  blocker in a parallel thread.
  Moved to phase 2 with its children (the migration child sits in the
  Regions thread); its block or slab case is a Span region (see GEN.302)
  and the table is DB.38 (2026-10-16 review).
  Prerequisites: DB.60, ADM.94, MAP.218.

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
  Built in part (2026-10-09): The blocker fade rule is built as part of
  MAP.149 (PR #1141, static/nearfield.js). Nothing draws or steps to
  neighbouring regions yet, so this item stays open for the drawing, the
  arrow buttons and keys, and the browser test. Lane 2 builds it with
  MAP.150, MAP.148 and MAP.155.
  Scope after the review (2026-10-16 review): stepping to a neighbouring
  region at the same level (URL, breadcrumb, keys; re-read what
  MAP.200's nav pad already covers) and the Sector and System Maps
  adopting the shared drawing; MAP.141 draws the faint context regions
  once.

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
  Scope after the review (2026-10-16 review): the Select mode
  (radiogroup, select = stars, 10 px and 22 px reach with a
  disambiguation list) and the Bookmark action; NAV.49 adds the Waypoint
  action when waypoints exist.
  Decided (Boss, 2026-10-16, 3:00 PM EDT, decision card): this item goes
  after the MAP.146 subitems and before the label items. Owner: Core 2.

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
  Note (2026-10-10, nebula research): the nebula part of this item is
  built by MAP.179 (scattered nebulae as markers from level 8) and
  MAP.176 (regions); the black-hole and habitable-world markers stay
  here.
  Note (2026-10-14, GitHub issue sweep):
  [#1234](https://github.com/dwhagar/planetGen/issues/1234) (Boss,
  2026-10-14 00:19Z): "Add to the galaxy markers (optional) for
  habitable worlds, specific classes, etc available in various star
  systems from zoom levels in the galactic map. Also icons for space
  stations, etc." So the markers are optional per kind (MAP.123), cover
  specific classes as well as habitable worlds, appear at the zoom
  levels where the map budget allows, and stations and vessels get icons
  once they exist (USR.10, POP.11).
  Split (2026-10-16, Boss 1:07 PM EDT: look at each item for smaller
  chunks and shared work): this item is now the umbrella and is done
  when its subitems are: MAP.219, MAP.220, MAP.221. MAP.132 is the
  marker layer for the Galaxy Map plus the first kinds, but its notes
  already give the nebula kind to MAP.179 and MAP.176 (and MAP.179's
  text says it draws 'the far-zoom markers of MAP.132' while having no
  prerequisite on it), and the 2026-10-14 note adds station and vessel
  icons that need USR.10 and POP.11 in 9.3, a later version.
  Scope after the review (2026-10-16 review): this item is the umbrella
  of the marker layer (child a: the kind registry, per-kind toggles, the
  zoom budget and picking, shipping the black-hole kind), markers for
  habitable worlds and chosen classes (child b) and the station and
  vessel icons (child c, in 9.3).
  Prerequisites: MAP.219, MAP.220, MAP.221.

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
  Screenshot (Boss, 2026-10-14 13:37Z, copy in
  /mnt/project-files/screenshots/issue-758/): the top of the Galaxy Map
  page. The header band (breadcrumb, the title "Galaxy Map" and the "0
  placed sectors" pill) fills the full width but only its left fifth is
  used; the "Go to" bar below it is the only control on its row; and
  the right-hand panel holds only the help text until something is
  selected. The empty header space to the right of the title and the
  empty panel are the places to put the selection basics and menu items.

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
  The flight on the Galaxy Map (double-click to the thing clicked, "Go
  to", Enter) is MAP.226; this item keeps the select-and-open behaviour
  on every map (2026-10-16).

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
  Built in part (2026-10-09): The faint-context rule is built as part of
  MAP.149 (PR #1141, static/nearfield.js). Nothing draws the faint
  neighbours or the sectors above and below yet, so this item stays open
  for the drawing and the selection of the sector below. Lane 2 builds
  it with MAP.150, MAP.148 and MAP.155.
  Owns drawing the faint context regions once for every map (2026-10-16
  review).

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
  Note (2026-10-10, nebula research): Boss asked on 2026-10-10 21:46Z
  that nebulae be visible on the Galaxy Map; the work is split into
  MAP.172 to MAP.181 (MAP.182 was folded into GEN.241; the dark-theme contrast fix is MAP.175, the dust
  colour MAP.181, the grouped regions MAP.176 and the individual clouds
  MAP.178). This item keeps the fuzzy-boundary shape work.

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
  Scope after the review (2026-10-16 review): the sprite rule for
  neighbour galaxies reads `galaxy_components` for the Clouds and
  Sagittarius and the file of GEN.157 for the rest, and says how it
  hands over to the points of MAP.207.
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
  MAP.159 (packed binary tiles; needs MAP.158; one cache stamp bump with
  MAP.151's tile keys; MAP.154, the nested lists, merged in #1303), MAP.151 (needs ADM.29 from Foundations lane 1) and MAP.152 (needs MAP.148
  and MAP.150, done in #1382); this umbrella closes last.
  Headline: 8.1 (infinite zoom: the fly-through Galaxy Map).
  Prerequisites: MAP.153, MAP.155, MAP.157, MAP.159, MAP.225,
  MAP.226, MAP.227, MAP.228. Related:
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
  Owner: Core 2 (Boss, 2026-10-16, 11:32 AM EDT), right after MAP.152; the 8.1 headline.
  What remains after MAP.152 merged (2026-10-16) (Boss 2:29 PM EDT): the
  wheel toward the cursor (MAP.225), the double-click flight (MAP.226),
  the container named from the camera position (MAP.227) and the camera
  in the address (MAP.228); each was only in the Done line above and now
  has its own item. This umbrella closes last.
  Owner: Core 2 (Boss, 2026-10-16, 2:43 PM EDT: after the 8.1 phase 1
  top 5, Core 2 takes the infinite zoom tree: MAP.146 and its subitems,
  then the label items).

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
  Prerequisites: none. Related: MAP.146, MAP.150 (done, #1382), MAP.141, MAP.122,
  ADM.29.

- [ ] **MAP.160 Quantise the Galaxy Map GPU buffers (deferred)**
  Source: docs/design/galaxy-map-wire-format.md (section 5, from section
  5.6 of map-ui-and-frontend-libraries.md). The 52 bytes a star in
  section 5.6 is the GPU buffer, not the download, so quantising it (to
  about 20 to 24 bytes) saves GPU memory and the re-upload on each tile
  arrival, not network bytes; not measured on a real GPU. Done: colour,
  scalars and flags in `setStars` are byte-quantised, and a client that
  appends the new tile's stars instead of rebuilding every star is
  considered with it. Decided (Boss, 2026-10-10 18:44Z, defaults
  approved; default defer until a real
  phone or GPU measurement shows the upload matters): build it?
  Detail (2026-10-09, galaxy-map-wire-format.md): also considered here:
  the client appends the new tile's stars instead of rebuilding and
  re-uploading every star on each tile arrival (3,200 stars is a 166 KB
  upload, 11,000 is 572 KB each time), and `tileStars` copies an object
  per generated star. Optional.
  Prerequisites: none. Related: MAP.147, MAP.159, MAP.109.
  Design: [docs/design/galaxy-map-wire-format.md](design/galaxy-map-wire-format.md)

- [ ] **MAP.167 Gravity map: a heat map of the gravitational field inside a sector (the seminal feature of 9.0)**
  Boss (2026-10-10 20:57Z): 'Since we have point mass vectors to
  calculate orbits, I would like to be able to show a gravitational map
  of a sector, probably using colors and calculating it per zone. I
  would like to be able to show the gravitational gradients of the
  vector field as kind of a heat map of the inside of the sector.' Boss
  made it the seminal feature of 9.0 (phase 3 roll-out), with the
  foundations in phase 2. Umbrella: it is done when a viewer opens a
  sector, switches on the gravity layer and sees the pull, the well
  depth or the tidal strength as coloured zones with a legend, and the
  same layer on the System Map; subitems GEN.198, TEST.131, PERF.72 and
  GEN.199 (phase 2), then API.24, MAP.168, MAP.169 and UX.94 (phase 3),
  with MAP.170 (Galaxy Map layer) later.
  Note (2026-10-10): Boss (22:14Z) asked that the size of each gradient
  measure be built from the preset minimum movement values of the
  orbital update system (THRESHOLDS_M in physics/position.py: galactic
  0.01 mpc, about 2 AU; system 0.01 AU; planetary 100,000 km). The
  umbrella carries the rule for all layers: gradient step and zone size
  come from those thresholds, and the legend names the scale in use. See
  docs/design/gravity-map.md.
  Headline: 9.1 (the gravity map).
  Prerequisites: GEN.199, MAP.168, MAP.169, UX.94, API.24. Related:
  GEN.105, NAV.6.
  Design: [docs/design/gravity-map.md](design/gravity-map.md)

- [ ] **MAP.168 The Sector Map gravity layer: coloured zones for pull, well depth and tidal strength**
  Done: the Sector Map Menu has a Gravity layer with three modes; each
  zone of the shown slab is a translucent colour on a log scale with a
  colour-blind-safe ramp (transparency by strength, no fill above 50%
  opacity as for MAP.131, so stars and routes stay readable); the slab
  buttons step through the sector, a transparent-volume view shows every
  slab, Lagrange and saddle points are marked and pickable, and the
  colour choice is kept in the URL. It uses the one mapping engine and
  the API.24 grid; a browser test checks the layer draws, switches mode
  and keeps the URL.
  Scope after the review (2026-10-16 review): the gravity layer's data
  and drawing (zones, three modes, Lagrange points); its Menu entry, URL
  key and legend come from the layer registry.
  Prerequisites: API.24, MAP.125, MAP.217. Related: MAP.167, MAP.131.

- [ ] **MAP.169 The System Map gravity layer: the orbital plane as a heat map with Lagrange points and Hill spheres**
  Done: the System Map has the same Gravity layer on the system's
  orbital plane, with the Lagrange points of each pair of bodies marked
  and each Hill sphere outlined, switching modes and keeping the choice
  in the URL; a browser test checks it draws.
  Scope after the review (2026-10-16 review): the System Map's gravity
  data and Hill spheres only; it is the same registered layer on a
  second map.
  Prerequisites: MAP.168, MAP.217, GEN.289. Related: MAP.167, GEN.109.

- [ ] **MAP.170 The Galaxy Map gravity layer: the galaxy potential and region aggregates, coarse and optional**
  Later (phase 3+): a coarse gravity layer on the Galaxy Map drawn from
  the galaxy's smooth potential (GEN.115) and the per-region mass
  aggregates (MAP.151), so the spiral arms, bulge and halo show as
  wells; only if the cost study (PERF.72) says it is cheap enough.
  Prerequisites: MAP.168, MAP.217. Related: MAP.167.

- [ ] **MAP.171 The Galaxy Map keys a massive phenomenon's visibility to its map luminosity (replaces "always lit")**
  Boss (2026-10-10 21:38Z): "inject into the phenomena part of that, use
  the following mass/luminosity to allow us to store a 'map visibility'
  value for phenomena that aren't bright but have large mass. The idea
  being that a neutron star of say 2 solar masses should be given
  placement on the map when a star's brightness at 2 solar masses would
  be visible. Same for black holes, if a super massive black hole has a
  mass well above the brightest star, then it's visibility on the map
  should be keyed to that. This allows us to put items on the map that
  aren't bright enough to be seen but are of importance (major
  gravitational locations are always important to navigation)." Today
  the client law (MAP.148, built in PR #1160) gives black holes, neutron
  stars and quasars a fixed magnitude (PHENOMENON_MAGNITUDE = -40 in
  starmagnitude.js), so they stay lit at any distance, and the
  400-per-tile cap and the server floors treat them separately. Done:
  the tile carries the map luminosity of each phenomenon (in the wire
  format of MAP.157/158/161 and the packed tiles of MAP.159, as one more
  plane or field), the client turns it into an apparent magnitude
  exactly as for a star, so a 2 Msun neutron star appears when a 2 Msun
  star (about 17 Lsun) would, and a supermassive black hole far above
  the brightest star stays visible across the galaxy; the server's
  per-tile choice ranks phenomena and stars together by luminosity or
  map luminosity, so the 400 cap and the generated_star_floor_sol floors
  treat them alike. Colours stay as MAP.164 (black holes purple, neutron
  stars dark blue) and sizes as MAP.165 (by mass). Overlap: MAP.152 owns
  the distance-cut tile choice, the floor removal and the cap
  reconciliation in the same query, so this item follows MAP.152 and
  uses its tile choice. Not a prerequisite of the 8.1 headline
  (MAP.146). Owner: Foundations lane 2, queued right after MAP.152.
  Prerequisite: MAP.152. Related: MAP.148,
  MAP.164, MAP.165, MAP.159, MAP.146.
  Status (2026-10-09): DB.24 done (Foundations lane 1, 2026-10-14,
  #1317, schema v83): nebula_field holds every kept class M and N cloud,
  written by plan from the seed. Deviation: rows carry no object id
  until a sector stores the cloud (an id needs a nebula shape, about 112
  ms each, roughly 30 CPU hours for all clouds at default scale); a save
  fills uid and is_stored (the flag is is_stored because 'stored' is
  reserved on MySQL 8), and plan links clouds a galaxy already holds.
  Key unstored clouds by (cell_i, cell_j, cell_k, n).

- [ ] **MAP.172 Nebula cell aggregates in the region pyramid**
  Boss (2026-10-10 21:46Z): nebulae must be visible on the Galaxy Map,
  by colour or by combined regions; Phase 2 for the feature, foundations
  in Phase 1. From the research (docs/design/nebula-map-visibility.md
  section 9, PR #1180; report
  /mnt/project-files/research/nebula-map-visibility/report.md): in the
  MAP.151 pyramid (1 kpc down to about 8 pc) each cell stores the nebula
  count, volume, centroid, bounding radius, family shares, dominant
  class and cover, updated when a nebula is stored. Done: the
  aggregates, and a Monte Carlo check that the stored cover is within 5%
  of the true cover. Owner: Foundations lane 1, with MAP.151.
  Prerequisite: none. Related: MAP.151 (done, PR #1539), MAP.159.
  Design: [docs/design/nebula-map-visibility.md](design/nebula-map-visibility.md)

- [ ] **MAP.173 The Galaxy Map tile carries a nebula layer: regions and singles, nested like the star lists**
  Boss (2026-10-10 21:46Z): nebulae must be visible on the Galaxy Map,
  by colour or by combined regions; Phase 2 for the feature, foundations
  in Phase 1. From the research (docs/design/nebula-map-visibility.md
  section 9, PR #1180; report
  /mnt/project-files/research/nebula-map-visibility/report.md): regions
  plus single nebulae per tile, nested like MAP.154 so a parent list is
  a subset of its child, in the JSON tile first and the packed tile
  (MAP.159) after; budget 600 regions per tile ranked by volume; one
  cache-stamp bump together with MAP.147 and MAP.151. Done: the layer,
  the budget and a test on a seeded galaxy. Owner: Foundations lane 2,
  with MAP.151 and MAP.159.
  Prerequisites: MAP.172. Related: MAP.147, MAP.151 (done, PR #1539),
  MAP.159.
  Tile-format changes are one chain with MAP.171 and MAP.195 (one
  cache-stamp bump), so this item moved into the Fly-through thread;
  MAP.172 and GEN.202 moved to phase 1 as foundations (2026-10-16
  review).
  Design: [docs/design/nebula-map-visibility.md](design/nebula-map-visibility.md)

- [ ] **MAP.174 Nebula cover in the stage-view cell statistics, for filled and unfilled cells**
  Boss (2026-10-10 21:46Z): nebulae must be visible on the Galaxy Map,
  by colour or by combined regions; Phase 2 for the feature, foundations
  in Phase 1. From the research (docs/design/nebula-map-visibility.md
  section 9, PR #1180; report
  /mnt/project-files/research/nebula-map-visibility/report.md): the cell
  statistics behind the stage view (and the Color by switch) gain the
  share of each cell inside a nebula, for filled cells from stored
  nebulae and for unfilled cells from the field. Done: the statistic,
  with a test that filled and unfilled cells agree on a seeded galaxy.
  Owner: Foundations lane 2.
  Prerequisites: DB.24, MAP.172. Related: MAP.131, MAP.151 (done, PR
  #1539).
  Narrowed (2026-10-16 review): exposes MAP.172's cover in the
  stage-view statistics payload (filled from stored nebulae, unfilled
  from the field) plus the filled-versus-unfilled agreement test.
  Design: [docs/design/nebula-map-visibility.md](design/nebula-map-visibility.md)

- [ ] **MAP.176 Nebula regions on the Galaxy Map: grouped soft sprites in the dominant family colour**
  Boss (2026-10-10 21:46Z): nebulae must be visible on the Galaxy Map,
  by colour or by combined regions; Phase 2 for the feature, foundations
  in Phase 1. From the research (docs/design/nebula-map-visibility.md
  section 9, PR #1180; report
  /mnt/project-files/research/nebula-map-visibility/report.md): Phase 2,
  visible feature. fixed-grid groups at 16 px cells drawn as a soft
  sprite in the dominant family colour with an outline, split by fade
  (the MAP.153 rule) as the camera zooms in, budget 600, per-kind
  toggles (MAP.123) and a legend block. Done: the layer, the toggles and
  the legend, with a browser test on a seeded galaxy. Default: 16 px
  cells and a 600 budget. Owner: Foundations lane 2.
  Note (2026-10-10, nebula density research, PR #1188): add a stress
  test with 3,000 candidate regions keeping the 600 with the largest
  volume; the region counts stay above the 600 budget from 8 kpc to 1
  kpc even after the density cut, so tile size is set by the budget, not
  by density.
  Prerequisite: MAP.173. Related: MAP.142, MAP.123, MAP.153.
  Design: [docs/design/nebula-map-visibility.md](design/nebula-map-visibility.md)

- [ ] **MAP.177 Color by "Nebula cover" on the Galaxy Map**
  Boss (2026-10-10 21:46Z): nebulae must be visible on the Galaxy Map,
  by colour or by combined regions; Phase 2 for the feature, foundations
  in Phase 1. From the research (docs/design/nebula-map-visibility.md
  section 9, PR #1180; report
  /mnt/project-files/research/nebula-map-visibility/report.md): Phase 2,
  visible feature. a new choice in the MAP.131 Color by switch: a
  single-hue 5-bin ramp (0, under 5, 5 to 15, 15 to 30 and 30% or more),
  legend "share of the block inside a nebula", hover shows the
  percentage, works on unfilled cells. Off by default. Done: the choice,
  legend and hover with a browser test. Owner: Foundations lane 2.
  Prerequisite: MAP.174. Related: MAP.131.
  Design: [docs/design/nebula-map-visibility.md](design/nebula-map-visibility.md)

- [ ] **MAP.178 Field nebula clouds drawn one by one from the nebula table**
  Boss (2026-10-10 21:46Z): nebulae must be visible on the Galaxy Map,
  by colour or by combined regions; Phase 2 for the feature, foundations
  in Phase 1. From the research (docs/design/nebula-map-visibility.md
  section 9, PR #1180; report
  /mnt/project-files/research/nebula-map-visibility/report.md): Phase 2,
  visible feature. individual clouds from `nebula_field` drawn as
  sprites on the 1 to 4 px ramp of MAP.155, then as meshes,
  de-duplicated with stored nebulae. Done: a browser test finds one of
  each family on a seeded galaxy. Owner: Foundations lane 2.
  Note (2026-10-10, nebula density research, PR #1188): the far sprite
  is drawn at the equivalent radius (0.49 of the bounding radius) with
  the 2 px floor and the MAP.155 ramp; the bounding radius stays for
  picking and for the mesh load.
  Note (2026-10-10, nebula density research, PR #1188): realistic cover
  is about 2% at the median, so the ramp bins become 0, under 0.5%, 0.5
  to 1.5%, 1.5 to 3% and 3% or more (or per-view percentiles).
  Prerequisite: MAP.173. Related: MAP.142.
  Design: [docs/design/nebula-map-visibility.md](design/nebula-map-visibility.md)

- [ ] **MAP.179 Scattered nebulae as markers from level 8**
  Boss (2026-10-10 21:46Z): nebulae must be visible on the Galaxy Map,
  by colour or by combined regions; Phase 2 for the feature, foundations
  in Phase 1. From the research (docs/design/nebula-map-visibility.md
  section 9, PR #1180; report
  /mnt/project-files/research/nebula-map-visibility/report.md): Phase 2,
  visible feature. planetary nebulae and remnants from the scatter shown
  as the far-zoom markers of MAP.132 from level 8. Done: the markers,
  with a browser test. Owner: Foundations lane 2.
  Draws the far-zoom markers of the marker layer (first child of
  MAP.132), so it needs it (2026-10-16 review); it keeps the
  planetary-nebula and remnant kind.
  Prerequisites: GEN.202, MAP.173, MAP.219. Related: MAP.132.
  Design: [docs/design/nebula-map-visibility.md](design/nebula-map-visibility.md)

- [ ] **MAP.180 Picking, hover and fly-to for nebula regions and cover**
  Boss (2026-10-10 21:46Z): nebulae must be visible on the Galaxy Map,
  by colour or by combined regions; Phase 2 for the feature, foundations
  in Phase 1. From the research (docs/design/nebula-map-visibility.md
  section 9, PR #1180; report
  /mnt/project-files/research/nebula-map-visibility/report.md): Phase 2,
  visible feature. hover and click on a region or on cover: a region
  click flies to the zoom at which it splits (no page change), with a
  keyboard path. Done: hover, click and keyboard paths with browser
  tests. Owner: Foundations lane 2.
  Prerequisite: MAP.176. Related: MAP.150 (done, #1382).
  Design: [docs/design/nebula-map-visibility.md](design/nebula-map-visibility.md)

- [ ] **MAP.181 Dust colour for the dark nebula family on both themes**
  Boss (2026-10-10 21:46Z): nebulae must be visible on the Galaxy Map,
  by colour or by combined regions; Phase 2 for the feature, foundations
  in Phase 1. From the research (docs/design/nebula-map-visibility.md
  section 9, PR #1180; report
  /mnt/project-files/research/nebula-map-visibility/report.md): Phase 2,
  visible feature. umber (about #a9805a at 0.45 opacity) with a light
  rim on the dark theme, chosen with the contrast checker of MAP.175,
  and a matching light-theme colour. Done: the tokens and the 3:1 test
  pass for the new colours. Owner: Foundations lane 2.
  Prerequisite: MAP.175. Related: MAP.142.
  Design: [docs/design/nebula-map-visibility.md](design/nebula-map-visibility.md)

- [ ] **MAP.183 Optional, default not built for the first release: near-view P and Q nebulae from the seeded field below about 250 pc**
  Boss (2026-10-10 22:01Z): review the cloud density against the real
  Milky Way and decide how to manage it so things stay visible. From the
  research (docs/design/nebula-density-vs-reality.md section 9, PR
  #1188; report
  /mnt/project-files/research/nebula-map-visibility/density-report.md):
  classes P and Q are implicit in the seed (291,000 clouds of 0.2 to 1.6
  pc) and are invisible beyond about 250 pc. Done if built: below about
  250 pc the near view draws them from the seeded field with a per-tile
  cell budget. Default: not built for the first release. Owner:
  Foundations lane 2 if Boss asks.
  Prerequisites: none. Related: DB.24, MAP.173.
  Design: [docs/design/nebula-density-vs-reality.md](design/nebula-density-vs-reality.md)

- [ ] **MAP.189 The Galaxy Map opens flown in to the created sectors**
  GitHub issue [#1232](https://github.com/dwhagar/planetGen/issues/1232)
  (Boss, 2026-10-13 23:16Z): one big reason for the infinite zoom is to
  restrict the field to the created sectors: "Ideally when galactic map
  is loaded, it should fly into the zoom where the user can see all of
  the sectors in a region." Boss: get one group (the neighbourhood fill)
  working first. Done: when the Galaxy Map loads and sectors have been
  created, the camera flies in to a view that frames the group of
  created sectors (their bounding box plus a margin) and the field shows
  only that region at the start; the user can still zoom out to the
  whole galaxy. A browser test loads a galaxy with one neighbourhood and
  finds the camera framed on it. Owner: Foundations lane 2, after
  MAP.152 (MAP.150 done, #1382); part of the fly-through headline.
  Scope after the review (2026-10-16 review): fly the camera in to the
  frame of the largest group of created sectors and show only that
  region at the start; the frame comes from the charted-regions library
  (MAP.190), one group being the first case.
  Prerequisite: MAP.152. Related: MAP.146, MAP.151 (done, PR #1539), MAP.190.

- [ ] **MAP.190 Several separate regions of created sectors: which one the Galaxy Map opens on**
  GitHub issue [#1232](https://github.com/dwhagar/planetGen/issues/1232)
  (Boss, 2026-10-13 23:16Z): "Multiple regions will have to be handled
  with an algorithm of some kind, not sure what yet." Done: when created
  sectors form more than one separate group, the map has a rule for
  which group it opens on and a way to jump between the groups (a list
  of regions). Decided (Boss, 2026-10-14 13:27Z): the map opens on the
  largest group; a user's own home (USR.11) overrides it.
  Owns the one server-side derivation of the regions of created sectors,
  exposed as `GET /api/galaxy/regions` (sectors of each region, bounding
  box, outline edge set) with one contiguity rule (face adjacency; edge
  or corner contact is not contiguous, the default recorded on UX.48);
  the Galaxy Map's "which region to open on" and "jump between regions"
  and UX.48 read it (2026-10-16 review).
  Prerequisite: MAP.189. Related: MAP.189, MAP.151 (done, PR #1539).

- [ ] **MAP.191 Galaxy Map default view: no density shading, only nebulae shade the stars behind them**
  GitHub issue [#1233](https://github.com/dwhagar/planetGen/issues/1233)
  (Boss, 2026-10-14 00:17Z): "Default view should in the galaxy map
  should be no density shading, only shading should be that of nebulae
  that should change the colors of the stars shining through the nebulae
  and the boundaries should be fussy and translucent... going to the
  infinite zoom model we don't need to see the densities and anything
  but the stars and markers for where on the map known space is." Done:
  the default Galaxy Map shows stars and markers for known space with no
  density, age or luminosity fill; the "Color by" density and age views
  stay as options; nebulae tint the colour of the stars seen through
  them and have soft translucent boundaries (MAP.142). Owner:
  Foundations lane 2, after the nebula map work (MAP.176; MAP.155 merged in #1253).
  Prerequisites: MAP.214, MAP.215. Related: MAP.131, MAP.142, MAP.176,
  MAP.146.
  Split (2026-10-16, Boss 1:07 PM EDT: look at each item for smaller
  chunks and shared work): this item is now the umbrella and is done
  when its subitems are: MAP.214, MAP.215. Its Done has two halves of
  different versions.

- [ ] **MAP.192 See the objects whose positions changed in the last orbital update**
  GitHub issue [#1231](https://github.com/dwhagar/planetGen/issues/1231)
  (Boss, 2026-10-13 16:41Z): "I would like to be able to see the
  positions that updated in the last orbital update." Decided (Boss,
  2026-10-14 13:25Z): an admin page lists the objects that moved in the
  last update, with a link to view them on the Galaxy or Sector Map; the
  list stays until the next update runs; each entry shows how far the
  object moved and in what direction, as a vector; a "ghost" of the
  previous position can be shown on the map on request. Done: that page,
  the vector and the ghost, from the update's own record.
  Split (2026-10-16, Boss 1:07 PM EDT: look at each item for smaller
  chunks and shared work): this item is now the umbrella and is done
  when its subitems are: DB.57, MAP.216. Nothing records what a run
  moved today (orbit_simulation_state holds one timestamp), so this item
  is a new table written inside the advance functions plus an admin page
  plus a map layer.
  Prerequisites: DB.57, MAP.216. Related: GEN.105, GEN.139, ADM.36.

- [ ] **MAP.193 Show every object's axis of rotation, direction of rotation and axial tilt**
  GitHub issue [#1227](https://github.com/dwhagar/planetGen/issues/1227)
  (Boss, 2026-10-12 23:52Z): "For all objects show their axis of
  rotation, direction of rotation, and the axial tilt. Also should be
  able to see the vectors of motion, but I think that's in another item
  already." Done: the System Map and the object pages draw the spin
  axis, the direction of rotation and the axial tilt for every star,
  planet, moon and compact object, each toggled in the layer menu; the
  motion vectors are drawn by the same layer menu (no existing item
  draws them: they come with the orbital-update work, GEN.105). Owner:
  Foundations lane 2, after the fly-through map.
  Scope after the review (2026-10-16 review): the spin-axis, direction
  and tilt drawing and the motion-vector drawing; the toggles are
  registered layers.
  Prerequisites: MAP.217. Related: GEN.105, GEN.139, ADM.36, MAP.169.

- [ ] **MAP.194 Nebulae and globular clusters open their own zoom window from the Galaxy Map**
  GitHub issue [#1204](https://github.com/dwhagar/planetGen/issues/1204)
  (Boss, 2026-10-10 23:34Z): "Groupings of items like this should get
  their own zoom window from the galactic view that shows the specific
  arrangements. Neither care how many sectors they go through." Done:
  choosing a nebula or a globular cluster on the Galaxy Map flies into a
  view of that object alone, showing its own arrangement (cloud shape,
  cluster stars) whatever sectors it crosses. Owner: Foundations lane 2,
  after the nebula map work (MAP.176 to MAP.180).
  Prerequisites: MAP.222, MAP.223. Related: MAP.176, MAP.178, MAP.180,
  GEN.209.
  Split (2026-10-16, Boss 1:07 PM EDT: look at each item for smaller
  chunks and shared work): this item is now the umbrella and is done
  when its subitems are: MAP.222, MAP.223. MAP.194 covers two kinds of
  object, nebulae and globular clusters, but only the nebula half is
  blocked by MAP.176 (in its own thread).
  Scope after the review (2026-10-16 review): the nebula window (child
  a) and the cluster window (child b, in phase 2, which needs the parent
  link of GEN.209).

- [ ] **MAP.195 Galaxy Map finer tiles leave out the stars a coarser tile already sent**
  Foundations lane 2 (2026-10-14): MAP.159 sends packed tiles without
  leaving out the stars a parent tile already sent, because the page
  would need the parent tile in hand for every child. Done: finer tiles
  leave those stars out (31 to 43% of the bright records at levels 3 to
  6, all of them at levels 7 to 12) and the page keeps parents available
  so a child never draws without them, with a test that a tile and its
  parents together hold exactly the stars the full list holds. Measure
  the saving against the cost of holding the parents before building it.
  Status (2026-10-09): MAP.159 merged (#1305, Foundations lane 2): tiles
  are sent packed at about 20 bytes a star, kept packed in the server
  cache and the browser's IndexedDB. Leaving out stars a parent tile
  already sent is NOT built; that is this item.
  Builds only if the measurement says so; part of the tile-format chain
  MAP.171, MAP.195, MAP.173 (one cache-stamp bump) (2026-10-16 review).
  Prerequisite: none. Related: MAP.159 (done, #1305).

- [ ] **MAP.207 The Galaxy Map draws satellites, clusters and wanderers as points beyond the disk outline**
  Step 9 of the whole-sphere train. Components (the Clouds, Sagittarius,
  dwarfs, clusters) and wanderers (hypervelocity stars, kicked neutron
  stars) are drawn as points where the outline ends, with no change to
  the shading of the sectors. Map code uses `galaxy/spherical.py` and
  `galaxy/frame.py`, never `geometry.py`. Done: a visitor flying out of
  the disk finds the Magellanic Clouds where the data puts them, and
  nothing is drawn for the empty halo cells. Owner: Interface. Nebulae hosted by pre-placed stars show through the nebula
  items (MAP.172 to MAP.181) from their scatter rows (GEN.241).
  Prerequisites: GEN.238, MAP.217.
  Design: [docs/design/intergalactic-and-halo-contents.md](design/intergalactic-and-halo-contents.md)

- [ ] **MAP.209 Opening a sector from the Galaxy Map flashes one view, reloads into another, and leaves the details panel empty (bug)**
  Boss (2026-10-15 9:42 PM EDT, on his own server): "there's some kind
  of bug transitioning to sector view from galactic view, for a moment
  it shows one view, then reloads to a different view, but does not show
  any details about what's on the screen to the right of the viewing
  window." Reproduce first: fly from the Galaxy Map into a filled sector
  (galaxystageview.js decides when a sector is open and flies the
  camera, galaxysector.js loads the sector's scene and joins its picker
  layers to the map's), watch for the reload (a full page load such as
  openSectorPage in galaxymap3d.js, or the scene rebuilt twice), and check whether it
  is the same defect as MAP.210. Find the root cause, fix it, and add a
  browser test that opens a sector from the map and asserts one stable
  view and a filled details panel (a star's name and facts) once the
  sector is open. Check whether TEST.163 (the  count is 0 on the nav
  route page, red on main) shares the cause. Done: the move from the
  galactic view into a sector shows one continuous view with no reload,
  and the panel to the right of the viewer shows the opened sector's
  details. Owner: Core 2, next after its running item, with MAP.208 (the
  same map and scene code, so one lane avoids conflicts; Boss,
  2026-10-15 9:42 PM EDT, via the coordinator).
  Status (Core 2, 2026-10-16, 2:21 AM EDT): not reproduced. On a
  generated-database site in headless Chromium with reduced motion off
  (so the camera flight plays; every older map browser test ran with
  reduced motion), opening a charted sector by address, by clicking its
  cell and by a wheel run from the whole galaxy each gave one document
  (no reload), one URL, the panel filled with the sector's details, a
  star click selecting it, and no page or console errors; requests
  peaked at 3 a second, so the page limit is not the cause. Tests only
  merged (#1776, src/tests/test_web_browser_sector_open.py); the item
  stays open. Reading the code ruled out a polling reload (the stamp is
  used only on tile fetch) and openSectorPage runs only on a pinned
  sector page. Needed from Boss: browser and device, whether he used the
  wheel, a click or a link, the sector name, any red text in the browser
  console (F12) when it happens, and whether it still happens on the
  rebuilt database (his ran code from before migration 0093 and
  MAP.208).
  Prerequisites: none. Related: MAP.210, MAP.66, MAP.208, TEST.163.

- [ ] **MAP.210 Nothing on the map can be selected, even in sector view (bug)**
  Boss (2026-10-15 9:42 PM EDT, on his own server): "I also cannot
  select anything on the map even in sector view, when I should be able
  to select stars within the sector." Stars inside an opened sector
  should be pickable (a click or hover goes through the map's picker,
  mappick.js, and fills the info panel, showInfo in galaxysector.js) and so
  should the other pickable objects. This may be the cause of the empty
  details panel in MAP.209, so one root cause is possible; reproduce
  first, find the root cause, and fix it. Add a browser test that opens
  a sector, clicks a star's projected position and checks that the
  selection ring and the details panel fill, and the same for a star in
  the Galaxy Map's own tiles. Done: in sector view and on the Galaxy Map
  a click on a star selects it, the details panel shows its facts, and
  the browser test catches a regression. Owner: Core 2, with MAP.209
  (Boss, 2026-10-15 9:42 PM EDT, via the coordinator).
  Status (Core 2, 2026-10-16, 2:21 AM EDT): not reproduced, see
  MAP.209's status note (same sector-open tests, #1776, found one
  document, a filled panel and working star selection); the item stays
  open and waits on Boss's answers listed there.
  Prerequisites: none. Related: MAP.209, MAP.66, MAP.208, TEST.163.

- [ ] **MAP.212 A region selection tool on the Galaxy Map and the Sector Map, for admin actions**
  Boss (2026-10-16, 1:07 PM EDT), headline of 8.6 (ADM.57): the admin
  has complete control; the world adjusts to what the admin wants, after
  the admin approves the generation plan. An admin draws or types a
  region (sphere by centre and radius, box, shell, band, wedge, or the
  sectors in view) on the Galaxy Map and the Sector Map, sees it
  outlined with a count of the systems and phenomena inside, and passes
  it to an admin action (a density edit, a regeneration). Done: the tool
  on both maps (browser test), the count from the region object, and the
  hand-off to the actions of ADM.73 and ADM.80.
  Unchanged in scope; it passes a GEN.302 region and shows the count,
  and the outline of the sectors in view comes from Span (2026-10-16
  review).
  Prerequisites: GEN.302, MAP.217. Related: ADM.57, NAV.45.
  Design: [docs/design/edit-everything.md](design/edit-everything.md)

- [ ] **MAP.213 Scene, maps and the JavaScript twin for trees: per-node zoom, log-radial overview and positions_at in JavaScript**
  Split from GEN.129 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). The System Map and scene draw two
  stars at most. Draw a tree: the scene places every member from the
  nodes, the map zooms per node and has a log-radial overview so a
  10,000 AU outer pair and a 0.1 AU inner pair both read, and static/
  gets the JavaScript twin of positions_at, held equal by the existing
  twin tests. Done: a triple and a seven-star system draw with every
  member at the right place and zoom level, the JavaScript and Python
  positions agree in the twin test, and the page lists the tree.
  Owner: Core 3 (Boss, 2026-10-16, 2:42 PM EDT: Core 3 takes the 8.1 new
  generate scatters line after GEN.280; the multi-star train follows its
  spine, 3:00 PM EDT). This replaces any earlier Core or Core 1 owner
  line on this item.
  Owner: Core 4 (handed over from Core 3 by Boss's "Split it", 2026-10-16,
  4:49 PM EDT; Core 3 agreed 4:53 PM EDT; waits on GEN.311 (Core 3, after DB.53)). Core 4 asks Core 3
  before touching stellar evolution constants in tuning.py. A line here
  that named Core 3 or Core 4 earlier is replaced by this one.
  Prerequisites: GEN.311. Related: GEN.129.

- [ ] **MAP.214 Galaxy Map default view: no density, age or luminosity shading, only stars and markers for known space**
  Split from MAP.191 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). Boss (GitHub issue #1233): going
  to the infinite zoom model the default view needs only the stars and
  the markers for where known space is. The default Galaxy Map shows
  stars and markers with no density, age or luminosity fill; the Color
  by density and age views stay as options. Done: the default view has
  no shading fill, Color by density and age still work and are off by
  default, and a browser test checks the default.
  Prerequisites: none. Related: MAP.191.

- [ ] **MAP.215 Nebulae tint the stars seen through them, with soft translucent boundaries**
  Split from MAP.191 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). In the default Galaxy Map view
  nebulae are the only shading: they change the colour of the stars
  shining through them and have fuzzy, translucent boundaries (MAP.142),
  read from the nebula layer of the tile. Done: a star behind a nebula
  is tinted by the nebula's family colour, the boundary fades, and the
  density views are unchanged.
  Prerequisites: MAP.176, MAP.142. Related: MAP.191.

- [ ] **MAP.216 Admin page of the objects that moved in the last update, with the movement vector and a map ghost of the old position**
  Split from MAP.192 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). An admin page lists the objects
  that moved in the last update with a link to each on the Galaxy or
  Sector Map; each entry shows how far the object moved and in what
  direction as a vector; the ghost of the previous position is drawn on
  the map on request. Done: the page, the vector and the ghost, read
  from child a's record; the list stays until the next update runs.
  Scope after the review (2026-10-16 review): the admin page of moved
  objects and the ghost data; the "show ghosts" toggle is a registered
  layer.
  Prerequisites: DB.57, MAP.217. Related: MAP.192.

- [ ] **MAP.217 Map layer interface: one registry the Galaxy stage view, the Sector Map and the System Map share for a toggled layer (menu entry, modes, URL key, legend, picker layer, fade rule)**
  Found in the split review (2026-10-16, Boss 1:07 PM EDT: look at each
  item for smaller chunks and shared work, 1:33 PM EDT: shared
  foundations). Today every map wires its own layer controls (MAP.123's
  kind toggles are written per map in galaxysector.js, sectorscene.js
  and systemmap.js; mappick.js gives only the picker). Add
  `static/maplayers.js`: a layer declares `{id, label, modes, urlKey,
  legend, pickerLayer, fadedWithSectors, draw(scene)}`; the Menu, the
  URL state, the legend block and the picker registration are produced
  from the registry on all three maps, and a layer says whether it is
  hidden by the blocker fade (MAP.149 `nearfield.js`; a course layer is
  not, per MAP.121's note). Existing toggles move onto it only when
  their map is touched, no sweep. Done: the module, one real layer
  ported on each of the three maps (the existing per-kind toggle), the
  URL round trip, and a browser test that a registered layer appears in
  each Menu, keeps its URL key and draws its legend.
  Prerequisites: none. Related: MAP.168, MAP.169, MAP.192, MAP.193,
  MAP.212.
  The label engine (map labels) registers as a layer kind when this
  registry lands (2026-10-16).
  Owner: Core 2 (Boss, 2026-10-16, 3:00 PM EDT: a global top blocker, 40
  items wait on it), after the 8.1 phase 1 top 5.

- [ ] **MAP.218 Galaxy Map density control: a readout for everyone, an editor for admins, with preview, confirm and undo**
  Split from MAP.119 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). On the Galaxy Map info panel
  guests and users see a static readout ('Expected density: 6.3 systems
  per 4 pc sector (x1.0)'); an admin also gets a number field, a
  log-scale range, 'Reset to 1.0' and quick steps (x0.5, x2), a live
  preview (the map re-colours with the draft multiplier; 'Expected
  systems in view', 'N filled sectors unaffected'), a confirm dialog
  listing scope, old and new value and counts, and an undo toast with a
  history list. The scope is the slab or block in view, passed to child
  b as a region. Done: the panel for both audiences, the admin edit
  calling `set_override`, and a browser test that edits, previews,
  confirms and undoes on a seeded galaxy.
  Prerequisites: ADM.94. Related: MAP.119.

- [ ] **MAP.219 Galaxy Map marker layer: a kind registry with per-kind toggles, the zoom budget and picking, shipping the black-hole kind**
  Split from MAP.132 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). Build the one layer the other
  marker kinds draw into: a registry (kind, glyph, colour token, minimum
  zoom, source), the per-kind toggles of MAP.123, the budget of MAP.116
  (markers show only at the zooms where it allows), a 2 px floor under
  which a sprite stays hidden, hover and pick with the same keyboard
  path as other map objects, and the legend rows. Ship it with the
  black-hole kind from the stored black holes. Done: the layer, the
  registry with black holes as its first kind, toggle and budget
  behaviour, and a browser test that finds a black-hole marker on a
  seeded sector and hides it when the kind is switched off.
  Prerequisites: none. Related: MAP.132.

- [ ] **MAP.220 Markers for habitable worlds and chosen classes, optional per kind**
  Split from MAP.132 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). Add the kinds that issue #1234
  asks for: habitable worlds (by the habitability index, GEN.84) and
  chosen object classes, each optional through its toggle and shown at
  the zooms the budget allows. Done: both kinds in the registry of child
  a, and a browser test that finds one of each on a seeded sector.
  Prerequisites: MAP.219. Related: MAP.132.

- [ ] **MAP.221 Station and vessel icons on the Galaxy Map**
  Split from MAP.132 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). Register the station and vessel
  kinds in the marker layer once user objects and scattered vessels
  exist (issue #1234: 'icons for space stations, etc.'). Done: the two
  kinds in the registry with their toggles, and a browser test that
  finds a user's ship and a scattered station.
  Prerequisites: MAP.219, USR.29, POP.11. Related: MAP.132.

- [ ] **MAP.222 Nebula zoom window: choosing a nebula on the Galaxy Map flies into a view of that nebula alone**
  Split from MAP.194 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). Choosing a nebula (a region of
  MAP.176, a cloud of MAP.178) flies into a view that shows that nebula
  only: its shape and member clouds or hosted stars, centred and framed
  on its bounds, whatever sectors it crosses, with a way back to the
  Galaxy Map. Build the view as a generic 'group window' (a bounding
  volume, a member list, a draw function per kind) so a second kind can
  plug in. Done: the window for nebulae, entered from the map and from
  the keyboard, and a browser test that opens it on a nebula spanning
  two sectors.
  Prerequisites: MAP.176. Related: MAP.194.

- [ ] **MAP.223 Globular cluster zoom window: choosing a cluster flies into a view of its own stars**
  Split from MAP.194 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). Plug the cluster kind into the
  group window of child a: the members are the stars whose parent is the
  cluster (GEN.209), drawn from the King profile where a sector is not
  yet filled. Done: choosing a cluster on the Galaxy Map opens its
  window, and a browser test checks the window on a seeded cluster that
  crosses two sectors.
  Prerequisites: MAP.222, MAP.224. Related: MAP.194.

- [ ] **MAP.224 Cluster pages and the Galaxy Map use the parent link**
  Split from GEN.209 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). The cluster page lists the stars
  that point to it (and 'not yet charted' for the unfilled part), a
  star's page names its cluster, and the Galaxy Map and the Sector Map
  label a cluster by its name, or by its designation while it has none.
  Done: the page and the labels, and a browser test on a seeded cluster
  before and after its first sector is filled.
  Prerequisites: GEN.359. Related: GEN.209.

- [ ] **MAP.225 The wheel zooms toward the point under the cursor, with clearance**
  Found 2026-10-16 (Boss, 2:29 PM EDT: "Do we need to add the lacking
  items to the TODO?", after Core 2 reported what MAP.146 still lacks
  once MAP.152 merged): this piece was only in the Done line of MAP.146
  and in fly-through-view-distance.md section 4.3, with no item of its
  own. Boss (2026-10-09 22:24Z): "have them based on the center of where
  the cursor is clicked. So that we always are drilling down exactly
  where the user wants." Today the wheel zooms about the shared centre
  and Shift-click chooses it (MAP.199). Done: the wheel changes the
  distance by `exp(-0.0025 * deltaY)` toward the point under the cursor,
  which stays fixed on screen; that point is, in order, the picked star,
  remnant or cell under the cursor, else the first sample along the
  cursor ray where the integrated column density passes a threshold,
  else the point at the current focus depth; each step moves by
  `min(requested, 0.5 * (d_nearest - d_clear))`, so the speed is
  logarithmic in the distance to the nearest pickable thing and the
  camera never lands on a star or inside a wall by accident. Shift-click
  stays as "turn about this point". A browser test zooms toward a chosen
  star and checks it stays under the cursor.
  Owner: Core 2 (Boss, 2026-10-16, 2:43 PM EDT: after the 8.1 phase 1
  top 5, Core 2 takes the infinite zoom tree: MAP.146 and its subitems,
  then the label items).
  From Core 2 (2026-10-16, 4:03 PM EDT, MAP.226 merged as #1859): a
  double click on empty space or a block stays two clicks; the "half way
  toward the point under the cursor" double click is not built and
  belongs here. A second click within 350 ms and 8 px is now a double
  click, so clicking a selected object again to clear it needs a pause.
  Prerequisites: MAP.152. Related: MAP.146, MAP.199, MAP.125.

- [ ] **MAP.227 The camera names its container: the header and the breadcrumb come from its position**
  Found 2026-10-16 (Boss, 2:29 PM EDT: "Do we need to add the lacking
  items to the TODO?", after Core 2 reported what MAP.146 still lacks
  once MAP.152 merged): this piece was only in the Done line of MAP.146
  and in fly-through-view-distance.md section 4.3, with no item of its
  own. Done: "I am in sector X" is worked out from the camera position
  (`geometry.sector_address_at` is closed form) and shown in the header,
  and the breadcrumb becomes a trail of containers (galaxy, region,
  block, sector, system) derived from the position, not from stored
  state; clicking a crumb flies out to that container. A browser test
  flies into a sector and checks the header and every crumb.
  Owner: Core 2 (Boss, 2026-10-16, 2:43 PM EDT: after the 8.1 phase 1
  top 5, Core 2 takes the infinite zoom tree: MAP.146 and its subitems,
  then the label items).
  Prerequisites: MAP.226. Related: MAP.146, MAP.139.

- [ ] **MAP.228 The camera in the address: the URL and bookmarks hold the camera's target, distance and orientation**
  Found 2026-10-16 (Boss, 2:29 PM EDT: "Do we need to add the lacking
  items to the TODO?", after Core 2 reported what MAP.146 still lacks
  once MAP.152 merged): this piece was only in the Done line of MAP.146
  and in fly-through-view-distance.md section 4.3, with no item of its
  own. Done: the Galaxy Map URL carries the camera (target, distance,
  orientation) and opening it puts the camera there; bookmarks and the
  charted-region bookmarks (UX.45, UX.48) store the same value; old
  block URLs and stage tokens stop working (no backward compatibility,
  Boss's rule). A browser test copies the address, reloads and checks
  the camera is where it was.
  Owner: Core 2 (Boss, 2026-10-16, 2:43 PM EDT: after the 8.1 phase 1
  top 5, Core 2 takes the infinite zoom tree: MAP.146 and its subitems,
  then the label items).
  Prerequisites: MAP.225. Related: MAP.146, UX.45, UX.48, MAP.189.

- [ ] **MAP.229 Map label engine: names placed by priority, with no overlaps, shown and hidden by zoom**
  Boss (2026-10-16, 2:30 PM EDT): "Oh do we have map labels anywhere in
  our TODO lineup?" Checked the open and finished items, the design
  notes and the code: no. Only three small label pieces exist: the
  compass and neighbour sprite labels of `sectorscene.js`, the course
  markers' text sprites in `galaxymap3d.js`, and the HTML labels laid
  over the system view (`systemview3d.js`); the design notes carry no
  label placement, collision or zoom rule. Done: `static/maplabels.js`,
  one module on the map engine: a label is text, a 3D anchor, a kind, a
  priority and a distance or pixel-size range; on every camera move the
  anchors are projected, those outside the view or below the size range
  are dropped, and the rest are placed by priority on a screen grid with
  a padding so none overlap, up to a budget (default 60), with a short
  fade (none under reduced motion) and hysteresis so a label does not
  flicker at a range edge. Text is drawn with a halo so it reads on both
  themes at 3:1 contrast. The Menu gets a "Labels" entry with one toggle
  per kind, kept in the URL. Tests: unit tests of the placement (no
  overlap, priority wins, stable between frames) and a browser test that
  labels appear and disappear with zoom; the timing with 20,000 stars in
  view is reported.
  Owner: Core 2 (Boss, 2026-10-16, 2:43 PM EDT: after the 8.1 phase 1
  top 5, Core 2 takes the infinite zoom tree: MAP.146 and its subitems,
  then the label items).
  Prerequisites: MAP.152. Related: MAP.146, MAP.116, MAP.123, MAP.217.

- [ ] **MAP.230 Names of stars, systems and sectors on the Galaxy Map**
  Boss (2026-10-16, 2:30 PM EDT): "Oh do we have map labels anywhere in
  our TODO lineup?" Checked the open and finished items, the design
  notes and the code: no. Only three small label pieces exist: the
  compass and neighbour sprite labels of `sectorscene.js`, the course
  markers' text sprites in `galaxymap3d.js`, and the HTML labels laid
  over the system view (`systemview3d.js`); the design notes carry no
  label placement, collision or zoom rule. Done: the brightest stars and
  systems in view show their names and the sector in view shows its name
  or address, ranked by apparent brightness and nearness to the camera,
  through the label engine; a star without a name shows its designation
  only when zoomed in far enough; a kind toggle for each. A browser test
  checks a named star is labelled at the zoom where it is among the
  brightest in view and that the label goes when it is not.
  Owner: Core 2 (Boss, 2026-10-16, 2:43 PM EDT: after the 8.1 phase 1
  top 5, Core 2 takes the infinite zoom tree: MAP.146 and its subitems,
  then the label items).
  Prerequisites: MAP.229. Related: MAP.146, MAP.203.

- [ ] **MAP.231 The Sector Map and the System Map use the label engine**
  Boss (2026-10-16, 2:30 PM EDT): "Oh do we have map labels anywhere in
  our TODO lineup?" Checked the open and finished items, the design
  notes and the code: no. Only three small label pieces exist: the
  compass and neighbour sprite labels of `sectorscene.js`, the course
  markers' text sprites in `galaxymap3d.js`, and the HTML labels laid
  over the system view (`systemview3d.js`); the design notes carry no
  label placement, collision or zoom rule. Done: the compass and
  neighbour labels of the Sector Map and the body labels of the System
  Map are placed by the same engine (collision, priority, zoom), with
  the same Labels toggle and the same halo; no second label code is
  left. Browser tests on both maps keep the labels they show today.
  Owner: Core 2 (Boss, 2026-10-16, 2:43 PM EDT: after the 8.1 phase 1
  top 5, Core 2 takes the infinite zoom tree: MAP.146 and its subitems,
  then the label items).
  Prerequisites: MAP.229. Related: MAP.146, MAP.229.

- [ ] **MAP.232 Names for the big features: spiral arms, bar, bulge, halo, satellites and clusters**
  Boss (2026-10-16, 2:30 PM EDT): "Oh do we have map labels anywhere in
  our TODO lineup?" Checked the open and finished items, the design
  notes and the code: no. Only three small label pieces exist: the
  compass and neighbour sprite labels of `sectorscene.js`, the course
  markers' text sprites in `galaxymap3d.js`, and the HTML labels laid
  over the system view (`systemview3d.js`); the design notes carry no
  label placement, collision or zoom rule. Done: through the label
  engine, the galaxy's large features show their names at the zooms
  where they are readable: the arms, the bar, the bulge, the halo and
  disc, the satellites and wanderers (MAP.207) and the clusters and
  dwarfs of `galaxy_components` (GEN.237), with the names from the
  galaxy model and the component table (a small data table where the
  model has none), each a kind with its own toggle. A browser test
  checks one of each kind.
  Prerequisites: MAP.229, MAP.207. Related: GEN.237, GEN.209, MAP.145.

- [ ] **MAP.233 Names for nebulae on the Galaxy Map**
  Boss (2026-10-16, 2:30 PM EDT): "Oh do we have map labels anywhere in
  our TODO lineup?" Checked the open and finished items, the design
  notes and the code: no. Only three small label pieces exist: the
  compass and neighbour sprite labels of `sectorscene.js`, the course
  markers' text sprites in `galaxymap3d.js`, and the HTML labels laid
  over the system view (`systemview3d.js`); the design notes carry no
  label placement, collision or zoom rule. Done: nebula regions and
  singles show their names (or class and designation) through the label
  engine at the zooms where the region is large enough to read, with a
  toggle; the label does not hide the nebula's own pick and hover. A
  browser test checks a seeded nebula is labelled and unlabelled with
  zoom.
  Prerequisites: MAP.229, MAP.176. Related: MAP.176, MAP.180.

- [ ] **MAP.234 Galaxy Map key: a Key section in the help dialog and an on-map key that names everything the map draws**
  Boss (2026-10-16, 7:16 PM EDT, timeline): "ON the galaxy map, I just
  realized there is no key, but what are the greenish colored spheroid
  shapes? They seem to be 2 intersecting elipses to define a space. What
  are those?" Found by the read-only "galaxy map shapes" thread: they
  are superbubbles (GEN.254, merged #1814), each drawn as three faint
  great circles (`static/galaxymap3d.js`, "Superbubbles (GEN.254)",
  colour 0x6fb3d8, opacity 0.14, always on, never picked, no hover, no
  legend entry). The help dialog (`GALAXY_MAP_HELP` in
  `src/planetgen/web/maps/galaxymap3d.py`, "Reading the map") covers
  only density blocks and stars; it does not mention superbubbles,
  nebula clouds, supernova remnants, territory balls, course lines or
  the made-sector rings. The only on-map legends today are the Color by
  gradient bar (MAP.131) and the course key (NAV.20). Done: a Key
  section in the help dialog and a small on-map key (collapsible) naming
  every drawn thing with its colour and shape: density blocks, stars,
  nebula clouds, remnants, superbubbles, territories, the course line
  and the made-sector rings; each entry appears only when that thing is
  drawn at the current view; a browser test that every kind drawn on the
  map has a key entry (so a new layer cannot be added without one). Not
  covered elsewhere: DOC.7 is the help page (it lists "nebula and
  territory overlays" but not superbubbles) and MAP.217 gives a legend
  per toggled layer, but superbubbles are not a toggled layer today. The
  superbubble toggle, hover name and radius are MAP.235. Owner: Core 2
  (the Galaxy Map and its help dialog).
  Prerequisites: none. Related: MAP.131, NAV.20, MAP.217, DOC.7,
  GEN.254.

- [ ] **MAP.235 Superbubbles as a registered Galaxy Map layer: a toggle, a hover name and radius, and a key row**
  Split from MAP.234 (2026-10-16, 7:21 PM EDT). Superbubbles (GEN.254)
  are drawn always on as faint great circles with no toggle, no picking
  and no hover. Done: they are a layer on the MAP.217 registry (menu
  entry, URL key, legend row from the registry), on by default, off on
  request; hovering one names it and gives its radius and centre; it is
  a picker layer so a click selects it. Needs the registry (MAP.217).
  The key text itself is MAP.234. Owner: Core 2.
  Prerequisites: MAP.217. Related: MAP.234, MAP.217, GEN.254.

- [ ] **MAP.236 Superbubbles need clearer bounds and slightly more transparency on the Galaxy Map**
  Boss (2026-10-16, 7:24 PM EDT, timeline): "the super bubbles need
  better definition to the bounds and be slightly more transparent."
  Today each superbubble (GEN.254) is three faint great circles
  (`static/galaxymap3d.js`, "Superbubbles (GEN.254)", colour 0x6fb3d8,
  opacity 0.14, always on), so the bound reads as two crossing ellipses
  and the inside has no edge. Done: the bound of each bubble is clearly
  defined (a brighter rim where the shell meets the line of sight, or a
  crisp outline of the true sphere silhouette, not three crossing
  rings), and the fill is slightly more transparent than today; both
  read in the dark and the light theme and at every zoom stage, and the
  look is drawn the same on the Sector Map where a bubble is shown. A
  browser screenshot test checks that the rim is brighter than the fill
  and the fill alpha is lower than 0.14. The key entry and toggle are
  MAP.234 and MAP.235; this item changes only the look.
  Owner: Core 2 (the Galaxy Map; proposed, Core 2 reads the plan).
  Prerequisites: none. Related: MAP.234, MAP.235, GEN.254.

- [ ] **MAP.237 The Galaxy Map's "Charted only" view is invisible again: charted space cannot be told from the rest (bug)**
  Boss (2026-10-16, 7:24 PM EDT, timeline): "the 'charted only' has gone
  back to being invisible from the galactic map, so we need some kind of
  way to see where they are." MAP.111 made "Charted only" dim blocks
  with nothing generated, dim stars outside charted sectors and outline
  each charted block ("Charted only" in the Menu; `setChartedOnly` in
  `static/galaxymap3d.js` sets the `chartedOnly` star uniform and calls
  `stageView.setChartedOnly`), and MAP.112 let empty blocks be hovered
  with it on. It was visible then and is not now, so something since
  then broke it or hid it (the infinite-zoom and tile work changed the
  stages). Reproduce on a database with a charted region, find the root
  cause, and fix it; do not add a second way to see charted space (the
  look itself is MAP.237). Done: with Charted only on, charted space is
  visibly different from the rest at every stage (whole galaxy, a window
  of blocks, a single sector window), with a browser test that compares
  a screenshot (or the drawn colours) with the toggle on and off at each
  stage and fails if they are the same.
  Owner: Core 2 (the Galaxy Map; proposed, Core 2 reads the plan).
  Prerequisites: none. Related: MAP.111, MAP.112, MAP.238, MAP.190,
  UX.48.

- [ ] **MAP.238 Charted regions shown by colour and opacity, not outlined blocks: drop the outlines everywhere charted space is marked**
  Boss (2026-10-16, 7:24 PM EDT, timeline): "I also don't like the
  outlined blocks for it, they are really ugly and make it hard to
  navigate and see where you are. So rather I want to highlight with
  color and opacity." Today charted space is marked by an outline around
  each charted block (MAP.111) and UX.48 plans an outlined region shape.
  Done: charted blocks and regions are highlighted by a fill colour and
  opacity that rise with the filled share (the amber of generated blocks
  already grows more solid with it), with no outline drawn around
  charted blocks or regions on the Galaxy Map or the Sector Map;
  "Charted only" (the Menu) raises the highlight and dims the rest; the
  hover outline on a choice stays (it marks what the pointer is on, not
  charted space); the region highlight of UX.48 uses the same fill. Both
  themes, every stage, and a browser test that no outline mesh is built
  for charted blocks and the fill differs between a charted and an empty
  block. Needs MAP.237 (the toggle must show something at all before it
  is restyled), and changes the look MAP.111 gave it.
  Owner: Core 2 (the Galaxy Map; proposed, Core 2 reads the plan).
  Prerequisites: MAP.237. Related: MAP.237, MAP.111, MAP.190, UX.48.

- [ ] **MAP.239 Raise the Galaxy Map tile root to the 262 kpc reach: `TILE_ROOT` in `galaxy/viewport.py` and `galaxymap3d.js`**
  From the GEN.372 precision audit (Core 1,
  docs/design/numeric-precision-audit.md section 5, 2026-10-16, 7:39 PM
  EDT), filed by the TODO thread; verdict: 64-bit is enough for every
  stored value. Boss decided by card (2026-10-16, 7:36 to 7:38 PM EDT):
  the map tile root rises to the 262 kpc reach. Done: `TILE_ROOT` in
  `src/planetgen/galaxy/viewport.py` and
  `src/html/static/galaxymap3d.js` is the reach, so a tile key covers
  every address in it; stale comments fixed (`web/lib/tilepack.py:30`,
  `galaxyprisms.js:475`); `test_tiles_cover_the_reach` checks that every
  cell at the reach edge falls in a tile. Owner: Core 2.
  Prerequisites: none. Related: GEN.372, MAP.146.

- [ ] **MAP.240 Browser buffers hold camera-relative coordinates: territory dots, superbubbles, made rings, block prisms, stage views, system view dots and belts**
  From the GEN.372 precision audit (Core 1,
  docs/design/numeric-precision-audit.md section 5, 2026-10-16, 7:39 PM
  EDT), filed by the TODO thread; verdict: 64-bit is enough for every
  stored value. Boss decided by card: fix camera-relative now. Float32
  attribute buffers in the browser hold absolute coordinates (float32
  spacing at 262 kpc is far above a parsec), so far from the origin the
  shapes jitter and snap. Done: the buffers of territory dots,
  superbubbles, made-sector rings, block prisms, stage views, and the
  system view's dots and belts are built relative to the camera (or a
  per-object origin) in double precision on the CPU and the model matrix
  carries the offset; a browser test reads each attribute range and
  fails if any value exceeds the safe float32 range (a set bound, in
  parsecs, stated in the test) at a camera 262 kpc from the origin.
  Owner: Core 2. Related: the superbubble look (MAP.236) and the key
  (MAP.234) touch the same buffers.
  Prerequisites: none. Related: GEN.372, MAP.236, MAP.234, MAP.146.

- [ ] **MAP.241 Stars on the star maps should glow brighter, not grow bigger: size stays small, glow and colour carry luminosity and temperature (bug)**
  Boss (2026-10-16, 7:43 PM EDT, timeline): "Bug, stars on the starmap
  should glow brighter but not be bigger if that makes any sense so they
  are more obviously stars with glows and colors relative to the star's
  properties." Which map: the message says the star map; the default is
  both the Galaxy Map and the Sector Map (both draw stars the same way),
  and the Galaxy Map comes first. What exists: Boss's earlier rule
  (MAP.51/MAP.208, in the header of `static/galaxymap3d.js`) is "very
  small with a big glow", sizes by the star's radius, brightness by
  luminosity and colour by temperature. Today each star is a point
  `POINT_PX` across far away and, from about a sector away, grows into a
  core sized by the star's radius (`STAR_CORE_PX` 1.1 to 4.5 px) inside
  a halo whose width and strength both grow with luminosity
  (`STAR_MIN_PX` 6 to `STAR_MAX_PX` 30 px), so a luminous star becomes a
  bigger disc, not a brighter point; the Sector Map has its own curve
  (MAP.87, `static/starlight.js`). Done: a star's on-screen footprint
  stays small and close to the same for every star (the size no longer
  grows with luminosity or radius beyond a narrow band); luminosity
  raises the glow's intensity (a brighter, hotter-looking core and halo
  that stays tight), colour follows the blackbody temperature with
  enough saturation to read, and a dim star reads dim and a bright star
  reads bright at the same size; the halo-blend rule stays (crowded arms
  glow in their colours, not white); one shared star look for both maps
  (`static/starpoint.js` or `starlight.js`, not two curves); both
  themes; a browser test that two stars of very different luminosity
  have the same sprite size within a stated tolerance but different peak
  glow, and that colour follows temperature. The star key entry of
  MAP.234 says what size, glow and colour mean. Overlap checked: MAP.208
  (built the grow-into-a-core rule that this changes), MAP.87 (the
  Sector Map curve), MAP.164 and MAP.165 (colour and size, built),
  MAP.148 (visibility law, unchanged), MAP.217 (layer registry, no
  change), and the density/age coloring design, which is for blocks, not
  stars. Owner: Core 2 (proposed).
  Prerequisites: none. Related: MAP.208, MAP.87, MAP.164, MAP.165,
  MAP.234, MAP.148.

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
  between the two. Needs user accounts (USR.1) and their storage.
  Question (decided below): until user accounts exist, should saving be
  admin-only, or
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

  - [ ] **NAV.18 Saved courses list: open, rename, delete and the Direct/Route switch, over the saved-list store**
    Done: a Save course button on the NAV result; a Courses list
    (beside Bookmarks, from `bookmarks.js` or a sibling module) with
    open, rename and delete; the viewer switches between Direct and
    Route on a saved course.
    Owns the Courses list over the saved-list store (open, rename,
    delete) and the viewer's switch between Direct and Route; the Save
    action is NAV.4's; this adds the `course` kind to the store
    (2026-10-16 review).
    Prerequisites: UX.129.

  - [ ] **NAV.19 Saved courses in the account (after USR.7)**
    Done: a `user_courses` table in the control database, owned by an
    account; courses saved per browser can be imported into the account
    once; the API lists and edits a user's own courses only.
  Decided (2026-10-10): Boss (19:12Z): save courses per browser until
  accounts exist (as the default taken).
  Scope after the review (2026-10-16 review): the Save course action
  only: the button on the NAV result and on the map, the name prompt,
  and filling the record's forms (direct line, route, later the adjusted
  path); its list, open, rename, delete and switch text moved to NAV.18.
    The account side for the `course` kind of the saved-list store: the
    course payload and its recheck on open (NAV.17, NAV.39) in the
    shared `user_saved` table of USR.7, with no second table and no
    second import path (2026-10-16
    review).

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
  Prerequisites: NAV.21, NAV.22, NAV.23, MAP.217.
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
    Prerequisites: MAP.217.

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
  comparison. Question (decided below): a Hill sphere needs an orbit
  around a
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
    Scope after the review (2026-10-16 review): the obstacle query only;
    it returns each object's radius and basis from the keep-out table
    (NAV.54).
    Prerequisites: NAV.54.

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
    Waits on NAV.62: the study of how the gravity of the stars a cell
    leaves out is added to the course (Boss, 2026-10-16, 5:01 PM EDT).
    Prerequisites: NAV.55, NAV.56, NAV.62.

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
    Prerequisites: GEN.289.

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
  Decided (2026-10-10): Boss (19:16Z, Accept all): the answers recorded
  here stand: a lone star or object uses the galactic Hill radius;
  inside a system a star uses the radiation radius.
    Scope after the review (2026-10-16 review): routing around fields,
    the course test at the three scales, using the keep-out table
    (NAV.54) and the planner (NAV.26); the keep-out radius for fields is
    NAV.54's.
    Prerequisites: NAV.25, NAV.26, NAV.54.

- [ ] **NAV.36 Unknown-space jumps drawn red and glowing**
  Boss (2026-10-02 01:53Z): "a jump through unknown space is marked in
  red and glows to draw attention to it." Done: a hop NAV.12 flags as
  an unknown-space jump is drawn red with a glow in UX.35's horizontal
  route strip, on the NAV page's own map (`lib/navmap.py`) and on the
  Galaxy Map course (with NAV.20), with a legend entry; with
  `prefers-reduced-motion` it stays red without the pulse; it reads in
  both themes, and the route list also labels it in text so it isn't
  shown by color alone. Prerequisites: UX.35, MAP.217.
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
  unfilled cells. Decided (Boss, 2026-10-10 18:44Z, defaults approved;
  default yes): list those as
  "not generated yet" rows.
  Region query (stage C, research #1659): the nearby page reads by the
  cell key through NAV.60 and the radius ladder and early-stop count of
  NAV.61, so it stays fast in the dense core and shows "300+" from the
  capped count.
  Owner: Core 4, "Core 4 (2026-10-16)" (Boss approved the lane at 4:31
  PM EDT, 2026-10-16; queue in the lane report, queue position 7).

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

- [ ] **NAV.54 One keep-out table for every kind of object: fields, supermassive holes and quasars, moons, stars in their own system, nebulae and remnants**
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
  bounding sphere and the gas is thin. A test for each. Decided (Boss,
  2026-10-10 18:44Z, defaults approved; defaults taken): nebulae and
  remnants pass through with a
  note (chord length, column density, extinction) and no hard keep-out,
  with an optional soft cost later; a star inside its own system keeps
  out its radiation radius `max(10 R_star, sqrt(L / (4 pi F_lim)))`,
  F_lim = 50 kW/m^2 (0.17 AU for the Sun, 165 AU for an O3 star), hard
  and exempt when an end lies inside it; belts inside a system pass
  through with a note; the galactic centre keep-out of about 1.7 pc is
  exempt for any course to a star inside it.
  Found in the split review (2026-10-16, Boss 1:07 PM EDT: look at each
  item for smaller chunks and shared work, 1:33 PM EDT: shared
  foundations). Owns the keep-out radii and sits first in the Steering
  thread, with the knobs of NAV.55 before it: asteroid fields keep out
  `radius_ly`, a supermassive hole or quasar its sphere of influence G M
  / sigma^2, a moon its Hill radius, a star inside its own system its
  radiation radius max(10 R_star, sqrt(L / (4 pi F_lim))), nebulae and
  remnants pass through with a corrected note, and the galactic centre
  is exempt; F_lim and the floors read from tuning. Done as stated
  above, one test per kind; it no longer needs NAV.51.
  Prerequisites: NAV.55.
  Design: [docs/design/course-avoidance.md](design/course-avoidance.md)

- [ ] **NAV.55 A tuning block for the keep-out knobs**
  Put the keep-out knobs in `tuning.py`: F_lim = 50 kW/m^2, the star
  floor of 10 R_star, the below-resolution fraction (0.002 of the leg)
  and the polish tolerances.
  First in the Steering thread: the polish tolerances and the
  below-resolution fraction are read by the planner (NAV.26) and the
  F_lim and star floor by the radii table (NAV.54) (2026-10-16 review).
  Prerequisites: none.
  Design: [docs/design/course-avoidance.md](design/course-avoidance.md)

- [ ] **NAV.56 Census of overlapping keep-out spheres in a generated galaxy**
  Count how many pairs of keep-out spheres overlap in a generated
  galaxy. The answer decides whether NAV.26's cluster fallback ever
  fires.
  Needs the whole table (NAV.54) and runs before the planner (NAV.26),
  so the cluster fallback is built or skipped on the answer (2026-10-16
  review).
  Prerequisites: NAV.54.
  Design: [docs/design/course-avoidance.md](design/course-avoidance.md)

- [ ] **NAV.57 The Intergalactic Frame in navigation-frames.md and `navigation.py`**
  A frame type and the rule for choosing it (stage 2 of GEN.9).
  Prerequisite: GEN.157.
  Design: [docs/design/multiple-galaxies.md](design/multiple-galaxies.md)

- [ ] **NAV.58 Edit a course: add, remove, reorder and drag its waypoints, then save the edit**
  Boss (2026-10-10 21:05Z): the 8.5 headline is courses with waypoints
  that are editable and saved. Waypoints today are picked once and
  plotted (NAV.49). Done: a course drawn on the map (any display) can be
  edited in place: add a waypoint by picking an object, remove one,
  reorder them, or drag a waypoint to another object; the course is
  re-planned after each edit (bending around gravity wells once NAV.6 is
  built) and the edit is saved over the saved course or as a new one
  (NAV.4); a test checks add, remove, reorder and save.
  Prerequisites: NAV.49, NAV.4. Related: NAV.6, NAV.18, NAV.59.

- [ ] **NAV.59 Headline of 8.5: courses that bend around gravity wells, with waypoints, editable and saved courses**
  Boss (2026-10-10 21:05Z) chose this as the seminal feature of 8.5 (the
  end of phase 2). Umbrella: it is done when a viewer plots a course
  through waypoints, sees it bend around the gravity wells and Hill
  spheres of the objects it passes (NAV.6 and its subitems), edits the
  waypoints (NAV.58), saves the course and opens it again on the map
  (NAV.4, NAV.23). It is the pull model of the gravity map (GEN.198,
  9.0) seen from the course planner: NAV.6 reuses the evaluator once it
  exists.
  Headline: 8.5 (courses that bend around gravity wells).
  Prerequisites: NAV.6, NAV.49, NAV.58, NAV.4, NAV.23. Related: GEN.198,
  MAP.167.

- [ ] **NAV.60 What is nearby reads by the cell key: one statement for every kind, including objects with no sector and unbuilt scatter rows**
  Step 3 (design note sections 6.4 and 9 item 3). `objects_within` in
  `db/near.py` (NAV.43, built) and every other reader of a sphere use
  the `VALUES` list of at most 256 key runs joined with `STRAIGHT_JOIN`
  and `FORCE INDEX` to each positioned table, then the exact distance
  test in the same statement, for every kind with a position: systems,
  phenomena, facilities, comets, asteroid fields, and the unbuilt
  scatter rows (the "not generated yet" rows come from the same
  statement on the scatter tables). Extended objects search `r + the
  widest placed radius` and measure to the edge. Objects whose key is
  `NULL` (outside the box) come from a `UNION ALL` branch. The
  sector-table path is removed, not kept beside it. The shape of the
  search over several tables (a shared `WITH rg(lo, hi) AS (VALUES ...)`
  reused by a `UNION ALL` of per-table selects, against one
  `located_object` side table) was not measured: Core measures the union
  first. Done: results equal the old path's on a seeded galaxy, a 50 pc
  sphere at the Sun's density takes tens of milliseconds, and a
  plan-shape test (`EXPLAIN FORMAT=JSON` names the covering index and
  the join order) fails if a MySQL upgrade chooses a full scan. Owner:
  Interface (2026-10-15; moved off Core by the lane queue, Core owns the
  star model or the key it builds on). Boss (2026-10-15 5:32 PM EDT):
  "if it's not ready in time before the spherical sector layout is, then
  we'll skip it and move to the next sequence until it's finished
  research and we can plan it." The research is finished and merged
  (#1659), so stage C is live.
  Prerequisite: GEN.276. Related: NAV.43, NAV.44, NAV.45,
  NAV.61, TEST.165.
  Owner: Core 2 (Boss, 2026-10-16, 3:00 PM EDT: the region query goes to
  Core 2), after the infinite zoom tree and its labels.
  Design: [docs/design/spherical-region-queries.md](design/spherical-region-queries.md)

- [ ] **NAV.61 Nearest-first by a radius ladder, and a dense-core guard with an early-stop count**
  Step 4 (design note sections 6.5 and 9 item 4). "The nearest K"
  (NAV.12's route hops and the ordering of the nearby page) asks for 2
  pc, then 4, 8, 16 and 32, stopping at the first sphere that holds K
  rows, so it is cheap wherever it is; counts shown capped ("300+") come
  from `SELECT COUNT(*) FROM (SELECT 1 ... LIMIT 10001)`, at most 11 ms
  even where the sphere holds a million systems; the API refuses, or
  only counts, a query whose early-stop count passes a limit (Core's
  call; 100,000 systems is about 70 ms of scan) with the text that says
  to shrink the radius. A 50 pc sphere at the core holds about 7.9
  million systems, so no uncapped top-200 runs there. Keyset paging
  `(distance, kind, id)` stays. Done: the ladder returns the nearest K
  at the core and at the Sun in milliseconds, the cap fires with a clear
  message at the core, and NAV.44's page shows "300+" from the
  early-stop count. Owner: Interface (2026-10-15; moved off Core by the
  lane queue, Core owns the star model or the key it builds on). Boss
  (2026-10-15 5:32 PM EDT): "if it's not ready in time before the
  spherical sector layout is, then we'll skip it and move to the next
  sequence until it's finished research and we can plan it." The
  research is finished and merged (#1659), so stage C is live.
  Owner: Core 4, "Core 4 (2026-10-16)" (Boss approved the lane at 4:31
  PM EDT, 2026-10-16; queue in the lane report); held until its blocker
  merges: NAV.60 (Core 2). Core 4 takes it when it unblocks.
  Prerequisite: NAV.60. Related: NAV.12, NAV.44, NAV.45.
  Design: [docs/design/spherical-region-queries.md](design/spherical-region-queries.md)

- [ ] **NAV.62 Gravity of the stars and objects a cell leaves out: work out how a course and the orbital update account for the unresolved mass, derived from the density model**
  Boss (2026-10-16, 5:01 PM EDT, on the cell cap of GEN.364): "in the
  core there are just too many stars to make all of them, no matter how
  scientifically accurate it would be. I would like to figure out how we
  can account for the gravity of items left out, but let's put that as a
  TODO for the courses that bend around gravity wells phase." A cell
  that expects more systems than its cap resolves the brightest and
  leaves the rest unmade; GEN.281 already derives the unresolved count
  and mass of a cell from its address (`galaxy/resolution.py`), from the
  same density model as every other cell. Done: a written study (a
  section in docs/design/realistic-galactic-core.md) that says how the
  pull of that unresolved mass is added: as a smooth field summed per
  cell from the derived unresolved mass and the cell centre, how it
  joins the galactic potential of GEN.285 and the keep-out spheres of
  NAV.54, whether a course that bends around wells (NAV.26) and the
  orbital update need it at every zoom or only inside the dense core,
  its cost per query, and a recommendation with the work it needs as new
  items. No fixed table: the mass comes from the density model and the
  address, never from a stored list. No build in this item. Owner: not
  assigned; the courses lane of 8.5.
  Core 1 note (2026-10-16, 6:13 PM EDT): belts and asteroid fields have
  no mass model, so DB.55 and GEN.337 store a stated zero mass for them;
  the study says whether a real model is needed (GEN.295 refines nebula
  mass).
  Prerequisites: NAV.54. Related: GEN.364, GEN.281, GEN.303, GEN.285,
  NAV.26, NAV.6.
  Design: [docs/design/realistic-galactic-core.md](design/realistic-galactic-core.md)

## GEN: Generation and physics

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
  Prerequisites: GEN.33, GEN.146.
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
  Prerequisites: GEN.350, GEN.351, GEN.352, GEN.353, GEN.354.
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
  has M-dwarf and G/K sub-types (design doc 5.1). Decided (Boss,
  2026-10-10 18:44Z, defaults approved; default R): rocky rogues of 10
  to 16 Earth masses become class R
  rather than a widened S (cap S at 13,500 km); a 1-bar "dune world" is
  a 10 percent variant of N.

  - [ ] **GEN.33 One class per PR, each with its tests**
    Suggested split: R and S (the commonest missing types) first; then
    U, W, X, Y and Z. Each adds the class to `PLANET_CLASSES`, the
    reference pages, the rogue flag, moon eligibility, and a
    distribution test over 1,000 generated systems.
    Plan (2026-10-07): Moved to phase 2 under GEN.90 (the class
    refactor), after the habitability score GEN.89.
    Research (2026-10-09, atmospheres-retention-and-classes.md): confirm
    the PR order: retention module first, then R, the S fix, U, W/X, Y,
    N, Z, the Hycean flag, the GEN.29 sweep.
  Split (2026-10-16, Boss 1:07 PM EDT: look at each item for smaller
  chunks and shared work): this item is now the umbrella and is done
  when its subitems are: GEN.350, GEN.351, GEN.352, GEN.353, GEN.354.
  GEN.33's own rule is 'one class per PR' and the item lists six
  remaining classes (R, U, W, X, Y, Z) with different dependencies: R
  and Y need the retention module (GEN.91's first PR), Z needs GEN.92's
  pre-main-sequence rule and so cannot come before GEN.92, which
  currently waits for all of GEN.28.
    Done (added (2026-10-16 review)): a parametrized class-contract test
    driven by `tuning.PLANET_CLASSES` (zone flags, weight, radius and
    mass ranges, atmosphere, moon blacklist, rogue flag, habitable list,
    reference page) plus a helper that samples N systems per zone, so
    each class PR adds a table row and not a hand-written test; the
    research PR order becomes the order of the children of GEN.28.

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
  Prerequisites: GEN.28, GEN.27, GEN.91, GEN.92, GEN.146.
  Plan (2026-10-07): Under GEN.90: the refactor is the sweep, so it
  stays with it in phase 2.
  Research (2026-10-09, atmospheres-retention-and-classes.md): known
  oddities to fold in: classes A and B hold 13 to 17 kPa at S = 4.4 and
  v_esc 3.7 km/s, about 50 times beyond the cosmic shoreline; N sits in
  the ecosphere at S about 1.05 (below runaway); T and I overlap; rogue
  S reaches 17,600 km (recommend class R at 10 to 16 Me, or cap S at
  13,500 km); K, L and P fail the Armstrong and 20 kPa limits (P needs
  the stage cap from GEN.92).
  Checked 2026-10-10 22:19Z for Boss's generation-first priority: this
  waits on GEN.27 and GEN.28 (the new planet classes) and sits under the
  GEN.90 refactor, so it cannot move up yet; it starts when they land.

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

  - [ ] **GEN.42 A block-first draw for the backfill (umbrella)**
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
    Split (2026-10-16, Boss 1:07 PM EDT: look at each item for smaller
    chunks and shared work): this item is now the umbrella and is done
    when its subitems are: GEN.308, GEN.309, TEST.182. After the
    research re-scope it holds a sampling algorithm with a statistical
    test harness, the per-block choice between two draw paths wired into
    the backfill, and an addition to another item's golden test.
  Prerequisites: GEN.308, GEN.309, TEST.182.

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
    Split (2026-10-16, Boss 1:07 PM EDT: look at each item for smaller
    chunks and shared work): this item is now the umbrella and is done
    when its subitems are: GEN.306, GEN.307, TEST.181. Two independent
    code changes in different modules (population seeds in
    population/model.py; the object-ID collision counter in db/store.py)
    and the test that proves both, all under one blocker of TEST.77,
    DB.17, ADM.61 and API.17.
    For a sector made under a directive the rule is seed, version,
    address, the directive digest and the attempt number
    (`generation/directives.py` seeds attempt a >= 1 from the
    directives' digest).
    Prerequisites: GEN.306, GEN.307, TEST.181.

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
    Prerequisites: OPS.57, OPS.58.

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
    first exits at once. Prerequisites: OPS.30.
    Research (2026-10-09, db-check-and-parity-repair.md): the parity
    export and checksum must not cover the columns the orbit update
    rewrites, or the daily run invalidates all of them; the run also
    refiles systems and phenomena into other sectors, so sector
    membership of placed objects changes daily.
    Research (2026-10-09, ops-scheduling-and-rotation.md): restructure:
    implement the logic once in `planetgen.cli.nightly`
    (`planetgen.cli.maintenance` already exists: maintenance mode on and
    off, UX.112) (steps,
    lock, `--if-due`, `--dry-run`, `--only STEP`, rotating log file,
    exit codes 0/1/2/75) and make `scripts/maintenance.sh` and `.ps1`
    thin launchers that find the right Python and run it as the web
    user. Take `locks/maintenance.lock` plus a check of
    `jobs.active_job()`; exit 75 when skipped; alert after N skipped
    days. Treat Redis as optional. Iterate every galaxy database. Run
    orbits in-process or as `planetgen.cli.orbits` (this item's text
    says `updateOrbits.py`, which no longer exists). Decided (Boss,
    2026-10-10 18:44Z, defaults approved; default: skip the day while a
    Generate job runs, exit 75,
    alert after 3 skipped days).
    Note (2026-10-09): The Windows parts are dropped by OPS.39 (Boss,
    2026-10-10 04:25Z).
    Needs OPS.30 (the lock helper) first and the thread order is OPS.30,
    OPS.16, OPS.17 (2026-10-16 review).
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
    Prerequisites: GEN.304, OPS.57, GEN.57, GEN.135.

  - [ ] **API.17 Remote generation reproduces what the server would make**
    Done: a remote run (API.12's download, API.13's generation without a
    database) with the same seed and release produces exactly what the
    server would for those sectors, checked by fingerprint (GEN.58); and
    API.8 can verify an upload by re-running a sample of its sectors on
    the server and comparing. Prerequisites: API.12, API.13, GEN.57,
    API.8.
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
    Needs API.8, because its Done changes the verification flow (API.8
    verifies an upload by re-running a sample of its sectors on the
    server and comparing) (2026-10-16 review).
    Design: [docs/design/reproducible-galaxies.md](design/reproducible-galaxies.md)

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
    for Q rows. Decided (Boss, 2026-10-10 18:44Z, defaults approved;
    default Option A): re-scope N,
    extend X, Hycean as an R flag and retire T, or two-character class
    codes.
    Split (2026-10-16, Boss 1:07 PM EDT: look at each item for smaller
    chunks and shared work): this item is now the umbrella and is done
    when its subitems are: GEN.355, GEN.356. Its first deliverable (the
    retention module) is a pure library that GEN.28's R and Y and
    GEN.145 need, but it is bundled with class work and sits after
    GEN.28 in the thread; the class work it describes (Hycean flag,
    extended X, retire T) is already inside GEN.28's R and X.
    Prerequisites: GEN.355, GEN.356.
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
    Decided (Boss, 2026-10-10 18:44Z, defaults approved; default yes):
    dose above 0.1 Sv/yr forces
    equipment tier 3 and above 1 Sv/yr tier 4.
    Prerequisites: GEN.350, GEN.351, GEN.352, GEN.353.
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
    radius get an eccentricity draw (up to about 0.3). Decided (Boss,
    2026-10-10 18:44Z, defaults approved; default yes): the
    swallowed-giant flag.
    Decided (2026-10-09): Boss (2026-10-09 18:37Z): the defaults are
    accepted, including the swallowed-giant flag. The orbit-expansion
    law for adiabatic mass loss is verified: a (M_star + M_planet) =
    constant with the eccentricity unchanged (Jeans 1924; Hadjidemetriou
    1963; Veras et al. 2011), valid when the mass-loss time is much
    longer than the orbital period and the loss is isotropic; impulsive
    loss uses vis-viva and unbinds an orbit when more than half the mass
    goes; engulfment overrides. Worked example: 1 AU round a 1 Msun star
    that drops to 0.6 Msun moves to 1/0.6 = 1.67 AU.
    Split (2026-10-16, Boss 1:07 PM EDT: look at each item for smaller
    chunks and shared work): this item is now the umbrella and is done
    when its subitems are: GEN.357, GEN.358 (DB.59 was found already built). Three deliverables,
    and the first two are needed a version earlier than the third:
    GEN.99 (8.1-4) must rebuild planets with the same `apply_cloud` (its
    own research note says so) and draws an age-capped population from
    the nebula's stored `age_myr`, which the design says 'waits on a
    schema column' (nebula-and-asteroid-field-classes.md, Applying the
    rules at generation, and Nebula volume backfill step 3).
    Prerequisites: GEN.357, GEN.358.
    Design: [docs/design/nebula-and-asteroid-field-classes.md](design/nebula-and-asteroid-field-classes.md)

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
  Rewrite to the remainder (2026-10-16 review): the populations and
  per-kind factors are built (GEN.234, GEN.236); what is left is the
  population-dependent rules for the type draws. Fold into GEN.269 or
  close once GEN.270 to GEN.275 land.
  Research supplied (Boss, 2026-10-16, 5:55 PM EDT, in the project
  chat): the two constants Core 4 could not look up, the close-binary
  fraction against [Fe/H] and the encounter-disruption limit for wide
  pairs with the Jacobi radius as the outer limit. Sources named by Boss:
  Badenes 2018, Moe, Kratter and Badenes 2019, El-Badry and Rix 2018,
  El-Badry, Rix and Heintz 2021, Kratter and Lodato 2016, Bate 2018,
  Heggie 1975, Bahcall, Hut and Tremaine 1985, Jiang and Tremaine 2010,
  Stegmann 2024. The reference list is unverified and Core 4 checks it.
  Filed as derivation rules of the type draw, not fixed values, in
  docs/design/multistar-and-compact-systems.md section 6.6: the close
  share follows the cell's [Fe/H] and the wide-pair limit is computed
  from the cell's density, dispersion, age and tidal field.
  Boss answered the GEN.103 cards (2026-10-16, 6:12 to 6:15 PM EDT):
  (1) the close-binary trend is the pooled one, read at the star's own
  [Fe/H]; (2) the metallicity effect applies to every star mass, with a
  built-in check that the numbers act as expected (filed as TEST.196); (3)
  wide pairs have no metallicity dependence; (4) the encounter limit is
  a sharp cut at a_1/2; (5) the fixed 10,000 AU wide-pair cap is
  replaced by the derived limit. The last card, whether the metallicity
  effect fades above about 2 solar masses, was answered "all masses":
  no fade, and the check is the safeguard.
  Design note landed (Core 4, 2026-10-16, 7:16 PM EDT, #1891):
  docs/design/placement-by-population.md. The build follows in its own PR
  and covers 3.1 (close pairs by [Fe/H], every mass) and 3.3 (the derived
  wide-pair limit). Section 3.2, massive multiples sitting in young
  regions through `prevalence`, is NOT in that PR: it is GEN.373, with the
  multi-star items (GEN.269, GEN.270 to GEN.275), and its cap on the rise
  is the open constant (a design default, not a value from a source).
  GEN.103 closes when the 3.1 and 3.3 build merges. The metallicity
  effect applies to every mass, so the numbers need a ceiling that is not
  a per-mass value: Core 4's recommendation (2026-10-16, 7:25 PM EDT): no new
  constant; apply the pooled [Fe/H] factor F/F(0) to the odds p/(1-p) of
  a close companion, not to the chance p, for every mass. It cannot pass
  1, and a star already likely to have a companion rises less. Core 4's
  check values: G dwarf 0.43, 0.36, 0.25, 0.12 at [Fe/H] -3.0, -1.0,
  -0.2, +0.5 against Moe, Kratter and Badenes 53+-12, 40+-6, 24+-4,
  10+-3% (all within 1 sigma); an O star at [Fe/H] -0.7 goes 0.45 to 0.58
  (+29%) against the "no change within about 30%" read for O and early-B
  stars in the Large Cloud and stays under Sana's 0.69+-0.09 (plain
  scaling gives 0.77). Boss chose the odds form on a card (2026-10-16, 7:19 PM EDT): Decided,
  the metallicity factor applies to the odds p/(1-p) of a close
  companion, for every mass. Core 4 builds it.
  Owner: Core 4, "Core 4 (2026-10-16)" (Boss approved the lane at 4:31
  PM EDT, 2026-10-16; queue in the lane report, queue position 5).

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
  Headline: 8.3 (the orbital update system).
  Prerequisites: GEN.115, GEN.109, GEN.110, GEN.289, GEN.290, GEN.291.
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
    (integrated in the potential once it has a vertical offset). Decided
    (Boss, 2026-10-10 18:44Z, defaults approved; defaults taken):
    constant 229.3 km/s whatever
    the galaxy's size; Sun distance 8.128 kpc as the test radius and 8.2
    kpc for placing the Sun in the density model (change
    `GALACTIC_CENTER_DISTANCE_LY` from 25,800 ly); NFW only for now.
    Split (2026-10-16, Boss 1:07 PM EDT: look at each item for smaller
    chunks and shared work): this item is now the umbrella and is done
    when its subitems are: GEN.338, GEN.339, GEN.340. Four deliverables
    with different readers and timing.
    Prerequisites: GEN.338, GEN.339, GEN.340.
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
    sees that star. Decided (Boss, 2026-10-10 18:44Z, defaults approved;
    defaults taken): full Hill
    spheres with cover lists rather than clipping to the sector stencil;
    a star cluster counts as its members; no warnings for unbound flybys
    or hierarchical members.
    Boss (2026-10-15, 11:06 PM EDT): the influence set also includes the
    planets and moons of the system by their point masses, so planets and
    moons can be captured or ejected (GEN.289 to GEN.291).
    Note from Core 1 (2026-10-16, 12:24 AM EDT, GEN.265 merged, #1766):
    the influence rule reads scatter_systems.point_masses_solar
    (db/scatter_systems.py), and member offsets come from
    physics/system_tree.py (stored node phases); the nightly update
    itself writes no member or node row.
    Split (2026-10-16, Boss 1:07 PM EDT: look at each item for smaller
    chunks and shared work): this item is now the umbrella and is done
    when its subitems are: GEN.341, GEN.342, GEN.343, GEN.344. It holds
    four deliverables that other items need at different points: a
    search (GEN.110, GEN.198, GEN.290 read the set), the acceleration
    and update, the perturbed pass that nobody else builds
    (capture-and-ejection-design.md section 3, work items 2 and 4:
    GEN.290 only decides when to switch a system into it and GEN.291
    only reads its result), and the warnings (GEN.111).
    Scope after the review (2026-10-16 review): the influence set, the
    kernel that feeds the stepper, the perturbed pass and the warnings;
    no integrator of its own for stars (the stepper is child c of
    GEN.115).
    Precision note (GEN.372 audit, Boss card 2026-10-16, 7:36 to 7:38 PM EDT):
    the long force sums here use compensated summation (Neumaier or fsum,
    or mpmath for the long integrations from the galactic centre); the audit
    decided they wait for this item, not before.
    Prerequisites: GEN.341, GEN.342, GEN.343, GEN.344.
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
    `IMF_BREAKS_SOL[0]` is 0.08) or clamp merger stars to 0.08. Decided
    (Boss, 2026-10-10 18:44Z, defaults approved; defaults taken): the
    physical table, with an
    `ALWAYS_DESTROY_TERRESTRIAL` setting for the literal rule (about a
    quarter of Earth-Earth hits at galactic speeds destroy both);
    hit-and-run keeps both bodies, eroded; the fusion check covers brown
    dwarf and star involvement.
    Split (2026-10-16, Boss 1:07 PM EDT: look at each item for smaller
    chunks and shared work): this item is now the umbrella and is done
    when its subitems are: GEN.345, GEN.346. A 600-line design (outcome
    table, giant merger, fusion check, merger-remnant star, event-made
    fields) behind one item.
    Detection comes from the pair closest-approach library; this item
    keeps the outcome rules (child a) and applying them (child b)
    (2026-10-16 review).
  Prerequisites keep the umbrellas; each umbrella closes with its
  children (2026-10-16 review).
    Prerequisites: GEN.345, GEN.346.
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
  Needs the Hill-warning child of GEN.109 (child d), not the whole
  umbrella (2026-10-16 review).
  Prerequisites: USR.3, GEN.344.
  Design: [docs/design/orbital-updates.md](design/orbital-updates.md)

- [ ] **GEN.129 Multi-star systems of up to twelve stars**
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
  yes/no; the table's mean is 1.38 stars per system. Decided (Boss,
  2026-10-10 18:44Z, defaults approved; defaults taken): faithful rarity
  (about 1 in 33,000 systems have
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
  Boss asked (2026-10-15) whether a black hole or neutron star can have
  a star system. Today a scattered or rolled one is a standalone row
  with no system (`run_sector.py:449`); a system is built around one
  only on request with `phenomenon --anchor-system`, and a remnant
  passed to `StarSystem` skips binary generation (`system.py:255`). The
  design note (section 5.12) points at this item for the rules: the same
  fractions let a Wolf-Rayet star's companions (GEN.243) and the
  remnants' companions (GEN.130) be drawn as trees.
  Decided (Boss, 2026-10-15 5:39 PM EDT): scatters and backscatters get
  no systems, planets or moons, so bright stars, massive stars, neutron
  stars, black holes and the other scattered objects are bare rows;
  planets and moons are made only when a sector is being filled with
  everything; a sector fill must represent every kind of object, and
  when a fill meets a phenomenon that can have a system (or makes one
  that can) it generates that system if the real-world probabilities
  give it one; nebulae need stars and phenomena, not planets. For this
  item: multi-star systems around a scattered object appear at sector
  fill only.
  Reworded (Boss, 2026-10-15 5:51 PM EDT;
  docs/design/multi-star-systems.md): companions of stars and star-like
  objects are generated at scatter time, planets, moons and belts at
  sector fill; a scatter row is a system with a barycentre, a total mass
  and stored members. This item's steps (a) to (d) are built as the
  member-tree generator GEN.261 (with the scatter columns GEN.260); what
  stays here is the scene and map work (e), the orbital update of built
  systems (f) and the planet zones and flux sums (g). Replaces the
  earlier "multi-star systems around a scattered object appear at sector
  fill only".
  Catalogue (Boss 2026-10-15 6:08 PM EDT;
  docs/design/multi-star-systems.md sections 13 to 16): 'up to seven'
  becomes 'up to twelve' (the cap is a setting, GEN.267, default 12).
  Add a trapezia sentence: young non-hierarchical groups under their
  decay age are the one exception to hierarchical trees, stored as a
  `group` node (GEN.273).
  Split (2026-10-16, Boss 1:07 PM EDT: look at each item for smaller
  chunks and shared work): this item is now the umbrella and is done
  when its subitems are: DB.53, GEN.311, GEN.312, MAP.213, GEN.313,
  GEN.314. GEN.129's text says its steps (a) to (d) are built as GEN.261
  and GEN.260, but those built only the member-tree generator and the
  scatter_members and scatter_nodes tables.
  Owner: Core 3 (Boss, 2026-10-16, 2:42 PM EDT: Core 3 takes the 8.1 new
  generate scatters line after GEN.280; the multi-star train follows its
  spine, 3:00 PM EDT). This replaces any earlier Core or Core 1 owner
  line on this item.
  Prerequisites: DB.53, GEN.311, GEN.312, MAP.213, GEN.313, GEN.314.

- [ ] **GEN.130 Compact-remnant systems at sector fill: companions and pulsar planets for neutron stars and black holes**
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
  question for Boss, default as stated). Decided (Boss, 2026-10-10
  18:44Z, defaults approved; default: only through GEN.100's scatter,
  which promotes a scattered
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
  multistar-and-compact-systems.md section 8. Decided (Boss, 2026-10-10
  18:44Z, defaults approved; defaults taken): a globular-cluster model
  does not exist, so the
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
  Boss asked (2026-10-15) whether a black hole or neutron star can have
  a star system. Today a scattered or rolled one is a standalone row
  with no system (`run_sector.py:449`); a system is built around one
  only on request with `phenomenon --anchor-system`, and a remnant
  passed to `StarSystem` skips binary generation (`system.py:255`). The
  design note (section 5.12) points at this item for the rules: this
  item decides the company a black hole or neutron star keeps (white
  dwarf, young neutron star none, millisecond pulsar 1% with planets,
  dormant black hole with a Gaia-type companion), so a scattered remnant
  can be promoted to a system. Nothing in the whole-sphere train waits
  on it.
  Decided (Boss, 2026-10-15 5:39 PM EDT): scatters and backscatters get
  no systems, planets or moons, so bright stars, massive stars, neutron
  stars, black holes and the other scattered objects are bare rows;
  planets and moons are made only when a sector is being filled with
  everything; a sector fill must represent every kind of object, and
  when a fill meets a phenomenon that can have a system (or makes one
  that can) it generates that system if the real-world probabilities
  give it one; nebulae need stars and phenomena, not planets. For this
  item: the company a black hole or neutron star keeps (companions,
  disk, planets) is decided when its sector fills, at the real-world
  fractions, and a fill that already holds such a phenomenon builds its
  system if it draws one; a scattered row never carries one, so
  promoting a remnant at scatter time (the earlier default) no longer
  applies. Boss has ordered the phenomena catalogue and fill rules
  (docs/design/phenomena-catalogue-and-fill.md, in research); items
  follow from it.
  Rescoped (catalogue section 8 item 3, Boss's rules of 2026-10-15):
  this is now the compact-remnant content of the system pass of the
  fill-stages item, not a scatter-time promotion. Slices (a) to (c)
  above become what the pass builds for neutron stars, millisecond
  pulsars and stellar black holes at the probabilities of catalogue
  section 5: the companion, the kick-and-survival check, the
  Holman-Wiegert zones, pulsar planets at 0.7%. Done: a filled sector
  holding a black hole or neutron star gives it a system at those odds,
  and a scattered one never has one before the fill. The white-dwarf and
  planetary-nebula-star rows belong to the fill-stages item. Owner: Core
  (2026-10-15).
  Decided (Boss, 2026-10-15 5:51 PM EDT), refining the 5:39 PM rule:
  scatters still carry no planets or moons, but the stellar companions
  of stars and star-like objects (stars, former stars, protostars,
  neutron stars, black holes), including neutron star pairs, black hole
  and neutron star pairs, black hole pairs and stars with neutron star
  or black hole companions, are generated at scatter time, because such
  systems matter for gravity; planets and moons still come only at
  sector fill. Multi-star research is under way (docs to follow; the
  coordinator relays its items). For this item: the companions of
  neutron stars and black holes (white dwarf, main-sequence or giant
  companions, neutron star and black hole partners) are decided when the
  object is scattered, at the observed fractions; pulsar planets stay at
  sector fill.
  Reworded (Boss, 2026-10-15 5:51 PM EDT;
  docs/design/multi-star-systems.md): companions of stars and star-like
  objects are generated at scatter time, planets, moons and belts at
  sector fill; a scatter row is a system with a barycentre, a total mass
  and stored members. The companion orbit laws of neutron stars and
  black holes (including their pairs) move into GEN.266 and the
  member-tree generator GEN.261; neutron star and black hole companions
  are drawn at scatter time (GEN.262, GEN.263), and the black hole
  companion fraction is 1.25%, not 10% to 30%. What stays here: the
  white dwarf and planetary-nebula-star rows of the system pass, pulsar
  planets at fill (0.7% for millisecond pulsars) and the radiation
  environment.
  Catalogue (Boss 2026-10-15 6:08 PM EDT;
  docs/design/multi-star-systems.md sections 13 to 16): slice (c), the
  accretion luminosity and the accreting compact pair, is built by
  GEN.272 (X-ray binaries, microquasars, spider pulsars), which reads
  GEN.269's interaction states.
  Built-on note (Core 1, 2026-10-16, 8:41 AM EDT, GEN.251 merged as
  #1815): the planet model gives a recycled pulsar up to about 10
  planets, and a planetary nebula's circumbinary disk is not modelled;
  both belong to this item. The system pass (generation/system_pass.py,
  tuning.SYSTEM_PASS_ODDS) already rolls a millisecond pulsar's planets
  at 0.7% on its own seed.
  Split (2026-10-16, Boss 1:07 PM EDT: look at each item for smaller
  chunks and shared work): this item is now the umbrella and is done
  when its subitems are: GEN.317, GEN.318. After the rescopes, GEN.130
  holds two unrelated rows of the system pass: the neutron-star side
  (pulsar planets and the radiation environment) and the white-dwarf and
  planetary-nebula-star side.
  Owner: Core 3 (Boss, 2026-10-16, 2:42 PM EDT: Core 3 takes the 8.1 new
  generate scatters line after GEN.280; the multi-star train follows its
  spine, 3:00 PM EDT). This replaces any earlier Core or Core 1 owner
  line on this item.
  Owner: Core 3, in its multi-star and cluster line (the 3:15 PM EDT
  lane queue); a line here that named Core 4 was the TODO thread's
  counting error, corrected 2026-10-16, 4:50 PM EDT.
  Prerequisites: GEN.317, GEN.318.

- [ ] **GEN.134 Tune the star populations to the observed star-formation profile by galactic radius**
  From GEN.133 (Boss's 09:20Z request; the analysis is `docs/design/star-types-by-galactic-radius.md`, PR #807), whose four proposals were left unbuilt. Done: (1) the young and intermediate populations are weighted by the observed star-formation profile (peak at 5 kpc, about -0.28 dex per kpc beyond, a dip inside 3 kpc, the Central Molecular Zone as its own small young region); (2) the young population's B share is lowered, or its weight cut, so local B stars come to about 0.04%; (3) the bulge gets a small young tail (about 10% under 5 Gyr, between the HST and microlensing figures); (4) a metallicity gradient is added only if planet occurrence is later tied to it, otherwise the note records why not. The changes are reproducible (GEN.56), and a test compares the star type shares at the core, mid radius and rim against the note's table. Boss to confirm which of the four proposals he wants before the build starts.
  Research (2026-10-09, globular-clusters.md): item 4 now has its
  reason: planet occurrence is tied to metallicity in clusters (see
  GEN.158).
  Scope after the review (2026-10-16 review): the star-formation profile
  by radius (peak at 5 kpc, -0.28 dex per kpc beyond, the dip inside 3
  kpc), the B share and the bulge's young tail, without the Central
  Molecular Zone, which GEN.282 owns.
  Owner: Core 4, "Core 4 (2026-10-16)" (Boss approved the lane at 4:31
  PM EDT, 2026-10-16; queue in the lane report, queue position 11; start
  after Core 3 merges GEN.283 and GEN.333 (tuning.py and the profile)).
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
  ships with a generator epoch bump (the Python 3.9 CI leg that blocked
  TEST.77 is gone, OPS.42 merged in #1291). Related: GEN.66 follow-up on pinning astropy's constant
  sets.
  The epoch counter it bumps is built by the first child of OPS.28;
  GEN.135 needs it so the order is not a circle.
  Prerequisites: OPS.56.
  Design: [docs/design/reproducible-galaxies.md](design/reproducible-galaxies.md)

- [ ] **GEN.136 Fingerprint tree: a stored leaf digest per sector and a shell-band tree of digests**
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
  the `tests/galaxy_fingerprint.py` stand-in. Decided (Boss, 2026-10-10
  18:44Z, defaults approved; default: round to 9 digits): keep exact
  float hashing as built, or
  round?
  Split (2026-10-16, Boss 1:07 PM EDT: look at each item for smaller
  chunks and shared work): this item is now the umbrella and is done
  when its subitems are: GEN.304, DB.49, GEN.305. It bundles the pure
  encoder that three other items need early (the battery digest in
  OPS.28, the golden test TEST.77, the order test in GEN.57) with two
  things nobody else waits for: a schema change and the save-time tree.
  Prerequisites: GEN.304, DB.49, GEN.305.
  Design: [docs/design/generation-determinism.md](design/generation-determinism.md)

- [ ] **GEN.139 Orbit-update thresholds: per-object epoch, path-length rule and what the 0.01 mpc applies to (GEN.106 built)**
  Follow-up to GEN.106: (a) struck 2026-10-16 (a per-object `epoch_unix`
  and `next_update_due` already exist on every moving table, schema v66,
  GEN.106); (b) define the movement rule on path length, not
  displacement (Phobos and Deimos can never move 100,000 km); (c) state
  whether 0.01 mpc applies to the galaxy-frame velocity (every star due
  every 16 days at 230 km/s, so "1 day = 1 day" runs touch about 1/15 of
  the stars, 6.2e7 writes a day per 1e9 stars) or to peculiar motion
  only; `next_update_due` can use the current speed; (d) use
  orbital-solvers-and-integrators.md section 5 as the justification; (e)
  "visible" could mean half a pixel at the deepest zoom where the object
  type is drawn (needs a map-scale value). Decided (Boss, 2026-10-10
  18:44Z, defaults approved; defaults taken): stars follow their
  analytic orbit and are rewritten
  at sector exit or when the influence set changes, planets and moons
  are computed from phase on read with thresholds only scheduling
  perturbation re-fits; the deepest star zoom is the 4 pc sector view.
  Moved first in the Orbital update core thread and a prerequisite of
  the stepper (child c of GEN.115) (2026-10-16 review). Restated: a
  planet or moon is placed from its stored vector and epoch by the
  closed-form step on read, and the thresholds only schedule
  perturbation re-fits; GEN.289 replaces the phase columns that the
  earlier wording relied on.
  Prerequisites: none.
  Design: [docs/design/orbital-solvers-and-integrators.md](design/orbital-solvers-and-integrators.md)

- [ ] **GEN.140 Orbital math guards the edge-case table adds (GEN.108 built)**
  Follow-up to GEN.108: adopt the edge-case table of
  orbital-solvers-and-integrators.md section 7 as the checklist. Add the
  guards not in the built item: Stumpff functions by series for abs(z) <
  1, the `fmod(t, P) / P` phase form, DOUBLE-only period and phase
  columns, radial orbits (h = 0), retrograde equinoctial singularities
  (use state vectors), the r = 0 and b = 0 guards and a sub-step cap for
  the galactic centre (the finite-acceleration test at (0, 0, 0) belongs
  to the first child of GEN.115). The corrected universal
  equation is in orbital-updates.md
  10.6.
  Moved to the Positions thread of phase 1, before GEN.289, and a
  prerequisite of the perturbed pass (child c of GEN.109), which runs
  into its guards (2026-10-16 review).
  Prerequisites: none.
  Design: [docs/design/orbital-solvers-and-integrators.md](design/orbital-solvers-and-integrators.md)

- [ ] **GEN.141 Faster Kepler solver (Mikkola or Markley) with brentq as fallback**
  Swap `solve_eccentric_anomalies` for Mikkola or Markley with `brentq`
  as fallback (20 to 30 times faster for comets). Low priority. When
  planets gain eccentricity (MAP.70 and MAP.62 are built), reuse the
  same non-iterative solver in `body_positions.py` and
  `orbitpositions.js` so the parity test stays simple.
  Math Research (2026-10-14): production code does not yet call the vector
  `solve_eccentric_anomalies` (1.36 us against 20.9 us for the scalar
  `brentq`); use it wherever updates are batched (GEN.212).
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
  sector-relative values. Decided (Boss, 2026-10-10 18:44Z, defaults
  approved; default yes).
  Narrowed (2026-10-16 review): draws the peculiar velocities at
  generation (thin and thick disc and bulge dispersions); the velocity
  and mass columns and the rogue sector-relative offsets come from the
  second and third children of GEN.294. Needs GEN.294:b, placed after
  them.
  Prerequisites: GEN.336.
  Design: [docs/design/collisions-and-mergers.md](design/collisions-and-mergers.md)

- [ ] **GEN.143 A collision_events table, an admin report and a test that runs the whole collision path**
  `collision_events` table, the admin report in design section 7 and a
  "last collisions" admin view, plus a test that seeds two rogues on a
  collision course through the whole update to the field, giant merger
  and star branches. Cross-link collisions-and-mergers.md from
  nebula-and-asteroid-field-classes.md and interstellar-object-rates.md.
  Split (2026-10-16, Boss 1:07 PM EDT: look at each item for smaller
  chunks and shared work): this item is now the umbrella and is done
  when its subitems are: DB.56, GEN.347. GEN.110's own Done says 'every
  event in the admin report', but the table and report are GEN.143 in
  the next phase, so GEN.110 would ship with no place to write.
  Scope after the review (2026-10-16 review): the admin report
  (collisions-and-mergers.md section 7), the last-collisions admin view
  and the end-to-end test through the field, giant-merger and star
  branches; the table is DB.56.
  Prerequisites: DB.56, GEN.347.
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
  Moved to the first place of the new Galactic gravity thread, ahead of
  GEN.115 child a, which edits the same constants and tests at the same
  radius (2026-10-16 review).
  Prerequisites: none.
  Design: [docs/design/galactic-potential.md](design/galactic-potential.md)

- [ ] **GEN.145 Class S atmosphere rule: S keeps air unless the shoreline ratio is over 30**
  Class S landed with `atmosphere = None` in the ecosphere and cold
  zones, which breaks the shoreline by a factor of 50 or more. S is
  airless only for shoreline ratio over 30 (T_eq above 650 to 1,200 K).
  It changes existing output, so it needs GEN.92's "unchanged until
  regenerated" guard.
  Needs the shoreline functions (child a of GEN.91) and the
  unchanged-until-regenerated guard of GEN.92 (2026-10-16 review).
  Prerequisites: GEN.355, GEN.92.
  Design: [docs/design/atmospheres-retention-and-classes.md](design/atmospheres-retention-and-classes.md)

- [ ] **GEN.146 Teff-dependent habitable zone from the Kopparapu table**
  Replace the fixed 1.1 and 0.53 flux limits in
  `orbits.calculate_habitable_zone` with the Kopparapu table (design doc
  4.1): the fixed limits put the outer edge 33 percent too near for a
  3,200 K star and at 1.37 AU against 1.68 AU for the Sun. It changes
  the zone of every star, so it needs reproducibility handling (the
  GEN.57 family) and the epoch bump.
  Does not need GEN.91 (nothing in the Kopparapu table touches
  atmosphere retention); placed right after GEN.33 and a prerequisite of
  GEN.27, the R, Y and Z classes of GEN.28 and the GEN.29 sweep, so the
  classes are tuned against the final habitable-zone edges (2026-10-16
  review).
  Prerequisites: none.
  Design: [docs/design/atmospheres-retention-and-classes.md](design/atmospheres-retention-and-classes.md)

- [ ] **GEN.157 The `neighbor_galaxies` table and a verified data file of about 25 real galaxies**
  One Alembic revision, no change to existing tables, and the data file
  checked against McConnachie 2012 (VizieR J/AJ/144/4) and NED. Hold
  absolute magnitude and distance, never a hand-keyed apparent magnitude
  (the draft's recalled values were off by 0.7 to 1.9 mag). The layout
  already has Mpc and Gpc units, so a galaxy's centre fits the 76-bit
  position ID; add type code 13 `galaxy` in `names/object_id.py` only if
  the neighbours use it as their `uid`.
  Scope after the review (2026-10-16 review): the `neighbor_galaxies`
  table and data file for M31 and its group, the other Local Group
  members and the galaxies beyond (M81, M82, NGC 253, Centaurus A, M83,
  IC 342, Maffei 1, Virgo as one object); the nine Milky Way satellites
  are in `galaxy_components` (GEN.237).
  Prerequisites: none.
  Design: [docs/design/multiple-galaxies.md](design/multiple-galaxies.md)

- [ ] **GEN.161 Bright-first fill for cluster sectors**
  Research (2026-10-09, globular-clusters.md, PR #824; handoff in
  /mnt/project-files/research/handoff/globular-clusters.md): giants,
  blue stragglers and horizontal-branch stars first, by luminosity level
  (GEN.44); the dim main sequence on demand. The cost reason: every
  system of a 157-cluster Milky Way is 1.8e8 systems, 10 to 34 days of
  single-worker time at 5 to 16.6 ms each.
  Whole-sphere train (docs/design/intergalactic-and-halo-contents.md,
  Boss 2026-10-15): step 5. No cluster star reaches the bright floor or
  the mass limit, so a cluster is entirely fill-time; this item is its
  fill order.
  Built-on note (Core 2, 2026-10-16, 10:10 AM EDT, GEN.160 merged as #1817):
  `planetgen plan` draws and stores the 157 globular clusters
  (galaxy_components); GalaxyShape has a `clusters` field (array density
  forms leave clusters out); this item places their stars.
  Owner: Core 3 (Boss, 2026-10-16, 2:42 PM EDT: Core 3 takes the 8.1 new
  generate scatters line after GEN.280; the multi-star train follows its
  spine, 3:00 PM EDT). This replaces any earlier Core or Core 1 owner
  line on this item.
  Owner: Core 3, in its multi-star and cluster line (the 3:15 PM EDT
  lane queue); a line here that named Core 4 was the TODO thread's
  counting error, corrected 2026-10-16, 4:50 PM EDT.
  Design: [docs/design/globular-clusters.md](design/globular-clusters.md)
  Prerequisites: GEN.324.

- [ ] **GEN.162 Planet cull and blue stragglers in clusters**
  Research (2026-10-09, globular-clusters.md, PR #824; handoff in
  /mnt/project-files/research/handoff/globular-clusters.md):
  giant-planet probability multiplied by 10^(2 [Fe/H]), capped at 1;
  hard cut: drop planets with a > a_h = G M_host / sigma(r)^2 (2 AU at
  15 km/s, 4.4 at 10, 18 at 5, 440 at 1), with sigma(r) from the King
  model. Blue straggler count N_BSS ~ M_core^0.4 (main-sequence stars of
  1 to 1.7 Msun in the core). Decided by Boss (2026-10-09 19:02Z,
  default taken): a hard cut at a_h rather than a smooth exponential.
  Decided (Boss, 2026-10-10 18:44Z, defaults approved; default taken):
  blue stragglers
  scale as M_core^0.4 while millisecond pulsars and X-ray binaries scale
  with the encounter rate (his text says linearly with Gamma for both).
  Whole-sphere train (docs/design/intergalactic-and-halo-contents.md,
  Boss 2026-10-15): step 7, built after the profile (GEN.235) and the
  metallicity value (GEN.158).
  Split (2026-10-16, Boss 1:07 PM EDT: look at each item for smaller
  chunks and shared work): this item is now the umbrella and is done
  when its subitems are: GEN.325, GEN.326. 'Planet cull and blue
  stragglers' are two independent deliverables in different code (the
  planet draw against the star population), and only the cull is waited
  on: GEN.284 needs the cull, while the blue stragglers need cluster
  stars to exist (the cluster population is 0 until GEN.237, Core 1's
  note).
  Owner: Core 3 (Boss, 2026-10-16, 2:42 PM EDT: Core 3 takes the 8.1 new
  generate scatters line after GEN.280; the multi-star train follows its
  spine, 3:00 PM EDT). This replaces any earlier Core or Core 1 owner
  line on this item.
  Owner: Core 3, in its multi-star and cluster line (the 3:15 PM EDT
  lane queue); a line here that named Core 4 was the TODO thread's
  counting error, corrected 2026-10-16, 4:50 PM EDT.
  Design: [docs/design/globular-clusters.md](design/globular-clusters.md)
  Prerequisites: GEN.325, GEN.326.

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
  Whole-sphere train (docs/design/intergalactic-and-halo-contents.md,
  Boss 2026-10-15): step 7.
  Owner: Core 3 (Boss, 2026-10-16, 2:42 PM EDT: Core 3 takes the 8.1 new
  generate scatters line after GEN.280; the multi-star train follows its
  spine, 3:00 PM EDT). This replaces any earlier Core or Core 1 owner
  line on this item.
  Owner: Core 3, in its multi-star and cluster line (the 3:15 PM EDT
  lane queue); a line here that named Core 4 was the TODO thread's
  counting error, corrected 2026-10-16, 4:50 PM EDT.
  Prerequisites: GEN.317.
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
  spherical r^-3.5 profile with a core of a few kpc. Decided (Boss,
  2026-10-10 18:44Z, defaults approved; default: 0.5 dex and the S_N
  ranges above): the mass-function
  width and S_N ranges.
  Whole-sphere train (docs/design/intergalactic-and-halo-contents.md,
  Boss 2026-10-15): step 7; also the source of seeded fictional
  satellites for other shapes.
  Design: [docs/design/globular-clusters.md](design/globular-clusters.md)

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

- [ ] **GEN.197 Orbital event record: the `record_event` writer and the object-ID rules for ejection, capture, merger and split**
  Source: docs/design/object-id-options.md section 0 (Boss decided
  2026-10-09 22:39Z: 80-bit birth-location ID). Left over from GEN.172
  and TEST.110 when Foundations lane 1 merged them (PR #1075). Done: the
  orbital-update events that create these cases apply the ID rules: an
  ejected planet keeps its ID and becomes a rogue-planet row; a merge
  keeps the heavier body's ID and retires the other; a split gives each
  fragment a new run-time serial; a deleted ID is never reused. The
  allocators already exist (store._runtime_uids and
  store._next_body_numbers). A test shows an ejected planet keeps its
  ID. Decided (Boss, 2026-10-10 18:44Z, defaults approved; default yes):
  an ejected planet keeps its
  ID? Built with the orbital-update items that create these events (no
  ejection, merger or split code exists yet).
  Moved ahead of the collision item (child b of GEN.110) and the unbind
  item (GEN.291), which both apply these ID rules, and it no longer
  needs GEN.143: the allocators `store._runtime_uids` and
  `_next_body_numbers` exist, and its test can use a hand-built event
  (2026-10-16 review).
  Widened (2026-10-16 review): owns the writer `record_event(...)`,
  which also writes the admin activity-log line, and the ID helpers (an
  ejected planet keeps its ID and becomes a rogue-planet row; a merge
  keeps the heavier body's ID and retires the other; a split gives each
  fragment a new run-time serial via `store._runtime_uids` and
  `store._next_body_numbers`; a deleted ID is never reused). Needs
  DB.56. GEN.110 writes through it instead of defining its own event
  rows, and GEN.291 turns its "event log" half into calls to it.
  Prerequisites: DB.56. Related: GEN.105, GEN.143, GEN.170.

- [ ] **GEN.198 The gravity field evaluator: pull, potential and tidal gradient at any point from the point masses and the galaxy potential**
  Boss (2026-10-10 20:57Z): the orbits come from point-mass vectors, so
  the same masses give a gravitational field. Done: one function
  returns, for any point and time, the acceleration vector, the
  potential and the tidal tensor (its largest eigenvalue and the size of
  the gradient) from the point masses in and around a sector (the
  influence set of GEN.109) plus the galaxy's smooth potential
  (GEN.115), with Plummer softening; near sources are summed exactly and
  distant ones as per-region aggregates (PERF.72 sets the split). It is
  a pure function of the stored vectors, so the same galaxy gives the
  same field, and course planning (NAV.6) does not use it: it plans
  around keep-out spheres.
  Note (2026-10-10): Boss (22:14Z) asked that the size of each gradient
  measure be built from the preset minimum movement values of the
  orbital update system (THRESHOLDS_M in physics/position.py: galactic
  0.01 mpc, about 2 AU; system 0.01 AU; planetary 100,000 km). The
  gradient step of the evaluator is the threshold of the scale being
  drawn (galactic for a sector, system for the System Map, planetary for
  a moon neighbourhood), read from the same constants, not a copy. See
  the new section in docs/design/gravity-map.md.
  Prerequisites: GEN.338, GEN.342. Related: GEN.105, NAV.6,
  TEST.131,
  PERF.72, MAP.167.
  Scope after the review (2026-10-16 review): the tidal tensor and its
  eigenvalue, the galaxy potential term (GEN.115, GEN.285), the near and
  far split with per-region aggregates and the gradient step from
  THRESHOLDS_M; pull and potential at a point come from the field kernel
  (PERF.88). Course planning (NAV.6) reuses nothing from it: NAV.6,
  NAV.25 and NAV.26 plan around keep-out spheres.

- [ ] **GEN.199 Per-zone gravity grid for a sector: stored in a cache and rebuilt when the orbits update**
  Done: for a sector, the field (GEN.198) is sampled at the centre of
  each zone (16 per edge by default, a setting), giving the pull, the
  potential and the tidal strength per zone, with Lagrange and saddle
  points flagged. The grid is derived data held in a cache like the
  Galaxy Map tiles (not in the galaxy tables, so no galaxy migration),
  is built on the work queue, is thrown away when the orbital update
  (GEN.105) moves its sources, and a sector with no neighbours generated
  yet is marked uncharted-edge. A test builds a small sector and checks
  the grid against the evaluator zone by zone.
  Note (2026-10-10): Boss (22:14Z) asked that the size of each gradient
  measure be built from the preset minimum movement values of the
  orbital update system (THRESHOLDS_M in physics/position.py: galactic
  0.01 mpc, about 2 AU; system 0.01 AU; planetary 100,000 km). A zone
  edge is a whole number of those steps (the default sector zone of 0.25
  pc is 25,000 galactic steps) and never finer than one step; a cached
  zone is thrown away only when a source passes its own threshold, the
  same event that sets its next update due time (GEN.106). PERF.72 sets
  the multiples.
  Prerequisites: GEN.198, PERF.72. Related: GEN.105, PERF.72, MAP.167.

- [ ] **GEN.202 Scattered nebulae (planetary nebulae and remnants) join the nebula index, deduplicated by object id**
  Boss (2026-10-10 21:46Z): nebulae must be visible on the Galaxy Map,
  by colour or by combined regions; Phase 2 for the feature, foundations
  in Phase 1. From the research (docs/design/nebula-map-visibility.md
  section 9, PR #1180; report
  /mnt/project-files/research/nebula-map-visibility/report.md):
  planetary nebulae and supernova remnants from `phenomenon_scatter` go
  into the same nebula index as the field and stored nebulae,
  de-duplicated by object id, using the scatter index that PERF.68 adds.
  Done: one query gives every nebula kind in a box with no object twice,
  and a test shows it. Owner: Foundations lane 1.
  Prerequisite: DB.24. Related: GEN.176, PERF.71.
  Design: [docs/design/nebula-map-visibility.md](design/nebula-map-visibility.md)

- [ ] **GEN.209 Globular clusters are parent objects for their stars and are named when their sector is filled**
  GitHub issue [#1204](https://github.com/dwhagar/planetGen/issues/1204)
  (Boss, 2026-10-10 23:34Z): "Globular clusters will be a parent object
  for the stars within, the stars within the cluster. Clusters will get
  their own unique name when the sector is filled, not before." Done: a
  cluster is a stored object that its stars point to as their parent, it
  has no name until its first sector is filled and then takes one unique
  name, and the cluster pages and the Galaxy Map use the parent link.
  Owner: Foundations lane 2, with the globular cluster group (GEN.159 to
  GEN.164).
  Split (2026-10-16, Boss 1:07 PM EDT: look at each item for smaller
  chunks and shared work): this item is now the umbrella and is done
  when its subitems are: DB.61, GEN.359, MAP.224. GEN.209 needs a schema
  change that its text does not name: `galaxy_components.name` is NOT
  NULL and unique today (db/schema.sql, v97) while the cluster must have
  no name until its first sector is filled, and no column links a star
  to its cluster (the only parent column is `parent_component_id`
  between components).
  Scope after the review (2026-10-16 review): the migration for unnamed
  clusters and the star-to-cluster parent link (child a, Core 1), the
  sector fill linking cluster stars and naming the cluster at its first
  fill (child b), and the cluster pages and the Galaxy Map using the
  link (child c).
  Prerequisites: DB.61, GEN.359, MAP.224. Related: GEN.160, GEN.164.

- [ ] **GEN.212 Galactic move in the orbital update as NumPy batches (3 to 4 times on the whole update)**
  Math Research (2026-10-14;
  docs/design/math-and-library-efficiency.md).
  `db/store.advance_galactic_positions` spends 25.7 us a row in Python;
  a NumPy batch takes 3.8 us and refiles the same addresses on 100,000
  of 100,000 rows. The whole update gains 3 to 4 times because the SQL
  per row costs about as much. Also use the vector
  `solve_eccentric_anomalies` wherever updates are batched (the scalar
  `brentq` is 20.9 us against 1.36 us; production code does not call the
  vector one yet). Done: both folded into the GEN.105 chain with
  identical results. Owner: the Orbital Updates builders, with GEN.105.
  Narrowed (2026-10-16 review): batching the SQL write path
  (executemany, bulk updates) and calling the vector
  `solve_eccentric_anomalies` for batched comets; the NumPy move of the
  stars is the stepper (child c of GEN.115). Moved after that child so
  the loop is rewritten once.
  Prerequisites: GEN.340. Related: GEN.105, GEN.141, PERF.88.
  Design: [docs/design/math-and-library-efficiency.md](design/math-and-library-efficiency.md)

- [ ] **GEN.218 One habitability describer, in the short wording**
  Audit: /mnt/project-files/wording/text-output-audit-2026-10-15.md.
  Owner lane: Foundations Wording (2026-10-15), which writes the code
  for this item. Fold physics/protection.paragraph and
  physics/habitability_explain into one module that returns the rating
  rows (domain, colour, score, sentence) in Boss's wording from GEN.216 (done, PR #1482)
  (his example lines are the spec): one number formatter on the ladders,
  one gas-name table, PHI-6 everywhere, dead code deleted
  (habitability_explain.factors and helpers unless a row sentence needs
  them). Keeps GEN.148's temperature behaviour.
  Prerequisites: GEN.214 (done, PR #1472).

- [ ] **GEN.219 The habitability levels page and tooltips are generated from the limits**
  Audit: /mnt/project-files/wording/text-output-audit-2026-10-15.md.
  Owner lane: Foundations Wording (2026-10-15), which writes the code
  for this item. web/lib/habitability_ref.py reads the protection_limits
  and habitability constants (adds the chaotropicity, oxygen-high and
  heat-shielded rows; oxygen is listed under Pressure but scored in Air,
  fix the page to match the code), a test fails when a number on the
  page differs from the constant, and the page, chips, tooltips,
  class_pages.py and docs name PHI-6.
  Prerequisites: GEN.214 (done, PR #1472).

- [ ] **GEN.220 "Habitable" means the index**
  Audit: /mnt/project-files/wording/text-output-audit-2026-10-15.md.
  Owner lane: Foundations Wording (2026-10-15), which writes the code
  for this item. The system overview counts worlds by their PHI-6 rating
  (default: potentially habitable = overall Green or Blue; Boss can move
  the line), the row chip and search facet say what they mean (the
  rating, or "Habitable class" where it is the class flag), and the
  class pages label the flag "Habitable class". Touches UX
  (UX.32 chip decision, UX.101 symbols): coordinate before editing
  systempage.py.
  Prerequisites: GEN.218, UX.120. Related: UX.32, UX.101.

- [ ] **GEN.221 Short plain prose for the fixed text**
  Audit: /mnt/project-files/wording/text-output-audit-2026-10-15.md.
  Owner lane: Foundations Wording (2026-10-15), which writes the code
  for this item. One pass, with Boss's example lines as the model: the
  planet and moon paragraph (class description, atmosphere, composition;
  the "of a mix of" run-ons), the gas giant depth sentences (keep the
  data, drop the jargon), the star age and lifespan sentence and the
  STAR_EVOLUTION notes (split into a sentence about lifespan and one
  about life), comet, belt, wide and close binary sentences, the planet
  class descriptions that claim life or water (tuning.PLANET_CLASSES,
  class pages), and the phenomenon prose. In-universe wording per UX.42
  on public pages. Wording moves into the describer as it is rewritten.
  The phenomenon prose is not part of this item (2026-10-16 review):
  UX.108 writes it as it moves each kind into the Document.

- [ ] **GEN.229 Repair or re-roll the wide binaries that already break the stability rule**
  GEN.228 (PR #1610) makes new wide binaries stable and adds `planetgen
  check-db --deep --wide-binaries`, which counts the existing unstable
  ones. Boss (2026-10-15 07:44Z) saw two in his last galaxy whose
  orbital paths overlap and would end in planetary collisions. Run the
  count and report how many break the rule; the repair or re-roll itself
  needs Boss's word before anything is changed. Owner: Core. Source:
  Boss's report 2026-10-15 07:44Z.
  Decided (Boss, 2026-10-16, 3:00 PM EDT, decision card): run the count
  and report it first; the repair or re-roll needs Boss's word. Owner:
  Core 1, in its bug phase.
  Prerequisites: none. Related: GEN.228 (done, PR #1610).

- [ ] **GEN.231 Generate a whole block on the shell ladder: the generate --block mode and the block buttons come back**
  The Galaxy Map moved to the shell, band and wedge ladder
  (243/81/27/9/3 shells deep, then sectors; PR #1612,
  docs/design/spherical-block-ladder.md), and with it the old block
  addresses went away. `planetgen generate --block` is refused with "the
  Galaxy Map's blocks are being rebuilt on shells" (`cli/generate.py`,
  `tests/test_gen_cli_errors.py`) and the block generate buttons
  (`blockGenerateButtons`) were removed from `generatebuttons.js`. Done:
  a block of the new ladder (its m, shell, band and wedge) is a region
  the generate command and the work queue can fill, nearest sectors
  first, with the same progress, cancel and resume as a neighbourhood
  fill; `--block` takes the new block address; the Galaxy Map's block
  menu offers Generate again; the refusal text and its test are gone.
  Found by Maps (2026-10-15). Owner: Core for the region and the
  command; the button is Interface, after Core. Order it after PERF.99
  so the queue it uses is the new one.
  Finding (Quality, 2026-10-15): the `--block` refusal message in
  planetgen/cli/generate.py (about line 726) still says "the Galaxy
  Map's blocks are being rebuilt on shells" although the map is rebuilt;
  the option is still refused until it is ported to the
  m.shell.band.wedge key. This item reworded the message to say what is
  true until the mode is built.
  Split (2026-10-16, Boss 1:07 PM EDT: look at each item for smaller
  chunks and shared work): this item is now the umbrella and is done
  when its subitems are: GEN.310, UX.127. Its own note says 'Core, then
  Interface for the buttons': a server and command change and a browser
  change owned by two lanes.
  Prerequisites: GEN.310, UX.127. Related: GEN.101, MAP.120.

- [ ] **GEN.238 Wanderers: hypervelocity stars reach out to the grid, and kicked neutron stars and black holes are an ejecta population**
  Part of step 6 of the whole-sphere train. Hypervelocity stars keep
  their pre-placement (direction uniform on the sphere, distance uniform
  to the reach, density falling as r^-2, speed 500 to 1,000 km/s) but
  lose the outline test; the code's 1e3 to 1e4 total is read as "out to
  100 kpc", which puts 84% of them beyond the disk. Kicked neutron stars
  (20% to 40% of the isolated ones) and black holes (about 2%) are an
  `ejecta` population: density N / (4 pi r^2 v T) out to v T, speeds 300
  to 500 km/s above escape, random direction. They are not pre-placed
  (about 1e7 rows inside 262 kpc, 1e-8 per cell): they join the cell's
  `lambda` for neutron stars and black holes in the profile and are
  drawn at fill. Positions carry epoch and velocity, so the nightly
  update (GEN.223) moves them like everything else. Done: placement is a
  pure function of the seed, none is redrawn for landing outside the
  outline, and the counts per radius match the law. Owner: Core
  (2026-10-15), after GEN.232 (a fast body's path integration runs for
  minutes until it is capped).
  Built-on note (Core 1, 2026-10-16, 12:03 PM EDT, GEN.236 merged as
  #1824, no migration): SectorProfile.kind_factors carries one factor
  per rare kind (supernova remnant = gas fraction; molecular cloud =
  cold gas times the young share against the solar circle's; hosted
  nebulae = cold gas alone; rogue planet and comet = planet occurrence
  at the place's [Fe/H]; bound remnants, planetary nebulae and brown
  dwarfs = 1; no magnetar factor, GEN.278 gates it); the fill and the
  scatters apply them; metal-poor stars cap their object count at
  ceil(16 x occurrence) below [Fe/H] -0.5 and thin gas giants.
  Owner: Core 3 (Boss, 2026-10-16, 2:42 PM EDT: Core 3 takes the 8.1 new
  generate scatters line after GEN.280; the multi-star train follows its
  spine, 3:00 PM EDT). This replaces any earlier Core or Core 1 owner
  line on this item.
  Prerequisites: GEN.235, GEN.323. Related: GEN.223, GEN.240.
  Design: [docs/design/intergalactic-and-halo-contents.md](design/intergalactic-and-halo-contents.md)

- [ ] **GEN.239 High-velocity clouds and the Magellanic Stream as placed gas clouds (later)**
  Step 10 of the whole-sphere train. High-velocity clouds are the one
  gas item in the halo worth rows: a few hundred, placed like the
  molecular cloud field (`nebula_field`) at 2 to 15 kpc, and the
  Magellanic Stream at about 55 kpc and beyond. The discrete cloud count
  is not known from the searches done; get it before building. Done: the
  clouds exist as placed gas objects and the Stream is a chain of them.
  Later, after the components. Owner: Core (2026-10-15).
  Prerequisite: GEN.237. Related: GEN.202, GEN.241.
  Design: [docs/design/intergalactic-and-halo-contents.md](design/intergalactic-and-halo-contents.md)

- [ ] **GEN.240 The three scatters fill the whole sphere with one algorithm: band groups, the phenomena pass on the object-first sampler, and a component pass**
  Step 6 of the whole-sphere train. The bright-star and mass passes keep
  their algorithm (one Poisson count per shell from a certified bound,
  bands in proportion, accept by true density over bound;
  `object_first.py`) but take the halo term in `band_majorants`, the new
  populations in `fractions`, and band groups of 64 bands beyond the
  outline: to 262 kpc that is about 2 minutes a pass instead of 2.25
  hours, with under 10% efficiency loss beyond 8 kpc. The phenomena pass
  walks every band of every shell today (`_band_rows`), which cannot
  reach 262 kpc; it moves to the same object-first sampler, with its
  rates, cut and below-cut draw at sector time unchanged. Components
  (clusters, Clouds, Sagittarius) get their own pass: total star count
  from the mass, positions by inverse transform of the profile table,
  each star kept or dropped by the same pass cut, and the sector fill
  adds the component's density to `lambda` so no pre-placed star is
  counted twice. At the defaults the passes place nothing in the halo,
  clusters or old satellites, a few thousand stars in the Clouds, and
  about 2.4 million extra black-hole rows outside the disk (4% on top of
  the disk's scatter). Done: the three passes cover the reach with the
  same code, a pre-placed star is never drawn again at fill, and the
  yields match the design note's table. Owner: Core (2026-10-15), after
  GEN.232 (the cap on path integration for fast bodies); the wanderers
  are GEN.238. Decided (Boss, 2026-10-15): the research note's
  recommendation stands: reach 262 kpc (65,535 shells), 96-bit object
  IDs, one rule set everywhere with per-sector parameters, nothing
  generated unless asked, and an asked sector always generated and
  saved.
  Decided (Boss, 2026-10-15 5:27 PM EDT): the scatter from the origin
  keeps the banded, layered approach it has now, because the bands of a
  shell have fixed slices (Boss's recollection; Core checks it against
  `object_first.py` and the band and wedge counts of `spherical.py`
  before relying on it, and reports if it is wrong). The band groups of
  this item stay inside that structure. Decided (Boss, 2026-10-15 5:27
  PM EDT): early alpha, so the current database schema is not a
  constraint. Redesign tables, columns and layouts as the work needs,
  with one migration or a fresh galaxy; nothing is kept for the sake of
  production data. Decided (Boss, 2026-10-15 5:39 PM EDT): scatters and
  backscatters get no systems, planets or moons, so bright stars,
  massive stars, neutron stars, black holes and the other scattered
  objects are bare rows; planets and moons are made only when a sector
  is being filled with everything; a sector fill must represent every
  kind of object, and when a fill meets a phenomenon that can have a
  system (or makes one that can) it generates that system if the
  real-world probabilities give it one; nebulae need stars and
  phenomena, not planets. For this item: the three scatters write bare
  rows and never build a system; the system is built when the sector
  fills.
  Decided (Boss, 2026-10-15 5:51 PM EDT), refining the 5:39 PM rule:
  scatters still carry no planets or moons, but the stellar companions
  of stars and star-like objects (stars, former stars, protostars,
  neutron stars, black holes), including neutron star pairs, black hole
  and neutron star pairs, black hole pairs and stars with neutron star
  or black hole companions, are generated at scatter time, because such
  systems matter for gravity; planets and moons still come only at
  sector fill. Multi-star research is under way (docs to follow; the
  coordinator relays its items). For this item: each scatter row carries
  its stellar companions (their masses and orbit), but never planets or
  moons; the band groups and the object-first sampler must place the
  companions with the primary so a pair is one gravitational source.
  Reworded (Boss, 2026-10-15 5:51 PM EDT;
  docs/design/multi-star-systems.md): companions of stars and star-like
  objects are generated at scatter time, planets, moons and belts at
  sector fill; a scatter row is a system with a barycentre, a total mass
  and stored members. Each scatter row's members are drawn in the
  scatter passes (GEN.262), the band groups and the object-first sampler
  must treat a row as one system (`1 + p_placed` bookkeeping), and the
  whole-sphere passes take the same member columns (GEN.260).
  Region query (stage C): scatter rows keep their (shell, band, wedge)
  address columns and gain the point columns and cell key (DB.35); the
  scatter pass keeps its bands and order and the key is computed from
  whatever position it writes.
  Prerequisites: GEN.322, GEN.323, GEN.324. Related: GEN.238, GEN.223.
  Built-on note (Core 1, 2026-10-16, 12:03 PM EDT, GEN.236 merged as
  #1824, no migration): SectorProfile.kind_factors carries one factor
  per rare kind (supernova remnant = gas fraction; molecular cloud =
  cold gas times the young share against the solar circle's; hosted
  nebulae = cold gas alone; rogue planet and comet = planet occurrence
  at the place's [Fe/H]; bound remnants, planetary nebulae and brown
  dwarfs = 1; no magnetar factor, GEN.278 gates it); the fill and the
  scatters apply them; metal-poor stars cap their object count at
  ceil(16 x occurrence) below [Fe/H] -0.5 and thin gas giants.
  Split (2026-10-16, Boss 1:07 PM EDT: look at each item for smaller
  chunks and shared work): this item is now the umbrella and is done
  when its subitems are: GEN.322, GEN.323, GEN.324. The item itself says
  'the three scatters' and lists three deliverables: band groups and the
  halo term for the bright-star and mass passes (about 2 minutes a pass
  instead of 2.25 hours), moving the phenomena pass off its walk of
  every band (_band_rows) onto the object-first sampler, and a new
  component pass.
  Owner: Core 3 (Boss, 2026-10-16, 2:42 PM EDT: Core 3 takes the 8.1 new
  generate scatters line after GEN.280; the multi-star train follows its
  spine, 3:00 PM EDT). This replaces any earlier Core or Core 1 owner
  line on this item.
  Design: [docs/design/intergalactic-and-halo-contents.md](design/intergalactic-and-halo-contents.md)

- [ ] **GEN.241 Every nebula gets the stars it needs as real systems made with the nebula, inside its bounds**
  Step 12 of the whole-sphere train (design note section 5.10); replaces
  the text-only hosts and the "do not build" default of MAP.182, which
  is folded in here (hosted nebulae for stars pre-placed by the bright
  and mass passes: the hosted roll is decided from the star's own seed,
  so it is the same whether the star was pre-placed or made during a
  fill, and the page and the map show the nebula from the scatter row).
  The rule: a nebula's required stars are real star systems made in the
  same step as the nebula and placed inside its bounds, and which kinds
  it needs depends on its class, not on where it is. A nebula that
  cannot get its stars (the sector is full) is not made. Per class:
  emission C the one O or early-B star at the centre; D and E the seed
  star plus the cluster's ionizing members (stars of 5.5 solar masses
  and up that are still alive: 42 on average for D, 549 for E) in a
  Plummer sphere of half-mass radius 1 pc (D) and 2 pc (E) clipped to
  the nebula radius; emission-reflection G the early-B star; reflection
  F the B3 to A9 star; planetary nebulae keep their white dwarf and gain
  the close companion the class table says is "often" there (20%);
  supernova remnants keep their neutron star or black hole at 60% of
  core-collapse remnants; Bok globule P one protostar in 30% of
  globules; star-forming core Q 1 to 4 protostars (Poisson mean 2, at
  least 1) clumped near the centre; giant molecular cloud M and dark
  cloud N none. Protostars need the pre-main-sequence state (GEN.242).
  Cluster members and protostars are derived from the nebula's own seed
  and placed by position, and each sector generates only the members
  whose positions fall in its cube, the way `nebula_field` hands each
  sector its clouds, so a nebula across sectors needs no coordination.
  Any code path that creates a nebula of a hosted class without a host
  calls one function, `build_nebula_with_contents`. The Type Ia remnant
  share becomes a profile field (1 where there is no young population);
  emission, reflection, Wolf-Rayet and core-collapse remnants scale with
  the young share of the profile's mix, planetary nebulae and Type Ia
  remnants with the old share, molecular clouds and protostars with
  `gas`; `nebula_field` must read the cell's `gas`, not the disk
  position alone. Open question for Boss (default: only the ionizing
  members are real, the rest of the cluster stays the sector's ordinary
  background density): filling the whole cluster is the dearer
  alternative, about 8 minutes of generation per giant complex, behind a
  setting. Cost at the default, if every sector were filled: about 7e6
  systems; an asked sector holding a massive complex gets about 270
  extra systems and 7 seconds. Done: an asked halo cell never shows an H
  II region or a pulsar-bearing remnant, a planetary nebula always has
  its white dwarf, every hosted nebula holds its stars inside its
  radius, a Magellanic cell can show all of them, and a pre-placed
  star's nebula shows from its scatter row. Owner: Core (2026-10-15);
  the map shows the nebulae through the existing nebula items (MAP.172
  to MAP.181). Decided (Boss, 2026-10-15): the research note's
  recommendation stands: reach 262 kpc (65,535 shells), 96-bit object
  IDs, one rule set everywhere with per-sector parameters, nothing
  generated unless asked, and an asked sector always generated and
  saved.
  Decided (Boss, 2026-10-15 5:39 PM EDT): scatters and backscatters get
  no systems, planets or moons, so bright stars, massive stars, neutron
  stars, black holes and the other scattered objects are bare rows;
  planets and moons are made only when a sector is being filled with
  everything; a sector fill must represent every kind of object, and
  when a fill meets a phenomenon that can have a system (or makes one
  that can) it generates that system if the real-world probabilities
  give it one; nebulae need stars and phenomena, not planets. For this
  item: the stars and phenomena a nebula needs are real objects made
  with it; they get planets and moons only when the nebula's sector is
  being filled, never in a scatter.
  Catalogue (Boss 2026-10-15 6:08 PM EDT;
  docs/design/multi-star-systems.md sections 13 to 16): the young
  non-hierarchical groups and trapezia of section 13.1 are built by
  GEN.273 on this item's cores.
  Prerequisites: GEN.319, GEN.320, GEN.321. Related: MAP.172, GEN.150,
  GEN.202, GEN.99. Built-on note (Core 1, 2026-10-16, 8:41 AM EDT,
  GEN.251 merged as #1815): the stage F required stars (the stars a
  nebula needs as real systems) are not built by GEN.251; fill stages
  have their own random streams (generation/fill_stages.py), so this
  item adds its stage there.
  Split (2026-10-16, Boss 1:07 PM EDT: look at each item for smaller
  chunks and shared work): this item is now the umbrella and is done
  when its subitems are: GEN.319, GEN.320, GEN.321. Step 12 mixes three
  deliverables of different size.
  Owner: Core 3 (Boss, 2026-10-16, 2:42 PM EDT: Core 3 takes the 8.1 new
  generate scatters line after GEN.280; the multi-star train follows its
  spine, 3:00 PM EDT). This replaces any earlier Core or Core 1 owner
  line on this item.
  Design: [docs/design/intergalactic-and-halo-contents.md](design/intergalactic-and-halo-contents.md)

- [ ] **GEN.242 A pre-main-sequence state for young stars, so protostars can exist as real stars**
  Step 11 of the whole-sphere train (design note section 5.10). The star
  model has no pre-main-sequence state: `evolve_star` goes from the main
  sequence on. Add one (age under 1 Myr for stars below 8 solar masses,
  drawn from the IMF's low end for protostars), shown as such on the
  star page and in the description layer (wording through the Interface
  lane's document model, UX.105 and UX.108). Done: a star younger than
  the limit comes out of `evolve_star` as a pre-main-sequence star with
  its own radius, luminosity and type text, the older states are
  unchanged, and the page says so. Owner: Core (2026-10-15). Decided
  (Boss, 2026-10-15): the research note's recommendation stands: reach
  262 kpc (65,535 shells), 96-bit object IDs, one rule set everywhere
  with per-sector parameters, nothing generated unless asked, and an
  asked sector always generated and saved.
  Systems around protostars (design note section 5.12, answering Boss's
  2026-10-15 question). A protostar is a state of an ordinary star, so
  it gets an ordinary star system with no special path. No planets yet:
  the planet draw already drops every non-belt body around a star
  younger than `PLANET_MIN_STAR_AGE_GY` = 0.01 Gy (`system.py:1024`), so
  belts are the most it keeps and nothing new is needed for that. Add a
  `disk` descriptor on the star instead (class 0, I or II by age under
  0.15, 0.5 and 3 Myr; mass 0.01 to 0.3 of the star, outer radius 10 to
  300 AU), with no bodies in it. Companions come from the normal binary
  machinery at the Perseus multiplicity: 0.57 for Class 0 and 0.23 for
  Class I, higher than the field's 0.26 to 0.45. Done also means: a
  protostar's system shows its disk and its companions, and holds no
  planets.
  Decided (Boss, 2026-10-15 5:39 PM EDT): scatters and backscatters get
  no systems, planets or moons, so bright stars, massive stars, neutron
  stars, black holes and the other scattered objects are bare rows;
  planets and moons are made only when a sector is being filled with
  everything; a sector fill must represent every kind of object, and
  when a fill meets a phenomenon that can have a system (or makes one
  that can) it generates that system if the real-world probabilities
  give it one; nebulae need stars and phenomena, not planets. For this
  item: a protostar's disk, companions and belts come when its sector
  fills; the scattered row carries none.
  Catalogue (Boss 2026-10-15 6:08 PM EDT;
  docs/design/multi-star-systems.md sections 13 to 16): the
  circumstellar and circumbinary disk descriptor and the young
  non-hierarchical groups of sections 13.1 and 13.5 build on this state
  (GEN.273, GEN.275).
  Finding (Core 2, 2026-10-15, #1705): the design note's
  pre-main-sequence time 31 Myr x m^-2.5 gives 10 Gyr at 0.1 solar
  masses and 1.2 Gyr at 0.2 (models: about 1 Gyr and a few hundred Myr).
  GEN.261 capped it at 1 Gyr (KELVIN_HELMHOLTZ_MAX_YR in
  generation/multistar.py); this item confirms or replaces the law.
  Scope after the review (2026-10-16 review): The pre-main-sequence
  state and the disk laws (class by age, mass and radius); it fills the
  disk descriptor columns that the first child of GEN.129 adds, and does
  not add them.
  Owner: Core 3 (Boss, 2026-10-16, 2:42 PM EDT: Core 3 takes the 8.1 new
  generate scatters line after GEN.280; the multi-star train follows its
  spine, 3:00 PM EDT). This replaces any earlier Core or Core 1 owner
  line on this item.
  Prerequisites: GEN.235, DB.53. Related: GEN.241, UX.108.
  Design: [docs/design/intergalactic-and-halo-contents.md](design/intergalactic-and-halo-contents.md)

- [ ] **GEN.243 Wolf-Rayet stars and ring nebulae: a W state for the last 0.3 Myr of very massive stars, and nebula class X**
  Step 13 of the whole-sphere train (design note section 5.11; GEN.113's
  Tier 2 item 3, never built because the star model has no W state). The
  star: a star of initial mass 25 solar masses or more (solar
  metallicity) is Wolf-Rayet during the last 0.3 Myr of its life
  (`WR_MIN_INITIAL_MASS_SOL = 25`, `WR_PHASE_MYR`, at most half the life
  for the heaviest) instead of a supergiant, at 40% to 60% of its
  initial mass, 30,000 to 200,000 K and 1e5 to 1e6 solar luminosities;
  subtype WN for the first 60% of the phase, WC for the next 30%, WO for
  the last 10% (the heaviest only), as spectral text `WN6`, `WC8` and
  Yerkes class `W`. The threshold rises as metallicity falls (about 35
  at half solar, 56 at a fifth), so it is a population field in the
  profile, not a constant. The bright and mass passes already place the
  massive stars and then draw the state through `evolve_star`, so
  Wolf-Rayet stars come out of them without a new rate; none exists in
  any population older than about 10 Myr (halo, clusters, old
  satellites, bulge). The nebula: new class `X` in `NEBULA_CLASSES`,
  family emission, radius 2 to 35 pc, a thin shell with a hollow
  interior on the existing `nebula_shape` mesh, host rule `("W", 0, 9,
  ("X",), 0.4)` in `NEBULA_HOST_RULES` and
  `NEBULA_CLASS_MAX_AGE_MYR["X"] = 5.0`; its one required star is the W
  star at the centre, a real system (GEN.241), and after the phase the
  star becomes the core-collapse remnants the generator already makes.
  The Wolf-Rayet row of `anomalies.md` becomes a pointer here. Done: the
  Milky Way disk holds 1,200 +- 300 Wolf-Rayet stars at the default
  scale, none in an old population, every class X nebula holds a W star,
  and the star and nebula pages and the map show them. Owner: Core
  (2026-10-15); the map marker (MAP.132) is Interface's. Decided (Boss,
  2026-10-15): the research note's recommendation stands: reach 262 kpc
  (65,535 shells), 96-bit object IDs, one rule set everywhere with
  per-sector parameters, nothing generated unless asked, and an asked
  sector always generated and saved.
  Systems around Wolf-Rayet stars (design note section 5.12, answering
  Boss's 2026-10-15 question). A W star is a state of an ordinary star,
  so it gets an ordinary system with no special path. No planets (the
  same age cut as above; no planet is known around one), belts allowed.
  Companions come from the normal binary machinery with the W star's own
  fractions: 50% for WN and 90% for WC (the surveyed 0.96 capped), the
  companion an O or B star or a compact object, the period drawn from
  the WN law (mostly short, like O stars) or the WC law (peak near 5,000
  days). After GEN.129 the same fractions apply to trees. A planetary
  nebula's central white dwarf is already an ordinary system, so planets
  there are possible and engulfment removes the inner ones. Done also
  means: a W star's system carries its companions at those fractions and
  no planets. Decided (Boss, 2026-10-15 5:39 PM EDT): scatters and
  backscatters get no systems, planets or moons, so bright stars,
  massive stars, neutron stars, black holes and the other scattered
  objects are bare rows; planets and moons are made only when a sector
  is being filled with everything; a sector fill must represent every
  kind of object, and when a fill meets a phenomenon that can have a
  system (or makes one that can) it generates that system if the
  real-world probabilities give it one; nebulae need stars and
  phenomena, not planets. For this item: a Wolf-Rayet star's companions
  and belts come when its sector fills; the pre-placed or scattered row
  carries none.
  Catalogue (Boss 2026-10-15 6:08 PM EDT;
  docs/design/multi-star-systems.md sections 13 to 16): the Wolf-Rayet
  binary row of section 13.2 (WN 52%, WC 90%) and the kind `wolf-rayet`
  in the star model are built by GEN.269 and GEN.270 on this item's W
  state.
  Owns the Wolf-Rayet companion fractions (the WN and WC companion
  rates); GEN.270 drops its Wolf-Rayet row (2026-10-16 review).
  Owner: Core 3 (Boss, 2026-10-16, 2:42 PM EDT: Core 3 takes the 8.1 new
  generate scatters line after GEN.280; the multi-star train follows its
  spine, 3:00 PM EDT). This replaces any earlier Core or Core 1 owner
  line on this item.
  Prerequisites: GEN.320. Related: GEN.113, MAP.132, MAP.176.
  Design: [docs/design/intergalactic-and-halo-contents.md](design/intergalactic-and-halo-contents.md)

- [ ] **GEN.244 Young open clusters and OB associations as components (later)**
  Design note 5.12. A component of the same table as the globular
  clusters (GEN.237): a Plummer profile, age under 1 Gyr, disk-bound,
  with its own age window and [Fe/H]; the hosted nebulae of GEN.241 are
  its youngest form. Done: a cell inside a young cluster holds its
  members at the cluster's age and the map can show the cluster as a
  point or a region. Later: not part of the whole-sphere train and
  nothing in it waits on this; plan after it. Owner: Core, later
  (2026-10-15).
  Catalogue (Boss 2026-10-15 6:08 PM EDT;
  docs/design/multi-star-systems.md sections 13 to 16): OB associations
  hold the trapezia of GEN.273 and the runaway origins of GEN.274.
  Owner: Core 3 (Boss, 2026-10-16 12:16 AM EDT, goal issue #1677;
  replaces any earlier owner named above). Core 1 still owns every
  migration revision: ask it before adding a table or column.
  Scatter note (Core 3, 2026-10-16,
  docs/design/later-phenomena-scatter-rules.md (section 8)): the
  component feeds `lambda` and `mix` only; no per-star scatter rows; the
  backfill reads the component term; trapezia and runaway origins come
  from the member tree.
  Owner: Core 3 (Boss, 2026-10-16, 2:42 PM EDT: Core 3 takes the 8.1 new
  generate scatters line after GEN.280; the multi-star train follows its
  spine, 3:00 PM EDT). This replaces any earlier Core or Core 1 owner
  line on this item.
  Prerequisites: GEN.241, GEN.321. Related: GEN.158.
  Design: [docs/design/intergalactic-and-halo-contents.md](design/intergalactic-and-halo-contents.md)

- [ ] **GEN.245 Debris disks around A and F stars (later)**
  Design note 5.12. A debris-disk flag (Vega type) on A and F stars as a
  descriptor on the star, no new rows; the protoplanetary `disk`
  descriptor of young stars belongs to GEN.242. Done: a share of A and F
  stars show a debris disk on the star page, with a radius and a dust
  mass. Later: not part of the whole-sphere train and nothing in it
  waits on this; plan after it. Owner: Core, later (2026-10-15).
  Owner: Core 3 (Boss, 2026-10-16 12:16 AM EDT, goal issue #1677;
  replaces any earlier owner named above). Core 1 still owns every
  migration revision: ask it before adding a table or column.
  Scatter note (Core 3, 2026-10-16,
  docs/design/later-phenomena-scatter-rules.md (section 8)): fill-only
  attribute; it shares the `stars` disk columns with GEN.242.
  Prerequisite: GEN.242.
  Design: [docs/design/intergalactic-and-halo-contents.md](design/intergalactic-and-halo-contents.md)

- [ ] **GEN.246 Luminous blue variables and hypergiants (later)**
  Design note 5.12. States of the massive-star model like Wolf-Rayet (a
  short window before it, with an ejecta nebula); rare, tens in the
  Galaxy. Done: a very massive star passes through the state before it
  becomes a Wolf-Rayet star or explodes, at the real count, with an
  ejecta nebula. Later: not part of the whole-sphere train and nothing
  in it waits on this; plan after it. Owner: Core, later (2026-10-15).
  Owner: Core 3 (Boss, 2026-10-16 12:16 AM EDT, goal issue #1677;
  replaces any earlier owner named above). Core 1 still owns every
  migration revision: ask it before adding a table or column.
  Scatter note (Core 3, 2026-10-16,
  docs/design/later-phenomena-scatter-rules.md (section 8)): states of
  stars the passes place; the ejecta nebula is hosted at fill; no
  scatter work.
  Prerequisite: GEN.243. Related: GEN.241.
  Design: [docs/design/intergalactic-and-halo-contents.md](design/intergalactic-and-halo-contents.md)

- [ ] **GEN.247 Accreting white dwarf binaries: novae, cataclysmic and symbiotic variables, Type Ia progenitors (later)**
  Design note 5.12. Subtypes of white-dwarf-plus-companion systems once
  the tree and exotic-system rules exist (the symbiotic one is already
  in `anomalies.md`); they follow the old-star share, so they appear in
  the halo and in clusters. X-ray binaries, microquasars and magnetars
  are already designed in `anomalies.md` and
  `multistar-and-compact-systems.md` 7.3 and follow the young share.
  Done: a white dwarf with a close companion can be a nova or
  cataclysmic variable at its rate, and the Type Ia remnant share of
  GEN.241 reads from them. Later: not part of the whole-sphere train and
  nothing in it waits on this; plan after it. Owner: Core, later
  (2026-10-15).
  Multi-star catalogue (Boss 2026-10-15 6:08 PM EDT): GEN.271 builds
  cataclysmic variables, AM CVn, symbiotic and barium stars and GEN.272
  the accreting compact binaries, so this item keeps only what they
  leave (novae as events stay out of scope). Decide when GEN.271 lands
  whether this item is closed.
  Scope after the review (2026-10-16 review): only the Type Ia
  progenitor share, read from the accreting pairs GEN.271 builds; novae
  are events and stay out of scope. Decide when GEN.271 lands whether
  this item closes.
  Prerequisites become the cataclysmic and symbiotic rows (GEN.271,
  GEN.272): the Type Ia progenitor share is read from them, and if it
  can be read without a new field this item closes (2026-10-16 review).
  Prerequisites: GEN.271, GEN.272, GEN.318. Related: GEN.241.
  Design: [docs/design/intergalactic-and-halo-contents.md](design/intergalactic-and-halo-contents.md)

- [ ] **GEN.248 Tidal streams as placed components of the halo (later)**
  Design note 5.12. Streams (Sagittarius and the thin ones) as an
  elongated density term in the component table; each is too small to
  alter a cell on its own. Done: a cell on a stream expects the extra
  stars and the component table lists the stream. Later: not part of the
  whole-sphere train and nothing in it waits on this; plan after it.
  Owner: Core, later (2026-10-15).
  Prerequisite: GEN.237. Related: GEN.239.
  Design: [docs/design/intergalactic-and-halo-contents.md](design/intergalactic-and-halo-contents.md)

- [ ] **GEN.249 Be stars, carbon stars and Mira variables as star attributes (later)**
  Design note 5.12. State or attribute on an ordinary star (decretion
  disks, carbon stars, Mira variables); no new kind, no new rows. Done:
  each attribute appears on the star at its rate and the star page says
  so. Later: not part of the whole-sphere train and nothing in it waits
  on this; plan after it. Owner: Core, later (2026-10-15).
  Multi-star catalogue (Boss 2026-10-15 6:08 PM EDT): GEN.270 builds Be
  stars (about 20% of B stars) with their `be` flag and decretion disk;
  this item keeps carbon stars and Mira variables. Decide when GEN.270
  lands whether the Be part is closed here.
  Scope after the review (2026-10-16 review): carbon stars and Mira
  variables as star attributes, reusing the `flags` column of the
  multi-star migration; the Be part is GEN.270's. Needs GEN.270 so the
  attribute pattern is built once.
  Prerequisites: GEN.270.
  Design: [docs/design/intergalactic-and-halo-contents.md](design/intergalactic-and-halo-contents.md)

- [ ] **GEN.252 Brown dwarfs as stars: an L, T and Y state for companions and primaries, with the brown-dwarf desert**
  Step 5 of the phenomena train (catalogue section 8 item 5). An L/T/Y
  state from the `ROGUE_BROWN_DWARF_*` cooling tables, so a brown dwarf
  can be a companion or primary of a system, with the brown-dwarf-desert
  rule for close companions. Done: a system can hold a brown dwarf star
  at the observed rate, shown as one on its page. Later; nothing in the
  whole-sphere train waits on it. Owner: Core, later (2026-10-15).
  Catalogue (Boss 2026-10-15 6:08 PM EDT;
  docs/design/multi-star-systems.md sections 13 to 16): the brown dwarf
  rows of section 13.5 (close companion under 0.6% inside 3 to 10 AU,
  wide 3.2%, pairs 15% of brown dwarfs, planetary-mass companions) are
  built by GEN.275.
  Owns the L, T and Y brown-dwarf star state and nothing else; the
  brown-dwarf companion rows and rates are child a of GEN.275
  (2026-10-16 review).
  Owner: Core 3 (Boss, 2026-10-16, 2:42 PM EDT: Core 3 takes the 8.1 new
  generate scatters line after GEN.280; the multi-star train follows its
  spine, 3:00 PM EDT). This replaces any earlier Core or Core 1 owner
  line on this item.
  Owner: Core 4 (handed over from Core 3 by Boss's "Split it", 2026-10-16,
  4:49 PM EDT; Core 3 agreed 4:53 PM EDT; ready: no blocker). Core 4 asks Core 3
  before touching stellar evolution constants in tuning.py. A line here
  that named Core 3 or Core 4 earlier is replaced by this one.
  Prerequisites: none. Related: GEN.129.
  Design: [docs/design/phenomena-catalogue-and-fill.md](design/phenomena-catalogue-and-fill.md)

- [ ] **GEN.253 White-dwarf pollution flag and cool white dwarf presets**
  Step 6 of the phenomena train (catalogue section 8 item 6). 25% to 50%
  of cool white dwarfs are polluted by accreted rocky debris [S]; add
  the flag and the presets, part of the white-dwarf row of GEN.130.
  Done: a cool white dwarf shows a pollution flag at that share and its
  system presets are the cool-white-dwarf ones. Later. Owner: Core,
  later (2026-10-15).
  Catalogue (Boss 2026-10-15 6:08 PM EDT;
  docs/design/multi-star-systems.md sections 13 to 16): the white dwarf
  rows of section 13.3 (close double white dwarfs, post-common-envelope
  pairs, cataclysmic variables, AM CVn, symbiotic and barium stars) are
  built by GEN.271, which extends this item.
  Moved to phase 3 of 8.1, before GEN.271, which extends it (the
  white-dwarf pollution flag and cool white-dwarf presets); it is small:
  a flag and presets in the white-dwarf row of GEN.130 (2026-10-16
  review).
  Owner: Core 3 (Boss, 2026-10-16, 2:42 PM EDT: Core 3 takes the 8.1 new
  generate scatters line after GEN.280; the multi-star train follows its
  spine, 3:00 PM EDT). This replaces any earlier Core or Core 1 owner
  line on this item.
  Owner: Core 4 (handed over from Core 3 by Boss's "Split it", 2026-10-16,
  4:49 PM EDT; Core 3 agreed 4:53 PM EDT; waits on GEN.318 (Core 3)). Core 4 asks Core 3
  before touching stellar evolution constants in tuning.py. A line here
  that named Core 3 or Core 4 earlier is replaced by this one.
  Prerequisites: GEN.318.
  Design: [docs/design/phenomena-catalogue-and-fill.md](design/phenomena-catalogue-and-fill.md)

- [ ] **GEN.255 Drawn jets: Herbig-Haro flows, microquasar jets and pulsar jets along their spin axes**
  Split from GEN.114 (catalogue section 8 item 7). GEN.104 gave every
  rotating object a spin axis; draw the jets of protostars in class Q
  cores (up to 1.5e5 in the Galaxy), microquasars and pulsars along it,
  on the object pages and the maps. Done: each of the three kinds shows
  a jet along its spin axis at its own rate. Later. Owner: Core then
  Interface, later (2026-10-15).
  Owner: Core 3 (Boss, 2026-10-16 12:16 AM EDT, goal issue #1677;
  replaces any earlier owner named above). Core 1 still owns every
  migration revision: ask it before adding a table or column.
  Scatter note (Core 3, 2026-10-16,
  docs/design/later-phenomena-scatter-rules.md (section 8)): derived at
  render; microquasar jets read the stored accreting member.
  Prerequisite: GEN.242. Related: GEN.104, GEN.247.
  Design: [docs/design/phenomena-catalogue-and-fill.md](design/phenomena-catalogue-and-fill.md)

- [ ] **GEN.256 A radiation hazard score derived from distance to hot stars, Wolf-Rayet stars, remnants and magnetars**
  Split from GEN.114 (catalogue section 8 item 7). A derived value, not
  an object: the dose at a point from the O and B stars, Wolf-Rayet
  stars, remnants and active magnetars near it (a giant flare of 1e46
  erg delivers 8.4e7 erg/cm2 at 1 pc), feeding the habitability
  radiation term; a flag within 10 pc of an active magnetar. Done: a
  system page shows its score and the terms behind it. Later. Owner:
  Core, later (2026-10-15).
  Owner: Core 3 (Boss, 2026-10-16 12:16 AM EDT, goal issue #1677;
  replaces any earlier owner named above). Core 1 still owns every
  migration revision: ask it before adding a table or column.
  Scatter note (Core 3, 2026-10-16,
  docs/design/later-phenomena-scatter-rules.md (section 8)): reads
  scatter rows through the stage C region query; it needs the magnetar
  rows (GEN.278) to see the 10 pc flag before the fill.
  Prerequisite: GEN.243. Related: GEN.113.
  Design: [docs/design/phenomena-catalogue-and-fill.md](design/phenomena-catalogue-and-fill.md)

- [ ] **GEN.257 A Fermi-bubble overlay on the Galaxy Map**
  Split from GEN.114 (catalogue section 8 item 7). The two lobes of
  about 50 degrees and 8 kpc above and below the nucleus, drawn as an
  overlay on the Galaxy Map, no objects. Done: the overlay can be shown
  and hidden and sits where the data puts it. Later. Owner: Interface,
  later (2026-10-15).
  Owner: Core 3 (Boss, 2026-10-16 12:16 AM EDT, goal issue #1677;
  replaces any earlier owner named above). Core 1 still owns every
  migration revision: ask it before adding a table or column.
  Prerequisites: GEN.233, MAP.217. Related: GEN.113.
  Design: [docs/design/phenomena-catalogue-and-fill.md](design/phenomena-catalogue-and-fill.md)

- [ ] **GEN.258 Protostars embedded in giant molecular clouds**
  Step 8 of the phenomena train (catalogue section 8 item 8). The class
  P and Q cores give 2.6e5 protostars against 1.4e6 expected; embedded
  Q-type groups inside M clouds make up the rest. Default: accept the 5x
  shortfall for the first release, as GEN.241 and GEN.242 stand. Done: a
  giant molecular cloud holds embedded Q-type groups of protostars and
  the count reaches the expectation, or Boss accepts the shortfall.
  Later. Owner: Core, later (2026-10-15).
  Owner: Core 3 (Boss, 2026-10-16, 2:42 PM EDT: Core 3 takes the 8.1 new
  generate scatters line after GEN.280; the multi-star train follows its
  spine, 3:00 PM EDT). This replaces any earlier Core or Core 1 owner
  line on this item.
  Owner: Core 4, "Core 4 (2026-10-16)" (Boss approved the lane at 4:31
  PM EDT, 2026-10-16; queue in the lane report); held until its blocker
  merges: GEN.242, GEN.320 (Core 3). Core 4 takes it when it unblocks.
  Prerequisites: GEN.242, GEN.320. Related: GEN.242.
  Design: [docs/design/phenomena-catalogue-and-fill.md](design/phenomena-catalogue-and-fill.md)

- [ ] **GEN.259 Intermediate-mass black hole share gives about 100,000 in the Galaxy against a handful of candidates (bug)**
  `BLACK_HOLE_INTERMEDIATE_MASS_CHANCE` = 0.1% of black holes gives
  about 1e5 intermediate-mass holes in the Galaxy, four orders above the
  few candidates and none confirmed in the field (catalogue section
  4.2). It was kept on purpose (Boss's "smattering") and is left alone
  until Boss says. Decided (Boss, 2026-10-16, 3:00 PM EDT, decision
  card): cut the chance so the Galaxy holds hundreds, as the
  observations suggest. Done: the chance is cut so the Galaxy holds
  hundreds of intermediate-mass holes, and its tests change with it.
  Owner:
  Core (2026-10-15).
  Owner: Core 3, in the 8.1 new generate scatters line (Boss,
  2026-10-16, 2:42 PM EDT; the cut was approved at 3:00 PM EDT).
  Prerequisites: none. Related: GEN.236.
  Design: [docs/design/phenomena-catalogue-and-fill.md](design/phenomena-catalogue-and-fill.md)

- [ ] **GEN.263 Promote every bound compact pair to the scatter, whatever its mass**
  Step 4 (design note 4.5 and 11 item 4). A neutron star never passes
  the 14 solar mass cut, so a double neutron star would be a fill-time
  object. Scatter every bound compact pair (neutron star with neutron
  star, with black hole, black hole with black hole) whatever its mass,
  about 3e5 extra rows, and remove its members from the fill's below-cut
  draws. Default on (Boss may switch it off). Done: neutron star pairs
  exist as scatter rows with their members and the below-cut draw no
  longer makes them twice. Owner: Core (2026-10-15). Boss (2026-10-15
  5:51 PM EDT): stellar companions of stars and star-like objects,
  including neutron star and black hole pairs, are generated at scatter
  time because they matter for gravity; planets, moons and belts stay at
  sector fill. Defaults taken (reversible): compact pairs promoted to
  the scatter, black hole companion fraction 1.25%, cap 12 stars (7
  before the 6:08 PM EDT catalogue), hierarchical trees only, a row's
  members deleted when its sector builds it, pair rates from the middle
  of the evidence.
  Owner: Core 3 (Boss, 2026-10-16, 2:42 PM EDT: Core 3 takes the 8.1 new
  generate scatters line after GEN.280; the multi-star train follows its
  spine, 3:00 PM EDT). This replaces any earlier Core or Core 1 owner
  line on this item.
  Owner: Core 3, in its multi-star and cluster line (the 3:15 PM EDT
  lane queue); a line here that named Core 4 was the TODO thread's
  counting error, corrected 2026-10-16, 4:50 PM EDT.
  Prerequisites: none. Related: GEN.130.
  Design: [docs/design/multi-star-systems.md](design/multi-star-systems.md)

- [ ] **GEN.264 The fill builds a system from its stored members; the system pass keeps only planets, moons and belts**
  Step 5 (design note 8 and 11 item 5). `StarSystem(member_tree=...)`
  skips its own companion draw; stage A and B build from the stored
  members and nodes (copied into `orbit_nodes`), stage C draws planets,
  moons, belts and disks only (no planets under 10 Myr,
  millisecond-pulsar planets at 0.7%, Holman-Wiegert zones from the
  stored nodes); stages D and G call `draw_member_tree` for every object
  they create; the build deletes the child rows. A failed build drops
  the outermost body, then the outermost member, never the placed
  object. Rescopes GEN.251. Done: the stars of a built system equal the
  stored members, the nodes equal the stored nodes, and the child rows
  are gone after the build. Owner: Core (2026-10-15). Boss (2026-10-15
  5:51 PM EDT): stellar companions of stars and star-like objects,
  including neutron star and black hole pairs, are generated at scatter
  time because they matter for gravity; planets, moons and belts stay at
  sector fill. Defaults taken (reversible): compact pairs promoted to
  the scatter, black hole companion fraction 1.25%, cap 12 stars (7
  before the 6:08 PM EDT catalogue), hierarchical trees only, a row's
  members deleted when its sector builds it, pair rates from the middle
  of the evidence.
  Prerequisites: GEN.315, GEN.316. Related: GEN.251, GEN.130. Built-on
  note (Core 1, 2026-10-16, 8:41 AM EDT, GEN.251 merged as #1815): the
  companion draw at fill for stages D and G is not built by GEN.251;
  companion odds stay in multistar.MEMBER_STATE_BY_PLACED (GEN.262),
  which holds planets and moons only.
  Split (2026-10-16, Boss 1:07 PM EDT: look at each item for smaller
  chunks and shared work): this item is now the umbrella and is done
  when its subitems are: GEN.315, GEN.316. Step 5 hides two deliverables
  and a missing prerequisite.
  Scope after the review (2026-10-16 review): The fill reads and writes
  the orbit tables; it adds no migration of its own.
  Owner: Core 3 (Boss, 2026-10-16, 2:42 PM EDT: Core 3 takes the 8.1 new
  generate scatters line after GEN.280; the multi-star train follows its
  spine, 3:00 PM EDT). This replaces any earlier Core or Core 1 owner
  line on this item.
  Owner: Core 3, in its multi-star and cluster line (the 3:15 PM EDT
  lane queue); a line here that named Core 4 was the TODO thread's
  counting error, corrected 2026-10-16, 4:50 PM EDT.
  Design: [docs/design/multi-star-systems.md](design/multi-star-systems.md)

- [ ] **GEN.267 A MAX_STARS setting and a --stars N directive, with a prevalence knob for five or more stars**
  Step 11 (design note 11 item 11; small). `MAX_STARS` (default 12) and
  a `--stars N` directive so a user may force N members or raise the
  cap, and a prevalence knob for the N of 5 or more columns, as GEN.129
  already has. Done: the setting and the directive exist and a forced N
  is honoured or refused with the reason. Owner: Core (2026-10-15). Boss
  (2026-10-15 5:51 PM EDT): stellar companions of stars and star-like
  objects, including neutron star and black hole pairs, are generated at
  scatter time because they matter for gravity; planets, moons and belts
  stay at sector fill. Defaults taken (reversible): compact pairs
  promoted to the scatter, black hole companion fraction 1.25%, cap 12
  stars (7 before the 6:08 PM EDT catalogue), hierarchical trees only, a
  row's members deleted when its sector builds it, pair rates from the
  middle of the evidence.
  Catalogue (Boss 2026-10-15 6:08 PM EDT;
  docs/design/multi-star-systems.md sections 13 to 16): the default cap
  is 12, not 7 (ten candidate systems with multiplicity up to 12; theory
  allows several hundred; the 1.1e6 to 1.1e2 systems of N = 8 to 12 in
  the default galaxy are in the tail table). A directive may still raise
  it.
  Owner: Core 3 (Boss, 2026-10-16, 2:42 PM EDT: Core 3 takes the 8.1 new
  generate scatters line after GEN.280; the multi-star train follows its
  spine, 3:00 PM EDT). This replaces any earlier Core or Core 1 owner
  line on this item.
  Owner: Core 4, "Core 4 (2026-10-16)" (Boss approved the lane at 4:31
  PM EDT, 2026-10-16; queue in the lane report, queue position 9; start
  after Core 3 merges UX.116 (Generate page)).
  Prerequisite: GEN.261. Related: GEN.129, GEN.180.
  Design: [docs/design/multi-star-systems.md](design/multi-star-systems.md)

- [ ] **GEN.269 Interaction states in the star model: a Roche classifier, `interaction` on every pair, member state and flags, rejuvenated age and the new kinds**
  Step 13 (design note sections 13 and 14, items 1 to 5). In
  `physics/stellar_evolution.py`, `Star` and `star_labels.py`: the
  Eggleton Roche-lobe radius and filling factor at the system's age; a
  classifier that gives each pair an `interaction` of detached,
  semi-detached, contact, accreting, post-common-envelope or
  merger-product (stored as `interaction` on `scatter_nodes` and
  `orbit_nodes`); the member flags `mass_reversed`, `mass_gainer`,
  `merger_product`, `be`, `fragile` and `theoretical` (`runaway_origin`
  is filled in by GEN.274); a rejuvenated effective age for a merger
  product, with the real system age kept; and the new kinds
  `stripped-helium`, `thorne-zytkow` and `wolf-rayet` (GEN.243 supplies
  the W evolution; this item adds the kind). The states come from
  geometry and evolution, and the tuning knobs of GEN.268 align their
  rates with the observed ones. Done: every pair the tree generator
  draws carries an `interaction` consistent with its radii (audit test
  13 of section 15), the new kinds exist in the star model and the
  labels, and a merger product has a rejuvenated age. Owner: Core
  (2026-10-15). Boss (2026-10-15 6:08 PM EDT) widened the multi-star
  plan to "every observed and theoretical type of
  multi-star/star-like-object system", built "into our generator
  throughout to match as close as possible to observed and theoretical
  probabilities (observed are best, theoretical are good too)".
  Findings (Core 2, 2026-10-15, #1705): (1) outer companions with P over
  1e5 d should draw their mass from the IMF independently (design note
  3.3); the member tree uses m1*U(0.1,1) for all, so this item or a
  follow-up changes it. (2) Defaults taken in GEN.261, one constant each
  and reversible: massive-primary period law mean 3.5, sigma 1.5;
  neutron-star progenitors 8 to 30 solar masses; direct collapse when
  the progenitor is 40 or more (unknown progenitor: hole at 20 or more);
  black-hole tertiary 30%, compact in 10% of those; the N=8 to 12 tail
  comes out of the 7-star column. (3) test_system_type_rates requires a
  generator only once its owner is in _LANDED_OWNERS (261 now); each
  later owner adds its number.
  Progress (Core 2, PR #1745): the generator side merged on 2026-10-15:
  physics/roche.py, the interaction classifier on every pair, the member
  flags, a rejuvenated merger-product age and the new kinds
  stripped-helium, thorne-zytkow and wolf-rayet with labels. Still open:
  storing interaction on scatter_nodes and orbit_nodes, which waits on
  GEN.260 (migration 0093); exotic-star draws and rate alignment stay
  with GEN.270.
  Prerequisites: GEN.260, DB.53. Related: GEN.268, GEN.243, GEN.242,
  GEN.270, GEN.271, GEN.272. Finding from Core 2 (2026-10-16, 4:58 AM
  EDT, GEN.158 merged as #1789): a binary's secondary shares its
  primary's [Fe/H] (stored column Star.metallicity_fe_h, migration
  0095), but the members of a multi-star group built by the member-tree
  generator do not: they must share it too, so the storage part of this
  item gives every member of a group the group's [Fe/H].
  Scope after the review (2026-10-16 review): Its storage half writes
  the interaction and flag columns that the first child of GEN.129
  provides; no revision of its own.
  Owner: Core 3 (Boss, 2026-10-16, 2:42 PM EDT: Core 3 takes the 8.1 new
  generate scatters line after GEN.280; the multi-star train follows its
  spine, 3:00 PM EDT). This replaces any earlier Core or Core 1 owner
  line on this item.
  Owner: Core 3, in its multi-star and cluster line (the 3:15 PM EDT
  lane queue); a line here that named Core 4 was the TODO thread's
  counting error, corrected 2026-10-16, 4:50 PM EDT.
  Design: [docs/design/multi-star-systems.md](design/multi-star-systems.md)

- [ ] **GEN.270 Close-binary types: contact, Algol, spectroscopic and eclipsing flags, tertiary shares, merger products and Be stars**
  Step 14 (design note sections 13.1 and 13.2). On GEN.269's classifier:
  contact binaries (0.2% of FGK dwarfs, a 60% tertiary share, 1% of O
  and B binaries over-contact), semi-detached Algols (0.1% to 0.2% of
  disk stars, `mass_reversed` donors), the spectroscopic and eclipsing
  flags derived from the geometry (a pair eclipses when the inclination
  is within asin((R1+R2)/a) of edge-on), field blue stragglers and
  merger products (0.5% of old and intermediate stars), massive-star
  merger products (8%) and mass gainers (20%), hot subdwarfs and
  stripped helium stars, Be stars (about 20% of B stars, 30% with a
  stripped companion), Wolf-Rayet binaries (WN 52%, WC 90%), the R
  Coronae Borealis and FK Comae flags, and the Thorne-Zytkow default
  (100 in the Galaxy, flagged theoretical). Contact and Algol rates are
  tuned with the tight-tail weight. Done: over a seeded fill each type's
  count is within four Poisson sigma of its default (audit test 10) and
  no type appears where its gate forbids it (audit test 11). Owner:
  Interface (2026-10-15; moved off Core by the lane queue, Core owns the
  star model or the key it builds on). Boss (2026-10-15 6:08 PM EDT)
  widened the multi-star plan to "every observed and theoretical type of
  multi-star/star-like-object system", built "into our generator
  throughout to match as close as possible to observed and theoretical
  probabilities (observed are best, theoretical are good too)".
  Scope after the review (2026-10-16 review): contact, Algol,
  spectroscopic and eclipsing flags, tertiary shares, merger products
  and Be stars, without the Wolf-Rayet binaries (GEN.243 owns the
  Wolf-Rayet companion fractions). GEN.249 narrows to carbon stars and
  Mira variables, reusing the `flags` column; the Be part is this item.
  Owner: Core 3 (Boss, 2026-10-16, 2:42 PM EDT: Core 3 takes the 8.1 new
  generate scatters line after GEN.280; the multi-star train follows its
  spine, 3:00 PM EDT). This replaces any earlier Core or Core 1 owner
  line on this item.
  Owner: Core 4 (handed over from Core 3 by Boss's "Split it", 2026-10-16,
  4:49 PM EDT; Core 3 agreed 4:53 PM EDT; waits on GEN.269 (Core 3)). Core 4 asks Core 3
  before touching stellar evolution constants in tuning.py. A line here
  that named Core 3 or Core 4 earlier is replaced by this one.
  Prerequisite: GEN.269. Related: GEN.243, GEN.246, GEN.162, GEN.268.
  Design: [docs/design/multi-star-systems.md](design/multi-star-systems.md)

- [ ] **GEN.271 White-dwarf binary types: close double white dwarfs, post-common-envelope pairs, cataclysmic variables, AM CVn, symbiotic and barium stars**
  Step 15 (design note section 13.3); extends GEN.253. A white dwarf
  with a companion 30% (22% of them wide); a close double white dwarf
  9.5% of white dwarfs (axis log-uniform 0.001 to 4 AU, the Peters
  rule); double white dwarf merger products (10% of white dwarfs, mass
  0.8 to 1.1 solar masses, flagged); post-common-envelope pairs (axis
  0.005 to 2 AU); cataclysmic variables (6.5e-4 of white dwarfs, periods
  1.3 to 10 h, 36% magnetic); AM CVn binaries (7.5e-5 of white dwarfs,
  periods 5 to 65 min, a helium donor of kind `stripped-helium`); barium
  and CH stars (1% of G and K giants) and symbiotic stars (5e-6 of
  giants); white dwarf plus brown dwarf. A Type Ia explosion is an event
  and out of scope. Done: counts within four Poisson sigma of the
  defaults (close double white dwarf 9.5% of white dwarfs, cataclysmic
  variable 3.4e-5 of systems, AM CVn 3.9e-6), and every accreting pair
  has `interaction` accreting. Owner: Interface (2026-10-15; moved off
  Core by the lane queue, Core owns the star model or the key it builds
  on). Boss (2026-10-15 6:08 PM EDT) widened the multi-star plan to
  "every observed and theoretical type of multi-star/star-like-object
  system", built "into our generator throughout to match as close as
  possible to observed and theoretical probabilities (observed are best,
  theoretical are good too)".
  Related to GEN.253, not dependent on it: the table here is
  self-contained and GEN.253's pollution flag and cool white-dwarf
  presets stay independent (2026-10-16 review).
  Owner: Core 3 (Boss, 2026-10-16, 2:42 PM EDT: Core 3 takes the 8.1 new
  generate scatters line after GEN.280; the multi-star train follows its
  spine, 3:00 PM EDT). This replaces any earlier Core or Core 1 owner
  line on this item.
  Owner: Core 4 (handed over from Core 3 by Boss's "Split it", 2026-10-16,
  4:49 PM EDT; Core 3 agreed 4:53 PM EDT; waits on GEN.269 (Core 3) and GEN.253). Core 4 asks Core 3
  before touching stellar evolution constants in tuning.py. A line here
  that named Core 3 or Core 4 earlier is replaced by this one.
  Prerequisites: GEN.269, GEN.253. Related: GEN.253, GEN.130, GEN.268.
  Design: [docs/design/multi-star-systems.md](design/multi-star-systems.md)

- [ ] **GEN.272 Accreting compact binaries: X-ray binaries, accretion luminosity, microquasars, spider pulsars and compact triples**
  Step 16 (design note section 13.4); extends slice (c) of GEN.130.
  High-mass X-ray binaries (the 0.3% compact companion of a living
  massive star, 40% of them accreting, Be/X-ray 60% of those); low-mass
  X-ray binaries (5e-5 of neutron stars and 1e-3 of black holes, 1%
  accreting); the microquasar flag with jets along the spin axis
  (GEN.255); millisecond-pulsar companion classes (white dwarf 80%,
  black widow 8%, redback 4%, neutron star 5%, other 3%) with the
  ablation flag; pulsars with main-sequence or Be companions; accretion
  luminosity L = eta Mdot c^2 capped at the Eddington limit; dormant
  black hole and neutron star companions (the Gaia kind); compact
  triples (30% of compact binaries, at most three compact members);
  magnetar binaries at 0 by default. Done: an accreting pair shows an
  X-ray luminosity at or under Eddington, counts are within four Poisson
  sigma of the defaults, and no Be/X-ray binary appears in the halo.
  Owner: Core (2026-10-15). Boss (2026-10-15 6:08 PM EDT) widened the
  multi-star plan to "every observed and theoretical type of
  multi-star/star-like-object system", built "into our generator
  throughout to match as close as possible to observed and theoretical
  probabilities (observed are best, theoretical are good too)".
  Prerequisite: GEN.269. Related: GEN.130, GEN.255, GEN.263,
  GEN.268.
  Owner: Core 3 (Boss, 2026-10-16, 2:42 PM EDT: Core 3 takes the 8.1 new
  generate scatters line after GEN.280; the multi-star train follows its
  spine, 3:00 PM EDT). This replaces any earlier Core or Core 1 owner
  line on this item.
  Owner: Core 4 (handed over from Core 3 by Boss's "Split it", 2026-10-16,
  4:49 PM EDT; Core 3 agreed 4:53 PM EDT; waits on GEN.269 (Core 3)). Core 4 asks Core 3
  before touching stellar evolution constants in tuning.py. A line here
  that named Core 3 or Core 4 earlier is replaced by this one.
  Design: [docs/design/multi-star-systems.md](design/multi-star-systems.md)

- [ ] **GEN.273 Young groups and trapezia: the `group` node, offset and velocity columns, a decay age and the remnant draw**
  Step 17 (design note sections 13.1 and 14 item 6); works with GEN.241
  and GEN.244. A node of kind `group` with N members, each with an
  offset (AU) and a velocity (km/s) at an epoch, and a `decay_age_yr`;
  internal motion is a straight line, there are no Kepler elements, and
  the nightly update does not integrate it. Default 5% of massive (O and
  early B) systems younger than 0.3 Myr, drawn at fill in star-forming
  cores and at scatter for the few mass-pass stars of that age; the
  sector page shows the group with its remaining life. A group older
  than its decay age is never drawn: the generator draws its result, the
  lightest member ejected as a runaway row (`runaway_origin` dynamical)
  and a tighter hierarchical remnant. Done: every stored group is
  younger than its `decay_age_yr` (audit test 14), a group older than
  that is replaced by its remnant, and members carry their offset and
  velocity columns. Owner: Core (2026-10-15). Boss (2026-10-15 6:08 PM
  EDT) widened the multi-star plan to "every observed and theoretical
  type of multi-star/star-like-object system", built "into our generator
  throughout to match as close as possible to observed and theoretical
  probabilities (observed are best, theoretical are good too)".
  Prerequisites: GEN.327, GEN.328. Related: GEN.241, GEN.244, GEN.274,
  UX.114, GEN.268.
  Split (2026-10-16, Boss 1:07 PM EDT: look at each item for smaller
  chunks and shared work): this item is now the umbrella and is done
  when its subitems are: GEN.327, GEN.328. GEN.273 says it 'works with
  GEN.241 and GEN.244' and draws groups 'at fill in star-forming cores',
  but the cores are made by GEN.241 in a parallel thread of the same
  phase, with no prerequisite and no legal way to add one.
  Scope after the review (2026-10-16 review): The group generator and
  the decay-age rule; the offset, velocity and decay columns come from
  the first child of GEN.129.
  Owner: Core 3 (Boss, 2026-10-16, 2:42 PM EDT: Core 3 takes the 8.1 new
  generate scatters line after GEN.280; the multi-star train follows its
  spine, 3:00 PM EDT). This replaces any earlier Core or Core 1 owner
  line on this item.
  Owner: Core 3, in its multi-star and cluster line (the 3:15 PM EDT
  lane queue); a line here that named Core 4 was the TODO thread's
  counting error, corrected 2026-10-16, 4:50 PM EDT.
  Design: [docs/design/multi-star-systems.md](design/multi-star-systems.md)

- [ ] **GEN.274 Cluster-hosted compact binaries and runaway bookkeeping: `runaway_origin` and the compact companion left by a supernova**
  Step 18 (design note sections 13.4 and 13.6). `runaway_origin` (none,
  supernova, dynamical) on a runaway, at the observed rates (10% to 30%
  of O stars and 5% to 10% of B stars; two-thirds supernova, one-third
  dynamical for the massive ones); for the supernova channel a compact
  companion left behind at the rate the survival rule gives; companions
  unbound by a collapse (83% of non-interacting massive binaries) are
  not stored, the survivor is an ordinary row with the runaway flag;
  hypervelocity stars stay single. Cluster-hosted compact binaries
  replace the field fractions: a black hole with a star 10%, a
  millisecond pulsar with a companion 80% with a higher neutron-star
  share, and dynamical black-hole pairs from the cluster density. Done:
  counts within the ranges (audit test 10), no runaway carries a missing
  origin, and a cluster sector draws its compact binaries at the cluster
  rates, not the field rates. Owner: Interface (2026-10-15; moved off
  Core by the lane queue, Core owns the star model or the key it builds
  on). Boss (2026-10-15 6:08 PM EDT) widened the multi-star plan to
  "every observed and theoretical type of multi-star/star-like-object
  system", built "into our generator throughout to match as close as
  possible to observed and theoretical probabilities (observed are best,
  theoretical are good too)".
  Prerequisites: GEN.329, GEN.330. Related: GEN.159, GEN.160, GEN.161,
  GEN.162, GEN.163, GEN.164, GEN.273.
  Split (2026-10-16, Boss 1:07 PM EDT: look at each item for smaller
  chunks and shared work): this item is now the umbrella and is done
  when its subitems are: GEN.329, GEN.330. The first half
  (runaway_origin on runaways at the observed O and B rates, the compact
  companion a supernova leaves, unbound companions not stored) is field
  bookkeeping needing only the member tree.
  Scope after the review (2026-10-16 review): The runaway bookkeeping;
  the column comes from the first child of GEN.129.
  Owner: Core 3 (Boss, 2026-10-16, 2:42 PM EDT: Core 3 takes the 8.1 new
  generate scatters line after GEN.280; the multi-star train follows its
  spine, 3:00 PM EDT). This replaces any earlier Core or Core 1 owner
  line on this item.
  Owner: Core 4 (handed over from Core 3 by Boss's "Split it", 2026-10-16,
  4:49 PM EDT; Core 3 agreed 4:53 PM EDT; waits on GEN.329 and GEN.330 (Core 3)). Core 4 asks Core 3
  before touching stellar evolution constants in tuning.py. A line here
  that named Core 3 or Core 4 earlier is replaced by this one.
  Design: [docs/design/multi-star-systems.md](design/multi-star-systems.md)

- [ ] **GEN.275 Substellar and circumbinary types: brown-dwarf and planetary-mass companions, circumbinary disks and planets**
  Step 19 (design note section 13.5). Brown dwarf close companions
  (under 0.6% inside 3 to 10 AU, the desert), wide companions (3.2% at
  28 to 1,590 AU), brown dwarf pairs (15% of brown dwarfs),
  planetary-mass companions (under 13 Jupiter masses, flagged
  `planetary-mass`), circumbinary planets (at least 10% of close
  eclipsing binaries, the Holman-Wiegert P-type zone from the node,
  planets at fill only), circumbinary disks (a descriptor on the node,
  no bodies), S-type planets in wide pairs and rogue planet pairs.
  GEN.94 (the feasibility study) is done; GEN.95 applies the nebula
  conditions. Done: counts within four Poisson sigma of the defaults
  (brown dwarf wide companion 3.2%), a circumbinary planet's orbit lies
  outside the stable zone, and a disk carries no bodies. Owner:
  Interface (2026-10-15; moved off Core by the lane queue, Core owns the
  star model or the key it builds on). Boss (2026-10-15 6:08 PM EDT)
  widened the multi-star plan to "every observed and theoretical type of
  multi-star/star-like-object system", built "into our generator
  throughout to match as close as possible to observed and theoretical
  probabilities (observed are best, theoretical are good too)".
  Split (2026-10-16, Boss 1:07 PM EDT: look at each item for smaller
  chunks and shared work): this item is now the umbrella and is done
  when its subitems are: GEN.331, GEN.332. Two different kinds of work
  under one title: substellar companions are scatter-time member rows
  drawn by the tree generator (brown dwarfs, planetary-mass companions),
  while circumbinary planets, S-type planets, circumbinary disks and
  rogue planet pairs are fill-time planet and disk work that needs the
  zones of a tree and the disk descriptor.
  Scope after the review (2026-10-16 review): The circumbinary disk law
  fills the same descriptor on a node.
  Owner: Core 3 (Boss, 2026-10-16, 2:42 PM EDT: Core 3 takes the 8.1 new
  generate scatters line after GEN.280; the multi-star train follows its
  spine, 3:00 PM EDT). This replaces any earlier Core or Core 1 owner
  line on this item.
  Owner: Core 4 (handed over from Core 3 by Boss's "Split it", 2026-10-16,
  4:49 PM EDT; Core 3 agreed 4:53 PM EDT; waits on GEN.331 and GEN.332 (Core 3)). Core 4 asks Core 3
  before touching stellar evolution constants in tuning.py. A line here
  that named Core 3 or Core 4 earlier is replaced by this one.
  Prerequisites: GEN.331, GEN.332. Related: GEN.252, GEN.95, GEN.242.
  Design: [docs/design/multi-star-systems.md](design/multi-star-systems.md)

- [ ] **GEN.276 The sphere-to-key-runs cover: `sphere_runs` over the existing cell key (at most 256 merged runs)**
  Step 2 of the region-query train (design note sections 6.1 and 6.3,
  item 2). `cell_key(points)` weaves the three cube numbers `floor((p +
  reach) / cell)` (21 bits each, x in the top bit of each triple) into
  one 63-bit unsigned integer, with the cell 0.25 pc and the reach
  262,144 pc by default; `sphere_runs(center, r, budget=256)` lays the
  sphere over the nested grid and returns at most 256 ascending, merged
  runs of keys, about 1.5 ms in numpy. Integers and numpy only: it must
  not import numba, because the web app never does. When the plan stores
  another reach the cell size follows, `2 x reach / 2^21`. The 21-bit
  constants (the five mask-and-shift steps) and the outside-the-box case
  are not yet tested: test them here. Done: a fuzz against brute force
  shows every point inside a random sphere lands in a run, runs never
  exceed the budget, points just inside and just outside the box edge
  are handled as the design says (outside gives no key), and the
  function is the one the SQL column of DB.34 is tested against. Owner:
  Core (2026-10-15). Boss (2026-10-15 5:32 PM EDT): "if it's not ready
  in time before the spherical sector layout is, then we'll skip it and
  move to the next sequence until it's finished research and we can plan
  it." The research is finished and merged (#1659), so stage C is live.
  Prerequisites: none. Related: DB.34, NAV.60, TEST.165, DOC.51,
  GEN.233.
  The cell key (`planetgen/galaxy/cellkey.py`: `cell_key`, `cell_keys`,
  `sql_cell_key`, tested) already exists; only the sphere-to-key-runs
  cover is built here, numpy only, fuzzed against brute force, with the
  outside-the-box edge cases (2026-10-16 review).
  Owner: Core 2 (Boss, 2026-10-16, 3:00 PM EDT: the region query goes to
  Core 2), after the infinite zoom tree and its labels.
  Design: [docs/design/spherical-region-queries.md](design/spherical-region-queries.md)

- [ ] **GEN.282 Galactic core gas and young terms: the central molecular zone, the nebula rate and young clusters (Arches, Quintuplet, the central parsec)**
  Build step of the realistic galactic core plan
  (docs/design/realistic-galactic-core.md, section 9; Boss's item,
  2026-10-16, 12:43 AM EDT; merged as #1728). Owner Core 3, goal issue
  #1677. Done: the central molecular zone in the gas model, the nebula
  rate, the nuclear disc's young share, and young-cluster components for
  the central parsec, the Arches and the Quintuplet through GEN.244; no
  object is placed by name, the generator makes objects of those kinds
  where the profile puts them.
  Owns the Central Molecular Zone and the nuclear disc's young share
  (2026-10-16 review); GEN.134 does not build them.
  Prerequisites: GEN.244. Related: GEN.280, GEN.246.
  Design: [docs/design/realistic-galactic-core.md](design/realistic-galactic-core.md)

- [ ] **GEN.284 Planet survival at the galactic core: a metallicity factor and a survival formula instead of the hard cut**
  Build step of the realistic galactic core plan
  (docs/design/realistic-galactic-core.md, section 9; Boss's item,
  2026-10-16, 12:43 AM EDT; merged as #1728). Owner Core 3, goal issue
  #1677. Done: planets in dense components use the survival formula and
  a metallicity factor, not the hard cut at the Hill radius; the
  survival table of the note is a test.
  Scope after the review (2026-10-16 review): feed the nucleus' density
  and dispersion (from GEN.280) into the one survival function that the
  first child of GEN.162 builds, and add the survival table as a test;
  no second cull.
  Owner: Core 4, "Core 4 (2026-10-16)" (Boss approved the lane at 4:31
  PM EDT, 2026-10-16; queue in the lane report); held until its blocker
  merges: GEN.325 (Core 3). Core 4 takes it when it unblocks.
  Prerequisites: GEN.325. Related: GEN.158.
  Design: [docs/design/realistic-galactic-core.md](design/realistic-galactic-core.md)

- [ ] **GEN.285 Nuclear mass components in the galactic potential for the orbital update**
  Build step of the realistic galactic core plan
  (docs/design/realistic-galactic-core.md, section 9; Boss's item,
  2026-10-16, 12:43 AM EDT; merged as #1728). Owner Core 3, goal issue
  #1677. Done: the nuclear cluster and the nuclear disc enter
  galactic-potential.py so orbits near the centre feel them.
  Finding (Core 3, 2026-10-16, 8:42 PM EDT): `galactic-potential.py`
  does not exist; main has only the old rotation curve in
  `galaxy/galactic_orbit.py`. The Hernquist, Miyamoto-Nagai and NFW
  module is GEN.338 (needs GEN.144), unbuilt and unowned, in 8.3. So
  this item needs GEN.338 and moves from 8.1 to 8.3 behind it; it no
  longer gates the 8.1 headline (GEN.301). Owner: Core 3 proposed for
  GEN.144, GEN.338 and this item when 8.3 starts (nothing before 8.1
  closes).
  Prerequisites: GEN.338. Related: GEN.280 (built), GEN.198.
  Design: [docs/design/realistic-galactic-core.md](design/realistic-galactic-core.md)

- [ ] **GEN.286 Horizontal-branch and AGB phases in the stellar evolution model**
  Boss (2026-10-16, 1:12 AM EDT): "why not build the model of
  horizontal-branch stars into the mix?" The stellar evolution model
  (src/planetgen/physics/stellar_evolution.py) has no horizontal-branch
  phase, which makes the halo scatter counts low. Plan: section 5.3 of
  docs/design/realistic-galactic-core.md (merged, #1735). Done:
  evolve_star gains core helium-burning phases (red clump or red
  horizontal branch, RR Lyrae, blue horizontal branch, extreme
  horizontal branch: 10 to 60 solar luminosities, 4,600 to 35,000 K, 0.1
  Gyr) and AGB phases (100 to 2,500 solar luminosities, 0.02 Gyr) for
  stars under 2.2 solar masses between the giant phase and the white
  dwarf; the colour is drawn from [Fe/H] and age; star_population's
  _bright_windows, _bright_table and mass-band tables mirror the new
  phases; Yerkes III with the type text, no column, no migration.
  Measured: halo scatter rows at 100, 500 and 1000 floors are 12
  million, 3.8 million and 310,000 now; the horizontal branch adds none,
  the AGB adds 2.6%, 4% and 28%. [Fe/H] comes from GEN.158 (Core 2) or,
  until it lands, from interim population defaults (halo -1.6 and -2.2,
  disc and bulge -0.1). Boss (2026-10-15 10:06 PM EDT): the
  horizontal-branch stars are to be fully part of the generation of a
  sector, the fill as well as the scatter, in the proper proportions,
  with the giant-branch spread of GEN.287 built in the same pipeline;
  "go ahead and build the fix into
  your pipe". Acceptance test (Core 3's figures): sampled old, halo and
  bulge populations put about 0.8 to 1.2% of their stars on the
  horizontal branch, the colour shares follow [Fe/H] (-2.2, -1.5, -0.1),
  halo, bulge and globular-cluster cells all show them, and scatter plus
  fill still sum to the untiered expected total. Owner Core 3 (Boss,
  1:12 AM EDT). Goal issue #1677.
  Prerequisites: none. Related: GEN.158, GEN.283, GEN.234, GEN.237,
  GEN.159, GEN.160, GEN.236.
  Design: [docs/design/realistic-galactic-core.md](design/realistic-galactic-core.md)

- [ ] **GEN.287 The stellar evolution model spreads giant-branch stars evenly in log luminosity, which probably overstates the bright giants by about a factor of two against real halo counts (bug)**
  Finding from Core 3 (2026-10-16, 1:14 AM EDT, kind: model accuracy,
  written in its galactic core note
  docs/design/realistic-galactic-core.md, no fix planned there): the
  giant-branch stars of the stellar evolution model are spread evenly in
  log luminosity, which probably overstates the bright giants by about a
  factor of two compared with real halo counts. Boss wants densities and
  star counts to match real observations. Done: the giant-branch
  luminosity distribution follows the time spent at each luminosity (the
  real luminosity function) and the halo giant counts are checked
  against published halo counts; built together with GEN.286 in the same
  pipeline, not queued separately (Boss, 2026-10-15 10:06 PM EDT: "go
  ahead and build the fix into your pipe"). Owner Core 3, goal issue
  #1677.
  Prerequisites: none. Related: GEN.286, GEN.283.

- [ ] **GEN.289 Planets and moons kept as the universal point-in-space object, with point mass and movement vector stored, so perturbations have somewhere to land**
  Finding (orbital-encounters research, 2026-10-15, kind: gap in
  GEN.109): today a planet or moon is a distance, a period and a phase
  on a circle (physics and db/store.py advance_orbital_phases:
  orbital_phase_deg = MOD(phase + turn, 360), position and velocity
  recomputed from the phase), so no pull from another body can change
  its orbit. GEN.109 does not say that planets and moons must first
  become state vectors (position and velocity relative to the host star,
  at an epoch) advanced by a real propagator (kepler.universal_step for
  the quiet case, a sub-stepped integrator while a perturber is inside
  the system) before its influence set can apply. Done: planets and
  moons are stored as vectors with parent and dynamical state (one
  column per component, no JSON), the phase columns are removed (Core 1
  adds the drop revision after the readers move), the pages, the API and
  the JavaScript position parity read the derived elements, the quiet
  case advances by `kepler.universal_step`, and a test shows a quiet
  system costs about what `advance_orbital_phases` costs now. Source:
  Boss's question of
  2026-10-15 at 10:53 PM EDT (if the orbital motion ran 100 years of
  movement in 1 second, what would the program account for when a black
  hole or neutron star got close to a planetary system, and would
  planets be ejected?) and its answer,
  /mnt/project-files/research/orbital-encounters/answer.md. Boss asked
  on 2026-10-15 (11:06 PM EDT) for these gaps to be filled; they join
  the Phase 2 and 3 workflow and no lane builds them yet. Placement:
  version 8.3, phase 1 (the Positions thread), before GEN.109. Boss's
  design
  (2026-10-15, 11:06 PM EDT): at each step in which orbital changes are
  calculated, a planet or moon is moved by the summed gravitational
  influence of point masses (the vector of the gravitational gradient of
  the space it occupies), the other planets and moons of its system
  counted by their point masses too, plus passing objects such as black
  holes and neutron stars, so planets can be captured or ejected; the
  orbital path comes out right and agrees with the
  pass-through-the-sphere logic (the sector-crossing guard of GEN.223
  and GEN.265). The planets and moons of a system enter the influence
  set as point masses next to the passing bodies, so GEN.289's vectors
  are what GEN.109 moves.
  Boss (2026-10-15, 11:10 PM EDT): "Use our universal point in space
  object for all objects and it should have a space for the point mass
  and the vector of current movement." So every planet and moon is a
  SpatialPosition3D (physics/position.py, GEN.74), the same object every
  other body uses: its position, its velocity vector (system frame,
  relative to the parent) and its point mass (mass_kg, mu following)
  travel with it, and the gravity sum, the encounter search and the
  integrator read and write that one object only. Stored one column per
  component (x, y, z, vx, vy, vz, mass, epoch, parent type and id,
  dynamical state), no JSON, in the point-column layout GEN.223 built; the
  orbital elements shown on the pages are derived from the vector, and
  the circular phase columns are replaced, not kept beside them (no
  compatibility layer). GEN.294 gives every other object kind the same
  object; the one migration for both is DB.55, whose revision number
  Core 1 names. Design:
  /mnt/project-files/research/orbital-encounters/capture-and-ejection-design.md.
  Prerequisites: DB.55. Related: GEN.105, GEN.109, GEN.110,
  GEN.223,
  GEN.294.
  Rewritten (2026-10-16 review): this item is the build, not a design
  note. The shared migration moved to the migration child of GEN.294
  (one Core 1 revision); the prerequisites are the first and the
  migration children of GEN.294, not GEN.115.
  Read path (added (2026-10-16 review)):
  `physics/body_positions.positions_at` and its browser twin
  `static/orbitpositions.js` read the stored state vectors and advance
  them by the universal-step solver for the quiet case, with the Python
  and JavaScript parity test kept; NAV.27, ADM.72 and MAP.169 need it.

- [ ] **GEN.290 An encounter search that cannot miss a pass at large steps (100 years per second)**
  Finding (orbital-encounters research, 2026-10-15, kind: gap in
  GEN.109): GEN.109 describes forces from nearby bodies but no search
  for close passes that stays safe when one step covers many years. A
  pass at speed v and distance b needs a slice well under b / v (a
  star-on-star pass at 40 km/s and 1,000 AU takes about 120 years), so a
  century in one step steps straight over it. Done: a design for
  stepping each object in slices short enough that a pass is never
  missed, or an analytic closest-approach test between the paths of
  nearby bodies over the whole step (a body with a mass at or above the
  host's, such as a black hole or a neutron star, first), so that a
  system is switched to the sub-stepped integrator only while a
  perturber is inside its Hill sphere or its planets' orbits and every
  other system keeps the cheap clockwork; the cost estimate on the
  benchmark galaxy; and a test with a seeded black hole passing a seeded
  system at a large step. Source: Boss's question of 2026-10-15 at 10:53
  PM EDT (if the orbital motion ran 100 years of movement in 1 second,
  what would the program account for when a black hole or neutron star
  got close to a planetary system, and would planets be ejected?) and
  its answer, /mnt/project-files/research/orbital-encounters/answer.md.
  Boss asked on 2026-10-15 (11:06 PM EDT) for these gaps to be filled;
  they join the Phase 2 and 3 workflow and no lane builds them yet.
  Placement: Phase 2, after GEN.109 and GEN.289. Passing bodies (black
  holes, neutron stars, any body at or above the host's mass) are
  searched first; GEN.109's Hill-sphere influence radius decides who
  joins a system's set. Boss's design (2026-10-15, 11:06 PM EDT): at
  each step in which orbital changes are calculated, a planet or moon is
  moved by the summed gravitational influence of point masses (the
  vector of the gravitational gradient of the space it occupies), the
  other planets and moons of its system counted by their point masses
  too, plus passing objects such as black holes and neutron stars, so
  planets can be captured or ejected; the orbital path comes out right
  and agrees with the pass-through-the-sphere logic (the sector-crossing
  guard of GEN.223 and GEN.265).
  Prerequisites: GEN.289, GEN.341, GEN.343, GEN.334. Related: GEN.105,
  GEN.109,
  GEN.110,
  GEN.115,
  GEN.223.
  Scope after the review (2026-10-16 review): deciding per system when
  to switch to the perturbed pass (passing bodies of host mass or more
  first, the system-radius filter from the first child of GEN.109), the
  cost estimate on the benchmark galaxy and the seeded black-hole test;
  the closest-approach maths comes from the library.

- [ ] **GEN.291 A bound planet becomes a rogue planet when an encounter unbinds it: the unbind and capture test and the event log**
  Finding (orbital-encounters research, 2026-10-15, kind: missing
  behaviour): nothing in the update decides whether a planet stays
  bound; rogue planets come only from generation, never from the update.
  The design table in docs/design/orbital-updates.md (an object flung
  out is kept and flagged, never deleted) is marked to come. Done: after
  every encounter the planets' specific energy and eccentricity relative
  to their star are recomputed (state_vectors.elements_from_state
  handles hyperbolas); a planet at or above zero energy, or past the
  system's Hill radius, is flagged unbound and becomes a rogue planet
  row (kept, never deleted), a captured body gets an orbit around its
  new host, and every ejection, capture, collision and warning goes to
  the admin activity log with the event time. Collisions use GEN.110.
  Tests: a black hole passing a seeded system ejects the outer planets
  first, and the totals before and after match the rule. Source: Boss's
  question of 2026-10-15 at 10:53 PM EDT (if the orbital motion ran 100
  years of movement in 1 second, what would the program account for when
  a black hole or neutron star got close to a planetary system, and
  would planets be ejected?) and its answer,
  /mnt/project-files/research/orbital-encounters/answer.md. Boss asked
  on 2026-10-15 (11:06 PM EDT) for these gaps to be filled; they join
  the Phase 2 and 3 workflow and no lane builds them yet. Placement:
  Phase 2, last step of the orbital thread, after GEN.289, GEN.290 and
  GEN.110, and before GEN.105 is done (its parent). Boss's design
  (2026-10-15, 11:06 PM EDT): at each step in which orbital changes are
  calculated, a planet or moon is moved by the summed gravitational
  influence of point masses (the vector of the gravitational gradient of
  the space it occupies), the other planets and moons of its system
  counted by their point masses too, plus passing objects such as black
  holes and neutron stars, so planets can be captured or ejected; the
  orbital path comes out right and agrees with the
  pass-through-the-sphere logic (the sector-crossing guard of GEN.223
  and GEN.265).
  Prerequisites: GEN.289, GEN.290, GEN.343, GEN.346, GEN.197. Related:
  GEN.105,
  GEN.109, GEN.110,
  GEN.115, GEN.223, GEN.289, GEN.290.
  The unbind and capture test and the rogue-planet conversion; its
  "event log" half is calls to `record_event` (ejection, capture)
  (2026-10-16 review).

- [ ] **GEN.292 Run the orbital update for a stated span in one go from the command line (elapsed years)**
  Finding (orbital-encounters research, 2026-10-15, kind: gap in
  GEN.105): GEN.105 records Boss's wish (2026-10-07) for an option that
  advances the galaxy by a stated extra span in one go, with the default
  one day of real time for one day of orbit. orbit_clock(conn,
  elapsed_years=...) in db/store.py already steps by a chosen number of
  years, but only tests call it, so the command line (planetgen.cli)
  still only runs to wall-clock time. Done: a documented option on the
  orbital update command (for example --advance-years or --advance-days)
  that runs one update over that span, shows the span in the progress
  and the activity log, refuses a span over the one billion year
  interval cap unless told otherwise, and says in the output which
  behaviours do not hold at that span until GEN.290 lands (clockwork
  phases stay correct at any span; nothing perturbs). Test: a
  one-century step equals the same century run as many short steps for
  the clockwork objects. Source: Boss's question of 2026-10-15 at 10:53
  PM EDT (if the orbital motion ran 100 years of movement in 1 second,
  what would the program account for when a black hole or neutron star
  got close to a planetary system, and would planets be ejected?) and
  its answer, /mnt/project-files/research/orbital-encounters/answer.md.
  Boss asked on 2026-10-15 (11:06 PM EDT) for these gaps to be filled;
  they join the Phase 2 and 3 workflow and no lane builds them yet.
  Placement: Phase 3, with the gravity map thread, once the integrator
  exists (the option is how the behaviour of GEN.289 to GEN.291 is
  exercised at 100 years per second). Boss's design (2026-10-15, 11:06
  PM EDT): at each step in which orbital changes are calculated, a
  planet or moon is moved by the summed gravitational influence of point
  masses (the vector of the gravitational gradient of the space it
  occupies), the other planets and moons of its system counted by their
  point masses too, plus passing objects such as black holes and neutron
  stars, so planets can be captured or ejected; the orbital path comes
  out right and agrees with the pass-through-the-sphere logic (the
  sector-crossing guard of GEN.223 and GEN.265).
  Prerequisite: GEN.291. Related: GEN.105, GEN.223, GEN.290.

- [ ] **GEN.293 Long-step rules for hypervelocity stars that leave the sector grid and for parabolic comets**
  Finding (orbital-encounters research, 2026-10-15, kind: undefined
  behaviour at long steps): hypervelocity stars fly in a straight line
  with no galactic pull (db/store.py, galaxy/straight_line.py), so a
  long enough step puts one outside every sector; what index.sector_at
  does for a position that has left the galaxy was not traced. A
  parabolic comet's anomaly grows without wrapping
  (advance_comet_orbits), so a long step puts it very far from its star
  and nothing ends the pass or marks the comet gone. Done: the rule for
  each case written and tested: a star past the outer limit of the
  sector grid is kept and flagged (never deleted) and found by the
  intergalactic frame (GEN.223, GEN.238) instead of failing a sector
  lookup, and a parabolic comet beyond a stated distance from its host
  is marked as having left the system and stops being advanced. Tests at
  a step of 100 years and of 1,000,000 years. Source: Boss's question of
  2026-10-15 at 10:53 PM EDT (if the orbital motion ran 100 years of
  movement in 1 second, what would the program account for when a black
  hole or neutron star got close to a planetary system, and would
  planets be ejected?) and its answer,
  /mnt/project-files/research/orbital-encounters/answer.md. Boss asked
  on 2026-10-15 (11:06 PM EDT) for these gaps to be filled; they join
  the Phase 2 and 3 workflow and no lane builds them yet. Placement:
  Phase 3, with GEN.292, after GEN.290 (the long-step rules sit on the
  encounter search). Boss's design (2026-10-15, 11:06 PM EDT): at each
  step in which orbital changes are calculated, a planet or moon is
  moved by the summed gravitational influence of point masses (the
  vector of the gravitational gradient of the space it occupies), the
  other planets and moons of its system counted by their point masses
  too, plus passing objects such as black holes and neutron stars, so
  planets can be captured or ejected; the orbital path comes out right
  and agrees with the pass-through-the-sphere logic (the sector-crossing
  guard of GEN.223 and GEN.265).
  Prerequisite: GEN.290. Related: GEN.105, GEN.223, GEN.232, GEN.238.

- [ ] **GEN.294 Every object in the galaxy is held as one SpatialPosition3D with its position, point mass and movement vector stored as columns**
  Boss asked (2026-10-15, 11:12 PM EDT, project chat): "do we have a
  TODO currently that says every object should have a SpacialPosition3D
  property? Everything that exists in the galaxy." Answer: no open item
  says it. GEN.74 (done, PR #665, "used by every object") built the
  class (physics/position.py) and gave it to sector system and
  phenomenon entries, stars, planets, moons and comets in memory;
  GEN.223 (done, #1757) gave bright_stars and phenomenon_scatter their
  point columns. The stored objects are not uniform: planets, moons,
  comets, facilities and star systems keep position, a velocity vector
  and an epoch; black holes, neutron stars, rogue planets, nebulae,
  supernova remnants, interstellar comets and asteroid fields keep a
  position, an epoch and only a scalar galactic_orbital_speed_kms (no
  velocity vector); the quasar and the asteroid belts keep no epoch or
  vector; stars take their place from the system; point mass is stored
  for stars, planets, moons, black holes, neutron stars and rogue
  planets but not for comets, interstellar comets, asteroid fields,
  belts, nebulae, remnants, quasars or facilities. Done: every object
  that exists in the galaxy (stars and their companions, planets, moons,
  comets, interstellar comets, asteroid belts and fields, rogue planets,
  black holes, neutron stars, quasars, nebulae, supernova remnants,
  facilities, scatter points and scatter members) is held as one
  SpatialPosition3D carrying its position, its point mass and its
  movement vector (position, velocity, mass_kg with mu, epoch_unix),
  stored one column per component with no JSON and a cell key, and read
  back into the same object; a scalar speed column is replaced by the
  vector, never kept beside it (no compatibility layer). A body with no
  real mass (a belt's or field's bulk, a nebula) stores its total mass
  or a stated zero so the gravity sum can include or skip it by rule; a
  facility takes the mass of what it sits on or its own stated mass. A
  test lists every uid table and fails when one lacks the point, the
  mass or the vector, so a new object kind cannot ship without them. The
  planet and moon dynamical state (parent, bound state, the replaced
  circular phase) belongs to GEN.289, which builds on this object; both
  go in the one migration DB.55, revision 0102 (Core 1, 2026-10-16, 2:44 PM EDT: 0099 and 0100 are Core 3's; DB.54 merged first as 0101, 6:22 PM EDT, so DB.55 is 0102), with the later drop revisions following in order; Core 1
  owns every revision. Placed in Phase 2 before GEN.109,
  because the influence set sums the point mass of every object kind. No
  lane has it yet (Boss, 11:07 PM EDT: hand to TODO, nothing built).
  Prerequisites: DB.55, GEN.336, GEN.337, GEN.372. Related: GEN.74,
  GEN.109, GEN.223, GEN.289,
  GEN.290.
  Split (2026-10-16, Boss 1:07 PM EDT: look at each item for smaller
  chunks and shared work): this item is now the umbrella and is done
  when its subitems are: GEN.335, DB.55, GEN.336, GEN.337. It touches
  about fifteen object kinds (eight tables keep only a scalar speed,
  five more lack mass, epoch or vector), with roughly 24 sites in
  db/store.py, 15 in schema.sql and a generator per phenomenon, and it
  hides a Core 1 migration that three other items (GEN.142, GEN.110,
  GEN.109) also wanted columns in.
  Scope after the review (2026-10-16 review): this item is the umbrella
  of one definition of the stored point columns and a test of every
  object table (a), one additive Core 1 revision for the 8.3 schema (m,
  shared with GEN.289, GEN.142, GEN.110 and GEN.109), and the vector and
  mass for compact and free objects (b) and for the rest (c).
  Owner: Core 1 (Boss, 2026-10-16, 2:42 PM EDT: Core 1 takes the global
  top 5 blockers, starting right after GEN.237 merges; estimates from
  Core 1: GEN.335 merged about 7:30 PM EDT, DB.55 about 10:30 PM EDT).
  Waits on GEN.372, the numeric precision audit, before it closes (Boss,
  2026-10-16, 6:54 PM EDT); the audit does not hold up DB.55, GEN.336 or
  GEN.337, which Core 1 starts first.

- [ ] **GEN.295 Nebulae as point groups whose shape follows the orbital update (design and research)**
  Boss asked (2026-10-16, 1:01 AM EDT, project chat): "should nebulae be
  a group of points that then as orbits update the nebula changes
  shape?" Answer given to him: yes in principle, it follows from GEN.294
  (every object one SpatialPosition3D with point mass and movement
  vector). Design and research only, nothing built, no lane. Today a
  nebula is one point (centre, one scalar galactic_orbital_speed_kms,
  epoch_unix, radius_ly as a bounding sphere) plus a shape fixed when it
  was made: 4 to 8 metaballs in nebula_shape_balls (centre and radius,
  no mass, no velocity) and the ellipsoid, warp, noise and isovalue
  columns (GEN.75, galaxy/nebula_shape.py); it moves rigidly and never
  changes shape. Done: a design note
  (docs/design/nebulae-as-point-groups.md, with its
  docs/design/README.md row) that answers: (1) the points: each a
  SpatialPosition3D with mass and vector (cloud parcels), how many per
  nebula by class (the 4 to 8 balls today, or tens for fidelity), and
  the storage cost per count (seven doubles and a key is about 80 bytes
  a point, so 8 points over the 163,000 field clouds of classes M and N
  is about 100 MB [C], to be checked against
  nebula-density-vs-reality.md); whether the small P and Q clouds stay
  single points; (2) how the shape and bounds are derived at each
  update: the metaball field summed from the points' masses and radii,
  the isovalue, the bounding radius_ly and the containment tests
  (inside_nebula_id, the stars placed inside by GEN.241, the map's soft
  sprites MAP.176/MAP.178) recomputed, and when a stored shape is
  rebuilt rather than carried; (3) what point-mass gravity cannot
  capture and what stands in for it: gas pressure and turbulence,
  magnetic fields, radiation pressure and photoionisation, shocks, star
  formation consuming the gas, stellar feedback, and the expansion of
  supernova remnants and planetary nebulae (analytic radius against age
  rather than gravity); (4) the timescale: a giant molecular cloud of 30
  pc crossing at about 3 km/s takes about 10 Myr [C], galactic shear
  (Oort A about 15 km/s/kpc) shapes a cloud over a few tens of Myr [C],
  so shape change is invisible at one day a step and shows only on long
  steps (GEN.292) such as 1 Myr, about three hours at 100 years a
  second; (5) a recommendation, including whether the points orbit under
  GEN.109's influence rule or only the galaxy's potential (GEN.115) and
  a cap per nebula so the nightly cost stays bounded (PERF.103's row ceiling counts every point). Source:
  Boss's question and the coordinator's answer. Phase 3, after GEN.294,
  GEN.109 and GEN.292.
  Prerequisites: GEN.294, GEN.109, GEN.292. Related: GEN.75, GEN.115,
  GEN.241, GEN.289, GEN.293, MAP.176, MAP.178.

- [ ] **GEN.296 The uncharted-star information panel shows the star's metallicity [Fe/H], which is already stored on bright_stars**
  Finding from Core 2 (2026-10-16, 4:58 AM EDT, GEN.158 merged as
  #1789): the sector map's star information panel now shows "Metallicity
  [Fe/H]", but the panel for an uncharted star (a bright star not yet
  built into a system) does not, although the value is stored on
  bright_stars (migration 0095). Done: the uncharted-star panel (and the
  uncharted phenomenon page of MAP.211 where a star is shown) shows the
  same Metallicity [Fe/H] line with the same format as the built-star
  panel, with a test. Owner: unassigned (small, fits Core 2 beside
  GEN.159).
  Prerequisites: none. Related: GEN.158, MAP.211.

- [ ] **GEN.297 Calibrate the cluster generator against the Harris file: its distributions guide the laws, and no real cluster is placed in any galaxy**
  Boss (2026-10-16, 4:11 PM EDT, in the TODO thread): "There should be
  no option to pick a Milky Way preset. The Harris file is designed to
  give you an idea of how to build the procedure to generate a realistic
  galaxy." Follow-up from Core 2 (2026-10-16, 7:52 AM EDT, GEN.159
  merged as #1810, migration 0097): the cluster table
  (galaxy_components), the King tables and solver and a seeded generator
  are built. Done: the Harris file (kept at
  /mnt/project-files/data/harris/mwgc.dat) is read by a reference reader
  used by the tests and the design note, never loaded into a galaxy; the
  generator's laws (the count from the halo and bulge mass, the mass
  function, the radial profiles of the halo and bulge groups, the
  bimodal [Fe/H], the King concentration, the age spread) are calculated
  or derived from published science and checked against the file's
  distributions statistically, not row for row; a test shows the
  generated population matches the file's medians and spreads within
  stated bounds; the design note records the file's version, licence and
  the mass-to-light assumption. There is no Milky Way option and no real
  cluster (47 Tuc, Omega Centauri, M13) appears by name. Owner: Core 3,
  in the cluster line. The file arrived 2026-10-16, 3:04 PM EDT. The
  file: Harris 1996 (2010 edition), December 2010 revision, 157 clusters
  in three parts. Part I has the name, RA and Dec, galactic longitude
  and latitude, the distance from the Sun and from the Galactic centre
  (R0 = 8.0 kpc) and X, Y, Z; Part II the [Fe/H], reddening, V magnitude
  of the horizontal branch, distance modulus, integrated V, M_V, the
  colours, spectral type and ellipticity; Part III the radial velocity,
  dispersion, King concentration c (a "c" marks a core-collapsed
  cluster), core and half-light radii in arcminutes, central surface
  brightness and luminosity density, and the relaxation times. It has no
  mass and no age: the reader takes mass from M_V with a stated
  mass-to-light ratio and says so; gaps in a column fall back to the
  generator's own laws. Licence: free to use if the original source is
  named; cite "Harris 1996 (2010 edition)" in the design note and the
  file header. GEN.237 merged (#1851, 2026-10-16, 3:15 PM EDT): the
  component table exists.
  Boss (2026-10-16, 4:14 PM EDT, in the TODO thread): "I want to
  reiterate that the idea is that we generate a scientifically accurate
  galaxy as best we can. We generate the sectors of that, again as
  accurately as we can, with what we can actually do processing-wise and
  what we actually know. Astrophysics-wise that can be implemented I
  don't want any fixed content. Everything should be procedurally
  generated based on the random seeds that we already use. Random
  numbers are used where we can. Things are derived everywhere we can or
  calculated everywhere we can. If needed you should ask me to do
  research if there are pieces missing that you can't look up yourself."
  Prerequisites: GEN.237. Related: GEN.160, GEN.237, GEN.162.

- [ ] **GEN.298 Remnants read the superbubble: the remnant size at build and the remnant factor use the bubble a supernova sits in**
  Follow-up from Core 3 (2026-10-16, 8:40 AM EDT, GEN.254 superbubble
  field merged as #1814, route GET /api/galaxy/superbubbles): the
  supernova remnant size at build (the density-dependent Sedov-Taylor
  law of GEN.151) and the remnant factor do not yet read the superbubble
  a star sits in, so a remnant inside a bubble is sized as if in the
  ambient gas. Also noted: clouds that fall inside a bubble are lost
  (the CLI step --redo-nebula-field apply removes stored ones). Done: a
  remnant inside a bubble uses the bubble's interior density in its size
  law, with a test over a seeded bubble, and the design note records the
  rule. Owner: unassigned (Core 3's area).
  Owner: Core 3 (Boss, 2026-10-16, 2:42 PM EDT: Core 3 takes the 8.1 new
  generate scatters line after GEN.280; the multi-star train follows its
  spine, 3:00 PM EDT). This replaces any earlier Core or Core 1 owner
  line on this item.
  Prerequisites: none. Related: GEN.254, GEN.151, GEN.256.

- [ ] **GEN.299 The Quantity class (the unit class): one value, every unit version always available, all conversions and the shared ladders**
  Boss (2026-10-16, 11:35 AM EDT, project chat): "Add new item to Phase
  3, a "unit" class, similar to what we do with the Point in space
  class, so that whatever kind of unit it is (length, mass, luminosity,
  etc...) the class itself makes sure that every version of that data is
  available all the time, and it does all the conversions, uses the
  shared ladders, etc"
  Design and build. A `Quantity` class (the unit class), the way
  SpatialPosition3D (physics/position.py, GEN.74, GEN.294) is the
  point-in-space class: one object per measured value, built from a
  value in any unit of its kind, with every other version of that value
  (the SI base value and each rung of the kind's ladder) available at
  every moment as an attribute, kept in step when the value changes, so
  no caller converts by hand and no copy of a unit version can go stale.
  The kinds are the ladder families the repository already shows
  (length, mass, luminosity, time, temperature, pressure, speed,
  density, gravity, energy and the rest of
  docs/design/units-and-number-formatting.md); the class reads the
  shared ladders (util/ladders.py, the one data table that
  scripts/build_unit_ladders.py also writes for the browser) for its
  units, conversion factors and the rung a value is shown in, and the
  browser gets the same class from the same table, so Python and
  JavaScript cannot disagree. It does the conversions between kinds that
  are really one dimension (a mass as a gravitational parameter,
  SpatialPosition3D's length_unit_m and mu_of and mass_of), refuses a
  conversion between kinds that are not (length to mass), takes a value
  deterministically (GEN.135) and formats through format_quantity.
  Default decisions (they stand unless Boss says otherwise): stored
  columns stay one column per quantity in the SI base unit, with no
  second column or JSON for a unit version; every other version is
  computed from the class on read and is never stored. Kinds are added
  to the ladder table, not declared on the class, so a new kind is one
  table row.
  Done: the class and its kinds, all conversions and the formatting
  algorithm in Python and in the browser from the one ladder table,
  SpatialPosition3D reading its mass and length through it, a property
  test that every kind round-trips through every rung and that the
  Python and JavaScript copies agree, and the design note updated.
  Prerequisites: GEN.294, GEN.135. Related: GEN.74, GEN.109, UX.107,
  UX.108.
  Design: [docs/design/units-and-number-formatting.md](design/units-and-number-formatting.md)

- [ ] **GEN.300 Every quantity in the code moves onto the Quantity class, and hand conversions are removed**
  Boss (2026-10-16, 11:35 AM EDT): the Quantity class (GEN.299) makes
  every version of a value available and does all the conversions, so
  the places that still convert by hand move onto it. Done: every call
  site that multiplies by a unit factor, calls format_quantity or keeps
  a value next to its unit (the describer, the map side panels, the CLI
  query lines, the class page ranges, the wide-binary, comet, nebula and
  quasar numbers that UX.107 sends through the ladders, the API output,
  the stored-column readers) takes a Quantity instead; no hand
  conversion factor remains in src outside the ladder table, and a test
  greps for one so a new one fails the suite (no compatibility layer,
  the old helpers are removed). Placed after GEN.299 and UX.107 because
  it is a sweep of call sites, a different timing and owner from the
  class itself.
  Split (2026-10-16, Boss 1:07 PM EDT: look at each item for smaller
  chunks and shared work): this item is now the umbrella and is done
  when its subitems are: GEN.348, UX.128, GEN.349. 'Every quantity in
  the code' is a sweep over describer, map panels, CLI lines, class
  pages, the API and the stored-column readers, plus a grep test; that
  is several PRs, in two languages, and the grep test can only be turned
  on at the end.
  Prerequisites: GEN.348, UX.128, GEN.349. Related: GEN.294, UX.108.

- [ ] **GEN.301 Headline of 8.1: the new generate scatters, accurate galactic generation and scattering including the halo**
  Boss (2026-10-16, 1:07 PM EDT): "I also want the new generate scatters
  as a headline of 8.1 as well as infinite zoom. Accurate galactic
  generation and scattering including the halo." Umbrella and the second
  headline of 8.1 (the first is the fly-through Galaxy Map, MAP.146).
  Done when the three scatters fill the whole sphere out to the 262 kpc
  reach from one rule set and a per-sector profile (GEN.240), with the
  satellites, clusters and wanderers placed as components (GEN.237,
  GEN.238, GEN.244), every nebula built with real stars (GEN.241 to
  GEN.243), a realistic galactic core with real density and a per-cell
  cap (GEN.280 to GEN.285), the halo luminosity floor tier (GEN.283) and
  its sliders (UX.116), the stellar evolution fixes that keep the bright
  giants honest (GEN.286, GEN.287, GEN.259), and the audits that show
  every kind of object is represented at its expected frequency
  (TEST.159, TEST.160).
  Headline: 8.1 (the new generate scatters: accurate galactic generation
  and scattering, including the halo).
  Prerequisites: GEN.238, GEN.240, GEN.241, GEN.242, GEN.243, GEN.244,
  GEN.259, GEN.280, GEN.282, GEN.283, GEN.284, GEN.286, GEN.287,
  GEN.298, UX.116, TEST.159, TEST.160.
  Design: [docs/design/intergalactic-and-halo-contents.md](design/intergalactic-and-halo-contents.md)

- [ ] **GEN.302 A region object: a sphere, box, shell, band, wedge or set of sectors, and the systems and phenomena inside it**
  Boss (2026-10-16, 1:07 PM EDT), headline of 8.6 (ADM.57): the admin
  has complete control; the world adjusts to what the admin wants, after
  the admin approves the generation plan. One library type for "an area
  the admin chose": its shape and parameters, a stable text form for
  plans and the density table, a membership test, and the list of stored
  objects and unbuilt scatter cells inside it, read through the cell key
  (GEN.276). Done: the type with tests for every shape, containment edge
  cases and the cover of a region by cell key runs.
  Scope after the review (2026-10-16 review): only the new shapes
  (sphere, box, set), the union and the stable text form, and the object
  and scatter membership by cell key; shell, band, wedge and ball come
  from `planetgen/galaxy/span.py` (`Span`) and
  `planetgen/galaxy/region.py` (`Region`).
  Prerequisites: GEN.276. Related: ADM.57, NAV.60.
  Design: [docs/design/edit-everything.md](design/edit-everything.md)

- [ ] **GEN.303 Per-cell system cap as a configuration option: `generation.cell_cap`, default 5,000, lowest 500, highest 50,000**
  Boss (2026-10-16, 1:47 PM EDT, in the GEN.280 thread): "cell cap
  should be a configuration option." The per-cell cap of systems
  (GEN.281, now the adaptive cap that GEN.280 builds) becomes one named
  option `generation.cell_cap` in the ADM.42 model: default 5,000,
  lowest 500, highest 50,000, with help text that says what it limits
  and that changing it affects only later fills, not stored cells. The
  fill code reads it at the start of a run instead of a constant, and
  the settings page of ADM.84 lists it with no extra screen. Core 3
  keeps the cap in one named setting so this item only wires it. The
  bounds are defaults for Boss to change. Done: the option, the read at
  run start, a test that a run with a lower cap stores fewer systems per
  cell, and the entry in docs/config.md.
  Stored in the user's `config.json` like every option (one
  configuration file, ADM.96).
  The cap it wires is the adaptive one of GEN.364 (Boss, 2026-10-16,
  1:46 PM EDT); the option sets that cap and nothing else. It waits on
  ADM.84 for the settings page and needs no new migration.
  Prerequisites: ADM.84. Related: GEN.281, GEN.283, ADM.43.

- [ ] **GEN.304 The canonical fingerprint encoder: floats to 9 significant digits, SHA-256, structural and float digests**
  Split from GEN.136 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). Replace the exact-float encoding
  of `planetgen/db/fingerprint.py` (GEN.58) and the
  `tests/galaxy_fingerprint.py` stand-in with one canonical encoder as a
  library: floats as `format(x, ".8e")`, integers exact (hex above
  2^53), no Unicode normalisation, sorted keys, generation order for
  lists, no ids or timestamps, a `pgfp1` format tag, SHA-256; separate
  structural and float digests over named sections (system, star,
  planet, moon, belt, comet, phenomenon, population). Why: OPS.28's
  battery digest, TEST.77's tiers and GEN.57's order test all hash
  generated values and must hash them one way, and a 1-ulp difference in
  exp or sin between machines must not show. Done: the encoder with
  golden-vector tests (the same value gives the same string on Python
  3.12 and 3.13; a 1-ulp change in a float leaves the digest unchanged
  and a change in the 9th digit changes it; the structural digest
  ignores float noise); `db/fingerprint.py` uses it; the stand-in is
  deleted and its tests moved.
  Prerequisites: GEN.135. Related: GEN.136.

- [ ] **GEN.305 Store a leaf digest at save time and build the tree: roots and region digests**
  Split from GEN.136 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). Compute the leaf digest at save
  time from the in-memory sector (not a re-read) with child a's encoder,
  store it, and update that sector's band and shell nodes (one
  1,000-child rehash instead of a scan); galaxy and shell roots on
  demand; a region digest that includes the sector count; `planetgen
  fingerprint` prints them. DB.9 keeps its own storage checksum (a
  separate thing). Done: tests: changing one sector changes exactly its
  band node, shell node and the root; a region digest changes when a
  sector is added; the digests do not depend on save order or worker
  count.
  Prerequisites: GEN.304, DB.49. Related: GEN.136.

- [ ] **GEN.306 Population seeds keyed on the object, not on row ids**
  Split from GEN.57 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work). `population/model.py` draws a
  species' form from `Stream(species_id * 7919)`, a planet's life from
  `Stream(planet_id)` and its technology from
  `Stream(f"tech:{homeworld_planet_id}")`; row ids come from per-worker
  id blocks, so the population changes with the worker count and the
  fill order. Key every one of these streams on the object's unique ID
  (uid) or the sector address plus the object's index. Changes stored
  values, so it ships with a generator epoch bump (OPS.28 child a).
  Done: the same galaxy gets the same species, ages and technology
  indices at 1 and 4 workers and in two fill orders; a test says so; no
  stream is seeded from a database id.
  Prerequisites: OPS.56. Related: GEN.57.

- [ ] **GEN.307 The object-ID collision number keyed on the lower address**
  Split from GEN.57 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work). `db/store._claim_object_ids` gives
  the next collision number to an ID already held by a stored row, so
  when two sectors' objects share an ID cell the number depends on which
  sector saved first. Make it depend on addresses: the object in the
  sector with the lower address keeps number 0, the other takes the
  next. Check first whether the 96-bit ID work (the later phases)
  already removes the collision; if it does, close this child as not
  needed. Changes stored values, so it ships with an epoch bump. Done:
  two sectors whose objects share an ID cell get the same numbers in
  either save order; a test says so.
  Prerequisites: OPS.56. Related: GEN.57.

- [ ] **GEN.308 A block-first sampler for the backfill: a Poisson candidate count from the block's density bound, then thinning**
  Split from GEN.42 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work). In `generation/bright_stars.py`, a
  function that draws a block's candidates as one Poisson count from the
  block's density bound and thins each to its own cell's rate, exact (no
  cell's probability altered; `MIN_RELATIVE_DENSITY` keeps every cell
  above zero), keyed so each cell's result depends only on the seed, the
  address and the mass band (GEN.39). Not wired to the backfill yet.
  Why: this is the part GEN.40 and PERF.35 depend on and it is testable
  alone. Done: the function and its tests: an assertion that every
  cell's rate is at most its block bound; the stratified multiplicity
  chi-square (half-decade strata of lambda, 3,000 or more replicates)
  run on the real draw; GEN.43's floor test as the lowest-lambda strata;
  the bound uses the `detmath` helpers.
  Prerequisites: none. Related: GEN.42.

- [ ] **GEN.309 The backfill picks the draw path per block and gives the same stars wherever it lands**
  Split from GEN.42 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work). In the backfill driver, choose per
  block: the candidate path (child a) when the expected candidates are
  under about 5 percent of the block's cells, the per-cell path
  otherwise. Done: the pass runs inside the backfill, the backfill
  visits fewer sectors, and a test of order, partition and sub-region
  identity shows the same stars for the same seed whichever way a region
  is reached or split; PERF.31 measures the SQL part at scale (2 ms per
  sector warm, a `sector_stats` row per visited sector) and the number
  goes in the PR.
  Prerequisites: GEN.308, PERF.108. Related: GEN.42.

- [ ] **GEN.310 Fill a whole block of the shell ladder: planetgen generate --block and the queue job**
  Split from GEN.231 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). A block of the new ladder (its m,
  shell, band and wedge, docs/design/spherical-block-ladder.md) is a
  region that the generate command and the work queue can fill, nearest
  sectors first, with the same progress and cancel as a neighbourhood
  fill (there is no resume); `--block` takes the new block address; the
  work-queue job route (`api/workqueue.py`) accepts a block so the map
  can start it; the refusal text in `cli/generate.py` and its test in
  `test_gen_cli_errors.py` are deleted. Done: `planetgen generate
  --block <address>` fills the block on the new queue and a cancel stops
  it; an API job for a block starts and reports progress; tests for the
  command and the route.
  Prerequisites: PERF.99. Related: GEN.231.

- [ ] **GEN.311 Built systems store and load their orbit tree: insert_star_system, load_star_system and serialization read and write orbit_nodes**
  Split from GEN.129 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). Once orbit_nodes exists, a
  StarSystem has to put its tree there and read it back. Write the nodes
  in store.insert_star_system, read them in load_star_system and the
  serialization, N=2 first (a plain binary is one node) with the old
  binary columns kept filled until every reader moves, and make
  positions_at read the nodes the way physics/system_tree.py already
  does for a scatter row. Done: a generated binary, a triple and a
  seven-star tree round-trip through the database with equal nodes,
  positions_at gives the same positions before and after, and the old
  columns still equal the node for N=2.
  Owner: Core 3 (Boss, 2026-10-16, 2:42 PM EDT: Core 3 takes the 8.1 new
  generate scatters line after GEN.280; the multi-star train follows its
  spine, 3:00 PM EDT). This replaces any earlier Core or Core 1 owner
  line on this item.
  Prerequisites: DB.53. Related: GEN.129.

- [ ] **GEN.312 Naming for trees: path codes (Aa, Ab, Ba) with branch words at wide splits, and planet stems such as Pikkita Ba I**
  Split from GEN.129 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). GEN.62 names a primary and one
  companion; a tree needs a name per member and per planet host. Extend
  it with path codes from the tree (Aa, Ab, Ba), branch words at wide
  splits drawn from the seeded stream in a fixed order, and planet stems
  that carry the host's code. Done: every member of a generated tree has
  a unique, deterministic name, two builds of one seed give the same
  names, and an old two-star system keeps its names.
  Owner: Core 3 (Boss, 2026-10-16, 2:42 PM EDT: Core 3 takes the 8.1 new
  generate scatters line after GEN.280; the multi-star train follows its
  spine, 3:00 PM EDT). This replaces any earlier Core or Core 1 owner
  line on this item.
  Prerequisites: none. Related: GEN.129.

- [ ] **GEN.313 Orbital update of built systems: advance_orbital_phases for nodes and the GEN.109 hooks**
  Split from GEN.129 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). The nightly update advances the
  mean anomaly of each built orbit. Make store.advance_orbital_phases
  advance the nodes of orbit_nodes (not the old binary columns) and let
  the Hill-sphere rule of GEN.109 read a system's total mass. Done:
  after an update every node's mean anomaly has advanced by 360 times
  the elapsed time over its period, a scatter row's nodes are still
  never written, and a binary behaves as before.
  Owner: Core 3 (Boss, 2026-10-16, 2:42 PM EDT: Core 3 takes the 8.1 new
  generate scatters line after GEN.280; the multi-star train follows its
  spine, 3:00 PM EDT). This replaces any earlier Core or Core 1 owner
  line on this item.
  Prerequisites: GEN.311. Related: GEN.129.

- [ ] **GEN.314 Planet zones in a hierarchy and flux sums: S-type and P-type windows per node (Holman-Wiegert) and the flux a planet receives from every member**
  Split from GEN.129 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). A planet's stable orbits in a
  tree: an S-type limit around each star (the sibling subtree as a point
  mass, lower limit the engulfment floor), a P-type window around each
  pair node (lower limit Holman-Wiegert, upper limit the S-type limit of
  the pair against its own sibling), all capped at the disk size, plus
  the summed flux of every member for the habitability inputs (design
  multistar-and-compact-systems.md 5.6 and 6.5). A pure function over an
  OrbitNode tree (physics/, numpy only), using the existing
  holman_wiegert functions of physics/orbits.py. This is the one place
  the zones are worked out; GEN.264:a, GEN.130 and GEN.275:b call it.
  Done: the windows for N=2 equal today's circumbinary and wide-pair
  limits, the table of 5.6 (S-type windows above 0.3 AU: 77% at N=2
  falling to 13% at N=7) is reproduced on a seeded draw, and a planet's
  flux sums over all members.
  The zone function is shared: GEN.264 stage C, GEN.130 (pulsar and
  white-dwarf rows) and GEN.275 (circumbinary planets) call it and
  compute no zones of their own.
  Owner: Core 3 (Boss, 2026-10-16, 2:42 PM EDT: Core 3 takes the 8.1 new
  generate scatters line after GEN.280; the multi-star train follows its
  spine, 3:00 PM EDT). This replaces any earlier Core or Core 1 owner
  line on this item.
  Prerequisites: none. Related: GEN.129.

- [ ] **GEN.315 A fill builds a system from its stored members and nodes; the system pass keeps only planets, moons, belts and disks**
  Split from GEN.264 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). StarSystem(member_tree=...) skips
  its own companion draw. Stages A and B build the system from the
  stored members (the placed object is member 0, the others are built
  from seed and initial mass at the system age) and copy the nodes into
  orbit_nodes; stage C draws planets, moons, belts and disks only (no
  planets under 10 Myr, millisecond-pulsar planets at 0.7%, zones from
  the stored nodes through GEN.129:f); the build deletes the row's
  scatter_members and scatter_nodes in its own transaction; a failed
  build drops the outermost body, then the outermost member, never the
  placed object. Done: the stars of a built system equal the stored
  members, its nodes equal the stored nodes, the child rows are gone,
  and the system pass never draws a companion for a scattered object.
  Owner: Core 3 (Boss, 2026-10-16, 2:42 PM EDT: Core 3 takes the 8.1 new
  generate scatters line after GEN.280; the multi-star train follows its
  spine, 3:00 PM EDT). This replaces any earlier Core or Core 1 owner
  line on this item.
  Prerequisites: GEN.311, GEN.314. Related: GEN.264.

- [ ] **GEN.316 Stages D and G draw a member tree for every object the fill creates, and the two-star classes are retired**
  Split from GEN.264 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). A system made at fill (stage D,
  the cell's other systems; stage G, below-cut remnants) must have
  companions at the same odds as a scatter row, from the same function.
  Call draw_member_tree for each, replace BinaryStarProxy,
  WideBinaryPair and BINARY_SYSTEM_PROBABILITY_BY_SPECTRAL_CLASS with
  the N table (B and A binaries raised to about 0.75 and 0.60, table
  mean 1.38 stars per system). Done: a fill-created system draws its
  members from draw_member_tree at the table's odds, no code path uses
  the old two-star draw, and the census of stars per system matches the
  table.
  Owner: Core 3 (Boss, 2026-10-16, 2:42 PM EDT: Core 3 takes the 8.1 new
  generate scatters line after GEN.280; the multi-star train follows its
  spine, 3:00 PM EDT). This replaces any earlier Core or Core 1 owner
  line on this item.
  Prerequisites: GEN.315. Related: GEN.264.

- [ ] **GEN.317 Pulsar planets and the radiation environment of neutron stars**
  Split from GEN.130 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). The system pass already rolls a
  millisecond pulsar's planets at 0.7% (tuning.SYSTEM_PASS_ODDS). Build
  what they are: three types (disc 25%, ablated companion 50%, captured
  circumbinary giant 25% present as a type but switched off until a
  cluster says otherwise, GEN.163), carbon rich with no atmosphere and
  radiation tier 3, the planet model's cap of about 10 planets for a
  recycled pulsar corrected, planets placed in the zones of GEN.129:f;
  and the radiation fields (spin-down power, X-ray luminosity) feeding
  GEN.86 and GEN.87. Done: a millisecond pulsar that rolls planets gets
  them in those types and shares, with the stated radiation tier, no
  young pulsar has any, and a test pins the 0.7% and the 25/50/25 split.
  Owner: Core 3 (Boss, 2026-10-16, 2:42 PM EDT: Core 3 takes the 8.1 new
  generate scatters line after GEN.280; the multi-star train follows its
  spine, 3:00 PM EDT). This replaces any earlier Core or Core 1 owner
  line on this item.
  Prerequisites: GEN.314. Related: GEN.130.

- [ ] **GEN.318 White-dwarf and planetary-nebula-star rows of the system pass: engulfment of inner planets and the planetary nebula's circumbinary disk**
  Split from GEN.130 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). The white dwarf and the planetary
  nebula's central star take their planets from the star model today and
  have no row in tuning.SYSTEM_PASS_ODDS, and a planetary nebula's
  circumbinary disk is not modelled (Core 1's note on GEN.251). Add the
  rows: planets around a white dwarf with engulfment removing the inner
  ones, the central star's close companion as already drawn at scatter
  (20%), and the circumbinary disk as the disk descriptor of GEN.129:a
  with no bodies. Done: a planetary nebula's white dwarf has its row, an
  engulfed planet is absent, a circumbinary disk is stored and shown,
  and the white-dwarf fractions of the catalogue hold within Poisson
  bounds.
  Owner: Core 3 (Boss, 2026-10-16, 2:42 PM EDT: Core 3 takes the 8.1 new
  generate scatters line after GEN.280; the multi-star train follows its
  spine, 3:00 PM EDT). This replaces any earlier Core or Core 1 owner
  line on this item.
  Prerequisites: DB.53, GEN.314. Related: GEN.130.

- [ ] **GEN.319 Nebula rates read the profile: young and old shares, the cell's gas and the Type Ia share as profile fields**
  Split from GEN.241 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). The class rates of the nebula
  field and the scatter must follow the place: emission, reflection,
  Wolf-Rayet and core-collapse remnants scale with the young share of
  the profile's mix, planetary nebulae and Type Ia remnants with the old
  share, molecular clouds and protostars with `gas`; the Type Ia remnant
  share becomes a profile field (1 where there is no young population);
  nebula_field reads the cell's gas, not the disk position alone. Done:
  an asked halo cell never shows an H II region or a pulsar-bearing
  remnant, a Magellanic cell can show every class, and the rates of the
  Milky Way default are unchanged inside the disk (a test compares
  them).
  Owner: Core 3 (Boss, 2026-10-16, 2:42 PM EDT: Core 3 takes the 8.1 new
  generate scatters line after GEN.280; the multi-star train follows its
  spine, 3:00 PM EDT). This replaces any earlier Core or Core 1 owner
  line on this item.
  Prerequisites: none. Related: GEN.241.

- [ ] **GEN.320 build_nebula_with_contents: every hosted nebula gets its stars as real systems inside its bounds (classes C, G, F, planetary, remnants, P and Q)**
  Split from GEN.241 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). One function,
  build_nebula_with_contents, that every path creating a hosted nebula
  without a host calls: emission C gets its O or early-B star at the
  centre, G the early-B star, F the B3 to A9 star, planetary nebulae
  keep their white dwarf and gain the close companion (20%), supernova
  remnants keep their neutron star or black hole at 60% of core-collapse
  remnants, Bok globules P one protostar in 30%, star-forming cores Q
  one to four protostars (Poisson mean 2, at least 1) clumped near the
  centre; M and N none. A nebula that cannot get its stars is not made.
  Added as stage F in generation/fill_stages.py on its own stream; the
  hosted roll for a pre-placed star is decided from the star's own seed
  (folds MAP.182). Done: every hosted nebula of these classes holds its
  stars inside its radius, a planetary nebula always has its white
  dwarf, and a pre-placed star's nebula shows from its scatter row.
  Needs the stored nebula age and ionizing rate (child a of GEN.95),
  which the member draw reads for the age cap (2026-10-16 review).
  Owner: Core 3 (Boss, 2026-10-16, 2:42 PM EDT: Core 3 takes the 8.1 new
  generate scatters line after GEN.280; the multi-star train follows its
  spine, 3:00 PM EDT). This replaces any earlier Core or Core 1 owner
  line on this item.
  Prerequisites: GEN.242. Related: GEN.241.

- [ ] **GEN.321 Emission regions D and E: the ionizing members of a young cluster in a Plummer sphere, cut by sector cube**
  Split from GEN.241 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). Classes D and E get the seed star
  plus the cluster's ionizing members (stars of 5.5 solar masses and up
  still alive: 42 on average for D, 549 for E) in a Plummer sphere of
  half-mass radius 1 pc (D) and 2 pc (E) clipped to the nebula radius,
  derived from the nebula's own seed and placed by position so each
  sector generates only the members whose positions fall in its cube.
  Use the Plummer profile and sampler of GEN.237, not a second one.
  Default as the parent says: only the ionizing members are real;
  filling the whole cluster is a setting. Done: every D or E region
  holds its members inside its radius, a nebula across sectors needs no
  coordination (two sectors give the same members), and an asked massive
  complex costs about 270 extra systems and 7 seconds.
  Owner: Core 3 (Boss, 2026-10-16, 2:42 PM EDT: Core 3 takes the 8.1 new
  generate scatters line after GEN.280; the multi-star train follows its
  spine, 3:00 PM EDT). This replaces any earlier Core or Core 1 owner
  line on this item.
  Prerequisites: GEN.320. Related: GEN.241.

- [ ] **GEN.322 Band groups and the halo term: the bright-star and mass passes cover the 262 kpc reach**
  Split from GEN.240 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). The two passes keep their
  algorithm (one Poisson count per shell from a certified bound, bands
  in proportion, accept by true density over bound) but take the halo
  term in band_majorants, the new populations in fractions, and band
  groups of 64 bands beyond the outline; rows carry their companions
  with the primary and a row counts as 1 + p_placed. Done: the passes
  cover the reach with the same code, inside the outline they yield
  exactly the rows they yielded before, beyond 8 kpc under 10%
  efficiency is lost, a pass to 262 kpc takes about 2 minutes, and a
  pre-placed star is never drawn again at fill.
  Owner: Core 3 (Boss, 2026-10-16, 2:42 PM EDT: Core 3 takes the 8.1 new
  generate scatters line after GEN.280; the multi-star train follows its
  spine, 3:00 PM EDT). This replaces any earlier Core or Core 1 owner
  line on this item.
  Prerequisites: GEN.237. Related: GEN.240.

- [ ] **GEN.323 The phenomena pass runs on the object-first sampler**
  Split from GEN.240 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). The phenomena pass walks every
  band of every shell (_band_rows), which cannot reach 262 kpc. Move it
  to the same object-first sampler as a, with its rates, cut and
  below-cut draw at sector time unchanged and the cell key written from
  whatever position it writes. Done: the pass covers the reach, a seeded
  galaxy inside the outline gets the same phenomenon rows as before (a
  test compares), and about 2.4 million extra black-hole rows appear
  outside the disk at the defaults.
  Owner: Core 3 (Boss, 2026-10-16, 2:42 PM EDT: Core 3 takes the 8.1 new
  generate scatters line after GEN.280; the multi-star train follows its
  spine, 3:00 PM EDT). This replaces any earlier Core or Core 1 owner
  line on this item.
  Prerequisites: GEN.322. Related: GEN.240.

- [ ] **GEN.324 A component pass: clusters, the Clouds and Sagittarius place their stars, and the fill adds the component's density to lambda**
  Split from GEN.240 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). Components get their own pass:
  total star count from the mass, positions by inverse transform of the
  profile table, each star kept or dropped by the same pass cut; the
  sector fill adds the component's density to lambda so no pre-placed
  star is counted twice. Done: at the defaults the pass places nothing
  in clusters and old satellites and a few thousand stars in the Clouds,
  no cluster star is counted twice, and the yields match the design
  note's table.
  Owner: Core 3 (Boss, 2026-10-16, 2:42 PM EDT: Core 3 takes the 8.1 new
  generate scatters line after GEN.280; the multi-star train follows its
  spine, 3:00 PM EDT). This replaces any earlier Core or Core 1 owner
  line on this item.
  Prerequisites: GEN.322. Related: GEN.240.

- [ ] **GEN.325 Planet cull in dense components by the survival formula and the metallicity factor**
  Split from GEN.162 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). Giant-planet probability is
  multiplied by 10^(2 [Fe/H]), capped at 1, and a planet survives only
  inside the largest orbit that outlives 12 Gyr of encounters: expected
  ionizing passes n pi a^2 v t (1 + 2 G M / (a v^2)), v = sqrt(6) times
  the dispersion, M = 0.5 solar masses
  (docs/design/realistic-galactic-core.md section 7, decision 4: the
  survival formula, not the hard cut at a_h). One function,
  planet_survival_limit_au(density, dispersion), that takes the place's
  density and dispersion; here the cluster's King n(r) and sigma(r) feed
  it. Done: the table of the cluster note (188, 11, 2.6 and 2.1 expected
  passes at its four rows) is a test, a giant planet in a metal-poor
  cluster is thinned by the factor, and the same function takes any
  density and dispersion (GEN.284 feeds it the nucleus).
  Owner: Core 3 (Boss, 2026-10-16, 2:42 PM EDT: Core 3 takes the 8.1 new
  generate scatters line after GEN.280; the multi-star train follows its
  spine, 3:00 PM EDT). This replaces any earlier Core or Core 1 owner
  line on this item.
  Prerequisites: none. Related: GEN.162.

- [ ] **GEN.326 Blue stragglers in clusters: N_BSS scales as the core mass to the 0.4**
  Split from GEN.162 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). Blue stragglers are main-sequence
  stars of 1 to 1.7 solar masses drawn within the core radius, their
  number N_BSS about M_core^0.4 (nearly independent of the encounter
  rate at fixed core mass; Boss's default taken). Done: a cluster's blue
  straggler count follows the law within Poisson bounds, they sit in the
  core, and none appear outside clusters.
  Owner: Core 3 (Boss, 2026-10-16, 2:42 PM EDT: Core 3 takes the 8.1 new
  generate scatters line after GEN.280; the multi-star train follows its
  spine, 3:00 PM EDT). This replaces any earlier Core or Core 1 owner
  line on this item.
  Prerequisites: GEN.161. Related: GEN.162.

- [ ] **GEN.327 The group node and the remnant draw: offset and velocity columns, a decay age, and the lightest member ejected**
  Split from GEN.273 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). A node of kind group with N
  members, each with an offset (AU) and a velocity (km/s) at an epoch,
  and decay_age_yr; internal motion is a straight line, no Kepler
  elements, and the nightly update does not integrate it
  (physics/system_tree.py already flies a group's members straight).
  Default 5% of massive (O and early B) systems younger than 0.3 Myr at
  scatter time for the few mass-pass stars of that age. A group older
  than its decay age is never drawn: the generator draws its result, the
  lightest member ejected as a runaway row (runaway_origin dynamical)
  and a tighter hierarchical remnant. Stored on scatter rows and built
  systems (the columns of GEN.129:a); the sector page shows a group with
  its remaining life. Done: every stored group is younger than its
  decay_age_yr (audit test 14), a group older than that is replaced by
  its remnant, and members carry offset and velocity.
  Owner: Core 3 (Boss, 2026-10-16, 2:42 PM EDT: Core 3 takes the 8.1 new
  generate scatters line after GEN.280; the multi-star train follows its
  spine, 3:00 PM EDT). This replaces any earlier Core or Core 1 owner
  line on this item.
  Prerequisites: DB.53. Related: GEN.273.

- [ ] **GEN.328 Groups at fill in star-forming cores and OB regions**
  Split from GEN.273 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). The same group draw, called when
  a fill builds a star-forming core or an emission region D or E (5% of
  massive systems younger than 0.3 Myr): the trapezia of the OB
  associations of GEN.244. Done: a seeded D or E region and a Q core
  hold groups at the default share, every one younger than its decay
  age, and a group across the cube boundary is made by one sector only.
  Owner: Core 3 (Boss, 2026-10-16, 2:42 PM EDT: Core 3 takes the 8.1 new
  generate scatters line after GEN.280; the multi-star train follows its
  spine, 3:00 PM EDT). This replaces any earlier Core or Core 1 owner
  line on this item.
  Prerequisites: GEN.327, GEN.321. Related: GEN.273.

- [ ] **GEN.329 Runaway bookkeeping: runaway_origin and the compact companion left by a supernova**
  Split from GEN.274 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). runaway_origin (none, supernova,
  dynamical) on a runaway at the observed rates (10% to 30% of O stars
  and 5% to 10% of B stars; two-thirds supernova, one-third dynamical
  for the massive ones); for the supernova channel a compact companion
  left behind at the rate the survival rule gives; companions unbound by
  a collapse (83% of non-interacting massive binaries) are not stored,
  the survivor is an ordinary row with the runaway flag; hypervelocity
  stars stay single. Done: runaway counts within the ranges (audit test
  10) and no runaway carries a missing origin.
  Owner: Core 3 (Boss, 2026-10-16, 2:42 PM EDT: Core 3 takes the 8.1 new
  generate scatters line after GEN.280; the multi-star train follows its
  spine, 3:00 PM EDT). This replaces any earlier Core or Core 1 owner
  line on this item.
  Prerequisites: DB.53. Related: GEN.274.

- [ ] **GEN.330 Cluster-hosted compact binaries: the cluster rates replace the field fractions**
  Split from GEN.274 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). In a cluster the compact-binary
  fractions follow the cluster, not the field: a black hole with a star
  10%, a millisecond pulsar with a companion 80% with a higher
  neutron-star share, and dynamical black-hole pairs from the cluster
  density; the cluster's soft binaries are cut by the survival rule of
  GEN.162:a applied to the pair. Done: a cluster sector draws its
  compact binaries at the cluster rates, a field sector at the field
  rates, and counts are within the ranges.
  Owner: Core 3 (Boss, 2026-10-16, 2:42 PM EDT: Core 3 takes the 8.1 new
  generate scatters line after GEN.280; the multi-star train follows its
  spine, 3:00 PM EDT). This replaces any earlier Core or Core 1 owner
  line on this item.
  Prerequisites: GEN.329, GEN.163, GEN.325. Related: GEN.274.

- [ ] **GEN.331 Substellar companions: brown-dwarf close (the desert) and wide companions, brown-dwarf pairs and planetary-mass companions**
  Split from GEN.275 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). Brown dwarf close companions
  (under 0.6% inside 3 to 10 AU, the desert), wide companions (3.2% at
  28 to 1,590 AU), brown dwarf pairs (15% of brown dwarfs) and
  planetary-mass companions (under 13 Jupiter masses, flagged
  planetary-mass), as member rows of the tree generator using GEN.252's
  state. Done: counts within four Poisson sigma of the defaults (brown
  dwarf wide companion 3.2%) and no brown dwarf in the desert beyond the
  stated rate.
  Owner: Core 3 (Boss, 2026-10-16, 2:42 PM EDT: Core 3 takes the 8.1 new
  generate scatters line after GEN.280; the multi-star train follows its
  spine, 3:00 PM EDT). This replaces any earlier Core or Core 1 owner
  line on this item.
  Prerequisites: GEN.252. Related: GEN.275.

- [ ] **GEN.332 Circumbinary and S-type planets, circumbinary disks and rogue planet pairs at fill**
  Split from GEN.275 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). Circumbinary planets (at least
  10% of close eclipsing binaries) in the P-type window of the node from
  GEN.129:f, planets at fill only; S-type planets in wide pairs;
  circumbinary disks as a descriptor on the node (GEN.129:a), no bodies;
  rogue planet pairs. GEN.95 applies the nebula conditions. Done: a
  circumbinary planet's orbit lies outside the stable zone, a disk
  carries no bodies, and counts are within four Poisson sigma of the
  defaults.
  Owner: Core 3 (Boss, 2026-10-16, 2:42 PM EDT: Core 3 takes the 8.1 new
  generate scatters line after GEN.280; the multi-star train follows its
  spine, 3:00 PM EDT). This replaces any earlier Core or Core 1 owner
  line on this item.
  Prerequisites: DB.53, GEN.314, GEN.315. Related: GEN.275.

- [ ] **GEN.334 Pair closest-approach library: sweep, pair windows, window cap and chunking**
  Found in the split review (2026-10-16, Boss 1:07 PM EDT: look at each
  item for smaller chunks and shared work, 1:33 PM EDT: shared
  foundations). `physics/closest_approach.py`: the closest-approach
  `sweep(dp, dv, R, dt)` of collisions-and-mergers.md section 5.2 (miss
  distance as a vector first, then the entry time), the pair-window rule
  of 5.3 (test X against Y's stored segment over [max(e_X_prev, e_Y),
  e_X_now]), the window cap tau_max and the one-sector reach cap,
  differences taken from sector-relative values, and chunking of a long
  step with a chunk length for galactic curvature (0.1 Myr). One
  function answers 'do these two paths come within R in this window',
  for R = R1 + R2 (collision) or a system radius (encounter). Done:
  tests for head-on hit and time, arrival after the window, hit at the
  window end, separating, already overlapping, tangent graze, miss by
  1e-7, zero relative speed, zero dt, and 20,000 random cases against
  dense time sampling with no mismatch.
  GEN.110 and GEN.290 each narrowed to what is theirs: GEN.110 to the
  outcome rules and applying them, GEN.290 to choosing per system when
  to switch to the perturbed pass, the cost estimate and the seeded
  black-hole test (2026-10-16 review).
  Prerequisites: none. Related: GEN.110, GEN.290.

- [ ] **GEN.336 Compact and free objects hold a velocity vector and a mass: stars, black holes, neutron stars, rogue planets, interstellar comets, asteroid fields, nebulae, remnants**
  Split from GEN.294 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). These eight tables keep only the
  scalar `galactic_orbital_speed_kms`. Done: generation writes the
  vector and the mass (a field's or nebula's bulk mass, or a stated zero
  that the gravity sum skips by rule); each table reads back into one
  SpatialPosition3D through child a; rogue planets also carry
  sector-relative offsets so collision tests difference small numbers
  (absolute doubles are good to about 14 km at 8 kpc); the scalar column
  is removed from the table, its model, the store, describe, edits,
  sector_paths and the generators in the same PR (a small drop revision
  from Core 1, no compatibility layer); the tables leave the allow-list.
  Owner: Core 1 (Boss, 2026-10-16, 2:42 PM EDT: Core 1 takes the global
  top 5 blockers, starting right after GEN.237 merges; estimates from
  Core 1: GEN.335 merged about 7:30 PM EDT, DB.55 about 10:30 PM EDT).
  Migration revision (Core 1, 2026-10-16, 6:13 PM EDT, renumbered 6:22
  PM EDT after DB.54 merged as 0101): 0103, after 0102 (DB.55).
  Prerequisites: DB.55. Related: GEN.294.

- [ ] **GEN.337 Epoch, vector and mass for the rest: belts, quasars, comets, facilities, scatter points and members**
  Split from GEN.294 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). Belts and the quasar keep no
  epoch or vector; comets, facilities, scatter points and members keep
  no mass (members also no epoch). Done: each stores an epoch, a vector
  (a stated zero for a body that does not move, such as the quasar) and
  a mass (its own, the mass of what a facility sits on, or a stated
  zero), reads back through child a's definition, the allow-list is
  empty and the coverage test now fails for any table missing a part.
  Owner: Core 1 (Boss, 2026-10-16, 2:42 PM EDT: Core 1 takes the global
  top 5 blockers, starting right after GEN.237 merges; estimates from
  Core 1: GEN.335 merged about 7:30 PM EDT, DB.55 about 10:30 PM EDT).
  Prerequisites: DB.55, GEN.336. Related: GEN.294.

- [ ] **GEN.338 The galactic potential as one pure module: Hernquist bulge, Miyamoto-Nagai disk and NFW halo, with settings and the rotation-curve tests**
  Split from GEN.115 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). One function of position returns
  acceleration, potential and circular speed from a Hernquist bulge, a
  Miyamoto-Nagai disk and an NFW halo (the flattened log halo only if
  its v0 and R_c are corrected and tested). Masses and scales are
  settings with the document's values at the default galaxy shape
  (lengths follow the disk scale length of 2,600 pc, bulge 1,580 pc,
  radius about 15.3 kpc; masses scale by s^p with p = 1); Plummer
  softening of 1 pc applies to the central black hole only; the
  acceleration is finite at (0, 0, 0). Why: GEN.109, GEN.198 and the
  stepper read it and no generator or table has to change for it. Done:
  the module and its settings with tests for a circular speed of 229.3
  +- 0.1 km/s at 8.128 kpc, 213 to 231 km/s on a 2 kpc grid from 4 to 20
  kpc, Kz(1.1 kpc) about 75 Msun/pc2 and local density about 0.10
  Msun/pc3 at 5 percent, a circular orbit keeping its radius to 1 pc
  over 5 Gyr at dt <= 1 Myr, Lz conserved to 1e-10, energy drift under
  1e-3 per Gyr, and a finite pull at the origin. No database or
  generator change.
  Drops the duplicate change of the Sun distance: GEN.144 (first in the
  Galactic gravity thread) names the 8.128 kpc test radius and the 8.2
  kpc placement value, and this module only reads them (2026-10-16
  review).
  Prerequisites: GEN.144. Related: GEN.115.

- [ ] **GEN.339 Generated and stored galactic velocities follow the potential: the old curve is replaced, dispersions are seeded, stored velocities are reseeded**
  Split from GEN.115 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). Replace
  `galactic_orbit.calculate_galactic_orbit` and its two constants
  (GALACTIC_ROTATION_FLAT_VELOCITY_KMS,
  GALACTIC_ROTATION_CORE_RADIUS_PC) with child a's circular speed; seed
  each object as v_c(R) along the tangent plus a peculiar part (Jeans
  asymmetric drift and the dispersion table of galactic-potential.md
  section 8; runaways keep `runaway_speed_kms`; a stand-alone facility
  takes v_c(R) with no dispersion); reseed or migrate the stored
  velocity vectors of stars, systems, phenomena and facilities. Why: the
  old curve gives 206 km/s at the Sun and would leave every star
  eccentric in the new potential. Done: a new galaxy's stars at 8.128
  kpc average 229.3 km/s with the stated dispersions, an old galaxy run
  through the reseed step gives the same distribution, and the old
  function and constants are gone.
  Prerequisites: GEN.338, GEN.294. Related: GEN.115.

- [ ] **GEN.340 One star stepper: kick-drift-kick in the potential, a step per object, on arrays**
  Split from GEN.115 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). Replace the rotate-about-Z scheme
  (`_galactic_turn`) in `advance_galactic_positions`, exact only for a
  planar circular orbit, with a kick-drift-kick step through child a's
  acceleration: step eta 0.05 of the pericentre crossing time, capped at
  5 Myr, block steps, a step cap that sets a flag, an extra-acceleration
  hook so GEN.109 adds the influence set to the same step, and the
  rotation kept for a facility with no vertical offset. Written on
  arrays from the start, which is the NumPy batching GEN.212 asks for.
  Done: a star at the Sun's radius keeps its radius to 1 pc over 5 Gyr,
  Lz is conserved to 1e-10, one century in one step agrees with a
  hundred one-year steps for stars without an encounter, and sector
  refiling is untouched.
  Prerequisites: GEN.338, GEN.294, GEN.139. Related: GEN.115.

- [ ] **GEN.341 Influence-set search: every object inside the largest nearby Hill sphere, found by sector bounds, with cover lists and a brute-force check**
  Split from GEN.109 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). Pure search, no motion. For a due
  object, find the objects whose point masses count: by distance from
  the sector's bounds (`geometry.enumerate_sectors_within_radius`, not
  the ring and slot stencil), the radius widened by the cell reach
  (about 4 pc; 21.6 percent of neighbours are missed otherwise),
  per-sector `max_mass_kg` and `max_mass_id` (columns from GEN.294:m)
  kept up to date, cover lists for objects whose Hill radius exceeds 1.5
  sector edges, `max_exact_members` (default 64) with one monopole per
  sector for the rest, a hierarchy as one point mass at its barycentre,
  a star cluster as its members, a lone rogue with an empty set. Done:
  tests compare the set against a brute force and `cKDTree` for sectors
  at the core, mid radius, the rim, the top and bottom layers and across
  the slot-0 seam, the exact pairwise Hill criterion, a negative control
  without the margin, the cap, and a seed-reproducibility test.
  Prerequisites: GEN.294, GEN.289. Related: GEN.109.

- [ ] **GEN.342 The acceleration kernel: the pull of the influence set plus the galactic gradient, as one function, and the velocity update from it**
  Split from GEN.109 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). One function returns the
  acceleration at a point from child a's set (pairwise point masses,
  Plummer 1 pc on the central black hole only, other masses unsoftened
  or 1e-3 to 1e-2 pc, no a_min filter) plus the galactic gradient of
  GEN.115's module, leaving out the object itself. For a star it is the
  extra term in GEN.115's stepper; for a planet or moon it works in the
  host's frame with the host's own pull subtracted. Trajectories change
  only when an update runs. GEN.198 later wraps this same kernel instead
  of summing masses a second time. Done: tests against a direct sum for
  random sets, momentum conserved for an isolated pair, a star's orbit
  unchanged to the stepper's tolerance when the set is empty, and finite
  everywhere.
  The force on a body comes from the shared field kernel (PERF.88); this
  item adds the influence-set term, the galactic gradient hook and the
  velocity update (2026-10-16 review).
  Prerequisites: GEN.341, GEN.338, GEN.340. Related: GEN.109.

- [ ] **GEN.343 The perturbed pass for a system: host, siblings, moons and outside perturbers integrated through an encounter**
  Split from GEN.109 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). capture-and-ejection-design.md
  section 3, Tier 2: an adaptive integrator (DOP853 or IAS15-class,
  REBOUND is not a dependency) through an encounter, in the host's frame
  with the host tide subtracted, moons with their planet, positions
  relative to the host, energy-based error control, the Roche limit and
  contact checked every sub-step. No item builds this today. Done: tests
  against a direct three-body integration, energy conserved to a stated
  tolerance in an isolated system, a head-on analytic flyby, far passes
  against the tidal approximation, and one 100-year step against a
  hundred one-year steps; a system with no perturber still takes the
  cheap clockwork.
  Prerequisites: GEN.342, GEN.289, GEN.140. Related: GEN.109.

- [ ] **GEN.344 Hill-radius warnings by the gravitational criterion, to the debug and activity logs**
  Split from GEN.109 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). A warning is raised only for a
  bound or capture candidate, a deflection of at least 1 degree, or the
  same system (collisions-and-mergers.md section 3.3), not for a plain
  Hill entry (about 200 a day on a full-size galaxy; plain entries go to
  debug). Hierarchy members are exempt from mutual warnings and
  captures, a flyby inside the outermost orbit's Hill zone is flagged
  'disrupts hierarchy', and an unbound flyby warns nothing. Done: a
  warning reaches the debug and activity logs with both ids and
  positions, a seeded pair inside each other's Hill sphere warns once
  and a hierarchy pair never does.
  Prerequisites: GEN.341, GEN.342. Related: GEN.109.

- [ ] **GEN.345 Collision outcomes as a pure function: regime table, giant merger, fusion check, conservation**
  Split from GEN.110 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). The section 5.1 outcome table of
  collisions-and-mergers.md as one function from two bodies (mass,
  radius, composition, impact speed, impact parameter) to an outcome and
  its products, seeded from the pair's ids and the window's start. A
  terrestrial body hitting a giant, brown dwarf or star is absorbed;
  hit-and-run and partial erosion leave both alive, eroded; two
  terrestrial bodies become an asteroid field of class U; two giants
  merge keeping momentum, with the new radius from
  `giant_radius_km(M_new)`; the fusion check (0.075 Msun, 78.6 Mjup, a
  setting) covers any substellar pair and yields a star flagged 'merger
  remnant' with no planets (the star generator takes 0.075 to 0.08 Msun
  or clamps to 0.08); `ALWAYS_DESTROY_TERRESTRIAL` is the setting for
  Boss's literal rule; the Roche factor 2.44 is for fluid bodies only
  (rocky 1.26). Done: a table-driven test for every row, mass and
  momentum conserved to a stated tolerance, the same input gives the
  same outcome, a borderline-within-3-percent flag; no database.
  Prerequisites: none. Related: GEN.110.

- [ ] **GEN.346 Applying a collision to the galaxy: new asteroid field, merged giant, merger-remnant star, retired bodies, events**
  Split from GEN.110 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). The update step that finds
  colliding rogues with the closest-approach library (differences taken
  from sector-relative values, the 5.3 pair windows and window cap,
  pairs that include a body already consumed in the run dropped, nothing
  deleted mid-loop) and applies child a's outcome: it creates the field
  (mass_kg, formed_at, origin_event_id and a velocity, columns from
  GEN.294:m, not a second migration), the merged giant or the star
  (reported to the admin), retires the destroyed bodies under the
  object-ID rules of GEN.197, and writes each event through GEN.143's
  record. Done: two seeded rogues on a collision course go through the
  update to a field, a merged giant and a star, each in the event
  record.
  Prerequisites: GEN.345, GEN.334, GEN.341, GEN.142, DB.56, GEN.197,
  DB.55. Related: GEN.110.

- [ ] **GEN.347 The last-collisions admin view and the test that runs the whole collision path**
  Split from GEN.143 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). The admin report of design
  section 7 and a 'last collisions' admin view over child a's table; the
  cross-links from nebula-and-asteroid-field-classes.md and
  interstellar-object-rates.md to collisions-and-mergers.md; a test that
  seeds two rogues on a collision course and runs the whole update
  through the field, giant-merger and star branches. Done: the view
  lists the three outcomes of the test and the test passes.
  Prerequisites: DB.56, GEN.346. Related: GEN.143.

- [ ] **GEN.348 Describer, CLI query lines and class page ranges take Quantity**
  Split from GEN.300 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). The Python text output: the
  describer, the CLI query lines, the class page ranges and the
  wide-binary, comet, nebula and quasar numbers that UX.107 sends
  through the ladders take a Quantity instead of multiplying by a unit
  factor or calling format_quantity beside a unit. Done: none of those
  call sites converts by hand, output is unchanged, the old helpers they
  used have no caller.
  Prerequisites: GEN.299. Related: GEN.300.

- [ ] **GEN.349 API output and stored-column readers take Quantity; the grep test fails on a hand-written conversion factor; the old helpers are removed**
  Split from GEN.300 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). The last call sites (the API
  output and the readers of stored columns), the removal of every helper
  the sweep made unused, and a test that greps src for a hand conversion
  factor outside the ladder table so a new one fails the suite. Done:
  the grep test passes and no compatibility layer is left.
  Prerequisites: GEN.348, UX.128. Related: GEN.300.

- [ ] **GEN.350 Class R, sub-Neptune, with the envelope-survival rule, the Hycean flag and the retirement of class T**
  Split from GEN.28 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work). Rock and ice core under a
  hydrogen-helium envelope, 1.6 to 3.2 Earth radii (cap 4), 3 to 10
  Earth masses, hot, ecosphere and cold; the envelope survives or is
  stripped by the energy-limited loss of the retention module; the
  Hycean window (1 to 10 Earth masses, 1.4 to 2.6 radii, T_eq 150 to 500
  K, density 1.5 to 3) is a flag on R; rocky rogues of 10 to 16 Earth
  masses become R; class T merges into R and is retired. Done: zone
  flags, weights, ranges, atmosphere, moon eligibility, rogue flag,
  reference page and the class-contract and 1,000-system tests; no old
  row or test expects T.
  Prerequisites: GEN.33, GEN.355, GEN.146. Related: GEN.28.

- [ ] **GEN.351 Class U, icy world**
  Split from GEN.28 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work). Water ice over rock, 500 to 3,000 km,
  cold (Ganymede, Callisto, Triton, Pluto, Eris). Data only. Done: zone
  flags, weights, ranges, atmosphere, moon eligibility, rogue flag,
  reference page and the class tests.
  Prerequisites: GEN.33. Related: GEN.28.

- [ ] **GEN.352 Classes W and X, subsurface ocean bodies**
  Split from GEN.28 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work). W: small, moon-only (not a rogue),
  about 50 to 500 km, the Enceladus analog. X: terrestrial sized, 500 to
  11,500 km, the Europa analog and larger, extended to take the
  cold-zone role GEN.91's option A gives it, using the shell and ocean
  depth rogue_surface.py already has. Done: both classes with zone
  flags, weights, ranges, atmosphere, moon eligibility (W moon-only),
  rogue flag (X), reference pages and the class tests.
  Prerequisites: GEN.33. Related: GEN.28.

- [ ] **GEN.353 Class Y, Titan-like world**
  Split from GEN.28 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work). Thick nitrogen atmosphere, methane
  rain, hydrocarbon lakes, 1,500 to 4,000 km, cold, allowed only where
  the shoreline ratio is under 0.7 and the vapour-pressure gate passes.
  Done: zone flags, weights, ranges, atmosphere, moon eligibility, rogue
  flag, reference page and the class tests.
  Prerequisites: GEN.355, GEN.146. Related: GEN.28.

- [ ] **GEN.354 Class Z, lifeless temperate world**
  Split from GEN.28 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work). Earth-sized rocky world with an
  atmosphere in the ecosphere and no life (Boss, 2026-10-03: an
  atmosphere that was never metabolised), with the M-dwarf and G/K
  sub-types of design section 5.1. It comes last because it needs the
  partial pressures and the pre-main-sequence rule of GEN.92. Done: zone
  flags, weights, ranges, atmosphere, moon eligibility, rogue flag,
  reference page and the class tests, and life stage none.
  Prerequisites: GEN.92, GEN.146. Related: GEN.28.

- [ ] **GEN.355 Atmosphere retention module with its test table**
  Split from GEN.91 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work).
  `planetgen/physics/atmosphere_retention.py`: the shoreline,
  Jeans-species, condensation and energy-limited functions of design
  section 3.6, XUV history and exposure from GEN.86's activity module
  (not the stand-in), no class change. It answers 'when is an atmosphere
  possible' in the docs: shoreline ratio below 30, species by lambda (30
  kept), the cold side by vapour pressure, H2 and water by
  energy-limited XUV. Done: unit tests on the 13-body table (Earth 0.17
  dense, Venus 0.43, Mars 1.77 tenuous, Mercury 53.7 none, Moon 81.6
  none, Titan 0.59) and the answer written in the design note.
  Prerequisites: GEN.33. Related: GEN.91.

- [ ] **GEN.356 Class N re-scoped as the hot, lifeless world; its life chemical and hidden timeline removed**
  Split from GEN.91 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work). N (the Venus analog) sits in the
  ecosphere at S about 1.05, below the runaway threshold, carries a life
  chemical and a hidden uncapped timeline. Re-scope it to the hot zone
  only, lifeless, up to 1.8 Earth radii, with no life data stored; a
  1-bar 'dune world' is a 10 percent variant of N. It changes existing
  output, so it uses GEN.92's 'unchanged until regenerated' guard. Done:
  N is generated only above the runaway threshold with no life chemical
  or timeline, the class page text matches, and existing galaxies are
  unchanged until regenerated.
  Prerequisites: GEN.355, GEN.92. Related: GEN.91.

- [ ] **GEN.357 Generation knows the cloud it sits in: surrounding_cloud extended, containers found early, StarSystem.apply_cloud**
  Split from GEN.95 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work). `surrounding_cloud` (set only at load
  today) gains extinction, radius, age, host stars and centre at
  generation; `run_sector` finds the containers of the sector box once
  (`store._placed_containers`, `NebulaShape.contains`); after the
  position is chosen `StarSystem.apply_cloud(cloud)` regenerates planets
  only when the cloud's rules change something and stamps the radiation
  tier, enrichment factor and composition flags. The same call serves
  GEN.99's regenerate in place. Done: a system inside a seeded nebula
  gets the cloud while it is generated, a system outside is
  byte-identical.
  Owner: Core 4, "Core 4 (2026-10-16)" (Boss approved the lane at 4:31
  PM EDT, 2026-10-16; queue in the lane report). READY: its blocker DB.59
  turned out to be built already (the nebula and remnant age columns have
  existed since v87 and v38; Core 2, 6:05 PM EDT), so Core 4 can take it
  now.
  Prerequisites: none. Related: GEN.95, DB.59.

- [ ] **GEN.358 GEN.94's rules and the four exotic-environment additions reach planet formation**
  Split from GEN.95 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work). Apply the rule table inside
  `apply_cloud`: the photoevaporation cut radius scales with host mass
  (200/50/10 AU is the solar-mass calibration); planetary nebulae H to L
  flag a circumbinary young disc only for a close-binary central star
  and never planets; an engulfed giant of over 5 Jupiter masses sets the
  'swallowed giant' flag; surviving planets near the engulfment radius
  get an eccentricity draw up to about 0.3; orbits expand by a (M_star +
  M_planet) = constant under adiabatic mass loss (impulsive loss by
  vis-viva, unbound above half the mass lost, engulfment overrides).
  Done: a test for each rule and the worked example (1 AU round 1 Msun
  dropping to 0.6 Msun moves to 1.67 AU).
  Prerequisites: GEN.357. Related: GEN.95.

- [ ] **GEN.359 Sector fill links cluster stars to their cluster and names the cluster at its first fill**
  Split from GEN.209 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). When a sector fill draws stars
  from a cluster's profile, each star gets the cluster link; the first
  fill that touches an unnamed cluster gives it one unique name from the
  naming codec, inside one transaction so two parallel workers cannot
  both name it; a rerun or a resumed fill finds the name already there.
  A catalogue cluster (GEN.297) keeps its catalogue name. Done: the link
  and the once-only name, with tests for a first fill, a second sector
  of the same cluster, a resumed fill and two workers filling at once.
  Prerequisites: DB.61. Related: GEN.209.

- [ ] **GEN.360 Asteroid fields and belts render from their uid: one seeded generator in a shared asteroidfield.js for dust, rocks and heroes**
  Split from GEN.112 (2026-10-16, Boss (2026-10-16, 3:00 PM EDT):
  "recommendations are otherwise approved" (the decision cards). The
  plan in docs/design/nebula-and-asteroid-field-classes.md, section
  "Asteroid fields and belts as rendered object systems", is approved;
  this is step 1 of its PR order). Today `systemview3d.js` seeds the
  belt from the database id, which changes when the database is rebuilt,
  and draws 700 points. Done: the render seed is `unit_seed(galaxy_seed,
  "asteroid-render", render_version/uid_hex)` computed on the server and
  passed in the view data; one `asteroidfield.js` generates dust, rocks
  and heroes from it (truncated power-law sizes, belt and field
  distributions, exaggerated sizes with a stated factor) and replaces
  `BELT_PARTICLES = 700`; the same uid gives the same bodies after a
  database rebuild. Tests: the first 1,000 sizes and positions hash the
  same for the same seed and differ for another uid; the sizes match the
  truncated Pareto (KS test); no rock overlaps; Kirkwood gaps stay
  empty.
  Prerequisites: none. Related: GEN.112, GEN.110.
  Design: [docs/design/nebula-and-asteroid-field-classes.md](design/nebula-and-asteroid-field-classes.md)

- [ ] **GEN.361 The asteroid field page has a 3D view**
  Split from GEN.112 (2026-10-16, Boss (2026-10-16, 3:00 PM EDT):
  "recommendations are otherwise approved" (the decision cards). The
  plan in docs/design/nebula-and-asteroid-field-classes.md, section
  "Asteroid fields and belts as rendered object systems", is approved;
  this is step 2 of its PR order). `phenomenonrender.py` lists
  `asteroid_field` in `_NO_VIEW`, so a field page has no 3D view. Done:
  the entry is removed and the field page shows its field through the
  shared generator of the first step, with the stated size factor and
  the real body count ("about 2e9 bodies larger than 1 km") on the page;
  a browser test opens a field page and checks the view draws.
  Prerequisites: GEN.360. Related: GEN.112.
  Design: [docs/design/nebula-and-asteroid-field-classes.md](design/nebula-and-asteroid-field-classes.md)

- [ ] **GEN.362 Rocks are instanced meshes in chunks with three levels of detail, within a frame budget**
  Split from GEN.112 (2026-10-16, Boss (2026-10-16, 3:00 PM EDT):
  "recommendations are otherwise approved" (the decision cards). The
  plan in docs/design/nebula-and-asteroid-field-classes.md, section
  "Asteroid fields and belts as rendered object systems", is approved;
  this is step 3 of its PR order). Done: rocks are per-chunk
  `InstancedMesh` objects (a 4x4x4 grid, so whole chunks cull) from 4 to
  6 noise-displaced variants at 80, 180 and 320 triangles, chosen by
  chunk distance; the dust, rock and hero layers appear and fade with
  the camera distance rules of the plan (sprite beyond 8 radii, dust
  inside 3, rocks inside 1, heroes inside 0.2); colours come from the
  composition family. A browser test measures 5,000 instances against a
  frame budget and reports the figure.
  Prerequisites: GEN.361. Related: GEN.112.
  Design: [docs/design/nebula-and-asteroid-field-classes.md](design/nebula-and-asteroid-field-classes.md)

- [ ] **GEN.363 Belt bodies move on the GPU: the vertex shader solves Kepler from per-instance elements**
  Split from GEN.112 (2026-10-16, Boss (2026-10-16, 3:00 PM EDT):
  "recommendations are otherwise approved" (the decision cards). The
  plan in docs/design/nebula-and-asteroid-field-classes.md, section
  "Asteroid fields and belts as rendered object systems", is approved;
  this is step 4 of its PR order). A belt is still one object for
  positional updates; only its drawing moves bodies. Done: the vertex
  shader solves Kepler from per-instance elements (a, e, i, Omega,
  omega, M0) and a time uniform, mean motion sqrt(mu/a^3); nothing is
  stored per body; the view's `update(years)` drives it; field bodies
  are static or spin rigidly and slowly from the seed; picking still
  selects the whole belt or field. A test compares the shader positions
  with the Python solver for a few bodies.
  Prerequisites: GEN.362. Related: GEN.112, GEN.139.
  Design: [docs/design/nebula-and-asteroid-field-classes.md](design/nebula-and-asteroid-field-classes.md)

- [ ] **GEN.365 Hot white dwarf and planetary-nebula central star radius disagrees with the radius that its luminosity and temperature imply (bug)**
  Found by a lane (2026-10-16), filed by the TODO thread. The radius
  stored for a hot white dwarf or the central star of a planetary nebula
  does not match the radius that L and T give (L = 4 pi R^2 sigma T^4).
  It was carried over from GEN.149. Done: the generator derives the
  radius from L and T for these stars or the luminosity from the radius,
  whichever is the model input, and a test checks the relation within 1%
  on a sample. Fails on main: yes. Frequency: every such star.
  Reproducer: draw 100 hot white dwarfs and compare R with the radius
  from L and T.
  Prerequisites: none. Related: GEN.149.

- [ ] **GEN.366 Check the reference satellite table against McConnachie 2012: it only tests that the derived population looks plausible**
  Found by a lane (2026-10-16), filed by the TODO thread. The table of
  Milky Way satellites (11 rows) was written from memory of McConnachie
  (2012) and is unverified. Done: every row (position, distance,
  luminosity, size) is checked against the published table, corrected
  where it differs, with the source named in a comment and a test that
  pins the checked values. GEN.237 (merged as #1851, 2026-10-16, 3:15 PM
  EDT) created the table; Core 1 reports that the Large and Small
  Magellanic Cloud, Sagittarius and eight dwarf galaxies (sizes, masses
  and distances) were recalled from memory and could not be checked.
  A galaxy planned before #1851 has no satellites until `planetgen plan`
  runs again on it.
  Folded into GEN.369 (Boss, 2026-10-16, 4:01 PM EDT, "Draw them"): the
  real rows are kept only as reference data for the plausibility test of
  GEN.369, with no Milky Way option (Boss, 2026-10-16, 4:11 PM EDT), so
  this check matters only for that test and waits on GEN.369.
  Prerequisites: GEN.369.

- [ ] **GEN.367 The Magellanic Clouds' gas is not in `kind_factors` and `nebula_field`: the Clouds hold stars but no nebulae**
  Found by a lane (2026-10-16), filed by the TODO thread. The Large and
  Small Magellanic Cloud components make stars but their gas is not in
  `kind_factors` or `nebula_field`, so the Clouds have no nebulae. Done:
  the Clouds get a gas model (H I and H II mass and the nebula kind mix)
  in the nebula field, with the numbers from the literature named, and a
  test that a Cloud sector holds nebulae at the stated rate. Later item;
  candidate for the deep-sky version.
  Prerequisites: GEN.237. Related: GEN.237, GEN.282.

- [ ] **GEN.368 Plan the appearance rates of everything the Milky Way supplies as fixed content: the count, place and size laws, with sources, so a random galaxy draws them from its seed**
  Boss (2026-10-16, 4:00 PM EDT, in the Core 1 thread; he tapped "Draw
  them" on the decision card at 4:01 PM EDT): asked whether the
  satellites are placed as in the Milky Way, said this is a random
  procedural galaxy so identical things make no sense, and wants the
  probabilities of things appearing worked out and planned. Done: one
  design section (in docs/design/intergalactic-and-halo-contents.md, a
  new part beside 5.5) that lists every place the generator copies the
  Milky Way as fixed content and replaces each with a stated draw: the
  satellite population (count follows the halo mass, distances a halo
  distribution, masses a satellite mass function), the globular cluster
  count and law, the galaxy shape (disk scale length and height, bulge
  share, arm count, halo mass), the central black hole mass and the
  nuclear cluster from their mass relations, the home position, and the
  rotation curve and enclosed mass. Each row gives its law, its
  parameters, the sources named, and the source it comes from; there is
  no Milky Way option (Boss, 2026-10-16, 4:11 PM EDT), so the real
  values appear only as reference data for plausibility tests.
  Boss reads and approves the table before any build. The fixed places
  found on 2026-10-16 (TODO thread survey): MILKY_WAY_SATELLITES (11
  copied rows, generation/satellites.py); the Harris clusters of GEN.297
  and CLUSTER_COUNT_MILKY_WAY = 157 (tuning.py); every shape default of
  `generate` (scale length 2,600 pc, height 300 pc, bulge 3.11, COBE bar
  1,580 pc); the central black hole (Sgr A*, 4.3e6) and nuclear cluster
  (GEN.280); GALAXY_SUN_RADIUS_PC = 8,200 and
  GALACTIC_CENTER_DISTANCE_LY = 25,800; MILKY_WAY_MASS (1.15e12) as the
  central mass of every Hill sphere (binary.py, star.py, keepout.py);
  the rotation curve (220 km/s, 3,000 pc; physics/constants.py). Owner:
  Core 1 (it built GEN.237 and knows the table), after the global top 5
  blockers unless Boss moves it up.
  Boss (2026-10-16, 4:04 PM EDT, in the Core 1 thread): "we don't need
  them to be the exact distances. What we need is for that to be
  calculated or derived or generated, preferably calculated and derived
  using real science." The table prefers a calculated or derived law
  over a drawn one for every row, with the real science named; a drawn
  value is allowed only where the science gives a distribution.
  Method (Boss, 2026-10-16, 4:04 PM EDT: derive, do not hand-set; the
  chain is Core 1's suggestion, each source to be checked in GEN.368):
  the galaxy's halo mass from its own mass model -> virial radius ->
  subhalo count from the subhalo mass function (dN/dM roughly M^-1.9,
  from the Aquarius and Via Lactea simulations) -> radial placement from
  a subhalo profile (less concentrated than the dark matter) -> stellar
  mass from the stellar-to-halo mass relation (abundance matching) with
  a low-mass cutoff for reionization -> size from the size-mass relation
  -> [Fe/H] from the mass-metallicity relation (Kirby 2013) -> age and
  profile from the mass and gas fraction, with tidal stripping near the
  disk. Each link names its source in the doc and the code. The real
  Milky Way rows serve only as a check that the derived population looks
  plausible.
  Boss (2026-10-16, 4:11 PM EDT, in the TODO thread): "There should be
  no option to pick a Milky Way preset. The Harris file is designed to
  give you an idea of how to build the procedure to generate a realistic
  galaxy." The Milky Way is a reference for building and checking the
  laws, never selectable content.
  Boss (2026-10-16, 4:14 PM EDT, in the TODO thread): "I want to
  reiterate that the idea is that we generate a scientifically accurate
  galaxy as best we can. We generate the sectors of that, again as
  accurately as we can, with what we can actually do processing-wise and
  what we actually know. Astrophysics-wise that can be implemented I
  don't want any fixed content. Everything should be procedurally
  generated based on the random seeds that we already use. Random
  numbers are used where we can. Things are derived everywhere we can or
  calculated everywhere we can. If needed you should ask me to do
  research if there are pieces missing that you can't look up yourself."
  This widens the survey from the Milky Way to all fixed content: any
  catalogue, table of named objects, preset, hand-placed object or
  hand-set value that stands for content rather than for a physical
  constant. Each row of the table is replaced by a derived or calculated
  law, with a seeded draw only where the science gives a distribution.
  The item also keeps a list of the pieces that cannot be looked up,
  each with what is needed and why, and hands it to Boss as research
  requests (the TODO thread puts them to him).
  Finding (read-only thread, 2026-10-16, 7:21 PM EDT, from Boss's 7:24
  PM EDT question "should the super massive black hole at the center of
  the galaxy always be at exactly 0,0,0?"; repo facts checked on main):
  the nucleus is pinned to (0,0,0) in three places: `store.GALACTIC_CENTER_PLACEMENT`
  / `insert_quasar` (snaps any placement to the origin),
  `phenomenon_scatter.nucleus_row` (position (0,0,0) in cell (0,0,0)) and
  `run_sector.add_galactic_nucleus` (position (-dist,0,0) in the sector
  frame); hypervelocity stars start from the origin too. Decided (Boss,
  2026-10-16, 7:27 PM EDT, card "Exactly 0,0,0"): the central black hole
  stays exactly on the origin, with no seeded offset or drift; the read-only
  thread's offset proposal is not built. The hole's mass is still drawn,
  not derived: `BLACK_HOLE_SUPERMASSIVE_MASS_RANGE_SOLAR` (1e6, 1e8)
  log-uniform (`tuning.py`), quasar 1e8 to 1e10 log-uniform, both
  independent of the bulge mass or the galaxy shape (real holes track
  bulge mass and dispersion, M-sigma, Haring and Rix; figures from memory,
  to be checked here), and the sphere-of-influence dispersion is a
  constant 100 km/s (`BLACK_HOLE_SUPERMASSIVE_VELOCITY_DISPERSION_KMS`,
  course-avoidance.md). Work for this item's central black hole row:
  derive the hole mass and the dispersion from the galaxy's bulge mass or
  dispersion by an M-bulge or M-sigma relation with its scatter as a
  seeded draw, instead of the independent log-uniform range; the active
  and quiescent split stays a draw. Kind: design gap, not a bug.
  Prerequisites: none. Related: GEN.237, GEN.369, GEN.297, GEN.370.
  Design: [docs/design/intergalactic-and-halo-contents.md](design/intergalactic-and-halo-contents.md)

- [ ] **GEN.369 Satellites drawn from the galaxy seed: count, distance, mass, size, profile, [Fe/H] and age are calculated from stated laws, with no Milky Way option**
  Boss (2026-10-16, 4:00 PM EDT, in the Core 1 thread; he tapped "Draw
  them" on the decision card at 4:01 PM EDT): asked whether the
  satellites are placed as in the Milky Way, said this is a random
  procedural galaxy so identical things make no sense, and wants the
  probabilities of things appearing worked out and planned. Decision:
  the 11 copied Milky Way satellites of GEN.237 (`MILKY_WAY_SATELLITES`,
  generation/satellites.py) are replaced by a drawn population from the
  galaxy seed: the count follows the halo mass, distances a halo
  distribution, masses a satellite mass function, size, profile
  (exponential for gas-rich, Plummer for dwarfs), [Fe/H] and age follow
  mass relations. Supersedes the "Real Milky Way satellites at their
  real places" open default in
  docs/design/intergalactic-and-halo-contents.md section 8, item 2.
  Done: `place_satellites` draws the population from
  `unit_seed(galaxy_seed, ...)` through the laws GEN.368 sets, the real
  rows are removed from the generator and kept only as reference data
  for the plausibility test of GEN.366 (Boss, 2026-10-16, 4:11 PM EDT:
  no option to pick a Milky Way preset), a test
  shows two seeds give different populations and one seed gives the same
  one, and the count scales with the halo mass. Owner: Core 1 (the code
  owner), after the global top 5 blockers unless Boss moves it up; the
  component pass GEN.324 places stars for whatever component rows exist,
  so it does not wait on this.
  Prerequisites: GEN.368, GEN.371. Related: GEN.237, GEN.366, GEN.324,
  GEN.297.
  Waits on GEN.371, the feasibility study of generating satellites
  through the scatter (Boss, 2026-10-16, 4:02 PM EDT); if the study
  finds it feasible, the drawn satellites are made by the scatter in the
  cells they occupy.
  Boss (2026-10-16, 4:04 PM EDT, in the Core 1 thread): "we don't need
  them to be the exact distances. What we need is for that to be
  calculated or derived or generated, preferably calculated and derived
  using real science." Acceptance: the satellites' count, distances and
  sizes are calculated from real science wherever a published relation
  exists (the subhalo mass function and its radial distribution for
  count and distance, the mass-size and mass-metallicity relations for
  size and [Fe/H]), drawn only where the science gives a distribution
  and not a value, and never copied or hand-picked; each law names its
  source in the code and in the GEN.368 table.
  Method (Boss, 2026-10-16, 4:04 PM EDT: derive, do not hand-set; the
  chain is Core 1's suggestion, each source to be checked in GEN.368):
  the galaxy's halo mass from its own mass model -> virial radius ->
  subhalo count from the subhalo mass function (dN/dM roughly M^-1.9,
  from the Aquarius and Via Lactea simulations) -> radial placement from
  a subhalo profile (less concentrated than the dark matter) -> stellar
  mass from the stellar-to-halo mass relation (abundance matching) with
  a low-mass cutoff for reionization -> size from the size-mass relation
  -> [Fe/H] from the mass-metallicity relation (Kirby 2013) -> age and
  profile from the mass and gas fraction, with tidal stripping near the
  disk. Each link names its source in the doc and the code. The real
  Milky Way rows serve only as a check that the derived population looks
  plausible.
  Boss (2026-10-16, 4:14 PM EDT, in the TODO thread): "I want to
  reiterate that the idea is that we generate a scientifically accurate
  galaxy as best we can. We generate the sectors of that, again as
  accurately as we can, with what we can actually do processing-wise and
  what we actually know. Astrophysics-wise that can be implemented I
  don't want any fixed content. Everything should be procedurally
  generated based on the random seeds that we already use. Random
  numbers are used where we can. Things are derived everywhere we can or
  calculated everywhere we can. If needed you should ask me to do
  research if there are pieces missing that you can't look up yourself."
  Design: [docs/design/intergalactic-and-halo-contents.md](design/intergalactic-and-halo-contents.md)

- [ ] **GEN.370 Hill spheres, keep-out spheres and the default tidal distance use the Milky Way's mass and Sol's distance for every galaxy: take the galaxy's own enclosed mass and the cell's own distance (bug)**
  Found by the TODO thread (2026-10-16, survey for Boss's 4:00 PM EDT
  request): `MILKY_WAY_MASS` (1.15e12 solar masses) is the central mass
  of every Hill sphere (`generation/binary.py`, `generation/star.py`,
  `galaxy/keepout.py`) and `GALACTIC_CENTER_DISTANCE_LY` (25,800) is the
  default distance when none is given, so a galaxy of another size or
  mass gets Milky Way tides, and a point-mass central mass overstates
  the enclosed mass inside the disk (the orbit function already avoids
  that, `galaxy/galactic_orbit.py`). Done: the central mass is the
  galaxy's mass enclosed within the star's radius from the density
  model, the default distance is removed in favour of the position the
  caller has, and a test shows a galaxy twice the size changes the Hill
  radius. Fails on main: no failure, a fixed assumption. Related:
  GEN.368, GEN.285.
  Prerequisites: none. Related: GEN.368, GEN.285, GEN.223.

- [ ] **GEN.371 Feasibility study: generate the satellite galaxies through the scatter, making only the region each one occupies**
  Boss (2026-10-16, 4:02 PM EDT, in the Core 1 thread): "we'll need a
  feasibility study. I would like to generate them in the scatter if
  possible, even if we just generate the region that they are in." Done:
  a written study (a section beside 5.5 of
  docs/design/intergalactic-and-halo-contents.md) that answers whether a
  drawn satellite can be made by the same scatters as the rest of the
  sphere, so its stars, gas and clusters come from the one rule set and
  a profile, and what it costs when only the cells the satellite
  occupies are generated. It covers the profile a satellite cell needs
  (its density law, [Fe/H], age and star-formation history), how the
  cell key and the 96-bit ID address a satellite's cells, how many cells
  and systems a dwarf, the Clouds and a Sagittarius-sized stream make at
  the default cell size, the run time per satellite by the measured cost
  per cell, how the component pass GEN.324 and the cluster chain feed
  it, and a recommendation with the work it needs as new items. No build
  in this item. Owner: Core 1, with GEN.368, after the global top 5
  blockers unless Boss moves it up.
  Boss (2026-10-16, 4:04 PM EDT, in the Core 1 thread): "we don't need
  them to be the exact distances. What we need is for that to be
  calculated or derived or generated, preferably calculated and derived
  using real science." The study also reports which of a satellite's
  properties can be calculated rather than drawn from the halo mass and
  the galaxy's history, and what the scatter needs to make its region
  from those.
  The study follows the same chain: the galaxy's halo mass from its own
  mass model -> virial radius -> subhalo count from the subhalo mass
  function (dN/dM roughly M^-1.9, from the Aquarius and Via Lactea
  simulations) -> radial placement from a subhalo profile (less
  concentrated than the dark matter) -> stellar mass from the
  stellar-to-halo mass relation (abundance matching) with a low-mass
  cutoff for reionization -> size from the size-mass relation -> [Fe/H]
  from the mass-metallicity relation (Kirby 2013) -> age and profile
  from the mass and gas fraction, with tidal stripping near the disk.
  Each link names its source in the doc and the code. The real Milky Way
  rows serve only as a check that the derived population looks
  plausible.
  Boss (2026-10-16, 4:14 PM EDT, in the TODO thread): "I want to
  reiterate that the idea is that we generate a scientifically accurate
  galaxy as best we can. We generate the sectors of that, again as
  accurately as we can, with what we can actually do processing-wise and
  what we actually know. Astrophysics-wise that can be implemented I
  don't want any fixed content. Everything should be procedurally
  generated based on the random seeds that we already use. Random
  numbers are used where we can. Things are derived everywhere we can or
  calculated everywhere we can. If needed you should ask me to do
  research if there are pieces missing that you can't look up yourself."
  Prerequisites: GEN.368. Related: GEN.368, GEN.369, GEN.324, GEN.237.
  Design: [docs/design/intergalactic-and-halo-contents.md](design/intergalactic-and-halo-contents.md)

- [ ] **GEN.372 Numeric precision audit at galactic scale: every float coordinate and sum, how far it gets from its origin, where the first error shows, and what needs a wider or split representation**
  Boss (2026-10-16, 6:54 PM and 6:55 PM EDT, project chat): are 64-bit
  numbers enough for positions at galactic scale, or do we need 128 or
  256 bits, and who should own all of our bugs and concerns about it.
  This is the one place they are collected. Existing coverage is
  scattered: docs/design/orbital-updates.md (positions as integer
  milliparsecs relative to the sector and kilometres relative to a
  system, and a float64 spacing line of 32.8 km at 8.128 kpc),
  docs/design/orbital-solvers-and-integrators.md (a float64 spacing
  table), docs/design/galactic-potential.md (cancellation and float32
  notes), the spherical sector design (a coarse integer address with a
  local offset), the field kernel of PERF.88 (float64 sums), TEST.192 (a
  cancellation in a test), and a GEN.142 note that says "about 14 km at
  8 kpc". These disagree: 2^-40 pc is about 28 km at 8.128 kpc, and at
  the 262 kpc reach the spacing is about 900 km below 262,144 pc and
  about 1,800 km above it. Done: a written audit
  (docs/design/numeric-precision-audit.md) with one table row for each
  float quantity the galaxy stores, computes or draws (the galactic
  position of every object kind, the sector-relative offset, system,
  planet and moon frames, velocities, epochs and phases, orbital
  elements, the field kernel and potential sums, the camera and render
  path including float32 in the browser, the JavaScript twins, the API
  wire formats, the 96-bit ID and cell key arithmetic), giving the
  largest distance from its origin it reaches (up to the 262 kpc reach),
  the spacing of a float64 there in kilometres, the smallest change that
  must survive (a moon's move, the 0.01 mpc threshold, a planet's
  orbit), and the first place the error shows. It says for each row
  whether float64 plus the coarse integer address and a local offset is
  enough, or where a wider sum is needed (a compensated or integer sum,
  or mpmath, for long integrations from the galactic centre), recommends
  the cases that need work as new items, and names a test for each (for
  example a differenced-position test at 262 kpc). No build in this
  item. Findings rule: any lane that meets a precision bug or concern
  sends it to the TODO thread, which attaches it here. Owner: Core 1 (Boss
  chose it by card, 2026-10-16, 6:59 PM EDT), who owns the point columns
  (DB.55, GEN.336, GEN.337, GEN.294) and the IDs, with Core 2 reading the
  camera and render rows.
  Audit written (Core 1, 2026-10-16, 7:39 PM EDT; docs/design/numeric-precision-audit.md,
  merging in its own PR): verdict, 64-bit is enough for every stored
  value, no 128-bit or 256-bit column. Boss decided four items by card
  (7:36 to 7:38 PM EDT): the nightly update keeps a local offset (document
  and test only, TEST.200); the map tile root rises to the 262 kpc reach
  (MAP.239); browser buffers become camera-relative now (MAP.240);
  compensated sums wait for GEN.109. Items filed from the audit's section
  5: MAP.239, MAP.240, GEN.374 (phase formula), DB.65 (cell key clamp),
  GEN.375 (field kernel local origin), GEN.376 (region haversine),
  TEST.200 (local offset limit), TEST.201 (block_stats order). GEN.372
  resolves when the audit merges.
  Finding attached (Core 2, 2026-10-16, 7:27 PM EDT): TEST.199, the
  halo subtraction in the density skeleton loses precision in the far tail
  (a relative gap of 2.4e-9 against 1e-9 at z = 36,706 pc).
  Prerequisites: none. Related: GEN.294, GEN.139, GEN.142, GEN.223,
  MAP.146, NAV.62, TEST.192, PERF.88, DB.55, TEST.199.
  Design: [docs/design/orbital-updates.md](design/orbital-updates.md)

- [ ] **GEN.373 Massive multiples sit in young regions: set `prevalence` from the cell's young share, with a derived cap on the rise**
  Split from GEN.103 (Core 4, 2026-10-16, 7:16 PM EDT; the coordinator
  asked it be an item of its own at 7:17 PM EDT).
  docs/design/placement-by-population.md section 3.2: high-mass stars
  are far more often in multiples, and the most massive sit in clusters
  and associations. The mass dependence is already in
  `BINARY_SYSTEM_PROBABILITY_BY_SPECTRAL_CLASS` and
  `STAR_COUNT_PROBABILITIES`; the place dependence is not. Done:
  `prevalence` of `draw_member_tree` is the cell's young share (the
  profile's `populations["young"].share`) over the solar circle's
  reference share (`profile._reference_young_share(shape)`), so it is 1
  where a cell looks like the solar circle, above 1 in arms and
  associations and 0 in the halo; the five-or-more columns scale with
  it. The cap on the rise is the one open constant (a design default,
  not a value from a source): under Boss's no-fixed-content rule it is
  derived, not typed in; the odds form Core 4 recommended for GEN.103
  (the factor applied to the odds, so nothing passes 1) is the model.
  Boss decides the cap on a card. Test: the massive-multiple census bins
  of TEST.164, and the metallicity check of TEST.196. Owner: Core 4,
  with the multi-star items (GEN.269 to GEN.275).
  Prerequisites: none. Related: GEN.103, GEN.269, GEN.270, TEST.164,
  TEST.196.
  Design: [docs/design/placement-by-population.md](design/placement-by-population.md)

- [ ] **GEN.374 Orbit phase formulas reduce first, and the JavaScript Kepler solver uses tolerance 1e-12: Python, SQL and JavaScript agree**
  From the GEN.372 precision audit (Core 1,
  docs/design/numeric-precision-audit.md section 5, 2026-10-16, 7:39 PM
  EDT), filed by the TODO thread; verdict: 64-bit is enough for every
  stored value. Done: the phase formula reduces the angle modulo one
  turn before multiplying (`body_positions` and `system_tree` in Python,
  `src/html/static/orbitpositions.js` in JavaScript); the JavaScript
  Kepler solver's tolerance is 1e-12; tests show that the Python, SQL
  and JavaScript phase functions agree within a stated bound at long
  times (millions of days from the epoch). Owner: Core 2, or whoever
  owns `orbitpositions.js`.
  Prerequisites: none. Related: GEN.372, GEN.109, GEN.313.

- [ ] **GEN.375 Field kernel subtracts a local origin before its loop in `physics/sector_path.py`**
  From the GEN.372 precision audit (Core 1,
  docs/design/numeric-precision-audit.md section 5, 2026-10-16, 7:39 PM
  EDT), filed by the TODO thread; verdict: 64-bit is enough for every
  stored value. Done: the field kernel of
  `src/planetgen/physics/sector_path.py` subtracts a local origin (the
  path's own start or the cell centre) from every position before the
  loop, so the sums run on small numbers; a test shows the result agrees
  with a high-precision reference within the float64 spacing at a sector
  262 kpc from the origin. Owner: Core 1, with GEN.294.
  Prerequisites: none. Related: GEN.372, GEN.294.

- [ ] **GEN.376 Region test uses the haversine form for angles under 1e-3 rad in `galaxy/region.py`**
  From the GEN.372 precision audit (Core 1,
  docs/design/numeric-precision-audit.md section 5, 2026-10-16, 7:39 PM
  EDT), filed by the TODO thread; verdict: 64-bit is enough for every
  stored value. Done: the angle in the region test of
  `src/planetgen/galaxy/region.py` is computed with the haversine form
  below 1e-3 rad (the arccos form loses all precision there); a test
  with two points a few parsecs apart at 262 kpc gives the right angle
  within a stated relative bound. Owner: Core 4 (region query).
  Prerequisites: none. Related: GEN.372, GEN.276, NAV.60.

- [ ] **GEN.377 Scatter limits by place for the neutron star and black hole scatter below the cut: the core and halo tiers apply to it, not one global mass cut**
  Found by Core 3 (2026-10-16, 8:34 PM EDT, after GEN.333 merged as
  #1904), filed by the TODO thread at its request as a later 8.1 item.
  GEN.333 applies the tiers (core, baseline, halo) to the star scatter;
  the phenomena scatter below the cut (neutron stars and black holes
  from stellar density) still keeps ONE global mass cut, so the core and
  halo tiers do not apply to it. Done: the neutron star and black hole
  scatter below the cut takes its mass limit by place (core, baseline,
  halo) the same way GEN.333 does for stars, with the tier chosen per
  band by the dominant lambda term; the three limits are the same
  settings the Generate page already shows (UX.116), with no extra
  setting; a test shows that a core band and a halo band cut the
  phenomena at their own tier limits and the baseline is unchanged.
  Changed behaviour recorded with GEN.333: the brightest white dwarf is
  now 50 solar luminosities (it was 100), so the luminosity floors and
  caps may go as low as 50. Owner: Core 3.
  Prerequisites: none. Related: GEN.333, GEN.283, UX.116, GEN.301.

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
  Folded into the backfill-as-queue-units child of PERF.98 (chunks of at
  least 10,000 cells, one enumeration, pre-warmed tables); it closes
  when that child ships (2026-10-16 review).

- [ ] **PERF.20 Short-term caching through the work queue and API**
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
  Split (2026-10-16, Boss 1:07 PM EDT: look at each item for smaller
  chunks and shared work): this item is now the umbrella and is done
  when its subitems are: PERF.109, PERF.110. Its plan is finished
  (design performance-eta-queue-and-caching.md 4.4 lists four build
  items) but the item still reads 'needs planning', and two of the four
  are already built: (a) the epoch counter is PERF.104's `epoch` tag
  (docs/design/cache-warm-up.md, status built) and (b) single-flight for
  tiles and stages is `web/lib/singleflight.py`, used in
  `web/lib/tilecache.py`.
  Prerequisites: PERF.109, PERF.110.

- [ ] **PERF.35 An interval or chunk ledger for untouched sectors once block-first backfill lands**
  Replace the one-`sector_stats`-row-per-visited-cell ledger of
  untouched sectors with an interval or chunk ledger once GEN.42's
  block-first draw lands; at 24 billion cells the ledger is 24 billion
  rows.
  Prerequisites: GEN.309.
  Design: [docs/design/sampling-backfill-and-resume.md](design/sampling-backfill-and-resume.md)

- [ ] **PERF.36 Memory and request guard: list at most about 50,000 candidate cells, stream the rest, and refuse enumerations above about 2 million cells in a web request**
  Stream the enumeration instead of materialising it; refuse to
  enumerate more than about 2 million cells in a web request. Add
  sampling-based estimation above about 20,000 sectors to the PERF.3
  estimate (200 samples land within 1% of the full 7.455 million-system
  core); the span and core modes need it.
  Decided (Boss, 2026-10-16, 3:00 PM EDT, decision card): both caps
  stand, the title's 50,000 for a listed set of candidate cells and the
  text's 2 million for what a web request may enumerate; they guard two
  different things.
  Prerequisites: none.
  Design: [docs/design/fill-order-curves-and-core.md](design/fill-order-curves-and-core.md)

- [ ] **PERF.41 Stamp the tile and page caches with a TILE_FORMAT constant instead of the version (optional)**
  A `TILE_FORMAT` constant in the stamp base in place of `__version__`,
  with a golden test that fails when the tile or page structure changes
  without a bump. Releases are 26 to 123 a day, so a development server
  never keeps a cache. Decided (Boss, 2026-10-10 18:44Z, defaults
  approved; default): decide after
  PERF.31 shows the cold rebuild cost; switch if a full rebuild takes
  more than a minute.
  Design: [docs/design/performance-eta-queue-and-caching.md](design/performance-eta-queue-and-caching.md)

- [ ] **PERF.47 The PERF.31 benchmark records the buffer pool, table sizes and worker start-up cost**
  Research (2026-10-09, generation-performance-study.md, PR #835;
  handoff in
  /mnt/project-files/research/handoff/generation-performance.md; from
  Boss's requests of 19:08Z and 19:21Z, generation being his slowest
  point): record `innodb_buffer_pool_size`, the table sizes and the
  worker start-up cost with every benchmark run, so a result can be
  compared with the next.
  Design: [docs/design/generation-performance-study.md](design/generation-performance-study.md)

- [ ] **PERF.48 Low priority: a numeric-only INSERT formatter or C driver for bright_stars and phenomenon_scatter**
  Research (2026-10-09, generation-performance-study.md, PR #835;
  handoff in
  /mnt/project-files/research/handoff/generation-performance.md; from
  Boss's requests of 19:08Z and 19:21Z, generation being his slowest
  point): 13.8 down to 10.1 microseconds a row. Low priority.
  Prerequisites: none.
  Design: [docs/design/generation-performance-study.md](design/generation-performance-study.md)

- [ ] **PERF.72 Research the cost of a sector's gravity grid and where to cut between exact and aggregated sources**
  Done: a note (docs/design/gravity-map.md, section 'How it is
  computed', extended with measurements) says what a grid of 16 zones
  per edge costs for a dense sector (core, bulge and arm) and a sparse
  one, with numpy and the per-region aggregates of MAP.151, picks the
  near/far distance and the zone count, gives the error of the aggregate
  against the exact sum, and says when a grid is recomputed after an
  orbital update (every step, or only for sectors whose sources moved by
  a set amount). Ends with the build items filed.
  Math Research (2026-10-14): the field is a sum over point masses at many
  grid points, the same shape as the path integrator's acceleration; cost it
  as a numba or NumPy kernel from the start (not measured in the note).
  Prerequisites: GEN.198. Related: GEN.198, GEN.199, MAP.151
  (done, PR #1539).
  Scope after the review (2026-10-16 review): the cost study of the grid
  on top of the kernel (measure the kernel and the aggregate split, pick
  the near and far distance and the zone count, state the error); no
  kernel code of its own.
  Design: [docs/design/gravity-map.md](design/gravity-map.md)

- [ ] **PERF.75 Keyset paging for the other data tables and the jump-by-value box (Systems name sort done)**
  From the PERF.71 research (docs/design/slow-reads-and-timeouts.md, PR
  #1175; Boss's ask of 2026-10-10 20:50Z; measured on 2,000,000
  synthetic systems, MariaDB 10.11): `OFFSET` makes page n cost n times
  as much (rows 0, 20,000, 500,000 and 1,500,000 take 0.004, 0.13, 5.8
  and 12.1 s idle); a keyset page costs 2 ms at row 1,500,000 (8.5 s
  before). Done: the table routes and the API list routes take
  `after=<sort key, id>` next to `offset=`, `static/datatable.js` pages
  forward by key and jumps by value (a letter or a sector; the scroll
  bar position becomes an estimate), the condition is written with `OR`
  and never as a row constructor (it took 5.0 s on MariaDB 10.11, no
  index range), and a test checks the plan on MariaDB and MySQL. Decided
  by default (the three research defaults, 2026-10-10 21:41Z; stands
  unless Boss objects): jump-by-value is accepted in place of
  jump-by-row-number on the big tables. Built (PR #1199): keyset paging
  for the Systems name sort only. Still open: the other sorts, tables and
  API list routes, and the jump-by-value box (a letter or a sector).
  Owner: Foundations lane 1.
  Prerequisites: none. Related: PERF.64, PERF.70, PERF.71.
  Design: [docs/design/slow-reads-and-timeouts.md](design/slow-reads-and-timeouts.md)

- [ ] **PERF.86 Redo the work queue: drop the task and subtask logic and use Redis simply**
  GitHub issue [#1229](https://github.com/dwhagar/planetGen/issues/1229)
  (Boss, 2026-10-13 02:41Z): "I was mistaken about how Redis worked and
  so I made our queuing system way too complicated. I need to look at
  what we're doing and junk all of the logic on what makes a task and
  subtask and all of that. It's going to be a lot of work too but I
  think I'm wasting a lot of efficiency." Done: a design note on what
  the queue does today, what it needs to do, and the simplest design on plain Redis queues,
  then the rebuild; PERF.40, PERF.78, ADM.15 (all folded into PERF.95,
  done, PR #1506), PERF.29 and PERF.30 waited for this plan and are
  folded in or dropped. Decided (Boss, 2026-10-14
  13:24Z): keep a "cancel now" that stops a job, progress bars and ETAs
  per stage, per-stage timing records (Boss wants stats on how long
  everything takes) and the job log with its download link; an admin
  starts or cancels a run, so there is no automatic resume after an
  interruption; several workers on one job are fine. CHANGED from the
  2026-10-10 decision: the scatter no longer stays one job. Each scatter
  becomes its own item in the queue, run one after another, and the same
  applies to the phases of anything. Boss asks the research to cover the
  best ways to use Redis in a situation like ours and to give him
  options before the rebuild starts. Owner: a research thread first
  (coordinator to assign), then Foundations lane 1.
  Prerequisites: none. Related: PERF.40 (done, PR #1506), PERF.78 (done, PR #1506),
  PERF.29, PERF.30, ADM.15 (done, PR #1506), PERF.24.
  Redis review (2026-10-14): research done (Research Lane 1, PR #1350,
  docs/design/redis-and-the-work-queue.md, with
  /mnt/project-files/research/redis-review/). Boss adopted its Option A
  and approved every default in it on 2026-10-14 21:08Z (see the
  Decided line below). Filed from the same review, independent of the
  rebuild: PERF.93 (bug, done, PR #1364), PERF.94 (done, PR #1383), OPS.48 (done, PR #1395), ADM.54 (done, PR #1401).
  Decided (Boss, 2026-10-14 21:08Z): Option A of
  the Redis review, every default approved. Build steps filed as PERF.95
  (queue spine), PERF.96 (run record and ETA library, done, PR #1553), PERF.97 (mass
  scatter end to end), PERF.98 (remaining stages), PERF.99 (pages and
  cancel), PERF.100 (delete old machinery, migration). PERF.86 is done
  when they are; it follows the dead-worker bug PERF.93 (done) and the client
  module PERF.94 (done) in Foundations lane 1. DOC.20 rewrites the notes after
  step 6.

- [ ] **PERF.91 Bulk random draws for the scatter layers (4 to 12 times on batches above 1,000 draws)**
  Math Research (2026-10-14;
  docs/design/math-and-library-efficiency.md). State is handed between
  `random.Random` and NumPy's MT19937 so draws and generator state are
  identical: no galaxy changes and no reseed. Only the default-scale
  scatter reaches batches this large (0.3 s of 16 s at benchmark scale),
  so this is low priority until the default-scale scatter is the
  bottleneck. Done: the scatter layers draw in bulk, a test proves same
  seed, same draws, same final state. Owner: Foundations lane 1, last of
  the math items.
  Prerequisites: none. Related: GEN.56, PERF.58.
  Design: [docs/design/math-and-library-efficiency.md](design/math-and-library-efficiency.md)

- [ ] **PERF.97 Work queue step 3: the mass scatter end to end on the new queue**
  Part of PERF.86. Decided by Boss (2026-10-14 21:08Z): Option A of
  docs/design/redis-and-the-work-queue.md (RQ and Redis stay; Redis is
  the only place a run lives; each scatter its own queue item run in
  order; counters plus an event stream for progress; cancel in about 0.5
  s; the MySQL job tree retired); Redis required for every run including
  one worker; workers on demand with an idle timeout and no always-on
  service; one Redis instance; no Redis password while it binds
  loopback; stage timing into the existing MySQL stage tables from the
  completing unit; the event stream is the live log, written to the jobs
  directory at run end with the last 20 kept. Each step leaves main
  working, one PR, Foundations lane 1. Section 5.1 and 5.11 step 3. One
  stage on the new spine, the mass scatter: static chunks (units per
  stage capped near 4,000), JSON arguments, the done-set barrier (SCARD
  equals the unit count, SET NX on an advanced flag so one unit
  advances), the run's seed independent of worker count, planetgen
  galaxy following the events, and --workers 1 using the same path. Each
  scatter layer set is its own queue item run in order (Boss,
  2026-10-14). Done when: a galaxy fill whose mass scatter runs this way
  gives the same galaxy for the same seed at any worker count and the
  old path still runs the other stages.
  Owner: Core 2 (Boss, 2026-10-16, 3:00 PM EDT: the old Core queue runs
  after the new top 5 work, routed by area): PERF.97 right after the 8.1
  phase 1 top 5, PERF.98 to PERF.100 after the infinite zoom tree.
  Prerequisite: none (PERF.96 done, PR #1553). Related: PERF.86.

- [ ] **PERF.98 Work queue step 4: the remaining stages and the sector and backfill units on the new queue**
  Part of PERF.86. Decided by Boss (2026-10-14 21:08Z): Option A of
  docs/design/redis-and-the-work-queue.md (RQ and Redis stay; Redis is
  the only place a run lives; each scatter its own queue item run in
  order; counters plus an event stream for progress; cancel in about 0.5
  s; the MySQL job tree retired); Redis required for every run including
  one worker; workers on demand with an idle timeout and no always-on
  service; one Redis instance; no Redis password while it binds
  loopback; stage timing into the existing MySQL stage tables from the
  completing unit; the event stream is the live log, written to the jobs
  directory at run end with the last 20 kept. Each step leaves main
  working, one PR, Foundations lane 1. Section 5.1 and 5.11 step 4. Move
  every other stage (skeleton, phenomena, luminosity, sectors, link,
  backfill, population, settle) and the sector and backfill units onto
  run records, units and barriers; one-shot stages are one unit. Stage
  timing goes into generation_stage_runs and generation_stats from the
  unit that completes the stage; a hard-killed worker is seen at once by
  the work-horse-killed handler, and every job gets a real job_timeout.
  Done when: every generation command runs on the new queue and the old
  executor is unused.
  Split (2026-10-16, Boss 1:07 PM EDT: look at each item for smaller
  chunks and shared work): this item is now the umbrella and is done
  when its subitems are: PERF.106, PERF.107, PERF.108. Step 4 moves
  eight stages, the sector units and the backfill units onto the new
  queue, plus stage timing, the killed-worker handler and job timeouts:
  it touches generation/run_plan.py (1,594 lines), run_galaxy.py
  (1,664), run_sector.py (884) and bright_stars.py, and every later
  queue item (PERF.99, PERF.100, GEN.231, DB.29's gate) waits for all of
  it.
  Owner: Core 2 (Boss, 2026-10-16, 3:00 PM EDT: the old Core queue runs
  after the new top 5 work, routed by area): PERF.97 right after the 8.1
  phase 1 top 5, PERF.98 to PERF.100 after the infinite zoom tree.
  Prerequisites: PERF.106, PERF.107, PERF.108. Related: PERF.86.

- [ ] **PERF.99 Work queue step 5: the Generate page, Queue page and job log read Redis; cancel now through the stop command**
  Part of PERF.86. Decided by Boss (2026-10-14 21:08Z): Option A of
  docs/design/redis-and-the-work-queue.md (RQ and Redis stay; Redis is
  the only place a run lives; each scatter its own queue item run in
  order; counters plus an event stream for progress; cancel in about 0.5
  s; the MySQL job tree retired); Redis required for every run including
  one worker; workers on demand with an idle timeout and no always-on
  service; one Redis instance; no Redis password while it binds
  loopback; stage timing into the existing MySQL stage tables from the
  completing unit; the event stream is the live log, written to the jobs
  directory at run end with the last 20 kept. Each step leaves main
  working, one PR, Foundations lane 1. Section 5.5 and 5.11 step 5. The
  Generate page, the Queue page and the job log read the run record, counters and event stream; cancel sets
  the run to cancelling, cancels the waiting units, sends the stop
  command to the started ones and clears planetgen:run:active (about 0.5
  s); the stale-run Clear action (section 5.6) is added; the event
  stream is the live log and is written to the jobs directory at run end
  for the download link, last 20 kept. Done when: the pages need no
  MySQL job tree and cancel is measured under a second on a 75-job run.
  The all-visitors banner is not in this step: it is UX.3 (UX),
  built on the run record after PERF.96 (decided by the
  coordinator, 2026-10-15).
  Owner: Core 2 (Boss, 2026-10-16, 3:00 PM EDT: the old Core queue runs
  after the new top 5 work, routed by area): PERF.97 right after the 8.1
  phase 1 top 5, PERF.98 to PERF.100 after the infinite zoom tree.
  Prerequisite: PERF.98. Related: PERF.86, UX.3.

- [ ] **PERF.100 Work queue step 6: delete the old machinery and the work_* tables (Alembic migration)**
  Part of PERF.86. Decided by Boss (2026-10-14 21:08Z): Option A of
  docs/design/redis-and-the-work-queue.md (RQ and Redis stay; Redis is
  the only place a run lives; each scatter its own queue item run in
  order; counters plus an event stream for progress; cancel in about 0.5
  s; the MySQL job tree retired); Redis required for every run including
  one worker; workers on demand with an idle timeout and no always-on
  service; one Redis instance; no Redis password while it binds
  loopback; stage timing into the existing MySQL stage tables from the
  completing unit; the event stream is the live log, written to the jobs
  directory at run end with the last 20 kept. Each step leaves main
  working, one PR, Foundations lane 1. Section 5.3, 5.7 and 5.11 step 6.
  Delete queue/work.py's WorkQueue, the job tree and the lease,
  RQExecutor, RQFuture and Channel, web/jobs.py's files and
  web/job_runner.py, queue/progress_file.py, submit_each's adaptive
  chunking, the in-process no-Redis fallback and GENERATION_QUEUE, with
  one Alembic revision dropping work_jobs, work_job_args, work_tasks and
  work_lease (Foundations lane 1 only adds migrations). PERF.29 and
  PERF.30 (resume) are dropped when this lands: Boss decided there is no
  automatic resume. Done when: the suite is green with the old code gone
  and the line count is recorded in the PR.
  Owner: Core 2 (Boss, 2026-10-16, 3:00 PM EDT: the old Core queue runs
  after the new top 5 work, routed by area): PERF.97 right after the 8.1
  phase 1 top 5, PERF.98 to PERF.100 after the infinite zoom tree.
  Prerequisites: DB.29. Related: PERF.86, PERF.29, PERF.30, DOC.20.

- [ ] **PERF.105 Slow page requests queue a job that builds the answer into Redis, with a status route (the mechanism)**
  Boss (2026-10-16, 12:21 PM EDT, project chat): "when a page requests
  something that's going to take a while (say listing all the black
  holes in the phenomena page) it should add a task to get that
  information to redis, have a progress bar in the search page so the
  user knows it is loading, that is one of the main issues, give it to
  Core 2's chain". Owner Core 2, in its queue after MAP.152 and MAP.146 (Boss, 12:23 PM EDT: MAP.152, then MAP.146, then this item). The mechanism. A
  page route that expects a slow answer (a full list of one phenomenon
  kind, a wide search, a large count) checks the tagged Redis page cache
  that PERF.104 built; on a miss it does not build the answer inside the
  request: it queues one RQ job that builds the answer, writes it to the
  same cache under the same tags (so later edits invalidate it the way
  PERF.104 already does), and returns at once with a job id. The job
  reports progress as counters and events in Redis (PERF.86, Option A:
  Redis is the only place a run lives), shows on /admin/queue beside the
  PERF.104 "Cache warm-up" jobs, and can be cancelled. A status route
  (JSON, polled by the page) returns state, the done and total counts
  and the ETA, and the finished answer's cache key. A route is marked
  slow by a cost rule in one place, not by a timer inside each page
  (default: the request is slow when its query would read more than a
  stated row budget, set in tuning, or when the cache miss took longer
  than a stated number of seconds the last time); two visitors asking
  the same thing share one job. If Redis is down the page falls back to
  building the answer in the request, as it does today, and says so.
  Done: the mechanism, the cost rule, the status route, the job on the
  queue page, and tests (a slow request queues exactly one job, a repeat
  request attaches to it, an edit invalidates the stored answer, Redis
  down falls back).
  Prerequisites: none. Related: PERF.86, PERF.99, PERF.104, PERF.75,
  UX.3.

- [ ] **PERF.106 The plan command on the run record, with stage timing, the killed-worker handler and job timeouts**
  Split from PERF.98 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). Move the stages of `planetgen
  plan` onto run records, units and barriers as PERF.97 did for the mass
  scatter: the skeleton as one unit, the phenomena and luminosity
  scatters as static chunks (each scatter its own queue item, run in
  order), and the other plan modes (`band`, `map-luminosity`,
  `habitability`, `nebula-field`, the redo list). This first step also
  adds what every later stage reuses: stage timing written into
  `generation_stage_runs` and `generation_stats` by the unit that
  completes the stage, `work_horse_killed_handler` marking the unit and
  the run at once, and a real `job_timeout` (a multiple of the stage's
  expected time) on every unit instead of -1. Done: `planetgen plan` in
  every mode runs on the new queue and gives the same galaxy for the
  same seed at any worker count; stage rows land in MySQL as today; a
  hard-killed horse fails the run at once; no unit has `job_timeout` -1;
  one test per mode.
  Prerequisites: PERF.97. Related: PERF.98.

- [ ] **PERF.107 The sector units and the stages around them: start, neighborhood, sectors, link, population and settle**
  Split from PERF.98 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). Move the `galaxy` command's
  sector stages (`start`, `neighborhood`, `sectors`) to units that are
  static chunks of sectors (about 4,000 units per stage at most, JSON
  arguments, the run's seed independent of worker count), and `link`,
  `population` and `settle` as one-shot units, all on the run record
  with the done-set barrier. Done: `planetgen galaxy` without the
  backfill runs on the new queue and gives the same sectors for the same
  seed at 1 and 4 workers; `--workers 1` uses the same path; test for
  each stage.
  Prerequisites: PERF.106. Related: PERF.98.

- [ ] **PERF.108 The backfill as queue units: chunks of at least 10,000 cells, one enumeration (folds in PERF.18)**
  Split from PERF.98 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). Move the `backfill` stage and its
  units onto the new queue: the unit is a chunk of the backfill's target
  sectors of at least about 10^4 cells (the per-cell draw is about 0.06
  ms, so smaller units are all overhead); the target set is enumerated
  once around the run's boundary and filtered by tier radius, not once
  per block (the 60-fold redundant list); the band and mass tables stay
  pre-warmed in the worker (`cli/worker.py`, PERF.42). This is PERF.18's
  remaining work, the old queue pointers of which
  (`stellarObjects/workQueue.py`, in-request backfills that no longer
  exist in `web/` or `api/`) are gone. Done: `planetgen galaxy` with the
  backfill runs wholly on the new queue and gives the same stars for the
  same seed at any worker count; the old executor is unused by any
  generation command; a test with a backfill of at least two chunks.
  Prerequisites: PERF.107. Related: PERF.98.

- [ ] **PERF.109 Cache the Generate page estimate**
  Split from PERF.20 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). The estimate (`--estimate-only`,
  up to 120 s) is repeatable and costly: cache it in front of the work,
  in the web or API process, keyed by (kind, arguments digest, epoch,
  speed-stats version), with the epoch from PERF.104's `epoch` tag and
  single-flight through `web/lib/singleflight.py`, so concurrent
  identical requests run one estimate (a queue round trip costs 2.4 s,
  so no RQ for cache reads). PERF.105's slow-route rule decides whether
  a miss is queued. Done: a second identical estimate request answers
  from the cache in milliseconds; a change of plan or epoch or new stage
  statistics misses; two concurrent identical requests run the estimate
  once; tests.
  Prerequisites: none. Related: PERF.20.

- [ ] **PERF.110 XFetch refresh for costly entries that expire by age**
  Split from PERF.20 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). For entries that expire by age
  and cost more than a stated build time (the stats page, search
  aggregates), recompute shortly before expiry with a probability that
  rises toward it (XFetch, beta 1), so a burst of visitors does not all
  miss at the expiry second. Done: a helper in the page cache with an
  injectable clock and random source, the stats page and the search
  aggregates on it, and a test that under 50 concurrent readers at the
  expiry time one recompute happens before it and none after.
  Prerequisites: none. Related: PERF.20.

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
  was an unverified estimate); after the cut it is about 2.7e5 rows.
  Decided (Boss, 2026-10-10 18:44Z, defaults approved; defaults taken):
  after a repair
  the orbit-updated columns come back as generated and the next orbit
  update carries on (a small daily snapshot of ids, phases and sector
  only if exact positions matter); G = 32, m = 2; repairing a bad
  clustered-index page may mean restarting the database server with
  `innodb_force_recovery=1`.
  Split (2026-10-16, Boss 1:07 PM EDT: look at each item for smaller
  chunks and shared work): this item is now the umbrella and is done
  when its subitems are: DB.43, DB.44, DB.45, DB.46, DB.47, DB.48. The
  largest item of the slice: a stdlib codec, a schema change, a checksum
  that every writer of exported columns must respect, a parity file kept
  current on every save, a repair command, and a table rebuild for an
  unreadable clustered-index page (nine research points in its text).
  Prerequisite GEN.57 dropped (2026-10-16 review): only the
  regenerate-from-seed fallback, DB.17, needs it; the parity half needs
  DB.8 only.
  Prerequisites: DB.43, DB.44, DB.45, DB.46, DB.47, DB.48.
  Design: [docs/design/reproducible-galaxies.md](design/reproducible-galaxies.md)

- [ ] **DB.16 Store the generator epoch and run id on each sector instead of four version text columns**
  DB.7 stored `version_key`, the version, python and platform as text on
  every sector; at 12 billion sectors those strings repeat. Store
  `generator_epoch` (SMALLINT) and `run_id` (foreign key to the existing
  `generation_runs`) instead, with the key text looked up from the run.
  Decided (Boss, 2026-10-10 18:44Z, defaults approved; default: epoch
  plus run id): or keep the full
  key text per sector as DB.7 shipped. ALGORITHM=INSTANT, idempotent
  migration, per db-check-and-parity-repair.md.
  Prerequisites: OPS.56.
  Design: [docs/design/reproducible-galaxies.md](design/reproducible-galaxies.md)

- [ ] **DB.17 Repair by regenerating a damaged sector from its seed when parity cannot rebuild it**
  The second half of DB.9: when the parity file cannot rebuild a sector,
  regenerate it from the galaxy seed, the sector's stored directives and
  the current galaxy's internal record (ADM.50), then replay its edits
  from the edit log.
  (Boss, 2026-10-09 20:42Z, dropped the pending-delta JSON and daily
  merge this item used to read; the seed is still used internally for
  repair. Boss, 2026-10-09 20:52Z: keep the repair from the seed, replaying
  the edit log.) Compare the result with the stored leaf digest. Decided
  (Boss, 2026-10-10 18:44Z, defaults approved; default
  yes): ship the parity half first, without GEN.57, GEN.58 and OPS.14.
  Note (2026-10-09): admin regenerate of a planet, moon or belt
  (`admin/edits.py`) draws from the process stream, so an edited sector
  cannot be rebuilt from its seed alone, which is why this repair
  replays the edit log; the web handlers were not checked (Research Lane
  1).
  Also needs GEN.179 (the stored directives it reads) and ADM.50 (the
  galaxy's internal record), both moved into this thread before it
  (2026-10-16 review).
  Prerequisites: GEN.57, DB.47, GEN.179.
  Design: [docs/design/db-check-and-parity-repair.md](design/db-check-and-parity-repair.md)

- [ ] **DB.29 update.sh and every database migration refuse to start while jobs are active, unless --yes is passed or the user confirms cancelling them**
  Boss, 2026-10-15 02:55Z: "Add TODO item, when migrating the database
  to a new version, if there are active jobs, ask the user to cancel the
  jobs first before changing anything." Widened by Boss, 02:56Z: "In
  fact, anytime update is run, it won't start without either --yes
  passed OR if the user specifically confirms they want to cancel all
  jobs prior to starting the upgrade or install." Done: before it
  changes anything, update.sh (every install and every upgrade, not only
  runs that migrate) and the migrate command check for active jobs
  (running or queued generation runs and other background jobs). With
  none, they proceed as before. With any, they list them and stop until
  either `--yes` was passed, or the user types an explicit confirmation
  that all active jobs may be cancelled; only then are the jobs
  cancelled, and then the install or upgrade starts. Boss, 02:59Z: "the
  update script needs to know if it's being run either from systemd or
  cron, and if it has, assume no unless --yes is passed. Or say 30
  seconds of waiting for input and it assumes no and cancels the
  update." So update.sh detects a run from systemd or cron (no person
  at the keyboard) and then assumes "no": it changes nothing and exits
  non-zero unless `--yes` was passed. When it does prompt, 30 seconds
  without an answer count as "no": it changes nothing and exits
  non-zero. Without `--yes` and without a terminal to ask on, it does
  the same.
  update.sh still never edits the user's Redis configuration (OPS.48).
  A test covers each path (no jobs, jobs and `--yes`, jobs and a typed
  confirmation, jobs and no answer within 30 seconds, jobs and a run
  from systemd or cron), and the check is the one place every
  migration goes through.
  Dependency: how the active jobs are read depends on the work queue
  steps. PERF.96 (run record and counters in Redis, done, PR #1553)
  supplies the library call that lists active runs, which the check
  uses; and PERF.99 supplies the cancel through the stop
  command, after which the check uses them. PERF.100 (which deletes the
  old work tables with a migration) must pass through this check too.
  Owner: Database, after DB.28 (done, PR #1476) and before PERF.100.
  Prerequisites: PERF.99. Related: PERF.95 (done, PR #1506), PERF.96
  (done, PR #1553), PERF.99, PERF.100, OPS.48, DB.18, DB.27 (done, PR
  #1446), DB.28 (done, PR #1476).
  Owner: Core 1, in its bug phase after the global top chain (Boss,
  2026-10-16, 2:43 PM EDT: "Core 1 will start working on bugs after";
  order set at 3:00 PM EDT: non-test bugs, then the database-version
  test failures, then the failures that hit full runs, then the one-off
  load flakes). Fix the cause; no skipping or quarantine.

- [ ] **DB.30 One migration for the per-species environment limits**
  Needed by POP.12 (done, PR #1554; Generate), same rule as DB.27 and
  DB.28: only Database adds migrations. Add the columns or
  the table that hold each species' own environment limits, one value
  per row of `LIMIT_ROWS` in `planetgen/physics/protection_limits.py`
  (pressure, oxygen, temperature, wet-bulb, dose and the other rows it
  holds), NULL for the human default, as laid out in
  docs/design/habitability-index.md section 9 ("a later per-species
  limit replaces a row in this table and nothing else"). Done: the
  migration (next after DB.29's if that lands first) with its old-schema
  fixture and tests, merged so the species limits of POP.12 (done, PR #1554) can be stored. Owner:
  Database, after DB.29.
  No longer waits on DB.29 (2026-10-16 review): DB.29 moved into the
  queue thread before PERF.100, and DB.30 is the urgent one for POP.12.
  Prerequisites: none. Related: POP.12 (done, PR #1554), POP.13, DB.27.

- [ ] **DB.31 Store structure, not prose, for life and flavour**
  Audit: /mnt/project-files/wording/text-output-audit-2026-10-15.md.
  Owner lane: Foundations Wording (2026-10-15), which writes the code
  for this item. Replace planet_evolutionary_paragraphs and
  moon_evolutionary_paragraphs text with the stage, the pace and the
  ages; flavor_text with a flavour key; drop hab_note (computed) and the
  run-config markdown column. Then life_stage_from_paragraphs and the
  string searches in decide_flavor_text go. Changes what generation
  stores, not the galaxy (draws are unchanged). Database
  owns the migration and agrees its timing; Foundations Wording writes
  the code against it. The phi4* column names stay (decided default); a
  rename to phi6* is not planned.
  Split (2026-10-16, Boss 1:07 PM EDT: look at each item for smaller
  chunks and shared work): this item is now the umbrella and is done
  when its subitems are: DB.50, DB.51, DB.52. It mixes two migrations
  (adding the structure columns, then dropping three kinds of prose
  columns and two tables) with a cross-cutting code change that another
  lane writes ('Database owns the migration and agrees its timing;
  Foundations Wording writes the code against it').
  Prerequisites: DB.50, DB.51, DB.52.

- [ ] **DB.33 MySQL 8.0.46 reflects the foreign keys fk_facilities_type and fk_facility_field_values_field as RESTRICT, so check-db calls a healthy database damaged (bug)**
  Reported by Quality (2026-10-15, 4:56 PM EDT run, relayed by the
  coordinator). On MySQL 8.0.46 the foreign keys `fk_facilities_type`
  and `fk_facility_field_values_field` reflect with `ondelete` RESTRICT
  while the models declare none, so `check-db` reports a healthy
  database as damaged and six tests fail (test_db_check, test_db_models,
  test_db_check_deep). Quality confirmed (2026-10-15, combined run) that
  it fails alone on a clean main, so no lane's change causes it; the six
  tests are test_a_healthy_database_passes,
  test_a_stale_system_count_is_a_warning_not_damage,
  test_the_command_exits_one_on_damage_and_zero_when_clean,
  test_models_match_a_database_built_from_schema_sql,
  test_a_changed_column_is_noticed and
  test_a_deep_check_off_a_terminal_needs_yes. Declaring
  ondelete='RESTRICT' alone would fix 8.0.46 but break MySQL 8.4 and
  MariaDB, which reflect none, so the compare must treat all three as
  equal. Fix: compare RESTRICT, NO ACTION and none as equal in
  `db/check.py` `_model_drift` and in `db/models.py`. Done: check-db
  passes on a healthy MySQL 8.0.46 database, the six tests pass, and a
  real drift in a foreign key's delete rule is still reported. Owner:
  Core (2026-10-15).
  Update (Core 1 and Core 2, 2026-10-15): the owner is now Core 2 (Core
  1 handed it over). Eight tests fail on MySQL 8.0.46, also
  test_the_check_writes_nothing (KeyError 'table_name':
  information_schema column names come back upper case, the test reads
  row["table_name"]) and test_a_deep_check_off_a_terminal_needs_yes;
  they fail on main alone, always. Regenerating models.py with
  scripts/generate_db_models.py on that server makes the drift failures
  vanish, so the fix may also be a regenerate plus the case-insensitive
  read.
  Quality (2026-10-15, #1697): the test_db_check column-name case fix
  (upper-case information_schema names) is merged; the RESTRICT drift
  itself is still open.
  Seen again (Quality, 2026-10-15): the six DB.33 tests failed again in
  the `changed` run, same cause (fails on main).
  Confirmed again by Core 3 (2026-10-16, 9:05 PM EDT, GEN.149 test run):
  the same six tests fail on clean origin/main (5b48f34) on every run (2
  of 2), with the same remove_fk and add_fk items for fk_facilities_type
  and fk_facility_field_values_field. Owner: Core 2 (in its queue
  after the plan items; Core 2 confirmed 9:17 PM EDT it is its to fix).
  Owner: Core 1, in its bug phase after the global top chain (Boss,
  2026-10-16, 2:43 PM EDT: "Core 1 will start working on bugs after";
  order set at 3:00 PM EDT: non-test bugs, then the database-version
  test failures, then the failures that hit full runs, then the one-off
  load flakes). Fix the cause; no skipping or quarantine.
  Prerequisites: none.

- [ ] **DB.36 Tables for change plans, their operations and the edit journal (migration only)**
  Boss (2026-10-16, 1:07 PM EDT), headline of 8.6 (ADM.57): the admin
  has complete control; the world adjusts to what the admin wants, after
  the admin approves the generation plan. One Alembic migration adds
  `change_plans` (who asked and in which role, the request, the state:
  draft, requested, proposed, approved, applied, undone, stale, rejected
  or withdrawn, the role that must approve it, who approved or rejected
  it and when with a note, the job that applied it; 10.0 adds only
  request comments and notices, DB.39, so the plan tables hold requests
  from the start and need no rebuild), `change_plan_ops` (plan, order,
  operation
  kind, target object ID, before and after values as JSON, status) and
  `edit_journal` (an append-only list of every applied operation: plan,
  object ID, aspect, before, after and the generator epoch, indexed by
  object ID and by sector). Core 1 owns the migration number. Done: the
  migration, the regenerated models, an old-schema fixture, an
  idempotent upgrade and a round-trip test; no code reads the tables
  yet.
  Prerequisites: none. Related: ADM.57, ADM.58, ADM.59.
  Design: [docs/design/edit-everything.md](design/edit-everything.md)

- [ ] **DB.37 Lock marks table (migration only)**
  Boss (2026-10-16, 1:07 PM EDT), headline of 8.6 (ADM.57): the admin
  has complete control; the world adjusts to what the admin wants, after
  the admin approves the generation plan. One Alembic migration adds
  `object_locks` (object ID, aspect or all, locked by, locked at,
  reason): a lock is a standing keep mark. Core 1 owns the migration
  number. Done: the migration, models and an idempotent upgrade with a
  round-trip test.
  Prerequisites: none. Related: ADM.57.
  Design: [docs/design/edit-everything.md](design/edit-everything.md)

- [ ] **DB.38 Region density overrides table (migration only)**
  Boss (2026-10-16, 1:07 PM EDT), headline of 8.6 (ADM.57): the admin
  has complete control; the world adjusts to what the admin wants, after
  the admin approves the generation plan. One Alembic migration adds
  `region_density_overrides` (region shape and parameters, multiplier
  from 0.1 to 10, set by, set at, undone at), append-only, for the
  density edits of MAP.119 and ADM.57. MAP.119's block and slab case is
  a region too. Core 1 owns the migration number. Done: the migration,
  models and an idempotent upgrade with a round-trip test.
  Prerequisites: none. Related: MAP.119, ADM.57.
  Design: [docs/design/edit-everything.md](design/edit-everything.md)

- [ ] **DB.39 Edit request support tables (migration only): request comments and notification records**
  Boss (2026-10-16, 1:43 PM EDT), headline of 10.0 (USR.25): three
  roles, User, Moderator and Admin; users submit edit requests,
  moderators approve user requests, admins approve moderator requests,
  and an admin's own edits need no approval (only the preview and
  confirm of ADM.66). Build on the accounts of 8.2 and the change plans
  of 8.6. The requester, the required approver role and the decision are
  columns of the change plan tables (DB.36), so a request is a plan. One
  Alembic migration adds `request_comments` (request, author, body,
  created) and `request_notifications` (request, recipient, kind, sent).
  Core 1 owns the migration number. Done: the migration, models and an
  idempotent upgrade with a round-trip test.
  Prerequisites: none. Related: USR.25, DB.36.
  Design: [docs/design/roles-and-requests.md](design/roles-and-requests.md)

- [ ] **DB.40 Parity data backup part: the parity file of the DB repair copied into the backup archive**
  Boss (2026-10-16, 1:59 PM and 2:00 PM EDT): "let's add backups of the
  control DB and backups of the parity data for DB repair." and "All
  those should be bz2 tarballed together." Boss (2:00 PM EDT): "we also
  have to build a way to use that data so build that out and let's make
  this targetted for v9.0 so incorporate the chain for these (restore /
  repair using backups) with that due date. Settings page is v8.1".
  Target 9.0 (the backups and repair chain, umbrella ADM.91); the
  settings page and the settings backups stay in 8.1 (ADM.83 to ADM.85).
  All backup pieces are parts of one `.tar.bz2` archive (ADM.87);
  retention is `settings.backups_kept` (default 20, 1 to 200). DB.9
  keeps parity for every sector in an append-only parity file outside
  the database (path in `config.json`). A damaged database is only
  repairable with a parity file that matches it, so the file is a backup
  part, written to `parity/`: the copy is taken under the parity
  writer's lock and cuts at the last complete record (the file is an
  append-only log of self-describing records with a header scan, so no
  index is stored); it can be switched off for a very large galaxy with
  the option `settings.backup_parity` (default on), and the manifest
  records its size and the sector count it covers. Done: the part
  registered with ADM.87, the lock-held copy under a concurrent save,
  and the size and coverage in the manifest (the restore is DB.41).
  Prerequisites: DB.46. Related: ADM.85, DB.8, DB.17.

- [ ] **DB.41 Parity data restore: put a backed-up parity file back and check it against the database**
  Boss (2:00 PM EDT): "we also have to build a way to use that data so
  build that out and let's make this targetted for v9.0 so incorporate
  the chain for these (restore / repair using backups) with that due
  date. Settings page is v8.1". Target 9.0 (the backups and repair
  chain, umbrella ADM.91); the settings page and the settings backups
  stay in 8.1 (ADM.83 to ADM.85). All backup pieces are parts of one
  `.tar.bz2` archive (ADM.87); retention is `settings.backups_kept`
  (default 20, 1 to 200). Using the parity part of an archive: restore
  writes the file back (to a temporary name first, then renamed), then
  runs the DB.8 check so every record that no longer matches the
  database is recomputed by the stale-record rule of DB.9 instead of
  being read as damage, and reports how many sectors are unprotected
  until their next save. `planetgen.cli.backup restore --part parity
  ARCHIVE` and the backup list use it. Done: a restore on a fixture
  galaxy followed by a passing check, a record from before an edit
  treated as stale and not as damage, and a refused restore while a save
  is writing parity.
  Prerequisites: DB.40, DB.8. Related: DB.9, ADM.85, ADM.91.

- [ ] **DB.42 Repair the galaxy from an archive: `repair-db --from-archive` uses the parity part of a backup when the live parity file is missing or damaged**
  Boss (2:00 PM EDT): "we also have to build a way to use that data so
  build that out and let's make this targetted for v9.0 so incorporate
  the chain for these (restore / repair using backups) with that due
  date. Settings page is v8.1". Target 9.0 (the backups and repair
  chain, umbrella ADM.91); the settings page and the settings backups
  stay in 8.1 (ADM.83 to ADM.85). All backup pieces are parts of one
  `.tar.bz2` archive (ADM.87); retention is `settings.backups_kept`
  (default 20, 1 to 200). The galaxy database itself is not in the
  archive (too large; it can be regenerated from the seed, DB.17, and
  repaired from parity, DB.9). `repair-db` of DB.9 reads the live parity
  file; with `--from-archive ARCHIVE` (or the choice on the repair page)
  it reads the archive's parity part instead when the live file is
  missing, damaged or newer than the damage. It applies the stale-record
  decision table against the live database, rebuilds every damaged
  sector it can, and lists the sectors it cannot repair (a stale record,
  more than m members of a group lost) for the seed regeneration of
  DB.17, with a report of what came from which source. Done: tests on a
  copy of a small galaxy: the live parity file deleted and the damage
  repaired from the archive, a group that lost too many members listed
  and handed to DB.17, and a stale archive recompute that does not
  damage a healthy sector.
  Prerequisites: DB.41, DB.9, DB.17. Related: DB.9, DB.17, ADM.91.

- [ ] **DB.43 Reed-Solomon codec in the standard library: planetgen/util/rscodec.py**
  Split from DB.9 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work). A GF(256) Reed-Solomon module with no
  new dependency: a Cauchy generator scaled so row 0 is plain XOR,
  `RS(k, m)` with `encode` (k equal-length blocks to m parity blocks),
  `decode` (any k of the k+m blocks back to the data) and `update` (new
  parity when one member changes, from the old and new bytes). The
  tested code is in docs/design/db-check-and-parity-repair.md 5.9. Why:
  the parity file (child d) and the repair command (child e) both stand
  on it, and it can be built and tested with no database. Done: the
  module and its tests: random erasure patterns of up to m blocks decode
  exactly for (k, m) of (16,1), (16,2), (32,2) and (8,3); `update`
  equals a fresh `encode`; unequal member lengths work through zero
  padding; the time to encode and to rebuild one 32-member group of 43
  KB blocks is recorded in the PR.
  Prerequisites: none. Related: DB.9.

- [ ] **DB.44 Migration only: the sector storage checksum and the parity dirty flag**
  Split from DB.9 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work). One Alembic revision adds to
  `sectors` a `storage_checksum` (BINARY(32), NULL) and a `parity_dirty`
  flag (TINYINT, NOT NULL, default 1, so every existing sector counts as
  not yet protected). Written with the DB.18 helpers (ALGORITHM=INSTANT,
  idempotent). Core 1 owns the migration number. Why: the checksum and
  the flag are the two columns the rest of DB.9 reads, and a migration
  is Core 1's alone. Done: the migration, the regenerated models, an
  old-schema fixture, an idempotent upgrade and a round-trip test; no
  code reads the columns yet.
  Prerequisites: DB.18. Related: DB.9.

- [ ] **DB.45 Storage checksum and the dirty flag: the export, the checksum and every writer**
  Split from DB.9 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work). Compute the sector's storage checksum
  as SHA-256 over a lossless, id-preserving dump of its stored rows,
  read back after the write, leaving out the columns the daily orbit
  update rewrites (phases, positions, velocities, reflex offsets, binary
  positions, comet anomalies, a placed object's `sector_id`; list in
  design section 5.2). Every writer of exported columns sets
  `parity_dirty` in its own transaction: the direct UPDATEs in
  `db/store.py` (names, positions, velocity), `api/edits.py` and
  `db/edits.py`; finding every writer is part of this item. DB.8's check
  reads the checksum and reports a dirty sector as "unprotected", not
  damaged. Why: the parity file can only be correct if a changed row is
  always noticed. Done: the export and checksum, every writer marked,
  and tests: each writer sets the flag; an orbit update does not change
  the checksum; a flipped stored value is reported as damage and a dirty
  sector as unprotected.
  Prerequisites: DB.44. Related: DB.9.

- [ ] **DB.46 The parity file: groups, records, a delta update on every save and the stale-record rule**
  Split from DB.9 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work). A parity file kept outside the
  database (path in `config.json`) with G = 32 members and m = 2 parity
  blocks per group (both in `config.json`), groups formed at save time
  from sectors of one export-size class that are far apart in sector id
  (the W = 256 rule, shrinking for small galaxies); an append-only log
  of self-describing records (magic, format version, group id, `seq`,
  per-member address, length and checksum, parity bytes, CRC32C), an
  in-memory index rebuilt by a header scan, compaction as a function.
  Each save updates the group's parity from the member's before-image
  inside the saving transaction (`P_j' = P_j xor c_ji x (old xor new)`),
  after checking the before-image against its old checksum, then clears
  `parity_dirty`. The stale-record decision table of design 5.6 is code.
  Why: this is the part that makes repair possible and it must cost
  milliseconds per save. Done: the file, the update on save, and tests:
  save, edit and delete keep parity equal to a fresh encode; a torn tail
  loses nothing; a stale record is recomputed, not treated as damage;
  overhead is about 7 percent on a 1,000-sector fixture.
  Prerequisites: DB.43, DB.45. Related: DB.9.

- [ ] **DB.47 planetgen repair-db: rebuild damaged sectors from the parity file**
  Split from DB.9 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work). `generate.py repair-db` (the
  `planetgen` command) uses DB.8's check to find damaged sectors,
  applies the decision table, rebuilds each sector from its group's
  parity (reading the other members), re-runs the check and lists
  anything it could not repair (for DB.17's seed fallback). After a
  repair the columns the orbit update rewrites come back as generated
  and the next update carries on (decided, 2026-10-10). Why: the
  user-facing half of DB.9. Done: a test damages stored rows in a copy
  of a small galaxy, repairs them and gets a passing check; a sector
  whose record is stale or whose group lost more than m members is
  listed, not guessed.
  Prerequisites: DB.46. Related: DB.9.

- [ ] **DB.48 repair-db: rebuild the table when a clustered-index page is unreadable**
  Split from DB.9 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work). When a page of the `sectors`
  clustered index cannot be read, rebuild the table: restart the
  database server with `innodb_force_recovery=1` (the command says so
  first and asks), dump the surviving rows, recreate the table, reload
  with foreign-key checks off, and insert the rebuilt sectors with their
  ids (design section 4). Why: single-sector rebuild in place does not
  work for a bad index page, and it is the one repair that needs a
  server restart, so it is a separate decision for the admin. Done: a
  test flips bytes in an `.ibd` page of a copy while the server is
  stopped, runs the repair and gets a passing check with every sector
  present.
  Prerequisites: DB.47. Related: DB.9.

- [ ] **DB.49 Migration only: the leaf digest column and the shell-band tree table**
  Split from GEN.136 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). One Alembic revision adds
  `sectors.leaf_digest` (BINARY(32), NULL) and a `sector_tree_nodes`
  table for the tree over the shell, band and wedge ladder: a band node
  (shell, band) over its wedge leaves, a shell node over its bands, the
  galaxy root on demand. Written with the DB.18 helpers; Core 1 owns the
  number. Done: the migration, models, an old-schema fixture, an
  idempotent upgrade and a round-trip test; no code reads them yet.
  Prerequisites: DB.18. Related: GEN.136.

- [ ] **DB.50 Migration only: structure columns for life and flavour**
  Split from DB.31 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work). One Alembic revision adds, on
  `planets` and `moons`, the structure that replaces the evolutionary
  paragraphs (the life stage, the pace and the ages; the PR names the
  columns) and a flavour key that replaces `flavor_text`. Nothing is
  dropped. Core 1 owns the migration number. Done: the migration,
  regenerated models, an old-schema fixture, an idempotent upgrade and a
  round-trip test; no code reads the columns yet.
  Prerequisites: none. Related: DB.31.

- [ ] **DB.51 Generation stores structure, not prose; the readers and the string searches go**
  Split from DB.31 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work). Generation writes the stage, pace,
  ages and flavour key; readers move to them: `db/query._life_stages`,
  `population/model.py`, `db/edits.py`, `db/fingerprint.py`;
  `life_stage_from_paragraphs` and the string searches in
  `decide_flavor_text` are deleted; `hab_note` is computed from the
  describer's rating rows (GEN.218) instead of stored and parsed by
  regex in `generation/planet.py`; the run-config
  `system_configs.markdown` column is no longer read. The fingerprint's
  content changes, so it bumps the generator epoch if OPS.28 child a has
  landed (a PR rule, not a build blocker). Draws are unchanged. Done:
  the same seed gives the same stage, pace, ages and flavour as before
  the change (a test compares them with the old paragraphs on a fixture
  galaxy); nothing reads the old columns or tables.
  Prerequisites: DB.50, GEN.218. Related: DB.31.

- [ ] **DB.52 Migration only: drop the prose columns and tables**
  Split from DB.31 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work). One Alembic revision drops
  `planet_evolutionary_paragraphs`, `moon_evolutionary_paragraphs`,
  `flavor_text` (keep `flavor_text_count` only if the flavour key logic
  still reads it), `hab_note` on planets and moons and
  `system_configs.markdown`. Core 1 owns the number. Done: the
  migration, regenerated models, the old-schema fixture updated, an
  idempotent upgrade and a round-trip test; the models and `schema.sql`
  have none of them.
  Prerequisites: DB.51. Related: DB.31.

- [ ] **DB.53 Built-system star storage: the orbit_nodes table, stars.parent_node_id, planets.host_node_id, the catalogue columns and the disk descriptor (one Core 1 revision)**
  Split from GEN.129 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). A built multi-star system has
  nowhere to keep its tree: scatter_nodes exists (migration 0093) but
  the built side, orbit_nodes, does not. This item is the one Core 1
  revision for everything the multi-star train stores on built systems,
  so GEN.264, GEN.269, GEN.273, GEN.274, GEN.242 and GEN.275 do not each
  ask for their own. Create orbit_nodes with the same names and units as
  scatter_nodes (system id, node_no, parent_node_no, SMALLINT slot,
  mass_solar, separation_km, eccentricity, inclination_deg,
  ascending_node_deg, arg_periapsis_deg, mean_anomaly_deg, epoch_unix,
  period_years, kind close, wide, compact or group, interaction, flags)
  so a fill copies column for column; add stars.parent_node_id and
  planets.host_node_id; add to the built star rows the member columns
  scatter_members already has (state, flags, offset_*_au,
  velocity_*_kms, decay_age_yr for the members of a group node); and add
  one disk descriptor (disk_kind, disk_mass_solar, disk_inner_au,
  disk_outer_au) on a star and on an orbit node, nullable, no bodies in
  it. A plain binary converts to one node (the migration converts
  existing galaxies, or the galaxy is rebuilt: Boss, 2026-10-15 5:27 PM
  EDT, early alpha). Done: the revision applies on a fresh and on an
  existing galaxy, check-db knows the new tables and columns, and a test
  pins column parity with scatter_nodes and scatter_members so a copy is
  column for column.
  Shared by GEN.242, GEN.264, GEN.269, GEN.273, GEN.274 and GEN.275,
  which each assumed columns on tables that do not exist yet; one Core 1
  revision for all of them.
  Owner: Core 1 (Boss, 2026-10-16, 2:42 PM EDT and 3:00 PM EDT: the
  global top blockers; after the GEN.335, DB.55, GEN.336, GEN.337,
  GEN.294 chain).
  Migration revision (Core 1, 2026-10-16, 6:13 PM EDT): 0104, after
  0103 (GEN.336).
  Prerequisites: none. Related: GEN.129.

- [ ] **DB.55 One additive migration for the 8.3 schema: point columns for every kind, heaviest object per sector, asteroid-field provenance**
  Split from GEN.294 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). Core 1 owns every revision, and
  four items of 8.3 each asked for columns on the same tables (GEN.294,
  GEN.142's velocities, GEN.110's mass_kg and origin on asteroid_fields,
  GEN.109's per-sector heaviest object). One revision, the next free
  after 0098, adds them nullable and backfilled: x, y, z
  (sector-relative, with the cell key where a table lacks it), vx, vy,
  vz, mass_kg, epoch_unix and next_update_due wherever missing for black
  holes, neutron stars, rogue planets, interstellar comets, asteroid
  fields, belts, nebulae, remnants, quasars, comets, facilities and
  scatter points and members (the vector backfilled from the scalar
  speed as a tangential vector of the same speed, the mass by the stated
  rule); `sectors.max_mass_kg` and `max_mass_id`;
  `asteroid_fields.formed_at` and `origin_event_id`; the planet and moon
  state columns GEN.289 needs (parent_type, parent_id, dynamical_state,
  vector, mass_kg). The scalar speed columns and the circular phase
  columns are dropped by the code items that stop reading them, each
  with its own small revision from Core 1. Why: no other lane then edits
  Alembic, and no two items add the same column. Done: the revision
  applies on a galaxy at the previous head and on an empty one,
  downgrades as the other revisions do, and the schema docs list the new
  columns.
  Owner: Core 1 (Boss, 2026-10-16, 2:42 PM EDT: Core 1 takes the global
  top 5 blockers, starting right after GEN.237 merges; estimates from
  Core 1: GEN.335 merged about 7:30 PM EDT, DB.55 about 10:30 PM EDT).
  Prerequisites: none. Related: GEN.294.

- [ ] **DB.56 The orbital_events table (migration only): collision outcomes, mergers, splits, ejections, captures, unbinds and warnings**
  Split from GEN.143 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). One row per event in
  `collision_events` with the columns of collisions-and-mergers.md
  section 7 (event id, run id, impact and detection time, outcome code,
  regime, seed string; each participant; the numbers; the conservation
  check; the products; context and flags including `numeric_warning`),
  and a `record_collision(conn, event)` writer. Core 1 places the
  revision. Done: the revision applies, an event written by the writer
  reads back whole, no view yet.
  Placed before child b of GEN.110, which writes its events here
  (2026-10-16 review).
  Scope after the review (2026-10-16 review): one table for every
  orbital event kind (collision outcome, merger, split, ejection,
  capture, unbind, warning) with id, run id, event time, detected time,
  kind, outcome code, participant ids, product ids, seed string and
  queried numbers as columns, plus `destroyed_by_event_id` and
  `origin_event_id` on the objects. Core 1 names the revision (in the
  window of the migration child of GEN.294 and GEN.289). The writer is
  GEN.197.
  Prerequisites: none. Related: GEN.143.

- [ ] **DB.57 The update records what it moved: a replace-each-run table of previous and new positions**
  Split from MAP.192 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). The advance step writes one row
  per object it moved (table, uid, previous and new position, epoch
  before and after) into a table that the next update replaces, bounded
  by the run's row ceiling (PERF.103). Core 1 places the revision. Done:
  after an update the table lists exactly the objects counted as moved,
  a second update replaces it, and a test on a small galaxy checks both.
  Prerequisites: GEN.294, GEN.212, GEN.340. Related: MAP.192.

- [ ] **DB.58 Control-schema revision for accounts: role, email and status, the single-Owner index, account_tokens and invite_uses**
  Split from USR.2 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work). USR.2, USR.4, USR.5 and USR.6 each
  need control-database tables; Core 1 places one revision (the next
  control version) instead of three. It adds role (user, admin, Owner),
  email and status (active, disabled) to the accounts table (Boss,
  2026-10-10: add the columns now, rename `admin_users` to `accounts`
  later), a database-level single-Owner unique index on a generated
  column, `account_tokens` (hashed token, purpose, account, created,
  expires, used_at: the shared table for USR.4, USR.5 and USR.6) and
  `invite_uses`; existing admins get role admin and the lowest id
  becomes the Owner. Done: the revision applies on a control database
  with several admins and on a fresh install (install creates exactly
  one Owner), a second Owner insert fails in MySQL, and the schema
  header comments are updated.
  Prerequisites: none. Related: USR.2.

- [ ] **DB.60 Migration: record the density multiplier used on each filled sector (migration only)**
  Split from MAP.119 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). One Alembic migration adds
  `density_multiplier_used` (decimal) and `override_id` to the
  filled-sector record, so a rebuilt galaxy can tell which override a
  sector was filled under (map-ui-and-frontend-libraries.md section 8).
  Core 1 owns the migration number. Done: the migration, models, an
  old-schema fixture, an idempotent upgrade and a round-trip test; no
  code writes the columns yet.
  The stored density multiplier is also added to the stored directive of
  GEN.179 (or to the DB.38 migration); say so in both items (2026-10-16
  review).
  Prerequisites: DB.38. Related: MAP.119.

- [ ] **DB.61 Migration: unnamed clusters and a star-to-cluster parent link (migration only)**
  Split from GEN.209 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). One Alembic migration makes
  `galaxy_components.name` nullable (the unique key stays over named
  rows) and adds a nullable cluster link, a reference to
  `galaxy_components.id`, on `stars` and on `bright_stars`, indexed, so
  a star can point to its cluster as its parent. Core 1 owns the
  migration number and settles the column names. Done: the migration,
  models, an old-schema fixture, an idempotent upgrade and a round-trip
  test; existing clusters keep their names and no code writes the link
  yet.
  Prerequisites: none. Related: GEN.209.

- [ ] **DB.62 Upload bookkeeping tables in the control database (migration only)**
  Split from API.11 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work). One Alembic migration on the control
  database (control v12 in api-design-standards.md section 7.3) adds
  `upload_runs` (id, owner, key id, galaxy database, state, created,
  last contact), `upload_batches` (run, seq, sha256, bytes, received
  at), `upload_units` (run, unit key, state open, sealed, verified,
  finalized or rejected, manifest hash, errors), `upload_claims` (run,
  sector address, state) and `upload_id_ranges` (run, table, first id,
  last id). Core 1 owns the migration number. Done: the migration,
  models, an old-schema fixture, an idempotent upgrade and a round-trip
  test; no code reads the tables yet.
  Prerequisites: none. Related: API.11.

- [ ] **DB.63 Table for ships and stations, user-owned or generated (migration only)**
  Split from USR.10 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work). One Alembic migration adds the table
  for vessels and stations: object ID (a new object-kind code, shared
  with the kind list of names/object_id.py), kind (ship or station),
  name, an owner that is a user or NULL for a generated one, the
  SpatialPosition3D columns of GEN.294 (position, point mass, movement
  vector), a course (the route plan, JSON) and creation data. Core 1
  owns the migration number and settles the columns. Done: the
  migration, models, an old-schema fixture, an idempotent upgrade and a
  round-trip test; no code reads the table yet.
  Prerequisites: GEN.294. Related: USR.10.

- [ ] **DB.64 Table for stored constellations (migration only)**
  Split from VIEW.3 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work). One Alembic migration adds the
  per-system constellation store: system uid, constellation index, name
  object ID (the `constellation` kind, VIEW.4), the member stars and the
  edges (JSON), and the seed and version it was built under. Core 1 owns
  the migration number. Done: the migration, models, an old-schema
  fixture, an idempotent upgrade and a round-trip test; no code writes
  it yet.
  Prerequisites: none. Related: VIEW.3.

- [ ] **DB.65 Clamp the cell key index at 2^21-1 in `cellkey.py` (`cell_index`, `cell_keys`) and in the SQL expression the migration generates**
  From the GEN.372 precision audit (Core 1,
  docs/design/numeric-precision-audit.md section 5, 2026-10-16, 7:39 PM
  EDT), filed by the TODO thread; verdict: 64-bit is enough for every
  stored value. Done: `cell_index` and `cell_keys` in
  `src/planetgen/.../cellkey.py` clamp the per-axis index at 2^21-1 so a
  coordinate past the reach cannot wrap into another cell's key, and the
  generated SQL expression does the same (one migration, which needs a
  revision from Core 1); a test pins that a point beyond the reach maps
  to the edge cell in both the Python and the SQL form. Small. Owner:
  Core 1 (needs a migration revision).
  Migration revision: 0105, after DB.53's 0104 (Core 1, 2026-10-16, 7:42 PM EDT;
  numbers are assigned at merge time against main).
  Prerequisites: none. Related: GEN.372, DB.55.

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

  Prerequisites: API.30, ADM.95.
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
    When the claims exist, applying a change plan (ADM.60) refuses a
    sector an upload run has claimed, as server-side runs do.
    Needs the upload bookkeeping tables (first child of API.11), which
    it uses (`upload_runs`, `upload_claims`, `upload_id_ranges`)
    (2026-10-16 review). Added to Done: applying a change plan refuses a
    sector an upload run has claimed, with a test.
    Prerequisites: DB.62.

  - [ ] **API.11 Staging tables**
    Done: staging storage (a galaxy schema migration) for received
    batches, flagged incomplete until a whole unit (a system, a sector)
    has arrived and passed API.8's checks, then copied into the real
    tables in one transaction.
    Research (2026-10-09, api-design-standards.md): re-scope from
    "staging storage (a galaxy schema migration)" to upload bookkeeping
    tables in the control database plus a disk spool (`uploads.dir`),
    unless Boss specifically wants rows staged in the galaxy schema
    (about 30 content tables, 65 migrations). Decided (Boss, 2026-10-10
    18:44Z, defaults approved; default: spool plus bookkeeping). Verify
    and finalize of a sealed
    unit run as an RQ job (200 if done within the short wait, else 202
    with a job id; default yes).
    Split (2026-10-16, Boss 1:07 PM EDT: look at each item for smaller
    chunks and shared work): this item is now the umbrella and is done
    when its subitems are: DB.62, API.29. API.11 mixes a
    control-database migration with the spool and unit-state code, and
    it has an order defect: API.10 (the first item, no prerequisite)
    already needs `upload_runs`, `upload_claims` and `upload_id_ranges`
    (api-design-standards.md section 7.2 and 7.3), which API.11 creates
    after it.
    Scope after the review (2026-10-16 review): this item is the
    umbrella of the upload bookkeeping tables in the control database
    (child a, Core 1) and the spool and unit-state code (child b).
    Prerequisites: DB.62, API.29.

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
  Split (2026-10-16, Boss 1:07 PM EDT: look at each item for smaller
  chunks and shared work): this item is now the umbrella and is done
  when its subitems are: API.30, ADM.95. API.3 is marked 'Parent' and
  its prerequisites are API.10, 11, 12, 13, 14, 8 and ADM.13, but its
  own Done holds two deliverables none of those builds: the client
  driver (`generate.py --remote URL --api-key KEY` for the same commands
  as a local run, which checks the version, downloads, generates,
  uploads in batches, resumes and reports) and the remote run showing up
  in the job tree and job page like a server-side run.
  Scope after the review (2026-10-16 review): the umbrella of its seven
  subitems and two children: `generate.py --remote URL --api-key KEY`
  for the same commands as a local run (child a), and a remote run shown
  as a job in the job tree and on the job page with its progress and
  report (child b).

- [ ] **API.4 API compatibility data in the docs**
  Boss (2026-10-01 19:28Z): "let's make API compatibility data and put
  that in the docs as another TODO item to add but not do yet." Done:
  the API has a version of its own, separate from the release number,
  and the docs carry a compatibility table: each API version, the
  release that introduced it, the routes and payload formats it
  changed, and which client versions (`generate.py` remote mode,
  API.3, and the web pages' `lib/apiclient.py`) can talk to which
  server versions. The table is updated in the same PR as any API
  change. Question (decided below): is the API version bumped by hand,
  or by a
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
  Decided (2026-10-10): Boss (19:12Z): the API version is bumped by
  hand, with a CI check against the OpenAPI file.
  Prerequisites: API.26.

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
  generating, not after. Built on API.4's compatibility data. Question
  (decided below): how many older versions the server keeps accepting.
  Research (2026-10-09, api-design-standards.md): answer "how many older
  versions": the current major, plus the previous major for one release
  cycle after a major bump, for upload payloads only (no read
  converters). Add `GET /api/version` (public, database-free, returns
  `api_version`, `min_client`, `max_client`, `release`, `version_key`,
  `galaxy_schema`) and the `PlanetGen-API-Version` header both ways;
  `generate.py --remote` checks before reading config. A refusal is `400
  api_version_unsupported` (or `410` past a sunset), not `426`; the
  check does not depend on `/api/health`. Add the version-key comparison
  (design doc 4.5). Decided (Boss, 2026-10-10 18:44Z, defaults approved;
  defaults taken): a remote
  client may differ in release within the accepted range, with a
  warning; a version-key mismatch (OS, architecture, Python) warns and
  does not refuse.
  Decided (2026-10-10): Boss (19:12Z): keep the current major plus the
  previous major for one release cycle (upload payloads only).

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
  Prerequisites: DB.58, USR.26.

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
  Question (decided below): which corrections the server makes on its
  own, and
  which reject the unit so the client regenerates it.
  Research (2026-10-09, api-design-standards.md): adopt the auto-correct
  versus reject rule (design doc 7.4) as the default answer: correct
  only a clashing name and pure derived columns; reject everything else;
  the server returns a rename map; a rejection does not release the
  claim.
  Decided (2026-10-10): Boss (19:16Z, Accept all): the server auto-fixes
  only a clashing name and purely derived columns, and rejects
  everything else so the client regenerates it.
  Narrowed (2026-10-16, Boss 1:33 PM EDT): the checks, the
  one-transaction finalize and the limited auto-fix are the unit
  finalizer ADM.82, built in 8.6 for change plans. This item is what
  remains: the batch and unit sync with checksums, the staged units, and
  calling the finalizer when a unit is sealed.
  Prerequisites: ADM.82.

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
  Narrowed (2026-10-16, Boss 1:33 PM EDT: build once what 8.6 and 8.9
  both need): the recipe schema, its validation and its runner are
  API.28, built in 8.6 with the regeneration specifications. This item
  is what remains: the route that accepts a recipe, the job it runs as,
  and the 400 and 422 JSON error bodies.
  Prerequisites: API.28.

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
  Needs the region object (GEN.302, built in 8.6): a recipe applied to a
  region uses the same shapes as the admin region tools.
  Prerequisites: API.18, GEN.302.

- [ ] **API.24 GET /api/sectors/<id>/gravity: the per-zone gravity grid as JSON**
  Done: the endpoint returns a sector's grid (GEN.199) for one mode
  (pull, well depth or tidal strength), as a compact array with the zone
  size, the origin, the unit and the Lagrange or saddle points, with the
  usual API scopes and logging (API.15), queued with the 202 behaviour
  when the grid has to be built; documented in the API docs and the
  compatibility data.
  Prerequisite: GEN.199. Related: MAP.167.

- [ ] **API.26 Bump API_VERSION for the shell, band and wedge sector address**
  Found by Quality after the switch (PR #1607). `API_VERSION` is still
  3, but the switch changed the sector address parameters of the API
  from ring, layer and slot to shell, band and wedge, which is a
  breaking change (`src/tests/test_api_version.py` says to bump it and
  add its line to the change list in `docs/api.md`). Owner: Core
  (2026-10-15). Done: `API_VERSION` is 4, `docs/api.md` lists the
  change, and the version header tests pass. Open question for Boss
  (default: bump to 4; there is no backward compatibility and Boss is
  the only user, so no 3 route is kept). Not a bug.
  Prerequisites: none. Related: API.4.

- [ ] **API.27 A read-only profile lookup for any address in the reach, so an uncharted sector can be previewed**
  Step 8 of the whole-sphere train. A route returns the profile of an
  address (expected systems, mix, gas, dark matter, kind factors)
  computed from the address alone, with nothing stored and no sector row
  needed, for the sector page (UX.113), the neighbourhood cost line and
  any client. It follows the API standards (error bodies, versioning,
  the per-key log of API.15). Done: the route answers in microseconds
  for any address inside the reach, refuses one beyond it with the
  grid's own message, is documented in docs/api.md and has tests. Owner:
  Interface asks Core for the profile function (GEN.235). Decided (Boss,
  2026-10-15): the research note's recommendation stands: reach 262 kpc
  (65,535 shells), 96-bit object IDs, one rule set everywhere with
  per-sector parameters, nothing generated unless asked, and an asked
  sector always generated and saved.
  Built-on note (Core 1, 2026-10-16, 3:17 AM EDT, GEN.235 merged as
  #1784): call galaxy.profile.sector_profile(address, skeleton) for the
  profile and GalaxyBounds.describe_miss for the "outside the outline"
  line.
  Owner: Core 4, "Core 4 (2026-10-16)" (Boss approved the lane at 4:31
  PM EDT, 2026-10-16; queue in the lane report, queue position 2).
  Prerequisite: GEN.235. Related: UX.113, API.15.
  Design: [docs/design/intergalactic-and-halo-contents.md](design/intergalactic-and-halo-contents.md)

- [ ] **API.28 Recipe schema and runner: a recipe fixes, ranges or randomises any aspect, and a library turns it into objects**
  Boss (2026-10-16, 1:07 PM EDT), headline of 8.6 (ADM.57): the admin
  has complete control; the world adjusts to what the admin wants, after
  the admin approves the generation plan. Shared with 8.9 (generate
  anywhere; Boss, 1:33 PM EDT: build once what both need). A recipe
  (API.18) and a regeneration specification (ADM.70) are the same thing:
  for each aspect of an object, a fixed value, a range or "draw it".
  Build the schema once as data on top of the aspect catalogue (ADM.62)
  and the keep-or-regenerate marks (ADM.64): a recipe for a system,
  planet, moon or phenomenon with nested children, validation that
  separates nonsense (envelope errors) from rejected values (kind and
  physics checks), the size caps of the API design note, and a runner
  that draws the object from it with the object's own seed stream. No
  route, no job, no HTTP in this item. Done:
  `planetgen/generation/recipe.py` with the schema, the validator and
  the runner, a round trip from marks to a recipe and back, and tests
  for fixed, ranged and random fields.
  Prerequisites: ADM.62, ADM.64. Related: API.18, API.19, ADM.70,
  ADM.57.
  Design: [docs/design/edit-everything.md](design/edit-everything.md)

- [ ] **API.29 Upload spool and unit states: received batches on disk, units assembled at seal, finalize as a queue job**
  Split from API.11 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work). Add the `uploads.dir` spool (gzip
  batches named by run and seq, bounded by the staging cap), the unit
  state machine (open, sealed, verified, finalized or rejected) over
  `upload_units`, assembly of a unit's JSON at seal, deletion of a
  unit's batches after it finalizes with its manifest kept, the
  incomplete flag until a whole unit has arrived, and verify and
  finalize of a sealed unit as an RQ job (200 inside the short wait,
  else 202 with a job id). It calls the finalizer through API.8. Done:
  the spool and the state machine with tests for each transition, a
  restart that keeps staged units, and the spool cap.
  Prerequisites: DB.62, API.10. Related: API.11.

- [ ] **API.30 generate.py --remote: the client driver for sector and galaxy runs**
  Split from API.3 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work). Add the remote mode to the CLI for
  the same options as a local run: check the API version first (API.5),
  download the state (API.12), generate with all workers into the outbox
  (API.13), upload the batches (API.14) with the backoff rules of
  api-design-standards.md section 6.7 (429 and 503 wait, 413 halves the
  batch, 409 on a claimed sector skips it), seal units, resume an
  interrupted run from the outbox and the server's list of held batches,
  poll for the result and print a report of what was added, skipped and
  renamed. First version: the `sector` and `galaxy` commands; the
  commands whose result is not a star system or a sector (`plan`, the
  bright-star scatter and backfill, `population`) refuse with a clear
  message until their upload units are defined (open question for Boss,
  default: sector and galaxy first). Done: a remote run against a test
  server fills the same sectors a local run does, survives a killed
  client and a restarted server, and a refused unit is reported.
  Prerequisites: API.5, API.12, API.13, API.14, API.8. Related: API.3.

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
  Part of, or linked from, the job page (ADM.10, PR #294). Question
  (decided below): should
  an upload with no contact for a long time be flagged as stale on the
  page?
  Research (2026-10-09, api-design-standards.md): the stale flag
  defaults to 3 days without contact and only flags; the page reads
  `upload_runs`, `upload_units` and `upload_claims`.
  Decided (2026-10-10): Boss (19:16Z, Accept all): flag an upload that
  has been silent for 3 days; flag only, never auto-clear.

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
  Moved to phase 3, Move and orbit, before ADM.71, and built on the
  shell and the plan engine: the action becomes the "add" plan operation
  that goes through approval; ADM.71 takes the placement functions from
  here and narrows to the move operation and the drag (2026-10-16
  review).
  Prerequisites: ADM.92, ADM.58, ADM.65, ADM.66, ADM.82.

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
  Prerequisites: MAP.120, ADM.77, ADM.78, ADM.79, ADM.80, ADM.81.

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
  Moved to phase 3, Move and orbit, before ADM.72 (2026-10-16 review):
  it builds the orbit-edit operation once (state vector to elements in
  any frame the position object supports, the escape and swept-path
  collision validation, the recompute, the unbound report) as a plan
  operation; ADM.72 adds the element form and the moons following on the
  same operation.
  Prerequisites: GEN.109, ADM.92, ADM.58, ADM.65, ADM.66.

- [ ] **ADM.43 A full configuration page under Admin**
  Boss (GitHub issue [#515](https://github.com/dwhagar/planetGen/issues/515), 2026-10-08 00:39Z): "A page under admin is
  needed to configure everything EXCEPT for the database information.
  All other settings should be accessed and changeable from there."
  Moved from 8.2 to 8.1 phase 1 (Boss, 2026-10-16, 1:43 PM EDT: "that
  needs to get bumped up because we're starting to desperately need
  that"), and broken into ADM.83 (the validated writer with backups and
  audit rows), ADM.84 (the page, the diff before save and the restart
  banner) and ADM.85 (backups, restore, import and export). Done when
  those three are: an Admin page generated from ADM.42's model lists
  every option with its help text, edits and validates it, saves the
  user's own `config.json`, and says which changes need a restart. The
  Owner tier (the
  recent-password prompt) is ADM.86 in 8.2, because it needs the Owner
  account of USR.2.
  Research (2026-10-09, settings-seo-and-accounts.md): re-scope to three
  tiers (admin, Owner, read-only host-level). `jobs.python`, `jobs.dir`,
  log paths, `secret_key`, `proxy_fix`, `admin_cookie_insecure`,
  `redis.url` and `api_base_url` are read-only in the browser: they are
  code-execution or request-forgery paths. Add: a diff before save,
  timestamped backups with restore, audit rows, per-field errors, a
  pending-restart banner from a startup snapshot, import and export, a
  recent-password prompt for the Owner tier. `set-permissions.sh` and
  the installers need no settings folder: the one file is `config.json`,
  already owned by the web user.
  Decided (Boss, 2026-10-10 18:44Z, defaults approved; default): show
  everything, edit in the browser only
  for the admin and Owner tiers.
  Read side built with ADM.42 (PR #912): `config.json` and a
  `settings.json` overlay are read (the overlay goes with ADM.96),
  nothing writes either yet; options marked `x-editable` false are
  command line only; `x-min_role` owner marks Owner-tier options;
  `util/appconfig.py` is gone (log paths are in `util/logpaths.py`);
  `docs/config.md` has a generated table (`python -m planetgen.cli.config
  docs`).
  Backups (Boss, 2026-10-16, 1:59 and 2:00 PM EDT) are part of this
  page: ADM.87 is the one archive mechanism (one bz2 tarball per backup;
  layout fixed for the settings, control and parity parts) and ADM.85 is
  the screen. Boss (2:00 PM EDT): the control database and parity parts,
  their restores and the repair from backups are targeted for 9.0
  (umbrella ADM.91); the settings page and the settings backups stay in
  8.1.
  One file (Boss, 2:25 and 2:26 PM EDT): all configuration lives in the
  user's own copy of `config.json.example`, named `config.json`; the
  settings chain edits only that copy, never the example, and updates
  never overwrite it (ADM.96). Secrets in it can be replaced but never
  read back, and diffs, audit rows, logs and exports show "(changed)".
  Owner: Settings page lane, "Settings page (2026-10-16)" (Boss,
  2026-10-16, 2:45 PM EDT: it finishes the whole settings chain before
  it goes idle).
  Prerequisites: ADM.85.

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
  Decided (Boss, 2026-10-10 18:44Z, defaults approved; default):
  generated system and sector pages
  are `noindex,follow`, top pages only in the sitemap, switchable with
  `seo.detail_pages`.
  Moved from 8.2 to 8.1 phase 1 with ADM.43, in the same thread right
  after it (Boss, 2026-10-16, 1:43 PM EDT).
  Its options are stored in the user's `config.json`, the one
  configuration file (ADM.96).
  Owner: Settings page lane, "Settings page (2026-10-16)" (Boss,
  2026-10-16, 2:45 PM EDT: it finishes the whole settings chain before
  it goes idle).
  Prerequisites: ADM.43.

- [ ] **ADM.52 A star system editor: plan a system as a 2D diagram or from the command line**
  GitHub issue [#1228](https://github.com/dwhagar/planetGen/issues/1228)
  (Boss, 2026-10-13 02:27Z): "We need a star system editor that the user
  if admin can plan the system graphically or via the command line.
  Graphically will be the order of objects as a 2D diagram, nice and
  simple. Drag and drop radial distance from the star. Positions will be
  randomized." Done: an admin editor shows a system as a simple 2D
  diagram of its objects in order, lets the admin drag each object's
  radial distance, then saves the system with random positions on the
  orbits; the same plan can be written as a command line recipe. Owner:
  Foundations lane 2, after the fly-through map.
  Prerequisites: ADM.92, ADM.58, ADM.66, ADM.82, API.28. Related:
  ADM.32, ADM.35, ADM.36, API.18.
  Moved to phase 2, Approval and history, between ADM.66 and ADM.67
  (2026-10-16 review): the drag edit saves through a plan
  (set-orbit-distance operations with random positions, through
  `save_system_edits`, which ADM.82 wraps), and the command line writes
  the same plan as a recipe (API.28). If it grows past a day, split the
  command line recipe into its own item after the editor.

- [ ] **ADM.55 Admin usability audit: walk every admin page and task, rank what is hard to use, and propose one layout**
  Boss (2026-10-16, 12:52 AM EDT): "also we need to run an interface
  update on everything once this is done so all the controls work
  properly, then run the ADM update chain since our admin pages need
  desperately to be streamlined and fixed, they aren't very usable."
  Done: walk every admin page and task (the admin start page, stats and
  its tables, generation and its jobs, sectors and system generator, the
  queue and its confirmations, facility types, maintenance and the edit
  actions on object pages) as an administrator would, one real task at a
  time (add a sector, regenerate, backfill, erase, edit an object, read
  a job, change a setting); for each record clicks needed, controls that
  fail or mislead, wording, missing feedback, missing or doubled pages,
  and the order of the menu; write a design note in docs/design (with
  its README row in the same PR) that ranks the defects and proposes one
  admin layout: one shell with one navigation, the same table, form,
  confirmation and job-log pattern everywhere. Queued after UX.115; the
  audit starts the ADM chain.
  Prerequisite: UX.115. Related: ADM.56, ADM.43, ADM.35, UX.102.

- [ ] **ADM.56 Streamline and fix the admin pages from the audit: one admin shell, one navigation and one pattern for tables, forms, confirmations and job logs**
  Boss (2026-10-16, 12:52 AM EDT): "also we need to run an interface
  update on everything once this is done so all the controls work
  properly, then run the ADM update chain since our admin pages need
  desperately to be streamlined and fixed, they aren't very usable."
  Done: build the layout the audit (ADM.55) proposes and fix every
  defect it ranks above its line, so the admin pages are usable end to
  end: one shell and navigation, the same table, form, confirmation and
  job-log pattern on every page, clear feedback and error messages, and
  nothing left that is doubled or dead. The new admin items of the chain
  (ADM.13, ADM.35, ADM.32, ADM.36, ADM.52) are built inside this shell;
  the settings page (ADM.43, ADM.44) is built earlier, in 8.1, in the
  admin layout of the day, and this item moves it into the shell.
  Prerequisites: ADM.92, ADM.93. Related: ADM.43, ADM.44, ADM.13,
  ADM.35, UX.102, UX.49.
  Split (2026-10-16, Boss 1:07 PM EDT: look at each item for smaller
  chunks and shared work): this item is now the umbrella and is done
  when its subitems are: ADM.92, ADM.93. ADM.56 is 'build the audit's
  layout and fix every defect it ranks above its line': its size is
  unknown until ADM.55 reports, and one PR cannot restyle and repair
  every admin page.
  Scope after the review (2026-10-16 review): this item is the umbrella
  of the shell (child a: one layout, one navigation and the shared
  table, form, confirmation and job-log components) and moving every
  existing admin page onto it with the fixes the audit ranks above its
  line (child b). The new admin items ADM.32, ADM.36, ADM.52 and POP.13
  build on the shell.

- [ ] **ADM.57 Headline of 8.6: edit everything, the admin can change anything and the world adjusts after the admin approves a generation plan**
  Boss (2026-10-16, 1:07 PM EDT): "I want the 'edit everything' to be a
  headline feature for 8.6, and by that I mean that the admin has
  complete control they can edit and change anything how they want. They
  can move stars in a sector, change orbits in a system, change
  densities for a sector or region, regenerate any star they want to any
  specification they want. The world will adjust to what the admin
  wants. That means if I want to change an A to an M planet the system
  will figure out how to rearrange the star system to make that happen
  but it will get approval for the new generation plan first. In regen
  screens the admin should also be able to select aspects they want to
  keep and what they want regenerated." Umbrella of 8.6, broken into
  items in the 8.6 version file (design:
  [docs/design/edit-everything.md](design/edit-everything.md)); it is done
  when its items are.
  Boss (2026-10-16, 1:49 PM EDT): an admin approves their own plan at
  the time, so an admin confirms rather than requests; there is no
  second approver and the plan preview stays. The roles of 10.0 (USR.25)
  add requests from users and moderators on the same plans.
  Headline: 8.6 (edit everything: the admin can change anything and the
  world adjusts, after approving a plan).
  Prerequisites: ADM.32, ADM.35, ADM.36, ADM.52, MAP.119, MAP.120,
  POP.13, ADM.68, ADM.71, ADM.72, ADM.74, ADM.77, ADM.78, ADM.79,
  ADM.80, ADM.81, TEST.175, TEST.176, TEST.177, DOC.54.

- [ ] **ADM.58 Change plan engine: operations, the effect diff and the plan states, as a library**
  Boss (2026-10-16, 1:07 PM EDT), headline of 8.6 (ADM.57): the admin
  has complete control; the world adjusts to what the admin wants, after
  the admin approves the generation plan. A plan is an ordered list of
  operations (move an object, set an aspect, regenerate with a keep set,
  add, remove) on stored objects. Build `planetgen/admin/plans.py`: the
  operation types, validation of each against the stored object (Hill
  spheres, spacing, bounds), the computed effect diff (every row that
  would change, with before and after, and every object touched
  indirectly), the state machine (draft, requested, proposed, approved,
  applied, undone, stale, rejected, withdrawn: `requested` is a plan
  someone below the approving role submitted, used by 10.0 and built
  here so no state is added later; an admin's own plan goes from draft
  to proposed and the same admin confirms it, with no second approver,
  Boss 2026-10-16, 1:49 PM EDT; every transition records who and in
  which role) and the stale check (a plan is stale when any
  object it touches changed after the plan was proposed). An operation
  serialises in the versioned envelope the remote upload uses
  (`StarSystem.to_dict()`, API.13), so a plan and an upload describe a
  system the same way. Pure Python over loaded objects, no database, so
  it can be tested on its own. Done: the library with unit tests for
  every operation kind, the diff and the stale check.
  Staleness (2026-10-16 review): "changed" means changed by an edit (the
  journal row or an edit-version stamp per object), not by the nightly
  orbital update; a pending move or orbit plan re-bases its before
  values after an update has run and the approval screen shows the
  drift. Add "orbital update movement" to the tests.
  Waits for 8.6 (Boss, 2026-10-16, 3:00 PM EDT: a global top blocker
  that is not started before its version); no lane yet.
  Prerequisites: none. Related: ADM.57, ADM.6, ADM.7, ADM.8.
  Design: [docs/design/edit-everything.md](design/edit-everything.md)

- [ ] **ADM.59 Change plan store: save, load and list plans, and mark them stale**
  Boss (2026-10-16, 1:07 PM EDT), headline of 8.6 (ADM.57): the admin
  has complete control; the world adjusts to what the admin wants, after
  the admin approves the generation plan. Read and write plans through
  the tables of DB.36: `save_plan`, `load_plan`, `list_plans` (by state,
  by object, by admin), and the stale marking when an object a proposed
  plan touches is edited by anything else. Done: the store with tests
  against MySQL, and a plan that survives a restart with its computed
  diff intact so the admin approves exactly what they were shown.
  Prerequisites: DB.36, ADM.58. Related: ADM.57.
  Design: [docs/design/edit-everything.md](design/edit-everything.md)

- [ ] **ADM.60 Apply, undo and audit a plan: one job, one activity-log entry per operation**
  Boss (2026-10-16, 1:07 PM EDT), headline of 8.6 (ADM.57): the admin
  has complete control; the world adjusts to what the admin wants, after
  the admin approves the generation plan. Each unit is written through
  the unit finalizer (ADM.82). An approved plan is applied, under the
  approver's identity, as a visible queue job (the journal and the
  activity log name the requester and the approver) (the Generate page
  job log): operations run in order inside
  one transaction per sector, each writes its journal rows and an
  activity-log entry, and a failure rolls the whole plan back and leaves
  it proposed with the reason. An applied plan can be undone as long as
  nothing newer touched the same objects: undo is the inverse plan built
  from the stored before values. Apply refuses while the update guard is
  on (maintenance, a running update, or the orbital update). Done:
  apply, undo and audit with tests, including a failed apply that
  changes nothing and an undo that restores every stored row.
  Prerequisites: ADM.59. Related: ADM.57, OPS.41, ADM.82.
  Design: [docs/design/edit-everything.md](design/edit-everything.md)

- [ ] **ADM.61 Edits survive repair and regeneration: replay the journal after a sector is rebuilt from its seed or from parity**
  Boss (2026-10-16, 1:07 PM EDT), headline of 8.6 (ADM.57): the admin
  has complete control; the world adjusts to what the admin wants, after
  the admin approves the generation plan. A sector rebuilt from its seed
  (DB.17) or repaired from a parity file (DB.9) would lose every admin
  edit, and GEN.57 says a sector depends only on the seed, the version
  and its address. Keep that rule: the generated sector is the base and
  the edit journal is replayed over it, in order, so the rebuilt sector
  equals the edited one. A journal row that no longer applies (the base
  changed with a generator epoch) is listed for the admin, not dropped.
  Done: replay after both rebuild paths, a test that a sector with edits
  rebuilds to the same stored rows, and the list of rows that cannot
  replay.
  Prerequisites: ADM.60, DB.17, DB.9, GEN.57. Related: DB.9, DB.17,
  GEN.57, GEN.59.
  Design: [docs/design/edit-everything.md](design/edit-everything.md)

- [ ] **ADM.62 Aspect catalogue for systems, stars, planets, moons and belts: every stored field and every generation input, with how they depend on each other**
  Boss (2026-10-16, 1:07 PM EDT), headline of 8.6 (ADM.57): the admin
  has complete control; the world adjusts to what the admin wants, after
  the admin approves the generation plan. An aspect is one thing about
  an object that can be kept or regenerated: for a star its class, mass,
  luminosity, age, metallicity, name, position and companions; for a
  planet its class, mass, radius, orbit, rotation, atmosphere,
  composition, moons, name and facilities; for a system the set of
  bodies and their order. Write the catalogue as data in
  `planetgen/admin/aspects.py`: each aspect, the stored columns it
  covers, the generation input that draws it, whether it is derived from
  other aspects (a planet's radius follows its class and mass) and what
  depends on it. Done: the catalogue with a test that every stored
  column of every body table belongs to exactly one aspect.
  Prerequisites: none. Related: ADM.57, ADM.35, ADM.6, ADM.7.
  Design: [docs/design/edit-everything.md](design/edit-everything.md)

- [ ] **ADM.63 Aspect catalogue for phenomena: black holes, neutron stars, nebulae, remnants, rogue planets, comets, asteroid fields and quasars**
  Boss (2026-10-16, 1:07 PM EDT), headline of 8.6 (ADM.57): the admin
  has complete control; the world adjusts to what the admin wants, after
  the admin approves the generation plan. The same catalogue for the
  phenomenon tables (`PHENOMENON_TABLES` in planetgen/db/edits.py),
  including the scatter rows of the whole-sphere build: each aspect, its
  columns, its generation input and its dependencies. Done: the
  catalogue extended, with the same every-column-belongs-to-one-aspect
  test.
  Prerequisites: ADM.62. Related: ADM.57.
  Design: [docs/design/edit-everything.md](design/edit-everything.md)

- [ ] **ADM.64 Keep-or-regenerate resolution: from the admin's marks to the exact set of aspects redrawn**
  Boss (2026-10-16, 1:07 PM EDT), headline of 8.6 (ADM.57): the admin
  has complete control; the world adjusts to what the admin wants, after
  the admin approves the generation plan. Every aspect is marked keep,
  regenerate (the default draw), set to a value, or draw within a range.
  Resolve the marks into what the generator redraws: a regenerated
  aspect redraws its derived aspects unless they are kept (kept radius
  with a regenerated class is a conflict to report), a kept aspect is
  fed to the generator as a fixed input, and the draw uses the object's
  own seed stream so the same marks give the same result. Done:
  `resolve(marks)` returning the redraw set, the fixed inputs and any
  conflicts, with tests for every dependency edge in the catalogue;
  `regenerate_planet` in planetgen/admin/edits.py gains the fixed-input
  path.
  Prerequisites: ADM.62. Related: ADM.57, ADM.35.
  Design: [docs/design/edit-everything.md](design/edit-everything.md)

- [ ] **ADM.65 Locks: keep marks that stay on an object or one of its aspects, shown on the pages and honoured by every plan**
  Boss (2026-10-16, 1:07 PM EDT), headline of 8.6 (ADM.57): the admin
  has complete control; the world adjusts to what the admin wants, after
  the admin approves the generation plan. An admin can lock an object or
  one aspect of it (a star's class, a planet's orbit). Every plan treats
  a locked aspect as kept and refuses an operation that would change it
  unless the admin lifts the lock inside the plan; the system page and
  the maps mark locked objects. The plan engine takes a lock predicate
  and this item registers it, so every plan refuses an operation on a
  locked aspect. Done: lock and unlock from the object page, the
  registered predicate, the marks, and a test that the plan engine
  refuses a change to a locked aspect (the regeneration, rearrangement
  and density items each test their own).
  Prerequisites: DB.37, ADM.62, ADM.58, MAP.217. Related: ADM.57.
  Design: [docs/design/edit-everything.md](design/edit-everything.md)

- [ ] **ADM.66 Plan approval screen: every operation with its before and after, approve, reject one, or cancel**
  Boss (2026-10-16, 1:07 PM EDT), headline of 8.6 (ADM.57): the admin
  has complete control; the world adjusts to what the admin wants, after
  the admin approves the generation plan. When a plan is proposed the
  admin sees it before anything changes: a summary, then each operation
  with the values before and after, the objects touched indirectly,
  warnings (a locked aspect, an unbound orbit), the number of rows and
  sectors it changes, and buttons to approve the whole plan, reject
  single operations (which re-validates the rest), or cancel. Done: the
  screen, per-operation reject with re-validation, approval stored with
  who and when, and a browser test. An admin confirms their own plan on
  this screen (no second approver, Boss 2026-10-16, 1:49 PM EDT); when a
  plan came from a request the same screen shows the requester, their
  note and the role that must approve, which is how 10.0 reuses it
  (USR.18).
  Prerequisites: ADM.60. Related: ADM.57.
  Design: [docs/design/edit-everything.md](design/edit-everything.md)

- [ ] **ADM.67 Plan diagram: the system before and after as the 2D diagram, and the sector with moved stars marked**
  Boss (2026-10-16, 1:07 PM EDT), headline of 8.6 (ADM.57): the admin
  has complete control; the world adjusts to what the admin wants, after
  the admin approves the generation plan. Reuse the system diagram of
  ADM.52 to draw a plan: the system before and after side by side with
  changed, added and removed bodies marked, and the sector view with
  moved systems drawn from old to new position. Done: both views on the
  approval screen with a browser test.
  Prerequisites: ADM.52, ADM.66, MAP.217. Related: ADM.57, ADM.52.
  Design: [docs/design/edit-everything.md](design/edit-everything.md)

- [ ] **ADM.68 Plan history: every plan with who asked, approved, applied and undid it, and an undo button**
  Boss (2026-10-16, 1:07 PM EDT), headline of 8.6 (ADM.57): the admin
  has complete control; the world adjusts to what the admin wants, after
  the admin approves the generation plan. An admin page lists plans by
  state, object, requester and approver, opens any plan as it was
  approved, and undoes
  an applied plan when nothing newer touched its objects (the undo
  itself is a plan that goes through approval). Done: the page, undo
  through approval, and a browser test.
  Prerequisites: ADM.66, ADM.60. Related: ADM.57.
  Design: [docs/design/edit-everything.md](design/edit-everything.md)

- [ ] **ADM.69 Keep-or-regenerate matrix: the one screen component every regen screen uses**
  Boss (2026-10-16, 1:07 PM EDT), headline of 8.6 (ADM.57): the admin
  has complete control; the world adjusts to what the admin wants, after
  the admin approves the generation plan. In every regen screen the
  admin picks what to keep and what to regenerate: a list of the
  object's aspects, each a keep, regenerate or set-to-value control,
  with presets (regenerate everything, keep everything except one
  aspect, keep only the name), locked aspects shown as kept and fixed,
  and the conflicts from the resolution step shown as they appear. Build
  it as one component for the object, system, sector and region screens.
  Done: the component, used by one screen to prove it, with a browser
  test.
  Prerequisites: ADM.64, ADM.65. Related: ADM.57, ADM.35.
  Design: [docs/design/edit-everything.md](design/edit-everything.md)

- [ ] **ADM.70 Specification form: "regenerate to this specification" built from the aspect catalogue and the system generator options**
  Boss (2026-10-16, 1:07 PM EDT), headline of 8.6 (ADM.57): the admin
  has complete control; the world adjusts to what the admin wants, after
  the admin approves the generation plan. The input half of a regen
  screen: the values the admin may set (star class and mass, age, planet
  count, a habitable world, a binary, the +/- options of the system
  command and the system file keys) generated from the catalogue so a
  new aspect appears in every form without code. Done: the form, a round
  trip to and from a recipe (API.28), and a browser test.
  Prerequisites: ADM.62, ADM.69. Related: ADM.57, ADM.35, API.28.
  Design: [docs/design/edit-everything.md](design/edit-everything.md)

- [ ] **ADM.71 Move a star system or phenomenon within its sector or to another sector, as a plan**
  Boss (2026-10-16, 1:07 PM EDT), headline of 8.6 (ADM.57): the admin
  has complete control; the world adjusts to what the admin wants, after
  the admin approves the generation plan. The admin moves a star system
  or a phenomenon by new coordinates (the three coordinate forms of
  ADM.32, relative to the sector or the galaxy centre), by dragging on
  the Sector Map, or to the emptiest spot: the plan checks Hill spheres
  and cover lists, changes the sector membership and the neighbours'
  cover lists, and never moves anything else unless the admin approves
  it. Done: the plan operation, the Sector Map drag, and tests for a
  move inside a sector, across a sector edge and a refused move.
  Scope after the review (2026-10-16 review): the move operation and the
  drag; the placement functions (emptiest point, the three coordinate
  forms, random outside every Hill sphere) come from ADM.32.
  Prerequisites: ADM.32, ADM.58, ADM.65. Related: ADM.57, ADM.32.
  Design: [docs/design/edit-everything.md](design/edit-everything.md)

- [ ] **ADM.72 Edit a body's orbit: distance, eccentricity, inclination and angles, with its moons following, as a plan**
  Boss (2026-10-16, 1:07 PM EDT), headline of 8.6 (ADM.57): the admin
  has complete control; the world adjusts to what the admin wants, after
  the admin approves the generation plan. The admin sets a planet's or
  belt's orbital elements (distance, eccentricity, inclination, node,
  argument, phase) from its page or the system diagram of ADM.52; moons
  follow their planet, the plan checks Hill spacing and stability and
  recomputes `next_update_due`, and the velocity-vector edit of ADM.36
  is the same operation in another form. Done: the plan operation, the
  page form, and tests for a planet with moons and for a refused
  crossing orbit.
  Prerequisites: ADM.36, ADM.52, ADM.58, ADM.65, GEN.289. Related:
  ADM.57,
  ADM.36, ADM.52.
  Scope after the review (2026-10-16 review): the element form and the
  moons following; the operation itself is ADM.36's.
  Design: [docs/design/edit-everything.md](design/edit-everything.md)

- [ ] **ADM.73 Density overrides for any region: store them, show them, and apply them to later fills**
  Boss (2026-10-16, 1:07 PM EDT), headline of 8.6 (ADM.57): the admin
  has complete control; the world adjusts to what the admin wants, after
  the admin approves the generation plan. Extend the density edit of
  MAP.119 from the slab or block in view to any region: the admin picks
  a region (the tool of MAP.212), sets a multiplier from 0.1 to 10,
  previews the effect, confirms, and can undo; the override table of
  DB.38 holds it and every later fill reads it. Done: the edit, preview
  and undo, the later-fill read, the logged entry, and tests that two
  overlapping regions combine in a stated order.
  Prerequisites: DB.38, GEN.302, MAP.212, ADM.94. Related: ADM.57,
  MAP.119.
  Scope after the review (2026-10-16 review): the region tool hand-off
  for any shape, the preview over an arbitrary region and the test that
  two overlapping regions combine; the store and the fill read are the
  second child of MAP.119 (which already stores
  `region_density_overrides`, DB.38).
  Design: [docs/design/edit-everything.md](design/edit-everything.md)

- [ ] **ADM.74 Adjust the world to a density change: a plan that adds or removes systems so a region matches its new density**
  Boss (2026-10-16, 1:07 PM EDT), headline of 8.6 (ADM.57): the admin
  has complete control; the world adjusts to what the admin wants, after
  the admin approves the generation plan. A density change on a region
  that already holds systems proposes the plan that makes the stored
  world match it: the systems to add (from the bright-star and fill
  backfill, with their draws) or to remove (the fewest, never a locked,
  named or populated one without the admin's approval), with counts per
  sector. Done: the plan generator, the approval of it through the
  screen of ADM.66, and tests that a raised and a lowered density each
  reach the target within a stated tolerance.
  Prerequisites: ADM.73, ADM.58, ADM.65. Related: ADM.57.
  Design: [docs/design/edit-everything.md](design/edit-everything.md)

- [ ] **ADM.75 Rearrange a system around a change to a planet, moon or belt: the plan that fewest changes make valid again**
  Boss (2026-10-16, 1:07 PM EDT), headline of 8.6 (ADM.57): the admin
  has complete control; the world adjusts to what the admin wants, after
  the admin approves the generation plan. The admin changes one body (a
  class A planet to class M, a mass, a distance) and the system works
  out how to rearrange itself: which other bodies move, resize,
  regenerate, appear or go, chosen to change as little as possible,
  honouring keep marks and locks, validated for spacing, Hill spheres,
  the habitable zone and moons. The result is a plan, never applied
  silently. Generalises `validation.stabilize_star_system`, which moves
  bodies and never removes them. Done: `rearrange(system, change,
  marks)` returning a plan, tests for the A to M case, a refused case
  that explains why, and a determinism test.
  Prerequisites: ADM.58, ADM.64, ADM.65. Related: ADM.57, ADM.6.
  Design: [docs/design/edit-everything.md](design/edit-everything.md)

- [ ] **ADM.76 Rearrange a system around a change to a star: class, mass, luminosity or companions**
  Boss (2026-10-16, 1:07 PM EDT), headline of 8.6 (ADM.57): the admin
  has complete control; the world adjusts to what the admin wants, after
  the admin approves the generation plan. The same planner for a star
  change: a new class or mass moves the habitable zone and the frost
  line, so bodies rescale or move with it; a companion added or removed
  in a multi-star system changes which orbits are stable. Extends
  `change_star` in planetgen/admin/edits.py from a single star to every
  kind of system. Done: the star cases of `rearrange`, tests for a
  single star, a binary and a three-star system, and a refused case that
  explains why.
  Prerequisites: ADM.75. Related: ADM.57, ADM.7.
  Design: [docs/design/edit-everything.md](design/edit-everything.md)

- [ ] **ADM.77 Regenerate a planet, moon or belt to a specification, with a keep set, rearranging the system if it must**
  Boss (2026-10-16, 1:07 PM EDT), headline of 8.6 (ADM.57): the admin
  has complete control; the world adjusts to what the admin wants, after
  the admin approves the generation plan. The regen screen of a body:
  the specification form and the matrix; the request becomes a plan (the
  body redrawn with the kept aspects fixed, plus whatever the rearranger
  changes) that goes through approval. Done: the screen for planets,
  moons and belts, the plan, and a browser test of the A to M example
  from the headline.
  Prerequisites: ADM.70, ADM.69, ADM.75, ADM.66. Related: ADM.57,
  ADM.35, ADM.8.
  Design: [docs/design/edit-everything.md](design/edit-everything.md)

- [ ] **ADM.78 Regenerate a star or a whole system to a specification, with a keep set**
  Boss (2026-10-16, 1:07 PM EDT), headline of 8.6 (ADM.57): the admin
  has complete control; the world adjusts to what the admin wants, after
  the admin approves the generation plan. The regen screen of a star or
  system: choose the specification, keep the aspects wanted (the name,
  the position, a planet), regenerate the rest; the rearranger handles a
  changed star. Done: the screen, the plan through approval, and a
  browser test that a kept name, position and chosen planet survive.
  Prerequisites: ADM.70, ADM.69, ADM.76, ADM.66. Related: ADM.57,
  ADM.35, ADM.7.
  Design: [docs/design/edit-everything.md](design/edit-everything.md)

- [ ] **ADM.79 Regenerate a phenomenon to a specification, with a keep set**
  Boss (2026-10-16, 1:07 PM EDT), headline of 8.6 (ADM.57): the admin
  has complete control; the world adjusts to what the admin wants, after
  the admin approves the generation plan. The regen screen for every
  phenomenon kind, from the phenomenon catalogue; a regenerated nebula
  keeps its stars, or redraws them, as marked. Done: the screen for each
  kind and a browser test for two of them.
  Prerequisites: ADM.63, ADM.70, ADM.69, ADM.66. Related: ADM.57,
  ADM.35.
  Design: [docs/design/edit-everything.md](design/edit-everything.md)

- [ ] **ADM.80 Regenerate a sector or a region keeping chosen aspects: locked objects, positions, names**
  Boss (2026-10-16, 1:07 PM EDT), headline of 8.6 (ADM.57): the admin
  has complete control; the world adjusts to what the admin wants, after
  the admin approves the generation plan. From the Galaxy Map,
  regenerate a sector, block, slab, wedge or any region (MAP.212) with
  the matrix: keep the stars' positions and names and redraw their
  planets, keep named systems, redraw the rest. The plan counts what
  changes per sector and goes through approval; a long run applies as a
  job. Done: the action, the plan, tests that locked and kept aspects
  survive, and a browser test.
  Prerequisites: ADM.69, MAP.212, ADM.74, ADM.66, ADM.65. Related:
  ADM.57, ADM.35, MAP.120.
  Design: [docs/design/edit-everything.md](design/edit-everything.md)

- [ ] **ADM.81 Edit any stored field of any object as a plan: one generic editor driven by the aspect catalogue**
  Boss (2026-10-16, 1:07 PM EDT), headline of 8.6 (ADM.57): the admin
  has complete control; the world adjusts to what the admin wants, after
  the admin approves the generation plan. The detail screen of every
  object kind offers each stored field as an edit through a plan, with
  the same approval, history and undo; the field list comes from the
  catalogue, so a new column appears without new screens. Done: the
  editor for every kind with a browser test for a star, a planet, a
  nebula and a black hole.
  Prerequisites: ADM.63, ADM.66, ADM.69. Related: ADM.57, ADM.35.
  Design: [docs/design/edit-everything.md](design/edit-everything.md)

- [ ] **ADM.82 Unit finalizer: validate a system or sector, write it in one transaction, correct only a clashing name and derived columns, or reject with reasons**
  Boss (2026-10-16, 1:07 PM EDT), headline of 8.6 (ADM.57): the admin
  has complete control; the world adjusts to what the admin wants, after
  the admin approves the generation plan. Shared with 8.9 (generate
  anywhere; Boss, 1:33 PM EDT: build once what both need). A change plan
  (ADM.60) and a remote upload (API.8) both end the same way: a unit (a
  star system or a sector) is checked and written whole or not at all.
  Build it once in `planetgen/admin/finalize.py`: the physics and shape
  checks of the existing validator (ADM.5,
  `planetgen/generation/validation.py`), the check that the unit stays
  inside its sector, the one-transaction write through `db/edits.py`,
  the auto-fix limited to a clashing name and purely derived columns
  (the rule Boss accepted for uploads, 2026-10-10), a rename map, and a
  rejection that lists every reason and changes nothing. Done: the
  library with tests for a valid unit, a name clash, a derived column
  fix and each kind of rejection, and `save_system_edits` goes through
  it.
  Prerequisites: none. Related: ADM.57, ADM.60, API.8.
  Design: [docs/design/edit-everything.md](design/edit-everything.md)

- [ ] **ADM.85 Settings backups, restore, import and export on the settings page**
  Boss (2026-10-16, 1:43 PM EDT): "I know there's a config page set
  somewhere, that needs to get bumped up because we're starting to
  desperately need that." Part of the configuration page of ADM.43
  (GitHub issue
  [#515](https://github.com/dwhagar/planetGen/issues/515)), moved from
  8.2 to 8.1 phase 1. The settings page is also the one screen for
  backups (ADM.87): it lists the archives (time, who, the parts each
  holds and their sizes), takes "Back up now" (every part available),
  downloads one (the Owner only, because an archive holds the whole
  `config.json` with its secrets; until the Owner exists in 8.2 the
  archives are owner-only files and an admin can export a version
  without secrets), shows the diff of an
  archive's settings part against the current settings, and restores the
  whole archive or one part after a confirmation and an audit row (on
  the web a restore of the settings part puts back only the options an
  admin may change; only the command line restores the whole file); the
  control database part (ADM.88, restore ADM.89) and the parity part
  (DB.40, restore DB.41) arrive in 9.0 and show up in the same list and
  restore through the same screen, with the guard of ADM.89, so this
  screen is built once for every part. Export and import are for the
  settings part alone: export
  downloads the editable options as JSON without secrets (a secret shows
  only "(changed)") or command-line-only options;
  import reads
  such a file through the same validation and diff as an edit. Done: the
  list, back up now, download, whole and one-part restore, export and
  import with a browser test and a test that an import with one bad
  value changes nothing.
  Owner: Settings page lane, "Settings page (2026-10-16)" (Boss,
  2026-10-16, 2:45 PM EDT: it finishes the whole settings chain before
  it goes idle).
  Prerequisites: ADM.84. Related: ADM.43.

- [ ] **ADM.86 Owner-tier settings: Owner-only options behind a recent-password prompt**
  Boss (2026-10-16, 1:43 PM EDT): "I know there's a config page set
  somewhere, that needs to get bumped up because we're starting to
  desperately need that." Part of the configuration page of ADM.43
  (GitHub issue
  [#515](https://github.com/dwhagar/planetGen/issues/515)), moved from
  8.2 to 8.1 phase 1. Options marked `x-min_role` owner (SMTP,
  `web.base_url`, the From address) are editable only by the Owner, and
  only after the Owner re-entered their password within the last 10
  minutes; admins see them disabled with the reason. This needs the
  Owner account that USR.2 creates, so it belongs to 8.2 and not to the
  umbrella ADM.43. Done: the prompt, the tier check in the writer, and
  tests that an admin and a stale-password Owner are both refused.
  Prerequisites: USR.26. Related: ADM.43, USR.3.

- [ ] **ADM.88 Control database backup part: a dump of the planetgen_control schema (admin logins, API keys, run records) in the backup archive**
  Boss (2026-10-16, 1:59 PM and 2:00 PM EDT): "let's add backups of the
  control DB and backups of the parity data for DB repair." and "All
  those should be bz2 tarballed together." Boss (2:00 PM EDT): "we also
  have to build a way to use that data so build that out and let's make
  this targetted for v9.0 so incorporate the chain for these (restore /
  repair using backups) with that due date. Settings page is v8.1".
  Target 9.0 (the backups and repair chain, umbrella ADM.91); the
  settings page and the settings backups stay in 8.1 (ADM.83 to ADM.85).
  All backup pieces are parts of one `.tar.bz2` archive (ADM.87);
  retention is `settings.backups_kept` (default 20, 1 to 200). The
  control database (`control_database`, default `planetgen_control`)
  holds the admin identities, API keys, naming key and run records;
  nothing backs it up. This item is the part that dumps it into the
  archive: `mysqldump --single-transaction` (the credentials go through
  a temporary 0600 option file, never the command line; no replication
  flags, which MySQL 8.4 and MariaDB 11.4 do not share), leaving out the
  sessions table, written to `control/`. The archive holds password
  hashes and key hashes, so only admins may download it, and from 8.2
  only the Owner (ADM.86). Done: the part registered with ADM.87, a
  missing `mysqldump` reported with the install step, a test that
  sessions are not in the archive, and a test that the dump loads into a
  scratch schema (the restore itself is ADM.89).
  Prerequisites: ADM.87. Related: ADM.43, ADM.85, USR.2.

- [ ] **ADM.89 Control database restore: load a backed-up planetgen_control into a scratch schema, check it, then replace the live schema, from the command and from the settings page**
  Boss (2:00 PM EDT): "we also have to build a way to use that data so
  build that out and let's make this targetted for v9.0 so incorporate
  the chain for these (restore / repair using backups) with that due
  date. Settings page is v8.1". Target 9.0 (the backups and repair
  chain, umbrella ADM.91); the settings page and the settings backups
  stay in 8.1 (ADM.83 to ADM.85). All backup pieces are parts of one
  `.tar.bz2` archive (ADM.87); retention is `settings.backups_kept`
  (default 20, 1 to 200). Using the control part of an archive: the
  restore loads the dump into a scratch schema, checks the tables, row
  counts and the schema version against the live code (an archive from a
  newer version is refused, an older one is migrated by the normal
  update), then replaces the control schema in one step, behind a
  confirmation page that names what will be lost (every change since the
  archive) and behind the update guard (OPS.41, OPS.53, OPS.54: refused
  while maintenance, an update or a run is active). It is logged, ends
  every session and makes the admin sign in again. `planetgen.cli.backup
  restore --part control ARCHIVE` and the "Restore" action of the backup
  list (ADM.85) use it, so the control part restores like the settings
  part. Done: the restore with tests for a good archive, a refused
  archive from a newer version, a refused restore while the guard is on,
  and a restore that fails in the scratch schema and leaves the live
  schema untouched.
  Prerequisites: ADM.88, ADM.85, OPS.41. Related: ADM.85, ADM.43,
  ADM.91.

- [ ] **ADM.90 Database health and repair page: run the check, see damaged and unprotected sectors, repair from the live parity or from an archive, restore parts, as a queue job**
  Boss (2:00 PM EDT): "we also have to build a way to use that data so
  build that out and let's make this targetted for v9.0 so incorporate
  the chain for these (restore / repair using backups) with that due
  date. Settings page is v8.1". Target 9.0 (the backups and repair
  chain, umbrella ADM.91); the settings page and the settings backups
  stay in 8.1 (ADM.83 to ADM.85). All backup pieces are parts of one
  `.tar.bz2` archive (ADM.87); retention is `settings.backups_kept`
  (default 20, 1 to 200). One admin page for everything in the chain:
  the DB.8 check result (damaged, unprotected and healthy sector counts,
  when it last ran), a "Run check" button that runs as a queue job, the
  list of damaged sectors with the source each repair would use (live
  parity, an archive chosen from the backup list, or the seed of DB.17),
  a "Repair" action that starts `repair-db` as a queue job with a
  progress bar and the report, and the restore actions of the backup
  list for the control and parity parts, all behind the update guard and
  a confirmation that names what will change. Done: the page and a
  browser test that runs a check on a damaged fixture galaxy, repairs it
  from an archive and shows a passing check.
  Prerequisites: DB.42, ADM.89, DB.41. Related: ADM.85, DB.8, ADM.91.

- [ ] **ADM.91 Backups and repair, target 9.0: one bz2 archive of the settings, the control database and the parity data, and the ways to use it (restore the control database, restore parity, repair the galaxy, the repair page)**
  Boss (2026-10-16, 1:59 PM and 2:00 PM EDT): "let's add backups of the
  control DB and backups of the parity data for DB repair." and "All
  those should be bz2 tarballed together." Boss (2:00 PM EDT): "we also
  have to build a way to use that data so build that out and let's make
  this targetted for v9.0 so incorporate the chain for these (restore /
  repair using backups) with that due date. Settings page is v8.1".
  Target 9.0 (the backups and repair chain, umbrella ADM.91); the
  settings page and the settings backups stay in 8.1 (ADM.83 to ADM.85).
  All backup pieces are parts of one `.tar.bz2` archive (ADM.87);
  retention is `settings.backups_kept` (default 20, 1 to 200). Umbrella
  of the chain; it is done when its items are. The archive mechanism and
  the settings part are in 8.1 (ADM.87, with ADM.83 to ADM.85); the
  control database and parity parts, their restores, the repair from an
  archive, the repair page, the end-to-end test and the guide are in
  9.0.
  Prerequisites: ADM.88, ADM.89, DB.40, DB.41, DB.42, ADM.90, TEST.180,
  DOC.56, OPS.55. Related: ADM.43, ADM.87, DB.9.

- [ ] **ADM.92 Admin shell: one layout, one navigation and the shared table, form, confirmation and job-log components**
  Split from ADM.56 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work). The admin pages have no common frame,
  so every new admin page invents its own table, form, confirmation and
  job log. ADM.55's design note proposes one layout; this item builds
  only the frame and the shared parts, not the page-by-page conversion.
  Build the admin base template and the one navigation (menu order from
  the audit), and the four shared parts as template macros or partials
  plus their CSS: the table, the form (Shoelace fields, UX.49), the
  confirmation page and the job log (the UX.102 logging standard).
  Convert the admin start page as the worked example and write the usage
  rules in the design note. Done: the shell and the four parts exist and
  are documented, the admin start page uses them, every item of the ADM
  chain can be built on them, and a browser test opens the converted
  page in both themes with no console error.
  Prerequisites: ADM.55. Related: ADM.56.

- [ ] **ADM.93 Move every existing admin page onto the shell and fix the defects the audit ranks above its line**
  Split from ADM.56 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work). With the shell in place, convert the
  remaining admin pages (stats and its tables, generation and its jobs,
  sectors and the system generator, the queue and its confirmations,
  facility types, maintenance, the edit actions on object pages) and fix
  each defect ADM.55 ranks above its line: failing or misleading
  controls, wording, missing feedback, doubled or dead pages, the menu
  order. ADM.55 files one item per page group when it reports, and this
  item is the umbrella over those; if it finds fewer than three groups
  it stays one item. Done: every admin page uses the shell and the
  shared parts, no defect above the audit's line is left open, and
  nothing is doubled or dead.
  Prerequisites: ADM.92. Related: ADM.56.

- [ ] **ADM.94 Density override store and the fill read: set, undo and combine overrides for any region, and have later fills use them**
  Split from MAP.119 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). Write the library under
  planetgen/admin: `set_override(region, multiplier, admin, reason)`,
  `undo_override(id)` (append-only, the undo is a reverse row),
  `multiplier_for(sector)` combining overlapping regions in a stated
  order (narrowest wins), the range 0.1 to 10 on a log scale (a typed
  value beyond it with a warning), an activity-log entry for each
  change, and the read in the sector fill so a later fill multiplies its
  `--density` and records `density_multiplier_used` and `override_id`.
  Only unfilled sectors are affected; the function reports how many
  sectors in a region are filled and unaffected. A block or slab is a
  region object (GEN.302), so there is no separate block-or-slab path.
  Done: the library over `region_density_overrides` (DB.38), the fill
  read, and tests for set, undo, two overlapping regions in a stated
  order, a filled sector left alone and the logged entry.
  Prerequisites: DB.38, DB.60, GEN.302. Related: MAP.119.

- [ ] **ADM.95 A remote run is a job: it shows in the job tree and the job page with its progress and report**
  Split from API.3 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work). An upload run appears in the job tree
  and on the job page (ADM.10, ADM.12) like a server-side run: who
  started it, its claimed sectors, units sealed, verified and finalized,
  rejected units with reasons, last contact, and the final report; the
  log follows the UX.102 standard. Done: a run created through the
  upload routes appears in the tree, the page updates as batches arrive,
  and a browser test checks the page for a running, a finished and a
  stalled run (the stalled one links to the incomplete uploads page,
  ADM.13).
  Prerequisites: API.10, API.14, API.8. Related: API.3.

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
  when Secure. Ideally before USR.2 so new accounts start on it. Decided
  (Boss, 2026-10-10 18:44Z, defaults approved; default yes): move from
  PBKDF2-600k to Argon2id (a
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
  Decided (Test Suite Rebuild thread, 2026-10-14 14:04Z, list in
  /mnt/project-files/test-suite-rebuild/test-fix-decisions.md): fold
  into the rebuild, no separate hunt; stage 1 keeps each failure's
  traceback and re-runs it alone, stage 3 gives each worker its own
  database; closed after 5 clean full runs after stage 3, or fixed from
  the captured traceback on recurrence. Owner: Quality (2026-10-15, took over from the Test Suite Rebuild thread).
  Recurrence (Quality, 2026-10-15, 4:56 PM EDT run): the test hung, not
  just failed, and looks like GEN.232's uncapped path integration of a
  fast star. Core fixes GEN.232; Quality re-checks this test after it
  and closes this item on 5 clean full runs as decided.
  Quality confirms (2026-10-15, combined run, MySQL 8.0.46): the test
  runs past 900 s and fails alone on main, in db/sector_paths.py settle
  and physics/sector_path.py integrate_path, _simplify and hermite; it
  points at GEN.232.
  Intermittent (Quality, 2026-10-15): it passed in the 6:13 PM EDT run
  on the same code that timed out alone at 4:20 PM EDT, so it is not
  deterministic; GEN.232 merged in #1660, so Quality re-checks it in the
  next full run.
  Owner: Core 1, in its bug phase after the global top chain (Boss,
  2026-10-16, 2:43 PM EDT: "Core 1 will start working on bugs after";
  order set at 3:00 PM EDT: non-test bugs, then the database-version
  test failures, then the failures that hit full runs, then the one-off
  load flakes). Fix the cause; no skipping or quarantine.
  Prerequisites: none.

- [ ] **TEST.116 test_bughunt_end_to_end stores M2V or M6V where it expects K2V on the MySQL 8.4 and MariaDB legs once (bug)**
  Reported by Bugfixes lane 1 (2026-10-10 03:57Z): the K2V check in
  test_bughunt_end_to_end failed once on MySQL 8.4 and once on MariaDB
  in run 37882710453 (the stored star was M2V or M6V instead of K2V); it
  did not reproduce on current main in 83 local runs. Done: the cause is
  found (a seed or ordering dependence in the test, or a real bug in
  star-type selection) and the test is made robust without loosening it,
  or the product bug is fixed. Decided (Boss, 2026-10-10 18:44Z,
  defaults approved; default: leave
  open until it recurs, then investigate with the failing run's data).
  Decided (Test Suite Rebuild thread, 2026-10-14 14:04Z, list in
  /mnt/project-files/test-suite-rebuild/test-fix-decisions.md): stays
  open for the first recurrence (with the captured seed): the test pins
  the star type, so a different stored type may be a real pin bug in the
  product, not test luck. Owner: Quality (2026-10-15, took over from the Test Suite Rebuild thread).
  Owner: Core 1, in its bug phase after the global top chain (Boss,
  2026-10-16, 2:43 PM EDT: "Core 1 will start working on bugs after";
  order set at 3:00 PM EDT: non-test bugs, then the database-version
  test failures, then the failures that hit full runs, then the one-off
  load flakes). Fix the cause; no skipping or quarantine.
  Prerequisites: none. Related: TEST.111, TEST.71.

- [ ] **TEST.126 The heavy WebGL browser-map tests time out when four workers run them together (bug)**
  Two System Map browser tests that failed under load (seen by Foundations lane 2 and the rebuild thread) were fixed on 2026-10-14 (#1367).
  Seen again on 2026-10-15 (Bugfixes lane 2, NAV.47 suite, load 8.9, call 40.6 s on
  gw0): src/tests/test_web_browser_controls.py::test_every_plain_button_does_something[web.sector];
  passed alone, main not checked.
  Seen again on 2026-10-14 12:44Z (Bugfixes lane 2, PERF.81 suite, 49 min):
  test_galaxy_map_free_camera_from_an_arc_down, ::test_galaxy_map_buttons (a
  playwright error), test_bookmarks_keep_a_course_pick_on_the_galaxy_map and
  ::test_galaxy_map_hover_lights_every_choice_while_picking_a_course, all
  passing on re-run. The three TEST.125 moon tests still fail alone and on main.
  Reported by Foundations lane 2 (2026-10-10 20:50Z): in a full suite
  under `python3 -m pytest -n 4 -q tests` (about 15,450 tests, 65-69
  min, 4 cores, MariaDB and Redis local) these pass alone and fail
  together, a pattern of SwiftShader WebGL tests timing out when four
  workers run them at once. MAP.148 run:
  tests/test_web_browser_maps.py::test_galaxy_map_free_camera_from_an_arc_down,
  ::test_galaxy_map_buttons (a playwright error),
  ::test_nav_ends_picked_from_bookmarks_on_every_page. MAP.149 run: the
  same three, plus
  tests/test_web_browser_fixture_maps.py::test_bookmarks_keep_a_course_pick_on_the_galaxy_map
  and
  ::test_galaxy_map_hover_lights_every_choice_while_picking_a_course,
  and
  tests/test_web_system_nav_row.py::test_the_navigate_menu_goes_where_the_wide_buttons_go[From
  here-a.nav-wide >> nth=0] (Bugfixes lane 1 saw the nav-row test the
  same way). A fixture-maps-only run under -n 4 (22 min) had 12 failures
  that all passed serially:
  test_bookmarks_keep_a_course_pick_on_the_galaxy_map,
  test_galaxy_map_hover_lights_every_choice_while_picking_a_course,
  test_galaxy_map_bookmark_star_saves_the_stage_and_the_menu_opens_it,
  test_galaxy_map_does_not_prefetch_on_a_data_saving_connection,
  test_galaxy_breadcrumb_always_matches_the_view_its_url_names,
  test_galaxy_map_draws_point_phenomena_that_link_to_their_pages,
  test_the_map_help_is_a_menu_item_that_opens_a_dialog,
  test_galaxy_map_point_phenomena_hover_and_offer_nav_links,
  test_galaxy_map_opens_zoomed_in_on_the_charted_space,
  test_galaxy_map_pick_mode_clicks_down_through_every_stage,
  test_the_slab_rail_shows_only_while_the_stage_has_slab_buttons and
  test_the_map_card_goes_wide_but_the_title_keeps_the_page_edge. Full
  tracebacks were not kept. Done: the cause is found (a time budget the
  tests do not scale with load, too many WebGL pages at once per worker,
  or a real wait that is missing) and the tests are made reliable under
  -n 4 (a shared limit on concurrent WebGL pages, waits on page state
  not fixed times, or a serial marker for the heaviest), without
  skipping or loosening any; a full -n 4 run shows none of them failing.
  Decided (Test Suite Rebuild thread, 2026-10-14 14:04Z, list in
  /mnt/project-files/test-suite-rebuild/test-fix-decisions.md): fold
  into the rebuild, no separate hunt; stage 1 runs the browser tests as
  a separate phase on two workers; closed after 5 clean full runs, and
  stage 4 gives any still timing out a state-based wait. Owner: Quality (2026-10-15, took over from the Test Suite Rebuild thread).
  Note (2026-10-10): Bugfixes lane 2 (21:06Z) adds
  tests/test_web_browser_fixture_maps.py::test_galaxy_map_color_by_switches_the_fill_with_a_legend_and_keeps_it_in_the_url:
  failed once in the NAV.9 run under load, passed alone.
  Sighting (Bugfixes lane 2, 2026-10-14, full suite of 50 min 43 s with
  UX.95): tests/test_web_browser_maps.py::test_galaxy_map_buttons (a
  playwright error) and ::test_nav_ends_picked_from_bookmarks_on_every_page
  failed once and pass alone, not checked on main. Same family; the
  rebuild's captured traceback and serial re-run tell load from real.
  Owner: Core 1, in its bug phase after the global top chain (Boss,
  2026-10-16, 2:43 PM EDT: "Core 1 will start working on bugs after";
  order set at 3:00 PM EDT: non-test bugs, then the database-version
  test failures, then the failures that hit full runs, then the one-off
  load flakes). Fix the cause; no skipping or quarantine.
  Prerequisites: none. Related: TEST.111, TEST.121, TEST.122.

- [ ] **TEST.127 test_nebula_shape_endpoint_serves_a_mesh failed once in a full parallel run (bug)**
  Reported by Foundations lane 2 (2026-10-10 20:50Z):
  tests/test_api.py::test_nebula_shape_endpoint_serves_a_mesh failed in
  a full `-n 4` run (the MAP.148 run, 19:27Z) with an assertion starting
  'assert 1' (the rest was cut from the log) and passed on a serial
  rerun; main was not checked. Bugfixes lane 1 has not seen it. Done:
  the next failure's full traceback is captured (or a loaded run
  reproduces it), the cause is found (a shared nebula or sector fixture,
  an id or ordering dependence, or a real mesh bug) and the test is made
  reliable without loosening it. Decided (Test Suite Rebuild thread, 2026-10-14 14:04Z, list in
  /mnt/project-files/test-suite-rebuild/test-fix-decisions.md): fold
  into the rebuild, no separate hunt; stage 1 keeps each failure's
  traceback and re-runs it alone, stage 3 gives each worker its own
  database; closed after 5 clean full runs after stage 3, or fixed from
  the captured traceback on recurrence. Owner: Quality (2026-10-15, took over from the Test Suite Rebuild thread).
  Owner: Core 1, in its bug phase after the global top chain (Boss,
  2026-10-16, 2:43 PM EDT: "Core 1 will start working on bugs after";
  order set at 3:00 PM EDT: non-test bugs, then the database-version
  test failures, then the failures that hit full runs, then the one-off
  load flakes). Fix the cause; no skipping or quarantine.
  Prerequisites: none. Related: TEST.111, TEST.116.

- [ ] **TEST.128 test_spatial_position_db and test_web_db_fields fail under a parallel run on one MySQL (bug)**
  Reported by Bugfixes lane 1 (2026-10-10 20:50Z): intermittent failures
  of test_spatial_position_db and test_web_db_fields (both under tests/)
  when the suite runs in parallel (`-n 4`) on one MySQL; both pass alone
  and none fail on a quiet main. The assertion messages were not saved.
  Done: the next failure's traceback is captured, the cause is found (a
  shared database name or schema, a leftover row, or timing under load;
  compare TEST.111 and TEST.123) and the tests are made reliable without
  skipping or loosening them. Decided (Test Suite Rebuild thread, 2026-10-14 14:04Z, list in
  /mnt/project-files/test-suite-rebuild/test-fix-decisions.md): fold
  into the rebuild, no separate hunt; stage 1 keeps each failure's
  traceback and re-runs it alone, stage 3 gives each worker its own
  database; closed after 5 clean full runs after stage 3, or fixed from
  the captured traceback on recurrence. Owner: Quality (2026-10-15, took over from the Test Suite Rebuild thread).
  Sighting (Foundations lane 1, 2026-10-14, full suite of 39 min 51 s
  with PERF.89): tests/test_work_queue_admin.py::
  test_a_paused_queue_node_stands_by_until_resumed[2] failed once and
  passes alone: the same parallel-run family on the shared database.
  Owner: Core 1, in its bug phase after the global top chain (Boss,
  2026-10-16, 2:43 PM EDT: "Core 1 will start working on bugs after";
  order set at 3:00 PM EDT: non-test bugs, then the database-version
  test failures, then the failures that hit full runs, then the one-off
  load flakes). Fix the cause; no skipping or quarantine.
  Prerequisites: none. Related: TEST.111, TEST.123.

- [ ] **TEST.129 Two orbit-update tests fail on the MySQL 8.4 leg: the stored last_updated_at rounds up (bug)**
  Reported by Foundations lane 1 (2026-10-10 20:55Z), from the CI
  run on main after PR #1161 (MySQL 8.4 leg; not caused by PERF.64):
  test_admin_script_cli.py::test_update_orbits_never_moves_backwards_when_the_clock_goes_back
  and
  test_edge_admin_scripts.py::test_update_orbits_first_run_sets_a_starting_point_then_advances.
  On MySQL 8.4 the stored last_updated_at rounds up, so it reads
  '0.000000 years in the future' and the elapsed time is -1 s; the
  second test expects 'years elapsed' and gets the clock-moved-back
  message. Likely cause: fractional-second rounding of NOW() in the
  orbit update (MariaDB truncates, MySQL 8.4 rounds). Done: the orbit
  update stores and compares a time that both servers treat the same way
  (truncate to the column's precision before comparing, or use the same
  clock on both sides), checked on MySQL 8.4 and MariaDB; a product bug
  if a real run could see a negative elapsed time. Decided (Test Suite Rebuild thread, 2026-10-14 14:04Z, list in
  /mnt/project-files/test-suite-rebuild/test-fix-decisions.md):
  implement, one small PR, number 4 in the thread's order; the stored
  last-updated time is compared the way MySQL 8.4 rounds it; the thread
  installs MySQL 8.4 to verify. Owner: Quality (2026-10-15, took over from the Test Suite Rebuild thread).
  Status (2026-10-09): Fixed by the Test Suite Rebuild thread
  (2026-10-14, #1313, suite green): the stored last-updated time is
  compared the way MySQL 8.4 rounds it. Retire this item once the
  scheduled MySQL 8.4 CI run on main is green.
  Owner: Core 1, in its bug phase after the global top chain (Boss,
  2026-10-16, 2:43 PM EDT: "Core 1 will start working on bugs after";
  order set at 3:00 PM EDT: non-test bugs, then the database-version
  test failures, then the failures that hit full runs, then the one-off
  load flakes). Fix the cause; no skipping or quarantine.
  Prerequisites: none. Related: TEST.116, TEST.126.

- [ ] **TEST.130 test_the_scene_positions_at_the_epoch_match_the_stored_ones fails when a random system has comets (bug)**
  Reported by Foundations lane 1 (2026-10-10 20:55Z), from the CI
  run on main after PR #1161:
  tests/test_body_positions.py::test_the_scene_positions_at_the_epoch_match_the_stored_ones
  fails because a comet's position at its epoch differs from the stored
  one by about 1e10 km, for a random system that has comets. It is
  intermittent because the system is random. Likely a real disagreement
  between the comet path in positions_at and the position stored for the
  comet. Done: the cause is found (the stored position or the scene path
  is wrong for comets, or the test compares different epochs or frames),
  the product or test is fixed so scene and stored positions agree for
  comets, and the test draws a system with a comet every time so it
  cannot pass by luck. Decided (Test Suite Rebuild thread, 2026-10-14 14:04Z, list in
  /mnt/project-files/test-suite-rebuild/test-fix-decisions.md):
  implement, one small PR, number 5 in the thread's order; the test
  draws a system with a comet every time. Owner: Quality (2026-10-15, took over from the Test Suite Rebuild thread).
  Status (2026-10-09): Test Suite Rebuild thread (2026-10-14, #1313):
  not reproducible; the test now builds comets by construction so it
  exercises the comet path every time. Stays open until the cause is
  found or five clean full runs pass.
  Owner: Core 1, in its bug phase after the global top chain (Boss,
  2026-10-16, 2:43 PM EDT: "Core 1 will start working on bugs after";
  order set at 3:00 PM EDT: non-test bugs, then the database-version
  test failures, then the failures that hit full runs, then the one-off
  load flakes). Fix the cause; no skipping or quarantine.
  Prerequisites: none. Related: TEST.123, GEN.122.

- [ ] **TEST.131 Tests for the gravity field evaluator: two bodies, Lagrange points, softening and the rotation curve**
  Done: the evaluator (GEN.198) matches Newton for two bodies to 1e-9
  relative, puts the Sun-Earth L1 and L2 points where the published
  distances are (about 1.5 million km) with a saddle in the tidal
  gradient, stays finite inside a softened source, and reproduces the
  galaxy potential's rotation curve (229.3 km/s at 8.128 kpc); the
  near/far split agrees with the exact sum to a stated tolerance on a
  random sector; the same seed gives the same field twice.
  Decided (Test Suite Rebuild thread, 2026-10-14 14:04Z, list in
  /mnt/project-files/test-suite-rebuild/test-fix-decisions.md): stays
  with whoever builds the gravity evaluator (GEN.198): they write these
  tests in the same PR as a normal gravity test file. Owner: that
  builder.
  Prerequisite: GEN.198. Related: GEN.198, GEN.115.

- [ ] **TEST.134 test_many_processes_starting_at_once_start_one_job fails under a parallel run (bug)**
  Bugfixes lane 1 (2026-10-10, full local suite 19:41Z to 20:54Z, -n 4):
  src/tests/test_web_jobs_damage.py::test_many_processes_starting_at_once_start_one_job
  fails when the suite runs in parallel and passes alone. Done: the
  cause is found (a shared Redis or database state with the other
  workers, or a timing limit) and the test is made to pass under -n 4
  without skipping, disabling or loosening it away from what it checks.
  Decided (Test Suite Rebuild thread, 2026-10-14 14:04Z, list in
  /mnt/project-files/test-suite-rebuild/test-fix-decisions.md): fold
  into the rebuild, no separate hunt; stage 1 keeps each failure's
  traceback and re-runs it alone, stage 3 gives each worker its own
  database; closed after 5 clean full runs after stage 3, or fixed from
  the captured traceback on recurrence. Owner: Quality (2026-10-15, took over from the Test Suite Rebuild thread).
  Owner: Core 1, in its bug phase after the global top chain (Boss,
  2026-10-16, 2:43 PM EDT: "Core 1 will start working on bugs after";
  order set at 3:00 PM EDT: non-test bugs, then the database-version
  test failures, then the failures that hit full runs, then the one-off
  load flakes). Fix the cause; no skipping or quarantine.
  Prerequisites: none. Related: TEST.126, TEST.127, TEST.128.

- [ ] **TEST.135 Big-galaxy query-budget test for the nebula reads**
  Boss (2026-10-10 21:46Z): nebulae must be visible on the Galaxy Map,
  by colour or by combined regions; Phase 2 for the feature, foundations
  in Phase 1. From the research (docs/design/nebula-map-visibility.md
  section 9, PR #1180; report
  /mnt/project-files/research/nebula-map-visibility/report.md): the
  nebula reads (field table, aggregates, tile layer, scatter join) are
  added to the big-galaxy query budget test of TEST.133, so none
  examines more rows than a tile needs. Done: the nebula reads are in
  the test and it passes on 2,000,000 systems. Decided (Test Suite Rebuild thread, 2026-10-14 14:04Z, list in
  /mnt/project-files/test-suite-rebuild/test-fix-decisions.md):
  implement after TEST.133 and DB.24, in the same opt-in `budget` tier;
  Foundations lane 1 drops it from its queue. Owner: Quality (2026-10-15, took over from the Test Suite Rebuild thread).
  Prerequisites: MAP.172, MAP.173, GEN.202. Related: PERF.68.
  Design: [docs/design/nebula-map-visibility.md](design/nebula-map-visibility.md)

- [ ] **TEST.137 test_a_neighborhood_from_the_generate_page[core] fails once under load (bug)**
  Reported by Bugfixes lane 2 (2026-10-10 22:47Z):
  tests/test_neighborhood_edge.py::test_a_neighborhood_from_the_generate_page[core]
  failed once in a full run (started 21:43Z, 57 min) under load and
  passed alone; the traceback was not kept. Done: the next failure's
  traceback is captured (or a loaded run reproduces it), the cause is
  found (a timeout, a shared table between workers, or a real defect in
  the neighborhood path) and the test is made reliable without loosening
  it. Decided (Test Suite Rebuild thread, 2026-10-14 14:04Z, list in
  /mnt/project-files/test-suite-rebuild/test-fix-decisions.md): fold
  into the rebuild, no separate hunt; stage 1 keeps each failure's
  traceback and re-runs it alone, stage 3 gives each worker its own
  database; closed after 5 clean full runs after stage 3, or fixed from
  the captured traceback on recurrence. Owner: Quality (2026-10-15, took over from the Test Suite Rebuild thread).
  Owner: Core 1, in its bug phase after the global top chain (Boss,
  2026-10-16, 2:43 PM EDT: "Core 1 will start working on bugs after";
  order set at 3:00 PM EDT: non-test bugs, then the database-version
  test failures, then the failures that hit full runs, then the one-off
  load flakes). Fix the cause; no skipping or quarantine.
  Prerequisites: none. Related: TEST.111, TEST.128, TEST.132 (done, PR #1578).

- [ ] **TEST.138 test_open_map_menus_hold_no_overlap[web.galaxy-1280] fails once in a full parallel run (bug)**
  Reported by Foundations lane 2 (2026-10-10 23:33Z):
  tests/test_web_browser_controls.py::test_open_map_menus_hold_no_overlap[web.galaxy-1280]
  failed in a full suite run (-n 4, 15,485 passed and 7 failed, 74 min,
  load average up to 14) and passed when re-run alone with 3 others (10
  min); not checked on main. Same load-timeout family as TEST.126. Done:
  the next failure's traceback is captured, the cause is found (a load
  timeout in the page wait, or a real overlap that appears only at that
  moment) and the test is made reliable without loosening it. Owner:
  Bugfixes lane 2, with TEST.126.
  Decided (Test Suite Rebuild thread, 2026-10-14 14:04Z, list in
  /mnt/project-files/test-suite-rebuild/test-fix-decisions.md): fold
  into the rebuild, no separate hunt; same family and fix as TEST.126.
  Owner: Quality (2026-10-15, took over from the Test Suite Rebuild thread).
  Recurrence (Quality, 2026-10-15, 4:56 PM EDT run):
  `test_open_map_menus_hold_no_overlap[web.sector-600]` also failed once
  under load; same family, same fix.
  Same family and fix as TEST.126 (`tests/browser_waits.py`): close with
  it, reopen only on a captured traceback (2026-10-16 review).
  Owner: Core 1, in its bug phase after the global top chain (Boss,
  2026-10-16, 2:43 PM EDT: "Core 1 will start working on bugs after";
  order set at 3:00 PM EDT: non-test bugs, then the database-version
  test failures, then the failures that hit full runs, then the one-off
  load flakes). Fix the cause; no skipping or quarantine.
  Prerequisites: none. Related: TEST.126, TEST.128, TEST.137.

- [ ] **TEST.140 test_a_body_is_reached_by_its_id_and_a_row_number_names_nothing fails once in a full parallel run (bug)**
  Reported by Bugfixes lane 2 (2026-10-14 12:44Z, PERF.81 run, -n 4, 49
  min):
  tests/test_api_object_ids.py::test_a_body_is_reached_by_its_id_and_a_row_number_names_nothing
  failed once in a full run and passed alone; the traceback was not
  kept. Seen again (Foundations lane 1, 2026-10-14 14:04Z full run, 27
  min): the failure showed a facility host id 003 where 009 was
  expected, and the test passed twice alone; unrelated to that lane's
  change. Done: the next failure's traceback is captured, the cause is
  found (a shared table between workers, a timeout, or a real defect in
  the ID lookup) and the test is made reliable without loosening it.
  Decided (Test Suite Rebuild thread, 2026-10-14 14:04Z, list in
  /mnt/project-files/test-suite-rebuild/test-fix-decisions.md): fold
  into the rebuild, no separate hunt; stage 1 keeps each failure's
  traceback and re-runs it alone, stage 3 gives each worker its own
  database; closed after 5 clean full runs after stage 3, or fixed from
  the captured traceback on recurrence. Owner: Quality (2026-10-15, took over from the Test Suite Rebuild thread).
  Owner: Core 1, in its bug phase after the global top chain (Boss,
  2026-10-16, 2:43 PM EDT: "Core 1 will start working on bugs after";
  order set at 3:00 PM EDT: non-test bugs, then the database-version
  test failures, then the failures that hit full runs, then the one-off
  load flakes). Fix the cause; no skipping or quarantine.
  Prerequisites: none. Related: TEST.128, TEST.137, TEST.138.

- [ ] **TEST.141 test_interrupting_a_parallel_galaxy_run_leaves_no_half_written_sector[15-True] fails once in a full parallel run (bug)**
  Reported by Bugfixes lane 2 (2026-10-14 12:44Z, PERF.81 run, -n 4, 49
  min):
  tests/test_work_queue_failures.py::test_interrupting_a_parallel_galaxy_run_leaves_no_half_written_sector[15-True]
  failed once in a full run and passed alone; the failure message was
  not kept. Done: the next failure's message is captured, the cause is
  found (an interrupt landing before a worker starts, a shared Redis
  queue between workers, or a real half-written sector) and the test is
  made reliable without loosening it. Decided (Test Suite Rebuild thread, 2026-10-14 14:04Z, list in
  /mnt/project-files/test-suite-rebuild/test-fix-decisions.md): fold
  into the rebuild, no separate hunt; stage 1 keeps each failure's
  traceback and re-runs it alone, stage 3 gives each worker its own
  database; closed after 5 clean full runs after stage 3, or fixed from
  the captured traceback on recurrence. Owner: Quality (2026-10-15, took over from the Test Suite Rebuild thread).
  Owner: Core 1, in its bug phase after the global top chain (Boss,
  2026-10-16, 2:43 PM EDT: "Core 1 will start working on bugs after";
  order set at 3:00 PM EDT: non-test bugs, then the database-version
  test failures, then the failures that hit full runs, then the one-off
  load flakes). Fix the cause; no skipping or quarantine.
  Prerequisites: none. Related: TEST.128, TEST.132 (done, PR #1578), TEST.134.

- [ ] **TEST.144 test_every_plain_button_does_something[web.sector] failed once in a full parallel run (bug)**
  Reported by Foundations lane 2 (2026-10-14 13:17Z, shared suite -n 4,
  15,573 passed and 9 failed, 54 min):
  tests/test_web_browser_controls.py::test_every_plain_button_does_something[web.sector]
  failed once under load and passes serially and on a re-run; not
  checked on main. Same load-timeout family as TEST.126 and TEST.138.
  Done: the next failure's traceback is captured, the cause is found (a
  load timeout in the page wait, or a button that really does nothing at
  that moment) and the test is made reliable without loosening it.
  Decided (Test Suite Rebuild thread, 2026-10-14 14:04Z, list in
  /mnt/project-files/test-suite-rebuild/test-fix-decisions.md): fold
  into the rebuild, no separate hunt; same family and fix as TEST.126.
  Owner: Quality (2026-10-15, took over from the Test Suite Rebuild thread).
  Same family and fix as TEST.126 (`tests/browser_waits.py`): close with
  it (2026-10-16 review).
  Owner: Core 1, in its bug phase after the global top chain (Boss,
  2026-10-16, 2:43 PM EDT: "Core 1 will start working on bugs after";
  order set at 3:00 PM EDT: non-test bugs, then the database-version
  test failures, then the failures that hit full runs, then the one-off
  load flakes). Fix the cause; no skipping or quarantine.
  Prerequisites: none. Related: TEST.126, TEST.138.

- [ ] **TEST.148 Investigate resuming a test run where it left off when the build is the same**
  Boss (2026-10-15 00:09Z): wants to be able to restart a test run where
  it left off so long as the build is the same, so a run that is stopped
  or times out does not start over (a full run is about 49 minutes).
  Boss also set the standing rule that every test run starts with a
  timeout of at least 90 minutes and is known to be done by its own
  completion notice or final summary line, not by polling (now in
  docs/testing.md). Done when the Test Suite Rebuild thread reports the
  options (for example a record of passed tests keyed by the git commit
  and a hash of the working tree, a resume flag in scripts/test.py that
  skips them, and what must invalidate a resume: any changed file,
  database or Redis state, the Python or library versions) with a
  recommendation, and Boss picks one; the build is then filed as its own
  items. The thread reported (2026-10-15, /mnt/project-files/test-suite-rebuild/resume-a-run-options.md): a per-test progress log, a build id hashed from the tree, uncommitted changes, lock files, test command, Python and database versions, and a `--resume` flag that skips tests that passed; about three small PRs. Options: A, progress log plus an explicit `--resume` (recommended; a resumed run says so in its last line and in the PR); B, automatic resume (not recommended for the merge gate); C, re-run only what failed. Boss picks. Owner: Quality (2026-10-15, took over from the Test Suite Rebuild thread), part of TEST.139 (done, PR #1605), not a lane.
  Decided (Boss, 2026-10-15, tapped Build
  it): option A, a per-test progress log plus an explicit
  `scripts/test.py <command> --resume`, same build id only; not
  automatic and not part of the merge gate. Build items TEST.150 (done, PR #1431),
  TEST.151 (done, PR #1484) and TEST.152 (done, PR #1520), owner the Test Suite Rebuild thread (TEST.133
  done, #1405).
  Quality, 2026-10-15: after a container restart `scripts/test.py
  changed --resume` refused ("packages changed since": the build id
  includes the installed packages), so a run interrupted by a restart in
  the cloud environment cannot be resumed; the browser phase has to be
  rerun alone with `scripts/test.py browser`. Decide whether the build
  id should keep the package list.
  Prerequisites: none. Related: TEST.139 (done, PR #1605).

- [ ] **TEST.155 Wording and index agree: property test**
  Audit: /mnt/project-files/wording/text-output-audit-2026-10-15.md.
  Owner lane: Foundations Wording (2026-10-15), which writes the code
  for this item. A property test over random worlds and the sample
  systems that fails when a domain's colour, score, the band its
  sentence names, the gear and the equipment label disagree; when a body
  has life text but the index gives it no life (check the 130 C case
  first); when an airless body shows Air as anything but Red (shared
  with GEN.215, done in PR #1478); and when any quantity in the describer is formatted
  with a raw :.Nf. Test items are owned by Test Suite Rebuild: that
  thread approves, Foundations Wording writes it.
  Prerequisite: GEN.218.

- [ ] **TEST.156 Rewrite the tests that name ring, layer and slot**
  The 79 affected test files are rewritten as each of the grid items
  lands; a MISSED line is a selection bug. Owner: Quality (2026-10-15, took over from the Test Suite Rebuild thread). Source: item 9 of section 6 of the plan. Boss decided the
  move (2026-10-15 06:39Z): "two meridian lines from pole to pole make a
  wedge (vert), shell (radius), and band (horizontal), together they
  identify a single 3D segment of the sphere". TOP PRIORITY in Phase 1,
  ahead of every lane queue.
  Quality, 2026-10-15 (#1697): the non-map step is done; the map tests
  that name ring, layer or slot remain.
  Design: [docs/design/spherical-sector-coordinates.md](design/spherical-sector-coordinates.md)

- [ ] **TEST.159 Whole-sphere tests: required stars per nebula kind, Wolf-Rayet counts, tiling at the reach, profile and empty-sector determinism, and the pre-placed share per population**
  Step 15 of the whole-sphere train. Tests: no emission, reflection or
  core-collapse remnant in a halo or old cluster cell; a planetary
  nebula always has its white dwarf; a core-collapse remnant has a
  compact core at the stated share; every hosted nebula holds its host
  star and every D or E region its members inside its radius; a
  Wolf-Rayet nebula holds a W star and the Milky Way count is 1,200 +-
  300 with none in old populations; tiling and the object ID round trip
  at the widest shell of the reach; the profile and an empty sector are
  deterministic; the pre-placed share of every population; no double
  count of a cluster star; the King solver against the textbook
  concentrations. Done: each test above exists and runs in the normal
  suite. Owner: Quality (all TEST items are Quality's), written beside
  each build step.
  Split (2026-10-16, Boss 1:07 PM EDT: look at each item for smaller
  chunks and shared work): this item is now the umbrella and is done
  when its subitems are: TEST.184, TEST.185. The item lists tests of two
  kinds of code.
  Prerequisites: TEST.184, TEST.185. Related: TEST.156.
  Design: [docs/design/intergalactic-and-halo-contents.md](design/intergalactic-and-halo-contents.md)

- [ ] **TEST.160 Representation audit: seeded sectors in each region class show every kind at its expected frequency, with systems and required contents checked**
  Step 4 of the phenomena train (catalogue section 7). The registry
  gives the expected count of each kind per sector for each region class
  (core, bar, thin disk arm and interarm, thick disk, halo,
  Magellanic-type satellite, globular cluster core, intergalactic
  floor); fill a fixed set of seeded sectors in each, count by stage,
  and compare: each kind within four sigma of its expectation, any kind
  expected 20 or more times appears at least once, none appears where
  its gate forbids it, the fraction with a system matches the system
  table, and every hosted nebula, planetary nebula, core-collapse
  remnant, D or E region and Wolf-Rayet nebula holds its required
  contents. Also an admin report. Audited, not forced: Boss's reading of
  "all types represented" is pending; default is audited frequency.
  Done: the audit runs in the suite in a few minutes and as a report.
  Owner: Quality (all TEST items are Quality's). Boss's rules
  (2026-10-15 5:39 PM EDT): scatters and backscatters carry no systems,
  planets or moons; nebulae need stars and phenomena, not planets;
  planets and moons only when a sector is filled with everything; a fill
  represents every kind and builds the systems phenomena can have at
  real-world odds.
  Decided (Boss, 2026-10-15 5:51 PM EDT): "all types represented" means
  audited frequency, not forced: accurate probabilities and no forced
  representation in a sector, except a kind that landed there in a
  scatter before generation. The open question is closed; the audit
  compares counts to expectation and never forces a kind in.
  Scope after the review (2026-10-16 review): the representation audit
  on the census harness (the new first item of the Whole-sphere audits
  thread); TEST.161, TEST.162 and TEST.164 add their checks to it, not
  counters of their own.
  Prerequisites: GEN.240, GEN.282, GEN.321,
  TEST.186.
  Related: TEST.159.
  Design: [docs/design/phenomena-catalogue-and-fill.md](design/phenomena-catalogue-and-fill.md)

- [ ] **TEST.161 Fill-stage and system-pass tests: stream independence, same seed same system, and the probability table**
  Tests beside the fill-stages item (catalogue sections 3 and 5): adding
  a kind to a stage moves no draw of another stage; the same object
  gives the same system in any sector fill order; each kind's system
  share matches the table within Poisson bounds; a scatter row carries
  no planet or moon before its sector fills; no room drops the outermost
  body first. Done: each test exists and runs in the normal suite.
  Owner: Quality (all TEST items are Quality's), written beside the
  build step.
  Owner: Core 4, "Core 4 (2026-10-16)" (Boss approved the lane at 4:31
  PM EDT, 2026-10-16; queue in the lane report); held until its blocker
  merges: TEST.186 (Core 3). Core 4 takes it when it unblocks.
  Prerequisites: TEST.186. Related: TEST.159.
  Design: [docs/design/phenomena-catalogue-and-fill.md](design/phenomena-catalogue-and-fill.md)

- [ ] **TEST.162 Multiplicity audit: companion states, row bookkeeping, stability, determinism, no planets in a scatter, fill reuse, motion and the census anchors**
  Step 9 (design note section 9, tests 1 to 9). Companion states within
  four Poisson sigma over a seeded scatter of 1e6 rows; object density
  before and after members; every stored tree passes the stability,
  Kozai and Peters screens; two scatters with the same seed store
  identical members; nothing in `scatter_members` or `scatter_nodes` is
  a planet, moon or belt; a built system equals its stored members and
  the child rows are gone; no member row is written by the nightly
  update; the double pulsar, V404 Cygni and the wide-companion examples;
  black hole pairs between the gravitational-wave stock and the Lamberts
  count. Done: each test exists and runs in the suite, and the
  representation audit report includes them. Owner: Quality (all TEST
  items are Quality's). Boss (2026-10-15 5:51 PM EDT): stellar
  companions of stars and star-like objects, including neutron star and
  black hole pairs, are generated at scatter time because they matter
  for gravity; planets, moons and belts stay at sector fill. Defaults
  taken (reversible): compact pairs promoted to the scatter, black hole
  companion fraction 1.25%, cap 12 stars (7 before the 6:08 PM EDT
  catalogue), hierarchical trees only, a row's members deleted when its
  sector builds it, pair rates from the middle of the evidence.
  The multiplicity checks as census checks on the harness (2026-10-16
  review).
  Owner: Core 4, "Core 4 (2026-10-16)" (Boss approved the lane at 4:31
  PM EDT, 2026-10-16; queue in the lane report); held until its blocker
  merges: GEN.315 and TEST.186 (Core 3). Core 4 takes it when it
  unblocks.
  Prerequisites: GEN.315, TEST.186. Related: TEST.160, TEST.161.
  Design: [docs/design/multi-star-systems.md](design/multi-star-systems.md)

- [ ] **TEST.164 Catalogue type audit: type counts, region gates, the N tail, interaction consistency and group decay**
  Step 20 (design note section 15, tests 10 to 14). Over a seeded fill
  of at least 1e6 systems per region class, each type's count lies
  within four Poisson sigma of its default (contact 0.2% of FGK, Algol
  0.1% to 0.2%, cataclysmic variable 3.4e-5 of systems, AM CVn 3.9e-6,
  double white dwarf 9.5% of white dwarfs, brown dwarf wide companion
  3.2%, runaway 10% to 30% of O stars); no type appears where its gate
  forbids it; the number of systems with N = 8 to 12 matches the tail
  table; every semi-detached pair has a donor at or over its Roche lobe,
  every contact pair both, and no accreting pair has a non-compact
  receiver unless it is a mass gainer; every stored group is younger
  than its decay age. Done: the five tests exist and run in the suite,
  and the representation audit report includes them. The docs of the
  catalogue are in DOC.50. Owner: Quality (2026-10-15). Boss (2026-10-15
  6:08 PM EDT) widened the multi-star plan to "every observed and
  theoretical type of multi-star/star-like-object system", built "into
  our generator throughout to match as close as possible to observed and
  theoretical probabilities (observed are best, theoretical are good
  too)".
  Prerequisites: GEN.269, GEN.270, GEN.271, GEN.272, GEN.275, GEN.327,
  GEN.329, TEST.186. Related: TEST.162, DOC.50, TEST.160.
  The type-count, gate, N-tail, interaction and group-decay checks as
  census checks on the harness (2026-10-16 review).
  Owner: Core 4, "Core 4 (2026-10-16)" (Boss approved the lane at 4:31
  PM EDT, 2026-10-16; queue in the lane report); held until its blocker
  merges: GEN.269-272, 275, 327, 329 and TEST.186. Core 4 takes it when
  it unblocks.
  Design: [docs/design/multi-star-systems.md](design/multi-star-systems.md)

- [ ] **TEST.165 Re-measure the region query on a real generated galaxy, warm and cold**
  Step 6 (design note sections 8 and 9 item 6). The benchmark scripts in
  /mnt/project-files/research/spherical-region-queries/ become
  `scripts/bench_region_queries.py` (seed, centres, the recommended
  statement against the plain box and the sector walk) and run once on a
  generated galaxy of a few million systems, warm and cold (restart with
  a small buffer pool and `O_DIRECT`), with real wide rows and ids from
  several generation runs; the numbers go into the note. Done: the
  script is in the repo, the note carries the real-galaxy figures
  against the research's 42 ms at 50 pc (10 million rows), and any
  figure more than twice the research's is reported to Core. Owner:
  Quality (2026-10-15). Boss (2026-10-15 5:32 PM EDT): "if it's not
  ready in time before the spherical sector layout is, then we'll skip
  it and move to the next sequence until it's finished research and we
  can plan it." The research is finished and merged (#1659), so stage C
  is live.
  Owner: Core 4, "Core 4 (2026-10-16)" (Boss approved the lane at 4:31
  PM EDT, 2026-10-16; queue in the lane report); held until its blocker
  merges: GEN.276 (Core 2). Core 4 takes it when it unblocks.
  Prerequisite: GEN.276. Related: NAV.60, DOC.51.
  Design: [docs/design/spherical-region-queries.md](design/spherical-region-queries.md)

- [ ] **TEST.167 test_random_start_mode_respects_min_start_density hit a 900 s timeout in a full run and passes alone (bug)**
  Reported by Interface (2026-10-15, 5:43 PM EDT gate run, relayed by
  the coordinator): `tests/test_generate*`
  `test_random_start_mode_respects_min_start_density` hit its 900 s
  timeout in a `scripts/test.py changed` run and passes alone. It may be
  the same slow fast-body path as TEST.111 (GEN.232 merged in #1660, so
  check whether it still hangs). Done: the cause is found and the test
  is made reliable without skipping or loosening it, or the product bug
  it hides is fixed. Owner: Quality (2026-10-15).
  Owner: Core 1, in its bug phase after the global top chain (Boss,
  2026-10-16, 2:43 PM EDT: "Core 1 will start working on bugs after";
  order set at 3:00 PM EDT: non-test bugs, then the database-version
  test failures, then the failures that hit full runs, then the one-off
  load flakes). Fix the cause; no skipping or quarantine.
  Prerequisites: none. Related: TEST.111, GEN.232.

- [ ] **TEST.168 Control sweep: a browser test drives every control on every page against a spherical galaxy with 96-bit IDs**
  Boss (2026-10-16, 12:52 AM EDT): "also we need to run an interface
  update on everything once this is done so all the controls work
  properly, then run the ADM update chain since our admin pages need
  desperately to be streamlined and fixed, they aren't very usable."
  Done: a browser test (Playwright, in the existing browser suite) opens
  every page of the audit list in UX.115 on a small generated spherical
  galaxy that has scatter rows, a satellite and a cluster, and presses
  or fills every control once, failing on a console error, a 4xx or 5xx
  response or a control that does nothing; it checks 96-bit IDs in
  inputs and URLs on each page that takes one. Quality runs it after
  UX.115's fixes. Queued after UX.115.
  Widened (2026-10-16 review): owns the shared browser-test watcher and
  sweep harness: `sweep(page, path, fill=...)` returning dead controls,
  console errors and bad responses (beside
  `src/tests/test_web_browser_controls.py`), a hook for a per-page task
  script and for the 96-bit ID input check, and the UX.115 page list.
  Console capture is hand-written today in test_web_a11y.py,
  test_web_browser_sector_open.py, test_web_browser_tables.py and
  test_web_browser_maps.py.
  Prerequisite: UX.115. Related: TEST.163, TEST.156, TEST.169.

- [ ] **TEST.169 Admin page browser tests: every admin page and task works end to end after the streamline**
  Boss (2026-10-16, 12:52 AM EDT): "also we need to run an interface
  update on everything once this is done so all the controls work
  properly, then run the ADM update chain since our admin pages need
  desperately to be streamlined and fixed, they aren't very usable."
  Done: a browser test per admin task of the audit (ADM.55): add a
  sector, regenerate, backfill, erase, edit an object, read a job and
  its log, change a setting, confirm and cancel a queue action; fails on
  a console error, a bad response or a control that does nothing. Runs
  after ADM.56 and again at the end of the chain.
  Scope after the review (2026-10-16 review): the per-task scripts for
  the admin tasks of ADM.55 (add a sector, regenerate, backfill, erase,
  edit an object, read a job and its log, change a setting, confirm and
  cancel a queue action), written on the harness of TEST.168; it adds
  nothing to the failure rules.
  Prerequisites: ADM.56, TEST.168. Related: TEST.168, ADM.55.

- [ ] **TEST.170 test_math_check::test_whole_run_is_fast and the Galaxy Map breadcrumb folding test fail under full-suite load and pass alone (bug)**
  Finding from Core 1 (2026-10-15, PR #1757, kind: load flake, fails on
  main: not seen, how often: once in the 89-minute changed-tests gate
  each): in the gate run, test_math_check::test_whole_run_is_fast failed
  (a timing assertion) and
  test_web_browser_fixture_maps::test_galaxy_breadcrumb_stays_on_one_line_folding_its_middle
  failed once and passed on rerun; both pass alone. Reproducer:
  scripts/test.py changed (the full parallel gate), then run each test
  alone. Done: each test is made robust to load without loosening what
  it proves (the speed test measures against a relative budget or the
  work done instead of a wall-clock limit; the breadcrumb test waits for
  layout to settle before measuring) or its cost is explained, and the
  same gate shows neither red. Owner: unassigned. Related: TEST.111,
  TEST.138, TEST.163, TEST.167.
  Scope after the 2026-10-16 review: only the
  `test_math_check::test_whole_run_is_fast` timing test (a relative
  budget); the breadcrumb folding test joins the TEST.126 family.
  Owner: Core 1, in its bug phase after the global top chain (Boss,
  2026-10-16, 2:43 PM EDT: "Core 1 will start working on bugs after";
  order set at 3:00 PM EDT: non-test bugs, then the database-version
  test failures, then the failures that hit full runs, then the one-off
  load flakes). Fix the cause; no skipping or quarantine.
  Prerequisites: none. Related: TEST.111, TEST.138, TEST.163, TEST.167.

- [ ] **TEST.171 test_gen_resume::test_sectors_a_forced_scatter_skipped_fill_correctly_afterwards failed once under gate load with no scattered system (0 > 0) and passes alone (bug)**
  Finding from Core 1 (2026-10-16, 3:17 AM EDT, GEN.235 gate, kind:
  flaky test, fails on main: not checked, how often: once in the 1:47 AM
  to 3:03 AM EDT gate under load, 9 of 9 passes alone):
  src/tests/test_gen_resume.py:299 asserted some_scattered > 0 and got 0
  > 0. Likely cause (Core 1): an unseeded galaxy seed leaves no bright
  star in the first wedges, so a seed fix is the probable remedy.
  Reproducer: scripts/test.py changed (the full parallel gate), then run
  the test alone in a loop with different seeds. Done: the test pins its
  seed (or builds a guaranteed scatter case) without loosening what it
  proves, and the same gate shows it green. Owner: unassigned.
  Owns pinning the random galaxy seed in the shared `_seed_galaxy` test
  helper (2026-10-16 review); TEST.173 shares the helper.
  Owner: Core 1, in its bug phase after the global top chain (Boss,
  2026-10-16, 2:43 PM EDT: "Core 1 will start working on bugs after";
  order set at 3:00 PM EDT: non-test bugs, then the database-version
  test failures, then the failures that hit full runs, then the one-off
  load flakes). Fix the cause; no skipping or quarantine.
  Prerequisites: none. Related: TEST.170, TEST.111.

- [ ] **TEST.172 test_galaxy_map_forward_to_current_jumps_to_the_newest_view failed once: clicking a Galaxy Map block with 0 sectors charted did nothing, and it passes alone (bug)**
  Finding from Core 3 (2026-10-16, 6:31 AM EDT, kind: browser test,
  intermittent, fails on main: not checked, how often: once in a
  12-minute run of the sample-galaxy tests, 469 passed, and 3 of 3
  passes alone):
  src/tests/test_web_browser_maps.py::test_galaxy_map_forward_to_current_jumps_to_the_newest_view
  clicked the Galaxy Map block "Block 15·0·30, longitude 111.3°-115.1°,
  14578-15550 pc from the core, 0 sectors charted" and nothing happened;
  the hover picked a rim block with no sectors charted. Possibly related
  to MAP.210 (nothing selectable on the Galaxy Map), which was not
  reproduced. Reproducer: run the sample-galaxy browser tests in a loop.
  Done: the cause is found (the test picks a block that cannot be
  entered, or the click really does nothing on an empty rim block) and
  the test or the page is fixed without loosening what the test proves.
  Owner: unassigned.
  Owner: Core 1, in its bug phase after the global top chain (Boss,
  2026-10-16, 2:43 PM EDT: "Core 1 will start working on bugs after";
  order set at 3:00 PM EDT: non-test bugs, then the database-version
  test failures, then the failures that hit full runs, then the one-off
  load flakes). Fix the cause; no skipping or quarantine.
  Second sighting (Core 3, 2026-10-16, 8:34 PM EDT, GEN.333 gate: full
  `changed`, 118 minutes, load 16): failed again under load and passed
  alone (call 28 s on a loaded worker); the same gate failed
  `test_galaxy_map_drill_down_by_clicks` (filed as TEST.202).
  Prerequisites: none. Related: MAP.210, TEST.170, TEST.171.

- [ ] **TEST.173 test_gen_bright_scatter_edges::test_a_fill_after_a_shell_failed_mid_scatter_builds_no_leftover_star failed once under gate load with an empty ring in _nearest (bug)**
  Finding from Core 1 (2026-10-16, 7:18 AM EDT, kind: flaky test, fails
  on main: unknown, how often: once in a loaded scripts/test.py changed
  run from 6:00 to 7:12 AM EDT with 4 workers, 13 of 13 passes alone on
  an idle machine):
  src/tests/test_gen_bright_scatter_edges.py::test_a_fill_after_a_shell_failed_mid_scatter_builds_no_leftover_star
  raised ValueError: min() iterable argument is empty in _nearest(...)
  because the ring_cells intersect given was empty: the sector the
  backfill reached drew no stars (the TEST.132 edge). Not touched by
  GEN.101 or GEN.251. Reproducer: loop the test under load. Done: the
  test pins its seed or builds a guaranteed star case, or _nearest
  handles an empty ring where that is a real edge, without loosening
  what the test proves, and the same gate shows it green. Owner:
  unassigned.
  No separate seed fix (2026-10-16 review): check that its nearest
  backfill cell draws a star under the seed TEST.171 pins.
  Owner: Core 1, in its bug phase after the global top chain (Boss,
  2026-10-16, 2:43 PM EDT: "Core 1 will start working on bugs after";
  order set at 3:00 PM EDT: non-test bugs, then the database-version
  test failures, then the failures that hit full runs, then the one-off
  load flakes). Fix the cause; no skipping or quarantine.
  Prerequisites: none. Related: TEST.132, TEST.171, TEST.172.

- [ ] **TEST.174 test_web_a11y::test_page_is_accessible_and_fits[web.sector-desktop-light] failed once under gate load and passes alone (bug)**
  Finding from Core 2 (2026-10-16, 7:48 AM EDT, kind: browser test, load
  flake, fails on main: not checked, how often: once in a GEN.159 gate,
  load 8.2, passes alone):
  src/tests/test_web_a11y.py::test_page_is_accessible_and_fits[web.sector-desktop-light]
  failed under load. The failure text was not captured. Reproducer:
  scripts/test.py changed (the full parallel gate), then run the test
  alone. Done: the failure text is captured on the next sighting, the
  cause is found (a layout or wait timing under load, or a real overflow
  on the sector page) and the test is made robust without loosening what
  it proves. Owner: unassigned. Related: TEST.170, TEST.172.
  Sighting (Core 3, 2026-10-16, 8:40 AM EDT):
  test_web_a11y.py::test_page_is_accessible_and_fits[web.sector-phone-light]
  failed once under load in the changed run: a request to /galaxy/tiles
  aborted (net::ERR), and it passes alone. Same family, a different
  parametrisation (phone-light).
  Scope after the 2026-10-16 review: capture the failure text on the
  next sighting; the phone-light sighting (an aborted `/galaxy/tiles`
  request) is a page-stat flake of the TEST.126 family.
  Owner: Core 1, in its bug phase after the global top chain (Boss,
  2026-10-16, 2:43 PM EDT: "Core 1 will start working on bugs after";
  order set at 3:00 PM EDT: non-test bugs, then the database-version
  test failures, then the failures that hit full runs, then the one-off
  load flakes). Fix the cause; no skipping or quarantine.
  Prerequisites: none. Related: TEST.170, TEST.172.

- [ ] **TEST.175 Edit everything end to end: the headline scenarios with approval and undo**
  Boss (2026-10-16, 1:07 PM EDT), headline of 8.6 (ADM.57): the admin
  has complete control; the world adjusts to what the admin wants, after
  the admin approves the generation plan. Browser and server tests of
  the headline cases as an admin would run them: change an A planet to
  class M and approve the rearranged system; move a star in a sector;
  change an orbit; change a region's density and approve the added and
  removed systems; regenerate a star to a specification keeping its name
  and a planet; then undo each. Done: one test per case, each asserting
  the stored rows before, after approval and after undo.
  Prerequisites: ADM.77, ADM.78, ADM.71, ADM.72, ADM.74, ADM.80,
  TEST.168.
  Related: ADM.57.
  Design: [docs/design/edit-everything.md](design/edit-everything.md)

- [ ] **TEST.176 Plan engine properties: any sequence of edits leaves every system valid, and undo restores every row**
  Boss (2026-10-16, 1:07 PM EDT), headline of 8.6 (ADM.57): the admin
  has complete control; the world adjusts to what the admin wants, after
  the admin approves the generation plan. A randomised test over many
  seeded systems: apply random plans (moves, orbit edits, class changes,
  star changes) and assert every system still validates and no locked
  aspect changed; undo them in reverse and assert every stored row
  equals the original. Done: the test, runs in the gate with a fixed
  seed set and in `full` with a larger one.
  Prerequisites: ADM.76, ADM.60, ADM.77. Related: ADM.57.
  Design: [docs/design/edit-everything.md](design/edit-everything.md)

- [ ] **TEST.177 Edited sectors survive repair, regeneration and the golden seed**
  Boss (2026-10-16, 1:07 PM EDT), headline of 8.6 (ADM.57): the admin
  has complete control; the world adjusts to what the admin wants, after
  the admin approves the generation plan. A golden-seed test with edits
  applied: rebuild the sector from its seed (DB.17) and from parity
  (DB.9) and assert the replayed journal gives the edited rows; assert a
  generator epoch change lists the rows that cannot replay instead of
  dropping them. Done: the test.
  Prerequisites: ADM.61, ADM.80. Related: ADM.57, TEST.77.
  Design: [docs/design/edit-everything.md](design/edit-everything.md)

- [ ] **TEST.178 Roles end to end: user, moderator and admin through the whole request chain**
  Boss (2026-10-16, 1:43 PM EDT), headline of 10.0 (USR.25): three
  roles, User, Moderator and Admin; users submit edit requests,
  moderators approve user requests, admins approve moderator requests,
  and an admin's own edits need no approval (only the preview and
  confirm of ADM.66). Build on the accounts of 8.2 and the change plans
  of 8.6. Browser and server tests: a user requests, a moderator
  approves and it applies; a moderator requests and only an admin can
  approve; an admin edits with preview and confirm and no second
  approver; a self-approval is refused; a role change takes effect at
  once; a hidden request stays hidden. Done: one test per case.
  Prerequisites: USR.18, USR.19, USR.20, USR.21, USR.22, USR.23, USR.24.
  Related: USR.25.
  Design: [docs/design/roles-and-requests.md](design/roles-and-requests.md)

- [ ] **TEST.179 Permission matrix test: every route and page checked for every role against the catalogue**
  Boss (2026-10-16, 1:43 PM EDT), headline of 10.0 (USR.25): three
  roles, User, Moderator and Admin; users submit edit requests,
  moderators approve user requests, admins approve moderator requests,
  and an admin's own edits need no approval (only the preview and
  confirm of ADM.66). Build on the accounts of 8.2 and the change plans
  of 8.6. A generated test walks the route map and the page list and
  asserts the response for user, moderator, admin and Owner matches the
  catalogue, so a new route without a permission fails. Done: the test
  in the gate.
  Prerequisites: USR.14. Related: USR.25.
  Design: [docs/design/roles-and-requests.md](design/roles-and-requests.md)

- [ ] **TEST.180 Backup and repair end to end: take a full archive, lose the control schema and the live parity file, restore and repair, and prove the galaxy and the logins come back**
  Boss (2:00 PM EDT): "we also have to build a way to use that data so
  build that out and let's make this targetted for v9.0 so incorporate
  the chain for these (restore / repair using backups) with that due
  date. Settings page is v8.1". Target 9.0 (the backups and repair
  chain, umbrella ADM.91); the settings page and the settings backups
  stay in 8.1 (ADM.83 to ADM.85). All backup pieces are parts of one
  `.tar.bz2` archive (ADM.87); retention is `settings.backups_kept`
  (default 20, 1 to 200). One scenario test on a small galaxy: take a
  full archive (settings, control, parity), damage sectors in the galaxy
  database, delete the live parity file, drop the control schema, then
  restore the control part, repair from the archive, run the check and
  sign in again; also an archive with a corrupted part refused, one-part
  restore leaving the others alone, retention keeping the newest full
  archive, and every restore refused while the update guard is on. Done:
  the test in the gate, with the run named in the PR.
  Also covers the nightly step (OPS.55): due and not due, and the
  pruning by kind.
  Prerequisites: ADM.89, DB.42, ADM.90, OPS.55. Related: ADM.91.

- [ ] **TEST.181 The order-independence test: the same sectors at 1 and 4 workers and in two orders**
  Split from GEN.57 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work). A test generates the same small set
  of sectors (one reached by a backfill, two sharing an ID cell) at 1
  and at 4 workers and in two orders and compares their fingerprints
  (GEN.58's encoder now, GEN.136:a when it lands); nearest-system links
  and population data stay out of the comparison (they are rebuilt from
  content). Why: this is the proof GEN.57 states, and TEST.77's tier B
  stands on it. Done: the test runs in the `changed` gate and fails when
  a or b is reverted.
  Prerequisites: GEN.306, GEN.307. Related: GEN.57.

- [ ] **TEST.182 Put one block-first backfill into the golden galaxy**
  Split from GEN.42 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work). TEST.77's golden galaxy runs a
  backfill whose blocks are too small to take the candidate path; force
  the candidate path for one block (a test hook) and update the golden
  through the `golden_update` tool, with the epoch bump. Done: the
  golden contains a block-first backfill and fails when the candidate
  path's draws change.
  Prerequisites: GEN.309, TEST.77. Related: GEN.42.

- [ ] **TEST.183 A test that fails when a new output path formats its own lines**
  Split from UX.102 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work). A guard test (an AST scan with an
  allow-list) that fails when code under `planetgen/` writes a log or
  report line with its own timestamp or divider instead of the layer.
  Done: the test, with the allow-list empty except for named exceptions
  with a reason.
  Prerequisites: UX.122, UX.123. Related: UX.102.

- [ ] **TEST.184 Whole-sphere tests, nebula and Wolf-Rayet side: contents per class, the Wolf-Rayet count, tiling, determinism and the King solver**
  Split from TEST.159 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). A planetary nebula always has its
  white dwarf; a core-collapse remnant has a compact core at the stated
  share; every hosted nebula holds its host star; a Wolf-Rayet nebula
  holds a W star and the Milky Way count is 1,200 +- 300 with none in
  old populations; tiling and the object ID round trip at the widest
  shell of the reach; the profile and an empty sector are deterministic;
  the King solver against the textbook concentrations. Done: each test
  exists and runs in the normal suite.
  Prerequisites: GEN.320, GEN.243, TEST.186. Related: TEST.159.

- [ ] **TEST.185 Whole-sphere tests, component side: pre-placed share per population, no double count of a cluster star, nothing hosted in a halo or old cluster cell, D and E members inside their radius**
  Split from TEST.159 (2026-10-16, Boss 1:07 PM EDT: look at each item
  for smaller chunks and shared work). No emission, reflection or
  core-collapse remnant in a halo or old cluster cell; every D or E
  region holds its members inside its radius; the pre-placed share of
  every population; no double count of a cluster star between the passes
  and the fill. Done: each test exists and runs in the normal suite.
  Prerequisites: GEN.240, GEN.321, TEST.186. Related: TEST.159.

- [ ] **TEST.186 Census audit harness: fill seeded sectors per region class, count by kind and stage, compare to expectation within k Poisson sigma, print one report**
  Found in the split review (2026-10-16, Boss 1:07 PM EDT: look at each
  item for smaller chunks and shared work, 1:33 PM EDT: shared
  foundations). TEST.160 builds the harness as a module
  (`src/tests/audit/census.py`): a fixed set of seeded sectors per
  region class, counts by stage and kind, `within_sigma(observed,
  expected, k=4)`, the 'expected 20 or more appears at least once' rule,
  and the report writer (the admin report of TEST.160). TEST.159,
  TEST.161, TEST.162 and TEST.164 add their checks to it; their own Done
  lines already say 'the representation audit report includes them'.
  Done: the module and TEST.160's audit on it; the other four add
  checks, not counters.
  Prerequisites: none. Related: TEST.159, TEST.160, TEST.161, TEST.162,
  TEST.164.

- [ ] **TEST.187 Browser-test failure reports are buried in the werkzeug request log: set that logger to WARNING in the test server**
  Found by a lane (2026-10-16), filed by the TODO thread. When a browser
  test fails, the pytest report ends with hundreds of lines of werkzeug
  request log ("GET /static/... 200"), and the assertion text sits above
  them. Done: the test server fixture sets the `werkzeug` logger to
  WARNING, a failing browser test report shows the assertion and the
  page console first, and one test proves the request log stays out of
  the captured output. Fails on main: yes (every failure). Frequency:
  every browser test failure. Reproducer: make any web.* browser test
  fail and read the report.
  Prerequisites: none.

- [ ] **TEST.188 The skimage marching-cubes call reads `.shape` in a way NumPy 2.5 deprecates: change the call and keep the warning out of the log**
  Found by a lane (2026-10-16), filed by the TODO thread. The nebula
  mesh tests emit a NumPy deprecation warning from the skimage
  marching-cubes path (`.shape` assignment on an array). Done: the call
  uses the supported form (`reshape`), the test run has no deprecation
  warning from this path, and the filterwarnings list holds no entry for
  it. Fails on main: no failure, a warning on every nebula mesh test.
  Frequency: every run. Reproducer: `python -W error::DeprecationWarning
  -m pytest -q src/tests -k nebula_shape`.
  Prerequisites: none.

- [ ] **TEST.189 A test that runs `migrate` leaves maintenance mode on: restore `state/maintenance.json` in the fixture (bug)**
  Found by a lane (2026-10-16), filed by the TODO thread. A test that
  runs the `migrate` command leaves
  `/var/lib/planetGen/state/maintenance.json` set, so later tests (and
  the web server on a developer machine) answer with the maintenance
  page. Done: the fixture that runs migrations restores the maintenance
  file to its earlier state on exit, whether the test passes or fails,
  and a test asserts the file is unchanged after a migration test. Fails
  on main: yes, when that test runs. Frequency: every run of the test.
  Reproducer: run the migrate test, then `cat
  /var/lib/planetGen/state/maintenance.json`.
  Prerequisites: none.

- [ ] **TEST.190 A failing session fixture makes every dependent test wait its full 900 s timeout: fail the dependents at once (bug)**
  Found by a lane (2026-10-16), filed by the TODO thread. When a
  session-scoped fixture (for example the browser-test server) fails to
  start, each test that needs it waits its whole 900 s timeout before
  failing, so one fixture failure costs hours of the run. Done: a
  fixture that fails marks its dependents as errors at once (pytest does
  that already for a raised error, so the cause is in how the fixture
  waits; find it and fix it), and a test shows a dead fixture ends the
  run in seconds. Fails on main: yes. Frequency: whenever a session
  fixture fails. Reproducer: make the server fixture raise after its
  start and time the run.
  Prerequisites: none. Related: TEST.167.

- [ ] **TEST.191 The plan tests write real settings JSON into `/var/lib/planetGen/jobs/galaxy-settings` and log "No naming key drawn ... Unknown database": give them a temporary directory and a test database**
  Found by a lane (2026-10-16), filed by the TODO thread. The `plan`
  tests write galaxy settings files into the real jobs directory on the
  machine and log "No naming key drawn ... Unknown database", so a
  developer machine collects stray files and the naming path runs
  against no database. Done: the plan tests use a temporary jobs
  directory and the test database, the log line is gone, and a test
  asserts the real directory is untouched. Fails on main: no failure, a
  side effect on every run. Frequency: every run. Reproducer: run the
  plan tests and list the jobs directory.
  Prerequisites: none.

- [ ] **TEST.193 A module fixture's `DROP DATABASE IF EXISTS` blocks for 853 s on a metadata lock that another connection holds, so a browser test hangs in teardown (bug)**
  Found by Core 2 (2026-10-16, 5:04 PM EDT), filed by the TODO thread.
  In the gate for PERF.88,
  `src/tests/test_web_browser_sector_reach.py::test_no_scene_is_loaded_while_the_camera_is_far`
  hung in teardown: the module fixture `site_db`
  (`src/tests/test_web_a11y.py:275`) ran `DROP DATABASE IF EXISTS` and
  waited 853 s on a metadata lock, in the gate's first run and again in
  its retry. A connection left open by something (the same kind of lock
  as the interstellar_comets one Boss hit on a migration) holds the
  lock. Done: the fixture closes or kills the connections that hold the
  database before the drop, or fails at once with the name of the holder
  instead of waiting; the teardown never waits more than a few seconds;
  and a test shows a held connection no longer stalls it. Fails on main:
  not checked, the PERF.88 diff touches no test fixture. Frequency: once
  in 3 gates. Reproducer: none; the same file passed alone (7 of 7)
  right after. Related to TEST.190 (a failing session fixture makes
  dependents wait their full timeout).
  Core 2 reproduced it twice in a row (2026-10-16, 5:52 PM EDT, the
  OPS.56 gate and its rerun alone): the holder is an idle connection
  (state Sleep, database `planetgen_a11y_*`, the same MySQL) that stays
  open while `DROP DATABASE planetgen_a11y_*` waits 853 s; the PERF.88
  gate hit the same shape. Reproducer: run `scripts/test.py changed`
  with the browser tests selected and read `show full processlist` in
  mysql during the teardown.
  Sighting (Core 2, 2026-10-16, 8:12 PM EDT, DB.18 gate: scripts/test.py
  changed, 468 passed, 1 failed, 2 errors): three browser tests hung
  about 15 minutes in setup or teardown:
  test_web_browser_sector_open::test_opening_a_sector_is_one_continuous_view_with_a_filled_panel,
  test_web_browser_sector_reach::test_no_scene_is_loaded_while_the_camera_is_far
  and test_web_browser_units::test_the_server_text_is_the_automatic_text
  (passes alone). Fails on main: yes (the same two failed in the OPS.56
  gate). Frequency: every full gate this evening, so it is no longer rare
  and costs each lane about 15 minutes a gate; one root cause is likely.
  Owner: Core 1 (Boss, 2026-10-16, 8:23 PM EDT, card "Core 1 next"): starts right after
  its current stack and the GEN.372 audit merge, ahead of its bug phase.
  TEST.193 and TEST.195 are likely one root cause, so Core 1 takes both together.
  Prerequisites: none.

- [ ] **TEST.194 `test_system_map_selection_drill_and_measure` failed once in the browser phase of a full run under load and passed alone: find the cause and make it robust (bug)**
  Found by the Settings lane (2026-10-16, 5:47 PM EDT), filed by the
  TODO thread.
  `src/tests/test_web_browser_maps.py::test_system_map_selection_drill_and_measure`
  failed once in the browser phase of a full `scripts/test.py changed`
  run (load 8.7, 2 browser workers) and passed when rerun alone (call
  122.9 s on gw1). Done: the cause is found (a wait that is too short
  under load, or a real race in the system map drill), the test waits on
  the condition it needs and not a fixed time, and it passes under the
  same load with 2 browser workers, with no quarantine. Fails on main:
  unknown. Frequency: 1 of 1 run under load. Reproducer: none beyond the
  full browser phase; try the file with two browser workers and a busy
  machine.
  Prerequisites: none.

- [ ] **TEST.195 `test_opening_a_sector_is_one_continuous_view_with_a_filled_panel` timed out at 30 s fetching `/sector/<id>/galaxy` inside a full gate run and passed alone (bug)**
  Found by Core 2 (2026-10-16, 5:52 PM EDT), filed by the TODO thread.
  `src/tests/test_web_browser_sector_open.py::test_opening_a_sector_is_one_continuous_view_with_a_filled_panel`:
  `page.request.get` of `/sector/<id>/galaxy` timed out at 30 s in the
  full gate run with 2 browser workers, and the test passed alone. Done:
  the cause is found (a slow first request that builds the whole galaxy
  scene under load, or a real stall in the route), and the test either
  waits on the condition it needs or the route is made fast enough, with
  no quarantine. Fails on main: unknown. Frequency: once so far.
  Reproducer: none beyond the full browser phase with 2 workers.
  Sighting (Core 2, 2026-10-16, 8:12 PM EDT, DB.18 gate: scripts/test.py
  changed, 468 passed, 1 failed, 2 errors): three browser tests hung
  about 15 minutes in setup or teardown:
  test_web_browser_sector_open::test_opening_a_sector_is_one_continuous_view_with_a_filled_panel,
  test_web_browser_sector_reach::test_no_scene_is_loaded_while_the_camera_is_far
  and test_web_browser_units::test_the_server_text_is_the_automatic_text
  (passes alone). Fails on main: yes (the same two failed in the OPS.56
  gate). Frequency: every full gate this evening, so it is no longer rare
  and costs each lane about 15 minutes a gate; one root cause is likely.
  Owner: Core 1 (Boss, 2026-10-16, 8:23 PM EDT, card "Core 1 next"): starts right after
  its current stack and the GEN.372 audit merge, ahead of its bug phase.
  TEST.193 and TEST.195 are likely one root cause, so Core 1 takes both together.
  Prerequisites: none.

- [ ] **TEST.196 Close-binary check: generated close-binary fractions by primary mass and [Fe/H] bin against the published ranges, with out-of-range bins flagged**
  Boss (2026-10-16, 6:15 PM EDT, on the GEN.103 cards): the metallicity
  effect on close binaries applies to every star mass, "but let's build
  into that a check to see if the numbers act as expected". Done: a
  check in the census harness (TEST.186, `src/tests/audit/census.py`)
  that fills seeded sectors across the [Fe/H] range, counts the
  close-companion share (inside about 10 AU) by primary-mass bin and
  [Fe/H] bin, compares each bin to the published range with the Poisson
  allowance (the Moe, Kratter and Badenes 2019 points, as Core 4 reads
  them: 53%, 40%, 24% and 10% at [Fe/H] -3.0, -1.0, -0.2 and +0.5 for
  solar-type stars, and whatever the sources give for other masses),
  flags every out-of-range bin in the audit report, shows that the
  wide-companion share does not move with [Fe/H], and fails the test run
  when a bin is out of range. The ranges are reference values for the
  check, never inputs to the draw. Owner: Core 4 (near TEST.186, which
  is Core 3's harness; ask Core 3 before editing census.py).
  Status (Core 4, 2026-10-16, 7:16 PM EDT): TEST.186 is not on main yet,
  so the check cannot be built. The note lists the expected ranges. The
  massive-star flatness bins fail with plain scaling of the metallicity
  factor; with the odds form Core 4 recommends in GEN.103 they are
  expected to pass at the edge. The test lands only with a green run (no
  skip, no quarantine). The 3.2 prevalence cap is GEN.373 and does not
  gate this test.
  Prerequisites: GEN.103, TEST.186. Related: GEN.103, TEST.186,
  TEST.161, TEST.164.

- [ ] **TEST.197 `test_page_is_accessible_and_fits[web.sector-phone-light]` failed under load and passed alone: find the cause and make it robust (bug)**
  Found by Core 3 (2026-10-16, 6:22 PM EDT), filed by the TODO thread.
  `src/tests/test_web_a11y.py::test_page_is_accessible_and_fits[web.sector-phone-light]`
  failed in a full `changed` run (78 minutes, load 9.6) on a branch that
  touches neither maps nor the sector page, and passed alone. Done: the
  cause is found (a timing assumption under load, or a real layout race
  on the phone sector page), the test waits on the condition it needs
  and not a fixed time, and it passes under the same load, with no
  quarantine. Fails on main: not checked. Frequency: once. Reproducer:
  none beyond the full suite at load about 9.
  Prerequisites: none.

- [ ] **TEST.198 `test_galaxy_map_does_not_prefetch_on_a_data_saving_connection` failed with a `/galaxy/tiles` net::ERR under load and passed alone: find the cause and make it robust (bug)**
  Found by Core 3 (2026-10-16, 6:22 PM EDT), filed by the TODO thread.
  `src/tests/test_web_browser_fixture_maps.py::test_galaxy_map_does_not_prefetch_on_a_data_saving_connection`
  failed with "request failed ... /galaxy/tiles?... net::ERR" in the
  same full run (load 9.6) and passed alone. Done: the cause is found
  (the test or the server drops a request when a page closes under load,
  or the test reads the request log too early), the test waits on the
  condition it needs, and it passes under the same load, with no
  quarantine. Fails on main: not checked. Frequency: once. Reproducer:
  none beyond the full suite at load about 9.
  Prerequisites: none.

- [ ] **TEST.199 `test_skeleton_bound_is_an_upper_bound_at_every_angle` fails on main for a far-tail case: the halo subtraction loses precision, 2.4e-9 against a 1e-9 tolerance (bug)**
  Found by Core 2 (2026-10-16, 7:27 PM EDT, in its ADM.50 gate), filed
  by the TODO thread. Kind: test, deterministic for the shrunk case.
  `src/tests/test_fuzz_galaxy_grid.py::test_skeleton_bound_is_an_upper_bound_at_every_angle`
  failed once in an 88 minute `changed` run, found by hypothesis, and
  fails on main. Failure text: `assert (1.4882030977187565e-12 -
  1.488203070457534e-12) <= ((2.7261222351049663e-20 * (1 + 1e-09)) +
  1e-300)`. Reproducer:
  `shape=build_galaxy_shape(disk_scale_length_pc=50.0,
  disk_scale_height_pc=259.0, bulge_scale_radius_pc=10.0,
  bulge_amplitude=0.0, arm_count=1, pitch_angle_rad=1.0,
  arm_amplitude=0.0)`, `r_cyl=0.0`, `z=36706.0`, `theta=0.0`: `k_norm *
  _raw_density` minus `halo_density` is 2.7261222417e-20 against a
  skeleton bound of 2.7261222351e-20, a relative gap of 2.4e-9 against
  the 1e-9 tolerance, so subtracting the halo from a sum two orders of
  magnitude larger loses precision in the far tail. Done: the cause is
  settled (a cancellation in the subtraction: compute the component
  density without subtracting the halo from a total, or state the
  tolerance the cancellation allows and why), the bound holds for this
  case and a test pins it, and hypothesis finds no further case in a
  full run; no quarantine, no loosened tolerance without the stated
  reason. It is also a precision finding, so it is attached to GEN.372
  (the precision audit). Related: TEST.192 (an earlier fuzz bound).
  Owner: Core 1 (proposed: its bug phase and GEN.372).
  Prerequisites: none. Related: GEN.372, TEST.192, GEN.301.

- [ ] **TEST.200 The nightly local offset limit is documented and pinned by `test_local_offset_after_a_rotation_is_within_one_ulp_of_the_absolute_point`**
  From the GEN.372 precision audit (Core 1,
  docs/design/numeric-precision-audit.md section 5, 2026-10-16, 7:39 PM
  EDT), filed by the TODO thread; verdict: 64-bit is enough for every
  stored value. Boss decided by card: the nightly update keeps a local
  offset, and this item documents it and tests it, with no code change
  to the update. Done: the limit (how far an object's local offset may
  grow before it is folded into the absolute point) is written in the
  orbital update's design doc, and
  `test_local_offset_after_a_rotation_is_within_one_ulp_of_the_absolute_point`
  pins that after the largest rotation the update applies the local
  offset is within one unit in the last place of the absolute point.
  Owner: Core 1.
  Prerequisites: none. Related: GEN.372, GEN.223, GEN.294.

- [ ] **TEST.201 `block_stats` totals do not depend on the order of the rows: a test**
  From the GEN.372 precision audit (Core 1,
  docs/design/numeric-precision-audit.md section 5, 2026-10-16, 7:39 PM
  EDT), filed by the TODO thread; verdict: 64-bit is enough for every
  stored value. Done: a test sums the block statistics (`block_stats`)
  of the same rows in several random orders and requires identical
  totals (or totals within a stated bound), failing if the sum is order
  dependent. Compensated sums for the long accumulations wait for
  GEN.109 (Boss card): see the note on GEN.109. Owner: Core 1.
  Prerequisites: none. Related: GEN.372, GEN.109.

- [ ] **TEST.202 `test_galaxy_map_drill_down_by_clicks` failed under load and passed alone: call 28 s to 52 s on a loaded worker (bug)**
  Found by Core 3 (2026-10-16, 8:34 PM EDT, in the GEN.333 gate: full
  `changed`, 118 minutes, load 16), filed by the TODO thread.
  `src/tests/test_web_browser_maps.py::test_galaxy_map_drill_down_by_clicks`
  failed under load and passed alone (call 52 s on a loaded worker); the
  same gate also failed
  `test_galaxy_map_forward_to_current_jumps_to_the_newest_view` (call 28
  s; TEST.172 carries that sighting). Fails on main: not checked.
  Frequency: once. Reproducer: none beyond a loaded full run. Done: the
  cause is found (a fixed wait or a click before the stage view is
  ready, or a real slow path in the drill-down), the test waits on the
  condition it needs, it passes under the same load, and no quarantine.
  Likely the same family as TEST.197 (a timing assumption under load)
  and the browser hangs of TEST.193 and TEST.195, so Core 1 reads them
  together. Owner: Core 1, bug phase.
  Prerequisites: none. Related: TEST.172, TEST.197, TEST.193, TEST.195.

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

  Headline: 8.2 (user accounts, with sign-up by an email invitation that
  only an admin can send).
  - [ ] **USR.2 Accounts with roles: user, moderator, admin and Owner**
    Boss: "Admin can
    then make a user admin or take away admin rights on everything but the
    primary 1st admin account generated at install, that'll have a
    designation of 'Owner' and no other account can override that." Done:
    one accounts table in the control database (with an email address)
    holding a role per account (user, moderator, admin, Owner; moderator
    is stored and valid from the start, has no permissions until 10.0
    and no page assigns it, Boss 2026-10-16, 1:43 PM EDT); exactly one
    Owner,
    which is the account install creates; an admin page where admins
    promote a user to admin or demote an admin to user, refused for the
    Owner; sessions and API keys that work for every role while the
    admin-only pages and API routes stay admin-only; each role change
    written to the audit log. Question (decided below): does
    `admin_users` become
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
    Decided (2026-10-10): Boss (19:16Z, Accept all): one accounts table
    with a role column (add the columns now, rename later); the lowest
    id is the Owner.
    Roles must fit 10.0 (roles and permissions, USR.25) without a
    rebuild: the role is one text column holding one of the four values;
    every admin-only page, route and command goes through one guard
    (`require_role`) so 10.0 swaps that single function for the
    permission check `can()` (USR.12, USR.14); no code outside that
    guard compares a role; the role change page of this item writes the
    audit row through one function that 10.0 reuses for the moderator
    role (USR.13).
    Split (2026-10-16, Boss 1:07 PM EDT: look at each item for smaller
    chunks and shared work): this item is now the umbrella and is done
    when its subitems are: DB.58, USR.26, USR.27. It mixes a
    control-database migration (role, email, status, single-Owner index
    and, per the later flows, the shared token and invite tables), the
    role model in sessions, API keys and guards (which USR.3, USR.4,
    USR.5, USR.7, USR.8 and API.6 actually wait for), an admin page and
    a console CLI.
    Prerequisites: DB.58, USR.26, USR.27.

  - [ ] **USR.3 SMTP settings in the admin config**
    Boss: "We'll use SMTP for
    email which means admin config needs SMTP settings." Done: SMTP host,
    port, security (STARTTLS or TLS), username, password and From address,
    set from an admin page (and the config file/installer), with a "send
    test email" button; one small mail module that every flow in USR.4
    to USR.6 uses, which logs failures and never shows the SMTP
    password. Question (decided below): is the SMTP password kept in the
    config
    file (like the database password) or in the control database, and is
    it encrypted there? Who can change SMTP settings: any admin, or only
    the Owner? What do invites and resets do when SMTP isn't configured
    (show the link to the admin to pass on by hand?)?
    Research (2026-10-09, settings-seo-and-accounts.md): answer as in
    B9: stdlib smtplib, Owner-only (SMTP, `base_url`, From address; open
    question for Boss, default yes), password in the user's
    `config.json` (replaceable, never read back), no SMTP means invites
    show the link only; sent from an RQ job.
    Decided (2026-10-10): Boss (19:16Z, Accept all): SMTP settings are
    Owner-only; the password lives in the user's `config.json`, not the
    database.
    The SMTP form sits on the settings page of 8.1 (ADM.84) and its
    Owner-only tier is ADM.86 (Boss, 2026-10-16, 1:43 PM EDT).
    Scope after the review (2026-10-16 review): adds the smtp.* options
    (Owner tier, the password in `config.json`) to the settings model
    and builds the mail module and the "send test email" action; the
    form is the generated settings page (ADM.84 and ADM.86), not a
    second form.
    Prerequisites: ADM.86, USR.26.

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
    an admin. Question (decided below): does an admin optionally type
    the invitee's
    email so the link is sent for them, or only copy the link? Is the
    invite page rate-limited, and is there a cap on open invites? Does
    a multi-use link record who used it?
    Research (2026-10-09, settings-seo-and-accounts.md): add the
    `invite_uses` table, atomic redeem, an open-invite cap of 50 and
    rate limits (B5).
    Decided (2026-10-10): Boss (19:16Z, Accept all): invites capped at
    50 open at once; the invitee email is optional.
    Scope after the review (2026-10-16 review): the invite page (use
    count, expiry choices, confirmation for infinite plus never, the
    open-invite cap of 50, revoke, the list) and registering through a
    valid link; tokens and `invite_uses` come from the account-token
    library and the first child of USR.2. Placed after USR.5 because its
    sign-up ends in the set-password email loop.
    Prerequisites: DB.58, USR.26, USR.3, USR.28, USR.5.

  - [ ] **USR.5 Email loop for setting and resetting passwords**
    Boss: "an
    email loop for password setting / resetting". Done: a new account
    sets its first password from an emailed link; "forgot password" on
    the login page emails a reset link; links are single-use, short-lived
    (stored hashed) and end the account's other sessions once used; the
    page never says whether an email address has an account; resets are
    rate-limited per address and per IP alongside the existing login
    backoff and SEC.1. Question (decided below): how long a reset link
    lasts
    (30 minutes? 1 hour?); does changing the email address also need an
    email confirmation to the old and new addresses; does the Owner's
    reset need anything extra?
    Research (2026-10-09, settings-seo-and-accounts.md): lifetimes:
    reset 30 minutes, first-password 24 hours, email change 1 hour; a
    notification email on password and email change; an Owner reset asks
    for the TOTP code when enrolled.
    Decided (2026-10-10): Boss (19:16Z, Accept all): reset links last 30
    minutes, a first password link 24 hours, an email change 1 hour with
    confirmation.
    Scope after the review (2026-10-16 review): the emailed set-password
    and forgot-password flows, per-address and per-IP limits, ending
    other sessions on use; tokens come from the account-token library.
    Placed before USR.4, which needs it.
    Prerequisites: USR.3, USR.28.

  - [ ] **USR.6 Owner transfer**
    Boss: "The owner CAN (with specific approval
    and double password confirmation and email loop confirmation) assign
    owner to someone else who then has to accept via email loop." Done:
    only the Owner can start a transfer, from a page that asks them to
    confirm the choice explicitly, enter their password twice, and then
    confirm from an email link; the chosen account then gets an email and
    must accept from its own link; only when both are done does Owner
    move; every step goes to the audit log. Question (decided below):
    what the old
    Owner becomes (admin?); does the new Owner have to be an admin
    already; how long the pending transfer lasts and whether the Owner
    can cancel it; what happens if the Owner loses their email or
    password (recovery from the server's command line?).
    Research (2026-10-09, settings-seo-and-accounts.md): the state
    machine of B10: the old Owner becomes admin; the target must be an
    admin with TOTP; it expires after 48 hours; the Owner can cancel;
    the confirm link lands on a neutral page because of
    `SameSite=Strict`.
    Decided (2026-10-10): Boss (19:16Z, Accept all): the old Owner
    becomes an admin; the target must be an admin with 2FA; a pending
    transfer expires in 48 hours.
    Scope after the review (2026-10-16 review): the Owner-transfer state
    machine and its pages; tokens come from the account-token library.
    Prerequisites: USR.5, USR.28.

  - [ ] **USR.7 A user-level interface with bookmarks**
    Boss: "a full user
    level interface to allow users to bookmark". Done: signed-in users
    (any role) get an account page and can bookmark sectors, systems,
    planets, phenomena and NAV courses, see them in a list, name them and
    remove them; bookmarks are stored per account in the control
    database. Question (decided below): which objects can be bookmarked,
    and can
    users also add notes? What else a user can do that an anonymous
    visitor can't (is the site still public to read, or sign-in only?)?
    Do bookmarks survive a galaxy regenerate (object ids change), and if
    not, what does a broken bookmark show? Per-browser bookmarks
    (`static/bookmarks.js`, MAP.23) shipped in PR #234; storing them in
    the database still needs decision 4 of the drill-down design and a
    migration.
    Research (2026-10-09, settings-seo-and-accounts.md): nothing new
    beyond B11 (forced sign-in would force `seo.indexing=off`).
    Decided (2026-10-10): Boss (19:16Z, Accept all): the site stays
    public to read; bookmarks are for signed-in users; notes come later.
    Scope after the review (2026-10-16 review): the account backend of
    the saved-list store (one `user_saved(account_id, kind, name,
    payload, galaxy_instance, created_at)` table in the control database
    for every kind (bookmark, course), the routes for a user's own rows
    only, the account page and the one-time import of
    the browser list); the list and rename UI is UX.47's.
    Prerequisites: DB.58, USR.26, UX.129.

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
    and admin. Decided (Boss, 2026-10-10 18:44Z, default approved): a
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
    `src/planetgen/web/system_page.py`. Decided (Boss, 2026-10-10
    18:44Z, defaults approved; default): user sessions last 30 days
    while admins stay at 12 hours (moved to the second child of USR.2, a
    role-aware session lifetime that applies at every sign-in; this item
    keeps the page route, the per-user limit and the concurrency cap).
    Prerequisites: USR.26.

- [ ] **USR.9 `seo.privacy_note` text, "download my data" and "delete my account" on the account page**
  A `seo.privacy_note` text setting shown on the account page, plus
  "download my data" and "delete my account" (B11). Not asked for by
  Boss; drop if unwanted.
  Needs USR.7 (the account page and the bookmarks) and the roles child
  of USR.2; "delete my account" is the delete-by-self rule USR.2 sets,
  and the data export reads each account list through the saved-list
  store so saved courses join it when NAV.19 lands (2026-10-16 review).
  Prerequisites: USR.7, USR.26.
  Design: [docs/design/settings-seo-and-accounts.md](design/settings-seo-and-accounts.md)

- [ ] **USR.10 Users add their own objects (ships and stations) and track their movement on the map**
  GitHub issue [#1226](https://github.com/dwhagar/planetGen/issues/1226)
  (Boss, 2026-10-12 23:27Z): "Allow users to add their own objects, not
  planets or anything but like a ship or a station or some such that
  they can track the movement of in the map." Done: a signed-in user can
  add a ship or station object, give it a position and movement, and see
  it and its track on the maps; the object belongs to the user. Decided
  (Boss, 2026-10-14 13:20Z): a user object's movement is a course from
  the route planner (this may change later); everyone can see every
  user's objects; no limit per user yet. Owner: not assigned.
  Split (2026-10-16, Boss 1:07 PM EDT: look at each item for smaller
  chunks and shared work): this item is now the umbrella and is done
  when its subitems are: DB.63, USR.29, USR.30. USR.10 holds a migration
  (a table for user objects), the library and API (an object with a
  position, a course from the route planner and an owner), and the pages
  and map display.
  Scope after the review (2026-10-16 review): the umbrella of the table
  for ships and stations (a, Core 1), the vessel and station objects
  with ownership (b) and the pages and map display (c); POP.11 needs
  only a and b.
  Prerequisites: DB.63, USR.29, USR.30. Related: USR.8, GEN.105, NAV.4.

- [ ] **USR.11 A "home" field on user accounts: where the Galaxy Map zooms to when it opens**
  GitHub issue #1232 (Boss, 2026-10-14 13:27Z): "users have the option
  to specify where to zoom into. Add a 'home' field to user accounts
  that will specify the location to zoom into." Done: a signed-in user
  can set a home (a star system, sector or position, picked from the map
  or the object page), it is stored on the account (Alembic migration,
  Foundations lane 1 owns migrations), and the Galaxy Map opens flown in
  to it instead of the largest group of created sectors (MAP.190); users
  without a home get the default. Owner: not assigned.
  Scope after the review (2026-10-16 review): the home field on the
  account and the override of the default; the default (the largest
  group) comes from the charted-regions library (MAP.190).
  Prerequisites: MAP.189, MAP.190. Related: MAP.189, MAP.190, USR.8.

- [ ] **USR.12 Permission catalogue and one `can()` check: user, moderator, admin and Owner as data**
  Boss (2026-10-16, 1:43 PM EDT), headline of 10.0 (USR.25): three
  roles, User, Moderator and Admin; users submit edit requests,
  moderators approve user requests, admins approve moderator requests,
  and an admin's own edits need no approval (only the preview and
  confirm of ADM.66). Build on the accounts of 8.2 and the change plans
  of 8.6. Name every permission the site checks (edit own account,
  request an edit, approve a user request, approve a moderator request,
  confirm and apply a plan, manage roles, change settings, see the audit
  log) and the role-to-permission map for user, moderator, admin and
  Owner, as data in `planetgen/admin/permissions.py` with one
  `can(account, permission, object=None)` that every route, page and CLI
  command will call. Owner keeps everything an admin has plus the
  Owner-only options of the settings page. Done: the catalogue and the
  check with tests for every role and permission, and a generated table
  in docs.
  Prerequisites: USR.2. Related: USR.2, ADM.86.
  Design: [docs/design/roles-and-requests.md](design/roles-and-requests.md)

- [ ] **USR.13 Moderator role: assign and remove it on the admin accounts page, with the audit entry**
  Boss (2026-10-16, 1:43 PM EDT), headline of 10.0 (USR.25): three
  roles, User, Moderator and Admin; users submit edit requests,
  moderators approve user requests, admins approve moderator requests,
  and an admin's own edits need no approval (only the preview and
  confirm of ADM.66). Build on the accounts of 8.2 and the change plans
  of 8.6. USR.2 stores the role in one column and already accepts
  "moderator". Admins promote a user to moderator and demote a moderator
  to user on the accounts page; each change is audited (who, whom, from,
  to); the Owner and admins are not affected. Done: the page actions,
  the audit rows and tests that a moderator holds exactly the moderator
  permissions.
  Prerequisites: USR.12, USR.27. Related: USR.2, USR.6.
  Design: [docs/design/roles-and-requests.md](design/roles-and-requests.md)

- [ ] **USR.14 Move every route, page and command from "is admin" checks to permission checks**
  Boss (2026-10-16, 1:43 PM EDT), headline of 10.0 (USR.25): three
  roles, User, Moderator and Admin; users submit edit requests,
  moderators approve user requests, admins approve moderator requests,
  and an admin's own edits need no approval (only the preview and
  confirm of ADM.66). Build on the accounts of 8.2 and the change plans
  of 8.6. Replace every direct role test (`role == "admin"`, admin-only
  decorators, template checks) in the API, the web pages and the CLI
  with `can()`, so a new role needs no edit to a route. Done: no raw
  role comparison is left in `src/planetgen` (a test scans for them),
  and the existing admin and user tests pass unchanged.
  Prerequisites: USR.12. Related: USR.2.
  Design: [docs/design/roles-and-requests.md](design/roles-and-requests.md)

- [ ] **USR.15 Submit an edit request: a user proposes an edit and the plan engine computes what it would change**
  Boss (2026-10-16, 1:43 PM EDT), headline of 10.0 (USR.25): three
  roles, User, Moderator and Admin; users submit edit requests,
  moderators approve user requests, admins approve moderator requests,
  and an admin's own edits need no approval (only the preview and
  confirm of ADM.66). Build on the accounts of 8.2 and the change plans
  of 8.6. A user (or a moderator) chooses an object and an edit from the
  operations their role may request (the catalogue lists them; a user
  starts with text aspects such as names and descriptions, and no
  structural changes), writes a short reason, and the plan engine
  (ADM.58) validates it and stores it as a plan in state requested, with
  the requester, their role and the role that must approve. Done: the
  submit flow for the object pages, validation that refuses an operation
  the role may not request, the stored plan with its computed diff, and
  tests for each role.
  Prerequisites: USR.12, ADM.58, ADM.59. Related: USR.25, ADM.57.
  Design: [docs/design/roles-and-requests.md](design/roles-and-requests.md)

- [ ] **USR.16 Request queue: who sees what**
  Boss (2026-10-16, 1:43 PM EDT), headline of 10.0 (USR.25): three
  roles, User, Moderator and Admin; users submit edit requests,
  moderators approve user requests, admins approve moderator requests,
  and an admin's own edits need no approval (only the preview and
  confirm of ADM.66). Build on the accounts of 8.2 and the change plans
  of 8.6. A page listing requests: a user sees their own; a moderator
  sees every user request and their own; an admin sees every request;
  filters by state, object and requester. Done: the queue page and API
  route with the visibility rules checked through `can()`, and tests
  that nobody sees a request the rules hide.
  Prerequisites: USR.15. Related: USR.25.
  Design: [docs/design/roles-and-requests.md](design/roles-and-requests.md)

- [ ] **USR.17 Approve, reject or send back a request, by role: moderators decide user requests, admins decide moderator requests**
  Boss (2026-10-16, 1:43 PM EDT), headline of 10.0 (USR.25): three
  roles, User, Moderator and Admin; users submit edit requests,
  moderators approve user requests, admins approve moderator requests,
  and an admin's own edits need no approval (only the preview and
  confirm of ADM.66). Build on the accounts of 8.2 and the change plans
  of 8.6. The decision rule: a user request is decided by a moderator or
  an admin, a moderator request by an admin, an admin needs none; nobody
  decides their own request. Approve applies the plan through ADM.60
  (the audit shows the requester and the approver); reject and send-back
  keep the reason. Done: the decision rules with tests for every role
  pair (including the refused self-approval), and the apply through the
  plan job.
  Prerequisites: USR.16, ADM.60. Related: USR.25, ADM.60.
  Design: [docs/design/roles-and-requests.md](design/roles-and-requests.md)

- [ ] **USR.18 Approver screen: the requester's note, the plan preview and the decision buttons**
  Boss (2026-10-16, 1:43 PM EDT), headline of 10.0 (USR.25): three
  roles, User, Moderator and Admin; users submit edit requests,
  moderators approve user requests, admins approve moderator requests,
  and an admin's own edits need no approval (only the preview and
  confirm of ADM.66). Build on the accounts of 8.2 and the change plans
  of 8.6. The page an approver opens: the reason, the object, the plan
  preview of ADM.66 (before and after, objects touched, warnings), the
  requester's earlier requests, and approve, reject and send-back with a
  note. Done: the screen reusing the approval screen, and a browser test
  for a moderator approving a user request and an admin approving a
  moderator request.
  Prerequisites: USR.17, ADM.66. Related: USR.25, ADM.66.
  Design: [docs/design/roles-and-requests.md](design/roles-and-requests.md)

- [ ] **USR.19 Comments on a request between the requester and the approver**
  Boss (2026-10-16, 1:43 PM EDT), headline of 10.0 (USR.25): three
  roles, User, Moderator and Admin; users submit edit requests,
  moderators approve user requests, admins approve moderator requests,
  and an admin's own edits need no approval (only the preview and
  confirm of ADM.66). Build on the accounts of 8.2 and the change plans
  of 8.6. A short thread on each request (DB.39): the approver asks, the
  requester answers; the requester can edit and resubmit, which
  recomputes the plan. Done: the thread on the request page, resubmit,
  and tests for ordering and permissions.
  Prerequisites: DB.39, USR.17. Related: USR.25.
  Design: [docs/design/roles-and-requests.md](design/roles-and-requests.md)

- [ ] **USR.20 "My requests" in the user interface: state of each request and withdraw**
  Boss (2026-10-16, 1:43 PM EDT), headline of 10.0 (USR.25): three
  roles, User, Moderator and Admin; users submit edit requests,
  moderators approve user requests, admins approve moderator requests,
  and an admin's own edits need no approval (only the preview and
  confirm of ADM.66). Build on the accounts of 8.2 and the change plans
  of 8.6. The account page of USR.7 lists the user's own requests with
  their state and lets them withdraw a pending one. Done: the page
  section with a browser test.
  Prerequisites: USR.16, USR.7. Related: USR.25, USR.7.
  Design: [docs/design/roles-and-requests.md](design/roles-and-requests.md)

- [ ] **USR.21 Request lifecycle: stale, withdrawn and expired requests**
  Boss (2026-10-16, 1:43 PM EDT), headline of 10.0 (USR.25): three
  roles, User, Moderator and Admin; users submit edit requests,
  moderators approve user requests, admins approve moderator requests,
  and an admin's own edits need no approval (only the preview and
  confirm of ADM.66). Build on the accounts of 8.2 and the change plans
  of 8.6. A request is stale when an object it touches changed after it
  was submitted (the stale check of ADM.59); a stale request goes back
  to the requester to resubmit; a pending request expires after 30 days.
  Done: the three transitions with tests, and the queue marks them.
  Prerequisites: USR.17. Related: USR.25, ADM.59.
  Design: [docs/design/roles-and-requests.md](design/roles-and-requests.md)

- [ ] **USR.22 Request limits: a per-user cap on open requests, a rate limit, and blocking a user from requesting**
  Boss (2026-10-16, 1:43 PM EDT), headline of 10.0 (USR.25): three
  roles, User, Moderator and Admin; users submit edit requests,
  moderators approve user requests, admins approve moderator requests,
  and an admin's own edits need no approval (only the preview and
  confirm of ADM.66). Build on the accounts of 8.2 and the change plans
  of 8.6. A user may hold at most 5 open requests and submit at most 20
  a day (settings on the settings page); a moderator or admin can stop
  an account from requesting, with the reason in the audit. Done: the
  limits and the block with tests, using the same pattern as the login
  rate limit.
  Prerequisites: USR.15. Related: USR.25, SEC.1.
  Design: [docs/design/roles-and-requests.md](design/roles-and-requests.md)

- [ ] **USR.23 Email notices: a request was submitted, decided or sent back**
  Boss (2026-10-16, 1:43 PM EDT), headline of 10.0 (USR.25): three
  roles, User, Moderator and Admin; users submit edit requests,
  moderators approve user requests, admins approve moderator requests,
  and an admin's own edits need no approval (only the preview and
  confirm of ADM.66). Build on the accounts of 8.2 and the change plans
  of 8.6. Approvers are emailed when a request needs them, and
  requesters when it is decided, through the queue-job mail module of
  USR.3; each account chooses on or off. Done: the notices as queue
  jobs, the setting on the account page, and a test with a captured
  mailbox.
  Prerequisites: USR.17, USR.3. Related: USR.25, USR.3.
  Design: [docs/design/roles-and-requests.md](design/roles-and-requests.md)

- [ ] **USR.24 Audit view for roles and requests: who changed a role, who asked, who decided**
  Boss (2026-10-16, 1:43 PM EDT), headline of 10.0 (USR.25): three
  roles, User, Moderator and Admin; users submit edit requests,
  moderators approve user requests, admins approve moderator requests,
  and an admin's own edits need no approval (only the preview and
  confirm of ADM.66). Build on the accounts of 8.2 and the change plans
  of 8.6. One admin page, and a section on the account page, listing
  every role change and every request with its decision, filterable by
  account, object and date, read from the audit log and the plan history
  of ADM.68. Done: the view with tests that each row names the requester
  and the approver.
  Prerequisites: USR.17, ADM.68. Related: USR.25, ADM.68.
  Design: [docs/design/roles-and-requests.md](design/roles-and-requests.md)

- [ ] **USR.25 Headline of 10.0: roles and permissions, user edit requests and moderators**
  Boss (2026-10-16, 1:43 PM EDT), headline of 10.0 (USR.25): three
  roles, User, Moderator and Admin; users submit edit requests,
  moderators approve user requests, admins approve moderator requests,
  and an admin's own edits need no approval (only the preview and
  confirm of ADM.66). Build on the accounts of 8.2 and the change plans
  of 8.6. Umbrella of 10.0, broken into items in the 10.0 version file;
  it is done when its items are. Headline: 10.0 (roles and permissions:
  user edit requests, moderators and admins).
  Prerequisites: USR.13, USR.14, DB.39, TEST.178, TEST.179, DOC.55.
  Related: USR.1, ADM.57.
  Design: [docs/design/roles-and-requests.md](design/roles-and-requests.md)

- [ ] **USR.26 Roles in sessions, API keys and route guards; the Owner rules; the accounts CLI**
  Split from USR.2 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work). Sessions and API keys carry the
  account's role; admin-only pages and API routes stay admin-only for
  every key and session; the Owner cannot be demoted, disabled or
  overridden, an admin cannot demote themselves, admins disable users
  and only the Owner disables admins, delete is by the account itself or
  the Owner, each role change is written to the audit log. User sessions
  last 30 days and admin sessions 12 hours (decided on USR.8, but it is
  a session rule). The console CLI `planetgen.cli.accounts
  --set-owner/--reset-password` is the recovery path. Done: a test per
  role over the guarded routes and keys, the Owner rules each refused,
  the lifetimes, and the CLI.
  Includes the role-aware session lifetime (users 30 days, admins 12
  hours; Boss, 2026-10-10 18:44Z, default approved), moved from USR.8
  (2026-10-16 review).
  Prerequisites: DB.58, SEC.32. Related: USR.2.

- [ ] **USR.27 Admin Users page: list accounts, promote, demote, disable and delete**
  Split from USR.2 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work). An admin page lists accounts with
  role and status; admins promote a user to admin or demote an admin to
  user (refused for the Owner), disable users, and the Owner also
  disables admins and deletes accounts, using child b's rules. Done: the
  page does each action, refuses what the rules refuse, and every change
  shows in the audit log.
  Prerequisites: USR.26. Related: USR.2.

- [ ] **USR.28 Account tokens: one issue-and-redeem library for invite, first-password, reset, email-change and owner-transfer links**
  Found in the split review (2026-10-16, Boss 1:07 PM EDT: look at each
  item for smaller chunks and shared work, 1:33 PM EDT: shared
  foundations). A small library over the `account_tokens` table (from
  USR.2:a): issue a random token stored hashed with a purpose (invite,
  first_password, reset, email_change, owner_transfer_old,
  owner_transfer_new), an account and an expiry; redeem atomically so a
  token works once and fails after expiry; revoke; the per-purpose
  lifetimes of Boss's decisions (reset 30 minutes, first password 24
  hours, email change 1 hour, owner transfer 48 hours); a hook for the
  rate limits the flows add. No flows and no pages. Done: tests for
  single use under two concurrent redeems, expiry, wrong purpose, hash
  only in the database, and the lifetimes.
  Prerequisites: DB.58. Related: USR.4, USR.5, USR.6.

- [ ] **USR.29 Vessel and station objects: add, move along a course and list, with ownership**
  Split from USR.10 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work). The library and API for the object of
  child a: a signed-in user adds a ship or station with a position and a
  course taken from the route planner (NAV.4; this may change later,
  Boss 2026-10-14), the object is owned by its user and visible to
  everyone (no limit per user yet), the orbital update moves it along
  its course, and each step is recorded as its track. The same function
  is what the population pass calls for a generated vessel (POP.11).
  Done: add, list, get and track routes with the usual API scopes,
  ownership checks, and tests for a course, a track and a refused edit
  by another user.
  Prerequisites: DB.63, NAV.4, GEN.105. Related: USR.10.

- [ ] **USR.30 Add a ship or station from a page and see it and its track on the maps**
  Split from USR.10 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work). The signed-in user's page to add an
  object, the object page, and the drawing of the object and its track
  on the Galaxy Map, the Sector Map and the System Map where the
  position falls. Done: the pages and the track drawing, and a browser
  test that adds a ship, sees it on the Sector Map and sees its track
  grow after an orbital update.
  Prerequisites: USR.29. Related: USR.10.

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
  the seed itself is not changed. The current galaxy's internal record
  (ADM.50, replacing ADM.18's settings file) keeps the epoch and
  `fp_spec`. Decided (Boss,
  2026-10-10 18:44Z, defaults approved; defaults taken): (1) is the same
  galaxy on Linux, macOS and
  Windows a goal, exact for integers, strings and structure and equal to
  9 significant digits for floats? Default yes; if not, drop the Windows
  and macOS CI additions and keep the epoch. (2) Is an integer epoch in
  the version record acceptable? Default yes. (3) Does a mismatch refuse
  a remote or rebuild run? Default: refuse on an epoch mismatch, warn on
  a battery-digest mismatch naming the first differing case.
  Note (2026-10-09): The Windows parts are dropped by OPS.39 (Boss,
  2026-10-10 04:25Z).
  Split (2026-10-16, Boss 1:07 PM EDT: look at each item for smaller
  chunks and shared work): this item is now the umbrella and is done
  when its subitems are: OPS.56, OPS.57, OPS.58. Three deliverables with
  three different clocks: the epoch counter and its PR rule (needs
  nothing), the battery digest (needs GEN.135 and the encoder), and the
  history row plus the galaxy record (a control-database migration and
  the update script).
  Prerequisites: OPS.57, OPS.58.
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
  Decided (Boss, 2026-10-10 18:44Z, defaults approved; default: print
  the manual command and end
  non-zero on a failed reload).
  Prerequisites: none.
  Design: [docs/design/ops-scheduling-and-rotation.md](design/ops-scheduling-and-rotation.md)

- [ ] **OPS.30 A lock helper for the maintenance run**
  `planetgen/util/locks.py` with `filelock` (pinned `==3.19.1` for
  Python 3.9, current otherwise), a sidecar info file, `timeout=0`
  helpers and `PermissionError` handling, plus a two-process test on
  Linux and macOS.
  Note (2026-10-09): The Windows parts are dropped by OPS.39 (Boss,
  2026-10-10 04:25Z).
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
  (DB.16) and in the current galaxy's internal record (ADM.50), shown on the admin status
  page beside the DB schema number, and returned by the API status
  response. The release version (MAJOR.REVISION.BUILD) stays as it is.
  Decided (Boss, 2026-10-09 21:02Z): the generator version and OPS.28's
  `generator_epoch` are the same number, bumped only when output changes
  for the same seed.
  Prerequisites: OPS.56.

- [ ] **OPS.55 Nightly full backup: a backup step in the daily maintenance run that takes every part into one archive, with a schedule setting and a nightly retention setting**
  Boss (2026-10-16, 2:06 PM EDT): "Yes plese." to a nightly full backup
  from the daily maintenance run. Target 9.0 (the backups and repair
  chain, umbrella ADM.91). The maintenance logic of OPS.16 (`--if-due`,
  `--only STEP`, lock, log, exit codes) gets a `backup` step that calls
  `planetgen.backup.take` of ADM.87 with every registered part
  (settings, control database ADM.88, parity data DB.40), so archives
  are taken without a click. The archive manifest carries a kind
  (`manual`, `settings-save`, `nightly`). The step is due when no
  nightly archive exists for the local day and the local hour has
  reached `settings.nightly_backup_hour` (0 to 23, default 4; the
  scheduler may call the script hourly), and it is switched off by
  `settings.nightly_backup` (default on). Nightly archives are pruned by
  their own count, `settings.nightly_backups_kept` (1 to 60, default 7),
  while `settings.backups_kept` (default 20, 1 to 200) keeps counting
  manual and settings-save archives, so a week of nightly archives never
  pushes out the archive an admin took by hand; the newest archive
  holding every part is never pruned. A failed step logs, exits non-zero
  like any other step and leaves the earlier archives alone; a run that
  finds the guard on (maintenance, an update, a run) skips and is
  retried when next due. The options go into the ADM.42 model and show
  on the settings page. Done: the step with tests for due and not due,
  off, a failure that leaves archives alone, pruning by kind, and a
  second run started during the first exiting at once; the daily run
  lists the step in its log.
  Prerequisites: OPS.16, ADM.88, DB.40. Related: OPS.16, ADM.87,
  ADM.91.

- [ ] **OPS.57 generator_battery_digest(): a database-free digest over 6 to 8 canned cases**
  Split from OPS.28 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work). A pure function in the package, no
  database and no network, over 6 to 8 tiny canned cases (a 3-system
  sector through `generate_sector` under a bound unit seed; one
  phenomenon; one nebula shape; one bright-star scatter on a tiny
  region; one population stream; one `unit_seed` known answer), each
  hashed with GEN.136's encoder into a structural and a float digest,
  returning the overall digest and a per-case digest (design
  generation-determinism.md 5.4, 3 to 5 s). Why: it is what OPS.15,
  OPS.14's mismatch check and TEST.77's tier A compare. Done: the
  function and a test that it is stable across two runs and two hash
  seeds, and that changing one canned draw changes exactly one case's
  digest.
  Prerequisites: GEN.135, GEN.304. Related: OPS.28.

- [ ] **OPS.58 The version history row and the galaxy's internal record carry the epoch, the battery digest and the environment**
  Split from OPS.28 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work). The `version_key_history` row
  (control schema v11, OPS.13) also holds `generator_epoch`, the battery
  digest, the lock-file hash and an environment JSON (libc, numpy,
  scipy, scikit-image versions): a control-database migration plus
  `planetgen.cli.version_history`. A row is added only when the key or a
  hash differs from the newest row (an unchanged update refreshes its
  last-seen date), the last 10 per galaxy are kept, and
  `requirements.lock` is pinned to LF in `.gitattributes`. The current
  galaxy's internal record (ADM.50) keeps the epoch and `fp_spec`.
  `version_check.galaxy_warning` names an epoch difference as the error
  and a battery-digest difference, with the first differing case, as a
  warning. Done: the migration with an old-schema fixture, the row, the
  record, the check, and tests for an unchanged update, a changed lock
  file and an epoch difference.
  Prerequisites: OPS.57. Related: OPS.28.

- [ ] **OPS.59 `generate galaxy --shell 0` does not stop on SIGTERM: the busy `worker --burst` child keeps running and the parent is left a zombie (bug)**
  Found by a lane (2026-10-16), filed by the TODO thread. Sending
  SIGTERM to the `generate galaxy --shell 0` command leaves its busy
  `worker --burst` child running and the parent a zombie, so a stopped
  run still holds Redis jobs and the database. Done: the command handles
  SIGTERM, forwards it to the worker child, waits for it and exits with
  the signal status; a test starts a run, sends SIGTERM and shows no
  process is left. Fails on main: yes. Frequency: every SIGTERM during a
  run. Reproducer: start the command, `kill -TERM` it, then `ps` for the
  child.
  Prerequisites: none.

## DOC: Documentation

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
  that every link inside the help section resolves. Decided (Boss,
  2026-10-10 18:44Z, defaults approved; default: Phase 2, after the
  features they describe have settled;
  Markdown source in docs/help/, built at update time, served at /help):
  other?
  Prerequisites: DOC.57, DOC.58. Related: DOC.7, DOC.8, DOC.9, DOC.10,
  DOC.11,
  DOC.12, DOC.13, DOC.14, DOC.15, DOC.16.
  Split (2026-10-16, Boss 1:07 PM EDT: look at each item for smaller
  chunks and shared work): this item is now the umbrella and is done
  when its subitems are: DOC.57, DOC.58. DOC.6 bundles the framework
  (template, index, '?' link, help search, version page, link test) with
  a coverage test that fails when any user-facing route lacks a help
  page.

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
  Prerequisites: DOC.57. Related: DOC.6.

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
  Prerequisites: DOC.57. Related: DOC.6.

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
  Prerequisites: DOC.57. Related: DOC.6.

- [ ] **DOC.10 The search and navigation help pages: search, nearby, the nav page and routes**
  Boss (2026-10-09 23:53Z): "Add to-do items to build static
  documentation pages for all features accessible through the web
  interface." Done: search (/search), nearby search (/nearby), the nav
  page (/nav), routes between systems and across sectors,
  within-N-parsec search, travel times (NAV.11, built), and how
  unknown space and asteroid fields affect a course. The page is a
  Markdown file in docs/help/ built into the Help section, linked from
  each page it describes, and updated by any later change to those
  features (a feature item is not done until its help text is).
  Prerequisites: DOC.57. Related: DOC.6.

- [ ] **DOC.11 The reference browser help pages: species, polities, object classes and phenomena**
  Boss (2026-10-09 23:53Z): "Add to-do items to build static
  documentation pages for all features accessible through the web
  interface." Done: /species, /polities, /classes and /phenomena with
  their detail pages: what each field means, how classes and codes are
  named, and how a phenomenon differs from a star system. The page is a
  Markdown file in docs/help/ built into the Help section, linked from
  each page it describes, and updated by any later change to those
  features (a feature item is not done until its help text is).
  Prerequisites: DOC.57. Related: DOC.6.

- [ ] **DOC.12 The account help pages: signing in, two-factor, the account page, API keys and bookmarks**
  Boss (2026-10-09 23:53Z): "Add to-do items to build static
  documentation pages for all features accessible through the web
  interface." Done: /login, the one-time code step, /account, two-factor
  setup, API keys and their scopes (API.9 when built), bookmarks (UX.47
  when built) and what is stored about a signed-in user. The page is a
  Markdown file in docs/help/ built into the Help section, linked from
  each page it describes, and updated by any later change to those
  features (a feature item is not done until its help text is).
  Prerequisites: DOC.57. Related: DOC.6.

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
  Prerequisites: DOC.57. Related: DOC.6.

- [ ] **DOC.14 The admin help pages: the queue, the stats page, settings, lockouts and the naming key**
  Boss (2026-10-09 23:53Z): "Add to-do items to build static
  documentation pages for all features accessible through the web
  interface." Done: /admin, /admin/queue and its tree, confirm and
  action pages, /admin/stats (naming key, lockouts; the galaxy settings panel goes with ADM.50),
  the worker count and the performance statistics (PERF.32), and which
  actions cannot be undone. The page is a Markdown file in docs/help/
  built into the Help section, linked from each page it describes, and
  updated by any later change to those features (a feature item is not
  done until its help text is).
  Prerequisites: DOC.57. Related: DOC.6.

- [ ] **DOC.15 A glossary and units help page: coordinates, scales, sector paths, object IDs, time and the in-universe wording**
  Boss (2026-10-09 23:53Z): "Add to-do items to build static
  documentation pages for all features accessible through the web
  interface." Done: the coordinate frames (galactic, sector, system),
  the unit ladder and when a unit switches (UX.22, UX.23, UX.78 (done, PR #1562), UX.81),
  sector addresses and paths, the object ID (GEN.170), the epoch and
  time scales, and the in-universe terms (UX.42, done), each with a worked
  example. The page is a Markdown file in docs/help/ built into the Help
  section, linked from each page it describes, and updated by any later
  change to those features (a feature item is not done until its help
  text is).
  Prerequisites: DOC.57. Related: DOC.6.

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
  Prerequisites: DOC.57. Related: DOC.6.

- [ ] **DOC.18 Last pass: re-run the documentation sweep when the other DOC items are done and Phase 1 is complete**
  Boss's instruction (2026-10-14 18:06Z): bring the docs current once
  Phase 1 is done. Items merged after the sweep (for example PERF.86,
  UX.102, OPS.45, the fly-through map MAP.151 to MAP.155, MAP.150, the
  logging standard) will have changed things the DOC items above cannot
  know about. Done: the sweep is run again on main, anything new is
  fixed in the same pass, README's headline feature list and the release
  headline are checked, and the result is recorded here.
  Prerequisites: DOC.19, DOC.20, DOC.24, DOC.47, DOC.48, DOC.49,
  DOC.50, DOC.51, DOC.52, DOC.53. Related: DOC.4.
  Runs last: re-run the sweep when the other 8.1 DOC items are done
  (2026-10-16 review).

- [ ] **DOC.19 After UX.102: document the log output standard for the CLI, the Generate page and admin logs**
  Requested by the coordinator from Boss's instruction of 2026-10-14
  18:06Z to compare the docs with merged work. UX.102 sets one logging
  standard (timestamps, format, dividers, step markers) for everything
  the user or an admin reads; once it is built, docs/cli.md,
  docs/html-interface.md (Generate page and job log) and docs/testing.md
  need a short section on it. Done: the standard is described once and
  the three docs link to it.
  Prerequisite: UX.102.

- [ ] **DOC.20 After PERF.86: rewrite the queue and scatter notes and the queue docstrings**
  Found by the full documentation sweep Boss asked for (2026-10-14
  18:06Z); every finding was checked against the code on main
  (8.0.1072), and the findings with file, line and evidence are in
  /mnt/project-files/reports/docs-sweep-2026-10-14/sweep-C.md, findings
  10 to 12, 14, 15 and part of 36. When each scatter becomes its own
  queue item run in turn (PERF.86), these describe the old model:
  generation-benchmark-report.md line 101 (scatter is still one job),
  performance-eta-queue-and-caching.md, work-queue-audit.md,
  sampling-backfill-and-resume.md section 7 (resume, which PERF.86
  decided against), library-migration.md section 3, architecture.md Flow
  4, and the docstrings of queue/work.py and queue/redisqueue.py. Done:
  all of them describe the queue as built.
  Prerequisite: PERF.86.

- [ ] **DOC.24 Bring Python 3.12, numba and astropy into the research notes and three source comments**
  Found by the full documentation sweep Boss asked for (2026-10-14
  18:06Z); every finding was checked against the code on main
  (8.0.1072), and the findings with file, line and evidence are in
  /mnt/project-files/reports/docs-sweep-2026-10-14/sweep-C.md, findings
  18 to 23. generation-determinism.md, reproducible-galaxies.md,
  units-and-number-formatting.md and library-migration.md still treat
  astropy as live and 3.9 as the floor; galactic-potential.md,
  orbital-solvers-and-integrators.md, db-check-and-parity-repair.md,
  api-design-standards.md and ops-scheduling-and-rotation.md need the
  same check; stale 3.9 comments remain in util/settings.py (line 39),
  queue/work.py (247) and cli/generate.py (1376). PERF.88 (the numba
  sector-path integrator) is still open, so say so. Done: the notes and
  comments match OPS.42 and GEN.211.
  Progress (PR #1592, Quality): the notes and two of the comments are
  done. One leftover: the Python 3.9 docstring and workaround
  `reject_lost_option_values` in cli/generate.py is dead on 3.12; Quality
  removes it after the spherical switch merges, because generate.py is in
  the switch.
  Quality, 2026-10-15 (#1698): DOC.21 is done; the
  `reject_lost_option_values` comment in planetgen/cli/generate.py is
  still open here.
  Prerequisites: none. Related: PERF.88.

- [ ] **DOC.45 Reference docs follow the description layer**
  Audit: /mnt/project-files/wording/text-output-audit-2026-10-15.md.
  Owner lane: Foundations Wording (2026-10-15), which writes the code
  for this item. Leftover sweep for docs/html-interface.md, api.md
  (/sections removed, /description added), cli.md (--format),
  system-file-format.md, example-systems.md,
  design/habitability-index.md (PHI-6 naming) and testing.md. Each build
  item updates the page it changes in its own PR; this is what is left,
  part of the DOC.18 to DOC.44 sweep at the end of Phase 1.
  Prerequisites: UX.105, GEN.218, GEN.219.

- [ ] **DOC.47 Rebuild docs/design/architecture.html from the new architecture.md**
  Found by Quality after DOC.27 and DOC.28 (PR #1628) rewrote the code
  map and flows in docs/architecture.md for the shell grid. The
  illustrated snapshot `docs/design/architecture.html` is still dated
  2026-10-07 and the README links it, so it shows the cylinder grid and
  the old code map. Done: the page is rebuilt from the new
  architecture.md (same layout, dated the day it is rebuilt), names
  shells, bands and wedges instead of rings, layers and slots, and the
  README link and every `architecture.html` reference in the docs still
  point at it. Owner: Quality (2026-10-15), after the Core and Maps
  changes settle, so the snapshot is not stale on the day it is made.
  Not a bug.
  Prerequisites: none. Related: DOC.27, DOC.28 (done, PR #1628).

- [ ] **DOC.48 Docs: the reach, the halo, the profile rules and the Wolf-Rayet and nebula-contents rules in the design notes**
  Step 14 of the whole-sphere train.
  `docs/design/galaxy-disk-density.md` says the stellar halo is not
  modelled beyond the floor, `docs/design/galaxy-coordinate-system.md`
  section 4 describes only the outline,
  `docs/design/globular-clusters.md` needs its status,
  `docs/design/interstellar-object-rates.md` needs the populations and
  the kind factors, and `docs/design/anomalies.md` points its Wolf-Rayet
  row at the design note. Done: all five describe the halo term, the
  reach versus the outline, the profile, the components, the kind
  factors, the nebula contents, the Wolf-Rayet state and class X, and
  the 96-bit layout as built. Owner: Quality, after the code.
  Built-on note (Core 1, 2026-10-16, 3:17 AM EDT, GEN.235 merged as
  #1784): section 5.9 of docs/design/intergalactic-and-halo-contents.md
  now has a "Built" paragraph and docs/cli.md is updated; this item
  carries the rest.
  Prerequisites: GEN.241, GEN.243.
  Design: [docs/design/intergalactic-and-halo-contents.md](design/intergalactic-and-halo-contents.md)

- [ ] **DOC.49 Docs: point the rate, anomaly, multi-star and nebula notes at the phenomena catalogue**
  Step 9 of the phenomena train. `interstellar-object-rates.md`,
  `anomalies.md`, `multistar-and-compact-systems.md` and
  `nebula-and-asteroid-field-classes.md` say where each kind now lives
  (`phenomena-catalogue-and-fill.md` and the registry). Done: each of
  the four names the catalogue section and the registry row for its
  kinds, and the Wolf-Rayet and compact-remnant rules agree with the
  build. Owner: Quality, after the code.
  Prerequisite: GEN.130. Related: DOC.48.
  Design: [docs/design/phenomena-catalogue-and-fill.md](design/phenomena-catalogue-and-fill.md)

- [ ] **DOC.50 Docs: the schema, the multi-star note and the catalogue say where the members live**
  Step 10 (design note 11 item 10). `database-schema.md` (the two new
  tables and the two new columns), `multistar-and-compact-systems.md`
  7.3 corrected (black hole companion 1.25%),
  `phenomena-catalogue-and-fill.md` sections 3 and 5 corrected, a
  pointer from `interstellar-object-rates.md`. Done: each says what the
  build does. Owner: Quality, after the code.
  Catalogue (Boss 2026-10-15 6:08 PM EDT;
  docs/design/multi-star-systems.md sections 13 to 16): also the schema
  doc for the new columns (`state`, `flags`, `interaction`, the SMALLINT
  `slot`, the `group` and `compact` node kinds, the offset, velocity and
  `decay_age_yr` columns), a pointer from
  `phenomena-catalogue-and-fill.md` section 4 to the catalogue sections
  13 to 16, and the cap of twelve. Runs after the catalogue code and
  audit.
  Prerequisites: GEN.264, TEST.164. Related:
  DOC.49.
  Design: [docs/design/multi-star-systems.md](design/multi-star-systems.md)

- [ ] **DOC.51 Docs: the cell key in the coordinate system and schema docs, the region note moves to built, the API states the cap and the ladder**
  Step 8 (design note section 9 item 8). `galaxy-coordinate-system.md`
  and `database-schema.md` describe the cell key, the column, the index
  and the guard; `spherical-region-queries.md` moves to "built" when
  DB.34, GEN.276 and the nearby item land; the API reference states the
  cap and the ladder. Done: each says what the build does. Owner:
  Quality (2026-10-15). Boss (2026-10-15 5:32 PM EDT): "if it's not
  ready in time before the spherical sector layout is, then we'll skip
  it and move to the next sequence until it's finished research and we
  can plan it." The research is finished and merged (#1659), so stage C
  is live.
  Prerequisites: GEN.276, NAV.60. Related: NAV.61, TEST.165.
  Design: [docs/design/spherical-region-queries.md](design/spherical-region-queries.md)

- [ ] **DOC.52 A docs index page: links to every HTML page under docs/, each with a description of what it contains**
  Boss (2026-10-15 7:01 PM EDT): "Add a TODO item in DOCS to build an
  index with links to pages with descriptions of what each page contains
  under docs to link the HTML together." Today the HTML pages under
  docs/ (the execution plan, the version plan, the TODO tree, the TODO
  reference, the architecture page, and any page added later) do not
  link to each other except by hand. Build `docs/index.html`: one entry
  per HTML page with its title, a one-or-two sentence description of
  what it contains and who it is for, and a link; the generated TODO
  pages and the execution plan gain a header link back to the index, and
  the index links the Markdown design notes and the GitHub issue list
  for the goals (docs/plan/issue-groups.md). Keep it generated where the
  page is generated (a small script reads each page's title and a
  description line, so a new page cannot be forgotten, with a test that
  fails when an HTML page under docs/ is missing from the index). Done:
  docs/index.html exists, every HTML page under docs/ is listed with a
  description and links back to it, and the test passes. Owner: Quality
  (2026-10-15).
  Prerequisites: none.

- [ ] **DOC.53 The tuning module header names the rule for a new signed constant: the fuzz tests accept only a dict or a positive "fall"**
  Finding from Core 2 (2026-10-16, 9:17 PM EDT, kind: documentation
  gap): the tuning fuzz tests reject negative top-level constants and
  unordered pairs, so a new signed constant must be a dict or a positive
  "fall". Done: one or two lines in the header of
  src/planetgen/tuning.py (and the tuning notes in docs if there are
  any) say so, so the next lane that adds a signed constant does not
  find out from a red test. Owner: any Core lane, with the next tuning
  change.
  Prerequisites: none. Related: DOC.52.

- [ ] **DOC.54 Admin help: how to edit, regenerate, lock, approve and undo**
  Boss (2026-10-16, 1:07 PM EDT), headline of 8.6 (ADM.57): the admin
  has complete control; the world adjusts to what the admin wants, after
  the admin approves the generation plan. A help page for the admin
  chain: the plan and its states, approving and rejecting operations,
  the keep-or-regenerate matrix, locks, density edits, undo and the
  history page, with a screenshot of the A to M example. Done: the page
  in the help section (DOC.6 template) and a link from every approval
  and regen screen.
  Needs the help framework (first child of DOC.6, built in 8.6 phase 2),
  so the headline can close (2026-10-16 review).
  Prerequisites: TEST.175, ADM.68, DOC.57. Related: ADM.57, DOC.6.
  Design: [docs/design/edit-everything.md](design/edit-everything.md)

- [ ] **DOC.55 Help: roles, edit requests and how an approval works**
  Boss (2026-10-16, 1:43 PM EDT), headline of 10.0 (USR.25): three
  roles, User, Moderator and Admin; users submit edit requests,
  moderators approve user requests, admins approve moderator requests,
  and an admin's own edits need no approval (only the preview and
  confirm of ADM.66). Build on the accounts of 8.2 and the change plans
  of 8.6. A help page for the roles, what each can request and approve,
  how to submit a request and follow it, and what moderators and admins
  see. Done: the page in the help section (DOC.6 template).
  Prerequisites: TEST.178. Related: USR.25, DOC.6.
  Design: [docs/design/roles-and-requests.md](design/roles-and-requests.md)

- [ ] **DOC.56 Backups, restore and repair: the deployment guide and the help page**
  Boss (2:00 PM EDT): "we also have to build a way to use that data so
  build that out and let's make this targetted for v9.0 so incorporate
  the chain for these (restore / repair using backups) with that due
  date. Settings page is v8.1". Target 9.0 (the backups and repair
  chain, umbrella ADM.91); the settings page and the settings backups
  stay in 8.1 (ADM.83 to ADM.85). All backup pieces are parts of one
  `.tar.bz2` archive (ADM.87); retention is `settings.backups_kept`
  (default 20, 1 to 200). `docs/deployment` and the help section
  (DOC.6): what a backup holds and does not (the galaxy database is not
  in it), where the archives live and their retention, how to take one,
  how to restore the settings, the control database or the parity data,
  how to repair the galaxy from parity and from an archive, what the
  guard refuses, and the `mysqldump` requirement. Done: the guide and
  the help page, linked from the settings page and the repair page.
  Prerequisites: TEST.180. Related: ADM.91, DOC.6.

- [ ] **DOC.57 Help section framework: page template, index, '?' links, help search and the link check**
  Split from DOC.6 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work). The web interface serves a Help
  section of static pages (no database, no JavaScript needed to read
  them): one shared template that matches the site shell, an index page
  that groups pages by feature, a '?' link in the header of every
  interface page that opens the matching help page (a route-to-page map
  in docs/help), a search box over the help text, a plain page of the
  current version, and the pages written as Markdown in docs/help/ and
  built into HTML by the update script (served at /help, as decided by
  Boss on 2026-10-10). Ship it with one placeholder page. Done: /help
  renders the index and a page in both themes, the '?' link works on
  every page, the update script builds docs/help, and a test checks
  every link inside the help section resolves.
  Template, index, the "?" link, search, the version page, the link
  check and the route-to-help coverage test with an allowlist of routes
  still without a page (DOC.7 to DOC.16 burn it down); a feature item is
  not done until its help text is, which needs the framework before 8.7
  (2026-10-16 review).
  Prerequisites: none. Related: DOC.6.

- [ ] **DOC.58 Help coverage test: every user-facing route has a help page and every help page has a route**
  Split from DOC.6 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work). A test lists every user-facing route
  (the route map the app already has) and fails when a route has no help
  page or a help page has no route; admin routes are covered by the
  admin pages. It runs in the gate. Done: the test, green once DOC.7 to
  DOC.16 and DOC.54 exist, with an allow-list file for routes that are
  deliberately undocumented (each entry needs a reason).
  Prerequisites: DOC.57, DOC.7, DOC.8, DOC.9, DOC.10, DOC.11, DOC.12,
  DOC.13, DOC.14, DOC.15, DOC.16, DOC.54. Related: DOC.6.
  The coverage test runs without the allowlist once every page exists
  (2026-10-16 review).

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
    time of day (a whole-sky chart does not need it). Decided (Boss,
    2026-10-10 18:44Z, defaults approved; defaults taken): exact stars
    for V <= 6.5 with tiers built by
    an angle-bin scatter and a background glow for the rest, the
    planet's page saying "sky incomplete" until its tiers ran; limiting
    magnitude 6.5 for every planet with horizon extinction and twilight
    only; physical colour from blackbody, desaturated for faint stars,
    with a "boost colour" option; galactic plane dust A_V 1.0 per kpc;
    the planet's north for the horizon is its spin pole.
    Research (2026-10-09, globular-clusters.md): the globular clusters
    of generated galaxies come from GEN.164.
    Split (2026-10-16, Boss 1:07 PM EDT: look at each item for smaller
    chunks and shared work): this item is now the umbrella and is done
    when its subitems are: VIEW.11, VIEW.12, VIEW.13, VIEW.14. VIEW.2 is
    the whole 'starmap seen from a planet' (photometry, which stars are
    visible, colour, light-time, dust, the band, orientation) and
    sky-view.md section 8 and VIEW.1's own 2026-10-09 research note
    already cut it into stages that 'each ship alone, in that order'.
    Scope after the review (2026-10-16 review): this item is the
    umbrella of the photometry library (a), the sky catalogue (b), the
    galactic band (c) and the orientation and horizon (d), cut along the
    stages of sky-view.md section 8.
    Prerequisites: VIEW.11, VIEW.12, VIEW.13, VIEW.14.

  - [ ] **VIEW.3 Render the view as a PNG, with constellations**
    Boss: "when
    it does that it will generate a PNG and it will generate
    constellations." Done: VIEW.2's view is drawn to a PNG (a sky
    projection, star size and colour by apparent brightness), and the
    brighter stars are grouped into constellations with lines and names
    from VIEW.4, stored so a planet keeps the same constellations each
    time. Blocked on VIEW.2's research pass. Question (decided below):
    whole-sky
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
    Decided (2026-10-10): Boss (19:16Z, Accept all): the whole-sky PNG
    comes first; seeded constellations are stored per system; the PNGs
    are cached.
    Split (2026-10-16, Boss 1:07 PM EDT: look at each item for smaller
    chunks and shared work): this item is now the umbrella and is done
    when its subitems are: VIEW.15, DB.64, VIEW.16, VIEW.17. VIEW.3 is a
    picture and a stored data set.
    Scope after the review (2026-10-16 review): the umbrella of the
    whole-sky PNG (a, needs only children a, b and c of VIEW.2), the
    constellation table (b, Core 1, in phase 1), the constellations (c)
    and the overlay (d).
    Prerequisites: VIEW.15, DB.64, VIEW.16, VIEW.17.

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
    VIEW.3. Question (decided below): what counts as a source list
    (licensing of
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
    Decided (2026-10-10): Boss (19:16Z, Accept all): constellation names
    use the name codec (the 2026-10-07 plan stands); the licensing,
    transliteration and uniqueness questions are moot.

- [ ] **VIEW.6 Settle the handedness of the generated galaxy before mapping real sky coordinates**
  The generated galaxy rotates counterclockwise about +Z and the real
  one clockwise about the north pole. Settle which frame sky charts and
  the neighbour galaxies use. Decided (Boss, 2026-10-10 18:44Z, defaults
  approved; default): keep the
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

- [ ] **VIEW.10 A HYG calibration test: a generated sky has about 9,000 stars at V <= 6.5**
  Pillow is already installed through scikit-image (pillow 11.3.0 below
  Python 3.10, 12.3.0 above) and no repo code imports it; declare it if
  Pillow is used (open question for Boss, default: SVG overlay, Pillow
  only for a self-contained PNG with text). A calibration test that
  reads HYG from a developer machine (not committed).
  Rewritten (2026-10-16 review): the "declare pillow" half is dropped
  (not declared; no repo code imports it; revisit only if VIEW.3 adds
  burned-in text); the calibration test moved to phase 2 after the
  photometry library and the sky tiers, which it needs.
  Prerequisites: VIEW.11, VIEW.7.
  Design: [docs/design/sky-view.md](design/sky-view.md)

- [ ] **VIEW.11 Photometry library: apparent magnitude, bolometric correction, blackbody colour and the extinction law**
  Split from VIEW.2 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work). Pure array functions in
  planetgen/physics/photometry.py with no database: the
  bolometric-correction table (replaced by VIEW.9),
  `apparent_magnitude`, the blackbody colour table (desaturated for
  faint stars, with the 'boost colour' option), and the extinction law
  applied to a magnitude and a colour (A_V 1.0 per kpc in the galactic
  plane, the horizon extinction left to child d). Done: the library with
  tests against the numbers of sky-view.md (Sun -0.07, A_B/A_V 1.3226,
  Sun B-V to 5,778 K), and the HYG calibration test of VIEW.10 run
  against it.
  Prerequisites: VIEW.1, VIEW.9. Related: VIEW.2.

- [ ] **VIEW.12 Sky catalogue: the stars an observer sees, with direction, brightness and light-time**
  Split from VIEW.2 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work). planetgen/galaxy/sky.py: for an
  observer (a planet, a moon or a position) arrays of position,
  luminosity, temperature and uid from `bright_stars`, generated stars
  and the sky tiers of VIEW.7, the vectorised first-order `sky_vectors`
  and the retarded-time step (age and luminosity at the time the light
  left; the position shift is under a tenth of a pixel), the 6.5
  magnitude cut, and the 'sky incomplete' state until the planet's tiers
  have run. Done: `sky_from(observer)` returning about 9,000 stars at V
  <= 6.5 for a seeded planet, with tests for the cut, the light-time
  step and the incomplete state.
  Prerequisites: VIEW.11, VIEW.7. Related: VIEW.2.

- [ ] **VIEW.13 Galactic band: the glow of the unresolved stars and dust along each line of sight**
  Split from VIEW.2 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work). planetgen/galaxy/sky_band.py: the
  ray-march over the galaxy density for the unresolved glow, cached by
  quantised position, with the dust model of VIEW.8 (the galactic law
  and `nebulae.extinction_av`) applied along the ray. Done: the band as
  an image array for an observer, compared against the resolved counts
  of child b, with a cache test.
  Prerequisites: VIEW.11, VIEW.8. Related: VIEW.2.

- [ ] **VIEW.14 Orientation and horizon: the planet's own sky, time of day and twilight**
  Split from VIEW.2 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work). The observer's frame from the
  planet's spin pole (north), axial tilt, rotation phase and latitude
  (GEN.104), the horizon frame (east is North x Up,
  navigation-frames.md), twilight and horizon extinction, and the
  planet's own star's glare; the chart's own l and b are labelled as
  measured in the generated frame (VIEW.6). Done:
  `horizon_view(observer, time)` with tests for a pole, the equator and
  a tilted planet.
  Prerequisites: VIEW.12, VIEW.13, VIEW.6. Related: VIEW.2.

- [ ] **VIEW.15 Whole-sky PNG: Hammer-Aitoff projection, written with numpy and zlib, cached and served by the API**
  Split from VIEW.3 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work). Draw the stars and band of VIEW.2 as
  a Hammer-Aitoff chart (star size and colour by apparent brightness,
  peak intensity and sigma from sky-view.md section 7, gamma-encoded
  after summing), write the PNG with numpy and `zlib` (no new
  dependency), cache it by system uid, epoch bucket, projection, size,
  limiting magnitude and version key, and serve `GET
  /api/systems/<id>/sky.png?proj=hammer&mag=6.5` with the usual scopes.
  Done: the route returns a cached PNG for a seeded system, a second
  request is a cache hit, and a test checks star positions against the
  projection.
  Prerequisites: VIEW.11, VIEW.12, VIEW.13. Related: VIEW.3.

- [ ] **VIEW.16 Constellations: the seeded partition over the brightest stars, stored per system, named by the codec, with Bayer labels**
  Split from VIEW.3 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work). The algorithm of sky-view.md section
  6 (seeds 20 degrees apart, one minimum spanning tree per region, 12
  stars each, over the 500 to 900 brightest stars), deterministic from
  the system's seed so a planet keeps its figures, stored through child
  b, named with the codec names of VIEW.4, and the Bayer letters
  assigned by apparent magnitude (derived, not stored). Done: the
  library and its stored result for a seeded system, with tests for
  determinism, stability when the limiting magnitude changes slightly,
  and unique names.
  Prerequisites: DB.64, VIEW.4, VIEW.12. Related: VIEW.3.

- [ ] **VIEW.17 Constellation overlay: SVG lines and names over the sky chart, and the planet-sky card**
  Split from VIEW.3 (2026-10-16, Boss 1:07 PM EDT: look at each item for
  smaller chunks and shared work). Draw the figures and their names as
  an SVG overlay on the PNG of child a (the PNG stays cacheable and the
  labels sharp), with a toggle, and register the planet-sky card with
  ADM.44's card renderer. Done: the overlay on the planet's page and the
  card, and a browser test that toggles the figures and finds a named
  constellation.
  Prerequisites: VIEW.15, VIEW.16. Related: VIEW.3.

## POP: Population and politics

Population, species and polities (POP.1 to POP.6) shipped in PRs #169 to
#184.

- [ ] **POP.11 Scatter space vessels and stations for the races in their territories**
  GitHub issue [#1226](https://github.com/dwhagar/planetGen/issues/1226)
  (Boss, 2026-10-12 23:27Z): "Also part of population and tech, have the
  ability to scatter space vessels and stations for appropriate races
  within their territories." Done: the population pass places vessels
  and stations for species whose tech level allows them, inside their
  territory, as stored objects like the user objects of USR.10. Which tech
  level allows vessels and which stations is set from the tech-level
  design (docs/design/Technological Assessment.md and the "Tech levels"
  section of population-and-politics.md, POP.8, done in PR #1422). Owner: not assigned.
  Prerequisites: DB.63, USR.29. Related: POP.7 (done, PR #1461), POP.10
  (done, PRs #1517 and #1524), USR.10.

- [ ] **POP.13 A species editor: create and edit species and their environment limits**
  Boss (2026-10-15 04:15Z): "Add a TODO item, I want the system to allow
  one to list all races that would be able to inhabit that environment.
  That means all species must also have that information. Add another
  TODO item for species editor." Done: an admin page lists the species
  and edits one or creates one: name, homeworld, life chemical, build,
  climate, size, the six tech indexes and tech level (`tech_source`
  becomes admin, the column POP.9 added for this), and the environment
  limits of POP.12 with the generated defaults shown beside each so an
  admin can see what was changed; saving checks the values (ranges, a
  homeworld that exists, a unique name) and refreshes the "could live
  here" lists; an admin-made species is marked as such and survives a
  reseed like other admin edits. Same shape as the facility type editor
  (POP.10, done) and the same admin patterns (ADM.9). Owner: UX
  for the page, Generate for the checks. Phase 2: nothing
  in Phase 1 needs it once POP.12 stores the data; Boss can pull it
  forward.
  Prerequisites: ADM.92. Related: POP.10 (done, PRs #1517 and #1524),
  POP.9, ADM.9.
  Moved to phase 2 with the shell as its prerequisite, so it is built on
  the shared table, form and confirmation and is not missed by the audit
  and the conversion (2026-10-16 review).

- [ ] **POP.14 A race is only placed where it can survive: a placement check and a consistency check over the data**
  Progress (2026-10-15): Generate built it and saved it on the branch
  `claude/pop14-held` (not merged, not started on main); the Core lane
  starts from that branch. It waits on DB.30 for storage.
  Boss (2026-10-15 04:18Z): "Also add another POP item that we ensure
  consistency that if we put a race on a planet it is within the limits
  of what can survive there." Done: (1) one check, built on the species
  limits and the scoring function of POP.12, that every way of putting a
  species on a planet or moon calls before it writes: generation
  (homeworld and any later settlement the scatter makes, POP.11), the
  species editor (POP.13), the API and any command or import; it refuses
  a placement outside the species' limits and says which domain fails,
  in plain in-universe wording; (2) a consistency check over the stored
  data (a report command and a test over a generated galaxy) that lists
  every species placed outside its limits; any it finds are fixed by
  regenerating the limits or the placement, with no shim and no
  tolerance flag. How it fits: the POP.12 default takes a species'
  limits from its homeworld, so every species passes its own homeworld
  by construction and this check matters for every other placement; an
  admin who edits limits in the species editor (POP.13) must not be able
  to leave a species outside its homeworld's conditions, so the editor
  runs this check on save and shows the failing domain. A placement that
  is only acceptable with equipment counts as outside the limits unless
  the placement says it is a habitat or vessel (default: only unaided
  survival counts for a planet or moon surface).
  Prerequisites: DB.30. Related: POP.13, POP.11, POP.10 (done, PRs #1517
  and #1524).
