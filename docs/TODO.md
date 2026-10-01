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

NAV, DB, API, OPS, DOC and POP have no open items today.

## Background

The long-term goal is a fully populated galaxy (every sector, system,
planet, moon and belt) stored in MySQL and browsed through a web
interface served from Apache. The original roadmap phases are all done:
object-graph serialization, the relational schema and migrations, CLI
tools writing to the database, lazy galaxy-scale generation from a
density skeleton, and the Flask API (`src/html/api/`) with server-rendered
pages served by the same app (`src/html/web/`): Galaxy, Sector and System
maps, search, NAV, admin auth and wiki publishing.

## Plan: what to do first

Boss (2026-10-01): "bugfixes and security are the two biggest concerns."
Bugs come first, then security; everything else waits behind them. The
research behind the security order is the design doc linked at the top
of the SEC section.

1. **Bugs, in this order:**
   1. MAP.45 (objects drawn outside the sector's wireframe), then
      MAP.46. MAP.45 goes first because it may mean stored positions
      are wrong; find out before the next full regenerate.
   2. Done: the bright stars on the Galaxy Map (MAP.47, MAP.48), the
      wedge lines past the galaxy's edge (MAP.43) and the generated
      systems that were hard to find (MAP.37).
   3. MAP.49: planet orbits drawn inside an asteroid belt on the System
      Map.
   4. MAP.50: names running off the edge of the map.
   5. The drill-down rework, built together: MAP.17 and MAP.19 (no free
      camera; big wedge, layer and region picks), with MAP.18 (dim
      everything but the hovered pick) and MAP.44 (lines kept to the
      zoomed block), then MAP.26 (open at the sector; the map's own Back
      and Forward). It needs Boss's answers to MAP.17's and MAP.19's
      open questions (what a slice and a region are) before it starts.
2. **Security, in this order** (login blocking first):
   1. SEC.20: log every failed and locked login with its address.
   2. SEC.1: the per-IP lockout in the control database, with SEC.21
      (the per-username backoff moved into the same table).
   3. SEC.23 (bug): a wrong current password on `/account` isn't
      counted.
   4. SEC.22: trusted-device cookie, so a lockout can't shut the real
      admin out.
   5. SEC.24 (common and breached password blocklist), then SEC.25
      (check the hashing cost, re-hash on login).
   6. SEC.26: two-factor sign-in for admins.
   7. SEC.27: the fail2ban recipe in the deployment docs.

   SEC.20 and SEC.23 are small and touch different files from the bugs,
   so they can run alongside step 1 if Boss wants.
3. **Database calls, the next update after bugs and security** (Boss,
   2026-10-01), in this order:
   1. PERF.12: check the schema once per process during generation.
   2. PERF.13: write each sector in batches.
   3. PERF.14: reserve a sector's names in bulk, safe with several
      writers at once (needed before any parallel generation).
   4. PERF.15 (fewer queries per web page), PERF.16 (full-text name
      search on whole words) and PERF.17 (a time limit on web
      statements). These touch the web read path, not generation, so
      they can run alongside 1 to 3.
4. **After that**, as before:
   - MAP.30 (slab list beside a 3:4 map) and bookmarks (MAP.23, which
     finishes the NAV page's map picks, MAP.22).
   - Sector Map stars as points of light (MAP.15).
   - The rest of generation at scale (PERF.1): the per-sector density
     stats (PERF.11) and speed records (PERF.10) feed the estimates
     (PERF.3, PERF.9); the work queue (PERF.8) comes before parallel
     generation (PERF.7).
   - Admin editing (ADM.1) starts with the validate module (ADM.5);
     user accounts (USR.1) start with roles (USR.2).
   - View from a planet (VIEW.1) waits on a research session with Boss,
     except the constellation names (VIEW.4).

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

Bugs Boss found on the Galaxy Map (2026-10-01): "bug fix, wedge lines
should not extend past the boundary of the galaxy. bug fix, bright stars
do not display past 500 seconds scale (1 block = 81 sectors across). bug
fix zooming reveals stars are being drawn but it takes a while to load,
another bugfix, it's really hard to find the generated star system on
the map, so everything not yet filled should be more transparent by a
lot with a much higher contrast." They were MAP.43, MAP.47,
MAP.48 and MAP.37, all fixed.

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
  doc; each subitem names its section. The nested ladder (MAP.28,
  `stellarObjects/galaxyDrill.py`), the stage contents API (MAP.29,
  `GET /api/galaxy/stage`), the stages (MAP.16,
  `static/galaxystages.js` and `static/galaxystageview.js`, with stage
  URLs `/galaxy?slab=`, `?at=`, `?sector=<designation>`), the Sector Map
  pick mode (MAP.21), the address bar (MAP.24, `/galaxy/locate`), the
  course overlay (MAP.27, `/galaxy?course=<from>,<to>`), generating from
  the map (MAP.20) and the "Show on Galaxy Map" links (MAP.25) have
  shipped. The old free camera stays behind the map's Free look button
  until MAP.17 removes it; Boss settled decision 2 (section 11) that
  way.

  - [x] **MAP.16 Drill-down stages**
    Done in 7.44.0 (PR #171); kept as the parent of its bug.
    - [ ] **MAP.17 (bug) No free camera; drill down from a top-down view by wedge, slice and block**

      Design: [docs/design/galaxy-drilldown-navigation.md](design/galaxy-drilldown-navigation.md), sections 5, 10 and 11

      Boss (2026-10-01): "bugfix, remove the
      ability to free form select a point to center on. Instead, we'll
      start at a top-down view. The user will select the wedge (quarter)
      aligned with the wedge lines, that we then zoom in on. From there
      the user select a slice of blocks that fit the current zoom level.
      It is selected by the mouse and it should be clearly highlighted the
      slice the user is going to click on and make sure the block size is
      such that the user can easily operate it via the correct method
      (responsive design, if a touch sized screen make sure the user can
      easily tap with a finger, but if it's a computer then the user is
      probably using a mouse). Once selected the user can choose any block
      in that slice. Once selected, the user zooms into that block, they
      can then see the next group of slices and the process repeats. This
      occurs until we get to the sector level. At no point can the user
      free rotate the map anymore." This settles decision 2 of the design
      doc (section 11, "Free camera"): no free look and no free rotation,
      not even the default's drag-rotate inside the 3D stages. Today the
      old free camera (click any point to center on it, drag to rotate)
      stays behind the map's Free look button (added with MAP.16, PR
      #171), and the 3D stages (1, 3, 5, 7, section 5.1) can be
      drag-rotated. Done:
      - The Free look button, free-camera picking and every drag-rotate
        are removed; the map can't be rotated at any stage.
      - The map opens top-down on the whole galaxy; the user picks a
        wedge (a quarter, its edges on the wedge lines) and the map zooms
        in on it.
      - At each level the user picks a slice of blocks sized to the
        current zoom, then any block in that slice, and the map zooms
        into that block and shows its slices; this repeats down to the
        sector level, where a click opens the sector.
      - The slice (and then the block) under the pointer is clearly
        highlighted before the click.
      - Blocks and slices are sized for the input: big enough to tap with
        a finger on touch screens (`pointer: coarse`), mouse-sized on a
        computer (the Responsive Web Design Standards' 44-48 px touch
        targets, under the UX section).
      - The design doc's sections 5 and 10-11 are updated to match, and
        breadcrumb, Back, stage URLs (`/galaxy?slab=`, `?at=`), bookmarks
        (MAP.23), the address bar (MAP.24) and the NAV course (MAP.27)
        keep working with the new steps.

      Builds on the fixed map bugs (wedge lines stop at the galaxy's
      edge, which matters more once wedges are what the user picks, and
      the filled blocks' contrast applies to the slices and blocks).
      Open questions: how the "wedge (quarter)" maps onto the design's
      nested ladder (243 -> 27 -> 3 ->
      1, section 3): is the first pick always one of four quarters, or
      one of the master wedges at that ring? What a "slice" is: a ring
      band, a layer (the current 3D stages pick a slab, a layer), or a
      row of blocks along the wedge? With no 3D view, how does the user
      pick a layer above or below the galactic plane (a side view, a
      layer list, or slices that run through the disk's thickness)? What
      replaces the free camera for MAP.20's neighborhood generate and
      MAP.22's NAV picking, which may want an arbitrary point?

    - [ ] **MAP.18 (bug) The block under the pointer is too hard to see from above**
      Boss (2026-10-01): "selection of blocks from top-down doesn't show
      highlight well enough, needs to be more obvious, when selecting
      everything but what your mouse is over should be dimmed". Today, in
      the top-down stages (2, 4, 6, 8), hovering a block only draws an
      outline around it (`galaxystageview.js`, `applyHover` ->
      `outlineBlock`), and nothing else changes; the 3D stages already
      fade the other slabs (`OTHER_SLAB_FADE`). Done: while the user is
      picking, everything except the block (or slice) under the pointer
      is dimmed, and the hovered one stands out clearly, in both themes,
      with the same look for a first tap on touch screens and for
      keyboard focus; the dimming clears when the pointer leaves the map.
      Goes with MAP.17 (the wedge, slice and block picks all need this
      highlight); unfilled blocks are now drawn much more transparent,
      which the dimming must still read against. Open
      questions: how dim the rest goes (say 25-40% opacity)? Does the
      hovered block also brighten or get a thicker outline? Does it
      animate (a short fade) or switch instantly, and does
      `prefers-reduced-motion` turn the fade off?

    - [ ] **MAP.19 (bug) Pick big wedges, not single blocks, between zooms**
      Boss (2026-10-01): "When selecting regions in the galactic map,
      selecting bigger arcs between zooms. Even if we can select a
      single block, that might be too small. So, have the user select
      large wedges and zoom into that wedge. Same plan: quadrant, layer,
      region, layer, region, ... , layer, sector. This might be a
      duplicate or need to enhance a bugfix." It refines MAP.17 (the
      wedge, slice and block drill-down) rather than duplicating it, and
      the two should be built together: MAP.17 says what to remove
      (free camera and rotation) and that picks go wedge, slice, block;
      this item sets the size of each pick. Today each top-down stage
      picks one block of the 243 -> 27 -> 3 -> 1 ladder (MAP.28), which
      can be a small target. Done: every pick is a large wedge (an arc
      of the ring band in view) that the map then zooms into, following
      the ladder quadrant, layer, region, layer, region, and so on, down
      to layer and sector; each pick stays large on screen at every
      stage (and finger-sized on touch, as MAP.17 asks); the
      highlight follows MAP.18. Open questions: how a "region" maps
      onto the nested ladder (one of the 3x3 = 9 children of a block, a
      wider arc across several blocks, or a new level between); whether
      the layer pick alternates with every region pick, as the ladder
      reads, or only where the disk is thick enough to have several
      layers; and whether MAP.17's "slice" is the same as this
      item's "layer" (if so, MAP.17's open question on what a slice
      is is settled).

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

  - [x] **MAP.25 "Show on Galaxy Map" links**
    Done (PR #188); kept as the parent of its bug.
    - [ ] **MAP.26 (bug) "Show on Galaxy Map" should open at the sector, and the map needs its own Back and Forward**

      Design: [docs/design/galaxy-drilldown-navigation.md](design/galaxy-drilldown-navigation.md), section 8

      Boss (2026-10-01): "when the
      user clicks "Show on Galaxy" It should be zoomed in to the sector
      level of the slice that we can see that sector, and back and forward
      buttons to travel ones own history on the map display". Today
      `/galaxy?sector=<designation>` opens the sector's stage 8
      (`galaxystages.parseStageQuery`), and the sector and system pages
      and search results link to it (MAP.25), but the NAV result's link
      opens the course overlay (`/galaxy?course=<from>,<to>`, MAP.27)
      over the galaxy. The map has no Back or Forward of its own; moving
      between stages relies on the breadcrumb and the browser's Back
      button. Done:
      - Every "Show on Galaxy Map" link (sector, system and search pages,
        the NAV result) opens the map zoomed in to the sector level of
        the slice holding that sector, with the sector in view and
        highlighted.
      - The map display has Back and Forward buttons that step through
        the user's own history on the map (each stage or position they
        visited), separate from but consistent with the browser's
        history and the stage URLs.

      This follows the wedge, slice and block drill-down of MAP.17
      (the "slice" here is MAP.17's slice), and builds on MAP.25's
      links. Open questions: for the NAV result, which sector is shown
      (the start, the destination, or a view that fits both, as
      MAP.27's course does today)? Is the map's history its own list or
      the browser's history (`history.pushState` per stage) with the
      buttons calling `history.back()`/`forward()`? Does it survive a
      page reload or a visit to a sector page and back? How far back
      does it go?

- [x] **MAP.3 A bigger Galaxy Map with controls underneath**
  Done in 7.55.0 (PR #178); kept as the parent of its subitem.
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
    change. Goes with MAP.17 (once the drill-down picks wedges and
    slices, the left-hand list holds whatever the current stage offers)
    and follows the Responsive Web Design Standards notes in the UX
    section (size classes, container queries). Open questions: is 3:4
    width to height (taller than wide) or height to width (4:3, wider
    than tall, close to today's 75vw height)? At what width does the
    list move to the left (the standards' expanded class, 840 px and up,
    or whenever a readable list column fits)? Does the rest of MAP.3's
    controls row stay under the map, or join the list on the left?

- [x] **MAP.5 Galaxy Map rework**
  Done (the pixel-sized mega-blocks plan Boss approved in the project's
  `galaxy-megablocks/report.md`); kept as the parent of its bugs. The
  map's code is `static/galaxyprisms.js`, `static/galaxymap3d.js`,
  `lib/galaxymap3d.py` and `stellarObjects/galaxyGeometry.py`.
  - [x] **MAP.36 One solid of blocks for filled and unfilled sectors**
    Done in 7.25.0 (PR #142).
  - [x] **MAP.42 Wedge lines from the center**
    Done in 7.9.0 (PR #120).

    - [ ] **MAP.44 (bug) Wedge lines and ring circles run far past a zoomed-in block**
      Boss (2026-10-01): "when zoomed in to a sector block it the radial
      lines should only be contained within the block and just a little
      beyond, and that should also be cylindar boundaries about the
      same." Today the wedge lines (`galaxymap3d.js`, `updateWedgeLevels`
      and `clipWedgeLevel`) are clipped only to the view ball around the
      camera's target, so when the drill-down (MAP.2) is zoomed in to a
      block they still run across the whole view, and the cylinder
      (ring) boundaries do the same. Done: zoomed in to a block, the
      wedge lines and the ring boundaries are drawn only across that
      block and a small margin past its edges, at every drill-down stage;
      zoomed out, they stop at the galaxy's edge as now.
      Open questions: how big "a little beyond" is (a fixed fraction of
      the block, such as 10-20%, or a few pixels on screen)? Do the lines
      fade out across the margin or stop sharply? Does the same apply to
      the 3D stages while MAP.17 still has them, or only to the
      top-down ones?

- [x] **MAP.11 Every kind of phenomenon on the Sector Map, clickable**
  Done in 5.51.0 (PR #81); kept as the parent of its bugs. The Sector Map
  is `html/lib/starmap.py` and `static/sectormap.js`.
  - [ ] **MAP.45 (bug) Rogue planets (and maybe other objects) drawn outside the sector's wireframe**
    First in the bug order (see the plan at the top): if stored
    positions are wrong, that needs to be known before the next full
    regenerate.
    Boss (2026-10-01): "rogue planets (and
    probably other objects) are shown outside the wireframe of the
    sector, so one of them is wrong. when you tackle this one, do a test
    where the rogue planets are a bright color and the background a dark
    color so you can see the distance. This is at the sector level".
    On the Sector Map (`html/lib/starmap.py`, `static/sectormap.js`) the
    wireframe is the sector's cylindrical grid cell (`_outline_data`),
    and the phenomena come from `queryDb.phenomena_near_sector`. Leads
    to check, not yet confirmed: that query takes every phenomenon whose
    sphere reaches the sector's bounding sphere (sized for the old cube,
    `edge_pc * sqrt(3) / 2`) plus every one generated with this
    sector "wherever it sits", so some outside points may be expected
    neighbors; or the phenomena's positions and the wireframe use
    different frames (the cell is rotated to the galaxy frame, and
    phenomena positions are galaxy-placed). Done: find which side is
    wrong (the wireframe, the object positions, or which objects are
    picked) and fix it, so every object generated in a sector is drawn
    inside its wireframe and anything shown from a neighboring sector
    reads as outside on purpose; check every phenomenon type and the
    star systems, not just rogue planets. As Boss asks, the fix includes
    a visual test that draws the rogue planets in a bright color on a
    dark background so the distance past the boundary is easy to see,
    plus an automated check that a sector's own objects fall inside its
    cell (`galaxyGeometry`'s `sector_address_at` giving back the
    sector's own address). Open questions: should nearby phenomena from
    other sectors still be drawn (they are on purpose today, per the
    map's hint "near this sector"), and if so, how are they told apart
    from the sector's own (dimmer, outside-only, or a toggle)? If the
    stored positions turn out wrong, do existing galaxies need a
    migration or a regenerate?
  - [ ] **MAP.46 (bug) Rogue planets are hard to find on the Sector Map**
    Boss (2026-10-01): "in sector view make sure rogue planets can be
    easily located." Today `static/sectormap.js` draws a rogue planet as
    a dim, dark-purple textured sphere (`roguePlanet`: core `#6b5a8a`
    fading to `#2a2438`, glow `#7d6aa8` at strength 0.8), "a dim,
    starless world lit only by its own internal heat", which nearly
    vanishes against the dark scene. Done: every rogue planet in a
    sector is easy to spot at the default zoom and when zoomed out, in
    both themes, without looking like a star; the sector page's list of
    its contents can point at each one on the map. Goes with MAP.45
    (rogue planets drawn outside the wireframe), whose bright-color test
    makes the same objects visible, and MAP.15's point-of-light style.
    Open questions: what makes them findable (a marker or ring around
    each, a brighter but still cool color, a label, a "highlight rogue
    planets" toggle, or a list that flies the camera to each)? Does the
    same apply to other dark objects (quiescent black holes, interstellar
    comets)?

- [x] **MAP.14 Bright stars on the Galaxy Map**
  Done in 7.42.0 (PR #160; it never had a number); kept as the parent of
  its bugs.

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

- [ ] **MAP.49 (bug) The System Map shows planet orbits inside an asteroid belt**
  Boss (2026-10-01): "the system view still is showing orbits of
  planets inside the orbit of an asteroid belt." Generation already
  keeps planets clear of belts (a planet after a belt keeps 5 Hill
  radii clear of its outer edge, and belts can't overlap the next
  planet out, per the CHANGELOG), so the first lead is the drawing:
  `html/lib/systemmap.py` places every orbit on a shared log scale
  (`_radial_px`), but `_belt_band_px` draws the belt's ring width
  linearly (`radius_px * (upper - lower) / distance`) and clamps it to
  9-40 px, so on the log scale the drawn band can reach past a
  neighboring planet's orbit. Second lead: systems stored before those
  generation fixes, which a regenerate would clear. Done: the belt's
  ring is drawn from its inner and outer edges mapped through the same
  log scale as the orbits (with a minimum visible width that never
  crosses a neighbor's orbit), on the System Map and in the 3D system
  view if it draws belts too; a test checks that no drawn planet orbit
  falls inside a drawn belt band for systems with belts; and a check
  over the stored data reports any planet whose orbit really is inside
  a belt (if there are any, that is a generation bug to fix too).

- [ ] **MAP.50 (bug) Names run off the edge of the map**
  Boss (2026-10-01): "names should not run off the screen edge." Taken
  to mean the name labels on the maps, the System Map first. Lead:
  `html/lib/systemmap.py` draws into a fixed 700 x 700 viewBox and
  places each label with `_label_sides_2d`, which tries below, above,
  right and left to avoid other labels and markers but never checks the
  candidate against the viewBox edges, so a name near the edge (an
  outer planet, a long name) is cut off; `static/systemmap.js`'s
  `layoutLabels` measures the real text the same way. Done: on the
  System Map (and its moon views), a label that would cross an edge
  picks a direction that fits, or is shifted inward, so every name is
  fully on screen at every size class; the Sector Map's and Galaxy
  Map's name labels are checked for the same problem and fixed if they
  have it; a test places labels for bodies near every edge and checks
  each label rectangle lies inside the view. Open question: when no
  direction fits, shift the label inward along the edge, or shorten it
  with an ellipsis and show the full name on hover? Default: shift it
  inward.

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

  - [ ] **PERF.6 Rate-limit SQL calls and make each call do more**
    Boss
    (2026-10-01): "ratelimiting calls to the sql database and seeing if
    we can investigate some way to make our DB calls more efficient, do
    more with less calls without impacting performance." The
    investigation is done (report of 2026-10-01: https://claude.ai/code/artifact/4111c1a8-0d63-4b6c-b109-1b8389538a13). Measured on
    a 40-system sector: 1.0 s generating, 1.9 s saving, through about
    98 single-row INSERTs and 17 SELECTs per system; one `executemany`
    wrote the same moon rows 2.8x faster than one INSERT each. The
    batching work is the subitems below. Decided defaults for the open
    questions:
    - Two separate limits. Generation: a cap on how many sector writers
      run at once (PERF.8's parents), at low priority, not a
      calls-per-second cap (batching already cuts its calls 30 to 50
      times). Web: Flask-Limiter's per-IP limits as today, plus a
      statement time limit (PERF.17).
    - Pools: each generation worker process keeps a pool of 1 or 2
      connections, so the total is the worker count plus the web's 10.
    - Benchmark: `generate.py sector --num-systems 40` five times before
      and after on the same machine; compare the median save time and
      the server's `Com_insert`/`Com_select` counts.

    - [ ] **PERF.12 Check the schema once per process during generation**
      `_db.get_connection(ensure_schema=True)` replays all of
      `schema.sql` (about 122 statements) on every checkout, and galaxy
      fill checks out twice per sector (`generate.py` `_fill_context`,
      `_db.save_sector`). Done: generation checks the schema once per
      process (`open_write`, or a per-process flag), and a filled sector
      sends no schema statements.

    - [ ] **PERF.13 Write each sector in batches**
      `insert_sector` and `insert_star_system` in `_db.py` write every
      system, star, planet, moon, belt, comet, composition row, life
      paragraph and spectrum value with its own INSERT. Done: a sector's
      rows are built in memory and each table is written with one
      `executemany` (multi-row INSERT) in the sector's transaction; row
      ids come from a small id-block table (one
      `UPDATE ... LAST_INSERT_ID(next + n)` per table per sector) so
      children know their parents' ids without relying on
      auto-increment order; `refresh_containment` writes its changes in
      one statement; and PERF.6's benchmark shows the save faster with
      the same rows written.

    - [ ] **PERF.14 Reserve a sector's names in bulk, safe with several writers at once**
      Four `generate.py sector` runs started together (2026-10-01) hit 3
      deadlocks on `system_name_registry`: `reserve_system_name`'s
      `SELECT ... FOR UPDATE` of a name not yet there takes a gap lock,
      and two writers then insert into each other's gap. 3 of the 4 runs
      failed (with error 1452 after the rollback), and nothing retries
      error 1213. Today only `web/jobs.py`'s one-job lock prevents this.
      Done: a sector reserves all its names in one
      `SELECT ... WHERE base_name IN (...)`, resolves collisions in
      Python and writes one multi-row upsert, locking names in sorted
      order; the sector's transaction is retried on 1213 and 1205; and
      the 4-process test runs clean. PERF.8's parallel writers depend on
      this.

    - [ ] **PERF.15 Fewer queries per web page**
      Done, in `queryDb.py`: search facets (12 queries per request,
      including `COUNT(*)` of every star, planet and moon) are cached
      until `galaxy_content_state` changes; `system_detail` loads moons
      once per system, not once per planet; `sector_detail` loads its
      stars in one query, not one per system; `list_systems` and
      `_search_result_systems` join `stars` once instead of three
      correlated subqueries per row, with an index on
      `stars (star_system_id, role)`; and `galaxy_placed_sectors` reads a
      stored system count (with PERF.11) instead of counting per sector.

    - [ ] **PERF.16 Search names without scanning every row**
      Search builds `LIKE '%term%'` (`_search_like_pattern`) and counts
      every match exactly (`_search_page`), a full scan of each body
      table that won't survive a large galaxy. Done: names are searched
      with a FULLTEXT index on whole words, and counts stop at the
      300-row result cap ("300+"). Boss (2026-10-01): "full text index,
      to match whole words", so a search no longer finds text in the
      middle of a word ("ara" doesn't find "Kemaral").

    - [ ] **PERF.17 A time limit on web database statements**
      Web connections have no statement timeout, so one runaway query
      holds one of mod_wsgi's 5 threads until Apache's 60 s request
      timeout. Done: the read-only pool's init command sets
      `max_statement_time` (MariaDB) or `MAX_EXECUTION_TIME` (MySQL),
      default 10 s, configurable in `config.json`; a timed-out query
      returns a clear error page.

  - [ ] **PERF.7 Parallelize sector and system generation, with stable progress bars**
    Boss (2026-10-01): "add a TODO item to parallelize
    sector and system generation and update the progress bars so that
    they stay stable, I want to keep the ETA until done and elapsed time
    and I know that'll require some customization of the status bar code
    as time estimates are to be calculated from a decaying average based
    on number of runs per second." This is the first user of PERF.8's
    work queue. Today `generate.py` fills sectors one after another in a
    single process, and `_generation_progress()` (rich `Progress`) shows
    elapsed time and rich's own ETA. Done: sector fill and the plan's
    bright-star scatter run through PERF.8's queue; the bars stay
    pinned at the bottom without flicker (as PERF.4 requires) even with
    many workers reporting at once; every bar keeps elapsed time and an
    ETA until done; and the ETA comes from a custom column that uses a
    decaying (exponentially weighted) average of tasks finished per
    second rather than rich's built-in estimate. PERF.3's measured stars
    per second and PERF.4's slow-layer bar use the same rate. How
    (report of 2026-10-01): each system's child task only builds the
    system (no database); its sector's parent places each system as it
    comes back (placement uses the built system's size, so it can't
    happen first) and writes the whole sector in one batched
    transaction (PERF.13, PERF.14); the bright-star scatter's tasks are
    one layer (or one ring batch of a dense layer) each, not one per
    star, since 63 million tasks would cost more than the work. Needs
    PERF.13 and PERF.14 first. Open questions: the decay constant (how
    fast the average forgets older runs); whether the rate is counted in
    systems, stars or sectors (sectors differ a lot in size, so systems
    per second may be steadier; the default is systems); and whether
    `progress.json` for the web jobs reports the same decayed rate so
    the Generate page and UX.3's banner show the same ETA (default:
    yes).

  - [ ] **PERF.8 A parallel background work queue in the API**
    Boss
    (2026-10-01): "I also want to do it in a specific way, ideally how it
    would work is we'd have a work queue... In fact new TODO item,
    implement a parallel background work cue in the API. So in the plan
    phase we'll parallelize the bright star generation for any selected
    level. Each star generation will be a task there. Then a sector will
    be one task, but that task will be running until it finished filling
    the sector. Each star system generated will be it's own task that the
    API will paralleize, the sending task will get a signal that the
    parallel task is done when it's finished and will be able to use that
    information to calculate how long it iwll take to complete a given
    series of tasks that got handed out. We'll also need limits so the
    system never uses more than 80% of the total CPU power and defers to
    other running processes in process scheduling." Boss (2026-10-01) on
    the scheduler: "we should see if we want the task scheduler to be a
    daemon or if it'll only spin up a loop when there are tasks, since
    99.99% of the time there won't be." Today background work
    is `web/jobs.py`'s one-job-at-a-time runner (an `active` lock,
    `state.json`, `progress.json`) launching `generate.py`, which does
    everything serially. Done:
    - A work queue with a pool of workers, run by an on-demand
      supervisor, not a daemon: queuing work starts the supervisor when
      none is alive, it runs while there are tasks, and it exits after
      about 60 s idle, so an idle server runs nothing.
    - Plan phase: the bright-star scatter for the selected level
      (PERF.5's band) runs as parallel tasks.
    - Sector fill: each sector is one parent task that stays running
      until its sector is full; each star system in it is its own child
      task that the queue runs in parallel.
    - When a child task finishes, its parent gets a signal, and the
      parent uses those signals to work out how long its handed-out
      tasks will take to finish (feeding PERF.7's ETA).
    - Limits: the workers never use more than 80% of total CPU, and they
      run at a lower scheduling priority (for example `nice` on Linux and
      macOS, below-normal priority on Windows) so other processes on the
      machine come first.

    Decided defaults (report of 2026-10-01: https://claude.ai/code/artifact/4111c1a8-0d63-4b6c-b109-1b8389538a13):
    - Processes, not threads (the GIL), through a `ProcessPoolExecutor`
      with the spawn start method so Linux, macOS and Windows behave the
      same.
    - 80%: `max(1, floor(0.8 x cores))` workers, one fewer when MySQL
      runs on the same machine, each at `os.nice(10)` or
      `BELOW_NORMAL_PRIORITY_CLASS`; optionally hold back new tasks
      while the load average is above 80% of the cores.
    - Not in the API process: mod_wsgi runs one process with 5 threads
      that serve every page and is recycled by Apache. The supervisor is
      its own detached process (`src/workQueue.py`), spawned the way
      `web/jobs.py` spawns `jobRunner.py` today; no systemd unit,
      launchd plist or Windows service to install. Command-line
      `generate.py` becomes the supervisor itself, in the foreground with
      its bars.
    - Queue storage: a `tasks` table in the control database (job,
      parent, kind, payload, state, attempts, timings, result) plus a
      supervisor lease row refreshed every 5 s; a lease older than 30 s
      counts as dead and the next enqueue or Generate page load starts a
      new supervisor. Claims are `UPDATE ... LIMIT n` (no `SKIP LOCKED`,
      which MariaDB 10.4 lacks). Only sector-level and scatter tasks go
      in the table, added a ring batch at a time; per-system child tasks
      live in the supervisor's memory.
    - Seeds: each child task seeds `random` from the job seed, its sector
      address and its request number, so results don't depend on worker
      count or finish order (positions stay OS-random, as today).
    - Names and positions: children don't touch the database; the
      parent places systems and reserves names in bulk (PERF.14).
    - Restart or cancel: a sector is one transaction, so a half-done
      sector rolls back and its task returns to the queue; cancelling a
      job stops new dispatches and ends its workers.
    - Reset, the skeleton build and schema work keep `web/jobs.py`'s
      one-at-a-time lock.

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

Design: [docs/design/login-brute-force-protection.md](design/login-brute-force-protection.md)

The login protection research of 2026-10-01 (what's built today, the
gaps found, standard methods compared, and the recommended design) is in
the design doc above; its section numbers are cited below. Order: SEC.20,
SEC.1 with SEC.21, SEC.23, SEC.22, SEC.24, SEC.25, SEC.26, SEC.27.

- [ ] **SEC.20 Log every failed and locked login with its address**

  Design: [docs/design/login-brute-force-protection.md](design/login-brute-force-protection.md), sections 1 and 3 (step 1)

  Today nothing records a failed or locked login: `admin_audit_log` only
  gets admin actions (`authz.py` calls `adminAuth.record_audit`), and a
  wrong password on the web `/login` form answers 200 with an error
  message (`web/admin_pages.py`, `login`), so Apache's access log can't
  tell it from a page view. Done: every failed login, every lock (per
  username today, per address with SEC.1) and every wrong current
  password (SEC.23) writes one log line with the time, client address
  (`request.remote_addr`) and username, and a row in `admin_audit_log`
  (`login.failed`, `login.locked`; the username as typed, capped in
  length; never the password); the `/login` form answers 401 for a wrong
  password and 429 for a lock; the admin stats page shows recent
  failures. The log line's format is fixed and documented, so SEC.27's
  fail2ban filter can match it. Open questions: the debug log (only
  written when `debug` is on) or a separate always-on auth log file, and
  where it lives (next to the debug log, mode 0640)? How long audit rows
  for failures are kept (a flood of failures shouldn't grow the table
  without bound: prune after 90 days, or cap the count)?

- [ ] **SEC.1 Lock out an IP address after failed logins**

  Design: [docs/design/login-brute-force-protection.md](design/login-brute-force-protection.md), section 3 (step 2)

  Boss: "3 failed
  login attempts triggers the script refusing to allow that IP address
  to login again for an increasing amount of time (Starts at 5 minutes,
  doubles every time for a max of 1 day). This should get integrated
  into the control database." This is separate from the per-username
  backoff already shipped (`src/html/api/loginbackoff.py`: 10 free
  failures per username, then 1 s doubling to 15 min, kept in memory
  per worker). Done: a table in the control database
  (`control_schema.sql`, alongside `admin_users`/`admin_sessions`)
  recording failures and lockouts per client IP (`request.remote_addr`,
  which already honors the `proxy_fix` setting); after 3 failures that
  IP gets 429 + `Retry-After` for 5 minutes, doubling on each further
  lockout up to 1 day; checked before the password, shared by all
  workers, kept across restarts, and logged as SEC.20 describes. The
  research adds these safeguards (default taken; Boss can change them):
  - `127.0.0.1`, `::1` and addresses in a new config allowlist are never
    locked, so the server's own admin can always get in.
  - When the site sits behind a proxy whose address isn't unwrapped
    (`proxy_fix` unset), every visitor shares one address and a lockout
    would block everyone for up to a day. If most logins come from one
    private or loopback address, the app warns at start-up and on the
    admin stats page.
  - A successful login clears that address's failure count; its
    doubling level decays (halves after a day with no lockout) rather
    than resetting, so one right guess doesn't wipe it.
  - An admin page and a command-line tool list current lockouts and
    lift one.
  - Covers every place a password is checked: `/login`, the API login,
    change-credentials (SEC.23), and later user logins, password resets
    and invite links (USR.4, USR.5) and the second step of SEC.26.

  Open questions: do both limits stay (per IP and per username)?
  Default: yes, see SEC.21. IPv6: lock the single address or its /64
  (one machine often holds a whole /64)? Default: the /64.

  - [ ] **SEC.21 Keep the per-username backoff in the control database too**
    Today `loginbackoff.py` keeps its counts in memory, so a restart
    clears every lock and each worker process counts separately (the
    deployment guides run one process today, so this only bites with
    more). Done: the per-username counts and locks live in SEC.1's
    table (or a sibling), shared by every worker and kept across
    restarts, with the same numbers (10 free failures, 1 s doubling to
    15 min) and the same rule that unknown usernames are counted like
    real ones; `LOGIN_BACKOFF_ENABLED` still turns it off for tests.

  - [ ] **SEC.22 A trusted-device cookie so lockouts can't shut out the real admin**

    Design: [docs/design/login-brute-force-protection.md](design/login-brute-force-protection.md), section 3 (step 4)

    Anyone who knows an admin's username can keep it locked by failing
    on purpose (15 minutes at a time today). OWASP's answer is a device
    cookie. Done: a successful login sets a long-lived, signed device
    cookie for that account (only its hash stored in the control
    database); a login from a browser holding a valid device cookie for
    that username skips the per-username lock (SEC.21), while the
    per-address lockout (SEC.1) and the per-IP rate limit still apply;
    changing credentials, or an admin action, revokes the account's
    device cookies. Open questions: how long a device cookie lasts
    (default 90 days); whether a device cookie also skips SEC.1's
    per-address lock (default: no).

- [ ] **SEC.23 (bug) Wrong current passwords on `/account` aren't counted**
  `POST /api/auth/change-credentials` re-checks the current password,
  but a wrong one isn't counted by the login backoff, and the web
  `/account` form only falls under the shared page limit (300 a minute;
  in-process API calls skip the API's default limits, `api/limiter.py`).
  Someone holding a stolen session cookie can guess the password there
  quickly. Done: a wrong current password counts as a failed login for
  that username (SEC.21) and address (SEC.1), is refused while either is
  locked, is logged (SEC.20), and the route gets the same explicit
  per-IP limit as login (`auth.LOGIN_RATE_LIMIT`). Can ship before SEC.1
  using the in-memory backoff.

- [ ] **SEC.24 Refuse common and breached passwords**

  Design: [docs/design/login-brute-force-protection.md](design/login-brute-force-protection.md), sections 2 and 3 (step 6)

  Today a password only has to be 12 characters and differ from the
  username (`adminAuth.MIN_PASSWORD_LENGTH`). NIST SP 800-63B-4 asks for
  a check against breached, common and expected passwords. Done: setting
  a password (change-credentials, the installer's first password, and
  later USR.5's reset) refuses one on a bundled offline list of common
  and breached passwords, plus a few site words ("planetgen", the
  username), with a message saying why; no composition rules are added.
  Open questions: minimum length 12 or NIST's 15 for a password used
  alone? Default: 12 until SEC.26 exists. Offline list only, or also Have
  I Been Pwned's k-anonymity API? Default: offline only. Which list and
  its license (for example the top 100,000 of a public corpus)?

- [ ] **SEC.25 Check the password hashing cost and re-hash on login**

  Design: [docs/design/login-brute-force-protection.md](design/login-brute-force-protection.md), section 3 (step 7)

  Passwords use werkzeug's default (scrypt N = 2^15, r = 8, p = 1 in
  werkzeug 3.1.9); OWASP suggests scrypt N = 2^17, argon2id, or
  PBKDF2-SHA256 at 600,000 rounds. Done: time a hash on a modest server;
  if well under about 250 ms, pick stronger settings in one place in
  `adminAuth.py` (minding memory: scrypt N = 2^17 needs 128 MiB per
  login, times the five worker threads); a successful login re-hashes a
  stored hash made with older settings. Open question: is a new
  dependency (argon2-cffi) acceptable, or stay with werkzeug's built-ins?
  Default: werkzeug's built-ins.

- [ ] **SEC.26 Two-factor sign-in (TOTP) for admins**

  Design: [docs/design/login-brute-force-protection.md](design/login-brute-force-protection.md), section 3 (step 8)

  The strongest single defense against a guessed or leaked password.
  Done: an admin can turn on a time-based one-time code (any
  authenticator app; QR code on the account page), gets single-use
  recovery codes (stored hashed), and then signs in with password plus
  code; the code step is throttled and logged like the password step;
  API keys are unaffected. Open questions: optional or required for
  admins (default: optional now, required for the Owner once USR.2
  exists)? Can an admin reset another admin's second factor, or only the
  command line? Is "remember this device for 30 days" (tied to SEC.22's
  cookie) allowed?

- [ ] **SEC.27 A fail2ban recipe in the deployment docs**

  Design: [docs/design/login-brute-force-protection.md](design/login-brute-force-protection.md), section 3 (step 9)

  Done: `docs/deployment/` gets an optional Linux section with a
  fail2ban filter matching SEC.20's log line and a jail that bans an
  address at the firewall after repeated failures (ignoring the server's
  own addresses and SEC.1's allowlist), and how to check and lift a ban;
  a short note on Apache-level options (mod_evasive, mod_security) and
  why they're optional. Nothing is installed by `install.sh`. Waits on
  SEC.20.

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
