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

NAV, DB, SEC, DOC and POP have no open items today.

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

1. **First:** PERF.5 (scatter bright stars in stages), done in PR #229.
2. **Galaxy Map and units:** done (UX.13, UX.14, MAP.15, MAP.30 and
   MAP.2 in PR #234, GEN.23 in PR #226).
3. **Generation estimates and progress:** done (PERF.3 and PERF.10 in
   PR #238, PERF.4 and PERF.9 in PR #258).
4. **Admin editing:** done (ADM.5 in PR #235, ADM.8 in PR #244, ADM.6
   and ADM.7 in PR #260).

Waiting behind those: PERF.11, UX.2, UX.3, GEN.9,
user accounts (USR.1,
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

- [ ] **MAP.57 The System Map writes NaN or infinite positions into its SVG (bug)**
  Found by the generation tests (2026-10-01): a body whose computed
  position is NaN or infinite is written straight into the System
  Map's SVG. Done: such a body is left out or drawn at a safe place
  with a note, the SVG never holds NaN or inf, and the strict xfail
  test for it passes.

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

- [ ] **MAP.59 Make it plain that a zoomed-in slab is a slab, not a wedge**
  Boss (2026-10-01 20:45Z): "We need to make it clearer, when we've
  zoomed into a specific slab, that we're viewing a specific slab and
  not a wedge. I'm not sure how to do that so do some research on that
  and then add the to-do items to make it happen."

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

- [ ] **MAP.78 Zooming into a wedge must show the whole wedge at every drill-down level (bug)**
  Boss (2026-10-01 21:11Z): "when the system zooms into a wedge, make
  sure it is the entire wedge as you drill down." Done: whenever the
  view zooms to a wedge (MAP.52) or to a slab or segment inside it
  (MAP.56), the zoom frames all of what was picked, with no part cropped
  by the map's edges or by the controls over it; the fit uses the map's
  actual size (MAP.53) and holds while the view rotates and when the
  window is resized. Under MAP.58's locked zoom, the locked level is
  this whole-wedge fit, not a closer one.

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
  focus), the help text says which, and a test pins it. Open question:
  which keys does Boss want?

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
    Done: from the corridor query of NAV.10, every
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
  Needs a galaxy schema migration for the position indexes.

  - [ ] **NAV.11 Travel times for the system-to-system route too**
    Today warp and fold times are shown only for the direct distance;
    the route shows only its length. Done: the route gets the same warp
    and fold tables, per hop and in total.

  - [ ] **NAV.12 A maximum hop length (open question)**
    The route's hops are the 6 nearest neighbours, with no limit on hop
    length, so a hop can be very long in a sparse region. Open question
    for Boss: should a route have a maximum hop (a ship's range)?
    Default: no limit, but the longest hop is shown.

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
  its goal of the same sectors at any worker count. Done: Boss decides
  whether generation should be reproducible from its seed; if yes,
  every draw comes from the seeded generator (per sector, so worker
  count doesn't matter) and a test checks that one seed gives the same
  sectors at any worker count. Open question for Boss: should a seed
  reproduce a galaxy?

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
    expected-versus-actual ratio to correct their estimates. Open
    questions: columns on `sectors` (a migration in the Database
    workstream) or a separate table; whether "actual" counts systems,
    stars, or both; what the decaying average is taken over (the ratio
    per density bucket, so it ties in with PERF.10, or one galaxy-wide
    figure); whether existing sectors are backfilled by a migration;
    and what happens to the stats when a sector is regenerated
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

## DB: Database and schema

DB.1 shipped in 7.35.0 (PR #152).

- [ ] **DB.2 Asteroid field and comet composition rows are written but never read (bug)**
  Found by the Database tests thread (TEST.11, PR #315, 2026-10-01):
  `asteroid_field_composition` and `interstellar_comet_composition` rows
  are saved, but nothing reads them back, so pages show the parent's
  `composition_summary` instead. Done: either the pages and the API read
  these rows, or generation stops writing them and a migration drops the
  tables. Open question: which of the two?

- [ ] **DB.3 resetDb while another process holds id blocks can duplicate primary keys (bug)**
  Found by the Database tests thread (PR #315, 2026-10-01): running
  `resetDb` from one process while another long-lived process (the web
  app or a generation worker) still holds cached id blocks lets the old
  process hand out ids the reset database gives out again, so inserts
  fail on duplicate primary keys. Done: a reset can't lead to reused
  ids, probably by not restarting ids at 1 after a reset (or by making
  holders drop their blocks), with a test that runs both processes.

- [ ] **DB.4 A database with an emptied schema_migrations table is treated as current (bug)**
  Found by the Database tests thread (TEST.9, PR #315, 2026-10-01): an
  old database whose `schema_migrations` table has been emptied is
  treated as up to date, so its migrations never run. Done: when the
  table is empty or missing on a database that has tables, the version
  is detected from the table shape (which tables and columns exist), the
  needed migrations run, and a test covers it.

- [ ] **DB.5 Several first connections to an empty database race to create the schema (bug)**
  Found by the parallel, population and navigation tests thread
  (2026-10-01), which worked around it in its tests: when several
  connections reach an empty database at once, each runs
  `_ensure_schema`, and one fails with IntegrityError 1062 "Duplicate
  entry '49' for key 'PRIMARY'" on `schema_migrations`; in a test run,
  that connection's teardown then left `DROP DATABASE` hanging. Real
  runs aren't hit today because the main process creates the schema
  before workers start. Done: schema creation is safe when several
  connections start at once (one creates it under a lock, the others
  wait and then see it current), with a test that opens several first
  connections together.

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

- [ ] **API.6 Admin-created user-level API keys that can read but not upload**
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
  out. Open question: does a user-level key belong to a user account
  (USR.1) or stand alone until user accounts exist?

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
