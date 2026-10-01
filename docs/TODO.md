# planetGen Roadmap

## How to use this document

This file tracks **open work only**. Finished work is not listed here:
`CHANGELOG.md` and git history are the record. For how things work today,
read the code and the reference docs (`src/stellarObjects/schema.sql` and
`docs/database-schema.md` for the schema, `docs/api.md` and
`docs/html-interface.md` for the web interface,
`docs/design/galaxy-coordinate-system.md` for galaxy geometry). When an
item ships, delete it here, renumber the rest, and describe the change
in the PR's `changes/` note.

Each open item says what's wrong (the observed symptom), where to start
looking, and what "done" means, so it can be picked up cold.

## Background

The long-term goal is a fully populated galaxy (every sector, system,
planet, moon and belt) stored in MySQL and browsed through a web
interface served from Apache. The original roadmap phases are all done:
object-graph serialization, the relational schema and migrations, CLI
tools writing to the database, lazy galaxy-scale generation from a
density skeleton, and the Flask API (`src/html/api/`) with server-rendered
pages served by the same app (`src/html/web/`): Galaxy, Sector and System
maps, search, NAV, admin auth and wiki publishing.

## Open items

Numbered in the order to work them, from the smallest change to the
largest, with new bug fixes put first. The numbers run across groups;
renumber when items are added or finished.

### Plan: what to do first

- **Extend the cache (8)**.
- **Galaxy navigation (70-79):** Boss's drill-down design of
   2026-10-01, specified in `docs/design/galaxy-drilldown-navigation.md`.
   70 (the nested block ladder), 71 (the stage contents API) and 72
   (the stages, which the map now opens on; the old free camera stays
   behind a Free look button) have shipped; 77 and 79 next.
- **Galaxy Map (12-19):** Boss approved the plan in the
   project's `galaxy-megablocks/report.md` (hybrid master-wedge
   slots, pixel-sized mega-blocks). 13-18 have shipped (pixel-sized blocks, the solid and its slice, filled and unfilled blocks with no marker dots, block info, smooth zooming, and the three.js decision in `html-interface.md`). 12 (the
   hybrid master-wedge slot rule) shipped in schema v35, and so have 19's
   follow-ups (opaque blocks drop the faces they share; translucent
   sorting showed no artifacts; block size stays per CSS pixel).
   Distance-based detail was left to the drill-down (see the note under
   "Galaxy navigation").
- **Features (25-36)** from the same notes: phenomena views (25), nebulae and
   remnants: placement, classes, containment and naming, plus asteroid
   field classes (27-31), the correlative update (32), navigation frames
   and speeds (33-34), and facilities (35-36). 28 and 31 (classes)
   shipped in schema v38, 29 (containment) in v39, 30 (naming) in v40,
   26 (nearest systems) in v41 and 35 (facilities) in v42; 32 (the
   correlative update moves everything) and 36 (placing facilities
   from the web) shipped with them, and so did 33-34 (navigation
   frames and speeds).
- Each change site in the code carries a `TODO(<area> #N)` comment
   naming its item here (areas: distances, system-list, site-header,
   search, phenomena, galaxy-map, sector-map, orbits, nav, facilities,
   security, physics, web-pages);
   grep for `TODO(` to see them all, or `TODO(galaxy-map` for one area.
- Items with a **Question for Boss** state the default taken; the work
   can start on that default.

### Galaxy Map and the sector standard (`static/galaxyprisms.js`, `static/galaxymap3d.js`, `lib/galaxymap3d.py`, `stellarObjects/galaxyGeometry.py`)

Today the map draws the analytic density as shrunk prisms, m sectors a
side (m a power of 3), sized by a volume budget that badly overestimates
the thin disk. The result is 70-290 px cubes with gaps, and the spiral
barely shows. The plan (report above, with renders) replaces that with a
continuous solid of mega-blocks sized from the screen's pixel scale.

Bugs Boss found on the map (2026-10-01): "bug fix, wedge lines should
not extend past the boundary of the galaxy. bug fix, bright stars do not
display past 500 seconds scale (1 block = 81 sectors across). bug fix
zooming reveals stars are being drawn but it takes a while to load,
another bugfix, it's really hard to find the generated star system on
the map, so everything not yet filled should be more transparent by a
lot with a much higher contrast."

96. [ ] **Bug: wedge lines run past the galaxy's edge.** Boss: "wedge
    lines should not extend past the boundary of the galaxy." Today
    `galaxymap3d.js` (`buildWedgeLines`) draws every master line of
    `galaxyprisms.wedgeLines` as a straight radial line out to
    `GALAXY_RADIUS * 1.02`, a fixed circle, so the lines carry on past
    the galaxy's real outline (which is not a circle at every bearing)
    and slightly past the radius itself. Done: each wedge line stops at
    the galaxy's boundary along its bearing, in the 3D view and the
    drill-down stages (items 70-72) alike. Open questions: which
    boundary counts (the outermost generated ring at that bearing, the
    density model's cutoff from `densityShape`, or the outermost layer
    extent from `galaxySkeleton.build_layer_extents`); and whether the
    bearing labels move in to the new line ends.

97. [ ] **Bug: bright stars vanish when zoomed out.** Boss: "bright
    stars do not display past 500 seconds scale (1 block = 81 sectors
    across)." Past that zoom the pre-placed bright stars (v43,
    `bright_stars`, drawn as points with a glow) stop showing; they
    should show at every zoom. Leads to check, not yet confirmed: each
    tile carries at most `queryDb.GALAXY_TILE_MAX_BRIGHT_STARS` (400)
    stars, and zoomed out the view may switch to tiles or a view
    radius that leaves stars out. Done: bright stars draw at every zoom
    out to the whole galaxy, thinned by luminosity if there are too many
    rather than disappearing. Open questions: "500 seconds" is taken as
    the scale readout at the zoom where blocks are 81 sectors a side;
    Boss to confirm which readout he meant. How many stars should the
    whole-galaxy view draw (the brightest N overall, or the brightest
    per tile)?

98. [ ] **Bug: stars take a while to appear after a zoom.** Boss:
    "zooming reveals stars are being drawn but it takes a while to
    load." After a zoom the bright stars (and the tiles they come in,
    `renderFromCache` in `galaxymap3d.js`) arrive late, so the view
    shows them popping in. Done: stars already loaded stay on screen
    through a zoom, the tiles for the new view load faster or ahead of
    time, and nothing visibly pops in. Open questions: where the time
    goes (the tile request, `galaxy_bright_stars_in_box`'s query, or
    rebuilding the points); whether to prefetch the next zoom level's
    tiles; and whether a separate, lighter star endpoint would help
    (ties in with item 90's fewer, bigger database calls).

99. [ ] **Bug: generated systems are hard to find on the map.** Boss:
    "it's really hard to find the generated star system on the map, so
    everything not yet filled should be more transparent by a lot with
    a much higher contrast." Done: blocks and sectors not yet filled
    draw much more transparent, and filled sectors stand out with much
    higher contrast against them, on the 3D map and on the drill-down
    stages (items 70-72). Open questions: how transparent the unfilled
    blocks go (and whether the density shape still reads at the galaxy
    scale); what "higher contrast" uses (a bright color, an outline, a
    glow like the bright stars); whether a block holding only a few
    filled sectors gets the filled look; and whether it follows the
    light and dark themes and keeps enough contrast in both.

### Galaxy navigation: the drill-down (`docs/design/galaxy-drilldown-navigation.md`)

Boss's design of 2026-10-01: the Galaxy Map becomes a drill-down. In 3D,
pick a slab; it is pulled out and shown from above; pick a block; its
contents fill the view as blocks 1/9 the size; repeat until single
sectors, where a click opens the sector. The ladder is 243 -> 27 -> 3 ->
1 sectors a side ("the bigger targets"), eight clicks from the galaxy to
a sector. Admins can generate a sector, a layer or a neighborhood (radius
asked in light-years) at the sector level, and the NAV page can pick its
start and destination on the map or in a sector. Everything below is
specified, with the math, in the design doc named in the heading; each
item names its section. 70 (the nested ladder,
`stellarObjects/galaxyDrill.py`), 71 (`GET /api/galaxy/stage`) and 72
(the stages, `static/galaxystages.js` and `static/galaxystageview.js`,
with stage URLs `/galaxy?slab=`, `?at=`, `?sector=<designation>`) have
shipped, and so have 77 (the address bar, `/galaxy/locate`) and 79 (the
course overlay, `/galaxy?course=<from>,<to>`); Web can do 63, 74, 75 and
78 now. The old free camera stays
behind the map's Free look button until Boss settles decision 2 (section
11). If it stays, its blocks could also grow with distance from the
camera (bigger blocks on the far side of the view ball, where the nested
ladder keeps the borders seamless); each stage draws one level, so the
stages themselves don't need it.

75. [ ] **NAV page picks on the map (section 9).** Boss: "from the nav
    menu select start and destination using either the text dropdowns as
    we have now or the galactic map interface to select. If it's within
    sector then it'll just use the sector interface." Done: beside each
    dropdown, "Pick on Galaxy Map" (`/galaxy?pick=...`, generated-only
    forced on, ending in 74's Sector Map pick mode), "Pick in this
    sector" once the other end is known, and a Bookmarks select. The
    Galaxy Map's side is in: `?pick=` shows the banner with Cancel back
    to NAV, keeps "Generated only" on, and a sector click opens that
    sector in pick mode.

76. [ ] **Bookmarks (section 8.2).** Done: a ☆ on the breadcrumb and
    info panels saves a stage, sector, system or phenomenon in
    `static/bookmarks.js` (per browser, up to 100, storage failures
    tolerated), with a map menu, Ctrl+1-9, rename and delete, and the
    entries offered by the NAV pickers. Shared bookmarks need Boss's
    decision 4 and a migration.

78. [ ] **"Show on Galaxy Map" links (section 8.1).** The sector page's
    link goes to the Quadrant table today. Done: sector, system and
    search pages link to `/galaxy?sector=<designation>`.

100. [ ] **Bug: no free camera; drill down from a top-down view by
    wedge, slice and block.** Boss (2026-10-01): "bugfix, remove the
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
    stays behind the map's Free look button (added with item 72, PR
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
      targets, under the Web interface section).
    - The design doc's sections 5 and 10-11 are updated to match, and
      breadcrumb, Back, stage URLs (`/galaxy?slab=`, `?at=`), bookmarks
      (item 76), the address bar (item 77) and the NAV course (item 79)
      keep working with the new steps.
    Ties in with the map bugs 96-99 (wedge lines stopping at the
    galaxy's edge matter more once wedges are what the user picks, and
    item 99's contrast applies to the slices and blocks). Open
    questions: how the "wedge (quarter)" maps onto the design's nested
    ladder (243 -> 27 -> 3 -> 1, section 3): is the first pick always
    one of four quarters, or one of the master wedges at that ring?
    What a "slice" is: a ring band, a layer (the current 3D stages pick
    a slab, a layer), or a row of blocks along the wedge? With no 3D
    view, how does the user pick a layer above or below the galactic
    plane (a side view, a layer list, or slices that run through the
    disk's thickness)? What replaces the free camera for item 73's
    neighborhood generate and item 75's NAV picking, which may want an
    arbitrary point?

101. [ ] **Bug: "Show on Galaxy Map" should open at the sector, and the
    map needs its own Back and Forward.** Boss (2026-10-01): "when the
    user clicks "Show on Galaxy" It should be zoomed in to the sector
    level of the slice that we can see that sector, and back and forward
    buttons to travel ones own history on the map display". Today
    `/galaxy?sector=<designation>` already opens the sector's stage 8
    (`galaxystages.parseStageQuery`), but not every "Show on Galaxy
    Map" link uses it: the NAV result's link opens the course overlay
    (`/galaxy?course=<from>,<to>`, item 79) over the galaxy, and the
    sector and system pages still link to the Quadrant table (item 78).
    The map has no Back or Forward of its own; moving between stages
    relies on the breadcrumb and the browser's Back button. Done:
    - Every "Show on Galaxy Map" link (sector, system and search pages,
      the NAV result) opens the map zoomed in to the sector level of
      the slice holding that sector, with the sector in view and
      highlighted.
    - The map display has Back and Forward buttons that step through
      the user's own history on the map (each stage or position they
      visited), separate from but consistent with the browser's
      history and the stage URLs.
    This follows the wedge, slice and block drill-down of item 100
    (the "slice" here is item 100's slice), and folds in item 78's links.
    Open questions: for the NAV result, which sector is shown (the
    start, the destination, or a view that fits both, as item 79's
    course does today)? Is the map's history its own list or the
    browser's history (`history.pushState` per stage) with the buttons
    calling `history.back()`/`forward()`? Does it survive a page reload
    or a visit to a sector page and back? How far back does it go?

### Web interface (`src/html/web/`, `src/html/static/`)

Web workstream. Boss's UX reference for both items below is
"Responsive Web Design Standards" (Boss's notes of 2026-10-01; a copy is
in the project's shared files at `ux-standards/`). The rules from it
that apply here: size classes by width (compact under 600 px, medium
600-839, expanded 840-1199, desktop 1200+), not by device; viewport media
queries only for the top-level frame and container queries for
everything inside it; on medium and wider screens a menu must not stretch
across the screen like a phone's; text columns capped at 45-75
characters, but a map or canvas may use the full width; touch targets at
least 44-48 px on coarse pointers (`pointer: coarse`), smaller is fine
for a mouse; spacing and type sized with `clamp()`.

62. [ ] **Menus sized to what they hold.** Boss (2026-10-01): "I want the
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

87. [ ] **Warn every visitor while a background job changes the
    galaxy.** Boss (2026-10-01): "a warning to all users on the UI when a
    task is running in the background which is modifying the starmap is
    going on with an ETA until it will be finished to the nearest hour
    rounding up." Today only the admin Generate page shows a running job
    (`web/jobs.py`: the `active` lock in the jobs directory, each job's
    `state.json`, and `progress.json` written by `generate.py` through
    `stellarObjects.progressFile`). Done: while a job that writes to the
    galaxy is running (plan, bright-star scatter, sector fill, reset, new
    galaxy, and later item 60's regenerate and item 73's slice and
    neighborhood generation), every page shows a banner to every
    visitor, signed in or not, saying the galaxy is being changed and
    when it should finish, as an ETA rounded up to the next whole hour
    (from `progress.json`, and from item 86's measured stars-per-second
    rate). Open questions: a banner at the top of every page, or only on
    the map and list pages? Does it also cover `generate.py` runs started
    from the command line, which don't take the web jobs lock today
    (they would need to write the same lock and progress file)? What does
    it say when there's no ETA yet?

102. [ ] **One meaningful-unit ladder for speeds.** Boss (2026-10-01):
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

103. [ ] **One meaningful-unit ladder for time periods.** Boss
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
    does); and whether elapsed-time and ETA displays for jobs (items 86-88)
    and the admin pages use the same function.

105. [ ] **Bug: put an object's data beside its 3D render when there's
    room.** Boss (2026-10-01): "if there is enough room next to the 3D
    render of an object, put the data segment next to the object." The
    3D renders (`static/bodyRendering.js`, used by the System Map,
    `static/systemmap.js`, and the Sector Map, `static/sectormap.js`)
    show an object's details in an info panel, which today can sit below
    the render even when the screen has space beside it. Done: when the
    space next to the render is wide enough, the data panel sits beside
    the object; when it isn't (phones, narrow windows), it stays below;
    the switch follows the Responsive Web Design Standards' size classes
    and container queries (see the notes above item 62), with no layout
    jump while the render loads. Open questions: which panels this
    covers (the System Map and Sector Map info panels, the object pages
    for planets, moons, stars and phenomena, or all of them)? What
    "enough room" means (a minimum width for the render plus a readable
    45-75 character text column)? Which side the panel goes on?

### Sector Map and generation (`static/sectormap.js`, `web/generate_page.py`, `generate.py`)

86. [ ] **Estimate size and time before bulk generation, and refuse
    what won't fit.** Boss (2026-10-01): "any directive to generate
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
    sectors, the admin Generate page's jobs, item 73's slice and
    neighborhood generation, item 60's sector regenerate, and the Galaxy
    Map's Generate buttons):
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

88. [ ] **A second progress bar for slow layers in the plan.** Boss
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
    to sector fill (`Sectors (ring … layer …)`), where item 86's rate is
    measured.

89. [ ] **Scatter bright stars in stages, one luminosity band at a
    time.** Boss (2026-10-01): "let's do a default of 100 solar
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
    down to 100 later is the kind of extra layer this item adds. Today the plan's scatter
    (`generate.py`, `--bright-star-min-luminosity`) clears `bright_stars`
    and redraws everything at or above the threshold, and refuses when
    any sector is already filled unless `--force` leaves those sectors
    out. Done: the galaxy stores its current star-fill level (the lowest
    luminosity already scattered); a scatter to a lower threshold keeps
    the existing bright stars and adds only stars from the new threshold
    up to (not including) the stored level, then lowers the stored
    level; asking for a level at or above the stored one does nothing
    (or says so); the Generate page and the CLI show the current level
    and offer "go down to N"; and item 86's size and time estimates cover
    just the new band. Open questions: where the level is stored (the
    control database, a `galaxy` row next to the skeleton, or derived
    from `MIN(luminosity_w)` in `bright_stars`)? What happens to sectors
    already filled when new, dimmer bright stars land in them: add the
    stars and build their systems in place, skip those sectors, or mark
    them for item 60's regenerate? Must the new band draw from the same
    random stream, so a 500-then-100 galaxy matches a straight-to-100
    one?

90. [ ] **Rate-limit SQL calls and make each call do more.** Boss
    (2026-10-01): "ratelimiting calls to the sql database and seeing if
    we can investigate some way to make our DB calls more efficient, do
    more with less calls without impacting performance." Today every
    database call goes through `stellarObjects/_db.py`'s `PooledDB`
    (one pool per config, `maxconnections=10`, `blocking=True`), with no
    limit on how often calls are made, and the Database thread measured
    about 60% of generation time going to per-row planet and moon saves.
    Done: an investigation first, written up before any code changes,
    that lists the hottest call sites (per-row inserts during sector
    fill, per-row reads on the web pages) and for each one whether it can
    be batched (`executemany`, multi-row `INSERT`, one save per system
    instead of per body, fewer round trips per page); then the batching
    changes that measure faster, and a rate limit on calls to the
    database. Open questions: what the rate limit protects against (the
    MySQL server being swamped by item 92's parallel workers, or web
    users hammering the API), and so whether it is a calls-per-second cap,
    a cap on concurrent connections, or both; whether it is one limit
    shared by generation and the web site or separate ones; how it
    interacts with the pool size once item 92 runs several workers at
    once (each process gets its own pool today); and what benchmark
    decides "without impacting performance" (a fixed test sector timed
    before and after?).

91. [ ] **Parallelize sector and system generation, with stable
    progress bars.** Boss (2026-10-01): "add a TODO item to parallelize
    sector and system generation and update the progress bars so that
    they stay stable, I want to keep the ETA until done and elapsed time
    and I know that'll require some customization of the status bar code
    as time estimates are to be calculated from a decaying average based
    on number of runs per second." This is the first user of item 92's
    work queue. Today `generate.py` fills sectors one after another in a
    single process, and `_generation_progress()` (rich `Progress`) shows
    elapsed time and rich's own ETA. Done: sector fill and the plan's
    bright-star scatter run through item 92's queue; the bars stay
    pinned at the bottom without flicker (as item 88 requires) even with
    many workers reporting at once; every bar keeps elapsed time and an
    ETA until done; and the ETA comes from a custom column that uses a
    decaying (exponentially weighted) average of tasks finished per
    second rather than rich's built-in estimate. Item 86's measured stars
    per second and item 88's slow-layer bar use the same rate. Open
    questions: the decay constant (how fast the average forgets older
    runs); whether the rate is counted in systems, stars or sectors
    (sectors differ a lot in size, so systems per second may be steadier);
    and whether `progress.json` for the web jobs reports the same decayed
    rate so the Generate page and item 87's banner show the same ETA.

92. [ ] **A parallel background work queue in the API.** Boss
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
    other running processes in process scheduling." Today background work
    is `web/jobs.py`'s one-job-at-a-time runner (an `active` lock,
    `state.json`, `progress.json`) launching `generate.py`, which does
    everything serially. Done:
    - A work queue the API owns, with a pool of workers.
    - Plan phase: each bright star generated for the selected level
      (item 89's band) is one task.
    - Sector fill: each sector is one parent task that stays running
      until its sector is full; each star system in it is its own child
      task that the queue runs in parallel.
    - When a child task finishes, its parent gets a signal, and the
      parent uses those signals to work out how long its handed-out
      tasks will take to finish (feeding item 91's ETA).
    - Limits: the workers never use more than 80% of total CPU, and they
      run at a lower scheduling priority (for example `nice` on Linux and
      macOS, below-normal priority on Windows) so other processes on the
      machine come first.
    Open questions: processes or threads (Python's GIL means CPU-bound
    generation needs processes, which then each need their own database
    connections, so the pool size and item 90's rate limit have to fit
    the worker count); how 80% is enforced (a worker count of 80% of the
    cores, or measuring load and throttling); whether the queue lives
    inside the API process or in a separate worker service the API talks
    to, and how command-line `generate.py` runs use it; how random seeds
    are handed to tasks so a galaxy comes out the same however many
    workers ran it and in whatever order they finished; how two systems
    in one sector avoid clashing on names and positions when they are
    built at the same time; and what happens to queued and half-done
    tasks when the server restarts or a job is cancelled.

93. [ ] **Weight the bright-star ETA by the shape of the galaxy.** Boss
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
    item 86 (the up-front time estimate uses the same per-layer
    weights), items 87 and 88 (the banner's ETA, and item 88's
    stars-remaining estimate for the current layer), item 89 (a staged
    scatter weights only the new luminosity band), and items 91 and 92
    (the decaying-average rate and the parallel tasks). Open questions:
    is the weight the expected star count alone, or does it also count
    rings and slots walked (an empty edge layer still costs some loop
    time)? How does the weighting combine with item 91's decaying
    average: the average measured in expected stars per second, or in
    layers per second and then scaled? Is the per-layer expected count
    worked out in a quick pre-pass at the start of every plan, or stored
    with the galaxy skeleton? Once item 92 runs layers in parallel and
    out of order, does the ETA add up the expected work still queued
    rather than following the layer order?

94. [ ] **Record generation speed across a log scale of densities.**
    Boss (2026-10-01): "record generation stats such as time per star,
    time per sector for a log scale of densities from 0.01 to the max
    expected density / actual density found. ... Both of these will be
    continued to be refined and calculated as long as the galaxy is in
    existence but as a decaying average." Today nothing records how long
    generation takes per star or per sector, and item 86's planned
    stars-per-second figure is a single number for the whole server.
    Done: density is split into log-scale buckets from 0.01 up to the
    highest density expected or found; every sector fill adds its time
    per star and time per sector to its density's bucket as a decaying
    average; the buckets keep updating for as long as the galaxy exists;
    and they can be read back by the tools below. These stats feed
    item 86 (time estimates before bulk generation, per bucket instead
    of one rate), item 87 (the banner's ETA), item 88 (the slow-layer
    bar's stars-remaining ETA), item 89 (the time for a new luminosity
    band), item 91 (the decaying-average rate) and item 93 (weighting
    the bright-star ETA by expected work per layer). Open questions: how
    many buckets and where their edges sit (per decade, half-decade?);
    whether the top edge is fixed from the density model's expected
    maximum or moves up when a denser sector is found; the decay
    constant (how fast old runs fade); whether the plan's bright-star
    scatter gets its own buckets (its cost per star differs from sector
    fill); whether it lives in the control database (survives a new
    galaxy, as item 86 suggests for stars per second) or the galaxy
    database (resets with it), and so what a regenerate or reset does to
    it; and whether item 92's parallel workers count wall time or CPU
    time per task.

95. [ ] **Store each sector's expected and actual density.** Boss
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
    after every fill. Item 94 places each sector in its density bucket
    with these numbers, and items 86, 89 and 93 use the
    expected-versus-actual ratio to correct their estimates. Open
    questions: columns on `sectors` (a migration in the Database
    workstream) or a separate table; whether "actual" counts systems,
    stars, or both; what the decaying average is taken over (the ratio
    per density bucket, so it ties in with item 94, or one galaxy-wide
    figure); whether existing sectors are backfilled by a migration;
    and what happens to the stats when a sector is regenerated
    (item 60) or the galaxy is reset.


104. [ ] **Bug: rogue planets (and maybe other objects) drawn outside
    the sector's wireframe.** Boss (2026-10-01): "rogue planets (and
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

106. [ ] **Give rogue planets a planet class, with a rogue flag in the
    class constants.** Boss (2026-10-01): "Rogue plants should get a
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
    planet's; and the class override and validation items (57-59) treat
    `"r"` as the rogue planet's zone. Open questions: which classes are
    allowed to be rogue (frozen, gas giant and barren classes are the
    obvious ones; does a class with life ever qualify)? Are the
    probabilities `PLANET_CLASS_PROBABILITIES` reweighted for rogues, or
    a separate rogue table? What happens to rogue planets already
    generated: a migration that assigns classes from their stored mass
    and type, or a regenerate? Does a rogue class change its rendering
    (it has no star to light it)?

107. [ ] **Bug: rogue planets are hard to find on the Sector Map.**
    Boss (2026-10-01): "in sector view make sure rogue planets can be
    easily located." Today `static/sectormap.js` draws a rogue planet as
    a dim, dark-purple textured sphere (`roguePlanet`: core `#6b5a8a`
    fading to `#2a2438`, glow `#7d6aa8` at strength 0.8), "a dim,
    starless world lit only by its own internal heat", which nearly
    vanishes against the dark scene. Done: every rogue planet in a
    sector is easy to spot at the default zoom and when zoomed out, in
    both themes, without looking like a star; the sector page's list of
    its contents can point at each one on the map. Goes with item 104
    (rogue planets drawn outside the wireframe), whose bright-color test
    makes the same objects visible, and item 108's point-of-light style.
    Open questions: what makes them findable (a marker or ring around
    each, a brighter but still cool color, a label, a "highlight rogue
    planets" toggle, or a list that flies the camera to each)? Does the
    same apply to other dark objects (quiescent black holes, interstellar
    comets)?

108. [ ] **Stars and glowing phenomena as points of light on the Sector
    Map.** Boss (2026-10-01): "make the stars in a sector more realistic
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
    yes; quiescent black holes and rogue planets, which item 107 must
    keep findable, probably not)?

### Star population (from the galaxy studies of 2026-09-30)

The random star model (mass from the Kroupa IMF, an age, then evolution:
`stellarObjects/stellarEvolution.py`), secondaries from a mass ratio at
the primary's age, planet/life gating by star age and engulfment, the
bright/dim sampling API (`stellarObjects/stellarPopulation.py`),
population densities (`galaxyDensity.population_densities`) and
`SpaceSector.add_preplaced_system` have shipped (bugs S1-S8 of the
project's `galaxy-studies/star-fix-spec.md`, and the Physics part of
`bright-star-preplacement-plan.md`), and so has the fill in
`generate.py` that uses them (population ages and pre-placed bright stars).

### Admin editing: overrides, delete and regenerate (Boss's notes of 2026-10-01)

Boss asked for these on 2026-10-01 (quoted where it matters). None is
designed yet; the open questions are listed in each item.

57. [ ] **Central validate module in `stellarObjects`.** Boss: "we should
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
    system and a star system, which generation and items 58-59 both
    call, with existing behavior unchanged. Prerequisite for 58 and 59.

58. [ ] **Admin override of a planet's or moon's class.** Boss: "it should
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
    item 57's functions, re-spacing orbits until it validates, and warn
    the admin when no stable layout exists. Open questions: does a forced
    change that can't be made stable still save (with the warning), or is
    it rolled back? May revalidation remove other bodies, or only move
    them?

59. [ ] **Admin override of a star.** Boss: "that will change the entire
    system but it will change the system to have as many objects as the
    original system had just their orbital positions will change, caveat
    there is if there are too many objects for the star (say a large star
    with a lot of objects changes to a small star that does not have
    orbital space, then it'll be truncated." Done: an admin control to
    change a system's star; the system keeps its planets, moons and belts
    (same count), orbits are re-spaced for the new star with item 57's
    validation, and outer objects are dropped when the new star lacks the
    room, telling the admin what was removed. Open questions: do the
    planets keep their classes, or are classes re-checked against the new
    star's zones (which could chain into item 58's revalidation)? Does
    this cover companion stars in multiple systems too?

60. [ ] **Delete and regenerate buttons on everything, sector down.**
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
    an audit-log entry, and item 57's validation after a single body is
    removed or regenerated. Open questions: does regenerating a sector
    keep manual overrides (items 58-59), renamed objects and placed
    facilities (item 35), or replace everything? Does regenerating keep
    the object's name? Does deleting a sector leave its slot unfilled
    (so it can be filled again) or mark it empty?

61. [ ] **Lock out an IP address after failed logins.** Boss: "3 failed
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
    workers, and logged to `admin_audit_log`. Open questions: does a
    successful login reset the doubling? Do both limits stay (per IP and
    per username), or does this replace the in-memory one? Should
    localhost or a configured allowlist be exempt so an admin can't lock
    themselves out, and should there be an admin "unlock" command?

### User accounts (Boss's notes of 2026-10-01)

Boss (2026-10-01): "a full user level interface to allow users to
bookmark this will, of course, require an email loop for password setting
/ resetting, invite only, so only an admin can invite a user which is done
by unique link". Today there are only admin accounts: `admin_users`,
`admin_sessions`, `admin_api_keys` and `admin_audit_log` in
`stellarObjects/control_schema.sql`, managed by `stellarObjects/adminAuth.py`.
Install seeds one admin with a random first password and
`must_change_credentials` (`bootstrap_control_schema`), and nothing in the
web interface or the CLI adds another account. There is no email support
and no saved-bookmark feature (pages only have bookmarkable URLs). Order:
64, then 65 and 66, then 67, 68 and 69. Login protection already exists
per username (`src/html/api/loginbackoff.py`, PR #131) and is planned per
IP address (item 61); both must cover user logins, password resets and
invite links too.

64. [ ] **Accounts with roles: user, admin and Owner.** Boss: "Admin can
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

65. [ ] **SMTP settings in the admin config.** Boss: "We'll use SMTP for
    email which means admin config needs SMTP settings." Done: SMTP host,
    port, security (STARTTLS or TLS), username, password and From address,
    set from an admin page (and the config file/installer), with a "send
    test email" button; one small mail module that every flow in 66-68
    uses, which logs failures and never shows the SMTP password. Open
    questions: is the SMTP password kept in the config file (like the
    database password) or in the control database, and is it encrypted
    there? Who can change SMTP settings: any admin, or only the Owner?
    What do invites and resets do when SMTP isn't configured (show the
    link to the admin to pass on by hand?)?

66. [ ] **Invite-only sign-up by unique link.** Boss: "only an admin can
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
    item 67's email loop sets their password; the account is a user, not
    an admin. Open questions: does an admin optionally type the invitee's
    email so the link is sent for them, or only copy the link? Is the
    invite page rate-limited, and is there a cap on open invites? Does
    a multi-use link record who used it?

67. [ ] **Email loop for setting and resetting passwords.** Boss: "an
    email loop for password setting / resetting". Done: a new account
    sets its first password from an emailed link; "forgot password" on
    the login page emails a reset link; links are single-use, short-lived
    (stored hashed) and end the account's other sessions once used; the
    page never says whether an email address has an account; resets are
    rate-limited per address and per IP alongside the existing login
    backoff and item 61. Open questions: how long a reset link lasts
    (30 minutes? 1 hour?); does changing the email address also need an
    email confirmation to the old and new addresses; does the Owner's
    reset need anything extra?

68. [ ] **Owner transfer.** Boss: "The owner CAN (with specific approval
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

69. [ ] **A user-level interface with bookmarks.** Boss: "a full user
    level interface to allow users to bookmark". Done: signed-in users
    (any role) get an account page and can bookmark sectors, systems,
    planets, phenomena and NAV courses, see them in a list, name them and
    remove them; bookmarks are stored per account in the control
    database. Open questions: which objects can be bookmarked, and can
    users also add notes? What else a user can do that an anonymous
    visitor can't (is the site still public to read, or sign-in only?)?
    Do bookmarks survive a galaxy regenerate (object ids change), and if
    not, what does a broken bookmark show?

### Project process (Boss's notes of 2026-10-01)

80. [ ] **Number TODO items by category, and build the version from
    them.** Boss (2026-10-01): "I want to renumber the TODO items in
    groups so like UI changes get 'UX.1' and API Changes at 'API.1' kind
    of thing, so that we can better track changes. We'll revamp
    everything so that tags are consistent through the documentation.
    We'll then build the build number (major feature set.revision.build)
    to be a composite of the change numbers for each category added up.
    (i.e. if we're on UX.4, API.8, and DB.12 we'd add those up to be
    4+8+12)." Nothing is renumbered yet. Today items carry one running
    number across groups (this file's "How to use this document" says to
    renumber when items are added or finished), code sites carry
    `TODO(<area> #N)` tags (`grep TODO(`), and the version (README badge,
    `src/stellarObjects/_version.py`, `CHANGELOG.md`; 7.37.0 as of this
    item) is bumped by `.github/workflows/stamp-version.yml` and
    `scripts/bump_version.py` from each merged PR's
    `changes/<name>.<patch|minor|major>.md` note. Done: every open item
    gets a category ID (`UX.1`, `API.1`, `DB.1`, ...); the same IDs are
    used in `TODO(...)` code tags, `changes/` notes, the changelog, PR
    titles and the design docs; and the version's third number is the
    sum of each category's counter. Open questions:
    - The category list and what each covers (for example UX, API, DB,
      MAP for the Galaxy/Sector/System maps, GEN for generation and
      physics, NAV, SEC for security, OPS for installers and hosting,
      DOC).
    - Is a category's counter the number of changes shipped in it, or the
      highest item ID? Does a finished item keep its ID (no more
      renumbering), so IDs are never reused?
    - How the post-merge Action counts: does each `changes/` note name its
      category and item ID (for example `ux-62.patch.md` or a front-matter
      line), and what happens to a PR that touches two categories or none
      (a pure bug fix)?
    - What "major feature set" and "revision" mean and who bumps them
      (still the `patch`/`minor`/`major` level of the note?). Does the
      build number reset when they go up? It can't, if it's a running sum
      of counters, so the version would only ever grow in its third
      place.
    - Do the 1-79 numbers already in commits, PRs and the changelog get a
      mapping table to the new IDs?

81. [ ] **A structural design document: the program's "circuitry".**
    Boss (2026-10-01): "build a structural design document of how the
    program works overall, bridging file names to what they contain and
    basically lays out the 'circuitry' of the program." Nothing like it
    exists today: `docs/` has reference docs per area (`api.md`,
    `database-schema.md`, `html-interface.md`, `config.md`,
    `testing.md`, ...) and `docs/design/` has topic designs, but nothing
    shows the whole. Done: one document (for example
    `docs/design/architecture.md`) that maps every top-level script,
    package and important module (`generate.py`, `src/stellarObjects/`,
    `src/html/api/`, `src/html/web/`, `src/html/lib/`, `src/html/static/`,
    installers, workflows) to what it holds, and traces the main flows
    through them: generating a galaxy, sector and system; storing and
    migrating the database; serving a page and a map; admin login and
    jobs; releases. Open questions: diagrams (Mermaid, which GitHub
    renders) or text only? How is it kept current (a CI check that every
    module is listed, or a rule that PRs update it)?

82. [ ] **Bring the design documents up to date, with the reasons.** Boss
    (2026-10-01): "clean up the design documents make sure they are all
    current, document how the program works the way it does and why and
    what choices were made that influenced each." Done: every file in
    `docs/design/` and `docs/analysis/` (and the reference docs in
    `docs/`) checked against the code, fixed or marked as historical; each
    says how that part works, why, and which choices and alternatives
    shaped it (for example the cylindrical sector grid, rendering systems
    from the database, three.js for the Galaxy Map, Flask-only site),
    drawing on the PRs and `CHANGELOG.md`. Open questions: do superseded
    designs get deleted or kept in an archive folder? Does this wait for
    item 80's category IDs so the docs are tagged once? Best done after
    item 81, which gives the map to hang them on.

### View from a planet (Boss's notes of 2026-10-01)

**Research first.** Boss: "view-from-planet will have to do calculations
on colors and A LOT Of stuff, so make special note of that, it will need
a full research pass." Before any code for items 83-84, Boss wants a
research session with him "into exactly how one would do that". Items 83
and 84 are blocked on it; item 85 is not.

83. [ ] **A starmap seen from a planet. RESEARCH WITH BOSS FIRST.** Boss:
    "Build a function to select a planet and generate an effective starmap
    from that planet based on all visible stars, this will have to include
    a lot, A LOT, of math so remind me to do research when we get there
    into exactly how one would do that, but it would have to account for
    where each star would have been at that light years back in time?"
    Done: pick a planet (or moon), and get every star visible from it
    with its direction and brightness as seen there, placed where it was
    when the light now arriving left it (light-travel time back along
    its galactic orbit; the correlative update now moves everything along
    galactic orbits, TODO 32, PR #157). The research pass covers at least:
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

84. [ ] **Render the view as a PNG, with constellations.** Boss: "when
    it does that it will generate a PNG and it will generate
    constellations." Done: item 83's view is drawn to a PNG (a sky
    projection, star size and colour by apparent brightness), and the
    brighter stars are grouped into constellations with lines and names
    from item 85, stored so a planet keeps the same constellations each
    time. Blocked on item 83's research pass. Open questions: whole-sky
    or a horizon view from a point on the surface? How are constellations
    chosen (bright-star patterns, by clustering, a set number per sky)?
    Are the PNGs cached on disk and served by the web interface, or made
    on request?

85. [ ] **Constellation names in the name generator.** Boss: "add to our
    name generator constellation name support based on constellation
    names throughout all known languages and then slice it up like we do
    for all our naming". Done: a constellation name list in
    `stellarObjects/names.py` gathered from constellation and star-group
    names across the world's languages and sky cultures (not just the 88
    IAU ones), and a constellation name generator that slices and
    recombines them into new names the same way stars, planets and
    sectors are named (`split_into_syllables` in `utils.py`, the
    prefix/suffix lists, the `offensive_words.txt` filter). Used by item
    84. Open questions: what counts as a source list (licensing of sky
    culture data such as Stellarium's), transliteration of non-Latin
    scripts, and whether names are unique per planet or galaxy-wide.

## Population and Politics

Items 51-54 shipped: schema v44 (`generate.py population`, the
`/api/species`, `/api/polities` and `/api/territories` endpoints; see
`docs/design/population-and-politics.md`), the Galaxy Map's Territories
overlay, and the pages (Species and a species page, Polities and a
polity page, "Dominant species" on a life world and "Territory of ..."
on an owned system), all hidden until population data exists.
