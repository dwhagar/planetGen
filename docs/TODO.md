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
   70-72 first, in order; it replaces the map's click-to-center and
   double-click zoom.
- **Galaxy Map (12-19):** Boss approved the plan in the
   project's `galaxy-megablocks/report.md` (hybrid master-wedge
   slots, pixel-sized mega-blocks). 13-18 have shipped (pixel-sized blocks, the solid and its slice, filled and unfilled blocks with no marker dots, block info, smooth zooming, and the three.js decision in `html-interface.md`). 12 (the
   hybrid master-wedge slot rule) shipped in schema v35. 19 is follow-ups.
- **Features (25-36)** from the same notes: phenomena views and stored nearest systems (25-26), nebulae and
   remnants: placement, classes, containment and naming, plus asteroid
   field classes (27-31), the correlative update (32), navigation frames
   and speeds (33-34), and facilities (35-36). 28 and 31 (classes)
   shipped in schema v38, 29 (containment) in v39, 30 (naming) in v40,
   26's storage in v41 and 35 (facilities) in v42; 32 (the
   correlative update moves everything) shipped with them, and so did
   33-34 (navigation frames and speeds).
- Each change site in the code carries a `TODO(<area> #N)` comment
   naming its item here (areas: distances, system-list, site-header,
   search, phenomena, galaxy-map, sector-map, orbits, nav, facilities,
   security, physics, web-pages);
   grep for `TODO(` to see them all, or `TODO(galaxy-map` for one area.
- Items with a **Question for Boss** state the default taken; the work
   can start on that default.

### Performance

8. [ ] **Add a cache so pages don't hit the database on every request.**
   Every page calls the Flask API through `html/lib/apiclient.py`,
   and every API route queries MySQL fresh, including results that rarely
   change (`/api/galaxy/sectors`, `/api/galaxy/shape`, sector and system
   detail). The 3D Galaxy Map already has one: its cube tiles are cached
   on disk by the web layer (`html/lib/tilecache.py`) and in the browser,
   and an edit refreshes only the tiles it touched (`GET /api/galaxy/
   changes`, from the v27 `modified_at` columns). Extend that pattern to
   the other pages, invalidating each page from its own rows'
   `modified_at`.

### Galaxy Map and the sector standard (`static/galaxyprisms.js`, `static/galaxymap3d.js`, `lib/galaxymap3d.py`, `stellarObjects/galaxyGeometry.py`)

Today the map draws the analytic density as shrunk prisms, m sectors a
side (m a power of 3), sized by a volume budget that badly overestimates
the thin disk. The result is 70-290 px cubes with gaps, and the spiral
barely shows. The plan (report above, with renders) replaces that with a
continuous solid of mega-blocks sized from the screen's pixel scale.

19. [ ] **Follow-ups (edge cases).**
   - Distance-based detail (bigger blocks farther from the camera), which
     the aligned wedges from 5 make seamless.
   - Order-independent transparency (weighted blended) if #15's sorting
     shows artifacts where translucent blocks intersect.
   - Phone performance at 390 px. Neighbouring full-size blocks share
     faces, so skipping a face whose neighbour exists would cut the
     vertex count; past that, one InstancedMesh per wedge-arc count.
   - DPR: `pcPerPixel` is per CSS pixel.


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
item names its section. Items 70-72 go in order (Galaxy Map thread); Web
can do 63 and 74 alongside, then 75 and 78 once 72 fixes the URLs.

70. [ ] **Nested ladder geometry (design doc section 3).** Today's
    `blockWedgeCount` picks each level's wedges on its own, so a child
    block sits inside one parent only sometimes (109 of 143 rings at
    243 -> 27, 1,162 of 1,286 at 27 -> 3). Done: the nested wedge rule
    (each level's wedge count a whole multiple of its parent's),
    sectors joining the level-3 block that holds their center, and
    `drillWedgeCount`, `drillParent`, `drillChildren`, `drillSlabs`,
    `drillBlockSectors`, `drillChainOf` in `static/galaxyprisms.js` with
    the same functions in `stellarObjects/galaxyGeometry.py`; a node
    parity test shows both agree and that a parent's sectors are exactly
    its children's.

71. [ ] **Stage contents API (section 7).** Done: `GET
    /api/galaxy/stage?at=m.I.s.S` returns one container's children with
    generated counts (generated sectors listed at m = 3), and the whole
    galaxy's level-243 blocks with no `at`; cached by `lib/tilecache.py`
    under the stamp, a change invalidating only its ancestor chain. No
    schema change.

72. [ ] **The drill-down stages (sections 4, 5, 8.1, 10).** Done: the
    eight stages on the Galaxy Map, with slab hover highlight, pull-out
    to a top-down view, the van Wijk-Nuij flight into a block, the slab
    strip, the breadcrumb with sibling menus, stage URLs
    (`/galaxy?slab=`, `?at=`, `?sector=`) with Back/Forward, keys,
    touch taps and the "Generated only" toggle; reduced motion cuts
    instead of animating. Click-to-center and double-click zoom go away
    (decision 2 in section 11 decides whether free look stays). Needs
    63's bigger map for room.

73. [ ] **Generate from the sector level (section 6).** Boss: "once
    we're down to a sector level we can tell a slice to generate all the
    sectors in that slice or click on a sector and generate it from the
    UI if you're admin", and "add a 'generate neighborhood' when at a
    sector selection level that will ask the radius in ly." Done, for an
    admin at stages 7-8: Generate this sector (today's `slot` mode),
    Generate this layer or slab (a new `generate.py galaxy --block
    m.I.s.S [--block-layer j]` mode plus a Generate page form), and
    Generate neighborhood with a light-year radius dialog (default 100
    ly, 13-652 ly, converted with `ly_to_pc`, an "up to about N sectors"
    estimate, a confirmation above 5,000), started without leaving the
    map and refreshed when the job ends. Web owns `generate.py` and the
    Generate page; Galaxy Map owns the buttons.

74. [ ] **Sector Map pick mode and Nav links (sections 9.1, 9.2).**
    Done: `/sectors/<id>?pick=from|to&...` shows a banner and a "Use as
    start/destination" button on a system or phenomenon, which lands on
    `/nav?from=...&to=...`; system and phenomenon pages and the Sector
    Map panel get "Nav from here" and "Nav to here".

75. [ ] **NAV page picks on the map (section 9).** Boss: "from the nav
    menu select start and destination using either the text dropdowns as
    we have now or the galactic map interface to select. If it's within
    sector then it'll just use the sector interface." Done: beside each
    dropdown, "Pick on Galaxy Map" (`/galaxy?pick=...`, generated-only
    forced on, ending in 74's Sector Map pick mode), "Pick in this
    sector" once the other end is known, and a Bookmarks select.

76. [ ] **Bookmarks (section 8.2).** Done: a ☆ on the breadcrumb and
    info panels saves a stage, sector, system or phenomenon in
    `static/bookmarks.js` (per browser, up to 100, storage failures
    tolerated), with a map menu, Ctrl+1-9, rename and delete, and the
    entries offered by the NAV pickers. Shared bookmarks need Boss's
    decision 4 and a migration.

77. [ ] **Address bar (section 9.3).** Done: a field over the breadcrumb
    takes a designation, `ring/layer/slot`, `x, y, z` pc or a name and
    flies to that sector's stage 8, in pick mode too.

78. [ ] **"Show on Galaxy Map" links (section 8.1).** The sector page's
    link goes to the Quadrant table today. Done: sector, system and
    search pages link to `/galaxy?sector=<designation>`.

79. [ ] **NAV course on the Galaxy Map (section 9.4).** Done: the NAV
    result's "Show on Galaxy Map" opens the smallest stage holding both
    endpoints with the course drawn. After 72 and 75.

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

63. [ ] **A bigger Galaxy Map with its controls underneath.** Boss
    (2026-10-01): "I also want the galaxy map box to be bigger, place the
    controls under it horizontally if possible, stacked if not, but use as
    much of the browser area as is reasonable to display the galaxy map."
    Today the map is a square capped at 36rem
    (`.galaxymap3d-panel .starmap-viewport` in `static/style.css`) inside
    the 72rem main column (`.app .app-main`), with the controls and info
    panel in a side column (`.starmap-side`, built in
    `lib/galaxymap3d.py`). Done: the Galaxy Map page lets the map use most
    of the browser window (wider than the 72rem column, and as tall as
    the window allows after the header, not forced square), the controls
    sit in a row under the map and wrap to a stack when the row doesn't
    fit, and the canvas resizes with the window (`static/galaxymap3d.js`
    must follow the new size; the Galaxy Map thread owns that file).
    Open question: does the block info panel go under the controls, or
    stay beside the map on wide screens?

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

### Phenomena (`lib/phenomenonmap.py`, `web/system_pages.py`, `web/sector_page.py`, `generate.py`)

25. [ ] **A view that suits each phenomenon.** Boss: "view for neutron
    stars should be not a 3D map but rather a rendered representation of
    the neutron star pulsing in rough time with its properties and show
    some way to show it's spinning from the render. Same for comets,
    rogue planets, etc, asteroid fields don't get a 3D render at all,
    black holes should get a 3D render representing their accretion
    disk around it. Similar for Quasars we should see a 3D rendering
    similar to what it might look like."
    - Today `render_phenomenon_map_panel` draws a flat SVG: a circle of
      `radius_ly` for nebulae, asteroid fields and remnants, a dot for
      everything else.
    - Neutron star: pulse in rough time with `spin_period_ms` (slowed to
      a visible rate, stated on screen), beams or a surface feature so
      the spin reads; non-pulsing ones just rotate.
    - Comet, rogue planet: a rendered body (the tail for a comet).
    - Asteroid field: no render.
    - Black hole: a three.js accretion disk (tilt, glow; intermediate
      and stellar sizes differ). Quasar: the disk plus jets when
      radio-loud.
    - `prefers-reduced-motion` gets a still frame; pages still read
      without JavaScript.
    - Nebulae and supernova remnants: Boss (2026-09-30) wants them
      generated and placed on the maps (#27); their own view keeps a
      map until a render is designed for them.

26. [ ] **Show phenomena's octant and everyone's nearest systems.**
    Storage shipped in schema v41 (2026-09-30): every placed phenomenon
    has `quadrant`, and `nearest_systems` holds the 3 nearest star
    systems to every placed system and phenomenon, searched across
    sector boundaries (`_db.refresh_nearest_systems`, filled at
    generation and refreshed by the correlative update, `updateOrbits.py`).
    `queryDb.phenomena_near_sector` returns `octant` and `nearest`,
    `sector_detail` returns `nearest` per system, and
    `queryDb.nearest_systems(conn, table, ids)` serves any page.
    - Left: show them on the sector page (`sector_page._contents`'s
      `octant`/`location`), system page and phenomenon page.

27. [ ] **Put nebulae and supernova remnants on the maps.** Generation
    shipped (2026-09-30): sectors now generate molecular clouds,
    planetary nebulae around their own hot white dwarf, H II regions
    around O and early-B stars and reflection nebulae around later B and
    A stars (`generate.add_star_hosted_nebulae`,
    `program_constants.NEBULA_HOST_RULES`), and a remnant's core drifts
    off-center by its birth kick. `queryDb.phenomena_near_sector` already
    lists every cloud that reaches a sector.
    - The Galaxy Map draws them (2026-10-01): each tile lists the clouds
      reaching into it (`queryDb.galaxy_clouds_in_box`).
    - Left: the Sector Map draws each cloud's extent (`sectormap.js`).

### Facilities (new)

36. [ ] **Place facilities from the web interface.** Boss: "The web
    interface should have a way to select within a star system where a
    facility goes in orbit around the star or around the planet, which
    will have calculated distances and orbital speeds the same way
    everything else does, based on the approximate mass." An admin form
    on the system page (`web/system_pages.system`), showing the
    calculated distance and speed before saving; facilities listed on
    the system page and drawn on the System Map; stand-alone ones on the
    sector page.
    - The database side shipped in schema v42 (2026-09-30): the
      `facilities` table, the rules in `program_constants.FACILITY_RULES`
      (`stellarObjects/facilities.py`), `_db.add_facility`, and the API:
      `POST /api/facilities`, `DELETE /api/facilities/<id>`,
      `GET /api/facilities/<id>`, `GET /api/systems/<id>/facilities`,
      `GET /api/sectors/<id>/facilities` and
      `GET /api/facilities/orbit?host_type=&host_id=&distance_km=` (the
      orbit to show before saving). See `docs/api.md`.
    - A colony makes its world inhabited: OR
      `queryDb.colonized_body_ids` into `queryDb._with_life_fields`.

### Web API (`src/html/api/routes.py`)

Low priority; nobody is waiting on these.

37. [ ] **The API can't create a system inside an existing sector.**
    `POST /api/systems` only creates standalone systems (`sector_id =
    NULL`, see `docs/api.md`). Attaching one to a sector needs the sector's
    placement and Hill-sphere separation logic (`SpaceSector.add_system`),
    which was left out of the write API to keep the admin-auth change
    small.

38. [ ] **The API can't edit a system's generated content.** `PATCH
    /api/systems/<id>` only renames. Changing stars/planets/moons/belts
    means `DELETE` then `POST` (regenerate). It may never need solving;
    kept here in case it does.

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

49. [ ] **System Map names never overlap.** Boss: "we need to make sure
    names on the system map clickable interface do not overlap."
    - Today `systemmap._label_sides_2d` places each label (4 directions,
      then a pushed "below"/"above" with a leader line, else dropped)
      against the others (plus seeded star-label rects) using an
      estimated width (`_label_half_width_px`: character count times a
      fixed width). Real text can run wider than the estimate, so
      labels can still collide.
    - Fix: make sure every star label and marker is in the collision
      set; measure the real text in the browser
      (`getBBox()` in `systemmap.js` after load and after each zoom
      step in `mapzoom.js`) and nudge or hide labels that still
      overlap, keeping the server placement as the no-script fallback.
    - Check every scene: single star, close and wide binaries, and the
      moon-centered scenes, at 390 px and 1280 px.

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

Exploratory ideas, not yet designed. Each needs a design pass before it
can be ordered against the work above.

51. [ ] Assign government ownership to star systems so that groups of
    systems form territories mapped in 3D space.
52. [ ] Flag worlds with life for generated names of their dominant
    species.
53. [ ] A database of spacefaring species.
54. [ ] Model younger and older civilizations: what differs with a
    society's age and how to store and present it.
