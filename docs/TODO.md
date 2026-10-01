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
   correlative update moves everything) shipped with them.
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

### Sector Map and generation (`static/sectormap.js`, `web/generate_page.py`, `generate.py`)

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

### Navigation and travel (`stellarObjects/navigation.py`, `queryDb.nav_between`, `web/nav_page.py`, `templates/nav.html`)

33. [ ] **Courses in "bearing mark mark" format on nested reference
    frames.** Boss: "Course projections should be in the format of:
    0-359 mark 0-359 with 0 mark 0 pointing toward the galactic core."
    Boss's design (summarized; the full text and pseudocode are in
    `docs/design/navigation-frames.md`):
    - North points toward the local dominant center of mass. Every
      local frame is a rigid transform of the absolute galactic
      Cartesian frame.
    - Galactic Standard Frame: origin the galactic core, +Z galactic
      north, +X a fixed zero meridian. Between sectors, bearing 000
      points at the core.
    - Sector Local Frame: origin the sector's barycenter (4 pc cell),
      North from the ship toward it, +Z the galactic +Z.
    - System Local Frame: origin the central star/barycenter, North
      from the ship toward it, +Z the star's net angular momentum
      (the ecliptic normal).
    - Math: D = target - ship; U = the plane's normal; N = (center -
      ship) with its U part removed, normalized; E = N x U. Bearing =
      atan2(D·E, D·N) in [0, 360), 000 = North, 090 = East. Mark =
      atan2(D·U, sqrt((D·N)² + (D·E)²)). Directly over the center pole,
      fall back to a fixed reference vector. Boss's `compute_course`
      pseudocode in the design doc is the reference.
    - Hand-offs: star to sector barycenter past the heliopause (~120
      AU); galactic frame when crossing a sector boundary (> 4 pc).
    - Today `navigation.course_between` returns azimuth/altitude on the
      galactic plane and `nav.html` shows them as separate rows.
    - **Questions for Boss:**
      - Marks from 0-359: the math gives -90 to +90. Default taken, from
        your own note: write mark as elevation mod 360, so 000-090 is up
        and 270-359 is down (270 = straight down), and nothing between
        091 and 269 appears. OK?
      - "0 mark 0 toward the galactic core" holds in the galactic frame;
        inside a system or sector, 0 points at the star or sector
        barycenter. Default taken: 0 mark 0 is toward the current
        frame's center.
      - The zero meridian "toward a reference quasar": here the quasar
        sits at the galactic center, so it can't set +X. Default taken:
        keep the galaxy's existing +X axis (ring slot 0).

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

### Installers and platforms (Boss's notes, 2026-09-30)

50. [ ] **PowerShell install and upgrade scripts, and bash scripts that
    also run on macOS.** Boss: "we need to write a powershell install and
    upgrade scripts as well as make sure our bash shell scripts will also
    work on macos as well as linux."
    - **Windows:** `install.ps1` and `update.ps1`, the counterparts of
      `install.sh` and `update.sh`: the same steps and the same prompts
      (check-only upgrades, the migrate-or-delete database prompt with
      its 30-second default, the migration progress bar), using a venv
      or the Windows Python launcher instead of apt, Windows services or
      Task Scheduler instead of systemd timers, and Apache on Windows
      (or IIS) paths and permissions (`icacls`) instead of `www-data`
      and `chown`.
    - **macOS:** every bash script (`install.sh`, `update.sh`,
      `scripts/deploy-common.sh`, `scripts/install-python-deps.sh`, and
      `examples/apache/*.sh`, `examples/maintenance/*.sh`) must run on
      macOS too. Known gaps: macOS ships bash 3.2, so
      `install-python-deps.sh`'s `declare -A` and `mapfile` fail there
      (require Homebrew bash, or rewrite them); apt is assumed (use
      Homebrew, or pip in a venv); systemd timers and `systemctl` (use a
      launchd plist); Debian Apache layout (`/etc/apache2`, `a2enmod`,
      `www-data`) versus Homebrew's (`/opt/homebrew/etc/httpd`, `_www`);
      logrotate (use newsyslog); and BSD versus GNU flags in `sed`,
      `stat`, `readlink`, `date` and `timeout` wherever they appear.
    - Keep the steps in step across the three platforms, so a change to
      one installer lands in all of them.
    - The Windows and macOS hosting guides (being written in `docs/` by
      the docs thread) describe the server setup; this item is only the
      scripts, and the guides should point at them once they exist.
    - Each script carries a `TODO(installers #50)` comment at its top.

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
