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
- **Galaxy Map (12-19):** Boss approved the plan in the
   project's `galaxy-megablocks/report.md` (hybrid master-wedge
   slots, pixel-sized mega-blocks). Work item 18 next (13-17 have shipped: pixel-sized blocks, the solid and its slice, filled and unfilled blocks with no marker dots, block info, and smooth zooming). 12 (the
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

18. [ ] **Keep three.js; record why.** It was checked on 2026-09-30:
   - Babylon.js is several MB, and deck.gl needs a bundler.
   - regl and raw WebGPU would mean rewriting picking, sprites and
     lighting by hand.
   - The CSP (`default-src 'self'`) and the no-build-step vendoring rule
     favor one vendored file.
   - The bottleneck is JavaScript listing work, not the renderer.

   three r186 already has InstancedMesh, BatchedMesh and a
   WebGPURenderer to move to later. Done means the rendering choice is
   written into `docs/html-interface.md`.

19. [ ] **Follow-ups (edge cases).**
   - Distance-based detail (bigger blocks farther from the camera), which
     the aligned wedges from 5 make seamless.
   - Order-independent transparency (weighted blended) if #15's sorting
     shows artifacts where translucent blocks intersect.
   - Phone performance at 390 px. Neighbouring full-size blocks share
     faces, so skipping a face whose neighbour exists would cut the
     vertex count; past that, one InstancedMesh per wedge-arc count.
   - DPR: `pcPerPixel` is per CSS pixel.


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
    - Left: the Sector Map draws each cloud's extent (`sectormap.js`) and
      the Galaxy Map shows them (`galaxymap3d.js`).

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
`bright-star-preplacement-plan.md`). What's left:

55. [ ] **Use population ages at sector fill (S7 call site).** In
    `generate.py`'s `generate_sector`, set each system's
    `SystemConfig.POPULATION` to
    `stellarPopulation.pick_population(galaxyDensity.population_densities(position_pc, shape))`
    for a galaxy-placed sector, so O/B stars and supergiants sit in the
    arms near the plane and the bulge has none (Database/Web own that
    file; the bright-star fill sets `MAX_STAR_LUMINOSITY_SOL` the same
    way).

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
