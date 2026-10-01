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
   slots, pixel-sized mega-blocks). Work items 15-18 in order (13, pixel-sized blocks, and 14, the solid and its slice, have shipped). 12 (the
   hybrid master-wedge slot rule) shipped in schema v35. 19 is follow-ups.
- **Features (23-36)** from the same notes: Galaxy Map generate buttons (24),
   phenomena views and stored nearest systems (25-26), nebulae and
   remnants: placement, classes, containment and naming, plus asteroid
   field classes (27-31), the correlative update (32), navigation frames
   and speeds (33-34), and facilities (35-36). 26, 27, 29, 30 and 35
   are schema changes; 28 and 31 (classes) shipped in schema v38.
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

15. [ ] **One solid of blocks for filled and unfilled sectors; no more
   marker dots.** Boss asked for this on 2026-09-30. The goal is to zoom in
   and out and find generated ("filled") sectors from the blocks alone.
   - **Remove the dots.** Drop the placed-sector sprites (the halo and
     core dots) and the planned-sector dots from `galaxymap3d.js`: the
     textures, `syncTier`, `withPinned`, marker scaling, and marker picking.
     Filled and unfilled sectors are both shown only through the
     continuous solid of blocks (#13, #14).
   - **Color by density.** Every block is colored by the density of the
     space it covers (`prismShade`). At m = 1, a filled sector is colored by its
     real system density (`placedDensityColor`'s scale).
   - **Opacity by how full a block is.**
     - Unfilled space is 20% (densest) to 50% (sparsest) transparent,
       scaled by density.
     - A block holding filled sectors grows more solid in proportion to
       its filled share (filled ÷ total sectors), reaching fully opaque when
       every sector is filled.
     - At large m the share is tiny (a handful of filled sectors among
       531,441), so give any filled content a minimum visible step, then
       scale it. A log of the count is one option.
     - Boss confirmed on 2026-09-30: the more filled sectors a block
       holds, the more solid it is, and fully solid once every sector is
       generated.
   - **Individual filled sectors appear only at sector zoom (m = 1).**
     Coarser, they show only through their block's opacity.
   - **Picking moves to blocks.**
     - Clicking a block shows its info (#16).
     - At m = 1, a filled block links to its sector page, and an unfilled
       one shows today's designation and CLI snippet.
     - Double-clicking a block with filled sectors zooms in toward them.
   - **What it needs:**
     - Per-block filled counts, counted in the browser from the tiles'
       placed lists, or totalled by the server for coarse views.
     - Translucent blocks drawn after opaque ones, sorted back to front.
     - Interior culling (#13) only where all neighbours are opaque.
     - A block's total sector count (`blockSectorCount`).

   Done means filled sectors can be found by zooming alone at every zoom,
   there are no marker sprites left, colors follow density, and the frame
   rate holds at the #13 block budget.

16. [ ] **Block info on click.** Clicking a block shows its sector ring,
   layer and slot ranges and its exact sector count (and how many are
   generated once tiles carry that). A click at m = 1 keeps today's sector
   panel (designation, CLI snippet, 8 corners).

17. [ ] **Smooth zooming: preload and prerender.** Today each zoom step
   rebuilds the whole prism set on the main thread, then waits on tiles.
   - Move block listing and geometry into a Web Worker. `galaxyprisms.js`
     already has no three.js import. Hand back transferable typed arrays.
   - Keep the geometry for the current m and the next finer and coarser
     m ready ahead of time, keyed by (m, slice, center cell).
   - Prefetch the tiles the next zoom step will need.
   - Animate wheel and button zoom over a few frames, and crossfade between
     block sizes instead of popping.

   Done means no dropped frames while zooming on a mid-range laptop, and
   no visible wait at any zoom step already visited.

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
   - Phone performance at 390 px.
   - DPR: `pcPerPixel` is per CSS pixel.
   - Reduced-motion users get instant zoom.


### Web interface (`src/html/web/`, `src/html/static/`)

### Sector Map and generation (`static/sectormap.js`, `web/generate_page.py`, `generate.py`)

24. [ ] **Generate buttons on the Galaxy Map's unfilled sectors.** Boss:
    "Sector map clicking on an unfilled sector should no longer give a command line
    but if admin is logged in then it should just add a button to
    generate that sector by itself or as a neighborhood or to generate
    the entire shell (not recommended), also let's add a 'generate
    column' option too." The Sector Map has them (generate.py's
    `--column`/`--shell`/`--slot --radius-pc`, the Generate page's
    column and shell modes, `sectormap.js generateButtons`). Left: the
    Galaxy Map's `showPlannedInfo`/`showCellInfo` in `galaxymap3d.js`
    still show the CLI snippet; give an admin the same four buttons
    (the page needs the same admin-only `generate` target
    `starmap.render_map_panel` gets). Visitors see the address and
    designation only.

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

26. [ ] **Octant and nearest systems for phenomena, stored.** Boss:
    "Phenomena in the sector list should also list the octant they are
    in. Also I'd like what stars are nearest to the phenomena at that
    moment and that should be pre-calculated and stored in the database,
    and in fact if we aren't already precalculate and store in the
    database the nearest 3 star systems to each, even if it crosses
    sector boundaries."
    - `sector_page._contents` sets phenomena's octant to None; use
      `spaceSector.classify_octant` (store it like
      `star_systems.quadrant`).
    - Only a text summary of a system's nearest 3 is stored today
      (`star_systems.location`, from `SpaceSector.nearest_neighbors`,
      within the sector only). Add a `nearest_systems` table (object
      kind and id, rank 1-3, neighbor system, distance) for every system
      and phenomenon, searched across sector boundaries, filled at
      generation and by #32. Schema change with a migration.
    - Show them on the sector page, system page and phenomenon page;
      `queryDb.phenomena_near_sector` returns them.

27. [ ] **Generate nebulae and supernova remnants, with the stars they
    need, and put them on the maps.** Boss (2026-09-30): "Nebulae and
    Remnants should be generated and placed on the map. Research if we
    need stars at the center of these or what kind of star, etc, so we
    can make them."
    - Rates come from `PHENOMENON_DENSITY_PC3` (v37; a per-volume rate suits objects this big).
    - Each class brings its central object (`NEBULA_CLASSES[...]["center"]`,
      table in `docs/design/nebula-and-asteroid-field-classes.md`): O/B stars for
      emission nebulae, a B or A star for reflection, one hot central
      star becoming a white dwarf for planetary, none for molecular
      clouds (protostars at most), a neutron star or black hole for
      core-collapse remnants and none for thermonuclear ones. The
      generator creates that star system inside the nebula, or places
      the nebula around a qualifying existing star.
    - Nebulae are up to 200 ly in radius, so one spans many 13 ly
      sectors: every sector it reaches lists it, the Sector Map draws
      its extent, and the Galaxy Map shows it.
    - Sites: `generate.generate_sector_phenomena`,
      `nebulaData.Nebula`, `supernovaRemnantData.SupernovaRemnant`,
      `queryDb.phenomena_near_sector`, `sectormap.js`, `galaxymap3d.js`.

29. [ ] **Record what sits inside a nebula.** Boss: "Add a DB field for
    if any stellar object (including systems) exist within a nebulae or
    similar (not asteroid fields, that wouldn't work) or stellar
    remnants if necessary."
    - A nullable `nebula_id` (the innermost containing nebula or
      supernova remnant) on `star_systems`, `rogue_planets`,
      `interstellar_comets`, `black_holes`, `neutron_stars`,
      `asteroid_fields`, `nebulae` (nesting) and stand-alone facilities
      (#35). Asteroid fields can sit inside a nebula but never contain
      anything. Schema change with a migration.
    - Containment is a 3D distance test against every nebula that
      reaches the object's sector, set at generation, when a later
      sector is generated inside an existing nebula, and by #32.
    - Show "inside <nebula>" on system, phenomenon and sector pages.
      Being inside also compresses a star's heliopause (down to ~0.2 AU
      in a dense cloud), which should feed system text, habitability and
      the navigation hand-off radius (#33).

30. [ ] **Names that follow one standard.** Boss: "Asteroid fields and
    comets should be named using a method that tells something about
    them by their name in letters and numbers in a standardized way.
    Nebulae should get names the same as star systems do, as do neutron
    stars, quasars, black holes, etc."
    - Nebulae, remnants, neutron stars, black holes, quasars and rogue
      planets go through the system-name registry
      (`_db.reserve_system_name`/`confirm_system_name`) instead of an
      unregistered `generate_phoneme_salad_name`.
    - Comets and asteroid fields get designations. Draft (IAU-style):
      `P/<system>-<n>` periodic comet, `C/<system>-<n>` long-period,
      `I/<sector designation>-<n>` interstellar comet, and
      `AF <class><size digit>-<sector designation>-<n>` asteroid field
      (`asteroid_fields.field_class`, schema v38). Examples are in the design
      doc. Needs a migration that renames existing rows.

### Correlative update (`src/updateOrbits.py`, `stellarObjects/_db.py`)

32. [ ] **Check and finish the correlative update.** Boss: "We need to
    double check our 'correlative update' to move everything (planets,
    moons, stars, phenomena, etc) in their orbital path, we had a script
    for it but let's make sure it's working for the current version.
    Each time this is run it should recalculate the nearest systems and
    store it in the database after it calculates the new galactic
    location."
    - `updateOrbits.py` advances phases through
      `_db.advance_orbital_phases`/`advance_comet_orbits`: planets,
      moons, stars, binaries, standalone phenomena's galactic phase and
      comets. No stale columns were found against schema v34, but it
      needs a real run against a v34 database.
    - It never moves anything's galactic position:
      `sectors.center_*_pc`, `star_systems.position_*_mpc`, phenomena
      `center_*_pc` and sector membership stay where they were. Move
      them along their galactic orbits, then recompute the nearest
      systems (#26).
    - Quasars aren't in the phenomena loop.
    - When an object's orbit carries it out of its sector it moves to
      the new sector. Boss (2026-09-30): "when a sector changes then we
      make sure the DB and all text is changed to point at the new
      sector location": `sector_id`, sector-relative positions, octant
      (`quadrant`), `star_systems.location`, containing nebula (#29),
      and any stored or rendered text naming the old sector.

### Facilities (new)

35. [ ] **Starbases, colonies and outposts in the database.** Boss: "I
    am going to have starbases, colonies, outposts, that kind of thing.
    Terrestrial facilities, orbital facilities, and stand-alone
    facilities (those parked in space)." Add them to the database and
    wire up adding a facility to a location, with these rules:
    - Gas giants can't have terrestrial facilities (orbital only).
    - Terrestrial worlds can have colonies (which automatically make
      the planet inhabited: OR it into `queryDb._with_life_fields`) and
      orbital facilities.
    - Asteroid belts can have asteroid facilities (an asteroid outpost
      or mining colony); asteroid fields can have asteroid outposts.
    - A star system can have an outpost in orbit around the star itself.
    - Stand-alone facilities are parked in space (a galactic position,
      like a phenomenon).
    - Orbits (around a star or a planet) get distance, period and speed
      from the host's approximate mass, the same way everything else
      does (`planetPhysics.calculate_orbital_period_years`,
      `utils.circular_orbital_speed_kms`), and move with #32.
    - Schema: one `facilities` table (kind, name, exactly one host,
      orbit columns), migration, API routes in `html/api/routes.py`.
    - Moons can host terrestrial and orbital facilities (Boss,
      2026-09-30).

36. [ ] **Place facilities from the web interface.** Boss: "The web
    interface should have a way to select within a star system where a
    facility goes in orbit around the star or around the planet, which
    will have calculated distances and orbital speeds the same way
    everything else does, based on the approximate mass." An admin form
    on the system page (`web/system_pages.system`), showing the
    calculated distance and speed before saving; facilities listed on
    the system page and drawn on the System Map; stand-alone ones on the
    sector page.

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
