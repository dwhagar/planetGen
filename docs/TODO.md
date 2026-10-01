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

- **Bug fixes (4-7)** come first, from Boss's notes of 2026-09-30. 4 is
   a small web change;
   5-7 change generation constants; the frequency research for 5 and 6
   is in `docs/design/interstellar-object-rates.md`.
- **Extend the cache (8)**.
- **Galaxy Map (12-19):** Boss approved the plan in the
   project's `galaxy-megablocks/report.md` (hybrid master-wedge
   slots, pixel-sized mega-blocks). Work items 17-18 in order (13-16 have shipped: pixel-sized blocks, the solid and its slice, filled and unfilled blocks with no marker dots, and block info). 12 (the
   hybrid master-wedge slot rule) shipped in schema v35. 19 is follow-ups.
- **Features (23-36)** from the same notes: generate buttons (23-24),
   phenomena views and stored nearest systems (25-26), nebulae and
   remnants: placement, classes, containment and naming, plus asteroid
   field classes (27-31), the correlative update (32), navigation frames
   and speeds (33-34), and facilities (35-36). 26, 27, 30 and 35
   are schema changes; 28 and 31 (classes) shipped in schema v38, 29
   (containment) in v39.
- **More pages (47)**: the Sector Map wireframe is small and can go in any time.
- Each change site in the code carries a `TODO(<area> #N)` comment
   naming its item here (areas: distances, system-list, site-header,
   search, phenomena, galaxy-map, sector-map, orbits, nav, facilities,
   security, physics, web-pages);
   grep for `TODO(` to see them all, or `TODO(galaxy-map` for one area.
- Items with a **Question for Boss** state the default taken; the work
   can start on that default.

### Bug fixes (Boss's notes, 2026-09-30)

4. [ ] **Tag search: collapsible groups and phenomena.** Boss: "Search
   by tag should have collapsible zones for each group of tags so it
   isn't overwhelming, add stellar phenomena types to its search list.
   I can't find any nebulae on the phenomena page and they aren't
   searchable."
   - Done: each "Browse by Tag" group is a `<details>`, open when one
     of its tags is active. Left: the phenomenon facet, after item 28.
   - `searchpage.TAG_FACETS`/`RESULT_PANELS` and
     `queryDb.SEARCH_TAG_FACETS`: add a phenomenon-type facet and a
     Phenomena result panel (none exists today).
   - Nebulae are missing because they're almost never generated, not
     because the page hides them: see `docs/design/interstellar-object-rates.md` (the v37 rates).

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

23. [ ] **Generate a column and a shell.** Boss defined them on
    2026-09-30: "For the shell I mean a ring through every layer and a
    column through the ring." Add two `generate.py galaxy` modes:
    - Column: every sector of one (ring, slot), all layers between
      `galaxy_column`'s `layer_index_min`/`max` for that ring.
    - Shell: every slot of one ring through every layer (a cylindrical
      shell). Far larger than today's ring batch (one ring at one
      layer), so it keeps the `--limit`/`--yes` guard.
    Today's modes are ring batch, one slot, local neighborhood and
    random start (`add_galaxy_arguments`, `run_galaxy`).

24. [ ] **Generate buttons on unfilled sectors.** Boss: "Sector map
    clicking on an unfilled sector should no longer give a command line
    but if admin is logged in then it should just add a button to
    generate that sector by itself or as a neighborhood or to generate
    the entire shell (not recommended), also let's add a 'generate
    column' option too so I can generate an entire column of sectors."
    - `sectormap.js showNeighborInfo`/`cliSnippet` and the Galaxy Map's
      `showPlannedInfo`/`showCellInfo`: drop the snippet; for an admin,
      show the four buttons (shell marked not recommended). Visitors see
      the address and designation only.
    - `generate_page.GALAXY_MODES`/`galaxy_argv`: add column and shell
      (#23); the buttons post a CSRF-protected form that
      starts the job with the address filled in.

55. [ ] **Star dots sized to real giants and white dwarfs.** From the
    star-type study (2026-09-30): `lib/starmap.py _star_dot_radius`
    caps at 14 px, so every giant and supergiant draws the same size.
    Once the star population fix adds real giants (10-200 solar radii)
    and white dwarfs (0.01), size the Sector Map dots on a log scale so
    a giant is visibly larger than a dwarf, and check that white dwarf
    and giant systems render on the system page and System Map.

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
      (`quadrant`), `star_systems.location`, containing nebula (`_db.refresh_containment`),
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

### More pages (Boss's notes, 2026-09-30)

47. [ ] **Draw the arc-segment wireframe on the Sector Map.** Boss: "Now
    that we have defined arc segments let's add a wireframe to the
    sector map." The server still sends the cell outline
    (`starmap._outline_data`, "outline" in the scene data, from
    `sector_cell_vertices_pc`); `sectormap.js` stopped drawing it in
    16d7eed as clutter. Bring it back as the sector's real arc segment:
    the inner and outer ring faces drawn as sampled arcs rather than the
    12 straight edges between 8 corners, thin and low-contrast so the
    stars stay the focus, in both themes. Consider faint outlines of
    the neighboring cells (ring, slot and layer boundaries) too.

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
