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

- **Bug fixes (1-7)** come first, from Boss's notes of 2026-09-30. 1
   (distance units) touches the most files; 2-4 are small web changes;
   5-7 change generation constants and need the frequency research
   first.
- **Extend the cache (8)**, then do the System Map route (9). The
   local-time change (22) is small and can go in any time.
- **Galaxy Map (10-21):** Boss approved the plan in the
   project's `galaxy-megablocks/report.md` (hybrid master-wedge
   slots, pixel-sized mega-blocks). Work items 10-18 in order. 10 and 11
   ship on today's prisms. 12 is the one data-deleting step and waits for
   Boss's go-ahead on the migration. 19-20 are follow-ups. 21 (wedge
   lines) can go in any time.
- **Features (23-31)** from the same notes: generate buttons (23-24),
   phenomena views and stored nearest systems (25-26), the correlative
   update (27), navigation frames and speeds (28-29), and facilities
   (30-31). 26 and 30 are schema changes.
- Each change site in the code carries a `TODO(<area> #N)` comment
   naming its item here (areas: distances, system-list, site-header,
   search, phenomena, galaxy-map, sector-map, orbits, nav, facilities);
   grep for `TODO(` to see them all, or `TODO(galaxy-map` for one area.
- Items with a **Question for Boss** state the default taken; the work
   can start on that default.

### Bug fixes (Boss's notes, 2026-09-30)

1. [ ] **Show every distance in its most meaningful unit.** Boss: "If
   it's under AU it's km, it goes: km < AU < mpc (milliparsec) < cpc
   (centiparsec) < ly < pc < kpc < Mpc < Gpc. All radii, all distances,
   all orbital distances, etc need to use this, so all distances should
   get passed through a helper function before displayed."
   - Add one helper (Python in `stellarObjects/utils.py`, re-exported by
     `html/lib/fmt.py`; a matching `static/distance.js`) that picks the
     largest unit the value is at least 1 of. Steps: 1 AU = 1.496e8 km,
     1 mpc = 206.3 AU, 1 cpc = 2,062.6 AU, 1 ly = 30.66 cpc, 1 pc = 3.26
     ly, then factors of 1,000.
   - Today every page has its own formatter. Sites, each marked
     `TODO(distances #1)`: `tabledisplay.format_body_distance` (moons
     always km, planets AU/km/ly), `format_star_radius` and
     `utils.format_length_km` (km only), `fmt.format_distance_ly`,
     `planetData.Planet.get_table_properties`,
     `asteroidData.AsteroidBelt.to_paragraph_list`,
     `systempage._belts_table_html` and `systemmap._belt_ring_svg` (raw
     km), `system_pages._ly`/`FIELD_SPECS`, `sector_page._contents`,
     `searchpage._km`, and in JS `systemmap.js formatDistanceKm`,
     `phenomenonmap.js formatSpan`, `sectormap.js formatLy`,
     `galaxymap3d.js formatPcLy`/`formatPc`. Also `views._sector_rows`,
     `galaxy_views._quadrant_summary_rows`, `starmap._cloud_data`,
     `navmap._scale_bar_html`, `routes._sector_wiki_content`, and the ly
     in `galaxy.html` and `nav.html`.
   - **Question for Boss:** on that ladder ly only ever shows between 1
     and 3.26 ly (above that it's pc). Is that intended, or should ly
     stay alongside pc (e.g. "4.2 pc (13.7 ly)")?

   Done means no page or text output formats a distance itself, and
   tests pin one value in each unit and each boundary.

2. [ ] **Planet list: one type chip, a habitable-moon chip, belt
   distances.** Boss: "if a planet is not terrestrial it is not
   habitable so no need to display both. Likewise, no need to say both
   terrestrial and habitable, but add a new field for all planets that
   appears only if one of the moons is habitable. Asteroid belts in the
   list should display with their distance from the star (see above) in
   a meaningful unit."
   - `systempage._planet_row_html` (list) and `_body_row` (table): show
     "Gas Giant", "Terrestrial" or "Habitable", never two of them.
   - New chip (e.g. "Habitable moon") when any moon is habitable;
     `queryDb._with_life_fields`/`system_detail` sets the flag.
   - `systempage._belt_row_html`: add the belt's distance (#1).

3. [ ] **A less dense top bar.** Boss: "Admin should be a menu dropdown
   with 'Admin', 'generate' and 'logout', if the search bar text entry
   is less than twice the size of the search button, don't display it.
   Galaxy, Sectors, and such should also appear in their own separate
   menu. In fact, settings like password, theme, and admin pages should
   all be under a 'settings' gear icon in upper right,
   galaxy/sectors/systems/phenomena/nav should be in a different menu
   only if there isn't enough room to comfortably print each as a
   button, and then the above laid out search bar logic."
   - `templates/base.html` (`account_links`, the header), `style.css`
     (`.site-header`, the 56rem/92rem collapses), `helpers.SECTIONS`,
     `theme.js` (the toggle moves into the gear menu).
   - Default taken: one gear menu holds Account (password), Theme,
     Admin, Generate, Stats and Logout; there is no separate Admin
     dropdown. A container query on the header can hide the search box
     and collapse the sections without JavaScript.
   - **Questions for Boss:** should Stats stay in the gear menu (it
     wasn't in the list)? For a visitor who isn't logged in, does the
     gear hold just Theme and Login?

4. [ ] **Tag search: collapsible groups and phenomena.** Boss: "Search
   by tag should have collapsible zones for each group of tags so it
   isn't overwhelming, add stellar phenomena types to its search list.
   I can't find any nebulae on the phenomena page and they aren't
   searchable."
   - `search.html` "Browse by Tag": each group becomes a
     `<details>`, open when one of its tags is active.
   - `searchpage.TAG_FACETS`/`RESULT_PANELS` and
     `queryDb.SEARCH_TAG_FACETS`: add a phenomenon-type facet and a
     Phenomena result panel (none exists today).
   - Nebulae are missing because they're almost never generated, not
     because the page hides them: see #5.

5. [ ] **Real-world frequencies for interstellar objects.** Boss: "The
   most common interstellar objects should be asteroid field and
   comets, look up actual stats for how common each stellar object is
   and adjust probability tables adding to constants if necessary so
   it's tweakable."
   - Today (`program_constants.PHENOMENON_RATE_PER_STAR_SYSTEM`, per
     star system, Poisson per sector in `generate.generate_sector_phenomena`):
     rogue planet 0.1, comet 0.05, asteroid field 0.05, neutron star
     5e-3, black hole 5e-4, nebula 2.5e-7, supernova remnant 1e-8. With a
     handful of systems per sector, a nebula is expected about once in
     millions of sectors.
   - Figures to check against sources before using (from memory, not
     yet verified): unbound interstellar objects like 'Oumuamua run
     about 0.1-0.2 per AU³ (Do, Tucker & Tonry 2018), so "a comet" here
     means a notable one, and the rate is a design choice; rogue planets
     may outnumber stars (Sumi et al. 2023 estimate ~20 per star, mostly
     small; Mroz et al. 2017 ≤0.25 Jupiter-mass per star); neutron stars
     ~1e9 and stellar black holes ~1e8 in the Milky Way; ~8,000 H II
     regions, 10,000-25,000 planetary nebulae and ~1,000-2,000 standing
     supernova remnants.
   - Nebulae and supernova remnants are light-years across, so a
     per-volume rate (scaled by galaxy density) suits them better than a
     per-system one; decide what rate makes them findable.
   - **Question for Boss:** the real rogue-planet count is larger than
     the star count, so "most common are asteroid fields and comets"
     means capping rogues well below reality. OK to rank the rates as
     comets > asteroid fields > rogue planets > neutron stars > black
     holes > nebulae > supernova remnants, all as constants?

6. [ ] **Rogue planets: terrestrial ones, and fewer overall.** Boss: "I
   see no terrestrial planets as rogue planets, is that intentional?
   Research and revise probabilities for appearance of rogue planets
   (there are a lot of them right now) and what type."
   - Not intentional. `roguePlanetData.RoguePlanet.__init__` draws mass
     linear-uniformly over 0.0005-10 Mjup, so only ~0.5% land under the
     0.05 Mjup gas-giant threshold. Microlensing says low-mass rogues
     outnumber giants; use a mass function (log-uniform or a power law)
     with its constants next to `ROGUE_PLANET_MASS_RANGE_JUPITER`.
   - The rate itself is set in #5.
   - Update `phenomenaPlausibility` (it recomputes the expected
     terrestrial share from the same draw).

7. [ ] **A supermassive black hole in every galaxy, and rare
   intermediate ones.** Boss: "add a supermassive black hole at the
   center (or near center) of each galaxy and a smattering (rare) of
   medium sized black holes."
   - `generate.add_galactic_nucleus` only places a quasar, 10% of the
     time; the other 90% have nothing at the center. Place a quiescent
     SMBH (Sagittarius A* is ~4.3e6 Msun) when the quasar roll fails.
     Its own table or a `black_holes` row with a supermassive class
     (schema change, migration).
   - Intermediate-mass black holes already exist as 2% of black-hole
     rolls (100-1,000 Msun, `compactRemnant.BlackHole.__init__`,
     `BLACK_HOLE_INTERMEDIATE_MASS_CHANCE`). Revisit the chance and the
     range (IMBHs span ~1e2-1e5 Msun). Default taken: keep them a black
     hole subtype with its own tweakable chance.

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

### System Map (`lib/systemmap.py`, `static/systemmap.js`)

9. [ ] **Draw the Measure distance path and route it around obstacles.**
   "Measure distance" ([5.46.32]) reports a straight-line distance plus,
   when the line crosses the scene's central body, a tangent-and-arc
   detour around that one body (`computeMeasurement`,
   `routeAroundCircle`). It should also draw the path on the map, and the
   route should avoid every body it would pass through (planets, moons,
   either star of a binary), keep a safe distance from stars rather than
   just clearing the surface, and not thread between the two stars of a
   close binary. Done means the drawn path and reported distance agree
   and both respect those clearances.

### Galaxy Map and the sector standard (`static/galaxyprisms.js`, `static/galaxymap3d.js`, `lib/galaxymap3d.py`, `stellarObjects/galaxyGeometry.py`)

Today the map draws the analytic density as shrunk prisms, m sectors a
side (m a power of 3), sized by a volume budget that badly overestimates
the thin disk. The result is 70-290 px cubes with gaps, and the spiral
barely shows. The plan (report above, with renders) replaces that with a
continuous solid of mega-blocks sized from the screen's pixel scale.

10. [ ] **Make the spiral arms stand out in the expected-density shading.**
   `prismIntensity` puts log density from 0.02 to 100 on one ramp, so the
   arms (a 1.4 / 0.6 contrast at the default `arm_amplitude`) span only
   about a tenth of it. Split each block's density into its azimuthal
   mean (bulge + disk, no arms) and the arm factor (density / mean), and
   let the arm factor drive about half the ramp. Mocked in
   `galaxy-megablocks/spiral-contrast-compare.png`. Done means the arms read clearly at
   full zoom-out and at 12 kpc, in both themes, and placed-sector dots
   still stand out on top.

11. [ ] **Scale readout in sectors, pc and ly.** Replace
   `updateScaleBar`'s "≈ N pc (reference)" with three lines:
   `1 px ≈ s sectors · pc · ly`, `1 block = m sectors across (m³) · pc ·
   ly`, and a 70 px bar in the same three units. Done means the readout
   updates on every zoom and resize, and is readable at 390 px.

12. [ ] **Hybrid master-wedge slot rule (next schema version).** Boss chose it on
   2026-09-30. There are 3 master wedges at the center, doubling (6, 12,
   ..., 1,536) once each would hold at least 8 slots. Each ring's slot
   count is the multiple of its zone's master count nearest
   `2*pi*(i + 1/2)`. Sector arcs stay within ±6% of the edge (±2% today),
   the total count is unchanged, and every master line runs from where it
   starts out to the edge.
   - Change `galaxyGeometry.ring_sector_count` and its mirror in
     `galaxyprisms.js` together.
   - Add a `ring_master_count` helper.
   - `_overlapping_slots` and `neighbor_addresses` become simple integer
     ratios across aligned boundaries.
   - Update `docs/design/galaxy-coordinate-system.md`.
   - Slot counts change in all but 15 of 3,856 rings, so the next migration
     deletes galaxy-placed sectors, systems and phenomena, as v32 and v33 did.
     `galaxy_layer` and `galaxy_column` are stored by ring and layer and stay
     valid.

   **Needs Boss's explicit OK before merging**, then `update.sh` and a
   regenerate. Done means tests pin the first rings (3, 9, 15, 21, 27,
   36, ...), check that every master line is a slot boundary in every
   ring outward, and check that arcs stay in 0.94-1.06.

13. [ ] **Mega-blocks sized from the pixel scale.**
   - Replace `sectorsPerPrism` and `prismsForView` with
     `blockSizeForScale`: the smallest power of 3 with `m * edge >=
     BLOCK_MIN_PX * pcPerPixel`, with BLOCK_MIN_PX = 4.
   - Block rings and layers are the sector grid scaled by m (odd m keeps
     layer 0 on the plane).
   - Block wedges are the innermost member ring's master wedges divided by
     a power of 2. This makes every block an exact set of whole sectors, so
     `groupSectorCount` stops binning by center angle.
   - Existence depends only on ring and layer, so list exposed (surface)
     blocks only, never the whole view volume. Keep a budget guard that
     steps m up by 3.

   Done means tests check m against the scale, counts within a few percent
   of m³, surface listing against brute force, and the budget at every zoom.

14. [ ] **Continuous blocks and a slice control.**
   - Draw blocks full size (fill 1, keep the thin face edges).
   - Add a Slice control, defaulting to "cut at the focus layer", with
     "whole solid" as the alternative. A solid only shows its terraced
     outside, and a zoomed-in camera sits inside it.
   - The near cut stays for when the camera is below the cut.

   Done means the arms show at full zoom-out, zoomed views look down on a
   continuous floor, and the control works by keyboard.

15. [ ] **One solid of blocks for filled and unfilled sectors; no more
   marker dots.** Boss asked for this on 2026-09-30. The goal is to zoom in
   and out and find generated ("filled") sectors from the blocks alone.
   - **Remove the dots.** Drop the placed-sector sprites (the halo and
     core dots) and the planned-sector dots from `galaxymap3d.js`: the
     textures, `syncTier`, `withPinned`, marker scaling, and marker picking.
     Filled and unfilled sectors are both shown only through the
     continuous solid of blocks (#13, #14).
   - **Color by density.** Every block is colored by the density of the
     space it covers (#10). At m = 1, a filled sector is colored by its
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
     - A block's total sector count (`groupSectorCount`, exact with #12).

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

20. [ ] **Remove the server's leftover density sampling.** The page draws
   density itself since the prisms landed, but `queryDb.galaxy_tiles`
   still accepts `density_key` and `galaxyViewport.density_points_for_tile`
   / `density_sample_points` still exist, as does `tilecache`'s `density`
   field. Done means they're gone with their tests, and the tile cache
   still works.

21. [ ] **Wedge lines from the center.** Boss: "Galaxy map should have
   meaningful wedge lines from the center to make navigation easier."
   Draw lines in the galactic plane from the core to the edge along
   the master-wedge boundaries (#12), labelled by bearing from the
   core so they match the course format (#28). A `LineSegments`
   overlay next to the content groups in `initGalaxyMap3d`, with a
   toggle. Can ship before #12 using today's ring-0 slot lines.

### Web interface (`src/html/web/`, `src/html/static/`)

22. [ ] **Show every timestamp in the viewer's own time zone.** Boss asked
    for this on 2026-09-30. Today, times are shown in whatever zone they
    were stored or formatted in.
    - Server-rendered pages. Examples:
      - `generate_page.py`'s job `created_text`, formatted with the
        server's `time.localtime`;
      - `admin_pages.py`'s API key `created_at` and the stats page's
        "Last system change", "Newest system" and "Last sector change";
      - `errors.py`'s error-page time.
    - API responses, such as `api/auth.py`'s key `created_at`,
      `last_used_at` and `revoked_at`. These should stay UTC ISO 8601 with
      an explicit offset, since clients convert them.
    - Job status in `static/generatejobs.js`, and any time the Galaxy Map
      shows.

    Approach: the server always emits UTC as `<time datetime="...Z">` with
    a UTC fallback text, and one small script (like `theme.js`) rewrites
    each one with `Intl.DateTimeFormat` in the browser's zone
    (`Intl.DateTimeFormat().resolvedOptions().timeZone`), also showing the
    zone's abbreviation. Pages still read correctly without JavaScript
    (UTC, labelled). MySQL `DATETIME` columns carry no zone, so first
    confirm the server session's `time_zone` is UTC or convert on read.
    Done means no page shows a bare, zone-less time, and tests pin the UTC
    markup.

### Sector Map and generation (`static/sectormap.js`, `web/generate_page.py`, `generate.py`)

23. [ ] **Generate a column.** Add a `generate.py galaxy` column mode:
    every sector of one (ring, slot) column, all layers between
    `galaxy_column`'s `layer_index_min`/`max` for that ring. Today's
    modes are ring batch (one ring at one layer), one slot, local
    neighborhood and random start (`add_galaxy_arguments`, `run_galaxy`).
    **Question for Boss:** by "the entire shell", do you mean today's
    ring batch (the whole ring at this sector's layer), or the whole
    ring through every layer (a cylindrical shell)? Default taken: the
    ring at this layer, which already exists.

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
    - `generate_page.GALAXY_MODES`/`galaxy_argv`: add column (and shell
      if it's new, #23); the buttons post a CSRF-protected form that
      starts the job with the address filled in.

### Phenomena (`lib/phenomenonmap.py`, `web/system_pages.py`, `web/sector_page.py`)

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
    - **Question for Boss:** nebulae and supernova remnants weren't
      named. Keep today's map for them?

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
      generation and by #27. Schema change with a migration.
    - Show them on the sector page, system page and phenomenon page;
      `queryDb.phenomena_near_sector` returns them.

### Correlative update (`src/updateOrbits.py`, `stellarObjects/_db.py`)

27. [ ] **Check and finish the correlative update.** Boss: "We need to
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
    - **Question for Boss:** when a system's galactic orbit carries it
      out of its sector, should it move to the new sector (sector
      contents change over time), or should sectors move with the
      galaxy's rotation so membership stays fixed?

### Navigation and travel (`stellarObjects/navigation.py`, `queryDb.nav_between`, `web/nav_page.py`, `templates/nav.html`)

28. [ ] **Courses in "bearing mark mark" format on nested reference
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

29. [ ] **Warp and fold speeds.** Replace `WARP_VELOCITY_EXPONENT`
    (plain w^(10/3)) and `WARP_FACTORS_FOR_NAV` with Boss's curves, and
    add fold travel times to the nav page:
    - Warp (w): speed in c = w^(10/3) + 1 / (1 + e^(-9.3575(w - 9.5)))
      × (198.9 / (10 - w)^0.75 + 1721.7 - w^(10/3)).
    - Dimensional fold (F): speed in c = 6F⁴ / (10 - F).
    - Keep every coefficient a named constant in `program_constants`.
    - Values (1 ly per 365.25 days at 1c; 1 kpc = 3,261.56 ly):

      | Warp | Speed (c) | ly/day | Days per ly | Days per kpc |
      |---:|---:|---:|---:|---:|
      | 1 | 1.0 | 0.003 | 365.25 | 1,191,286 |
      | 2 | 10.1 | 0.028 | 36.24 | 118,191 |
      | 4 | 101.6 | 0.278 | 3.60 | 11,726 |
      | 8 | 1,024.0 | 2.804 | 0.357 | 1,163 |
      | 9 | 1,520.1 | 4.162 | 0.240 | 784 |
      | 9.5 | 1,936.0 | 5.301 | 0.189 | 615 |
      | 9.9 | 2,822.7 | 7.728 | 0.129 | 422 |
      | 9.995 | 12,201.9 | 33.41 | 0.030 | 98 |

      | Fold | Speed (c) | ly/day | Days per ly | Days per kpc |
      |---:|---:|---:|---:|---:|
      | 4 | 256.0 | 0.701 | 1.43 | 4,653 |
      | 5 | 750.0 | 2.053 | 0.487 | 1,588 |
      | 6 | 1,944.0 | 5.322 | 0.188 | 613 |
      | 6.5 | 3,060.1 | 8.378 | 0.119 | 389 |
      | 7 | 4,802.0 | 13.15 | 0.076 | 248 |
      | 7.5 | 7,593.8 | 20.79 | 0.048 | 157 |
      | 8 | 12,288.0 | 33.64 | 0.030 | 97 |
      | 8.5 | 20,880.2 | 57.17 | 0.017 | 57 |

    Done means tests pin these values, and the nav page lists warp and
    fold travel times.

### Facilities (new)

30. [ ] **Starbases, colonies and outposts in the database.** Boss: "I
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
      `utils.circular_orbital_speed_kms`), and move with #27.
    - Schema: one `facilities` table (kind, name, exactly one host,
      orbit columns), migration, API routes in `html/api/routes.py`.
    - **Question for Boss:** can moons host facilities? Default taken:
      a terrestrial moon follows the terrestrial-world rules and a
      moon can have orbital facilities.

31. [ ] **Place facilities from the web interface.** Boss: "The web
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

32. [ ] **The API can't create a system inside an existing sector.**
    `POST /api/systems` only creates standalone systems (`sector_id =
    NULL`, see `docs/api.md`). Attaching one to a sector needs the sector's
    placement and Hill-sphere separation logic (`SpaceSector.add_system`),
    which was left out of the write API to keep the admin-auth change
    small.

33. [ ] **The API can't edit a system's generated content.** `PATCH
    /api/systems/<id>` only renames. Changing stars/planets/moons/belts
    means `DELETE` then `POST` (regenerate). It may never need solving;
    kept here in case it does.

## Population and Politics

Exploratory ideas, not yet designed. Each needs a design pass before it
can be ordered against the work above.

34. [ ] Assign government ownership to star systems so that groups of
    systems form territories mapped in 3D space.
35. [ ] Flag worlds with life for generated names of their dominant
    species.
36. [ ] A database of spacefaring species.
37. [ ] Model younger and older civilizations: what differs with a
    society's age and how to store and present it.
