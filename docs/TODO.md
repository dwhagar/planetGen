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

- **The known generation bugs (40-45) come first**; their tests are
   already written.
- **Bug fixes (2-7)** come next, from Boss's notes of 2026-09-30. 2-4
   are small web changes;
   5-7 change generation constants; the frequency research for 5 and 6
   is in `docs/design/interstellar-object-rates.md`.
- **Extend the cache (8)**, then do the System Map route (9). The
   local-time change (22) is small and can go in any time.
- **Galaxy Map (12-19):** Boss approved the plan in the
   project's `galaxy-megablocks/report.md` (hybrid master-wedge
   slots, pixel-sized mega-blocks). Work items 12-18 in order (14, the solid and its slice, shipped early at Boss's request). 12 is the
   one data-deleting step. 19 is follow-ups.
- **Features (23-36)** from the same notes: generate buttons (23-24),
   phenomena views and stored nearest systems (25-26), nebulae and
   remnants: placement, classes, containment and naming, plus asteroid
   field classes (27-31), the correlative update (32), navigation frames
   and speeds (33-34), and facilities (35-36). 26, 27, 28, 29, 30 and 35
   are schema changes.
- **More pages (46-49)**: the full systems list, the Sector Map
   wireframe, the system page layout and non-overlapping System Map
   names are small and can go in any time.
- **Installers (50)**: PowerShell install and upgrade scripts, and the
   bash scripts made to run on macOS too.
- Each change site in the code carries a `TODO(<area> #N)` comment
   naming its item here (areas: distances, system-list, site-header,
   search, phenomena, galaxy-map, sector-map, orbits, nav, facilities,
   security, physics, web-pages, installers);
   grep for `TODO(` to see them all, or `TODO(galaxy-map` for one area.
- Items with a **Question for Boss** state the default taken; the work
   can start on that default.

### Bug fixes (Boss's notes, 2026-09-30)

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
   - `systempage._belt_row_html`: add the belt's distance (`fmt.format_distance_km`).

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
   - Boss confirmed on 2026-09-30: a visitor who isn't logged in gets
     the pages (Galaxy, Sectors, Systems, Phenomena, Nav) and, under the
     gear, Theme and search. No Stats unless logged in. ("Admion" was a
     typo for Admin.)
   - Default taken: for a logged-in admin the gear menu holds Account
     (password), Theme, Admin, Generate, Stats and Logout; there is no
     separate Admin dropdown. A container query on the header can hide
     the search box and collapse the sections without JavaScript.
   - **Question for Boss:** a visitor still needs a way to log in; a
     Login entry in the gear menu is the default.

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
   - Boss supplied the research on 2026-09-30 and asked to reorganize the
     rates from it. Numbers, scaling functions, checks, measured cost and
     references: `docs/design/interstellar-object-rates.md`.
   - Today (`program_constants.PHENOMENON_RATE_PER_STAR_SYSTEM`, per
     star system, Poisson per sector in `generate.generate_sector_phenomena`):
     rogue planet 0.1, comet 0.05, asteroid field 0.05, neutron star
     5e-3, black hole 5e-4, nebula 2.5e-7, supernova remnant 1e-8.
   - The model: each type has a local density `n_i` (per pc³) at the
     solar neighborhood's `n_*0 = 0.14` stars per pc³, scaled by the
     sector's own stellar density: `n_i = n_i0 * n_* / 0.14`. That is a
     fixed rate per star, so the existing per-sector Poisson draw stays;
     only the rates change.

     | Type | n_i0 (pc⁻³) | Per star | Scales with |
     |---|---|---|---|
     | Comets and debris | 1e12 (1e11-1e14) | ~7e12 | n_* |
     | Terrestrial rogue planets | 0.7 (0.5-1.4) | 5 (2-10) | n_* |
     | Jupiter-mass rogue planets | 0.035 | ≤ 0.25 | n_* |
     | Rogue brown dwarfs | 0.03 | ~0.21 | n_* |
     | Runaway stars | 2.1e-3 | 0.015 | n_* |
     | Isolated neutron stars | 1e-3 | ~7e-3 | n_* |
     | Isolated black holes | 1e-4 | ~7e-4 | n_* |
     | Giant molecular clouds | 5e-6 | n/a | ρ_gas^1.4, filling 1-2% of arms |
     | Planetary nebulae | 3e-8 | ~2e-7 | n_* |
     | Supernova remnants | 1e-8 | n/a | n_* · ρ_gas |
     | Hypervelocity stars | 1e-10 at 8 kpc (research says 5e-9) | n/a | r_GC⁻² |
     | Isolated asteroid fields | ~0 | 0 | none (they disperse) |

   - Constants to add in `program_constants` (replacing
     `PHENOMENON_RATE_PER_STAR_SYSTEM`), each with its source in a
     comment: `REFERENCE_STELLAR_DENSITY_PC3 = 0.14`;
     `PHENOMENON_DENSITY_PC3` keyed by type (`"rogue-planet"` is the sum
     of the item 6 bins, `"brown-dwarf"`, `"runaway-star"`,
     `"hypervelocity-star"`, `"neutron-star"`, `"black-hole"`,
     `"molecular-cloud"`, `"planetary-nebula"`, `"supernova-remnant"`,
     `"comet"`, `"asteroid-field"`); `GMC_ARM_FILLING_FACTOR = 0.015`;
     `GMC_GAS_DENSITY_EXPONENT = 1.4`; `HVS_REFERENCE_RADIUS_PC = 8000`;
     `INTERSTELLAR_DEBRIS_DENSITY_PC3 = 1e12`; and one
     `PHENOMENON_RATE_SCALE` per type (default 1.0) to dial any type up
     or down without touching the research value.
   - Generation cost (measured, see the design doc): at these rates a
     local-density 4 pc sector gets about 50-60 rogue rows on top of
     roughly 650 body rows, about 8-10% more rows and under 1% more
     generation time. Cost is not the problem; the sector's look is.
     Rogues would outnumber systems about five to one in the Contents
     table and on the Sector Map.
   - Recommended way to apply them:
     - Rogue planets, brown dwarfs, neutron stars, black holes: rows,
       per-star rate times the sector's star count, Poisson, no cap.
     - Comets: the real density (~6e13 per sector) can't be rows. Show
       it as a computed sector figure ("about 10^13 interstellar comets
       and planetesimals"), and keep the `comet` rows as notable comets
       at a design rate (default 0.007 per pc³, which is today's 0.05
       per star); `INTERSTELLAR_DEBRIS_DENSITY_PC3` feeds the figure.
     - Isolated asteroid fields: rate 0 (Boss's research: a free
       asteroid field disperses in 1e6-1e7 years). Keep the type and
       table for hand-made fields and facilities (#31, #35); none get
       generated.
     - Runaway and hypervelocity stars: a flag and speed on an ordinary
       generated system (1.5% of systems; HVS by `r_GC⁻²`), not a new
       phenomenon table.
     - Brown dwarfs: a new rogue-body kind (a rogue planet row above
       13 Mjup, or its own table); schema change.
     - Nebulae and remnants: real point rates make them essentially
       never appear (a planetary nebula ~2e-6 per sector). Molecular
       clouds place by filling factor (about 1-2% of arm sectors sit
       inside one), with #27 generating the stars inside; planetary
       nebulae and remnants keep their real rates.
   - The research's comet separation column, hypervelocity density and
     remnant density don't add up; see "Checks" in the design doc.
   - **Question for Boss:** show rogue planets at the full research rate
     (about 45 per sector, five to one against systems), or scale them
     down for readability with `PHENOMENON_RATE_SCALE["rogue-planet"]`
     (0.1 gives about 5 per sector)? Default: full rate, with rogues
     grouped into one collapsible row in the sector's Contents table.
   - Regenerate the galaxy after this ships to see the new mix.

6. [ ] **Rogue planets: terrestrial ones, and fewer giants.** Boss: "I
   see no terrestrial planets as rogue planets, is that intentional?
   Research and revise probabilities for appearance of rogue planets
   (there are a lot of them right now) and what type."
   - Not intentional. `roguePlanetData.RoguePlanet.__init__` draws mass
     linear-uniformly over 0.0005-10 Mjup, so only ~0.5% land under the
     0.05 Mjup gas-giant threshold. Microlensing says low-mass rogues
     outnumber giants several to one; terrestrial rogues must stay
     present.
   - Boss's research (2026-09-30): terrestrial rogues 2-10 per star
     (Johnson et al. 2020; OGLE-2016-BLG-1928, 0.3-2 M⊕, Mróz et al.
     2020); Jupiter-mass rogues at most 0.25 per star (Mróz et al. 2017,
     superseding Sumi et al. 2011's 1.8). Sub-Neptunes are expected to
     dominate with terrestrials (Barclay et al. 2023) but have no rate.
   - Draw a mass bin by its per-star rate, then a log-uniform mass inside
     the bin. Constant to add next to `ROGUE_PLANET_MASS_RANGE_JUPITER`
     (which it replaces): `ROGUE_PLANET_MASS_BINS`, each bin
     `(min_mass_earth, max_mass_earth, per_star_rate)`:
     - terrestrial 0.1-2 M⊕, 5 per star;
     - sub-Neptune 2-20 M⊕, 1 per star (default, not from the research);
     - Saturn-class 20 M⊕-1 Mjup, 0.25 per star (default, not from the
       research);
     - Jupiter-mass 1-13 Mjup, 0.25 per star.
     Above 13 Mjup is a brown dwarf (#5). A single power law can't fit
     both ends (the design doc shows why).
   - The bins' summed per-star rate is the rogue-planet rate in #5.
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
    - Rates come from #5 (a per-volume rate suits objects this big).
    - Each class brings its central object (table in
      `docs/design/nebula-and-asteroid-field-classes.md`): O/B stars for
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

28. [ ] **Class nebulae and remnants A-Z by what's in them.** Boss:
    "We also want Nebulae and Stellar remnants to be classed like
    planets ... and we'll need to come up with What's IN the Nebulae.
    ... develop a letter-class system similar to planets (A to Z) based
    on contents of the nebulae."
    - Draft classes A-W (I and O unused, X-Z reserved) are in
      `docs/design/nebula-and-asteroid-field-classes.md`, built from
      Boss's reference document: diffuse (A-B), H II (C-E), reflection
      (F-G), planetary (H-L), molecular (M-Q), supernova remnants (R-W).
    - Add `NEBULA_CLASSES` next to `PLANET_CLASSES` in
      `program_constants` (description, contents, radius, nH,
      temperature, extinction, central-object rule, frequency),
      replacing `NEBULA_TYPES` and `SUPERNOVA_REMNANT_MORPHOLOGIES`.
    - Store the class and contents (dominant species, density,
      temperature, extinction) on `nebulae` and `supernova_remnants`;
      show them on the phenomenon page and as a search facet (#4).
    - **Question for Boss:** does "stellar remnants" also mean the
      compact objects (white dwarfs, neutron stars, black holes)? The
      draft classes supernova remnants only.

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
      (class from #31). Examples and the size rule are in the design
      doc. Needs a migration that renames existing rows.

31. [ ] **Class asteroid fields A-Z.** Boss: "Asteroid fields should also
    have classes (A to Z) based on composition and density and size."
    Draft in the design doc: the letter comes from composition
    (carbonaceous, stony, metallic, icy, basaltic, mixed, dust,
    collisional family) and density (today's sparse/typical/dense), and
    size is the digit in the designation (#30). `asteroidFieldData`,
    `ASTEROID_FIELD_*` constants, `asteroid_fields` table, phenomenon
    page. **Question for Boss:** should size be part of the letter
    instead?

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

### Known generation bugs (strict xfail tests)

Each has a test marked `xfail(strict=True)` that starts passing, and so
fails the run, once the bug is fixed; remove the marker in the same PR.
Each has a `TODO(physics #N)` comment where the fix goes.

40. [ ] **Moons orbit outside their planet's Hill sphere.**
    `planetPhysics.generate_moons` sets `high_orbit` to 5 Hill radii
    (the comment says 1/5). Test:
    `test_fuzz_system_generation.py::test_moons_orbit_inside_their_parents_hill_sphere`.
41. [ ] **Moons can orbit inside their planet.** `generate_moons`'
    `low_orbit` ignores the planet's radius. Test:
    `test_moons_orbit_outside_their_parents_body`.
42. [ ] **A close binary's planets can orbit inside the binary.**
    `StarSystem._generate_planets` has no floor at the stars' separation
    (about 5% of close binaries). Test:
    `test_circumbinary_bodies_orbit_outside_the_binary`.
43. [ ] **A planet's Hill sphere can overlap the belt inside it.**
    `StarSystem.validate_system` keeps a planet only 0.05 AU past a belt,
    but a belt after a planet must clear 5 Hill radii. Test:
    `test_planet_hill_sphere_clears_the_belt_inside_it`.

### More pages (Boss's notes, 2026-09-30)

46. [ ] **A paginated list of every system on the Systems page.** Boss:
    "the systems page should have a paginated list of all systems."
    Today `/systems` (`views.systems`, `_systems_panel`) lists only
    standalone systems (`sector_id="none"`). List every system, 50 rows
    a page through the shared pager (`html/lib/pagination.py`), with
    its sector and octant; keep the standalone list as its own panel or
    a filter. `apiclient.get_systems` without `sector_id` already pages
    all systems.

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

48. [ ] **System page: one ordered list of everything in orbit, with
    expandable moons, stars table first.** Boss: "on the star system
    page we list planets and moons twice, have moons expandable under
    the initial planet list so we can get rid of the 2nd table at the
    end of the page. Move the star table to the top of the tables under
    the clickable map interface."
    - Moons: in `systempage._planet_row_html` a planet's moons sit inside
      its `<details>`, after its whole description. Give them their own
      expandable group right under the planet's row (e.g. a nested
      "N moons" `<details>`), still working without script.
    - Remove the Planets & Moons table (`bodies_html` /
      `_planets_table_html`, rendered last in `system.html`). Its columns
      (class, type, zone, distance, period, gravity) move into each
      row's compact stats so nothing is lost; item 2's chip rules apply
      there.
    - Move `stars_html` (the stars table) up to be the first table,
      directly under `map_html`.
    - Belts and comets (Boss, 2026-09-30): "Asteroid belts and comets
      should be placed in the interactive list of objects in orbit
      around the star in their relative order from the star." Remove
      the Asteroid Belts and Comets tables too. `_orbiting_rows_html`
      already orders belts with planets by `orbital_index`, but appends
      every comet at the end; sort comets in by distance instead.
      Default taken: a comet sorts by its semi-major axis
      (`perihelion_distance_km / (1 - eccentricity)`), and parabolic
      ones (no finite axis) go last by perihelion. Each belt and comet
      row shows its distance (`fmt.format_distance_km`).

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
