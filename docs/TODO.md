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

- **Extend the cache (1)**, then do the System Map route (2). The
   local-time change (14) is small and can go in any time.
- **Galaxy Map (3-13):** Boss approved the plan in the
   project's `galaxy-megablocks/report.md` (hybrid master-wedge
   slots, pixel-sized mega-blocks). Work items 3-11 in order. 3 and 4 ship
   on today's prisms. 5 is the one data-deleting step and waits for Boss's
   go-ahead on the migration. 12-13 are follow-ups. Each change site
   in the code carries a `TODO(galaxy-map #N)` comment naming its item
   here; grep for `TODO(galaxy-map` to see them all.

### Performance

1. [ ] **Add a cache so pages don't hit the database on every request.**
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

2. [ ] **Draw the Measure distance path and route it around obstacles.**
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

3. [ ] **Make the spiral arms stand out in the expected-density shading.**
   `prismIntensity` puts log density from 0.02 to 100 on one ramp, so the
   arms (a 1.4 / 0.6 contrast at the default `arm_amplitude`) span only
   about a tenth of it. Split each block's density into its azimuthal
   mean (bulge + disk, no arms) and the arm factor (density / mean), and
   let the arm factor drive about half the ramp. Mocked in
   `galaxy-megablocks/spiral-contrast-compare.png`. Done means the arms read clearly at
   full zoom-out and at 12 kpc, in both themes, and placed-sector dots
   still stand out on top.

4. [ ] **Scale readout in sectors, pc and ly.** Replace
   `updateScaleBar`'s "≈ N pc (reference)" with three lines:
   `1 px ≈ s sectors · pc · ly`, `1 block = m sectors across (m³) · pc ·
   ly`, and a 70 px bar in the same three units. Done means the readout
   updates on every zoom and resize, and is readable at 390 px.

5. [ ] **Hybrid master-wedge slot rule (next schema version).** Boss chose it on
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

6. [ ] **Mega-blocks sized from the pixel scale.**
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

7. [ ] **Continuous blocks and a slice control.**
   - Draw blocks full size (fill 1, keep the thin face edges).
   - Add a Slice control, defaulting to "cut at the focus layer", with
     "whole solid" as the alternative. A solid only shows its terraced
     outside, and a zoomed-in camera sits inside it.
   - The near cut stays for when the camera is below the cut.

   Done means the arms show at full zoom-out, zoomed views look down on a
   continuous floor, and the control works by keyboard.

8. [ ] **One solid of blocks for filled and unfilled sectors; no more
   marker dots.** Boss asked for this on 2026-09-30. The goal is to zoom in
   and out and find generated ("filled") sectors from the blocks alone.
   - **Remove the dots.** Drop the placed-sector sprites (the halo and
     core dots) and the planned-sector dots from `galaxymap3d.js`: the
     textures, `syncTier`, `withPinned`, marker scaling, and marker picking.
     Filled and unfilled sectors are both shown only through the
     continuous solid of blocks (#6, #7).
   - **Color by density.** Every block is colored by the density of the
     space it covers (#3). At m = 1, a filled sector is colored by its
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
     - Boss's wording was "decrease the opacity by a factor proportional to
       the number of filled sectors". This item reads it as "less
       see-through", so filled regions stand out. Confirm before building.
   - **Individual filled sectors appear only at sector zoom (m = 1).**
     Coarser, they show only through their block's opacity.
   - **Picking moves to blocks.**
     - Clicking a block shows its info (#9).
     - At m = 1, a filled block links to its sector page, and an unfilled
       one shows today's designation and CLI snippet.
     - Double-clicking a block with filled sectors zooms in toward them.
   - **What it needs:**
     - Per-block filled counts, counted in the browser from the tiles'
       placed lists, or totalled by the server for coarse views.
     - Translucent blocks drawn after opaque ones, sorted back to front.
     - Interior culling (#6) only where all neighbours are opaque.
     - A block's total sector count (`groupSectorCount`, exact with #5).

   Done means filled sectors can be found by zooming alone at every zoom,
   there are no marker sprites left, colors follow density, and the frame
   rate holds at the #6 block budget.

9. [ ] **Block info on click.** Clicking a block shows its sector ring,
   layer and slot ranges and its exact sector count (and how many are
   generated once tiles carry that). A click at m = 1 keeps today's sector
   panel (designation, CLI snippet, 8 corners).

10. [ ] **Smooth zooming: preload and prerender.** Today each zoom step
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

11. [ ] **Keep three.js; record why.** It was checked on 2026-09-30:
   - Babylon.js is several MB, and deck.gl needs a bundler.
   - regl and raw WebGPU would mean rewriting picking, sprites and
     lighting by hand.
   - The CSP (`default-src 'self'`) and the no-build-step vendoring rule
     favor one vendored file.
   - The bottleneck is JavaScript listing work, not the renderer.

   three r186 already has InstancedMesh, BatchedMesh and a
   WebGPURenderer to move to later. Done means the rendering choice is
   written into `docs/html-interface.md`.

12. [ ] **Follow-ups (edge cases).**
   - Distance-based detail (bigger blocks farther from the camera), which
     the aligned wedges from 5 make seamless.
   - Order-independent transparency (weighted blended) if #8's sorting
     shows artifacts where translucent blocks intersect.
   - Phone performance at 390 px.
   - DPR: `pcPerPixel` is per CSS pixel.
   - Reduced-motion users get instant zoom.

13. [ ] **Remove the server's leftover density sampling.** The page draws
   density itself since the prisms landed, but `queryDb.galaxy_tiles`
   still accepts `density_key` and `galaxyViewport.density_points_for_tile`
   / `density_sample_points` still exist, as does `tilecache`'s `density`
   field. Done means they're gone with their tests, and the tile cache
   still works.

### Web interface (`src/html/web/`, `src/html/static/`)

14. [ ] **Show every timestamp in the viewer's own time zone.** Boss asked
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

### Web API (`src/html/api/routes.py`)

Low priority; nobody is waiting on these.

15. [ ] **The API can't create a system inside an existing sector.**
    `POST /api/systems` only creates standalone systems (`sector_id =
    NULL`, see `docs/api.md`). Attaching one to a sector needs the sector's
    placement and Hill-sphere separation logic (`SpaceSector.add_system`),
    which was left out of the write API to keep the admin-auth change
    small.

16. [ ] **The API can't edit a system's generated content.** `PATCH
    /api/systems/<id>` only renames. Changing stars/planets/moons/belts
    means `DELETE` then `POST` (regenerate). It may never need solving;
    kept here in case it does.

## Population and Politics

Exploratory ideas, not yet designed. Each needs a design pass before it
can be ordered against the work above.

17. [ ] Assign government ownership to star systems so that groups of
    systems form territories mapped in 3D space.
18. [ ] Flag worlds with life for generated names of their dominant
    species.
19. [ ] A database of spacefaring species.
20. [ ] Model younger and older civilizations: what differs with a
    society's age and how to store and present it.
