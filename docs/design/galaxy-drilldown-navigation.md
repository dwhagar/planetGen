# Galaxy Map drill-down navigation

Boss's design for getting around the Galaxy Map, recorded 2026-10-01.

**Status (2026-10-01, checked against 7.58.2):** mostly built.

| Piece | Section | TODO | Built in |
|---|---|---|---|
| Block ladder (`stellarObjects/galaxyDrill.py`, `drill*` in `static/galaxyprisms.js`) | 3 | MAP.28 | 7.41.2, PR #160 |
| Stage contents API (`GET /api/galaxy/stage`, the site's cached `/galaxy/stage`) | 7 | MAP.29 | 7.41.3, PR #160 |
| The stages (`static/galaxystages.js`, `static/galaxystageview.js`), stage URLs, breadcrumb, keys, touch, "Generated only"; the old camera behind Free look | 4, 5, 8.1 | MAP.16 | 7.44.0, PR #171 |
| Address bar (`/galaxy/locate`) | 9.3 | MAP.24 | 7.50.0, PR #172 |
| Course on the map (`/galaxy?course=<from>,<to>`) | 9.4 | MAP.27 | 7.52.0, PR #176 |
| Map Generate buttons and the light-year radius dialog | 6 | MAP.20 (map side) | 7.53.0, PR #177 |
| Sector Map pick mode and Nav from/to links | 9.1, 9.2 | MAP.21 | 7.58.0, PR #178 |
| Generate this layer or slab (`generate.py galaxy --block`) | 6 | MAP.20 | not built |
| NAV page pickers | 9 | MAP.22 | not built |
| Bookmarks | 8.2 | MAP.23 | not built |
| "Show on Galaxy Map" links with `?sector=` from sector, system and search pages | 8.1 | MAP.25 | not built (`?sector=` itself works) |

Two later requests change this design: MAP.17 (no free camera and no
drag-rotate at all; start top-down and pick a wedge, then a slice, then a
block), which settles decision 2 against Free look, and MAP.26 ("Show on
Galaxy Map" opens at the sector level, and the map gets its own Back and
Forward). Neither is built. The open work items are in `docs/TODO.md`,
and each one points back to a section here.

## 1. What Boss asked for

Boss's words, from the project thread on 2026-10-01:

> What if we select in sections? As in we click on a mega-block and it
> zooms in to just that megablock's contents to have smaller megablocks
> fill that space. From there, the mouse can select a layer in megablocks
> that then is pulled out into a slice that the user can select another
> megablock inside and it'll zoom it in to that contents, and repeat
> until we have sectors and clicking on one goes to that sector. In fact,
> do the slice of the 3D to zoom in and select a block be the first
> stage. Mouse will highlight the slice it's going to make and from the
> slice highlight the block it's going to zoom into. Then we can add
> bookmarks and integrate this with the navigation system, so we can
> from the nav menu select start and destination using either the text
> dropdowns as we have now or the galactic map interface to select. If
> it's within sector then it'll just use the sector interface.

> Also let's add UI so that once we're down to a sector level we can tell
> a slice to generate all the sectors in that slice or click on a sector
> and generate it from the UI if you're admin

> we'll use the bigger targets, also add a 'generate neighborhood' when
> at a sector selection level that will ask the radius in ly.

"Bigger targets" picks the 243 → 27 → 3 → 1 ladder (section 3) over
81 → 9 → 1.

## 2. Background: the research paper

Boss supplied a paper, "Mapping Volumetric Galactic Topologies to Planar
Interfaces" (an HCI report on cylindrical galaxies). Its main claims, and
how they apply here:

| Paper idea | Used here? |
| --- | --- |
| Discrete scale tiers instead of continuous scroll zoom | Yes: the drill-down stages are the tiers, picked by clicking instead of by the wheel |
| Breadcrumb with clickable levels; an address bar for direct entry | Yes (sections 5.4, 9.3) |
| Animated scale-space flight (van Wijk and Nuij, 2003) | Yes, for every move into a block (section 5.3) |
| Top-down polar plan view; elevation (r, z) slice | Yes: every slab is shown top-down, and the slab strip picks the height (section 5.2) |
| Bookmarks and history | Yes (sections 5.4, 8) |
| A 2D canvas replacing 3D | No. The 3D solid of blocks stays for the slab-picking stages. |
| Unwrapping each cylindrical shell into a strip | No. A ring is one sector wide, so its strip is a thin ribbon with a seam, and it cuts off the radial neighbors. |
| Morton (Z-order) hashing | No. Every address here is already found by closed-form math, and the designation hex is already a reversible packed key. |
| Fisheye lenses | No. They would bend the distances and bearings that NAV reports. |

## 3. The block ladder

### 3.1 Notation

- `e`: the sector edge, `edge_pc` (4 pc by default).
- A sector is `(i, j, k)`: ring `i`, layer `j`, slot `k`. Ring `i` has
  `N_i = ring_sector_count(i)` slots, and slot `k` spans
  `theta in [2*pi*k/N_i, 2*pi*(k+1)/N_i)`.
- A block is `(m, I, s, S)`: `m` sectors a side (a power of 3), block
  ring `I`, wedge `s`, slab `S`. Block ring `I` has `W_m(I)` wedges.
- Bounds of block `(m, I, s, S)`:
  - `r in [I*m*e, (I+1)*m*e)`
  - `theta in [2*pi*s/W_m(I), 2*pi*(s+1)/W_m(I))`
  - `z in [(S - 1/2)*m*e, (S + 1/2)*m*e)`
- Member rings: `I*m` to `I*m + m - 1`. Member layers: `S*m - (m-1)/2`
  to `S*m + (m-1)/2`, so slab 0 is centered on the plane, the same way
  layer 0 is. These are `blockSectorRanges` in `static/galaxyprisms.js`.

### 3.2 The ladder: 243, 27, 3, 1

Each step into a block divides its size by 9, except the last, which
divides by 3. Numbers below are for the default galaxy (edge ring 3,855,
635 layers from -317 to 317, about 12.3 billion sector cells). They were
computed from the code on 2026-10-01.

| Level | Block size | Block rings | Slabs | Wedges per ring |
| --- | --- | --- | --- | --- |
| 243 | 972 pc (3,170 ly) | 16 | 3 (-1 to 1) | 3 at the core, 96 at the edge |
| 27 | 108 pc (352 ly) | 143 | 25 (-12 to 12) | 3 at the core, 1,024 at the edge |
| 3 | 12 pc (39 ly) | 1,286 | 213 (-106 to 106) | 3 at the core, 7,680 at the edge |
| 1 (sectors) | 4 pc (13 ly) | 3,856 | 635 layers | `N_i` |

Slab 0 at level 243 holds layers -121 to 121 (within 486 pc of the
plane), which is nearly every star. Slabs 1 and -1 hold the thick disk
and bulge out to layer 317.

### 3.3 Rings and layers nest exactly

With `f = m / m'` (9 or 3) the child size `m'`:

- Rings: block ring `I` at size `m` covers sector rings `I*m` to
  `I*m + m - 1`. These are exactly child block rings `I*f` to
  `I*f + f - 1`, so the parent ring of child ring `I'` is `floor(I'/f)`.
- Layers: the children of slab `S` are child slabs `S*f - (f-1)/2` to
  `S*f + (f-1)/2` (9 or 3 of them), and the parent slab of child slab
  `S'` is `round(S'/f)`. Proof: the lowest child slab's lowest layer is
  `(S*f - (f-1)/2)*m' - (m'-1)/2 = S*m - (m-1)/2`, the parent's lowest
  layer. The top layer follows by symmetry. `f` and `m'` are odd, so the
  halves are whole numbers.

### 3.4 Wedges must be made to nest

Today's `blockWedgeCount(ring, m)` picks each level's wedge count on its
own (aligned to the master lines where it can, within 1.2x of
`round(2*pi*(I + 1/2))`). Child wedges then only sometimes sit inside one
parent wedge. Counted on the default galaxy, a child ring's wedges divide
evenly into its parent's in 109 of 143 rings at 243 → 27, 1,162 of 1,286
at 27 → 3, and 1,468 of 3,856 at 3 → 1. A drill-down needs each child to
sit wholly inside one parent, so the ladder uses a **nested wedge rule**:

```
T(i)    = max(3, round(2*pi*(i + 1/2)))            ring i's target count
W_243(I) = blockWedgeCount(I, 243)                   unchanged
W_27(i)  = W_243(p) * max(1, round(T(i) / W_243(p))),  p = floor(i / 9)
W_3(i)   = W_27(p)  * max(1, round(T(i) / W_27(p))),   p = floor(i / 9)
```

Every child count is a whole multiple of its parent's, so every parent
wedge line is also a child wedge line. With `q = W_child(i) / W_parent(p)`:

- Children of wedge `s` in child ring `i`: wedges `s*q` to `s*q + q - 1`.
- Parent of child wedge `s'`: `floor(s' / q)`.

The cost is a little arc error: a wedge's centerline arc stays within
0.93-1.07 of a block edge at level 27, and within 0.95-1.06 at level 3.
Level 243 keeps today's rule (0.85-1.18). Each parent holds 69-96
children per child slab at 243 → 27 and 74-90 at 27 → 3, times 9 child
slabs.

The new rule is used by the drill-down stages. If the free camera stays
(decision 2), its pixel-sized blocks should switch to the same rule, so
that both views draw the same block borders.

### 3.5 Sectors inside a level-3 block

Sector slots don't line up with level-3 wedges (only 1 of 1,286 rings
does), so sectors join the level-3 block that holds their center, the
rule `blockSlotRange` already uses. In ring `i` of block
`(3, J, s, S)` with `W = W_3(J)` and `N = N_i`:

```
first slot = ceil((2*s*N - W) / (2*W))
last slot  = ceil((2*(s+1)*N - W) / (2*W)) - 1
```

The members are rings `3J` to `3J + 2` and layers `3S - 1` to `3S + 1`.
That is 1 to 5 slots per ring (5 only near the core) and about 9
sectors per layer (8.2 to 10 on average).

The parent of sector `(i, j, k)`, using integer math only:

```
J = floor(i / 3)
S = round(j / 3)
s = floor((2*k + 1) * W_3(J) / (2 * N_i))
```

The rest of the chain follows from section 3.4. Each sector then has
exactly one block at every level, and a block's sectors are exactly the
sectors of its children.

### 3.6 What exists and what counts

- A block is drawn when any member sector is allowed by the stored
  outline (`galaxy_layer` and `galaxy_column`). In the browser,
  `blockExists`/`sectorAllowed` already compute this from the analytic
  bound.
- A slab is drawn when any of its blocks is.
- A block's total is its allowed sector count (`blockSectorCount`,
  adapted to the nested rule). Its generated count comes from the
  server (section 7). The look stays as today: glass when nothing is
  generated, more solid as the generated share grows, fully solid at
  100%.

### 3.7 Code

- `static/galaxyprisms.js`: `drillWedgeCount(level, ring)`,
  `drillParent(block)`, `drillChildren(block)` (grouped by child slab),
  `drillSlabs(block)`, `drillBlockSectors(block, layer)` and
  `drillChainOf(ring, layer, slot)`. Pure, and runs under node like the
  rest of the file.
- `stellarObjects/galaxyDrill.py`: the same functions in Python, for
  the API (section 7) and generation (section 6), plus
  `format_drill_key`/`parse_drill_key` for the `m.ring.wedge.slab` keys
  (`formatDrillKey`/`parseDrillKey` on the page). Built 2026-10-01
  (PR #160). Generation (section 6) does not use it yet.
- `tests/test_galaxydrill.py` runs both on every block ring and on
  sampled sectors and checks that they agree, and that a parent's
  sectors are exactly the sectors of its children.

## 4. The stages

Eight clicks reach a sector: slab, block, slab, block, slab, block,
slab, sector.

| Stage | View | Holds (default galaxy) | Click does |
| --- | --- | --- | --- |
| 1 Galaxy | 3D | Level-243 blocks, 3 slabs | Pull out a slab |
| 2 Galaxy slab | Top-down | One level-243 slab, up to 818 blocks | Fly into a block |
| 3 Block (972 pc) | 3D | 9 slabs of 69-96 level-27 blocks | Pull out a slab |
| 4 Block slab | Top-down | 69-96 level-27 blocks | Fly into a block |
| 5 Block (108 pc) | 3D | 9 slabs of 74-90 level-3 blocks | Pull out a slab |
| 6 Block slab | Top-down | 74-90 level-3 blocks | Fly into a block |
| 7 Block (12 pc) | 3D | 3 layers of about 9 sectors | Pull out a layer |
| 8 Sector layer | Top-down | About 9 sectors (up to 15) | Open the sector (section 5.5) |

Stages 7 and 8 are the sector selection level: single sectors are what
you see and pick. The admin Generate tools (section 6) appear there.

## 5. How the stages behave

### 5.1 3D stages (1, 3, 5, 7): pick a slab

- Only the current container is drawn: the galaxy, or one block's
  children. Every block of every slab is drawn (a few hundred to about
  900), with no surface-only listing, so the 60,000 budget never
  matters.
- **Hover:** a raycast finds the block under the pointer. Its whole slab
  is highlighted with an accent outline and full opacity, and the other
  slabs drop to about 35% opacity. The slab strip (below) highlights the
  same row, and a tooltip names the slab, for example "Slab -7: layers
  -22 to -20, 78-90 pc below the plane, 12 of 2,187 sectors generated".
- **Click:** the slab is pulled out. It rises by 0.6 of its own
  thickness while the others fade to 0, and the camera's polar angle
  goes from where it is to 0 (straight down). Both run over 600 ms with
  ease-in-out. With `prefers-reduced-motion` it cuts straight to the top
  view.
- **Drag** still rotates the container. The wheel zooms only between
  fitting the container and 1.5 times closer. Neither changes the stage.
- **Slab strip:** a vertical strip beside the map with one row per slab,
  top slab first. Each row shows the slab's layer range and a bar of its
  generated share. Hovering a row highlights the slab, and clicking it
  pulls it out. This is the paper's elevation panel, used as the picker,
  and it gives keyboard and screen-reader users the same choice.

### 5.2 Top-down stages (2, 4, 6, 8): pick a block

- **Camera:** looks straight down (-z). The view fits the slab's
  footprint: camera height `h = 1.1 * R_fit / tan(FOV/2)`, where
  `R_fit` is the largest distance from the footprint's center to its
  corners, with the narrower canvas side used for the FOV. Stage 2 keeps
  bearing 000 (+x) pointing right. Stages 4, 6 and 8 turn the view so
  the block's outward radial direction points up (core at the bottom),
  and a small compass arrow shows which way bearing 000 lies
  (decision 6).
- Wedge lines of the parent and grandparent levels stay drawn faintly as
  context (the paper's "multiscale residue").
- **Hover:** outlines the block under the pointer and shows a tooltip:
  address, bearing range, distance from the core, and generated over
  total sectors.
- **Click:** flies into the block (section 5.3). Its children fade in
  over the last 200 ms, and the next stage is a 3D view of them.
- Stage 2 is up to 818 blocks across a 32-block-wide disk, so the wheel
  zooms (up to 4 times) and drag pans inside it. Stages 4 to 8 fit on
  the canvas and don't zoom.

### 5.3 The flight into a block

The camera flies with van Wijk and Nuij's optimal zoom-and-pan path, in
the plane of the slab. Here `w` is the visible width,
`w = 2 * h * tan(FOV/2)`, and `u` is the distance from the start center
`c0` to the end center `c1`, with `rho = 1.4`.

```
b0 = (w1^2 - w0^2 + rho^4 * u^2) / (2 * w0 * rho^2 * u)
b1 = (w1^2 - w0^2 - rho^4 * u^2) / (2 * w1 * rho^2 * u)
r0 = ln(-b0 + sqrt(b0^2 + 1))
r1 = ln(-b1 + sqrt(b1^2 + 1))
S  = (r1 - r0) / rho                           path length

w(s) = w0 * cosh(r0) / cosh(rho*s + r0)
d(s) = (w0 / rho^2) * (cosh(r0) * tanh(rho*s + r0) - sinh(r0))
center(s) = c0 + (c1 - c0) * d(s) / u
```

When `u` is close to 0, `S = |ln(w1/w0)| / rho` and
`w(s) = w0 * exp(±rho*s)`. The flight takes `T = clamp(S / 1.1, 0.5, 1.6)`
seconds, with `s` advanced in proportion to time. During a flight, the
camera's tilt also eases from 0 to the 3D stage's default 35 degrees
from vertical. With `prefers-reduced-motion` the camera cuts straight to
the end. The same path is used for the address bar's jumps
(section 9.3) and for Back.

### 5.4 Breadcrumb, Back and keys

- A breadcrumb above the map shows the chain, for example:
  `Galaxy › Slab 0 › Block 7·14 › Slab -1 › Block 63·115 › Slab -7 ›
  Block 568·1036 › Layer -20 › Sector 1,705·-20·3,225` (a real chain
  under the nested rule). Blocks read
  `ring·wedge`, and the last crumb is the address `ring·layer·slot`.
  Each crumb's tooltip gives its bearing range and distance from the
  core. Clicking a crumb returns to that stage.
- Each crumb has a ▾ menu of its siblings (the other slabs of that
  container, or the neighboring blocks), so a sideways move needs no
  climb.
- Keys: ← → move to the previous or next wedge (block) at the same
  stage, ↑ ↓ the next slab up or down (or the next ring out or in, on
  top-down stages), Enter acts on the highlighted item, Esc or Backspace
  goes up one stage, and Home returns to the galaxy.
- **Touch:** the first tap highlights a slab or block, and a second tap
  on it confirms. Pinch zooms stage 2.
- A **"Generated only"** toggle dims and disables blocks with nothing
  generated, so existing content is easy to follow.

### 5.5 The sector at the end

At stage 8 (or by clicking a sector in stage 7's 3D view):

- A generated sector opens its sector page (the Sector Map).
- A sector that isn't generated shows today's cell panel (address,
  designation, coordinates, 8 corners). Admins also get the Generate
  tools (section 6).

## 6. Generating from the map (admin)

Only an admin with current credentials sees these, as with PR #151 (the
server puts the `generate` target in the scene data). They appear at
stages 7 and 8.

| Button | What runs | Size |
| --- | --- | --- |
| Generate this sector | `generate.py galaxy --ring i --layer j --slot k` (today's `slot` mode) | 1 |
| Generate this layer (stage 8) or slab (stage 7, the highlighted one) | New `--block 3.J.s.S --block-layer j`: the block's ungenerated allowed sectors on that layer (section 3.5) | Up to about 15 |
| Generate neighborhood… | `slot` mode with `--radius-pc`, after the radius dialog below | Grows with the radius cubed |

The `--block m.I.s.S` mode (with optional `--block-layer`) also accepts
larger blocks, so a "Generate this block" button for stage 5 or 3 is a
later option (decision 3). It obeys `generationLimits.py` like every
other mode.

**Neighborhood dialog:** it opens when the button is pressed and is
centered on the selected sector, generated or not.

- One number field, "Radius (light-years)". The default is 100 ly
  (today's 30.7 pc). The minimum is 13 ly (one sector edge). The maximum
  is `MAX_GENERATE_RADIUS_PC` converted to ly: 200 pc, about 652 ly.
- Conversion: `radius_pc = radius_ly / 3.26156`, rounded to 0.1 pc and
  sent as `slot_radius_pc` (the same `ly_to_pc` as
  `stellarObjects.utils`).
- A live estimate under the field: "about N sectors", with
  `N ≈ (4/3) * pi * (radius_pc / e)^3`, rounded to two significant
  figures. Sectors the outline doesn't allow are skipped by generation,
  so the line says "up to". About 1,900 at 100 ly, and about 15,000 at
  200 ly.
- Above 5,000 sectors, the Start button asks for a confirmation that it
  may run for a long time.

**As built (7.53.0):** the buttons and the radius dialog work as above
(default 100 ly, 13 to 652 ly, the estimate, and the confirmation past
5,000 sectors); the radius goes as `slot_radius_pc`. Generate this layer
or slab is not built: it waits for the `--block` mode. The progress line
below is not built either; the form still posts and follows the redirect
to the job page.

**After starting (planned):** the form posts to the admin Generate page as today.
Instead of following the redirect, the map sends it with `fetch` and
reads the job URL. The panel then shows a progress line (polling the
same status endpoint `generatejobs.js` uses) with a link to the job
page. When the job ends, the stage's contents are fetched again: the
stamp changes, and only the touched chain is invalidated. The fallback,
without the `fetch` path, is today's redirect to the job page.

## 7. Data: the stage contents API

`GET /api/galaxy/stage?at=m.I.s.S` returns one container's children with
generated counts. With no `at`, it returns the galaxy's level-243
blocks.

```json
{"at": "27.63.115.-1", "child_m": 3,
 "children": [{"ring": 568, "wedge": 1036, "slab": -7, "generated": 12}, ...],
 "sectors": null}
```

At `m = 3` the children are sectors, and `sectors` lists each generated
one: `{ring, layer, slot, id, name, system_count}`.

- **Query (as built, `queryDb.galaxy_stage`):** for a block, one
  `SELECT ... FROM sectors WHERE layer_index BETWEEN ? AND ? AND (...)`,
  where the bracket holds one `(ring_index = ? AND ring_slot_index BETWEEN
  ? AND ?)` clause per member ring, with each ring's slot range worked out
  first by section 3.5's integer formula (`_block_slot_range`). The rows
  are then grouped by child block in Python. At `m = 3` it also counts each
  sector's systems. Without `at`, it groups every generated sector by
  `(ring_index, FLOOR((layer_index + 121) / 243), ring_slot_index)` and maps
  each group to its level-243 block in Python (`drill_chain_of`). That
  returns no more rows than there are generated sectors. Children with no
  generated sectors are left out of `children`.
- **Caching:** the site's `/galaxy/stage` reads through
  `lib/tilecache.fetch_stage`, a disk cache under the galaxy stamp, the
  same way as tiles. A changed sector deletes only its ancestor chain's
  stages, which `/api/galaxy/changes` lists as `stages`.
- **Scale note:** if generated sectors ever reach the millions, the
  no-`at` summary should come from a per-level count table kept up to
  date on generation (a Database thread item at that point). No schema
  change is needed now.
- Totals (allowed sectors per block) are computed in the browser from
  the analytic bound, as today, so the API returns only generated
  counts.

## 8. URLs, history and bookmarks

### 8.1 Stage URLs

| Stage | URL |
| --- | --- |
| 1 | `/galaxy` |
| 2 | `/galaxy?slab=0` (a level-243 slab) |
| 3, 5, 7 | `/galaxy?at=243.7.14.0` (`m.ring.wedge.slab`) |
| 4, 6, 8 | `/galaxy?at=243.7.14.0&slab=-1` (a slab of that block's children, numbered at the child level) |

- Every stage change does `history.pushState`, and `popstate` restores
  the stage with the flight played in reverse. A reload or a shared link
  opens that stage directly, with the breadcrumb rebuilt by walking
  `drillParent`.
- The page validates `at` and `slab` (`galaxystages.parseStageQuery`:
  the block must exist and the child slab must be inside it). Anything
  invalid opens stage 1, with a notice saying why. The server does not
  check them; it only serves `/galaxy/stage` for a valid `at`.
- `?sector=<designation>` opens the stage 8 that holds that sector, with
  it selected (built in 7.44.0). Sector, system and search pages are to
  link this way (MAP.25, not built; the sector page still links to the
  Quadrant table).

### 8.2 Bookmarks

Not built yet (MAP.23).

- A ☆ button on the breadcrumb and on each info panel saves the current
  stage, sector, system or phenomenon:
  `{name, kind: "stage"|"sector"|"system"|<phenomenon type>, value, created}`,
  where `value` is a stage URL, a designation or a NAV endpoint
  (`system:12`).
- Stored per browser in `localStorage["planetgen.bookmarks.<db name>"]`,
  up to 100 entries, with every read and write wrapped so that a browser
  with storage blocked just shows no bookmarks. Shared bookmarks in the
  database are a separate option (decision 4).
- A Bookmarks menu sits on the map. Ctrl+1 to Ctrl+9 open the first
  nine. Bookmarks can be renamed and deleted from the menu.
- The NAV pickers list system and phenomenon bookmarks (section 9). A
  sector bookmark opens the system picker for that sector.
- One module, `static/bookmarks.js`, is shared by the map, the Sector
  Map and the NAV page.

## 9. NAV integration

The NAV page keeps its sector-then-system dropdowns. Each endpoint gains
the following (the NAV page part, MAP.22, is not built yet; 9.1 to 9.4
are):

1. **Pick on Galaxy Map**, which links to `/galaxy?pick=from&to=system:40`
   (or `pick=to&from=...`), carrying the other endpoint along.
   - The map shows a banner, "Choosing a destination · Cancel" (Cancel
     returns to `/nav` with the endpoint already chosen).
   - In pick mode "Generated only" is forced on, since NAV endpoints are
     systems and phenomena, which exist only in generated sectors.
   - At stage 8, clicking a sector goes to
     `/sector/<id>?pick=to&from=system:12`.
2. **Pick in this sector** appears once the other endpoint is known. It
   opens that endpoint's own sector straight in pick mode, so a pick
   inside one sector uses only the sector interface.
3. **Bookmarks:** a select listing the system and phenomenon bookmarks.

### 9.1 Sector Map pick mode

Built in 7.58.0 (`web/sector_page.py`, `_pick_mode`).
`/sector/<id>?pick=to&from=...` shows the same banner. Clicking a
system or phenomenon adds a **Use as destination** (or start) button to
its info panel. That button links to
`/nav?from=system:12&to=system:40`, which lands on the plotted course.
The endpoint strings are `nav_page.endpoint(kind, id)`.

### 9.2 Nav links from everywhere

System and phenomenon pages, and the Sector Map's info panel, get **Nav
from here** and **Nav to here**. The NAV page already accepts `to`
without `from`. Built: the system and phenomenon pages had "Navigate
from/to here" already, and the Sector Map's panel gained them in 7.58.0.

### 9.3 Address bar

A field over the breadcrumb takes a designation, `312/-3/1042` or
`ring 312 layer -3 slot 1042`, `x, y, z` in pc, or a sector or system
name. It flies to that sector's stage 8 (section 5.3).

As built (7.50.0): the page parses addresses itself
(`galaxystages.parseAddress`). A name goes to the site's `/galaxy/locate`,
which calls `GET /api/galaxy/locate` (`queryDb.galaxy_locate`: sectors and
systems named like the text, each with its sector address), not
`/api/search` as first planned. A name with several matches lists them to
pick from, and anything that can't be a sector says why. Pick mode on the
Galaxy Map is not built, so it does not work there yet.

### 9.4 The course on the map

The NAV result offers **Show on Galaxy Map**. That opens the smallest
stage holding both endpoints, with a line between them and the endpoints
marked.

As built (7.52.0): the link is `/galaxy?course=<from>,<to>`;
`web/galaxy_views.py` asks `nav_page.galaxy_course` for the waypoints in
galaxy-frame parsecs, and the map draws the course through its stops,
each ringed and the two ends named. A course that stays inside one sector
opens that sector instead. MAP.26 asks for every "Show on Galaxy Map"
link, this one included, to open at the sector level.

## 10. What changes on today's map

- Click-to-center and double-click zoom are replaced by the stages. The
  free continuous zoom stays only if decision 2 keeps a "Free look"
  toggle. (As built: the map opens on the stages and keeps Free look as
  a button; MAP.17 would remove it.)
- Wedge lines, density shading, the filled-share look, the sector/pc/ly
  scale readout and the info panel all stay.
- Each stage draws at most about 900 blocks, so the phone-performance
  follow-ups in MAP.1 matter only for free look. (MAP.1 shipped in
  7.42.1.)
- MAP.3 (bigger map, controls underneath) should land first or
  alongside: the breadcrumb, address bar and slab strip need the room.
  (It shipped in 7.55.0, after the stages.)
- Added since: the Territories overlay (7.54.0, population-and-politics.md)
  and the NAV course line (section 9.4) draw over the map.

## 11. Decisions for Boss

Each has a default, and work can start on it.

1. **Ladder:** 243 → 27 → 3 → 1. *(Decided by Boss, 2026-10-01: "the
   bigger targets".)*
2. **Free camera.** Default: drag-rotate inside the 3D stages only, with
   no free fly. Alternative: keep today's free zoom as a "Free look"
   toggle. *(As built in 7.44.0: the stages drag-rotate and Free look is
   kept as a button until Boss decides. Boss's MAP.17 of 2026-10-01
   goes further than the default: no Free look and no drag-rotate at any
   stage, starting top-down with a wedge pick. Not built yet.)*
3. **Bigger generate buttons.** Default: sector, layer or slab, and
   neighborhood at stages 7-8 only. Option: "Generate this block" at
   stage 5 (up to about 19,000 sectors) behind a confirmation.
4. **Bookmarks.** Default: per browser. Alternative: shared in the
   database (needs a migration).
5. **NAV picking.** Default: a round trip through the map and Sector Map
   pages back to NAV. Alternative: the map embedded in the NAV page.
6. **Orientation of block slabs.** Default: outward up (core at the
   bottom), with a compass arrow. Alternative: bearing 000 to the right
   everywhere, as in stage 2.
7. **Neighborhood radius default.** Default: 100 ly, the same as today's
   button. *(Built with this default in 7.53.0.)*

## 12. Build order and owners

| TODO | Piece | Owner | Depends on |
| --- | --- | --- | --- |
| MAP.28 | Nested ladder geometry, JS and Python, with a parity test (section 3). **Built, 7.41.2, PR #160** | Galaxy Map | none |
| MAP.29 | Stage contents API and caching (section 7). **Built, 7.41.3, PR #160** | Galaxy Map | MAP.28 |
| MAP.16 | The stages: views, hover, pull-out, flight, breadcrumb, URLs, keys, touch (sections 4, 5, 8.1). **Built, 7.44.0, PR #171** | Galaxy Map | MAP.28, MAP.29, MAP.3 |
| MAP.20 | Admin generation at the sector level: `--block` mode, Generate page form, map buttons, radius dialog, progress (section 6). **Map buttons and radius dialog built, 7.53.0, PR #177**; `--block`, the form and progress open | Web (`generate.py`, Generate page) and Galaxy Map (buttons) | MAP.28 (Python), MAP.16 |
| MAP.21 | Sector Map pick mode and Nav from/to links (sections 9.1, 9.2). **Built, 7.58.0, PR #178** | Web | the URL formats only |
| MAP.22 | NAV page: Pick on Galaxy Map, Pick in this sector, Bookmarks (section 9) | Web | MAP.16, MAP.21 |
| MAP.23 | Bookmarks (section 8.2) | Galaxy Map (module, map menu) and Web (NAV, Sector Map) | MAP.16 |
| MAP.24 | Address bar (section 9.3). **Built, 7.50.0, PR #172** | Galaxy Map | MAP.16 |
| MAP.25 | "Show on Galaxy Map" links with `?sector=` (section 8.1) | Web | MAP.16's URL format |
| MAP.27 | Course on the Galaxy Map (section 9.4). **Built, 7.52.0, PR #176** | Galaxy Map and Web | MAP.16, MAP.22 |
| MAP.17 | No free camera; drill down top-down by wedge, slice and block (bug against MAP.16; decision 2) | Galaxy Map | MAP.16 |
| MAP.26 | "Show on Galaxy Map" opens at the sector level; the map's own Back and Forward (bug; folds in MAP.25) | Galaxy Map and Web | MAP.17, MAP.25 |

MAP.3 (the bigger map) shipped in 7.55.0. MAP.27 shipped before MAP.22, using
the NAV result's link rather than the NAV pickers. What is left is MAP.20's
`--block` mode, MAP.22, MAP.23 and MAP.25, and the two bugs MAP.17 and MAP.26, which
change sections 4, 5 and 8.1 once Boss's open questions in them are
answered.

## 13. Sources

- "Mapping Volumetric Galactic Topologies to Planar Interfaces", the
  report Boss shared in the project on 2026-10-01 (not in the repo).
- J. J. van Wijk and W. A. A. Nuij, "Smooth and efficient zooming and
  panning", IEEE Transactions on Visualization and Computer Graphics
  9(4), 2003, pp. 424-431. The source of the formulas in section 5.3.
- This repo: `static/galaxyprisms.js` (`blockWedgeCount`,
  `blockSlotRange`, `blockSectorCount`), `stellarObjects/galaxyGeometry.py`
  (`ring_master_count`, `ring_sector_count`), `stellarObjects/generationLimits.py`,
  `web/generate_page.py`, `web/nav_page.py`.

## 14. Why it works this way

- **Discrete stages instead of free zoom.** Boss asked to select in
  sections: click a block, zoom into just its contents, and repeat down to a
  sector (section 1). The research paper Boss shared supports discrete scale
  tiers, a breadcrumb and animated flights (section 2).
- **243, 27, 3, 1.** Boss chose "the bigger targets" over 81, 9, 1
  (decision 1). One result, noted in section 5.1, is that each stage draws
  at most about 900 blocks, so the free camera's 60,000-block budget never
  matters.
- **A nested wedge rule.** With each level's wedge count picked on its own,
  child wedges sat inside one parent only in some rings (109 of 143 rings
  at 243 to 27, 1,468 of 3,856 at 3 to 1). A drill-down needs every child
  inside exactly one parent, so each child count is a whole multiple of its
  parent's. The cost is a slightly wider arc error (section 3.4).
- **Sectors join the level-3 block holding their center.** Sector slots do
  not line up with level-3 wedges (only 1 of 1,286 rings does), so exact
  nesting is impossible at the last step; the center rule is the one
  `blockSlotRange` already used.
- **Generated counts from the server, totals in the browser.** The
  allowed-sector totals follow from the stored outline and are cheap to
  compute on the page; only generated counts need the database. The
  stage query reads at most one row per generated sector.
- **A Free look button, for now.** Decision 2 was still open when the
  stages shipped, so the old camera stayed behind a button rather than
  being removed (7.44.0). MAP.17 is Boss's answer: remove it.
- **Addresses parsed on the page, names on the server.** A designation,
  ring/layer/slot or coordinates map to a sector with pure math the page
  already has, so only a name needs the database (`/galaxy/locate`).
- **Same math in Python and JavaScript, with a parity test.** The page and
  the server must agree on every block key; `test_galaxydrill.py` runs both.
- **three.js stays.** The map keeps the vendored three.js build; the reasons
  are in `docs/html-interface.md` ("Why the maps use three.js") and in
  `design-decisions.md`.

### Ideas from the paper that were rejected

Section 2's table records them: a 2D canvas replacing 3D (the 3D solid stays
for picking slabs), unwrapping each cylindrical shell into a strip (a ring
is one sector wide, so the strip is a thin ribbon that hides radial
neighbors), Morton hashing (every address is already closed form and the
designation is already a reversible key) and fisheye lenses (they would
bend the distances and bearings NAV reports).
