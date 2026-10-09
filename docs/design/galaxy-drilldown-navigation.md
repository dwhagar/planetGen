# Galaxy Map drill-down navigation

Boss's design for getting around the Galaxy Map, recorded 2026-10-01.

**Status (2026-10-02):** built. The arc pick (MAP.85) has replaced the
galaxy's first pick (section 4); section 15's later steps (slab, then
segment, MAP.56 onward) are still planned.

| Piece | Section | TODO | Built in |
|---|---|---|---|
| Block ladder (`planetgen/galaxy/drill.py`, `drill*` in `static/galaxyprisms.js`) | 3 | MAP.28 | 7.41.2, PR #160 |
| Stage contents API (`GET /api/galaxy/stage`, the site's cached `/galaxy/stage`) | 7 | MAP.29 | 7.41.3, PR #160 |
| The stages (`static/galaxystages.js`, `static/galaxystageview.js`), stage URLs, breadcrumb, keys, touch, "Charted only" | 4, 5, 8.1 | MAP.16 | 7.44.0, PR #171 |
| Top-down only, quarter, layer and arc picks, dimmed hover, wedge lines kept to the view, opening at a sector's layer, the map's Back and Forward | 4, 5, 8.1, 10, 11 | MAP.17, MAP.18, MAP.19, MAP.44, MAP.26 | 7.73.177 |
| Address bar (`/galaxy/locate`) | 9.3 | MAP.24 | 7.50.0, PR #172 |
| Course on the map (`/galaxy?course=<from>,<to>`) | 9.4 | MAP.27 | 7.52.0, PR #176 |
| Map Generate buttons and the light-year radius dialog | 6 | MAP.20 (map side) | 7.53.0, PR #177 |
| Sector Map pick mode and Nav from/to links | 9.1, 9.2 | MAP.21 | 7.58.0, PR #178 |
| Generate this layer or slab (`planetgen galaxy --block`, the Generate page's block form) | 6 | MAP.20 | 7.60.174 |
| NAV page pickers: Pick on Galaxy Map, Pick in this sector, the Bookmarks select | 9 | MAP.22 | 7.94.253 (the picks earlier) |
| Bookmarks (`static/bookmarks.js`): the breadcrumb ☆, page buttons, the map menu, keys 1-9 (MAP.81) | 8.2 | MAP.23 | 7.94.253 |
| "Show on Galaxy Map" links with `?sector=` from sector, system and search pages | 8.1 | MAP.25 | 7.64.174 |
| The whole galaxy in 3D with no grid lines, the arc pick in place of the quarter, outlines along real wedge lines, one scale line, the controls' Menu | 4, 5, 8.1, 10, 15 | MAP.85, MAP.52, MAP.60, MAP.55 | PR #369 |

Later requests changed this design, and the sections below describe the
map as it now is: MAP.17 (no free camera and no rotation at all; start
top-down and pick a quarter, then a slice), MAP.19 (every pick is a big
arc, not a single block), MAP.18 (everything but the hovered pick is
dimmed), MAP.44 (wedge lines kept to the part in view) and MAP.26 ("Show
on Galaxy Map" opens at the sector's own layer, and the map has its own
Back and Forward). Boss settled what a slice and a region are on
2026-10-01 ("Layer + arc", section 11, decision 8). The open work items
are in `docs/TODO.md`, and each one points back to a section here.

MAP.85 (Boss, 2026-10-01 23:53Z; arc size decided 2026-10-02 01:53Z)
then made the whole galaxy a 3D disk the visitor can turn, with no
sector, block or wedge lines drawn on it, and replaced the quarter with
an arc (section 4); MAP.55 put every control but Back, Forward, Up,
Reset and Bookmarks into a Menu and dropped the Wedges button and the
wedge lines; MAP.60 cut the scale readout to one line.

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
- `planetgen/galaxy/drill.py`: the same functions in Python, for
  the API (section 7) and generation (section 6), plus
  `format_drill_key`/`parse_drill_key` for the `m.ring.wedge.slab` keys
  (`formatDrillKey`/`parseDrillKey` on the page). Built 2026-10-01
  (PR #160). Generation (section 6) does not use it yet.
- `tests/test_galaxydrill.py` runs both on every block ring and on
  sampled sectors and checks that they agree, and that a parent's
  sectors are exactly the sectors of its children.

## 4. The stages

The whole galaxy is a 3D disk, opened straight down so the spiral arms
show (MAP.97), that the visitor can turn any way, and zoom in to about
twice as close but no further out than its fit (MAP.85, MAP.58's range);
it has no sector, block or wedge lines, so the stars and the spiral show
through. Below it each zoom step flies to a preset for what it shows
(MAP.97, Boss 2026-10-02): a block of several slabs (an arc, an entered
block, the cube of sectors) at an isometric slant (Boss, 2026-10-01,
"turn the angle to be isometric so that the user can select layers
easier"), so layers can be clicked on the map as well as picked with the
slab buttons, and a slab straight down. Every view can be turned any way
by any amount, through edge-on and under the plane (MAP.96), moved and
zoomed within a
short range, to make layers, blocks and sectors easier to pick (Boss,
2026-10-01, section 5.1); each step still opens on its own view. A stage is a container (the galaxy or
one block of the ladder, section 3) and the picks made inside it so far.
The picks go:

- **Arc** (MAP.85), at the galaxy only: a piece of the disk about 40°
  of bearing by a third of its radius (inner, middle or outer), running
  the whole height of the disk. An arc holds the level-243 blocks of
  its third of the block rings whose middle bearing falls in one of
  eight 45° bins: 45° is the width nearest 40° whose lines are wedge
  lines of every block ring past the core (their wedge counts are all
  multiples of 8), so an arc's sides run straight from its inner edge
  to its outer one; the core's three 120° blocks join the arc their
  middle falls in. That makes 24 arcs. Hovering one outlines it, top and
  bottom with its corners joined, and its neighbors faintly; nothing
  else is outlined. Older links' quarters (`q<n>`) still open.
- **Layer** (MAP.17's "slice"): a layer of the disk. While the view
  holds more than three slabs of the container's children, the choice
  is between the lowest, middle and highest third of them; with three
  or fewer, a single slab. Layers are picked from the strip beside the
  map, since from above one layer hides the others.
- **Region** (MAP.19's "arc"): one cell of a 3 by 3 grid over the view,
  three ring bands by three arcs of the view's bearing span (fewer where
  the view has fewer rings or blocks across), so each arc is about a
  third of the wedge in view. A block belongs to the cell its middle
  falls in. The map zooms into the region picked.
- An arc, quarter or region reaches as far as its blocks do, not just its
  even share of the bearings (Boss, 2026-10-01): near the core a block
  sits between meridians that run all the way in and can be wider than
  the share (quarter 0°–90° holds the 0°–120° core wedge), so the
  highlight, the label and the zoom all follow the real wedge.

They alternate: after the arc, a layer, then a region, then a layer
and so on (galaxystages.nextPickKind). A layer comes next whenever the
last pick wasn't a layer and the view still spans more than one slab;
otherwise a region while the view is more than one column of blocks
wide; otherwise a layer. A choice with only one option is taken without
asking, and when the view is down to a single block it is entered: its
children become the view and the picks start again inside it
(galaxystages.settleStage). In the end the view is the sectors of one
layer of a level-3 block, about nine of them, and the last pick is a
sector.

A sector deep in the disk takes about a dozen picks, each a big target.
For example (a real chain under the default galaxy): Galaxy, Quarter
90°–180°, Slab 0, Arc, Arc, then inside block 7·14 Slabs -1 to 1, Arc,
Slab -1, then inside block 63·115 Slabs -7 to -5, Arc, Slab -7, then
inside block 568·1036 Layer -20, and the sector.

Two shortcuts (Boss, 2026-10-01):

- **A block one sector tall shows its sectors.** A level-27 block whose
  sectors all lie in one layer (at the disk's top and bottom faces) skips
  its level-3 blocks: entering it shows its sectors at once, narrowed by
  arcs (galaxystages.thinSectors). Their generated counts come from each
  level-3 child that holds any.
- **The last cube is picked in 3D.** Once the view is the sectors of a
  level-3 block across several layers (27 at most), it is shown from a
  fixed 55° slant with its layers touching (no space between blocks or
  layers anywhere on the map, Boss 2026-10-01). The next pick is still a
  layer, so hovering any sector lights and outlines its whole layer
  (the others fade) and a click picks that layer, the same one the strip
  marks (MAP.91); the layer then shows its sectors one by one. The cube
  can be turned like any view below the quarters.

The admin Generate tools (section 6) appear once the view is one layer
of a level-3 block.

## 5. How the stages behave

### 5.1 The view

- Only the blocks in view are drawn: the container's children, narrowed
  by the picks so far. That is at most a few hundred blocks.
- **Camera:** straight down (-z). The view fits the footprint of the
  blocks in view: camera height `h = 1.1 * R_fit / tan(FOV/2)`, where
  `R_fit` is the largest distance from the footprint's center to its
  corners, with the narrower canvas side used for the FOV. The whole
  galaxy has galactic north up and bearing 000 to the right; every other
  view is turned so the middle of its bearing span points up (the core
  toward the bottom, decision 6).
- **Turning and moving** (below the galaxy and its quarters): drag turns
  the view (tilt up to 80° from straight down), right-drag or
  Shift-drag moves it (its middle stays within 1.5 fits of the stage's
  own), and the wheel or a pinch zooms from an eighth of the stage's fit
  to 2.5 times it. A drag never picks. Where the view can be turned, a
  layer can also be clicked on the map (the block under the pointer
  picks its layer). Reset view flies back to the stage's own view, and
  every step to another stage (into a pick, Up, Back, a crumb) opens on
  that stage's own view. At the galaxy and its quarters the wheel
  scrolls the page.
- **Only the wedge in view shows** (MAP.44; Boss, 2026-10-01: "only show
  the wedge"): the wedge lines, stars and clouds are kept to the part of
  the galaxy in view (its rings, bearings and layers) and stop there
  sharply. Over the
  whole galaxy they run to its edge (MAP.43) and carry their bearing
  labels. Ring boundaries are the faces of the blocks in view, so they
  never reach past it.

### 5.2 Picking

- **Hover** (MAP.18): the choice under the pointer (a quarter or an
  arc, with all its blocks) stays as it is, every other choice drops to
  25% opacity, and an accent outline goes round its area. A tooltip
  names it, for example "Arc 30°–60°, 5832–10692 pc from the core, 0 of
  338,081,040 sectors generated". The dimming clears when the pointer
  leaves the map. It switches at once, with no fade, so there is nothing
  for `prefers-reduced-motion` to turn off.
- **Click:** takes the choice, and the map flies into it (section 5.3).
  Inside a block, a click on a black hole, neutron star, quasar or small
  cloud shows it instead. Stars never take the click (MAP.101): a click
  on one picks what's under it, since in a dense sector they would cover
  it, and a star's details are on its sector's page.
- **The slab buttons (MAP.30, now MAP.54 and MAP.76):** while a layer
  is to be picked, one button per slab to the right of the map (in one
  column below it on a phone), each with its slab's name and generated
  count, and a leader line from the button to its slab on the map,
  redrawn whenever the view turns, zooms, pans or resizes. The buttons
  are always in slab-number order (MAP.110), running the same way as the
  slabs on screen (the highest first, or the lowest first once the view
  turns under the plane), and each line ends no higher than the one
  above it so the lines don't cross (on a phone the lines run up lanes
  along the map's right edge, nested, so that column runs from the slab
  lowest on screen down); a slab off the map
  gets a line ending in an arrow at the map's edge. Hovering or focusing
  a button lights its slab and line, clicking it takes the slab, and
  hovering a slab on the map lights its button and line. The lines are
  faint until lit. At other stages the box says which layers the view
  holds. Hovering any cube
  of a slab on the map lights the whole slab and outlines it round its
  full height, at every level including the 3 by 3 by 3 view (MAP.91,
  Boss 2026-10-02). This is the paper's
  elevation panel, and it gives keyboard and screen-reader users the
  same choice. The map beside it is 4:3 and no taller than the window
  (1:1 on a phone).
- Choices are big: four quarters, one slab at a time, at most nine
  arcs, so each is at least about a ninth of the map, finger-sized on
  a phone.
- **Touch:** the first tap highlights a choice and a second tap on it
  takes it.
- **Keys:** ← → move through the choices, ↑ ↓ to the next one along the
  same arc (or the next layer up or down), Enter takes the highlighted
  one, Esc or Backspace goes one step back out, and Home returns to the
  galaxy.
- A **"Charted only"** toggle (MAP.111) dims blocks with nothing
  generated (they stay pickable, MAP.112), dims the stars outside charted sectors
  to a fifth of their opacity, and outlines each charted block, so
  existing content is easy to follow.

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
view also turns to the new view's bearing; it stays looking straight
down. With `prefers-reduced-motion` the camera cuts straight to
the end. The same path is used for the address bar's jumps
(section 9.3) and for Back.

### 5.4 Breadcrumb, Back and Forward

- A breadcrumb above the map shows the chain, one crumb per step:
  `Galaxy › Quarter 90°–180° › Slab 0 › Arc 120°–150° › … › Block
  568·1036 › Layer -20 › Sector 1,705·-20·3,225`. Blocks read
  `ring·wedge`, and the last crumb is the address `ring·layer·slot`.
  Clicking a crumb returns to that step. The breadcrumb is always one
  line (MAP.93): when the steps don't fit, it shows the first, a "…"
  button whose menu lists the steps it hides, and as many of the last
  steps as fit, then the current one, measured again whenever the
  line's width changes; the ☆ stays on the line. At phone width
  (under 600 px) the line gives way to a round Steps button between
  Back and Forward whose menu lists every step, the current one marked
  (MAP.94); Reset stays beside the arrows.
- **Back** and **Forward** buttons beside the map (MAP.26) step through
  the stages visited on this map. They use the browser's own history:
  every stage change does `history.pushState`, with the map's own index
  in the entry, so the buttons and the browser's Back agree, a reload
  keeps the stage, and a visit to a sector's page and back returns to it.
  Back is off at the first stage of the visit and Forward once nothing
  is ahead.
- **Up** goes one step back out (the same as Esc), and **Reset** starts
  over at the whole galaxy (Home). Back, Forward, Up, Reset and
  Bookmarks are the only buttons in view; Reset view, Charted only
  and Territories are in the Menu beside them (MAP.55).

### 5.5 The sector at the end

Once the view is one layer of sectors:

- A generated sector opens in place (MAP.66): its stars, phenomena and
  bodies are drawn where the sector sits, the camera flies to fit it, and
  hover, tooltip, ring, info panel, ☆ and NAV links work as on the sector
  page. This is the last stage; it has its own URL,
  `/galaxy?sector=<designation>&open=1`, and Up closes it (back to
  `?sector=<designation>`, the sector still selected). A click on a
  neighboring sector steps sideways. The sector's own page, `/sector/<id>`
  (tables, text, edit tools), embeds the same engine locked to that sector
  (MAP.68): the Sector Map there is this map, with no steps or history of
  its own, and a click on a neighboring sector opens that sector's page.
- A star of an open sector can be opened in place too (MAP.125): its info
  panel's "Open system here" loads the system's scene
  (`/api/systems/<id>/scene`, `galaxysystem.js`), draws its stars, planets,
  moons, belts and comets at the star's place in the sector (compressed
  scale, or true scale from the info panel), and flies the camera to fit
  it. Its bodies are picked like any other object and "Fly to" takes the
  camera down to one and follows it. URL `&system=<id>` after `open=1`
  (and `&object=<ref>` for the body picked); Up closes the system (back to
  the sector), Back and Forward step through. The wheel carries the zoom
  across: zoomed right in on a sector, it opens the star picked (else the
  one nearest the middle of the view); zoomed right out of a system, it
  returns to the sector.
- A sector that isn't generated shows today's cell panel (address,
  designation, coordinates, 8 corners). Admins also get the Generate
  tools (section 6).

## 6. Generating from the map (admin)

Only an admin with current credentials sees these, as with PR #151 (the
server puts the `generate` target in the scene data). They appear
inside a level-3 block.

| Button | What runs | Size |
| --- | --- | --- |
| Generate this sector | `planetgen galaxy --ring i --layer j --slot k` (today's `slot` mode) | 1 |
| Generate this layer (once one layer of sectors is shown) | New `--block 3.J.s.S --block-layer j`: the block's ungenerated allowed sectors on that layer (section 3.5) | Up to about 15 |
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
  `planetgen.physics.units`).
- A live estimate under the field: "about N sectors", with
  `N ≈ (4/3) * pi * (radius_pc / e)^3`, rounded to two significant
  figures. Sectors the outline doesn't allow are skipped by generation,
  so the line says "up to". About 1,900 at 100 ly, and about 15,000 at
  200 ly.
- Above 5,000 sectors, the Start button asks for a confirmation that it
  may run for a long time.

**As built (7.53.0):** the buttons and the radius dialog work as above
(default 100 ly, 13 to 652 ly, the estimate, and the confirmation past
5,000 sectors); the radius goes as `slot_radius_pc`. Since GEN.23 the dialog's default
is 39 ly (about 12 pc, `DEFAULT_GENERATE_RADIUS_PC`). Since 7.60.174 a
level-3 block also offers "Generate this layer" and "Generate this
block" (`blockGenerateButtons` in `static/generatebuttons.js`), posting
the Generate page's block mode (`--block`, `--block-layer`). The
progress line below is not built; the form still posts and follows the
redirect to the job page.

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
| The galaxy | `/galaxy` |
| Picks at the galaxy | `/galaxy?p=a1.90,L0,r4` |
| A block, and picks inside it | `/galaxy?at=243.7.14.0&p=L-4~-2,r4` (`m.ring.wedge.slab`) |
| A sector, selected | `/galaxy?sector=<designation>` |
| A sector, opened in place | `/galaxy?sector=<designation>&open=1` (MAP.66) |
| A sector's page | `/sector/<id>`: the same engine locked to the sector, no query of its own (MAP.68) |
| A system, opened in place | `/galaxy?sector=<designation>&open=1&system=<id>` (MAP.125) |
| A body of the open system | `/galaxy?sector=<designation>&open=1&system=<id>&object=<ref>` (MAP.125; ref like `planet:12`) |

- A pick reads `a<band>.<bearing>` (the arc of band 0, 1 or 2, inner
  first, starting at that bearing in degrees), `q<n>` (quarter n, from
  older links),
  `r<n>` (region n of the 3 by 3 grid, band by arc, inner band first),
  `L<s>` (slab s of the container's children) or `L<lo>~<hi>` (slabs lo
  to hi).
- Every stage change does `history.pushState` (section 5.4). A reload or
  a shared link opens that stage directly, with the breadcrumb rebuilt
  by walking from the galaxy.
- The page checks `at` and every pick against the galaxy
  (`galaxystages.resolveStage`). Anything invalid opens the galaxy, with
  a notice saying why. An older link's `?slab=s` reads as a pick of that
  one slab. The server does not check them; it only serves
  `/galaxy/stage` for a valid `at`.
- `?sector=<designation>` (MAP.26) opens the sector level of the slice
  holding the sector: its level-3 block with the sector's own layer
  picked, the sector among its neighbours, highlighted and shown in the
  info panel. The sector, system and search pages link this way
  (MAP.25), and the NAV course still fits both ends (section 9.4).
- `hide=<kind>,<kind>` (MAP.79) lists the kinds of object left off an
  opened sector's map (`star`, `nebula`, `supernovaRemnant`,
  `asteroidField`, `blackHole`, `neutronStar`, `quasar`, `roguePlanet`,
  `interstellarComet`, `neighbor`). It rides on every stage URL and on a
  sector page's, and the Menu's "Show on the map" buttons change it in
  place (no history entry). Unknown names are ignored.
- `mark=<kind>,<kind>` (MAP.123) lists the phenomenon kinds drawn larger
  and brighter (the same names as `hide`, but never `star`); the Menu's
  "Highlight" buttons change it in place. A bookmark of a stage view keeps
  `hide`, `mark`, `stars`, `lum` and `color`.
- One scheme (MAP.67): a map's position is always a query of `/galaxy`
  (`at`, `p`, `sector`, `open`) and every move is a history entry, so
  Back, Forward, a reload and a bookmark behave the same at every level
  down to the open sector. The sector page is the one place the map has no
  query of its own, because the page's URL already names the sector and the
  map is locked to it; a pick kept for NAV (`pick`, `from`, `to`) rides on
  every URL either way.

### 8.2 Bookmarks

Built (MAP.23), with decision 4's default: per browser, no migration.

- A ☆ button on the map's breadcrumb (`galaxystageview.js`) saves the
  selected sector, or else the current stage; a ☆ Bookmark button on the
  system, phenomenon and sector pages (`templates/partials/bookmark.html`,
  its entry from `web.helpers.bookmark`) saves that page. An entry is
  `{name, kind: "stage"|"sector"|"system"|<phenomenon type>, value, created}`,
  where `value` is a stage URL (path and query), a sector designation (a
  sector with no galaxy address uses its page URL) or a NAV endpoint
  (`system:12`, `nebula:3`, from `nav_page.endpoint`), plus `url`, the
  page it opens, and `sectorId` for a generated sector. Kind and value
  identify it, so the map's ☆ and the sector page's button find the same
  sector. A saved place shows ★ ("Bookmarked"); pressing it again asks
  to remove it.
- Stored per browser in `localStorage["planetgen.bookmarks.<db name>"]`,
  up to 100 entries, with every read and write wrapped so that a browser
  with storage blocked just shows no bookmarks (and saving says so). The
  database name reaches the script as `data-bookmark-db` on the button,
  the menu or the NAV form. Shared bookmarks in the database remain the
  option decision 4 leaves for later.
- A Bookmarks menu (a `<details>`) sits in the map's controls row and
  opens in the row's flow, a full row wide, so it fits a phone. Each
  entry is a link with Rename (in place: Enter saves, Escape cancels) and
  Delete; Escape closes the menu. The keys 1 to 9 open the first nine
  while the map (its canvas or a control in its panel) has the focus,
  but not from a text box or select (MAP.81: they were Ctrl+1 to Ctrl+9,
  which Chrome and Firefox on Windows and Linux keep for switching tabs).
- The NAV page's Bookmarks select lists system and phenomenon bookmarks
  (section 9). A sector bookmark opens the system picker for that sector,
  and a saved map view opens the Galaxy Map there to pick on. Once both
  ends are set, the course page offers the same select for a new start
  or a new destination.
- Bookmarks keep a course pick (NAV.40): while a start or destination is
  being chosen, the Galaxy Map's menu and its 1 to 9 keys, and a
  Bookmarks menu under the sector page's pick banner, carry the pick
  (`data-pick`, `data-keep-name`, `data-keep-value`, `data-nav-url`). A
  system or phenomenon bookmark sets that end on the NAV page, and a
  sector page, a sector on the map or a saved map view opens in pick
  mode, all keeping the other end. The Galaxy Map keeps the pick on its
  own stage URLs too.
- One module, `static/bookmarks.js`, serves the map, the system,
  phenomenon and sector pages and the NAV page; importing it wires
  whatever each page has, with no inline script. The Sector Map's own
  info panel has no ☆ yet; its page's button covers the sector.

## 9. NAV integration

The NAV page keeps its sector-then-system dropdowns. Each endpoint gains
the following (all built: MAP.22, with sections 9.1 to 9.4):

1. **Pick on Galaxy Map**, which links to `/galaxy?pick=from&to=system:40`
   (or `pick=to&from=...`), carrying the other endpoint along.
   - The map shows a banner, "Choosing a destination · Cancel" (Cancel
     returns to `/nav` with the endpoint already chosen).
   - Every control works as outside pick mode, "Charted only" included
     (NAV.32). NAV endpoints are systems and phenomena, which exist only
     in generated sectors, so hover still lights every arc and choice as
     outside pick mode (NAV.31); one holding nothing generated says so
     and can't be taken.
   - At the sector level, clicking a sector goes to
     `/sector/<id>?pick=to&from=system:12`.
2. **Pick in this sector** appears once the other endpoint is known. It
   opens that endpoint's own sector straight in pick mode, so a pick
   inside one sector uses only the sector interface.
3. **Bookmarks:** a select listing the system and phenomenon bookmarks
   (`nav_page._bookmark_pick`, filled by `static/bookmarks.js`), beside
   the map picks at each step, hidden while there are none it can use.
   "Start Here" (or "End Here") goes to `/nav` with that endpoint set
   and the other one kept; a sector bookmark goes to its system picker
   (`from_sector` / `to_sector`).

### 9.1 Sector Map pick mode

Built in 7.58.0 (`web/sector_page.py`, `_pick_mode`).
`/sector/<id>?pick=to&from=...` shows the same banner. Clicking a
system or phenomenon adds an **End Here** (or **Start Here**) button to
its info panel. With the other end already chosen the button links to
`/nav?from=system:12&to=system:40`, which lands on the plotted course.
With none chosen (a pick begun by the NAV page with one end empty, or by
**Start Here** / **End Here** on a map with no pick on) the button keeps
the user on the map at the same view and zoom: `static/navpick.js` holds
the pick, the URL becomes `?pick=to&from=system:12` in place, the banner
names the end already chosen, and the user zooms out, in and across to
find the other end with the same picker (NAV.29, NAV.33). Cancel on such
a pick clears it in place. The endpoint strings are
`nav_page.endpoint(kind, id)`.

### 9.2 Nav links from everywhere

System and phenomenon pages, and the Sector Map's info panel, get **Nav
from here** and **Nav to here** (the pages), and the maps' panels
**Start Here** and **End Here** (NAV.29, 9.1). The NAV page already
accepts `to` without `from`.

### 9.3 Address bar

A field over the breadcrumb takes a designation, `312/-3/1042` or
`ring 312 layer -3 slot 1042`, `x, y, z` in pc, or a sector or system
name. It flies to the sector level of the slice holding that sector
(section 8.1).

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

- Click-to-center, double-click zoom, the +/− buttons, Slice and Free
  look are gone (MAP.17); the stages are the only way between places.
  Drag-rotate, panning and wheel zoom are back everywhere (section 5.1;
  the whole galaxy since MAP.85). The buttons beside the map are Back,
  Forward, Up, Reset and Bookmarks, with Reset view, Charted only and
  (with polities) Territories in the Menu (MAP.55).
- Density shading, the filled-share look and the info panel stay. The
  wedge lines and the Wedges button are gone (MAP.85), and the scale
  readout is one line, a bar and its length in sectors, pc and ly
  (MAP.60).
- Each stage draws at most a few hundred blocks.
- The Territories overlay (7.54.0, population-and-politics.md) and the
  NAV course line (section 9.4) draw over the map.

## 11. Decisions for Boss

Each has a default, and work can start on it.

1. **Ladder:** 243 → 27 → 3 → 1. *(Decided by Boss, 2026-10-01: "the
   bigger targets".)*
2. **Free camera.** *(Decided by Boss, 2026-10-01. MAP.17 first took
   away all rotation; later that day Boss allowed it again below the top:
   "Once zoomed into an arc or a block, the user can again freely rotate
   and move around the render ... The only time the user cannot freely
   rotate is when at the full galaxy or quarter galaxy zoom levels."
   Built that way: the galaxy and its quarters are locked top-down, every
   step below can be turned, moved and zoomed, and each step opens on its
   own view.)*
3. **Bigger generate buttons.** Default: sector, layer or slab, and
   neighborhood at stages 7-8 only. Option: "Generate this block" at
   stage 5 (up to about 19,000 sectors) behind a confirmation.
4. **Bookmarks.** Default: per browser. Alternative: shared in the
   database (needs a migration).
5. **NAV picking.** Default: a round trip through the map and Sector Map
   pages back to NAV. Alternative: the map embedded in the NAV page.
6. **Orientation of block slabs.** Default: outward up (core at the
   bottom), with a compass arrow. Alternative: bearing 000 to the right
   everywhere, as over the whole galaxy. *(Built as the default, without
   the compass arrow.)*
7. **Neighborhood radius default.** Default: 100 ly, the same as today's
   button. *(Built with this default in 7.53.0. Lowered to 39 ly, about
   12 pc, by GEN.23: the 100 ly around every generated sector now gets
   only its bright stars.)*
8. **What a slice and a region are** (MAP.17, MAP.19). *(Decided by
   Boss, 2026-10-01: "Layer + arc". A slice is a layer of the disk,
   picked from a side strip or list since nothing rotates; a region is
   an arc of the ring band in view, about a third of the current wedge,
   and the map zooms into it. The ladder is quarter, layer, region,
   layer, region, ..., sector, and the slice is the layer. Built as
   section 4 says. Defaults taken: a layer pick is a third of the slabs
   while there are more than three, the region grid is 3 rings by 3
   arcs, hover dims the rest to 25% at once, wedge lines go 15% past the
   view and stop sharply, and the NAV result still fits both ends of
   the course.)*

## 12. Build order and owners

| TODO | Piece | Owner | Depends on |
| --- | --- | --- | --- |
| MAP.28 | Nested ladder geometry, JS and Python, with a parity test (section 3). **Built, 7.41.2, PR #160** | Galaxy Map | none |
| MAP.29 | Stage contents API and caching (section 7). **Built, 7.41.3, PR #160** | Galaxy Map | MAP.28 |
| MAP.16 | The stages: views, hover, pull-out, flight, breadcrumb, URLs, keys, touch (sections 4, 5, 8.1). **Built, 7.44.0, PR #171** | Galaxy Map | MAP.28, MAP.29, MAP.3 |
| MAP.20 | Admin generation at the sector level: `--block` mode, Generate page form, map buttons, radius dialog, progress (section 6). **Map buttons and radius dialog built, 7.53.0, PR #177; `--block` and the form built, 7.60.174** | Web (`planetgen`, Generate page) and Galaxy Map (buttons) | MAP.28 (Python), MAP.16 |
| MAP.21 | Sector Map pick mode and Nav from/to links (sections 9.1, 9.2). **Built, 7.58.0, PR #178** | Web | the URL formats only |
| MAP.22 | NAV page: Pick on Galaxy Map, Pick in this sector, Bookmarks (section 9). **Built** | Web | MAP.16, MAP.21 |
| MAP.23 | Bookmarks (section 8.2). **Built** | Galaxy Map (module, map menu) and Web (NAV, Sector Map) | MAP.16 |
| MAP.24 | Address bar (section 9.3). **Built, 7.50.0, PR #172** | Galaxy Map | MAP.16 |
| MAP.25 | "Show on Galaxy Map" links with `?sector=` (section 8.1). **Built, 7.64.174** | Web | MAP.16's URL format |
| MAP.27 | Course on the Galaxy Map (section 9.4). **Built, 7.52.0, PR #176** | Galaxy Map and Web | MAP.16, MAP.22 |
| MAP.17 | No free camera; drill down top-down by wedge, slice and block (bug against MAP.16; decision 2). **Built, 7.73.177** | Galaxy Map | MAP.16 |
| MAP.26 | "Show on Galaxy Map" opens at the sector level; the map's own Back and Forward (bug; folds in MAP.25). **Built, 7.73.177** | Galaxy Map and Web | MAP.17, MAP.25 |

MAP.3 (the bigger map) shipped in 7.55.0. MAP.27 shipped before MAP.22, using
the NAV result's link rather than the NAV pickers. Everything in this
table is built: MAP.20's `--block` mode in 7.60.174, MAP.25 in 7.64.174,
MAP.17 and MAP.26 in 7.73.177, and MAP.22 and MAP.23 in 7.94.253
(bookmarks per browser, decision 4's default). The arc pick shipped in
PR #369 (section 15); the open work is what follows it, section 15.3
onward.

## 15. The arc pick (MAP.85)

Boss asked for it on 2026-10-01 (23:53Z) and moved it into phase 0 on
2026-10-02 (01:46Z); the arc size below is his decision of 01:53Z.
15.1 and 15.2 are built (PR #369, with MAP.52, MAP.55 and MAP.60) and
describe the map as it is. 15.3 is built too (MAP.56, MAP.54 with
MAP.76, MAP.77), and so are the turn about the middle of each view
(MAP.53) and framing every level whole, fitted to the map's size
(MAP.78). Still planned: the zoom lock below the arc (MAP.58), the ghost
of the other slabs (MAP.59) and the colours of 15.4.

> Redo the galactic selection, so that the galaxy map is 3D, we can
> manipulate it. The user doesn't select an entire wedge, just a large
> arc, then zoom in to select slab, and go from there. ... Don't show any
> sector lines at so that the star map really comes through and the
> spiral pattern. When the user mouses over an arc they can select (each
> arc goes from top to bottom so slab of that we'll do in the next part)
> then they can see it's boundaries and the boundaries of the other
> segments.

"Arc" means this first pick; a "region" is one cell of the 3 by 3 pick
below it (MAP.19, section 4) until MAP.56 replaces that.

### 15.1 The whole galaxy

- The whole galaxy opens as a 3D disk seen straight down, galactic north
  up the screen (MAP.97; it opened tilted 35 degrees before), that the
  user can turn any way (MAP.96). It is drawn as its
  stars, spiral structure and block shading only: no sector, block or
  wedge lines, and no Wedges button.
- Zoom on the whole galaxy is a short manual range: in to about twice
  the fitted view, out until the whole galaxy fits. Locking user zoom
  below the arc, with each pick's staged zoom fitting the view, is
  planned (MAP.58).
- The scale readout is one line: the bar and its length (MAP.60).
- Back, Forward, Up, Reset and Bookmarks stay beside the map; Reset
  view, Charted only and Territories are in a Menu (MAP.55). The
  buttons are text labels and carry a `data-icon` hook, so UX.28's icons
  can replace the words later.

### 15.2 The arc

- An arc is a piece of the disk bounded by bearing and by distance from
  the center, running the full height of the disk.
- Size: 45 degrees of bearing by a third of the disk's radius (inner,
  middle or outer), so the disk has 24 arcs (`ARCS_PER_TURN` = 8 and
  `ARC_BANDS` = 3 in `static/galaxystages.js`). Boss asked for about 40
  degrees (MAP.52's width); 45 is the nearest width whose sides are
  wedge lines of every block ring past the core (their wedge counts are
  all multiples of 8), so an arc's sides run straight from its inner
  edge to its outer one. A block joins the arc its middle bearing falls
  in; the core's three 120-degree blocks do the same.
- Hovering outlines the arc under the cursor along its blocks' own
  sides, top and bottom with its corners joined, and the neighboring
  arcs' outlines faintly. Nothing else is outlined.
- Clicking flies to the arc from the isometric slant, its middle
  bearing up the screen, fitted round its blocks. Its URL token is
  `a<band>.<bearing>` (for example `a1.90`); older links' quarters
  (`q<n>`) still open. Every level is framed whole, fitted to the map's
  actual width and height and centred on its picture, refitted when the
  map is resized and while the view turns (MAP.78); the view turns about
  the middle of what it shows, by drag, a one-finger touch drag, or
  Shift and the arrow keys (MAP.53).

### 15.3 Below the arc

The ladder is arc, slab, segment, then slab and segment again inside
each smaller block, down to a sector (MAP.56, built; the 3 by 3 region
pick is gone, and an older link holding one opens at the stage before
it). A segment's URL token is `s<ring>.<wedge>`.

- **Slab:** a height band of the arc, picked with one button per slab and
  a leader line from each button to its slab (MAP.54, layout MAP.76),
  replacing the slab slider. The arc view draws only the boundaries
  between slabs (MAP.77). Hovering the map lights and outlines the
  whole slab under the pointer, never a single cube (MAP.91).
- **Segment:** one drill block of the next level inside the slab (a
  sector inside a level-3 block), picked on the zoomed slab with the same
  hover highlight; picking it enters that block. The rest of
  the arc stays as a faint ghost so a slab doesn't read as a wedge
  (MAP.59).
- On navigation screens the hover highlight is the same arc highlight
  (NAV.31), and pick mode runs on the shared picker (NAV.15).

### 15.4 Colour instead of lines

With the lines gone, colour shows the structure (Boss, 2026-10-08):
a filled sector is a translucent fill. Its **opacity** follows the star
density (denser is more solid, never fully opaque), its **hue** the mean
stellar age (blue young, slate at the 4.5 Gy disk average, amber old) and
its **brightness** the summed luminosity of its stars on a log scale. A
block takes the same three from its sectors: ages averaged weighted by
star count, luminosities summed, densities as systems per generated
sector. The stage API serves those raw sums (`stats`, DB.14) and the page
works the colour out (MAP.128, MAP.129), ranking each block against the
others in view. Nothing is baked into the database, and the earlier rule
(the average colour of a sector's stars, which came out red nearly
everywhere) is retired (MAP.130). At sector zoom the map shows almost
every star of the sector (MAP.80, phase 1).

### 15.5 Order

| ID | Piece | Phase | Needs |
|---|---|---|---|
| MAP.63, MAP.64 | Shared map helpers (`static/mapcore.js`), one camera and input controller (`static/mapcontrol.js`) | 0 | **Built, PR #351**, with TEST.70, MAP.87, MAP.82 to MAP.84, NAV.30 and MAP.81 |
| MAP.60 | One scale line | 0 | **Built, PR #369** |
| MAP.55 | Buttons into a menu; back, forward, up, reset and bookmark showing | 0 | **Built, PR #369**, text labels with a `data-icon` hook for UX.28 |
| MAP.85 + MAP.52 | The arc pick, one PR | 0 | **Built, PR #369** |
| MAP.86, MAP.80 | Sector colours; most stars at sector zoom | 1 | MAP.85, PERF.11 (MAP.86); MAP.90, MAP.87 (MAP.80) |
| MAP.56, MAP.53, MAP.58, MAP.78, MAP.76, MAP.54, MAP.75, MAP.59, MAP.77 | Slab and segment ladder, rotate and fit, zoom limits, leader lines, mini map, ghost | 2 | MAP.85, in that order, one thread |
| MAP.65, NAV.15, NAV.31 | One picking layer, pick mode on it, the arc highlight while picking | 2 | MAP.56 (MAP.65); NAV.13, NAV.14, MAP.65 (NAV.15); NAV.15 (NAV.31) |
| MAP.66 to MAP.68 | The sector as the last stage, one engine for all three maps | 3 | MAP.65 |

## 13. Sources

- "Mapping Volumetric Galactic Topologies to Planar Interfaces", the
  report Boss shared in the project on 2026-10-01 (not in the repo).
- J. J. van Wijk and W. A. A. Nuij, "Smooth and efficient zooming and
  panning", IEEE Transactions on Visualization and Computer Graphics
  9(4), 2003, pp. 424-431. The source of the formulas in section 5.3.
- This repo: `static/galaxyprisms.js` (`blockWedgeCount`,
  `blockSlotRange`, `blockSectorCount`), `planetgen/galaxy/geometry.py`
  (`ring_master_count`, `ring_sector_count`), `planetgen/generation/limits.py`,
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
  being removed (7.44.0). MAP.17 was Boss's answer, and it is gone.
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

## 16. Planned changes (2026-10-03 and 2026-10-07)

- **Charted only** (MAP.111, done): the "Generated only" toggle was renamed;
  it dims stars outside charted sectors as well as blocks, outlines the
  charted sectors in each wedge, and never blocks picking.
- **Free navigation**: every slab, wedge and block can be picked whether
  or not anything is in it, the slab buttons never cover the pick area,
  and the breadcrumb is drawn from the map's own state on every change.
- **Select mode**: Galaxy (blocks and sectors, the default) or Star, where
  stars and other objects pick, hover and offer Bookmark and Waypoint.
- **Filters**: show, hide or highlight star classes and phenomenon kinds,
  and a luminosity floor, on the Galaxy and Sector Maps.
- **Less on the Galaxy Map**: no comets, rogue planets or asteroid fields
  above the sector level, and a screen-density budget per tile that gives
  sparse sectors more room than dense ones.
- **Streaming**: tiles with level of detail, picking through a BVH, and
  camera-relative rendering (library-migration.md).
- **Opening view**: with nothing in the URL, the map fits all charted
  space.
- **Neighbours on the Sector Map**: the surrounding sectors drawn dimmed,
  fading those between the camera and the sector.
- **Infinite zoom** (the end state): one 3D interface from the galaxy down
  to a moon.
