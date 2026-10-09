# Fly-through Galaxy Map: view distance, see-through near field, free camera (MAP.146 rebuilt)

Research and design, 2026-10-09, from Boss's message of 22:56Z. No code changes. The rewritten MAP.146 text and the proposed split are in sections 7 and 8, for the TODO thread to file.

Boss's goal: one continuous 3D world from the whole galaxy to a star system. Scroll wheel zooms, double-click flies to a clickable thing, the observer can be inside a block or sector and sees its contents without looking around what is in the way, the region being looked at is rendered well and the surroundings are visible but quiet, and only stars an observer could see are drawn.

## 1. Answer in brief

1. **Stars: one smooth law, not tile steps.** Draw a star at opacity `a = smoothstep((m_lim - m) / 1.5)`, where `m = M + 5 log10(d / 10 pc)` is its apparent magnitude from the camera and `m_lim` is chosen each frame so about N stars (default 20,000) are drawn. The inverse-square "floor luminosity at distance d" is the same law written the other way: `L_min(d) = L_min(d0) (d/d0)^2`, a straight line of slope 2 on a log-log plot. Today's rule is already this law, but frozen into steps (section 2.2).
2. **Near field: three cheap, standard tools together.** (a) A depth fade, so anything closer than a fraction of the focus distance dissolves. (b) A soft-edged see-through tube from the camera to the focus ("ghosting"), so what is in the way thins out as it gets in the way. (c) The container the camera is inside is drawn from the inside, as faint edges, with its contents at full strength. These replace the one-off blocker formula in MAP.121.
3. **What is visible is what is clickable.** The picker ignores anything whose current opacity is under 0.35 and anything nearer than the fade's start distance, using the same function the shader uses.
4. **Camera: free, with the observer able to be inside.** Wheel zooms toward the point under the cursor, speed proportional to the distance to the nearest object (so it can never overshoot into one). Double-click flies to the thing clicked on a van Wijk and Nuij path. The stage-by-stage presets and the arc, slab and segment picks stop being the way you navigate.
5. **Data under the camera: a distance-cut of the existing octree tiles plus MAP.146's aligned cells.** Nearer tiles are finer. A 60 degree forward view needs about 140 to 280 tiles instead of today's 27 (section 4.5; empty tiles not counted).
6. **Focus and context are one number.** Each region gets a prominence from 0 to 1 (the thing looked at, the container, the neighbours). Prominence scales opacity and lowers `m_lim` locally (context shows only stars 2.5 magnitudes brighter), so the sector in view is rich and its surroundings are present but cheap.
7. **MAP.146's earlier recommendation stays as its data layer:** exact-centred frame on aligned cells. For continuous flight, prefer the 3-ary cell pyramid (section 6).

## 2. What the map does today (checked in the repo)

### 2.1 Camera and stages
- `galaxystageview.js`: stage presets (galaxy, arc, slab, segment, block, sector, system). Each flies to a preset camera. Zoom runs from `MIN_ZOOM` 1/8 to `MAX_ZOOM` 2.5 times the stage's fit; the galaxy only from fit to twice as close (`GALAXY_MIN_ZOOM` 0.5). A system opened in place goes to `SYSTEM_MIN_ZOOM` 1e-9. A wheel move from one stage to the next takes one gesture with a 900 ms cool-down (`ACROSS_COOLDOWN_MS`, MAP.125).
- `galaxymap3d.js`: an orbit camera: a target and a distance back along a quaternion. The camera is always outside the thing shown. The near plane follows `min(BASE_NEAR, dist * 1e-3)` and the depth buffer is logarithmic, so a close near plane costs no precision.
- Wheel factor is `0.0025` per pixel of wheel travel (`WHEEL_ZOOM_PER_PX`), so a 100 px notch scales the distance by about 0.78.

### 2.2 Stars: the rule is a limiting magnitude in steps
`queryDb.generated_star_floor_sol(level)` lists in a tile only generated stars at least `0.004 L_sun * (tile_edge / 32 pc)^2`. That is a limiting magnitude of **13.35 at a distance equal to the tile edge**, at every level (computed: level 11, 32 pc, floor 0.004 L_sun, M 10.82; level 7, 512 pc, floor 1.02 L_sun, M 4.80; level 3, 8.2 kpc, floor 262 L_sun, M -1.22; each is 13.35 at its own edge). The floor quadruples (1.5 magnitudes) per level, and a view uses one level, so the limit steps by 1.5 magnitudes at each zoom level and is measured against the tile, not against the camera. Boss's gradient is this law made continuous and tied to the observer.
- Brightness on screen: `starBright` comes from log luminosity only (`STAR_LOG_LUMINOSITY` -4 to 6), with the Sector Map's boost (`starlight.js`, MAP.87). It does not change with distance: a supergiant 10 kpc away and an identical one 10 pc away draw the same.
- Stars are drawn with `depthTest: false`, at a fixed 6 to 30 pixels (`STAR_MIN_PX`, `STAR_MAX_PX`), so they are never hidden by blocks and never occluded by anything, including something in front of the camera.
- Budgets: 27 tiles per view of one level (`GALAXY_VIEW_MAX_TILES`), plus up to 8 finest tiles around a target; 150 to 850 generated stars per tile (`GALAXY_TILE_STAR_BUDGET`), 4,000 in a finest tile; 70,000 stars per view overall (`GALAXY_VIEW_MAX_STARS`, MAP.109).

### 2.3 What is in front of the camera
- Block fills use a `fade` uniform and `depthWrite: !translucent`; there is no near fade on blocks. Clouds do have one (`CLOUD_NEAR_FADE` [1.2, 2.0] in cloud radii, `smoothstep`), and CLOUD_PICK_MAX_PX 150 stops a cloud being picked once it fills the screen.
- Picking (`mappick.js`): layers in priority order; a mesh layer raycasts and takes the nearest hit with a meaningful entry; points are picked on screen within a pixel reach. A layer marked `occludes` blocks what is behind it. Nothing checks how visible a hit is, and a ray starting inside a container hits that container's own faces. That is the "clickable what is right in front of the camera, and the object the camera is in" problem.
- MAP.121 and MAP.141 already specify a blocker fade and faint neighbours (design doc `map-ui-and-frontend-libraries.md` section 5.4): neighbour alpha `a0 exp(-d / lambda)`, floor 0.08, ceiling 0.3, and a fade by the distance of a neighbour's centre from the camera-to-focus line. Good, but centre-based, per region, and not a general near field.

## 3. Algorithms researched

| Problem | Technique | Source / basis | Fit here |
|---|---|---|---|
| Which stars to show | Limiting magnitude from apparent magnitude (Pogson: 5 magnitudes = a factor of 100 in flux); `d_vis = 10^((m_lim - M)/5 + 1)` pc | standard photometry | Replaces tile steps. Table in section 4.1 |
| How many to show | Octree LOD where a node is drawn when it subtends enough solid angle at the camera (threshold theta), with a cap on stars held in memory | [Gaia Sky level-of-detail catalogs](https://gaia.ari.uni-heidelberg.de/gaiasky/docs/3.4.2/LOD-catalogs.html) and [data streaming](https://gaia.ari.uni-heidelberg.de/gaiasky/docs/3.0.1/Data-streaming.html) | Same shape as our octree tiles; section 4.5 |
| Dissolve what is close | Depth fade (smoothstep on view depth) | common practice; our clouds already do it | Section 5.1 |
| See into a focus without moving | Context-preserving rendering: reduce opacity where it hides the region of interest, using distance to the eye and accumulated opacity; "ghosting" (keep an occluder faintly) versus "cutaway" (remove it) | [Bruckner et al., Illustrative context-preserving volume rendering (2005)](https://www.cg.tuwien.ac.at/research/publications/2005/bruckner-2005-ICV/bruckner-2005-ICV-Paper.pdf) and the [2006 follow-up on ghosting and cutaways](https://www.cg.tuwien.ac.at/research/publications/2006/bruckner-2006-ICE/bruckner-2006-ICE-Paper.pdf); the Feiner and Seligmann cutaway paper is the usual origin, but I could not retrieve it | Ghosting is better than cutaway for spatial sense; section 5.2 |
| Order of translucent layers | Sorted back to front with `depthWrite: false` (few layers); weighted blended OIT (additive weighted average, weight falling with depth, needs two render targets and float textures; tends to look too transparent); stochastic or dithered (screen-door) transparency (order independent, noisy) | [McGuire and Bavoil, JCGT 2013](https://www.jcgt.org/published/0002/02/09/paper.pdf); [Enderton et al., Stochastic Transparency](https://www.cse.chalmers.se/~d00sint/StochasticTransparency_I3D2010.pdf) | Sorted for the few hundred translucent prisms we draw; dithered discard for the near fade; weighted blended only if overdraw artifacts show |
| Zoom toward a point without overshoot | Point-of-interest movement: approach the target at logarithmic speed, so each step covers a fixed fraction of the remaining distance | [Mackinlay, Card and Robertson, SIGGRAPH 1990](https://www.cs.ubc.ca/~tmm/courses/533-09/slides/nav.pdf) (summarised in these course slides; I did not read the original) | Wheel and double-click approach |
| Fly to a target smoothly across scales | Van Wijk and Nuij's pan-and-zoom path: second-order smooth, minimises perceived optical flow, takes the same time to cross a ratio of scales whatever the distance | van Wijk and Nuij, InfoVis 2003 (not retrieved; described in [Reach and North, Smooth, Efficient, and Interruptible Zooming and Panning](https://arxiv.org/pdf/1801.09358), which extends it) | Already planned in `galaxy-drilldown-navigation.md` 5.3 |
| Precision from galaxy to AU | Camera-relative rendering (positions minus the camera in double precision on the CPU, per-tile origins) and a logarithmic depth buffer | in place for the galaxy (MAP.102); section 4.4 for the numbers | Needed below about 10 pc |

Not verified and so not relied on: the specific limiting-magnitude or crowding rules of Celestia, Stellarium and Space Engine.

## 4. Design

### 4.1 The visibility law for stars (A in Boss's message)
- **Per star:** `m = M + 5 log10(d / 10)` with `M = 4.83 - 2.5 log10(L / L_sun)` and `d` the camera distance in parsecs. Opacity `a = smoothstep(0, 1, (m_lim - m) / 1.5)`. 1.5 magnitudes is the factor 4 that today's floors step by, so the ramp replaces each step by a ramp of the same width.
- **Brightness:** draw the point's strength from apparent flux, not luminosity. A star's glow and core alpha follow `10^(-0.4 (m - m_ref))` through the same tone curve the Sector Map uses (`starlight.js`), so the two maps agree and distance dims a star as it should. Close stars are tone-mapped (`F / (F + F_half)`) so one star next to the camera cannot flood the screen; their disc gets its true angular size once it is more than 1 pixel across.
- **Choosing `m_lim`:** each frame, build a histogram (64 bins) of `m` over the stars in the loaded tiles (at most 70,000, so one pass) and take the magnitude at which the count reaches `N_target`. Defaults: 20,000 on desktop, 8,000 on a phone (unmeasured starting values). This needs no model of the star population because it uses the stars actually loaded, and it automatically shows fewer faint stars in a packed core and more in thin space.
- **Distance a star can be seen** (`d_vis = 10^((m_lim - M)/5 + 1)` pc):

| star | M | m_lim 6.5 (naked eye) | 9 | 12 | 15 |
|---|---|---|---|---|---|
| red dwarf | 12 | 0.8 pc | 2.5 pc | 10 pc | 40 pc |
| Sun | 4.83 | 22 pc | 68 pc | 272 pc | 1.1 kpc |
| A star | 1 | 126 pc | 398 pc | 1.6 kpc | 6.3 kpc |
| giant | -1 | 316 pc | 1.0 kpc | 4.0 kpc | 16 kpc |
| bright giant | -5 | 2.0 kpc | 6.3 kpc | 25 kpc | 100 kpc |
| supergiant | -8 | 7.9 kpc | 25 kpc | 100 kpc | 400 kpc |

  Supergiants stay visible across the galaxy (it is 15 kpc in radius), which is right; a Sun-like star disappears past 1 kpc even at 15.
- **Server side, nothing to change at first.** The tile floors stay as a superset (a tile never lists less than any observer at a sensible distance would want). Once tiles are chosen by distance (4.5) the floor for a tile follows its distance and the two rules coincide.

### 4.2 Prominence: render what is looked at well (B, last paragraph)
Every region (cell, block, sector) has a prominence `p` in 0..1:
- **Looked at** (`p = 1`): contains the focus target, or its centre is within about 10 degrees of the view axis (`w = 1 - smoothstep(10, 25 degrees, angle)`).
- **Container the camera is inside** (`p = 0.8`).
- **Neighbours and the rest** (`p = clamp(exp(-d / lambda), 0.08, 0.3)`, the formula already in design doc 5.4), with `lambda` about one region width.

`p` scales a region's fill opacity and shifts its stars' limit: `m_lim_eff = m_lim - 2.5 (1 - p)` (a factor of 10 in flux). A context sector shows only its bright stars, with no per-star labels or hit tests, so it is cheap. The sector in view shows everything down to `m_lim`.

### 4.3 The camera (C)
- **Wheel:** `dist <- dist * exp(-0.0025 * deltaY)` as now, but toward the point under the cursor, which stays fixed on screen. That point is, in order: the picked star, remnant or cell under the cursor; else the first sample along the cursor ray where the model's integrated column density passes a threshold (the page already evaluates the density in `galaxyprisms.js`); else the point at the current focus depth.
- **Clearance:** each step moves by `min(requested, 0.5 * (d_nearest - d_clear))` where `d_nearest` is the distance to the nearest pickable thing, so speed is logarithmic in the distance to it (Mackinlay) and the camera never lands on a star or inside a wall by accident. `d_clear` is a small multiple of the thing's drawn radius.
- **Double-click:** target is the thing under the pointer (a star, sector, cell or system). Fly on the van Wijk and Nuij path, with the end distance chosen so the target fills about 60% of the view height. This is MAP.140's go-to, extended: a double-click on a block or cell flies in; a double-click on empty space flies a fixed fraction (half) of the way toward the point under the cursor. Offer the same flights from a "Go to" button and Enter, as MAP.140's notes require.
- **Observer inside:** "I am in sector X" comes from the camera position (`geometry.sector_address_at` is closed form), shown in the header and the breadcrumb, which becomes a trail of containers (galaxy, ..., block, sector, system) derived from the position, not from state. Clicking a crumb flies out to it.
- **URLs and bookmarks** hold the camera (target, distance, orientation). Old block URLs and stage tokens stop working (no backward compatibility, Boss's rule).
- **What goes:** the arc, slab and segment picks as the way to move, their presets and the 900 ms stage-hop cool-down. The slab strip can stay as an optional section-plane tool (clip everything above a layer) for looking at one slice; that is a view tool, not a stage.

### 4.4 Scale and precision
Float32 positions at distance D from the origin resolve `D / 2^24`: 0.003 pc (600 AU) at 49 kpc, 12 AU at 1 kpc, 0.12 AU at 10 pc, 0.012 AU at 1 pc. Galaxy-scale positions are fine for stars as points; anything closer than about 100 pc needs positions relative to the camera, computed in double precision per tile (a tile origin as a double, vertex data as float32 offsets, and a per-tile uniform `tileOrigin - camera`). The dynamic range from a 4 pc sector to 1 AU is 825,000, which a log depth buffer and a near plane of `1e-3 * distance` (both in place) cover. A system's disc becomes drawable when it subtends about 20 pixels: for a 50 AU system at 1000 px and 50 degrees of view that is about 3,000 AU (0.014 pc).

### 4.5 Which data the camera pulls (distance-cut)
Today one level covers the whole view (27 tiles). For a camera inside the scene, refine the octree by distance: split a node while its edge is more than `k` times its distance from the camera (and it meets the galaxy and the view cone). Counts below are upper bounds (empty tiles are not removed), for a 60 degree forward view:

| camera | k = 2 | k = 1 | k = 0.5 |
|---|---|---|---|
| Sun radius, in the plane | 128 | 140 | 452 |
| 2 kpc out, 100 pc above the plane | 112 | 272 | 634 |
| bulge, ring 100 | 128 | 204 | 560 |
| galaxy edge, 14 kpc | 128 | 280 | 664 |

All-round free look at `k = 1`: 684 to 936 tiles. At roughly 300 stars per tile, `k = 1` forward is 40,000 to 85,000 stars, so per-tile budgets must shrink with the number of tiles to stay under the 70,000-star cap; the `m_lim` histogram then trims it on screen. A request holds at most 128 tiles, so this is 2 to 8 requests, cached by key. A node's floor luminosity (edge squared) is the magnitude law at constant angular size, so choosing levels by distance and the magnitude law agree. Phase this: first the client-side law on the tiles already fetched, then the distance-cut.

### 4.6 The near field (B, middle paragraphs)
1. **Depth fade (always on).** Multiply a translucent fill's or an outline's alpha by `smoothstep(z0, z1, viewDepth)` with `z0 = 0.15 D` and `z1 = 0.5 D`, where `D` is the camera-to-focus distance (so the hidden zone scales with how close the user is looking). A camera inside a block loses the block's near faces by the same rule. The existing cloud fade `[1.2, 2.0]` radii is the precedent; for the 3-ary pyramid, apply it per cell at its own size.
2. **Focus tube (ghosting).** With camera `C`, focus `F` and fragment `P`: `t = ((P - C) . (F - C)) / |F - C|^2` and `rho = |P - C - t (F - C)|`. For `0 < t < 1` (in front of the focus) multiply alpha by `1 - (1 - a_min) * (1 - smoothstep(R_t, 1.6 R_t, rho))` with `R_t = max(r_focus, 0.12 D)` and `a_min = 0.08`. One dot product and one smoothstep per fragment; two uniforms. The tube can widen toward the camera into a cone for a cutaway feel. This reduces to the MAP.121 blocker rule for a whole-region fade and replaces it. Use ghosting, not removal: context stays faintly visible.
3. **Inside a container.** Draw its faces as back faces at low alpha (or edges only) and everything in it at full strength; do not raycast its own walls.
4. **Order.** Translucent prisms are sorted back to front with `depthWrite: false` (a few hundred; MAP.121 already says this). Dithered discard (a 4 x 4 or blue-noise threshold on the fade alpha) makes the near fade independent of order. If overlapping translucent fills still look wrong, add weighted blended OIT for the fills only, since stars are additive anyway.
5. **Picking agrees with drawing.** Compute `visibleAlpha(point) = fillAlpha * nearFade * tubeFade * prominence` on the CPU with the same constants; a hit counts only if it is at least 0.35, and nothing within `z0` is hit. Stars use the screen-space picker and are skipped if their opacity `a` is under 0.35.

### 4.7 Hand-off between scales
Cross-fade, with hysteresis so nothing flickers: a sector's own scene (all stars) fades in when the camera is within 1.5 edges of its centre and is full by 1 edge; it unloads past 3 edges. A system's disc fades in at about 20 pixels (section 4.4) and unloads at 10. Blocks fade out through the near fade as the camera enters them. One scene graph throughout (MAP.121's one engine).

## 5. Cost (what I can say without a browser)
- Shader: near fade, tube and magnitude ramp are a few multiply-adds and two `smoothstep`s per fragment or vertex; the tube needs three uniforms. Expected negligible next to overdraw. **Not measured.**
- Histogram for `m_lim`: one pass over up to 70,000 stars per frame. **Not measured**; if it is too slow, run it every few frames or on camera move only.
- Tiles: section 4.5 (upper bounds, 140 to 280 forward at `k = 1`, instead of 27).
- Wire: 52 bytes per star today; MAP.147 is measuring the format, and 100,000 stars at 52 bytes is 5.2 MB, so per-tile budgets must drop as tiles multiply.

## 6. Cells and regions (MAP.146's original content, kept as the data layer)
From `research/drilldown-region-sizes/report.md`: frame exactly centred on the click (the ADM.29/30 round region as layer, ring and wrap-aware slot ranges), data from the aligned cells covering it, read from the existing per-block cache, size as cells across. For continuous flight prefer the **3-ary pyramid** (1, 3, 9, 27, 81, 243, 729): a level switch changes cell size by 3 instead of 9, which a cross-fade hides, and the view can hold 5 to 15 cells across (up to 3,375 prisms) without a gap. Cost: cells are less square (centreline arc 0.84 to 1.25 of an edge against 0.94 to 1.06 at 9-ary). Per-level aggregates (counts, age, luminosity) must exist for every level so colour and density need no scan; with a moving camera the stage query that scans sector rows cannot serve. Ratio 2 was not measured.

## 7. MAP.146 rebuilt (text for the TODO thread)

Title: **MAP.146 Fly through the galaxy: scroll-zoom, double-click flight, distance-based visibility and a see-through near field**

Boss (2026-10-09 22:24Z): "I want the zoom drill down to be less specific, instead of set wedges, blocks, and slabs pre-determined, have them based on the center of where the cursor is clicked. So that we always are drilling down exactly where the user wants. We'll have to convert block/slab measures to ranges of layers/shells/slots for filling on demand." Boss (22:49Z and 22:56Z): the end goal is to scroll-zoom and fly smoothly through the interface to a location; the system does not do view distance well and renders and keeps clickable what is right in front of the camera and the object the camera is in; stars should be shown by a smooth distance-and-brightness gradient so only what the observer could see is drawn; inside a block or sector the observer should see its contents without looking around what is in the way, the nearer things becoming more transparent as they approach; double-click zooms or flies to a clickable thing, the wheel zooms in and out, the user is always in the 3D galaxy even when viewing a sector; the sectors looked at are rendered well and the surroundings stay visible but quiet; smooth from the galaxy down to a star system.

Done: the Galaxy Map has one free camera from the whole galaxy to a star system. The wheel zooms toward the point under the cursor and double-click flies to the thing clicked. A star is drawn at an opacity that follows its apparent magnitude from the camera, against a limit chosen to hold about 20,000 stars on screen, so only what an observer could see is drawn. Things nearer the camera than a fraction of the focus distance dissolve, a soft see-through tube thins what stands between the camera and the focus, and the container the camera is in is drawn from the inside. Only what is visible enough can be picked. The region looked at, and the container, are rendered at full strength and the rest faintly. The camera position names the container (the breadcrumb is derived from it) and is what URLs and bookmarks hold. Region data is read in aligned cells and described as layer, ring and slot-arc ranges for stats, fills and backfills. The arc, slab and segment picks stop being the way to move (a slab strip can remain as an optional section plane). Old stage URLs stop working.

Split proposed (the TODO thread assigns IDs and orders them by the dependencies below):
1. The star visibility law: apparent-magnitude opacity, flux-based brightness and the on-screen limit from a histogram (client side first), then tiles chosen by distance. Replaces the step floors in `generated_star_floor_sol` as the rule and folds MAP.116's budget table into it.
2. The near field: depth fade, focus tube, inside-a-container drawing, prominence, and pick agreement. One shared function; folds MAP.121's blocker fade and MAP.141's faint context into it.
3. The free camera: wheel zoom to the cursor, clearance, double-click flight (MAP.140's go-to extended), observer inside with the container named from position, camera URLs and bookmarks.
4. The data layer: exact-centred frame, aligned cells, per-level aggregates, ranges with slot wrap for ADM.29/30 and MAP.120, tile and stage cache keys (decide with MAP.147's wire format).
5. Scale hand-offs: galaxy, sector and system cross-fade with hysteresis, per-tile camera-relative origins below about 100 pc.

Order: 2 and 1 (client side) can start now and improve today's map; 3 needs 2; 4 and 5 follow. Related: MAP.121, MAP.122, MAP.125, MAP.131, MAP.134, MAP.140, MAP.141, MAP.147, ADM.29, ADM.30, GEN.101, GEN.126, NAV.13, NAV.14, MAP.59, MAP.75. Design: this report, and `docs/design/galaxy-drilldown-navigation.md` (to be marked superseded in part).

Open questions for Boss, with defaults:
- Stars on screen: default 20,000 (8,000 on a phone) with a 1.5 magnitude ramp; tune after a first build.
- Cell pyramid: default 3-ary, accepting cells 0.84 to 1.25 of an edge across; or keep 9-ary and accept a 2.45x gap in sizes.
- Context strength: default context regions at opacity 0.08 to 0.3 and 2.5 magnitudes shallower.
- Slab strip: default kept as an optional section plane, not a stage.
- Retiring the arc, slab and segment picks (MAP.85, MAP.56, MAP.17 flow): default yes.

## 8. Phases
1. Client-only, no schema change: items 1 (law on fetched tiles) and 2 (near field, picking agreement). Improves the map as it is.
2. Free camera and the observer inside (3), with the container named from position.
3. Distance-cut tiles and scale hand-offs (1b, 5); aligned data layer and per-level aggregates (4) alongside.

## 8a. Dependencies on other research (added 2026-10-09 23:10Z)

- **Research Lane 1, smoother star counts while zooming** (`docs/design/zoom-star-visibility.md`, merged #878). Its rank birth radius `R_b` gives each star an opacity that depends on the camera radius only; it is reversible and independent of the network. My magnitude law (4.1) depends on the observer's distance, so the two are one product per star: `a = a_rank(R) * a_mag(d) * a_near`. Lane 1's rule is the zoom-driven stage, and mine adds the free camera. Where both set a limit, the dimmer of the two wins, so there is no extra star budget. Build order: Lane 1 stage 1 first (client-only, one shader), then my law as a second factor in the same shader. Its stage 2 (server lists nested, one key at every level) is also what my distance-cut tiles need, since a child tile must contain the parent's stars in its box. Its tile keys are written against "a list of stars per tile, most luminous first", so they survive the MAP.146 change of tile shape. Its decisions for Boss (fade rate W, dense sector stars only at sector zoom, point objects from level 8, backfill shells) overlap my open questions on the star count and ramp; Boss should answer them together.
- **Lane 1's two cautions on combining the rules** (message 23:04Z, agreed): (1) Calibrate `m_lim` so `a_mag` is about 1 for a star at the target distance at any camera radius, and let it only dim stars much farther than the target and ghost the near field. Otherwise the two rules thin the same stars twice and the on-screen count falls below what the rank rule intends (Lane 1's scheme B alone left worst steps of +390% to +777% because tile caps, not a floor, cut the lists). (2) If `R` is redefined for a free camera (for example, distance to the nearest sector instead of to the target), keep it continuous in the camera position. `R` is the only input that makes the picture reversible, and a jump in `R` is a pop for every star at once.
- **Research Lane 3, MAP.147 wire format** (asked 23:08Z, no reply yet). My design adds to each tile a stable per-star id and, for the 3-ary pyramid, per-level aggregates (count, summed luminosity, mean age). It also changes the tile keys (MAP.146 item 1b). The payload numbers MAP.147 measures hold for the tile edges 16 pc to 65,536 pc, but any format chosen has to carry the id and the aggregate. Item 5 (aggregates) should wait for MAP.147.
- **Object IDs** (`research/object-ids/object-id-options.md`): the per-star id I assume is the one that note decides.

## 9. What I could not measure
- Frame time and shader cost (no browser or GPU here), the cost of the histogram, and how it feels.
- Tile counts ignore empty tiles; the stars per tile (about 300) is the existing budget order, not a count from a filled galaxy. 20,000 and 8,000 are starting values.
- Whether the 0.35 visibility pick threshold, the tube radius and the fade distances feel right. All constants are first guesses to tune.
- The original Mackinlay, Card and Robertson, Feiner and Seligmann, and van Wijk and Nuij texts: I used the course slides, abstracts and a follow-up paper for the first and last, and did not retrieve the second. Other systems' (Celestia, Stellarium, Space Engine) rules were not verified.
- How far camera-relative rendering already goes below the sector scale (MAP.102 marks it done; I read the star path, not the system stage).
- Ratio-2 cell squareness.

Scripts for the numbers: `research/drilldown-region-sizes/scripts/` and `cut.py` (tile counts) in this folder.
