# Making nebulae visible on the Galaxy Map

Boss (2026-10-10 21:46Z): "I also want nebula should be visible, even if they just change the color of the sector shading on the galactic map, we create them, I want to see them somehow, maybe by color maybe by combined region (if like 5 nebulas are close enough together to appear like 1 on a bigger map we could do that perhaps). ... let's add this to phase 2, with foundations in phase 1."

Informs: MAP.132 (overlay markers), MAP.142 (fuzzy boundaries), MAP.151 (region data layer and its per-level aggregates), MAP.153 to MAP.155 (rank fade, nested lists, other objects fade in), MAP.157 to MAP.161 (wire format), MAP.131 and MAP.143 (Color by), MAP.123 (per-kind toggles), GEN.47 (the molecular cloud field), GEN.99 and GEN.150 (hosted nebulae). Builds on [nebula-and-asteroid-field-classes.md](nebula-and-asteroid-field-classes.md), [zoom-star-visibility.md](zoom-star-visibility.md), [galaxy-map-wire-format.md](galaxy-map-wire-format.md) and [slow-reads-and-timeouts.md](slow-reads-and-timeouts.md).

Status: research, 2026-10-10; nothing here is built. Evidence tags: [S] seen in the repo, [C] computed or measured in this research, [R] recalled and not confirmed (listed at the end).

## Summary

- **Four separate reasons you cannot see them, and each is real.** (1) Most nebulae are not rows anywhere: the molecular-cloud field is *derived* from the galaxy seed and stored only when a generated sector reaches a cloud, and emission and reflection nebulae are created only when a system with an O or B star is built [S]. The default galaxy's field holds about 700,000 clouds [C]; the quarter-scale test galaxy had 7 stored nebulae in 44 sectors [C]. (2) A cloud is a sprite hidden under 2 px (`CLOUD_MIN_PX`) [S]; the field's median cloud is 2.7 pc across in radius, and at the opening view one pixel is 31 pc, so none draws until the camera is within about 12 kpc for the biggest and 4 kpc for the common dark cloud [C]. (3) Dark clouds are painted `#1c1c24` at 0.91 opacity on a transparent canvas [S], which on a dark page is 1.17:1 contrast (MAP.142's own figure). (4) A tile lists at most 200 clouds, largest first [S].
- **Recommended: three scales that hand over to each other, all keyed to pixel size.** Far: a **nebula cover** shading (the share of each block's volume inside a nebula), as a new "Color by" choice. Middle: **regions**, one soft sprite for each group of nebulae that would be under about 16 px apart, splitting into their members as you zoom in. Near: **individual** nebulae, drawn as today (sprite, then shape mesh). The same fixed grid carries all three, so panning never reshuffles a group.
- **Foundations (Phase 1) are data, not drawing:** a thin stored `nebula_field` table made once from the seed (about 700,000 rows, about 50 MB [C by estimate]), per-level cell aggregates in MAP.151's pyramid, the nebula layer in the tile wire format, and a contrast fix for the dark cloud colour. **The visible feature (Phase 2)** is the region layer, the cover shading, individual field clouds, picking and the legend.
- **Grouping by itself does not make the counts small.** Scale-relative cells at 16 px give 144 regions at 30 kpc and 574 at 15 kpc, but 3,238 at 8 kpc and 6,122 at 2 kpc in a dense arm [C]. So regions also need a budget (about 600) and the rank-birth fade MAP.153 already defines for stars; the cover shading carries what is not drawn.
- **Reading all of this from the database is cheap only if it is stored and indexed.** Computing the field on demand costs 0.07 to 0.55 ms per 50 pc cell [C], which is 1 minute for one level-3 tile; so precompute once, at plan time, and read it with a primary-key range, as [slow-reads-and-timeouts.md](slow-reads-and-timeouts.md) asks of every page.
- **One real dependency, not a question:** the field is probably 10 to 40 times too full (the GEN.47 note's own rate check). Cover shading at today's rate is about 11% nearly everywhere in an arm and would carry little signal; at the recommended rate (about 1e-7 pc^-3) clouds become rare and the individual and region layers matter most. The design works at either rate; the rate is Boss's separate decision.
- Fourteen build items, with their phase, are in section 9. Defaults in section 10.

## 1. What nebula data exists

[S] unless marked.

| Source | What | Where it lives | Created when | How many [C] |
|---|---|---|---|---|
| Molecular cloud field, classes M, N, P, Q (dark) | center, radius, class, shape parameters | derived from the galaxy seed, 50 pc cells (`galaxy/nebula_field.py`); a row in `nebulae` once a sector reaches it | stored when the first generated sector whose sphere reaches it is saved (`_db.insert_sector`) | about 700,000 over the default galaxy: 14.6% class M (radius 35 / 81 / 137 ly at p10 / median / p90), 43.8% N (5.5 / 26 / 45 ly), 27.5% P (0.6 / 1.6 / 2.7 ly), 14.1% Q (0.2 / 0.5 / 0.9 ly) |
| Emission and reflection nebulae, classes C to G | around an O star (always), B0 to B2 (half), later B and A (rarely) | `nebulae` row with the star's system | when the system is generated (`add_star_hosted_nebulae`, `NEBULA_HOST_RULES`) | none until built; the pre-placed bright stars hold about 4,200 O and 10,000 B0 to B2 stars at quarter scale, so about 270,000 and 640,000 at default scale |
| Planetary nebulae, classes H to L | around a young white dwarf | `phenomenon_scatter` kind `planetary-nebula`, and `nebulae` once built | scattered galaxy-wide; built with the sector | 445 scatter rows at quarter scale (about 28,000 at default) [R scaling] |
| Supernova remnants, classes R to W | shells | `phenomenon_scatter` kind `supernova-remnant`, `supernova_remnants` once built | same | 306 scatter rows at quarter scale (about 20,000 at default) [R scaling]; `supernova_remnants` empty in the test galaxy |
| Diffuse gas, classes A and B | not generated | | | 0 |

Per record the map receives: `type`, `id`, `name`, `descriptor` (the family), `class`, `radius_pc`, `x/y/z` (`galaxy_clouds_in_box`). A shape (`nebula_shape`: axes, warp, iso level) is served separately per nebula and fetched when the nebula is big enough on screen (`NEBULA_MESH_MIN_PX`). Families and look: diffuse `#e3a6c8`, emission `#ff6f91`, reflection `#6fa8ff`, planetary `#5be8c9`, dark `#1c1c24`, remnants orange with a blue rim (`NEBULA_LOOKS`, `cloudTexture`).

The scatter rows (planetary nebulae, remnants) are not drawn at all today: the tile's scattered points are only black holes, neutron stars and quasars (`SCATTERED_POINT_CLASSES`).

## 2. Why they cannot be seen today

1. **Not stored.** `galaxy_clouds_in_box` reads the `nebulae` and `supernova_remnants` tables. A galaxy that has generated only a small part of itself has almost none of its field there (7 nebulae in the 44-sector test galaxy [C]). The field, the O/B nebulae and the scatter nebulae are all implicit.
2. **Too small to draw.** A sprite shows only when its radius is at least 2 px (`CLOUD_MIN_PX`). At camera radius R the map shows about `0.001 R` parsecs per pixel (900 px high, 50 degree field) [C], so a cloud of radius r appears at about R = 480 r:

   | Cloud | Median radius | First drawn (camera radius) |
   |---|---|---|
   | Class M | 25 pc | about 12 kpc |
   | Class N | 8 pc | about 4 kpc |
   | Class P | 0.5 pc | about 240 pc |
   | Class Q | 0.15 pc | about 70 pc |

   The opening view (34 kpc, 31 pc a pixel) draws none [C].
3. **Colour.** The dark family is the field, which makes it two thirds of all nebulae by number; its fill is near-black on a transparent canvas, so on the dark theme it is invisible, and MAP.142 already counts a 1.17:1 contrast against a 3:1 target.
4. **The cap.** 200 per tile, largest first, then nothing for small ones. In the sample dense arm, 15,049 clouds are 2 px or bigger at a 4 kpc camera and the four view tiles show 800 [C].

## 3. The design: three scales, one grid

A nebula's visibility is a **size law**, not a brightness law: it is drawn when its angular size reaches a threshold, the way a star is drawn when its apparent magnitude reaches one (MAP.148). That gives each scale a clean hand-over.

| Camera radius (pc per pixel) | What is drawn | Why |
|---|---|---|
| Far (at or above about 8 kpc, 8 pc per pixel and up) | **Cover shading** in the block fill, and **region sprites** for the largest groups | individual clouds are under 2 px; the cover says where gas is |
| Middle (about 100 pc to 8 kpc) | **Regions**, splitting into **individual** nebulae as they pass 2 to 4 px | the member nebulae sit within a few pixels of each other |
| Near (under about 100 pc) | **Individual** nebulae: sprite from 2 px (with MAP.155's 1 to 4 px opacity ramp), shape mesh from `NEBULA_MESH_MIN_PX`; the focused one as the 3D field | the existing behaviour, now for field clouds that are not yet stored |

### 3.1 Cover shading (Boss's "change the colour of the sector shading")

The cover of a block is the share of its volume inside some nebula. It is a number per cell of the pyramid MAP.151 builds ("per-level aggregates"), so it can colour every cell the map draws, including **planned and unfilled cells**, because it comes from the seed and not from generated sectors. It is offered as a new choice in MAP.131's **Color by** switch, "Nebula cover", beside Density, Mean age, Luminosity, Star count (and MAP.143's Habitable worlds). That is how it avoids clashing with the existing sector colouring: the default (age hue, density opacity, luminosity brightness) is untouched, and "Nebula cover" replaces it for as long as it is chosen, the way every other statistic does. Sector space on the Sector Map is not tinted, as before.

Ramp and legend: a single hue that survives greyscale, 5 bins (0, under 5%, 5 to 15%, 15 to 30%, 30% and up), the dust brown of section 4 from pale to deep, with the legend line "share of the block inside a nebula". Measured on the sample arm: the cover is 11% median and 17% at p90 per 256 pc block at today's rate, and per 64 pc block 69% hold none and p90 is 40% [C]; so use the 256 pc cell and larger for shading and let the regions and individuals show the small structure. Across the galaxy it follows the gas, which is the arms and the inner disc (expected clouds by ring: 111,000 inside 2.5 kpc, 236,000 at 2.5 to 5, 211,000 at 5 to 7.5, 99,000 at 7.5 to 10, 33,000 at 10 to 12.5, 11,000 beyond) [C].

### 3.2 Regions (Boss's "combined region")

**Rule.** At a given zoom, take the pyramid level whose cell edge is about 16 pixels across (MAP.151's 3-ary levels, so 0.84 to 1.25 of that). Every nebula whose center is in a cell belongs to that cell's group. A group of one is drawn as itself; a group of two or more is drawn as one **region**: a soft sprite at the volume-weighted centroid, radius reaching the outer edge of its members (and never past the cell), coloured by the family with the most volume, with the member count in its tooltip. A region is drawn when it is 2 px across, like any cloud.

**Why a fixed grid and not "friends within 5 px".** A distance-linkage cluster chains: neighbours of neighbours join until a whole arm is one blob, and the grouping depends on the camera, so panning reshuffles it. A fixed grid cell gives the same groups at the same zoom for everyone, nests exactly (a parent's group is the union of its children's), and is the same structure MAP.154 uses for stars. It is also what the cell pyramid stores, so there is no clustering computation at request time.

**Splitting as you zoom in.** The parent region fades out and its child regions (or members) fade in over one octave of zoom, the same cross-fade MAP.153 gives a star between its tile ranks; nothing pops. Birth radius of a group is the camera radius at which its cell edge passes 16 px.

**Measured on the sample arm** (45,061 field clouds in a 4 kpc box 6 kpc from the core; 900 px view) [C]:

| Camera radius (pc) | pc per px | Individual clouds at 2 px or more | After today's 200 per tile | Regions at 16 px (cell edge) | Largest region |
|---|---|---|---|---|---|
| 30,000 | 31 | 0 | 0 | 144 (512 pc) | 14 px |
| 15,000 | 16 | 2,425 | 800 | 574 (256 pc) | 20 px |
| 8,000 | 8.3 | 4,683 | 800 | 3,238 (128 pc) | 22 px |
| 4,000 | 4.1 | 15,049 | 800 | 10,981 (64 pc) | 23 px |
| 2,000 | 2.1 | 20,700 | 1,191 | 6,122 (32 pc) | 29 px |
| 1,000 | 1.0 | 14,404 | 1,600 | 685 (16 pc) + 13,713 single | 47 px |
| 500 | 0.5 | 4,416 | 1,600 | 25 + 4,395 single | 69 px |
| 250 | 0.3 | 1,405 | 1,204 | 1 + 1,405 single | 12 px |

(The 8 px and 32 px versions give about 3 times more and 4 times fewer regions; 16 px is the middle.) The sample is a dense case; the 45,061 clouds are the view sphere's content, more than the frustum shows. Three conclusions:

- A region is never large: 14 to 69 px at 16 px cells, so it reads as a soft blob, not a wall.
- Counts stay in the thousands from 8 kpc to 1 kpc, so **a budget is required**: draw the largest 600 groups by volume, the rest folded into the cover. Use MAP.153's rank-birth rule with the group's volume as the rank, and the 600 is the same kind of number as the 400 stars per tile.
- Below about 250 pc the clouds are single again, so the near layer is the individual one, as today.

### 3.3 Individual nebulae

Unchanged in look. What changes is that they come from the field index, not only from stored `nebulae` rows, so a field cloud is drawn before its sector is generated. A cloud that is stored is the same object (its id is the stored object id; `cloud_object_id` gives a field cloud its id from the seed alone, GEN.176), so the map shows it once, and its page opens when it is picked. Emission, reflection, planetary and remnant objects appear as they are generated; the scatter's planetary nebulae and remnants are listed as points from level 8 like MAP.155's other point objects (section 9 item 11).

## 4. Colour and legend

The sector fill's colours carry age (hue), density (transparency) and luminosity (brightness); Color by modes use a purple, teal, yellow ramp. Nebulae are an **overlay layer of translucent sprites**, so their colours need to differ from the fill's hues and stay legible on top of it:

| Kind | Colour (existing) | Change |
|---|---|---|
| Emission | `#ff6f91` | keep |
| Reflection | `#6fa8ff` | keep |
| Planetary | `#5be8c9` | keep |
| Diffuse | `#e3a6c8` | keep |
| Remnant | orange with a pale blue rim | keep |
| **Dark** | `#1c1c24` at 0.91 | **replace with a dust colour that clears 3:1 on both themes** (proposal: a warm umber, about `#a9805a` at 0.45, with a thin light rim `#d9b48a` on the dark theme; keep the dark fill on the light theme). Checked with the project's contrast test, not by eye. |
| Cover shading | none yet | the same umber, pale to deep, 5 bins |

A region takes its dominant family's colour with a thin outline to separate it from a single nebula. The legend is a small "Nebulae" block under the existing Color by legend with one swatch per kind and a toggle each (MAP.123), plus the line "grouped when close together; zoom in to separate". Colour is never the only cue (MAP.142 asks for a non-colour cue per class): regions get the outline, dark clouds the rim, remnants their ring.

## 5. Data, storage and the wire format

### 5.1 Foundations

- **`nebula_field`, a thin table made once.** One row per field cloud: object id, cell index, center (mpc), radius, class, a flag set when its `nebulae` row exists. Made at plan time from the seed by iterating the 50 pc cells (3.4 million cells in the thin disc, 0.07 ms in sparse cells and 0.55 ms in dense ones, so about 4 to 30 minutes on one core and a few minutes over the RQ workers [C estimate]). About 700,000 rows at about 70 bytes: roughly 50 MB with indexes [C estimate]. It is immutable: the seed fixes it, so there is no refresh and no write contention with a fill. A stored `nebulae` row and its field row share the object id.
- **Cell aggregates in the MAP.151 pyramid.** For every cell from 1 kpc down to about 8 pc: count, summed volume, volume-weighted centroid, bounding radius, volume share per family, dominant class, **cover** (volume of member spheres inside the cell over the cell's volume; overlaps double-counted and capped at 1, which at 11% fill is within a few percent of the true union). From `nebula_field` by one `GROUP BY` per level; the levels from 1 kpc to 64 pc total well under a million rows (about 700 cells at 1 kpc, 11,000 at 256 pc, 170,000 at 64 pc over the thin disc [C estimate]). Stored non-field nebulae (emission, reflection, planetary, remnants) are added to the cells they sit in when stored.
- **The nebula layer in the tile.** Each tile gains `nebulae: {regions: [...], singles: [...]}`, nested the way MAP.154 nests stars (a child omits what its parent already sent); records are tile-relative integers in the packed format of MAP.159 (about 14 bytes a region: position as three 16-bit values, radius and volume as log bytes, member count, family shares as four 4-bit values, class letter), and ordinary JSON until then (about 120 bytes). 600 regions are 72 KB of JSON and about 9 KB packed, against 30 KB for today's 200 clouds. The cover value rides in the stage view's cell stats beside density and luminosity.
- **Reads.** A primary-key range on `(level, ix, iy, iz)` for aggregates and `(cell index, x)` for singles; both examine only the rows returned, and a query-budget test (slow-reads item 11) fails when they do not. No scan of the nebula tables, no sort of a large set.

### 5.2 Why not compute the field on demand

Drawing a level-3 tile (8 kpc across, 300 pc thick) by calling `cell_clouds` over its cells is 307,000 cells, about one minute [C]. The coarse levels are exactly where regions and cover live, so they cannot be derived per request. Precomputing is also the only way the cover can colour unfilled cells.

## 6. Cost on big galaxies

- Build once: 4 to 30 minutes on one core (a few minutes in parallel), about 50 MB of rows and under a million aggregate rows [C estimate].
- Per tile: a handful of primary-key lookups, no scan; at the 600-region budget about 72 KB JSON or 9 KB packed.
- Browser: 600 sprites plus the existing meshes; the cloud sprite path already handles several hundred. The cover shading adds no objects, only a colour on cells.
- A fill does not touch `nebula_field` or the aggregates (they are seed-fixed); a stored nebula added by a fill updates one aggregate cell per level. No effect on the statement limit.
- Cache stamp: the nebula layer changes the tile format, so it bumps the stamp with MAP.147 and MAP.151 (one bump for all of them).

## 7. Picking and hover

- **Individual** (unchanged): `cloudAtRay` picks the smallest visible cloud under the pointer that is at most 150 px across; inside its core (half the radius) a click opens its page, up close clicks go to the sectors inside it.
- **Region:** hover outlines it and shows "5 nebulae, largest Orion-class giant molecular cloud, dark 80% / emission 20%"; click flies to a camera radius at which the region splits (MAP.150's double-click flight), not to a page, because a region has no page. A region under the pointer yields to an individual nebula inside it once that is 2 px.
- **Cover:** hover on a cell shows the percentage (the same hover the other Color by modes use).
- **Touch:** tap shows the tooltip; a second tap flies. **Keyboard:** the legend's per-kind toggles and a "Nebulae" entry in the focus list.
- A region or cloud never blocks picking a sector: the pick order is the existing one (sector first at sector zoom, cloud first at galaxy zoom).

## 8. How it behaves under a fill, and what the visitor sees while it loads

The nebula layer is seed-fixed and read by key, so a fill neither slows it nor stales it (except one aggregate cell when a stored nebula appears). A late tile only delays nebulae that are still at zero opacity under the rank fade (the visitor sees fewer, not different, nebulae); the cover shading is in the stage view, which is already fetched once. If the nebula part of a tile times out, the tile is served without it and refilled (slow-reads section 3.3, per-piece budgets).

## 9. Build items (phase and order)

Numbers are placeholders; the TODO thread assigns IDs.

**Phase 1: foundations**

1. **Nebula field table.** `nebula_field` made at plan time from the seed (immutable, one row per field cloud with its object id and a stored flag); an RQ job over the cell grid; a test that its rows equal `clouds_reaching` for sampled sectors and that stored nebulae match by id. Schema migration with the previous-version fixture.
2. **Nebula cell aggregates.** Per-level count, volume, centroid, bounding radius, family shares, dominant class and cover for MAP.151's cells, from item 1 and from stored nebulae; incremental update when a nebula is stored; a test that cover agrees with a Monte Carlo estimate within 5%.
3. **Nebula layer in the tile format.** `regions` and `singles` per tile, nested like MAP.154, in JSON first and the packed form with MAP.159; the budget (600 regions) and rank (volume); cache stamp bumped with MAP.147 and MAP.151.
4. **Nebula cover in the stage-view cell stats**, next to density and luminosity, for filled and unfilled cells.
5. **Scatter nebulae into the same index.** Planetary nebulae and remnants from `phenomenon_scatter` listed with their stored counterparts, deduplicated by object id; an index for the query (see slow-reads item 4).
6. **Dark nebula contrast (bug).** The dark family's fill is near-black on a dark page; add per-theme tokens and a contrast test (3:1), the part of MAP.142 that does not need the shape work.
7. **Big-galaxy query test for the nebula reads** (rides slow-reads item 11).

**Phase 2: the visible feature**

8. **Region layer.** Fixed-grid groups at 16 px, soft sprite in the dominant family colour with an outline, split by fade (MAP.153 rule) as you zoom, budget 600, per-kind toggles (MAP.123), legend block.
9. **Color by: Nebula cover.** The new choice in MAP.131's switch, single-hue 5-bin ramp and legend, hover percentage, applied to unfilled cells too.
10. **Field clouds as individuals.** Draw nebulae from `nebula_field` as sprites (1 to 4 px ramp, MAP.155) and meshes, deduplicated with stored ones by id; a browser test finds one of each family on a seeded galaxy.
11. **Scatter nebulae as markers** from level 8 (MAP.132's far-zoom nebula markers), planetary nebulae and remnants.
12. **Picking, hover and flight** for regions and cover (section 7); a keyboard path.
13. **Dust colour** for the dark family on both themes, chosen with the contrast checker (replaces item 6's temporary value if needed).
14. **Optional: hosted H II and reflection nebulae for pre-placed O and B stars**, drawn before their sector exists, from the same seeded draw the sector fill will use (GEN.99's conditional redraw keyed by star seed and nebula id) so the two agree. This would show about 590,000 more nebulae at default scale [C scaling]; default is **not to build it** until the O and B star counts and GEN.150's radii are settled.

## 10. Defaults and the one open dependency

- **Store the field, do not derive it per request.** Default: store (section 5.2).
- **16 px cells, budget 600 regions.** Default as stated; tunable constants.
- **Cover shading is a Color by choice, off by default.** Default: off; the regions and individuals are on.
- **Dust colour for the dark family** replaces the near-black. Default: yes (it is a bug).
- **Hosted H II for pre-placed stars: no** (item 14) until the star counts and radii are settled.
- **The field's rate** (GEN.47 note: 5e-6 pc^-3 is 10 to 40 times the catalogues) is a separate decision of Boss's and stays one. The design needs nothing from it; at a lower rate the cover shading shows sparse gas and the region layer shows hundreds, not tens of thousands, of groups.

## Limits of this study

- Counts are from the in-memory field model on the default shape (`nebula_field.cell_clouds`, `gas_factor`), 250 pc by 25 pc integration grid for the totals and one 4 kpc sample box for grouping; arms and the inner disc differ from the sample.
- The sample view counts the clouds inside a sphere of 1.6 times the camera radius, which is larger than the frustum; real counts in view are lower, so the budget bites later than the table shows.
- Pixel sizes assume a 900 px high canvas and the 50 degree field of the page [S]; a phone is smaller and meets thresholds at higher camera radii.
- Row and byte sizes of the new tables are estimates from column widths, not a built table.
- The O and B and scatter counts are quarter-scale test-galaxy figures scaled by the 64-times volume the generation-performance study gives [R].
- The colour proposal is not tested on the real themes; the contrast checker settles it.

## Evidence notes

- [S] `galaxy/nebula_field.py`, `db/query.py` (`galaxy_clouds_in_box`, `GALAXY_TILE_MAX_CLOUDS`, `SCATTERED_POINT_CLASSES`), `static/galaxymap3d.js` (`NEBULA_LOOKS`, `cloudTexture`, `CLOUD_MIN_PX`, `cloudAtRay`, `CLOUD_PICK_MAX_PX`, `renderer.setClearColor`), `static/galaxyblocks.js` (`COLOR_MODES`, `STAT_RAMP`), `tuning.py` (`NEBULA_CLASSES`, `NEBULA_HOST_RULES`, `PHENOMENON_DENSITY_PC3`), TODO items MAP.131, MAP.132, MAP.142, MAP.143, MAP.151, MAP.153 to MAP.155, MAP.171.
- [C] scripts and outputs in `/mnt/project-files/research/nebula-map-visibility/` (`field_counts.py`, `cluster_by_zoom.py`).
- [R] the 64-times scaling of quarter-scale counts; the dust colour values; Kennicutt and Evans and the catalogue counts are from the GEN.47 note, not re-checked.
