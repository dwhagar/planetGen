# Smoother star and object visibility while zooming the Galaxy Map

How stars (and the other objects drawn as points) come and go as the camera zooms, why it happens in bursts today, and a way to make it a steady trickle. Boss asked (2026-10-09 22:30Z): "research different ways to step up and down the # of visible stars smoother than what we have now, I want a very smooth transition where stars and objects are slowly added as one zooms in."

Informs: MAP.147 (the wire format; nothing here depends on its answer), MAP.146 (it changes the tile keys; the rule below is written against "a list of stars per tile, most luminous first", whatever the tile shapes are). Builds on MAP.48, MAP.51, MAP.80, MAP.87, MAP.109, MAP.116 and MAP.123 (all done).

Status: research, 2026-10-09; the recommendation is mine. **Superseded in part by MAP.152 (2026-10-15):** tiles are now chosen by distance and a star's fade comes from its tile's share, not from a birth radius by rank (see `fly-through-view-distance.md`, "The fade under the cut"); the rank rule below describes what MAP.153 built before that. Stage 1 (section 5, MAP.153) was built: `static/starfade.js` holds the arithmetic, `galaxymap3d.js` ranks each tile's stars and the star shader applies the fade; stage 2 (section 8, MAP.154, extended by MAP.250) and stage 3 (MAP.155) are built too, see their "Built" and "As built" paragraphs below. Evidence tags: [S] seen in the repo's code or data, [C] computed in this research, [R] recalled and not confirmed (listed at the end).

## Summary

- **Today the picture changes at 11 camera radii and at no others.** Between them the set of drawn stars is frozen. On the scratch galaxy, 78 of 84 consecutive zoom steps of 9% changed nothing in one spot and 75 of 84 in another; at the others a burst of new stars appeared at once, up to 8.8 times the number on screen before (table 2) [C]. The burst then fades in over 300 ms, a time fade, not a zoom fade, and a star that leaves (zooming out) vanishes without any fade [S].
- **The cause is the tile ladder.** The map picks one tile level from the camera radius (`tileLevelForRadius`); each level has its own luminosity floor (four times fainter per level), its own cap per tile and its own list of point objects, so a new level swaps the whole star set in one frame. A second step is the "detail" tiles: every star within 8 pc of the target joins at once when the view radius drops under 160 pc [S].
- **Recommended: give every star a birth radius and draw its opacity as a smooth function of the camera radius.** The birth radius comes from the star's rank in its tile's list, which the server already sorts most luminous first, so the first stage needs no server or schema change. Across the tile levels the opacity cross-fades from the parent tile's rank to the star's own. On the same tile data the worst single step falls from +883% to +59% (dense spot) and from +775% to +74% (thin spot), the typical step is 3 to 4%, and nothing drops by more than 6% [C].
- **It is symmetric.** The opacity depends only on the camera radius, so zooming out removes stars as smoothly as zooming in adds them, a bookmarked view always draws the same picture, and the network no longer decides when stars appear (a late tile only delays stars that are still invisible).
- **Other objects** (black holes, neutron stars, quasars, nebula sprites, sector blocks) have their own hard switches. The same rule covers the point objects; clouds need a pixel-size ramp; blocks are a separate level-of-detail question (section 6).
- **Decisions for Boss** are in section 7: how fast stars should come in, whether a dense sector may show all its stars only at sector zoom, and whether the 100 ly backfill tiers (GEN.30) should become a smooth falloff.

## 1. What happens today

### 1.1 The ladder

*MAP.239: the tile root is now the grid's whole reach, 524,288 pc, three levels above the 65,536 pc root this ladder was written for. Every tile level below is three higher in the keys (a 16 pc tile is level 15, a 256 pc tile level 11), with the same edge and the same budgets; the three new coarsest levels (0 to 2) take the coarsest level's bright-star budget and list no generated stars.*

The camera orbit radius `R` sets a view radius `1.6 R`; the tile level is `floor(log2(65536 / view radius))` clamped to 0 to 12, so the level changes every time `R` halves [S: `galaxymap3d.js` `tileLevelForRadius`, `neededTiles`; `FETCH_RADIUS_FACTOR` 1.6]. Each level has fixed limits [S: `db/query.py` `generated_star_floor_sol`, `GALAXY_TILE_STAR_BUDGET`, `GALAXY_TILE_MAX_BRIGHT_STARS`, `GALAXY_TILE_POINT_BUDGET`]:

| Level | Tile edge (pc) | Camera radius served (pc) | Generated-star floor (L_sun) | Generated stars per tile | Bright stars per tile | Point objects per tile |
|---|---|---|---|---|---|---|
| 2 and coarser | 16,384 and up | 5,120 and up | none | 0 | 400 | 0 |
| 3 | 8,192 | 2,560 to 5,120 | 262 | 150 | 400 | 0 |
| 4 | 4,096 | 1,280 to 2,560 | 66 | 200 | 400 | 0 |
| 5 | 2,048 | 640 to 1,280 | 16 | 250 | 400 | 0 |
| 6 | 1,024 | 320 to 640 | 4.1 | 300 | 400 | 0 |
| 7 | 512 | 160 to 320 | 1.0 | 350 | 400 | 0 |
| 8 | 256 | 80 to 160 | 0.26 | 400 | 400 | 0 |
| 9 | 128 | 40 to 80 | 0.064 | 500 | 400 | 0 |
| 10 | 64 | 20 to 40 | 0.016 | 650 | 400 | 40 |
| 11 | 32 | 10 to 20 | 0.004 | 850 | 400 | 100 |
| 12 | 16 | 5 to 10 | none (every star) | 4,000 | 400 | 200 |

(Table 1. Radius ranges [C] from the formulas; the sweep in section 1.3 changed level at exactly these radii. The generated-star floor column and the flat 400 bright stars a tile are the values before MAP.152: MAP.152 as built, 2026-10-15, fly-through-view-distance.md: bright stars `GALAXY_TILE_BRIGHT_BUDGET` (300 for the coarsest tiles, falling to 40 at level 8 and finer, 300 in the finest) and the server floors (`generated_star_floor_sol`) are gone, which settles the old flat 400 bright stars a tile.)

### 1.2 Seven switches

1. **Level change** swaps floor (÷4), budget, per-sector allowance and the tile set at once [S].
2. **Bright-star cap.** Each tile lists its 400 most luminous pre-placed stars, an equal share from each population (`_stratified`) [S]. A tile one level finer covers an eighth of the volume and lists 400 again, so the cut falls by roughly the factor that brings 8 times as many stars into reach [C: the scratch galaxy lists 3,200 bright stars in 8 tiles at every level from 2 to 7].
3. **Detail tiles.** At a view radius of 160 pc or less (`DETAIL_VIEW_RADIUS_PC`) the finest tiles within 8 pc of the target are added, listing every star [S]. In the dense sector tested, 671 stars joined in one step at `R` = 102 pc [C].
4. **Point objects** (black holes, neutron stars, quasars) are listed only from level 10 (`POINT_PHENOMENON_MIN_LEVEL`) [S]. *As built (MAP.155): from level 8, ranked like the stars; see the note below.*
5. **No generated stars at level 2 and coarser** (the floor passes 1,000 L_sun) [S].
6. **The fade is by time.** A star new to the drawn set rises over `STAR_FADE_IN_MS` = 300 ms (`starBorn`); a star already on screen never fades again, and a removed star is simply not in the next geometry [S]. So a zoom shows nothing, then a burst, whatever the zoom speed.
7. **Cloud sprites** switch on at 2 px (`CLOUD_MIN_PX`), a hard edge [S]. *As built (MAP.155): a ramp; see the note below.*

**As built (MAP.155).** The three switches above that were still hard edges now ramp, all in `starfade.js` and the shaders in `galaxymap3d.js`:

- Point objects are listed from level 8 (160 pc tiles, `POINT_PHENOMENON_MIN_LEVEL`, with tile budgets 30, 40, 40, 100 at levels 8 to 11). The server sends them most luminous first, so `tileRanks` ranks them like stars under the key the page draws them with (`black_hole12`) and they follow the same birth radius, level glide and parent-list rules (`fadeOfStar`). Tiles coarser than level 8 list only the biggest scattered classes, which stay whole (`POINT_FADE`). A test holds the page's `POINT_RANK_MIN_LEVEL` equal to the server's constant.
- A cloud sprite's opacity ramps between 1 and 4 px of radius (`cloudSizeShare`, `CLOUD_RAMP_PX`) times the near-field fade, replacing the cut at 2 px. It can be aimed at only once its share reaches `PICK_MIN_VISIBLE` (0.35, about 2.3 px), so picking agrees with drawing.
- A star arriving with the zoom is a point of light (`POINT_PX`, `static/starpoint.js`) that grows to its size as its zoom share rises (twinned in the star vertex shader), instead of appearing at full size and only fading in. MAP.208: every star is a point of light until the camera is within `POINT_FAR_PC` of it and fully grown at `NEAR_SECTOR_PC`, one constant the Galaxy Map and the sector scene share.

Two more effects are not zoom steps but read as pops: a tile that has not arrived yet leaves its area on the previous level's stars (`carriedStars`, MAP.48) until the fetch lands, and the fetch is debounced 300 ms (`FETCH_DEBOUNCE_MS`) [S]. A cached level switches in the frame the radius crosses the boundary.

### 1.3 Measured steps

Method [C]: `research/scripts/zoomlod/steps_now.py` and `pops_now.py` ask for the tiles the page asks for at each camera radius (steps of 2^(1/8), 9%, from 6,000 pc to 4 pc), read them with the same queries `galaxy_tiles` uses, and keep the stars inside a sphere of radius `R` tan 25 deg around the target (a stand-in for the 50 degree frustum). For each pair of neighbouring steps it counts stars present in both, new at the closer step, and gone. The scratch galaxy is the quarter-scale test galaxy (44 filled sectors, 420,840 bright stars), so absolute counts are smaller than the default galaxy's; the shape of the steps does not depend on that.

Table 2: stars appearing in one 9% zoom step (every step with a change; all other steps changed nothing).

| Camera radius (pc) | Level | Dense spot: kept, new | Thin spot: kept, new |
|---|---|---|---|
| 5,502 | 2 to 3 | 3,066, +28 | 2,498, +26 |
| 2,751 | 3 to 4 | 2,283, +21 | 389, +6 |
| 1,375 | 4 to 5 | 1,048, +38 | 60, +9 |
| 688 | 5 to 6 | 203, +40 | 14, +81 |
| 344 | 6 to 7 | 50, +52 | 20, +155 |
| 172 | 7 to 8 | 21, +114 | 35, +76 |
| 102 | detail tiles | 76, +671 | none |
| 86 | 8 to 9 | 737, +29 | 21, 0 |
| 43 | 9 to 10 | 718, +1 | 5, 0 |

The dense spot is a filled sector of 547 systems near the galactic core of the scratch galaxy; the thin spot is a point with bright stars only. The largest burst is 8.8 times the stars kept (dense) and 7.8 times (thin) [C]. In the dense spot 68% of every star that appears over the whole zoom appears in the one step at 102 pc.

## 2. What "smooth" should mean

Five tests, all measurable on the tile data without a browser:

1. **No step is large.** The most stars added in one 9% zoom step, as a share of those on screen. Target: well under 100%, typical under 10%.
2. **Nothing is removed abruptly.** Stars leaving on a zoom in (not counting those leaving the frame) and arriving on a zoom out also change gradually.
3. **Reversible.** The picture depends on the camera radius only, not on the path or on timing.
4. **The network is invisible.** A tile arriving late must not change what is on screen at that moment.
5. **Bounded.** The star count a view loads stays under the MAP.109 cap (70,000, since MAP.249 45,000: `GALAXY_VIEW_MAX_STARS`, 2026-10-18) and the per-tile budgets stay as they are (they fall with distance since MAP.152, see Table 1).

## 3. Ways to do it

| Scheme | How | Smooth? | Cost and risks |
|---|---|---|---|
| A. Longer time fade | `STAR_FADE_IN_MS` from 300 to, say, 1,500 | No. The burst is still one event; it only lasts longer, and it continues to run after the zoom stops. | Trivial. Rejected. |
| B. A luminosity floor drawn as a line | Hide stars under `t(R) = 4e-5 R^2` L_sun, with a ramp over 0.6 dex; the per-level floors already are this line at the zoomed-in end of each level, so no server change | Only where a tile's list is cut by its floor. Where the 400 or per-sector caps cut it (everywhere at galaxy scale) the burst stays: worst step +390% dense, +777% thin [C] | Small, client only. A useful part of D, not enough alone. |
| C. Cross-fade whole tile levels by zoom | A star a level adds fades in with the zoom over one octave, all together | Better in a thin region, no better in a dense one: every star a level adds shares one fade curve, so a clump of equal stars all fade together (worst step +1,066% in the dense spot, table 3) | Small, client only. Needs the parent tile in memory (it is, the page prefetches one level). |
| **D. Rank birth radius** (recommended) | Each star's rank `r` in its tile's list gives it a radius `R_b = R* 2^W (N0/r)^(1/3)` at which it starts to appear; it fades in over `W` octaves of zoom. Across levels the opacity cross-fades from the parent tile's rank to the star's own | Yes: the density of drawn stars follows the density of listed stars, no tile edge in time or space shows | Client only for the first stage; section 4 |
| E. Hash thinning | Each star gets a fixed random number `u` (for example from `CRC32(id)`); it shows when `u` is under a fraction that grows smoothly with zoom | Yes, but it ignores brightness: a bright star may stay hidden while a dim one shows | Right for objects with no natural rank (rows of the same kind). Use as a tiebreak inside D. |
| F. Spacing rule | A star shows when the typical gap between stars at least as important as it, in its neighbourhood, is at least `s` pixels | Yes, and adapts to density with no budgets | Needs a density estimate per star; D is this rule with the density estimated from the list rank and the tile volume |
| G. Auto-exposure | Pick the cut each frame so about `N` stars sit in the frustum, smoothed in time (the Stellarium approach [R]) | Yes, but the count is constant, so nothing is "slowly added" | Needs all candidate stars loaded; time-dependent |
| H. Progressive streaming | Fetch the next level's stars in several pages, brightest first | It smooths the arrival, not the choice of stars | Larger change to the API; MAP.147 decides the format |
| I. Size and brightness easing | A new star starts one pixel and faint and grows to its size with the opacity | Adds to B to D: removes the "dot appears" look | Small shader change |
| J. Nested lists on the server | Make every parent list a subset of each child's list inside the child's box (same key, budgets that never fall with the level) | Needed for D to be exact; today the bright lists are stratified by population and the per-sector allowance differs by level | Server change, stage 2 |
| K. Apparent magnitude from the camera | Show a star when `L / d^2` (distance `d` from the camera, not the target) is above a threshold | Yes, and foreground stars appear before background ones | Heavier per frame; nests with D if the rank rule supplies the cap |

### 3.1 Rule D in detail

For a star of rank `r` (1 = most luminous) in the list of a tile of level `l`, edge `E`, define `R* = E / 3.2` (the closest camera radius the level serves, Table 1) and a reference size `N0` = 400. The star's birth radius is

    R_b(l, r) = R* * 2^W * (N0 / r)^(1/3)

and its opacity at camera radius `R` is `smoothstep(log2(R_b / R) / W)`, zero at `R >= R_b` and one at `R <= R_b / 2^W`. The exponent 1/3 is the one at which a tile's visible count follows the volume: halving `R` raises the visible count in a tile eightfold, and a tile one level finer holds an eighth of the volume, so the number of stars drawn per view stays level instead of jumping at the boundary. A faster growth (a smaller exponent) breaks the nesting between levels, so the number of stars in view is raised by raising the budgets of the finer levels (Table 1 already does), not by the exponent.

The level cross-fade: within the octave in which a tile level is used, `g` runs from 0 at the zoomed-out end (`R = E/1.6`) to 1 at the zoomed-in end (`R = E/3.2`), smoothed. A star's opacity is `(1 - g) a_parent + g a_own`, where `a_parent` is the opacity its rank in the parent tile's list would give (zero when the parent did not list it) and `a_own` the one its rank in this tile gives. At the boundary between two levels the second formula equals the first by construction, so nothing jumps whatever the lists contain; where the lists are nested (stage 2), `a_parent` is also what the zoomed-out view showed.

The detail tiles (8 pc around the target, listed in full) need no special case: they join whenever the view radius is 160 pc or less, their stars start with opacity 0 under the same formula, and the sector's red dwarfs come in over the zoom from about 35 pc to 8 pc, instead of all at 102 pc.

## 4. Prototype on real tile data

`research/scripts/zoomlod/pops_rank.py` (shared folder) runs scheme D, and `pops_smooth.py` scheme B, on the same tile data, over the same 84 steps. It measures the sum of opacities in the common region of two steps; "rise" is the sum of increases, "fall" the sum of decreases.

Table 3: the largest and typical rise per 9% zoom step, as a share of the opacity on screen before the step (steps with at least 20 stars on screen) [C].

| Scheme | Dense spot: median, 90th, 99th, worst | Thin spot: median, 90th, 99th, worst | Worst fall |
|---|---|---|---|
| Today (stars, instant) | 0, 0, 543, 883 | 0, 0, 217, 775 | not measured: the set changes at once |
| B. Floor line, ramp 0.6 dex | worst 390 | worst 777 | not measured |
| C. Level cross-fade (W 1) | 0, 10, 61, 1066 | 1, 29, 52, 74 | 0 |
| **D. Rank birth radius, W 1, N0 400** | 4, 17, 48, 59 | 3, 29, 60, 74 | 6 |
| D with each level's own cap as N | 3, 19, 48, 59 | 3, 29, 60, 74 | 6 |
| D with W 0.5 | 6, 20, 48, 54 | 5, 29, 46, 67 | 12 |
| D with W 2 | 3, 16, 45, 61 | 1, 29, 54, 74 | 2 |

Scheme B alone: worst step +390% dense and +777% thin. It removes the burst only where the floor cut the list; at galaxy scale the 400 cap cuts it. Scheme C, the plain level cross-fade, is as smooth as D in the thin spot and no better than today in a dense one, because the stars a level adds all get the same curve. A shorter fade (W 0.5) raises the worst fall to 12%; a longer one (W 2) lowers it to 2% and delays arrival.

Reading the table: D leaves a steady 3% to 30% per step in the middle of the zoom, which is the visible-density growth the scheme is built to produce (a fixed spot gains 30% per 9% of zoom when the count per view is level). The 59% to 74% worst steps sit just after a level change, where the cross-fade starts; a slower cross-fade (a smaller slope for `g`) would flatten them further and costs nothing but a slightly later arrival.

What D shows: at 33 pc the dense sector draws about 135 stars instead of 714, and all of them from about 8 pc [C]. That is a change of look, not an error. Whether the sector's stars should wait for sector zoom is Boss's decision (section 7).

## 5. Recommended plan

**Stage 1, client only (`galaxymap3d.js`).** In `tileStars`/`setStars`, add three per-star attributes: the birth radius from the star's own tile, the birth radius it had in the parent tile (a large sentinel when absent, or when the parent is not cached), and the tile edge. Add the camera radius as a uniform. The vertex shader computes `g` and the opacity every frame (a few arithmetic operations per star, no geometry rebuild while zooming, the cost the 300 ms fade has today). `STAR_FADE_IN_MS` stays only for stars that arrive from a fetch with opacity above zero. Keep `carriedStars`. Browser tests: a new hook returning the opacity sum at a given radius lets a test assert the step bound without pixel diffs.

**Stage 2, server.** Make the lists nested: the bright lists use one key for every level (a population weight in the key instead of the equal-share picking), and the per-level budgets never fall (they did not before MAP.152; they fall with distance since, see Table 1) so a child list always contains the parent's stars inside its box. The lists are not strictly nested today (the equal share per population and the per-sector allowance both break it); the prototype already runs on those lists, and where a child lists a star its parent did not, the star simply fades in across the child's octave, so D degrades gracefully. Stage 2 makes that exact. This stage is the one that changes the cache stamp.

*Built (MAP.154).* The bright lists rank by `luminosity x BRIGHT_STAR_POPULATION_WEIGHT[population]` (young 1, intermediate 10, old and bulge 400, bringing each population's top to about 1e6) at every level, with no height bands and no equal share; every level's budget was `GALAXY_TILE_MAX_BRIGHT_STARS`, so it never fell (superseded by MAP.152 as built, 2026-10-15, fly-through-view-distance.md: bright stars `GALAXY_TILE_BRIGHT_BUDGET` (300 for the coarsest tiles, falling to 40 at level 8 and finer, 300 in the finest) and the server floors (`generated_star_floor_sol`) are gone, which settles the old flat 400 bright stars a tile). The ranking key is the same at every level, so a child list holds every star of its parent's list inside its box when both lists use one limit (`test_a_child_tile_lists_every_star_of_its_parent_inside_it` passed at 7:31 AM EDT 2026-10-21 with the same limit of 6, 12 and 40 for every box of its chain). With the falling budgets of MAP.152 (300 to 40) a child that lists fewer stars than its parent is not guaranteed to hold all of the parent's stars inside its box; that case is to be re-checked by a test that passes each level its own budget. Decided (Boss, contradiction sweep card, 7:31 AM EDT 2026-10-21): the sentence says what the test shows and no more. The galaxy-wide sample that answers the biggest tiles uses the same key, so for those tiles nesting holds as far as the sample reaches. The generated-star lists (per-sector allowance) are not part of this change.

*Built (MAP.250).* The halo and the satellite galaxies are the lowest priority of every star budget. `LOWEST_PRIORITY_TIER` (`planetgen/db/query.py`) gives the `halo` population tier 1 and `satellite`, `satellite_young` and `satellite_intermediate` tier 2; every other population (the thin and thick disc, the bulge, the nuclear cluster and disc, the globular clusters) is tier 0. `_top_by_key` sorts by tier first, then by the MAP.154 key, so a list that cannot hold every star drops the satellite galaxies' stars first, then the halo's, and only then a star of the disc or the bulge, however bright the halo or satellite star is. All the bright-star lists go through it (`galaxy_bright_stars_in_box` on all three of its paths, `galaxy_brightest_stars`, and so the tiles, the finest tiles' 300-star detail lists, the galaxy-wide sample and a nebula's surroundings). The order is the same at every tile level, so the lists stay nested. The generated-star lists have no tier: their sectors are only generated inside the galaxy's outline, and a generated star stores no population. The view cap (70,000 stars when this was written, 45,000 since MAP.249) is the sum of the tiles' own caps, so it needs no rule of its own; a sector page shows every star its sector holds and cuts none.

**Stage 3, other objects.** Point objects: list them from level 8 (not 10) using the same rank rule, so a black hole fades in with the stars around it. Clouds: ramp opacity between 1 and 4 px (`CLOUD_MIN_PX` is a hard edge at 2 px). Star size: scheme I. Sector blocks and the fills have their own per-tile level of detail (the mega-block plan); the same "opacity as a function of camera radius" works there.

**Not recommended:** A (it hides nothing), G (constant count), and K until the perspective cost has been measured.

## 6. Things D leaves alone

- **Filters.** The class buttons and the dimmest-star slider (MAP.123) multiply on top; they are chosen by the viewer, not by zoom.
- **Charted only, wedge clip.** Unchanged.
- **Backfill tiers (GEN.30, replaced by GEN.187's four mass rings).** Around a generated sector the generator fills stars down to 100, 250, 500 and 750 L_sun in four shells at 10, 25, 50 and 100 ly [S: the retired `tuning.BRIGHT_STAR_BACKFILL_TIERS`]. On the map that is four concentric steps in the star density, independent of zoom. The generator already draws in eight bands per decade of luminosity (`canonical_bands`), so smaller tiers cost nothing extra to store; the map could also blend each shell edge by distance. Boss set these tiers; this note only records the effect.

## 7. Open questions for Boss

1. **How fast should stars come in?** Default taken: each star fades over one halving of the camera radius (`W` = 1) and a spot gains about 30% more stars per 9% of zoom. A larger `W` is slower and softer.
2. **May a dense sector show all its stars only at sector zoom?** Default taken: yes, the stars appear in rank order between 35 pc and 8 pc of view radius, bright ones first. The alternative is to bring the detail tiles in earlier, which costs more stars loaded.
3. **Point objects from level 8?** Default taken: yes (stage 3), so black holes and neutron stars appear at 160 pc instead of 40 pc.
4. **Smooth the backfill shells?** Default taken: leave GEN.30 alone.

## 8. Evidence notes

- [R] Stellarium, Celestia and Gaia Sky choose a limiting magnitude from the field of view or distance and fade objects near the limit; I did not re-read their code.
- [R] Potree-style point clouds order points in each octree node so that any prefix is a usable sample; I did not re-read it.
- [R] Dithered cross-fade between levels of detail is common in game engines.
- [C] The tile data is from the scratch galaxy of the generation-performance study (quarter scale). The default galaxy has 8 times the stars per layer count and larger caps will bind more tiles, which only widens today's bursts.
- [C] The "screen" is a sphere of radius `R` tan 25 deg around the target; the real frustum is a cone, so counts are an upper bound and the steps are the same.
- A side-by-side demo (synthetic field, same tile rules, both schemes, a stars-on-screen chart): `/mnt/project-files/research/zoomlod-demo/index.html`. In its chart the rank rule's worst step per 3.3% of zoom is +18%, against 5.7 times (dense) and 4.3 times (thin) today [C].
- Scripts: `/mnt/project-files/research/scripts/zoomlod/` (`steps_now.py`, `pops_now.py`, `pops_smooth.py`, `pops_rank.py`); results `/mnt/project-files/research/genperf-results/zoomlod_*.json`.

## 9. Stage 1 as built (MAP.153)

- `static/starfade.js`: `birthRadius`, `zoomShare`, `levelGlide`, `starZoomOpacity` (W = 1, N0 = 400, a constant for every level). The star shader repeats `zoomShare` and the glide in GLSL; the unit test (`tests/js/starfade.test.mjs`) works on the module, the browser test on the page's hook.
- `renderFromCache` gives every drawn star a vec4 `starFade` = (own birth radius, parent birth radius, tile edge, detail). The bright and generated lists are ranked separately, in the order the server sends them. A tile's parent counts only when it is cached in memory; the signature includes that, so the picture is rebuilt when a parent arrives.
- A star the parent tile listed inside a child's box that the child left out stays drawn (own radius 0) and fades out across the octave, so the lists need not be nested yet (stage 2, MAP.154, makes them so).
- Detail tiles (view radius 160 pc or less, 8 pc around the target) use their own rank alone with the level-12 edge; a star that is in both a detail tile and the regular tile takes the regular tile's fade. Stars carried over while a tile loads keep the fade they had; a coarser tile's stars carry its level's.
- Black holes, neutron stars and quasars are never ranked: their birth radius is past any camera radius (stage 3 ranks them).
- `STAR_FADE_IN_MS` still fades a star that newly joins the drawn set; since it also starts at zero opacity from zoom, the zoom fade leads.
- Test hook: `canvas.galaxyStarOpacity(radius)` returns the sum of the drawn stars' zoom opacities as if the camera orbited at `radius`, the count and the number shown, and the real camera radius.

## 10. Cut order: halo and satellite-galaxy stars go first (MAP.250, Decided)

Boss (2026-10-17, 6:48 PM EDT): "when choosing stars to show on any map, halo and satellite galaxies are the lowest priority and should always be cut first before anything else." Every place that spends a star budget (the tile lists and their per-tile budgets, the nested lists of stage 2, the on-screen limit, the view star cap (70,000, 45,000 since MAP.249), the stage views, the Sector Map) ranks a star of the halo or of a satellite galaxy below every other star whatever its luminosity; the rank in sections 1 to 9 (luminosity times the population weight) orders the stars inside a tier. The tier is part of the one key every level uses, so the nested lists stay subsets. Halo and satellite galaxies are two classes (separate switches in 8.5, `docs/design/intergalactic-and-halo-contents.md` section 10); the order between them is Decided (Boss, 6:57 PM EDT): satellite-galaxy stars are cut first, then halo stars, then everything else. The tile cache (PERF.129 to PERF.134) stores tiles built by this rule and its stamp bumps when the rule lands. The galaxy map lag Boss reported at the same time is MAP.249.

**Amendment, Decided (Boss, 7:45 PM and 7:46 PM EDT, MAP.251, MAP.252).** The full order, first to last, is: (1) satellite-galaxy stars, (2) halo stars, then (3) two ramps at the same time, phenomena from small to large by mass (nebulae are left out of the mass ranking) and stars from small dim ones to large bright ones, so the most massive and the brightest go last. Distance is part of the order: the farther an object is from the viewer, the sooner it is removed. Decided (Boss, 7:46 PM EDT, "That is exactly what I mean, thank you"): distance moves an object's place on its ramp toward the front the farther away it is,, the brightest and most massive near objects still go last, and nebulae, which have no mass rank, are cut with the phenomena ramp by their visible size (default, Boss may overrule). The server tile build knows no viewer (tiles are a pure function of key and database, cached by key and stamp), so the distance rule lives in the per-level tile budgets (`GALAXY_TILE_BRIGHT_BUDGET`, `GALAXY_TILE_STAR_BUDGET`, `GALAXY_TILE_POINT_BUDGET` falling with distance, nested), in the cut's reach (`viewport.fetch_cut_tiles`, `distance_cut_tiles`) and in the client's per-frame limit (`starfade.js`, `starmagnitude.js`), never in the tile query. Dropping is spread evenly over the viewing space (MAP.251: a cut removes the same fraction in every region, the order above decides who goes first in each region), and the dim of a view that shows charted sectors only follows the same order (MAP.252). MAP.250 is steps 1 and 2 of this order.

## 11. Even thinning inside the cut order (MAP.251)

A list that cannot hold every object is cut by one function, `_even_pick` in `src/planetgen/db/query.py`, and every tile list that spends a budget reads it (the pre-placed bright stars through `_top_by_key`, the placed point phenomena, the scattered points). The order, first kept to last: the tier (MAP.250, the halo and the satellite galaxies last), then the band of the rank key, widest first (`EVEN_DROP_BAND_DEX`, half a power of ten of the map key; for a point the size share in tenths), then a fixed number every object carries (`_thin_number`, a hash of its kind and id); the best tenth of the slots (at least one, `EVEN_DROP_ANCHOR_SHARE`) goes to the plain best whatever the band, so the single brightest object is never the one a thinned band leaves out. So the most massive and the brightest are the last to go, and the one band the cut falls in loses the same share of its objects in every region, instead of its dimmest end being truncated wherever that lies. The number belongs to the object, not the tile, so a longer list keeps every object a shorter one had and the lists stay nested wherever the database hands over the same candidates; the database reads `EVEN_DROP_POOL` times the limit of each ranking so the cut band has candidates to thin. What stays as it was: the generated stars are already shared out by sector (`SECTOR_ALLOWANCE_FUDGE`), the galaxy-wide sample (`galaxy_brightest_stars`) is the plain brightest, and the distance rule is the per-level budgets of section 10, not a thing the tile query knows. A star within half a power of ten of another can be dropped while the other is shown; that is the price of an even cut.

## 12. The dim of "Charted only" follows the cut order (MAP.252)

"Charted only" (MAP.111) dims what lies outside charted space. It dims by the cut order of section 10: every uncharted object carries a share of its light it keeps (`unchartedKeep` in `static/starpoint.js`): a satellite-galaxy star 0.03, a halo star 0.08, then two ramps at the same time from 0.12 to 0.5, a star by its luminosity (the share of `STAR_LOG_LUMINOSITY`) and a phenomenon by its mass (`massRampShare`, a solar mass to 10^8), so the smallest, dimmest object is dimmed hardest and the most massive and the brightest the least. A charted object keeps all of its light. The tier reaches the page in bits 1 and 2 of the packed tile's flag byte (`tilepack.py`, the same byte as `charted`; the halo is 1, a satellite galaxy 2, set by `db/query.py`'s `LOWEST_PRIORITY_TIER`). Distance (`unchartedFarShare`, the shader's `unchartedFar`): the farther an object is from the camera than the camera's focus, the less it keeps, from nothing extra at the focus's distance to `UNCHARTED_FAR_KEEP` (0.35) of its share at four times that far, log in between, so the far side of a view dims sooner than the near. `canvas.galaxyUnchartedKeeps()` reads the shares for tests. A scattered phenomenon is uncharted until its sector is built (`charted` in the tile's points); a nebula is not dimmed here (it has no mass rank; Boss's default was to cut it with the phenomena ramp by its visible size, which a later step can add).

## 13. The same order on the screen, and the views that spend no star budget (MAP.251)

The page's own limit follows the cut order. The limiting magnitude of section 10 (`static/starmagnitude.js`) is already a cut by apparent magnitude, so the brightest stars and the nearest are the last to go and the farther an object is, the sooner it goes, the same in every part of the view. What it lacked was the tier: a halo star or a star of a satellite galaxy was judged by its light alone. Each star now carries a penalty in magnitudes (`TIER_MAGNITUDE_PENALTY`: 3 for the halo, 6 for a satellite galaxy, a factor of about 16 and 250 in light) that is added to its apparent magnitude both where the page counts the stars on screen against the target (`histogramLimit`) and in the shader (`starTierMag`, for the opacity ramp and the brightness by flux), so when more stars are on screen than the target the limit takes the satellite-galaxy stars first, then the halo's, then the rest by apparent magnitude. The calibration limit still reads the unpenalised absolute magnitudes, so nothing near the focus is thinned more. With no more stars than the target everything shows whole, tiers included.

The stage views and the Sector Map spend no star budget and are left as they are. A stage of the drill-down draws the stored counts and sums of its blocks (`galaxy_stage`, `galaxy_block_stats`), not stars; the Sector Map draws every system of the one sector it shows. So the rule has nothing to cut there, and the lists they open from read the same ranked tile lists as the Galaxy Map.

