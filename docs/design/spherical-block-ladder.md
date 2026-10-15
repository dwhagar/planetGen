# Spherical sectors: the Galaxy Map's block ladder

Status: built, 2026-10-15 (MAP.205, MAP.206), after the coordinate switch (DB.32, GEN.227, UX.111). Replaces the cylinder ladder of [drilldown-region-sizes.md](drilldown-region-sizes.md) (MAP.151, 243/81/27/9/3 with ring, wedge and slab). The grid it sits on is [spherical-sector-coordinates.md](spherical-sector-coordinates.md) (`galaxy/spherical.py`); the meridian study that fed it is [spherical-ladder-meridians.md](spherical-ladder-meridians.md).

Evidence tags: [C] computed by the code in `galaxy/drill.py` over the real grid, [B] Boss's words, [S] read in the repository.

## 1. What a block is

[B, 06:58Z] A block holds the cells whose **centre** it holds; blocks are not cut on cell edges. So the ladder picks counts, and the grid's own wedge counts (multiples of 6 with primes up to 7) are left alone.

A block is `(m, shell, band, wedge)`, `m` in **243, 81, 27, 9, 3**; a block of 1 is a sector `(1, shell, band, wedge)` as `spherical.py` names it. The ladder keeps the 3-ary step of MAP.151b. The words keep the grid's meaning: "shell" is radius, "band" is latitude, "wedge" is longitude; "ring", "layer", "slab" and "slot" mean nothing here.

- **Block shell** `s = shell // m`: block shell `s` holds sector shells `s*m` to `s*m + m - 1`; block shell 0 also holds the core. Its centreline radius is `R = s*m + (m-1)/2` edges. Radius nests exactly: a block shell sits inside block shell `s // 3` of the next size up.
- **Rows (bands of blocks):** block shell `s` is cut by colatitude into `drill_band_count(m, s)` equal steps, always an **odd** number, so the galactic plane runs through the middle of the equatorial row. Rows are numbered like sector bands: signed, 0 on the plane, positive north.
- **Wedges of a row:** the row is cut by meridians into `drill_wedge_count(m, s, band)` equal steps from +X counterclockwise.

## 2. The counts

Only the top level is free; every lower level is a whole multiple of its parent, so a parent's edge is always a child's edge.

- **243 rows:** the odd integer nearest `pi * R / m` (a row's arc is then about `m` edges). **243 wedges:** `max(1, round(2 pi R sin(theta) / m))`, `theta` the row's centre colatitude.
- **81, 27, 9, 3 rows:** the parent block shell's (`s // 3`) count times the odd integer nearest `T / parent` with `T = pi R / m`. A child row therefore lies wholly inside one parent row (`child row // (child / parent rows)`), and every parent latitude edge is a child edge.
- **81, 27, 9, 3 wedges:** the wedge count of the parent row that holds the child's row times `max(1, round(T / parent))`, `T = 2 pi R sin(theta) / m`. A child wedge lies wholly inside one parent wedge (`wedge // (child / parent wedges)`), so every parent meridian is a child meridian.

The six meridians of the grid's cap rims are not special to blocks (a block's wedge count is not a multiple of 6); blocks nest on each other, and sectors join a block by centre, so the grid's wedge counts stay what they are.

Measured over every block shell of the galaxy (shells up to 3,764, block shells to 15 at 243 and to 1,254 at 3), arc of a block relative to its target of `m` edges [C]:

| m | Blocks, whole sphere | Row arc, 1st / 50th / 99th pct | Wedge arc, 1st / 50th / 99th pct | Worst row, wedge |
|---|---|---|---|---|
| 243 | 17,165 | 0.94 / 1.00 / 1.12 | 0.96 / 1.00 / 1.07 | 0.94 to 1.12, 0.95 to 1.16 |
| 81 | 0.44 million | 0.87 / 0.99 / 1.15 | 0.86 / 1.00 / 1.20 | 0.73 to 1.27, 0.76 to 1.33 |
| 27 | 11.5 million | 0.87 / 1.00 / 1.17 | 0.85 / 1.00 / 1.20 | 0.73 to 1.99, 0.76 to 1.38 |
| 9 | 310 million | 0.86 / 1.00 / 1.17 | 0.86 / 1.00 / 1.15 | 0.71 to 1.98, 0.46 to 1.38 |
| 3 | 8.3 billion | 0.86 / 1.00 / 1.17 | 0.86 / 1.00 / 1.15 | 0.70 to 1.95, 0.46 to 1.40 |

(The counts are of the whole sphere; the outline keeps the thin disk's part. The few worst cases are blocks next to the core and the poles.) Sectors per level-3 block, over shells 900 to 902 (1.06 million blocks): 18 to 48, median 29, against a nominal 27 [C].

## 3. Sectors join a block by centre

A sector `(shell, band, wedge)` joins the level-3 block:

- shell: `shell // 3`;
- row: the row holding its band's centre colatitude (the middle of the band; the plane for band 0; the pole for a cap);
- wedge: the wedge holding its centre longitude, in integers: `((2 w + 1) * block wedges) // (2 * band wedges)`, so a centre exactly on a meridian goes to the higher wedge (the cylinder ladder's rule).

The **core** joins block shell 0, row 0 (the plane), wedge 0 at every level; a **cap** joins wedge 0 of the polar row (its centre has no longitude). A block's sectors are found in closed form (`drill_block_sectors`): per sector shell, the bands whose centre colatitude is in the row (the band centres are monotone, so a bisection), and per band the wedges whose centre longitude is in the block wedge (integer first and last). Tests walk every sector of sample shells and check the chain, the child lists and the counts agree.

## 4. Keys, stages and windows

- **Key:** `m.shell.band.wedge`, e.g. `81.12.-3.40` (a sector is `1.shell.band.wedge`). `parse_drill_key` checks the size, the shell, that the row exists and the wedge is in the row's count.
- **A stage** of block `B` (m > 1) lists its children one size down grouped by child shell, innermost first; a level-3 block's children are its sectors. The galaxy home lists the level-243 blocks the outline reaches.
- **A window** of `g` cells across (3, 5, 7, 9, 11) is **a ball**: every block of the centre's size whose centre lies within `g * m / 2` edges of the centre block's centre. It is found in closed form per block shell (the sphere meets the ball in a cap of angular radius `psi`, so the rows are those within `psi` of the centre row, and per row the wedges within a longitude half-width), never by walking the ball. Today's cylinder window is a round region too (Region in `region.py`); the shell window is the same shape with different coordinates.
- **Block bounds** for drawing are analytic: radius `(s*m - 1/2) e` to `((s+1)*m - 1/2) e` (block shell 0 starts at 0), colatitude `row * pi / rows` to `(row+1) * pi / rows`, longitude `wedge * 2 pi / wedges` to `(wedge+1) * 2 pi / wedges`. A block is a shell sector, not a box; the page draws its outline from these numbers.

## 5. What is stored

`galaxy_block_stats` keeps one row per block that holds a generated sector, at every level above a sector, keyed `(m, shell_index, band_index, wedge_index)` (renamed from `ring_index, wedge, slab` by the DB.32 migration, which also rebuilds the rows because the blocks are different blocks). Sums are unchanged: `sector_count`, `systems`, `expected_systems`, `stars`, `age_sum`, `luminosity_sol`; `ladder_step` stays (3). `db/block_stats.py` keeps them current as sectors are saved and deleted, through `drill_chain_of`.

## 6. Port to the page

`static/galaxyprisms.js` copies the count rules step for step (`Math.floor(x + 0.5)` for round-half-up, the same operand order); `tests/test_galaxydrill.py` and the `drillChain` fixture in `tests/test_js_unit.py` hold the two equal over a sweep of shells and rows. Float agreement matters only at exact `.5` ties of the count formulas, which the sweep covers.

## 7. Decisions taken (change them if Boss disagrees)

| Question | Taken | Alternative |
|---|---|---|
| Ladder sizes | 243/81/27/9/3 as MAP.151b | Fewer or more levels (counts do not depend on the choice) |
| Block rows odd | Yes, the plane runs through the equatorial row's middle | Plane on a row boundary |
| Window shape | Ball of radius `g*m/2` edges around the centre block | Cube in x, y, z |
| Core and caps | Core: row 0, wedge 0 at every level; caps: wedge 0 of the polar row | Own blocks |
| Polar rows | Their wedge count follows the same rule, so it falls to 1 at the pole | A fixed minimum of 6 |

## 8. Work, in order

1. `drill.py` and its tests (done on the branch).
2. `window.py` on the ball, with the Generate region interface where it helps.
3. `db/block_stats.py` and `query.py` stage and window reads keyed on shell, band, wedge (the columns arrive with DB.32).
4. `galaxyprisms.js`, `galaxystages.js`, then the drawing (MAP.206): `galaxystageview.js`, `galaxymap3d.js`, picking, tile cuts, the slab control becoming a band control.
5. The hand-offs of MAP.152 (galaxy, sector and system hand-offs with hysteresis) build on the new stages and are not started before 4.
