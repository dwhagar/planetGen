# Named regions (GEN.825, GEN.827, DB.107, API.35, DOC.82)

A **region** is a 9-block of the Galaxy Map's block ladder (`spherical-block-ladder.md`) that holds at least one generated
sector: about 780 sectors, 36 pc across at the default 4 pc sector edge. A region exists, and shows on the galactic map, as soon
as any sector in its cell is generated; ungenerated space has no region, and no sector line may join two unrelated places. This
page is the naming, geometry and storage rules as built; the map, the pages and the search that use them are listed at the end
as not built.

## Names

`galaxy/named_regions.py` (pure: it reads no database) names a 9-block "Word Word Region", for example "Kanbee Get Region". The
name is the gated phoneme codec (`names/gated_phoneme_codec.py`, kind `region` in `names/naming_key.KINDS`, under the galaxy's
naming key) spelling seven hex digits as two words: digits 0 to 3 are the first word, digits 4 to 6 the second. No name comes
from a list, so none comes from the Milky Way or the Solar System.

- Digits 0 to 3 are the **family hash**: the hash of the cell's 81-block (`family_of`) under the key. All the regions of one
  81-block share the first word, and the 81-block alone is "Word Reach" (`family_name`, "Kanbee Reach"). An 81-block holds about
  729 9-blocks (681 to 1,430 where the grid squeezes), not 27: the ladder steps 81, 27, 9.
- Digits 4 to 6 are the cell's **ordinal** among the 9-blocks of its family (`ordinal_in_family`, three hex digits, up to 4,095),
  so two regions of one family never share a name.
- **Neighbour rule.** A family hash that equals the unsalted hash of a neighbouring 81-block with a smaller key is salted
  (":1", ":2", ...) until it differs, so the 81-blocks around one never share a first word and of two neighbours the one with
  the larger key changes. Neighbours are the 81-blocks holding the points one block's width from the anchor along the 26
  directions of the cube (`neighbour_families`). A tie with a neighbour that is itself salted is possible in principle (about 1 in
  65,000 per pair); the audit check below finds it.
- The name depends on the cell, the key and the sector edge only, never on what has been generated, so a name never moves when
  sectors are added or deleted. A new naming key renames every region (the search table of API.36 is refilled then).

Two blocks 324 pc apart sharing a first word is harmless; neighbours doing so is what the salt prevents.

## Anchor, framing sphere, extent

- **Anchor** (`anchor_pc`): the centre of the cell's analytic bounds (`drill_block_bounds`).
- **Framing sphere** (`framing_radius_pc`): radius `min(12 pc, inscribed radius)`. 12 pc is `REGION_RADIUS_PC`, the sphere a standard
  fill of one anchor fills (`tuning.DEFAULT_GENERATE_RADIUS_PC`, about 109 sectors). The inscribed radius is the distance from the
  anchor to the nearest face. Measured at the default edge over the first 1,500 block shells: 14.6 billion 9-blocks, 2,190 with
  an inscribed radius under 12 pc (next to the polar axis) and 336 under 8 pc, none under 7.8 pc (`scripts/measure_region_cells.py`).
  So every cell has a framing sphere of at least 7.8 pc, no cell is dropped, and fewer than 1 cell in 6 million is framed by less
  than 12 pc.
- **Extent** (`region_extent`, GEN.827): the bounding sphere of the **largest face-connected clump** of the cell's generated
  sectors (neighbours by `spherical.neighbor_addresses`; ties by the smallest address): the clump's centres plus a sector's
  circumradius (sqrt(3)/2 edges), the radius within 5 percent of the smallest sphere that holds them, never under the framing
  radius (a lone sector gets the 12 pc sphere about itself). The other clumps count toward the share charted but are not inside
  the sphere. Why the largest clump and not one sphere over all: in a model of a cell of 9 by 9 by 9 sectors, two ad hoc fills of
  12 pc radius leave 2 or more clumps in 41 percent of cells, and one sphere over them would cover a gap of empty space (over 1.5
  times the largest clump's sphere) in 41 percent; three fills, 45 percent and 21 percent (`scripts/measure_region_clumps.py`).
  A standard anchor fill makes one clump per cell.
- **Share charted** is `sector_count / capacity`: the generated sectors over the sectors the cell's outline can hold. A standard
  fill is about 14 percent. Regions rank by it.

## Storage (DB.107)

Table `region_stats` (schema 161, migration 0161; `database-schema.md`): one row per 9-block holding a generated sector, key
`(shell_index, band_index, wedge_index)` (level 9 implied), with `sector_count`, `capacity`, `clump_count`, and the extent sphere
`extent_x_pc`, `extent_y_pc`, `extent_z_pc`, `extent_radius_pc` (galaxy frame, parsecs). `db/region_stats.py` keeps it, called by
`db/block_stats.py` whenever it puts right the blocks above a sector, so a sector save, a deletion and `planetgen fill-block-totals`
leave the same rows whatever the order. A cell with no generated sector has no row. Readers: `region_of(conn, position_pc)`,
`region_row`, `charted_percent`, `regions`, `labels` (names and families of the stored regions under the active key).

## On the wire (API.35)

A charted region's node carries `name` and `parent` (its family, "Word Reach"); a 9-block stage names it and lists its children
the same way (`child_m` 9); only charted regions, only with a naming key. `TILE_FORMAT` is 5. See
`../interface/galaxy-map-wire-format.md` section 0.1.

## The audit check: region names stay unique among neighbours (DOC.82)

`tests/audit/region_names.py` (`findings(cells, key)`), run by `tests/test_audit_region_names.py` and by the release audits
(`reports/audit-playbook-2026-10-22.md`, part 1.3): for each region and each of its 26 neighbouring 9-blocks, no two share a full
name, no two in different 81-blocks share a first word, and no two in one 81-block share an ordinal. A finding is a line naming both
cells. An audit run builds the cells from `region_stats.regions` (or `cells_around(point)`) and the key in force and expects an
empty list; any finding is High (a map label that reads the same for two places beside each other).

## Not built yet

Map labels by zoom (MAP.314), the infobox region and parent links (MAP.315), "in Word Word Region" on sector and object pages
(UX.231; object pages say nothing of regions until then, `../interface/object-pages.md`), the search table and "Go to" box (API.36),
the `mass_sum` column (DB.108) and the one-anchor fill that makes a region (GEN.826).
