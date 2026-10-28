# tests/audit/region_names.py

"""
The region-name audit (DOC.82): "region names stay unique among neighbours". A region is a 9-block of the block ladder
(`galaxy/named_regions.py`), named "Word Word Region" under the galaxy's naming key. The names never come from a list, so
nothing but this check keeps two regions that sit side by side from reading the same on the map. For each region in `cells`
and each of its neighbours (the 9-blocks that hold the points one block's width from its anchor along the 26 directions of
the cube, diagonals included) the check finds:

- two regions with the same full name;
- two regions of different 81-blocks whose first word is the same (the family word is salted against exactly this, so a hit
  means the salt rule broke);
- two regions of one 81-block with the same ordinal (a cell given another's name).

A finding is one line naming the two cells and the words. No findings is a pass; an audit run reads the cells of the
regions the galaxy has (`db.region_stats.regions`) or a sample of cells around chosen points (`cells_around`).
"""

import math

from planetgen import tuning
from planetgen.galaxy import drill, named_regions, spherical

LEVEL_INDEX = drill.DRILL_LEVELS.index(named_regions.REGION_LEVEL)


def cell_at(point_pc, edge_pc=tuning.DEFAULT_SECTOR_EDGE_PC):
    """The 9-block that holds the galaxy-frame point `point_pc`, or `None` outside the addressable galaxy."""
    if math.sqrt(sum(c * c for c in point_pc)) > spherical.MAX_ADDRESSABLE_SHELL * edge_pc:
        return None
    return drill.drill_chain_of(*spherical.sector_address_at(tuple(point_pc), edge_pc))[LEVEL_INDEX]


def neighbour_cells(cell, edge_pc=tuning.DEFAULT_SECTOR_EDGE_PC):
    """The 9-blocks around `cell` (see the module text), sorted, the cell's own excluded."""
    x, y, z = named_regions.anchor_pc(cell, edge_pc)
    width = named_regions.REGION_LEVEL * edge_pc
    found = set()
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for dz in (-1, 0, 1):
                if dx == dy == dz == 0:
                    continue
                other = cell_at((x + dx * width, y + dy * width, z + dz * width), edge_pc)
                if other is not None:
                    found.add(other)
    found.discard(cell)
    return sorted(found, key=tuple)


def cells_around(point_pc, radius_cells=2, edge_pc=tuning.DEFAULT_SECTOR_EDGE_PC):
    """The 9-blocks in the cube of `2 * radius_cells + 1` block widths about `point_pc`, sorted: a sample of cells to audit."""
    x, y, z = point_pc
    width = named_regions.REGION_LEVEL * edge_pc
    found = set()
    for i in range(-radius_cells, radius_cells + 1):
        for j in range(-radius_cells, radius_cells + 1):
            for k in range(-radius_cells, radius_cells + 1):
                cell = cell_at((x + i * width, y + j * width, z + k * width), edge_pc)
                if cell is not None:
                    found.add(cell)
    return sorted(found, key=tuple)


def findings(cells, key, edge_pc=tuning.DEFAULT_SECTOR_EDGE_PC):
    """The lines that say where the names of `cells` and of their neighbours clash under naming key `key` (empty: none do)."""
    lines = []
    seen = set()
    for cell in cells:
        name = named_regions.region_name(cell, key, edge_pc)
        family = named_regions.family_of(cell)
        for other in neighbour_cells(cell, edge_pc):
            pair = tuple(sorted((tuple(cell), tuple(other))))
            if pair in seen:
                continue
            seen.add(pair)
            other_name = named_regions.region_name(other, key, edge_pc)
            other_family = named_regions.family_of(other)
            if name == other_name:
                lines.append(f"same name '{name}': {tuple(cell)} and {tuple(other)}")
            elif other_family != family and name.split(" ")[0] == other_name.split(" ")[0]:
                lines.append(f"same family word '{name.split(' ')[0]}': {tuple(cell)} ({tuple(family)}) and "
                             f"{tuple(other)} ({tuple(other_family)})")
            elif other_family == family and (named_regions.ordinal_in_family(cell)
                                             == named_regions.ordinal_in_family(other)):
                lines.append(f"same ordinal in {tuple(family)}: {tuple(cell)} and {tuple(other)}")
    return lines
