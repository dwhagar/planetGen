# planetgen/db/field_clouds.py

"""
The field clouds of one Galaxy Map node (MAP.178): the molecular clouds the galaxy's seed places (`nebula_field`, classes M and N) that no sector
has stored yet, listed beside the stored nebulae (`query.galaxy_clouds_in_box`) so the map draws a cloud before its sector is generated.

A field cloud has no row in `nebulae` until a sector stores it, and then its `nebula_field` row is `is_stored` and the stored nebula is the one
listed, so a cloud is never drawn twice. It has no object ID either (the page of a field cloud does not exist), so a cloud is `{"type": "nebula",
"id": "field-<i>.<j>.<k>.<n>", "field": True, ...}`: the name its row is keyed by. The family of M and N is `dark`.

A node lists them only while it is small: the node's reach plus the biggest cloud's radius is at most `MAX_REACH_PC`. A wider node shows the nebula
regions (`db.nodes.node_nebula_regions`), which hold these clouds, and its field would be tens of thousands of rows for a handful the map could draw.
The rows are read by the cell key (`db.region`), so only the rows near the node are read.
"""

import math

from planetgen import tuning
from planetgen.db import region
from planetgen.galaxy import nebula_field
from planetgen.physics.units import ly_to_pc

MAX_REACH_PC = 700.0
"""float: The widest a node's circumscribed sphere (half its diagonal plus the biggest cloud's radius) may be for it to list field clouds, parsecs.
At the default sector edge a block of 81 sectors reaches about 440 pc by this and fits; one of 243 reaches about 1,100 pc and does not."""

_COLUMNS = "p.cell_i, p.cell_j, p.cell_k, p.n, p.nebula_class, p.radius_ly, p.center_x_pc, p.center_y_pc, p.center_z_pc"
_POSITION = ("p.center_x_pc", "p.center_y_pc", "p.center_z_pc")


def cloud_id(row):
    """The id a field cloud is listed under, from the key of its `nebula_field` row."""
    return f"field-{int(row['cell_i'])}.{int(row['cell_j'])}.{int(row['cell_k'])}.{int(row['n'])}"


def field_clouds_in_box(conn, lo, hi, max_clouds):
    """
    The field clouds not yet stored whose sphere reaches into the box `[lo, hi)`, the largest first.

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        lo (tuple): `(x, y, z)` lower corner, parsecs.
        hi (tuple): `(x, y, z)` upper corner, parsecs.
        max_clouds (int): The most to list.

    Returns:
        list[dict]: `type` (`"nebula"`), `id` (`cloud_id`), `name` (`None`), `descriptor` (the family of the class), `class`, `radius_pc`, `x`/`y`/`z`
        and `field` (`True`), largest first, ties by the row's key; `[]` for a box wider than `MAX_REACH_PC`.
    """
    centre = tuple((a + b) / 2.0 for a, b in zip(lo, hi))
    reach = 0.5 * math.dist(lo, hi) + nebula_field.MAX_CLOUD_RADIUS_PC
    if reach > MAX_REACH_PC:
        return []
    clouds = []
    for row in region.region_rows(conn, "nebula_field", _COLUMNS, _POSITION, centre, reach, where="p.is_stored = 0"):
        radius_pc = ly_to_pc(float(row["radius_ly"]))
        position = (row["center_x_pc"], row["center_y_pc"], row["center_z_pc"])
        gap = math.sqrt(sum(max(lo[i] - position[i], 0.0, position[i] - hi[i]) ** 2 for i in range(3)))
        if radius_pc <= 0 or gap > radius_pc:
            continue
        clouds.append({
            "type": "nebula", "id": cloud_id(row), "name": None, "descriptor": tuning.NEBULA_CLASSES[row["nebula_class"]]["family"],
            "class": row["nebula_class"], "radius_pc": radius_pc,
            "x": position[0], "y": position[1], "z": position[2], "field": True,
        })
    clouds.sort(key=lambda cloud: (-cloud["radius_pc"], cloud["id"]))
    return clouds[:max_clouds]
