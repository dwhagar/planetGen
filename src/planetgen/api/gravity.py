# planetgen/api/gravity.py

"""
`GET /api/sectors/<id>/gravity`: the per-zone gravity grid of a sector as JSON (API.24; docs/design/orbits/gravity-map.md)
=========================================================================================================================

The grid job (`db.gravity_grid`, MAP.387) works the *local* field of a sector out and stores it; this module reads the file of the sector's current stamp
without importing numba (`db.gravity_stamp` and `galaxy.gravity_cache` are the compiled-code-free halves), adds the smooth galaxy at the zone centers when
the total field is asked for (the default, Boss 2026-11-05; `galaxy.gravity_field`), and answers one mode of it as a compact array. A sector whose grid is
not cached, or is stale, is built by a queued job (`queue.api_jobs.build_gravity_grid`): the route answers `202` with the job's `status_url` and the caller asks
again when it is done, like every long job of the API.
"""

import os

import numpy as np

from planetgen.db import gravity_stamp
from planetgen.galaxy import gravity_cache, gravity_field

SUBDIRECTORY = "gravity"
"""str: The folder of the web cache root (`web.lib.tilecache.cache_dir`) the grid files are kept in."""
FIELDS = ("total", "local")
"""tuple: What `field=` may name."""
SIGNIFICANT_DIGITS = 4
"""int: The digits kept in each number of the answer (the grid is float32; a color scale needs three)."""


def directory(root):
    """The folder the grid files live in under the cache root `root`, or `None` when the cache is off."""
    return os.path.join(root, SUBDIRECTORY) if root else None


def current_stamp(conn, sector_id, pot):
    """The stamp of the sector's grid now. Raises `ValueError` for a sector that is not stored."""
    target = gravity_stamp.sector_row(conn, sector_id)
    return gravity_stamp.source_stamp(conn, gravity_stamp.sectors_in_reach(conn, target), pot=pot)


def read_grid(conn, sector_id, pot, folder):
    """`(values, stamp)`: the sector's cached local grid `(zones^3, 10)` float32 of the current stamp, or `None` for the values when it is not there."""
    stamp = current_stamp(conn, sector_id, pot)
    return gravity_cache.read(folder, sector_id, gravity_stamp.ZONES, stamp), stamp


def _round(values):
    return [float(f"{v:.{SIGNIFICANT_DIGITS}g}") for v in values]


def payload(conn, sector_id, pot, local, stamp, mode, total=True):
    """
    The answer for `mode` of the grid `local` of the sector (see the module doc and `docs/reference/api.md`): `sector_id` is the caller's; the rest is
    `mode`, `field` ("total" or "local"), `unit`, `zones`, `zone_edge_pc`, `edge_pc`, `center_pc` and `origin_pc` (the corner of zone `[0, 0, 0]`, galaxy axes),
    `order` (the zone `[i, j, k]` is `values[(i * zones + j) * zones + k]`: `i` along x, the slowest), `values`, `minimum` and `maximum`, `saddles` (the zones
    named `[i, j, k]` with their `position_pc`) and the grid's `stamp`.
    """
    target = gravity_stamp.sector_row(conn, sector_id)
    zones = gravity_stamp.ZONES
    edge_pc = target["edge_mpc"] / 1000.0
    center = (float(target["center_x_pc"]), float(target["center_y_pc"]), float(target["center_z_pc"]))
    points = gravity_stamp.zone_points(center, edge_pc, zones)
    result = gravity_field.field(local, points, pot, total)
    values = gravity_field.scalar(result, mode)
    spots = gravity_field.saddles(result, zones)
    return {
        "mode": mode, "field": "total" if total else "local", "unit": gravity_field.UNITS[mode], "zones": zones,
        "zone_edge_pc": edge_pc / zones, "edge_pc": edge_pc, "center_pc": list(center),
        "origin_pc": [c - edge_pc / 2.0 for c in center], "order": "x-major",
        "values": _round(values), "minimum": float(f"{np.min(values):.{SIGNIFICANT_DIGITS}g}"), "maximum": float(f"{np.max(values):.{SIGNIFICANT_DIGITS}g}"),
        "saddles": [{"zone": list(zone), "position_pc": [float(v) for v in points[(zone[0] * zones + zone[1]) * zones + zone[2]]]} for zone in spots],
        "stamp": stamp,
    }
