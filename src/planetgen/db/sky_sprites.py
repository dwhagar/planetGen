# planetgen/db/sky_sprites.py

"""
The nebulae and supernova remnants of a sphere, as arrays (VIEW.31, the data of System sky S10 and S13)
=======================================================================================================

The sky a system sees (`reports/system-sky-from-the-system-2026-11-04.md`) draws a nebula or a supernova remnant as a sprite, not a point. This is
the one read of the stored ones that the sky needs: `objects_in_sphere(conn, center, radius_pc)` returns every placed nebula (`nebulae`) and supernova
remnant (`supernova_remnants`) whose center lies inside the sphere, on `db/region.py` (the cell key, at most 256 runs, a covering index), as the plain
arrays `galaxy/sky_sprites.py` turns into sprites. The exact distance is measured here with `detmath`; a row outside the radius is dropped.

A nebula's family is its `nebula_type` (the class's family in `tuning.NEBULA_FAMILIES`); a remnant's is `"supernova-remnant"`, or `"plerion"` for a
plerion (`morphology`). Both come sorted by kind and then ID, so the same database always gives the same arrays. Reads only.
"""

import math
from typing import NamedTuple

import numpy as np

from planetgen.db import region
from planetgen.galaxy import cellkey
from planetgen.util import detmath

NEBULA, REMNANT = "nebula", "supernova_remnant"
"""str: The two kinds a row can be (`api.ids.OBJECT_TABLES` kinds, so the caller can print the IDs)."""

TABLES = {NEBULA: "nebulae", REMNANT: "supernova_remnants"}
"""dict: The table of each kind."""

_POSITION = ("p.center_x_pc", "p.center_y_pc", "p.center_z_pc")
_NEBULA_COLUMNS = ("p.id, p.nebula_type AS family, p.center_x_pc, p.center_y_pc, p.center_z_pc, p.velocity_x_kms, p.velocity_y_kms, p.velocity_z_kms, "
                   "p.epoch_unix, p.radius_ly, p.density_cm3, p.extinction_av")
_REMNANT_COLUMNS = ("p.id, p.morphology, p.center_x_pc, p.center_y_pc, p.center_z_pc, p.velocity_x_kms, p.velocity_y_kms, p.velocity_z_kms, "
                    "p.epoch_unix, p.radius_ly, p.density_cm3, p.extinction_av, p.age_years")


class SpriteObjects(NamedTuple):
    """
    The nebulae and remnants of a sphere, one row each.

    Attributes:
        kind: `(n,)` `NEBULA` or `REMNANT`.
        ident: `(n,)` the row ID.
        family: `(n,)` the family `physics.nebula_sprites.surface_brightness` takes.
        position_pc: `(n, 3)` galaxy-frame parsecs at `epoch_unix`.
        epoch_unix: `(n,)` seconds, `nan` where none is stored.
        velocity_kms: `(n, 3)` km/s, `nan` where nothing is stored.
        radius_ly, density_cm3, extinction_av: `(n,)` the stored contents.
        age_years: `(n,)` a remnant's age at its epoch, `nan` for a nebula.
    """

    kind: np.ndarray
    ident: np.ndarray
    family: np.ndarray
    position_pc: np.ndarray
    epoch_unix: np.ndarray
    velocity_kms: np.ndarray
    radius_ly: np.ndarray
    density_cm3: np.ndarray
    extinction_av: np.ndarray
    age_years: np.ndarray


def _number(value, empty=math.nan):
    return empty if value is None else float(value)


def empty():
    """A `SpriteObjects` of no rows."""
    return SpriteObjects(np.empty(0, dtype=object), np.empty(0, dtype=np.int64), np.empty(0, dtype=object), np.empty((0, 3)), np.empty(0),
                         np.empty((0, 3)), np.empty(0), np.empty(0), np.empty(0), np.empty(0))


def _rows(conn, kind, center, radius_pc, runs):
    table = TABLES[kind]
    columns = _NEBULA_COLUMNS if kind == NEBULA else _REMNANT_COLUMNS
    found = []
    for row in region.region_rows(conn, table, columns, _POSITION, center, radius_pc, "p.radius_ly > 0", (), runs):
        point = (row["center_x_pc"], row["center_y_pc"], row["center_z_pc"])
        if detmath.dist(center, point) > radius_pc:
            continue
        family = row["family"] if kind == NEBULA else ("plerion" if row["morphology"] == "plerion" else "supernova-remnant")
        found.append((kind, int(row["id"]), family, point, _number(row["epoch_unix"]),
                      (_number(row["velocity_x_kms"]), _number(row["velocity_y_kms"]), _number(row["velocity_z_kms"])),
                      float(row["radius_ly"]), _number(row["density_cm3"]), _number(row["extinction_av"], 0.0),
                      _number(row["age_years"]) if kind == REMNANT else math.nan))
    return found


def objects_in_sphere(conn, center, radius_pc):
    """
    Every placed nebula and supernova remnant whose center is within `radius_pc` of `center` (galaxy-frame parsecs), as a `SpriteObjects`.

    Raises:
        region.RegionError: A point or a radius that cannot be searched.
    """
    runs = cellkey.sphere_runs(center, radius_pc)
    found = []
    for kind in (NEBULA, REMNANT):
        found.extend(sorted(_rows(conn, kind, center, radius_pc, runs), key=lambda item: item[1]))
    if not found:
        return empty()
    return SpriteObjects(
        np.array([item[0] for item in found], dtype=object), np.array([item[1] for item in found], dtype=np.int64),
        np.array([item[2] for item in found], dtype=object), np.array([item[3] for item in found], dtype=float).reshape(-1, 3),
        np.array([item[4] for item in found], dtype=float), np.array([item[5] for item in found], dtype=float).reshape(-1, 3),
        np.array([item[6] for item in found], dtype=float), np.array([item[7] for item in found], dtype=float),
        np.array([item[8] for item in found], dtype=float), np.array([item[9] for item in found], dtype=float))
