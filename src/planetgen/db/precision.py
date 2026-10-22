# planetgen/db/precision.py

"""
The precision table (DB.77, UN10): each length column declares its role, its reference frame and how far from that
frame's origin it must hold, and the table works out the smallest body it can place.

A double's spacing grows with the number (1 AU 3.1e-5 m, 10,000 AU 0.25 m, 1 pc 4 m, 8 kpc 32,768 m), so a position
stored in a column's unit keeps a metre only within a distance that depends on the unit. This is a limit of the
floating-point type, not of the unit class (docs/design/unified-unit-quantity.md, section 6); the table makes it a
number a test can read.

* `RESOLUTION_M` and `FLOOR_BODY_M`: a body is placed well when the column's spacing at its farthest reach is at most
  a tenth of the body (a 10 m body keeps a metre).
* `spacing_m(value, unit)` and `metre_reach_m(unit)`: the spacing of a double at a value in a unit, in metres, and the
  largest distance at which a unit keeps a metre.
* `REACHES`: the largest distance each kind of position must hold, in metres, each with the source of its number.
* `PRECISION`: for each length column of `column_units`, a `Precision`: its `role`, its `frame` and the key of its
  `reach`.
* `floor_body_m(table)`: the smallest body the table's position columns place, the spacing at their reach times
  ten; `REQUIRED_FLOOR_M` holds the tables that must place the project's 10 m floor body, and
  `tests/test_precision.py` fails when one cannot.

Roles:

* `position`: a coordinate in a frame (`galaxy`: from the galactic centre, `sector`: from the sector's centre,
  `system`: from the object's star or parent).
* `separation`: a distance between two bodies of a system (an orbit's size, a limit about a star); a coordinate-free
  length whose reach is the system's.
* `size`: the extent of the object itself (a radius, a depth, a scale length); it keeps its relative precision at any
  size, so it declares no reach.
* `relative`: a length in units of the object's own size (a nebula ball's centre in nebula radii).
* `derived`: a magnitude worked out from a position (the galactic radius of a point, the distance of a listed
  neighbour); it is not stored for precision.
"""

import dataclasses
import math

from planetgen import tuning
from planetgen.db import column_units, point_columns
from planetgen.physics import constants

RESOLUTION_M = 1.0
"""float: The metre a well-placed body keeps."""

FLOOR_BODY_M = 10.0
"""float: The smallest body the project places (GEN.550's floor), metres."""

PLACEMENT_FRACTION = 10.0
"""float: A body is placed well when its size is this many times the column's spacing."""

ROLES = ("position", "separation", "size", "relative", "derived")


def spacing_m(value, unit):
    """The spacing of a double at `value` of `unit`, in metres: the smallest step the column can store there."""
    return math.ulp(abs(value)) * column_units.unit_info(unit).factor


def metre_reach_m(unit, resolution_m=RESOLUTION_M):
    """
    The largest distance in metres a column of `unit` holds to `resolution_m`.

    A double of magnitude `v` steps by `2 ** (floor(log2 v) - 52)`, so the step stays within `resolution_m / factor`
    units while `floor(log2 v)` does not pass `52 + log2(resolution_m / factor)`; the reach is the first power of two
    above that, in metres.
    """
    factor = column_units.unit_info(unit).factor
    exponent = math.floor(52 + math.log2(resolution_m / factor))
    return 2.0 ** (exponent + 1) * factor


@dataclasses.dataclass(frozen=True)
class Reach:
    """
    How far from its frame's origin a kind of position must hold.

    Attributes:
        metres (float): The distance.
        source (str): Where the number comes from.
    """

    metres: float
    source: str


def _edge_maximum_pc():
    from planetgen.galaxy import dials
    return next(lever.maximum for lever in dials.LEVERS if lever.name == "edge")


REACHES = {
    "bound-orbit": Reach(
        10_000.0 * constants.AU_TO_M,
        "the bound-orbit limit of the design note (docs/design/unified-unit-quantity.md, section 6): the multistar "
        "FRAGILE_AXIS_AU, 10,000 AU; farther than that a body is parent-relative in a wider frame"),
    "stability-ceiling": Reach(
        2.2e6 * constants.AU_TO_M,
        "the Hill sphere of a 150 solar-mass star at 8 kpc (generation.star.calculate_hill_sphere, checked by "
        "tests/test_precision.py): the farthest a comet's closed orbit reaches"),
    "sector": Reach(
        math.sqrt(3.0) / 2.0 * tuning.DEFAULT_SECTOR_EDGE_PC * constants.PARSEC_M,
        "the half-diagonal of a cell of the sector grid (tuning.DEFAULT_SECTOR_EDGE_PC)"),
    "galaxy-edge": Reach(
        _edge_maximum_pc() * constants.PARSEC_M,
        "the largest Edge the galaxy dials allow (galaxy.dials, the grid's reach)"),
}
"""dict: Reach key -> `Reach`."""


@dataclasses.dataclass(frozen=True)
class Precision:
    """
    What one length column declares.

    Attributes:
        role (str): One of `ROLES`.
        frame (str or None): `galaxy`, `sector` or `system` for a position; `system` for a separation; else `None`.
        reach (str or None): A key of `REACHES` for a position or a separation, else `None`.
    """

    role: str
    frame: str | None = None
    reach: str | None = None

    def __post_init__(self):
        if self.role not in ROLES:
            raise ValueError(f"unknown role {self.role!r}")
        if (self.role in ("position", "separation")) != (self.reach is not None):
            raise ValueError("a position or a separation declares its reach, nothing else does")
        if self.reach is not None and self.reach not in REACHES:
            raise ValueError(f"unknown reach {self.reach!r}")


_SYSTEM_REACH = {"planets": "bound-orbit", "moons": "bound-orbit", "asteroid_belts": "bound-orbit",
                 "stars": "bound-orbit", "scatter_members": "bound-orbit", "comets": "stability-ceiling"}
"""dict: Table -> the reach of its `system` frame position."""

_FRAME_REACH = {"galaxy": "galaxy-edge", "sector": "sector"}

_SEPARATION_WORDS = ("distance", "separation", "periapsis", "apoapsis", "perihelion", "orbit_distance", "reflex_offset",
                     "wobble", "min_orbit", "habitable_zone", "lower_limit", "upper_limit", "a_crit",
                     "binary_primary_position", "binary_secondary_position", "binary_mutual_position")
_DERIVED = {("nearest_systems", "distance_pc"), ("system_owners", "distance_ly")}
_RELATIVE_UNITS = ("nebula-radius",)
_OTHER_POSITIONS = {
    ("sector_path_knots", "x_pc"): "galaxy", ("sector_path_knots", "y_pc"): "galaxy",
    ("sector_path_knots", "z_pc"): "galaxy",
    ("nebula_block_stats", "centroid_x_pc"): "galaxy", ("nebula_block_stats", "centroid_y_pc"): "galaxy",
    ("nebula_block_stats", "centroid_z_pc"): "galaxy",
}
"""dict: (table, column) -> frame, the positions that are not a table's point (`point_columns.LAYOUTS`)."""


def _classify(table, column, unit):
    """The `Precision` of a length column, from the point layouts and the column's name."""
    layout = point_columns.LAYOUTS.get(table)
    if layout is not None and column in layout.position:
        reach = _SYSTEM_REACH[table] if layout.frame == "system" else _FRAME_REACH[layout.frame]
        return Precision("position", layout.frame, reach)
    if layout is not None and layout.galaxy_position and column in layout.galaxy_position:
        return Precision("position", "galaxy", _FRAME_REACH["galaxy"])
    if (table, column) in _OTHER_POSITIONS:
        return Precision("position", _OTHER_POSITIONS[table, column], _FRAME_REACH["galaxy"])
    if unit in _RELATIVE_UNITS:
        return Precision("relative")
    if (table, column) in _DERIVED or column == "galactic_radius_pc":
        return Precision("derived")
    if any(word in column for word in _SEPARATION_WORDS):
        return Precision("separation", "system", _SYSTEM_REACH.get(table, "bound-orbit"))
    return Precision("size")


PRECISION = {
    (table, column): _classify(table, column, unit) for (table, column), unit in column_units.COLUMNS.items()
    if column_units.unit_info(unit).kind == "length"
}
"""dict: (table, column) -> `Precision`, for every length column of the column registry."""


def reach_m(table, column):
    """The distance in metres a position or separation column must hold, `None` for another role."""
    precision = PRECISION[table, column]
    return None if precision.reach is None else REACHES[precision.reach].metres


def column_spacing_m(table, column):
    """The spacing in metres of the column at its reach, `None` for a column with none."""
    metres = reach_m(table, column)
    if metres is None:
        return None
    unit = column_units.unit_of(table, column)
    return spacing_m(metres / column_units.unit_info(unit).factor, unit)


def floor_body_m(table):
    """
    The smallest body the table's position and separation columns place well, in metres: the largest spacing at the
    reach of any of them, times `PLACEMENT_FRACTION`. `None` for a table with no such column.
    """
    spacings = [column_spacing_m(t, c) for (t, c), p in PRECISION.items() if t == table and p.reach is not None]
    return PLACEMENT_FRACTION * max(spacings) if spacings else None


REQUIRED_FLOOR_M = {"planets": FLOOR_BODY_M, "moons": FLOOR_BODY_M}
"""dict: Table -> the body its columns must place well. Planets and moons hold the project's 10 m floor body inside
the bound-orbit limit. A comet at the stability ceiling, a wanderer in a galaxy frame and a star in a sector are
placed at the precision the table computes, a lower one (the figures are in docs/design/unified-unit-quantity.md,
section 6.1)."""
