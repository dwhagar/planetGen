# planetgen/db/point_columns.py

"""
The stored point of an object, defined once (GEN.335; the first of the four items of GEN.294).

Every object in the galaxy is one `physics.position.SpatialPosition3D`: a position, a point mass and a movement vector
at an epoch. The tables spell those parts in their own ways (`position_x_mpc` in the scatter, `center_x_pc` for a black
hole, `position_x_km` for a planet, a mass in kilograms or in solar masses), which is how the kinds drifted apart. This
module is the one place that knows the spellings:

* `LAYOUTS`: for each table that holds every part, a `PointLayout` naming its columns, their frame and units. The SQL of
  the tables' reads and writes takes its column lists from the layout (`select`, `assignments`, `columns`), and
  `spatial_from_row` / `values_from_spatial` turn a row into a `SpatialPosition3D` and back.
* `NOT_BODIES`: the tables that carry a point but hold no object of their own (a sector, a density component).

The tests (`test_point_columns.py`) list every table that holds an object and fail when one is on neither, so
a new kind of object cannot ship without a point, a mass, a movement vector and an epoch.

The parts and how a table is seen to have them are in `detected_parts`.
"""

import dataclasses
import re

from planetgen.db import column_units
from planetgen.physics import constants
from planetgen.physics.position import SpatialPosition3D

FRAMES = ("galaxy", "sector", "system")
"""tuple: What a stored position is measured from: the galactic centre, the sector's centre, or the object's star."""

PARTS = ("point", "mass", "vector", "epoch")
"""tuple: The parts every object must hold: a position, a point mass, a movement vector and the epoch they hold at."""


@dataclasses.dataclass(frozen=True)
class PointLayout:
    """
    The columns of one table's point.

    Attributes:
        table (str): The table.
        frame (str): One of `FRAMES`: what `position` is measured from.
        position (tuple): The three position columns, x, y, z.
        velocity (tuple): The three velocity columns, km/s, in the galactic axes (relative to the star for a `system`
            frame table).
        mass (str or None): The mass column, or `None` when the mass is in another table (`mass_from`).
        mass_from (tuple or None): `(table, key, column)`: the mass is the sum of `column` over that table's rows whose
            `key` is this row's id (a star system's is its stars').
        epoch (str): The column holding when the position and velocity hold, Unix seconds.
        due (str or None): The column holding when the object next moves (`None` when its owner's does).
        galaxy_position (tuple or None): The galactic position beside a `sector` frame position, parsecs.
    """

    table: str
    frame: str
    position: tuple
    velocity: tuple
    mass: str | None
    epoch: str
    due: str | None
    mass_from: tuple | None = None
    galaxy_position: tuple | None = None

    def __post_init__(self):
        if self.frame not in FRAMES:
            raise ValueError(f"unknown frame {self.frame!r}")
        if (self.mass is None) == (self.mass_from is None):
            raise ValueError("a layout holds its mass in a column or takes it from another table, not both")

    @property
    def position_unit(self):
        """The unit the position is stored in, a length unit of the unit registry, read from the column registry."""
        return self._unit_of(self.table, self.position[0], "length")

    @property
    def mass_unit(self):
        """The unit the mass is stored in, a mass unit of the unit registry: the mass column's, or the column the
        mass is summed from."""
        if self.mass_from is not None:
            source, _key, column = self.mass_from
            return self._unit_of(source, column, "mass")
        return self._unit_of(self.table, self.mass, "mass")

    @staticmethod
    def _unit_of(table, column, kind):
        unit = column_units.unit_of(table, column)
        if column_units.unit_info(unit).kind != kind:
            raise ValueError(f"{table}.{column} is stored in {unit!r}, not a {kind} unit")
        return unit

    def columns(self):
        """Every column the layout names, in the order position, galaxy position, velocity, mass, epoch, due."""
        return (*self.position, *(self.galaxy_position or ()), *self.velocity, *((self.mass,) if self.mass else ()),
                self.epoch, *((self.due,) if self.due else ()))

    def select(self, *parts, alias=""):
        """A SELECT list of the named parts (`position`, `galaxy_position`, `velocity`, `mass`, `epoch`, `due`), each
        column prefixed by `alias.` when one is given."""
        prefix = f"{alias}." if alias else ""
        return ", ".join(f"{prefix}{name}" for name in self.named(*parts))

    def assignments(self, *parts):
        """An UPDATE's SET list for the named parts: `x = ?, y = ?, ...`."""
        return ", ".join(f"{name} = ?" for name in self.named(*parts))

    def named(self, *parts):
        """The columns of the named parts, in the order named."""
        out = []
        for part in parts:
            value = getattr(self, part)
            if value is None:
                raise ValueError(f"{self.table} has no {part} column")
            out.extend(value if isinstance(value, tuple) else (value,))
        return tuple(out)

    def values(self, row, *parts):
        """The values of the named parts in a row (a mapping by column name), in the order of `named`."""
        return tuple(row[name] for name in self.named(*parts))


_SCATTER_POINT = dict(
    position=("position_x_mpc", "position_y_mpc", "position_z_mpc"),
    velocity=("velocity_x_kms", "velocity_y_kms", "velocity_z_kms"), epoch="epoch_unix", due="next_update_due")
_BODY_POINT = dict(
    frame="system", position=("position_x_km", "position_y_km", "position_z_km"),
    velocity=("velocity_x_kms", "velocity_y_kms", "velocity_z_kms"), mass="mass_kg",
    epoch="epoch_unix", due="next_update_due")

_FREE_POINT = dict(
    frame="galaxy", position=("center_x_pc", "center_y_pc", "center_z_pc"),
    velocity=("velocity_x_kms", "velocity_y_kms", "velocity_z_kms"), epoch="epoch_unix", due="next_update_due")

LAYOUTS = {
    "planets": PointLayout(table="planets", **_BODY_POINT),
    "moons": PointLayout(table="moons", **_BODY_POINT),
    "bright_stars": PointLayout(table="bright_stars", frame="galaxy", mass="mass_kg", **_SCATTER_POINT),
    "phenomenon_scatter": PointLayout(
        table="phenomenon_scatter", frame="galaxy", mass="mass_solar", **_SCATTER_POINT),
    "star_systems": PointLayout(
        table="star_systems", frame="sector", mass=None, mass_from=("stars", "star_system_id", "mass_kg"),
        galaxy_position=("galaxy_x_pc", "galaxy_y_pc", "galaxy_z_pc"), **_SCATTER_POINT),
    "black_holes": PointLayout(table="black_holes", mass="mass_solar", **_FREE_POINT),
    "neutron_stars": PointLayout(table="neutron_stars", mass="mass_solar", **_FREE_POINT),
    "rogue_planets": PointLayout(table="rogue_planets", mass="mass_kg", **_FREE_POINT),
    "interstellar_comets": PointLayout(table="interstellar_comets", mass="mass_kg", **_FREE_POINT),
    "asteroid_fields": PointLayout(table="asteroid_fields", mass="mass_kg", **_FREE_POINT),
    "nebulae": PointLayout(table="nebulae", mass="mass_kg", **_FREE_POINT),
    "supernova_remnants": PointLayout(table="supernova_remnants", mass="mass_kg", **_FREE_POINT),
    "quasars": PointLayout(table="quasars", mass="black_hole_mass_solar", **_FREE_POINT),
    "facilities": PointLayout(table="facilities", mass="mass_kg", **_FREE_POINT),
    "vessels": PointLayout(table="vessels", mass="mass_kg", **{**_FREE_POINT, "due": None}),
    "comets": PointLayout(table="comets", **_BODY_POINT),
    "asteroid_belts": PointLayout(table="asteroid_belts", **_BODY_POINT),
    "stars": PointLayout(
        table="stars", frame="system", position=("offset_x_au", "offset_y_au", "offset_z_au"),
        velocity=("velocity_x_kms", "velocity_y_kms", "velocity_z_kms"), mass="mass_kg",
        epoch="epoch_unix", due=None),
    "scatter_members": PointLayout(
        table="scatter_members", frame="system", position=("offset_x_au", "offset_y_au", "offset_z_au"),
        velocity=("velocity_x_kms", "velocity_y_kms", "velocity_z_kms"), mass="mass_solar",
        epoch="epoch_unix", due=None),
}
"""dict: Table -> the `PointLayout` of every table that holds all four parts. A `system` frame position is measured
from the object's star (a scatter member's offset, from its scatter point)."""

NOT_BODIES = {
    "sectors": "a cell of the grid, not an object in it",
    "galaxy_components": "a density component of the galaxy (a profile about a centre), drawn from when a sector is made",
    "object_index": "the list index: a row copied from each object's own table, never a source of one",
    "galaxy_containers": "a region of space drawn from the galaxy seed (a cluster's tidal sphere, a bubble's cavity): a centre and an extent, not an object in it",
    "nebula_field": "the index of cloud sites the nebula draws from: a centre and a radius, not an object",
    "orbital_event_bodies": "a body's state as an orbital event found it, kept as a record of the event: the body's own row "
                            "stays the source",
}
"""dict: Table -> why it carries a point (or a uid) and holds no object."""

BODY_PARENT_COLUMNS = ("star_system_id", "star_id", "planet_id", "moon_id", "asteroid_belt_id", "asteroid_field_id",
                       "parent_id")
"""tuple: The columns that tie a row to the system or the body it is in. A table with one and a point stores the point
parent-relative (frame `system`, TEST.241), except the tables of `ABSOLUTE_WITH_A_PARENT`."""

ABSOLUTE_WITH_A_PARENT = {
    "bright_stars": "the scatter point a system is built from: it is the system's own position (its row points at the "
                    "system it grows into), not a body inside it",
    "black_holes": "a remnant is placed free, at its own galactic point; `star_id` only names the star it replaced. "
                   "A remnant that joins a system takes a parent-relative state with the 8.4 revision",
    "neutron_stars": "the same as black_holes",
    "rogue_planets": "a free planet: its parent link is always `free` and it keeps its own galactic point; one that "
                     "joins a system takes a parent-relative state with the 8.4 revision",
    "interstellar_comets": "the same as rogue_planets",
    "facilities": "a facility keeps a free galactic point until the 8.4 revision gives it a parent-relative "
                  "state; no 8.1 item stores a position for one",
}
"""dict: Table -> why a table with a parent link stores an absolute position. The rule (TEST.241, section 5 point 3 of
docs/design/architecture/unified-system-barycentre.md): nothing inside a system stores an absolute sector or galactic position; it
stores the state parent-relative (about the parent's central object, in the galactic axes, with its epoch), and the
barycentre is the system's row point. A new entry needs a reason and a way out."""

_POINT_STEMS = (("center_{}_pc",), ("position_{}_km",), ("position_{}_mpc",), ("offset_{}_au",), ("galaxy_{}_pc",))
_MASS_COLUMNS = ("mass_kg", "mass_solar", "black_hole_mass_solar")


def detected_parts(columns):
    """
    The parts of `PARTS` a table holds, seen from its column names: a point is a complete x, y, z group
    (`center_*_pc`, `position_*_km`, `position_*_mpc`, `offset_*_au` or `galaxy_*_pc`), a mass is `mass_kg`,
    `mass_solar` or `black_hole_mass_solar`, a vector is `velocity_*_kms`, an epoch is `epoch_unix`.
    """
    names = set(columns)
    found = set()
    if any(all(stem[0].format(axis) in names for axis in "xyz") for stem in _POINT_STEMS):
        found.add("point")
    if names & set(_MASS_COLUMNS):
        found.add("mass")
    if all(f"velocity_{axis}_kms" in names for axis in "xyz"):
        found.add("vector")
    if "epoch_unix" in names:
        found.add("epoch")
    return found


def holds_a_point(columns):
    """Whether a table carries a position at all (any complete x, y, z group)."""
    return "point" in detected_parts(columns)


def spatial_from_row(layout, row, sector_center_m=(0.0, 0.0, 0.0), star_center_m=None, star_velocity_ms=(0.0, 0.0, 0.0),
                     mass_kg=None, sector_edge_pc=None):
    """
    The row's point as a `SpatialPosition3D`, in metres (so `length_unit` "m"), m/s and kg.

    Args:
        layout (PointLayout): The table's layout.
        row (mapping): The row, by column name.
        sector_center_m (tuple): The sector's centre, galactic metres: the origin of a `sector` frame position.
        star_center_m (tuple or None): The nearest star's centre, galactic metres: the origin of a `system` frame
            position, which has none without it.
        star_velocity_ms (tuple): The nearest star's velocity, galactic axes, m/s: what a `system` frame velocity is
            measured from.
        mass_kg (float or None): The mass in kilograms for a layout that takes it from another table.
        sector_edge_pc (float or None): The sector grid's edge, so the position knows its address.

    Returns:
        SpatialPosition3D
    """
    unit = column_units.unit_info(layout.position_unit).factor
    here = tuple(row[name] * unit for name in layout.position)
    velocity = tuple(row[name] * 1000.0 for name in layout.velocity)
    if layout.frame == "galaxy":
        galactic, relative = here, False
    elif layout.frame == "sector":
        galactic, relative = tuple(c + h for c, h in zip(sector_center_m, here)), False
    else:
        if star_center_m is None:
            raise ValueError(f"a {layout.table} position is measured from its star: pass star_center_m")
        galactic, relative = tuple(c + h for c, h in zip(star_center_m, here)), True
    if layout.mass is not None:
        value = row[layout.mass]
        mass_kg = None if value is None else value * column_units.unit_info(layout.mass_unit).factor
    point = SpatialPosition3D(
        galactic, sector_center_m, star_center_galactic=star_center_m if relative else None, mass_kg=mass_kg,
        star_velocity_galactic=star_velocity_ms, epoch_unix=row[layout.epoch], sector_edge_pc=sector_edge_pc)
    # Kept exactly as stored: the galactic sum above loses the small number to the large one (a planet 1.5e8 km from a
    # star 6 kpc out is good to ten kilometres as a galactic position).
    if layout.frame == "system":
        point.set_system_cartesian(*here)
        point.set_velocity_cartesian(*velocity, frame="system")
    else:
        if layout.frame == "sector":
            point.set_sector_cartesian(*here)
        point.set_velocity_cartesian(*velocity)
    return point


def values_from_spatial(layout, point):
    """
    The columns of `layout` that a `SpatialPosition3D` holds, as `{column: value}` in the stored units: the inverse of
    `spatial_from_row`. A `mass_from` layout's mass is not included (it is the sum of other rows).
    """
    unit = column_units.unit_info(layout.position_unit).factor
    if layout.frame == "galaxy":
        here = point.get_coordinates("galactic", "cartesian")
        velocity = point.get_velocity_vector("galactic")
    elif layout.frame == "sector":
        here = point.get_coordinates("sector", "cartesian")
        velocity = point.get_velocity_vector("galactic")
    else:
        here = point.get_coordinates("system", "cartesian")
        velocity = point.get_velocity_vector("system")
        if here is None or velocity is None:
            raise ValueError(f"a {layout.table} point needs a star reference to be stored")
    out = dict(zip(layout.position, (c / unit for c in here)))
    if layout.galaxy_position:
        out.update(zip(layout.galaxy_position, (c / constants.PARSEC_M for c in point.get_coordinates("galactic", "cartesian"))))
    out.update(zip(layout.velocity, (v / 1000.0 for v in velocity)))
    if layout.mass is not None:
        mass = point.mass_kg
        out[layout.mass] = None if mass is None else mass / column_units.unit_info(layout.mass_unit).factor
    out[layout.epoch] = point.epoch_unix
    return out
