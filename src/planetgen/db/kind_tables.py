# planetgen/db/kind_tables.py

"""
The registry of kind tables (DB.74): the tables that hold one kind of object, and what each one is, in one place.

The same list of kind tables used to be written out at 22 sites in 9 files (`rogue_planets` alone): the tables with an
object ID, the ones a sector files, the ones that take ids from the id blocks, the ones with row timestamps, the ones a
nebula can contain, and so on. Each of those lists is now a group of this registry, read through `tables(group)`, and
the web's phenomenon type and the generator's type of a phenomenon read their table from here too.

A row is a table, the reference kind of its objects (`galaxy.objectref`, `api.ids`), the generator's type for a
standalone phenomenon, and the groups it is in. The groups:

* `uid`: the objects carry a 96-bit `uid` (`api.ids.OBJECT_TABLES`).
* `id_block`: new rows take their `id` from `id_blocks` rather than AUTO_INCREMENT (`store.ID_BLOCK_TABLES`).
* `timestamped`: the table has `created_at` and `modified_at`.
* `sector_filed`: the rows are filed in a sector by `sector_id`; `system_filed`: by `star_system_id`; `star_filed`:
  a row may belong to a star by `star_id` with no `sector_id` of its own.
* `phenomenon`: a standalone phenomenon with a galaxy-frame centre, a `quadrant` and a generator type.
* `containable`: the rows can sit inside a nebula or a supernova remnant.
* `path`: the rows get a path through the sectors (`sector_paths`).

The list holds data only and imports nothing, so the store, the readers and the object references can all read it. The
order of a group is the order of `KIND_TABLES` unless the reader asks for a named order from `ORDERS`: the sites that
listed the tables by hand each had an order of their own (a progress counter counts the placed phenomena in one, the
ids list puts the facility among the system kinds), and a reader keeps the order it had. `tests/test_kind_tables.py`
checks every other place that names a kind table (the schema, the object references, the phenomenon queries) against
the registry, and every named order against its group.
"""

import dataclasses


@dataclasses.dataclass(frozen=True)
class KindTable:
    """
    One row of the registry.

    Attributes:
        table (str): The table.
        kind (str or None): The reference kind of its objects (`system`, `rogue_planet`, ...), or `None` for a table
            that holds no object of its own (`sectors` has the kind `sector`, but its id is a plain integer).
        generator_type (str or None): For a standalone phenomenon, its `tuning.PHENOMENON_TYPE_CHOICES` value.
        groups (frozenset): The groups the table is in.
    """

    table: str
    kind: str | None
    generator_type: str | None
    groups: frozenset


def _row(table, kind=None, generator_type=None, *groups):
    return KindTable(table, kind, generator_type, frozenset(groups))


_PHENOMENON = ("phenomenon", "id_block", "timestamped", "sector_filed")

KIND_TABLES = (
    _row("sectors", "sector", None, "id_block", "timestamped"),
    _row("system_configs", None, None, "id_block"),
    _row("star_systems", "system", None, "uid", "id_block", "timestamped", "sector_filed", "containable", "path"),
    _row("stars", "star", None, "uid", "id_block", "system_filed"),
    _row("planets", "planet", None, "uid", "id_block", "system_filed"),
    _row("moons", "moon", None, "uid", "id_block", "system_filed"),
    _row("asteroid_belts", "belt", None, "uid", "id_block", "system_filed"),
    _row("comets", "comet", None, "uid", "id_block", "system_filed"),
    _row("black_holes", "black_hole", "black-hole", "uid", "star_filed", "containable", *_PHENOMENON),
    _row("neutron_stars", "neutron_star", "neutron-star", "uid", "star_filed", "containable", *_PHENOMENON),
    _row("nebulae", "nebula", "nebula", "uid", "containable", *_PHENOMENON),
    _row("supernova_remnants", "supernova_remnant", "supernova-remnant", "uid", *_PHENOMENON),
    _row("quasars", "quasar", "quasar", "uid", *_PHENOMENON),
    _row("rogue_planets", "rogue_planet", "rogue-planet", "uid", "containable", "path", *_PHENOMENON),
    _row("interstellar_comets", "interstellar_comet", "comet", "uid", "containable", "path", *_PHENOMENON),
    _row("asteroid_fields", "asteroid_field", "asteroid-field", "uid", "containable", *_PHENOMENON),
    _row("facilities", "facility", None, "uid", "sector_filed", "system_filed"),
    _row("vessels", "vessel", None, "uid"),
)
"""tuple: Every kind table, in the order its groups list them."""

GROUPS = frozenset(group for row in KIND_TABLES for group in row.groups)
"""frozenset: Every group name."""

ORDERS = {
    "placed": ("black_holes", "neutron_stars", "nebulae", "supernova_remnants", "rogue_planets",
               "interstellar_comets", "asteroid_fields", "quasars"),
    "timestamped": ("sectors", "star_systems", "black_holes", "neutron_stars", "nebulae", "supernova_remnants",
                    "rogue_planets", "interstellar_comets", "asteroid_fields", "quasars"),
    "containable": ("star_systems", "rogue_planets", "interstellar_comets", "black_holes", "neutron_stars",
                    "asteroid_fields", "nebulae"),
    "slow": ("nebulae", "asteroid_fields", "black_holes", "neutron_stars", "supernova_remnants", "rogue_planets",
             "interstellar_comets", "quasars"),
    "ids": ("star_systems", "stars", "planets", "moons", "asteroid_belts", "comets", "facilities", "vessels", "nebulae",
            "asteroid_fields", "black_holes", "neutron_stars", "supernova_remnants", "rogue_planets",
            "interstellar_comets", "quasars"),
    "phenomenon_kinds": ("nebulae", "asteroid_fields", "black_holes", "neutron_stars", "supernova_remnants",
                         "rogue_planets", "interstellar_comets", "quasars"),
    "object_ids": ("rogue_planets", "black_holes", "neutron_stars", "nebulae", "supernova_remnants", "quasars",
                   "interstellar_comets", "asteroid_fields"),
}
"""dict: Name -> the tables of a group in the order one reader lists them. A reader that names one gets its group in
that order; every table of the group must be in it (`tests/test_kind_tables.py` checks the membership)."""


def _in_order(rows, order):
    if order is None:
        return list(rows)
    ranked = ORDERS[order]
    return sorted(rows, key=lambda row: ranked.index(row.table))


def tables(group, order=None):
    """The tables in `group`, in registry order or in the named `order`; `ValueError` for a group no row has."""
    if group not in GROUPS:
        raise ValueError(f"unknown kind-table group {group!r}")
    return tuple(row.table for row in _in_order((row for row in KIND_TABLES if group in row.groups), order))


def kind_tables(group, order=None):
    """`{kind: table}` for the rows of `group`."""
    return {row.kind: row.table for row in _in_order((row for row in KIND_TABLES if group in row.groups), order)}


def generator_type_tables(order=None):
    """`{generator type: table}` for the standalone phenomena."""
    return {row.generator_type: row.table for row in _in_order((row for row in KIND_TABLES if row.generator_type), order)}


def generator_types(order=None):
    """`{kind: generator type}` for the standalone phenomena."""
    return {row.kind: row.generator_type for row in _in_order((row for row in KIND_TABLES if row.generator_type), order)}
