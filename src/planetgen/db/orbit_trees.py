# planetgen/db/orbit_trees.py

"""
The orbit tree of a built star system in the database (GEN.311, GEN.269's storage half).

A built system's stars hang in a tree of orbits, one `orbit_nodes` row per node, with the names and units of
`scatter_nodes` (`db/scatter_systems.py`) so a fill copies a tree column for column and `physics.system_tree` places the
members the same way. A plain binary is one node: the relative orbit of the secondary about the primary, in the orbit
the pair was generated with. The old binary columns of `star_systems` stay filled until every reader has moved, and for
a pair they equal the node.

`binary_nodes` makes the node of a generated pair (with the `interaction` its Roche geometry gives), `write_nodes` and
`read_nodes` store and load a tree of `OrbitNode`s, and `member_row_values` is what each `stars` row says about its place
in the tree (`parent_node_id`, offset from the system's centre). `nodes_from_member_tree` takes the tree
`generation.multistar.draw_member_tree` drew.
"""

from collections import namedtuple

from planetgen.generation import multistar
from planetgen.physics import constants, system_tree
from planetgen.physics.system_tree import OrbitNode

AU_KM = constants.AU_TO_KM
SOLAR_MASS_KG = constants.SOLAR_MASS_TO_KG

ROOT_NODE_NO = 0
"""int: The root node of a pair, as migration 0104 numbered it."""

COLUMNS = ("node_no", "parent_node_no", "slot", "mass_solar", "separation_km", "eccentricity", "inclination_deg",
           "ascending_node_deg", "arg_periapsis_deg", "mean_anomaly_deg", "epoch_unix", "period_years", "kind",
           "interaction", "flags")
"""tuple: The columns of `orbit_nodes` after its key `star_system_id`, in the order of `OrbitNode`'s fields."""

_PairMember = namedtuple("_PairMember", ("kind", "radius_au", "mass", "initial"))

_COMPACT_CLASSES = {"VII": "white-dwarf", "NS": "neutron-star", "BH": "black-hole"}
"""dict: The Yerkes class of a remnant -> the member kind the Roche classifier knows it by."""


def _pair_member(star):
    kind = _COMPACT_CLASSES.get(star.yerkes_class, "main-sequence")
    mass_sol = star.mass / SOLAR_MASS_KG
    initial = getattr(star, "initial_mass_sol", None) or mass_sol
    return _PairMember(kind, star.radius / AU_KM, mass_sol, initial)


def pair_interaction(primary, secondary, separation_au, eccentricity):
    """The `interaction` of a generated pair of `Star`s from the Roche geometry of their radii at the pair's closest
    approach (`multistar.classify_pair`); `detached` for stars that stay inside their lobes."""
    interaction, _gainers = multistar.classify_pair(
        _pair_member(primary), _pair_member(secondary), separation_au, eccentricity)
    return interaction


def binary_nodes(star_system, epoch_unix):
    """
    The one node of a generated binary (empty for a single star), from the pair the system was generated with: its mass,
    separation, eccentricity, plane, phase and period, the phase held at `epoch_unix`. A close pair's orbit is circular,
    a wide pair keeps its own eccentricity; the argument of periapsis is zero (the phase is the mean anomaly).
    """
    kind = star_system.binary_type
    if kind is None:
        return ()
    if kind == "close":
        pair = star_system.star
        separation_au, eccentricity = pair.binary_separation_au, 0.0
        inclination, node, phase = (pair.binary_mutual_orbital_inclination_deg,
                                    pair.binary_mutual_orbital_ascending_node_deg, pair.binary_mutual_orbital_phase_deg)
        period = pair.binary_mutual_orbital_period_years
    else:
        pair = star_system.wide_binary
        separation_au, eccentricity = pair.separation_au, pair.eccentricity
        inclination, node, phase, period = pair.inclination_deg, pair.ascending_node_deg, pair.phase_deg, pair.period_years
    mass_solar = sum(star.mass for star in star_system.stars) / SOLAR_MASS_KG
    return (OrbitNode(
        node_no=ROOT_NODE_NO, parent_node_no=None, slot=0, mass_solar=mass_solar, separation_km=separation_au * AU_KM,
        eccentricity=eccentricity, inclination_deg=inclination, ascending_node_deg=node, arg_periapsis_deg=0.0,
        mean_anomaly_deg=phase, epoch_unix=float(epoch_unix), period_years=period, kind=kind,
        interaction=pair_interaction(star_system.primary_star, star_system.secondary_star, separation_au, eccentricity),
        flags=0),)


def nodes_from_member_tree(tree, epoch_unix):
    """The `OrbitNode`s of a `multistar.MemberTree`, their phases held at `epoch_unix`."""
    return tuple(OrbitNode(
        node_no=n.node_no, parent_node_no=n.parent_node_no, slot=n.slot, mass_solar=n.mass_solar,
        separation_km=n.separation_km, eccentricity=n.eccentricity, inclination_deg=n.inclination_deg,
        ascending_node_deg=n.ascending_node_deg, arg_periapsis_deg=n.arg_periapsis_deg,
        mean_anomaly_deg=n.mean_anomaly_deg, epoch_unix=float(epoch_unix), period_years=n.period_years, kind=n.kind,
        interaction=n.interaction, flags=n.flags) for n in tree.nodes)


def write_nodes(conn, star_system_id, nodes, replace=True):
    """Stores `nodes` as the system's tree, replacing any it had (`replace=False` for a system just inserted, which has
    none: the rows then go through `Connection.batched` with the rest of the sector)."""
    if replace:
        conn.execute("DELETE FROM orbit_nodes WHERE star_system_id = ?", (star_system_id,))
    sql = f"INSERT INTO orbit_nodes (star_system_id, {', '.join(COLUMNS)}) VALUES (?, {', '.join('?' * len(COLUMNS))})"
    for node in nodes:
        conn.execute(sql, (star_system_id, *(getattr(node, column) for column in COLUMNS)))


def read_nodes(conn, star_system_id):
    """The system's stored tree as a tuple of `OrbitNode`s, by node number (empty for a single star)."""
    rows = conn.execute(
        f"SELECT {', '.join(COLUMNS)} FROM orbit_nodes WHERE star_system_id = ? ORDER BY node_no",
        (star_system_id,)).fetchall()
    return tuple(OrbitNode(**{column: row[column] for column in COLUMNS}) for row in rows)


def star_places(star_system):
    """
    Where each star of a system sits in its tree, in `star_system.stars` order: `(parent_node_id, offset_au, slot)`
    with `parent_node_id` and `slot` `None` and a zero offset for a single star. A pair the fill generated hangs on the
    root node, the primary in slot 0 and the secondary in slot 1, at the stored binary position (the primary opposite the
    secondary in proportion to the masses). A system with `member_stars` (GEN.630) says where each star hangs in
    `star_slots`, and every offset is the star's place in the whole tree at the nodes' epoch.
    """
    slots = star_system.star_slots
    if slots is not None:
        offsets = system_tree.star_offsets_au([star.mass / SOLAR_MASS_KG for star in star_system.stars], slots,
                                              star_system.orbit_nodes, star_system.orbit_nodes[0].epoch_unix)
        return [(parent, offset, slot) for offset, (parent, slot) in zip(offsets, slots)]
    if star_system.binary_type is None:
        return [(None, (0.0, 0.0, 0.0), None)]
    if star_system.binary_type == "close":
        pair = star_system.star
        return [(ROOT_NODE_NO, (pair.binary_primary_position_x, pair.binary_primary_position_y,
                                pair.binary_primary_position_z), 0),
                (ROOT_NODE_NO, (pair.binary_secondary_position_x, pair.binary_secondary_position_y,
                                pair.binary_secondary_position_z), 1)]
    pair = star_system.wide_binary
    return [(ROOT_NODE_NO, (pair.primary_position_x_au, pair.primary_position_y_au, pair.primary_position_z_au), 0),
            (ROOT_NODE_NO, (pair.secondary_position_x_au, pair.secondary_position_y_au, pair.secondary_position_z_au), 1)]
