# planetgen/db/scatter_systems.py

"""
The members of a scattered multi-star system in the database (GEN.260 tables, GEN.265 motion).

A scatter row (`bright_stars`, `phenomenon_scatter`) with companions is one point in space, the system's barycentre:
its position, velocity and epoch columns move with the nightly update like any other scatter point
(`store.advance_scatter_points`), `system_mass_solar` is the mass it carries, and `member_count` says it has members.
The members and their orbits are `scatter_members` and `scatter_nodes`, which the nightly update never writes: a
member's place is computed on read from the node phases (`physics.system_tree`). A member carries the epoch its tree
was drawn at (`epoch_unix`, GEN.337): when its offset and velocity, where it has them, hold.

`store_member_tree` writes a drawn tree (`generation.multistar.MemberTree`) under a row, `delete_member_tree` takes it
away again (the build of a row does that in its own transaction), `read_member_tree` reads it back, and
`member_positions_pc` returns where each member is, in the galaxy frame. A fill that builds the rows reads whole trees
with `read_stored_trees` (every column of the member and node rows, GEN.315) and removes the built rows' children with
`delete_member_trees`.
"""

import collections
import dataclasses

from planetgen.physics import constants
from planetgen.physics.system_tree import OrbitNode, TreeMember, member_offsets_au

SOURCES = {"bright_stars": 1, "phenomenon_scatter": 2}
"""dict: Scatter table -> the `source` code of its members and nodes."""

_MEMBER_COLUMNS = (
    "source", "source_id", "member_no", "parent_node_no", "slot", "kind", "state", "flags", "mass_solar",
    "initial_mass_solar", "seed", "epoch_unix", "offset_x_au", "offset_y_au", "offset_z_au", "velocity_x_kms",
    "velocity_y_kms", "velocity_z_kms", "decay_age_yr", "radius_au", "radius_max_au", "mass_max_solar")
_READ_ROWS = 1000
_ORBIT_NODE_FIELDS = tuple(field.name for field in dataclasses.fields(OrbitNode))
_NODE_COLUMNS = (
    "source", "source_id", "node_no", "parent_node_no", "slot", "mass_solar", "separation_km", "eccentricity",
    "inclination_deg", "ascending_node_deg", "arg_periapsis_deg", "mean_anomaly_deg", "epoch_unix", "period_years",
    "kind", "interaction", "flags")


def _source(table):
    try:
        return SOURCES[table]
    except KeyError:
        raise ValueError(f"{table!r} is not a scatter table ({', '.join(SOURCES)})") from None


def delete_member_tree(conn, table, row_id):
    """Removes a row's members and nodes (there is no foreign key, so whoever builds or rebuilds the row does this)."""
    source = _source(table)
    for child in ("scatter_members", "scatter_nodes"):
        conn.execute(f"DELETE FROM {child} WHERE source = ? AND source_id = ?", (source, row_id))


def _member_values(source, row_id, tree, epoch_unix):
    """The `scatter_members` rows of a tree; a member of a `group` node (GEN.327) also has its offset, velocity and the
    group's decay age, every other member has NULL for the seven. Every member also keeps its `radius_au` (GEN.952), and the largest
    radius and mass it has had (`radius_max_au`, `mass_max_solar`, GEN.1009)."""
    return [(source, row_id, m.member_no, m.parent_node_no, m.slot, m.kind, m.state, m.flags, m.mass_solar,
             m.initial_mass_solar, m.seed, epoch_unix, *(m.offset_au or (None,) * 3), *(m.velocity_kms or (None,) * 3),
             m.decay_age_yr, m.radius_au, m.radius_max_au, m.mass_max_solar) for m in tree.members]


def _node_values(source, row_id, tree, epoch_unix):
    return [(source, row_id, n.node_no, n.parent_node_no, n.slot, n.mass_solar, n.separation_km, n.eccentricity,
             n.inclination_deg, n.ascending_node_deg, n.arg_periapsis_deg, n.mean_anomaly_deg, epoch_unix,
             n.period_years, n.kind, n.interaction, n.flags) for n in tree.nodes]


def _insert_members(conn, members):
    conn.executemany(
        f"INSERT INTO scatter_members ({', '.join(_MEMBER_COLUMNS)}) VALUES ({', '.join('?' * len(_MEMBER_COLUMNS))})",
        members)


def _insert_nodes(conn, nodes):
    conn.executemany(
        f"INSERT INTO scatter_nodes ({', '.join(_NODE_COLUMNS)}) VALUES ({', '.join('?' * len(_NODE_COLUMNS))})", nodes)


def store_member_tree(conn, table, row_id, tree, epoch_unix):
    """
    Stores a drawn tree under a scatter row, replacing any it had: the members and nodes, and the row's
    `member_count` and `system_mass_solar`. A single object (no nodes) stores no rows and clears both columns.

    Args:
        conn (Connection): An open, read-write connection.
        table (str): `bright_stars` or `phenomenon_scatter`.
        row_id (int): The row's id.
        tree (MemberTree): What `generation.multistar.draw_member_tree` returned.
        epoch_unix (float): The time the nodes' mean anomalies hold at.
    """
    source = _source(table)
    delete_member_tree(conn, table, row_id)
    if not tree.nodes:
        conn.execute(f"UPDATE {table} SET member_count = NULL, system_mass_solar = NULL WHERE id = ?", (row_id,))
        return
    _insert_members(conn, _member_values(source, row_id, tree, epoch_unix))
    _insert_nodes(conn, _node_values(source, row_id, tree, epoch_unix))
    conn.execute(f"UPDATE {table} SET member_count = ?, system_mass_solar = ? WHERE id = ?",
                 (tree.member_count, tree.system_mass_solar, row_id))


def _double(value):
    """`value` as SQL text that MySQL reads as a DOUBLE (an exponent makes it one), the same number the driver's own
    escaping sends, so a literal row and a parameter row store the same bits."""
    text = repr(float(value))
    return text if "e" in text else text + "e0"


def _update_counts(conn, table, counts):
    """Sets `member_count` and `system_mass_solar` of many rows with one statement (PERF.158): `counts` is
    `(member_count, system_mass_solar, row_id)` each, the numbers only, so the text is built here."""
    values = ",".join(f"ROW({int(row_id)},{int(count)},{_double(mass)})" for count, mass, row_id in counts)
    conn.execute_literal(
        f"UPDATE {table} t JOIN (VALUES {values}) AS v(id, n, m) ON t.id = v.id"
        " SET t.member_count = v.n, t.system_mass_solar = v.m")


def store_member_trees(conn, table, items, epoch_unix, batch_rows=2000):
    """
    Stores the trees of many freshly written scatter rows at once (`store_member_tree` for a batch): the members and
    nodes of every tree with nodes, and their rows' `member_count` and `system_mass_solar`, in a few statements (one
    `UPDATE` for the counts of the whole batch, PERF.158).
    The rows are new, so nothing is deleted first.

    Args:
        conn (Connection): An open, read-write connection.
        table (str): `bright_stars` or `phenomenon_scatter`.
        items (iterable): `(row_id, MemberTree)`; a single object (no nodes) is left as its row is.
        epoch_unix (float): The time the nodes' mean anomalies hold at.
        batch_rows (int): How many rows' members and nodes go to a statement.

    Returns:
        int: The rows that got members.
    """
    source = _source(table)
    stored = 0
    members, nodes, counts = [], [], []

    def flush():
        if counts:
            _insert_members(conn, members)
            _insert_nodes(conn, nodes)
            _update_counts(conn, table, counts)
        members.clear()
        nodes.clear()
        counts.clear()

    for row_id, tree in items:
        if not tree.nodes:
            continue
        members.extend(_member_values(source, row_id, tree, epoch_unix))
        nodes.extend(_node_values(source, row_id, tree, epoch_unix))
        counts.append((tree.member_count, tree.system_mass_solar, row_id))
        stored += 1
        if len(counts) >= batch_rows:
            flush()
    flush()
    return stored


StoredMember = collections.namedtuple(
    "StoredMember", "member_no parent_node_no slot kind state flags mass_solar initial_mass_solar seed epoch_unix radius_au")
StoredMember.__doc__ = """One stored member with the columns a build needs: its kind and state, flag bits, present and birth
mass, the seed its star is built from, and the radius the tree drew it with (`generation.multistar.Member` less the largest
radius and mass, which only the history rule reads); a row made before v184 has no radius (None)."""

StoredTree = collections.namedtuple("StoredTree", "members nodes")
StoredTree.__doc__ = """A row's stored tree: `members` (`StoredMember`s by member number, the placed object first) and
`nodes` (`OrbitNode`s by node number, with their interaction and flags), both empty for a single object."""


def read_stored_trees(conn, table, row_ids):
    """
    The stored trees of the rows `row_ids` of `table` as `{row_id: StoredTree}`, one per row that has members (a row
    with none is absent). The members and nodes are read with two queries however many rows there are.
    """
    source = _source(table)
    ids = sorted(set(row_ids))
    members = collections.defaultdict(list)
    nodes = collections.defaultdict(list)
    for start in range(0, len(ids), _READ_ROWS):
        chunk = ids[start:start + _READ_ROWS]
        marks = ", ".join("?" * len(chunk))
        for r in conn.execute(
                "SELECT source_id, member_no, parent_node_no, slot, kind, state, flags, mass_solar, initial_mass_solar,"
                f" seed, epoch_unix, radius_au FROM scatter_members WHERE source = ? AND source_id IN ({marks})"
                " ORDER BY source_id, member_no", (source, *chunk)).fetchall():
            members[r["source_id"]].append(StoredMember(
                r["member_no"], r["parent_node_no"], r["slot"], r["kind"], r["state"], r["flags"], r["mass_solar"],
                r["initial_mass_solar"], r["seed"], r["epoch_unix"], r["radius_au"]))
        for r in conn.execute(
                "SELECT source_id, node_no, parent_node_no, slot, mass_solar, separation_km, eccentricity,"
                " inclination_deg, ascending_node_deg, arg_periapsis_deg, mean_anomaly_deg, epoch_unix, period_years,"
                f" kind, interaction, flags FROM scatter_nodes WHERE source = ? AND source_id IN ({marks})"
                " ORDER BY source_id, node_no", (source, *chunk)).fetchall():
            nodes[r["source_id"]].append(OrbitNode(**{name: r[name] for name in _ORBIT_NODE_FIELDS}))
    return {row_id: StoredTree(tuple(found), tuple(nodes[row_id])) for row_id, found in members.items()}


def delete_member_trees(conn, table, row_ids):
    """Removes the members and nodes of the rows `row_ids` of `table` (the rows a fill has just built)."""
    source = _source(table)
    ids = sorted(set(row_ids))
    for start in range(0, len(ids), _READ_ROWS):
        chunk = ids[start:start + _READ_ROWS]
        for child in ("scatter_members", "scatter_nodes"):
            conn.execute(f"DELETE FROM {child} WHERE source = ? AND source_id IN ({', '.join('?' * len(chunk))})",
                         (source, *chunk))


DELETE_BATCH = 100_000
"""int: How many member or node rows one `DELETE` of `delete_all_member_trees` removes before it commits and goes on."""


def delete_all_member_trees(conn, table):
    """
    Removes every member and node of `table`'s rows: the scatter table was emptied. A scatter leaves millions of them, and one
    `DELETE` over that many rows was minutes of a single statement (PERF.180): when the other scatter has none, the two tables
    are truncated, and when it has some, the rows go in batches of `DELETE_BATCH` (the key leads with the source, so each batch
    is a short range), each committed on its own.
    """
    source = _source(table)
    for child in ("scatter_members", "scatter_nodes"):
        if conn.execute(f"SELECT 1 FROM {child} WHERE source <> ? LIMIT 1", (source,)).fetchone() is None:
            conn.execute(f"TRUNCATE TABLE {child}")
            continue
        while conn.execute(f"DELETE FROM {child} WHERE source = ? LIMIT {DELETE_BATCH}", (source,)).rowcount:
            conn.commit()


def delete_member_trees_where(conn, table, where_sql, args=()):
    """
    Removes the members and nodes of the rows of `table` that `where_sql` selects, to be called before those rows are
    deleted (there is no foreign key to do it). `where_sql` is a condition over the table's columns written with the
    alias `t` (`"t.luminosity_w < ? AND t.star_system_id IS NULL"`).
    """
    source = _source(table)
    for child in ("scatter_members", "scatter_nodes"):
        conn.execute(f"DELETE c FROM {child} c JOIN {table} t ON c.source = ? AND c.source_id = t.id WHERE {where_sql}",
                     (source, *args))


def read_member_tree(conn, table, row_id):
    """A row's stored tree as `(members, nodes)`: lists of `TreeMember` and `OrbitNode`, both empty for a single object."""
    source = _source(table)
    members = [
        TreeMember(
            member_no=r["member_no"], parent_node_no=r["parent_node_no"], slot=r["slot"], mass_solar=r["mass_solar"],
            offset_au=(None if r["offset_x_au"] is None else (r["offset_x_au"], r["offset_y_au"], r["offset_z_au"])),
            velocity_kms=(None if r["velocity_x_kms"] is None
                          else (r["velocity_x_kms"], r["velocity_y_kms"], r["velocity_z_kms"])))
        for r in conn.execute(
            "SELECT member_no, parent_node_no, slot, mass_solar, offset_x_au, offset_y_au, offset_z_au,"
            " velocity_x_kms, velocity_y_kms, velocity_z_kms FROM scatter_members"
            " WHERE source = ? AND source_id = ? ORDER BY member_no", (source, row_id)).fetchall()]
    nodes = [
        OrbitNode(**{name: r[name] for name in (
            "node_no", "parent_node_no", "slot", "mass_solar", "separation_km", "eccentricity", "inclination_deg",
            "ascending_node_deg", "arg_periapsis_deg", "mean_anomaly_deg", "epoch_unix", "period_years", "kind")})
        for r in conn.execute(
            "SELECT node_no, parent_node_no, slot, mass_solar, separation_km, eccentricity, inclination_deg,"
            " ascending_node_deg, arg_periapsis_deg, mean_anomaly_deg, epoch_unix, period_years, kind"
            " FROM scatter_nodes WHERE source = ? AND source_id = ? ORDER BY node_no", (source, row_id)).fetchall()]
    return members, nodes


def groups_past_decay(conn):
    """
    The ids of the `bright_stars` rows whose stored young group is no younger than its `decay_age_yr` (audit 14, GEN.327):
    a group is stored only while its star is younger than that, an older one is drawn as its remnant, so a healthy scatter
    returns an empty list. (Only a star has a group: a scattered phenomenon never does.)
    """
    return [r["id"] for r in conn.execute(
        "SELECT DISTINCT t.id FROM scatter_members m JOIN bright_stars t ON m.source = ? AND m.source_id = t.id"
        " WHERE m.decay_age_yr IS NOT NULL AND t.age_gy * 1e9 >= m.decay_age_yr ORDER BY t.id",
        (SOURCES["bright_stars"],)).fetchall()]


def member_positions_pc(conn, table, row_id, unix):
    """
    Where each member of a scatter row's system is at `unix`, in the galaxy frame, in parsecs, as
    `{member_no: (x, y, z)}`: the barycentre (the row's stored point, as of its own epoch) plus each member's offset
    from the tree. A row with no members is one object at its point, member 0.

    Args:
        conn (Connection): An open connection.
        table (str): `bright_stars` or `phenomenon_scatter`.
        row_id (int): The row's id.
        unix (float): The time to place the members at.

    Raises:
        LookupError: If there is no such row.
    """
    row = conn.execute(
        f"SELECT position_x_mpc, position_y_mpc, position_z_mpc FROM {table} WHERE id = ?", (row_id,)).fetchone()
    if row is None:
        raise LookupError(f"{table} has no row {row_id}")
    centre = (row["position_x_mpc"] / 1000.0, row["position_y_mpc"] / 1000.0, row["position_z_mpc"] / 1000.0)
    members, nodes = read_member_tree(conn, table, row_id)
    if not members:
        return {0: centre}
    offsets = member_offsets_au(members, nodes, unix)
    per_pc = constants.AU_PER_PARSEC
    return {number: (centre[0] + o[0] / per_pc, centre[1] + o[1] / per_pc, centre[2] + o[2] / per_pc)
            for number, o in offsets.items()}


def point_masses_solar(conn, table, row_ids):
    """
    The mass each scatter point carries, in solar masses, as `{row_id: mass}`: the system's mass where it has members
    (`system_mass_solar`), else the object's own. This is what the Hill-sphere influence rule (GEN.109) reads.
    """
    _source(table)
    ids = list(row_ids)
    if not ids:
        return {}
    own = "mass_kg / ?" if table == "bright_stars" else "mass_solar"
    args = ((constants.SOLAR_MASS_TO_KG,) if table == "bright_stars" else ()) + tuple(ids)
    rows = conn.execute(
        f"SELECT id, COALESCE(system_mass_solar, {own}) AS mass FROM {table} WHERE id IN ({', '.join('?' * len(ids))})",
        args).fetchall()
    return {r["id"]: r["mass"] for r in rows}
