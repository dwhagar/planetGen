# tests/test_census_multiplicity_db.py

"""
TEST.162, the database half of the multiplicity audit (`tests.audit.multiplicity` is the half that draws trees by the million): the
scatter plans store their rows' members and nodes, and the whole stored census is read back and tested.

- row bookkeeping: every row's `member_count` and `system_mass_solar` are its stored members' count and mass, its nodes are one fewer than
  its members, and no member or node belongs to a row that is not there;
- no planet, moon or belt is stored: the kind of every stored member is a member kind, and the scatter tables have no body columns;
- two plans of the same galaxy store identical members;
- a built system equals its stored members (the stars, their masses and every node of the tree) and the child rows of what was built
  are gone, for every row of the sectors with the most members;
- the nightly update moves the barycentres and writes no member or node row.
"""

import collections

import pytest

from planetgen.db import orbit_trees, scatter_systems, store
from planetgen.generation import multistar, run_galaxy, run_plan, stored_trees
from planetgen.galaxy.spherical import sector_position_pc
from planetgen.physics import constants
from tests.test_bright_star_scatter import EDGE_PC, _scatter_args as star_args, _seed_galaxy as star_galaxy
from tests.test_phenomenon_scatter import _scatter_args as phenomenon_args, _seed_galaxy as phenomenon_galaxy

SOLAR = constants.SOLAR_MASS_TO_KG
BUILT_SECTORS = 1
"""int: How many of the sectors holding the most stored members the fill reuse is checked on (a sector fill builds a system around every star)."""


def _scalar(conn, sql, args=()):
    return conn.execute(sql, args).fetchone()["n"]


def _snapshot(conn):
    """Every stored member and node, as they are."""
    return (conn.execute("SELECT * FROM scatter_members ORDER BY source, source_id, member_no").fetchall(),
            conn.execute("SELECT * FROM scatter_nodes ORDER BY source, source_id, node_no").fetchall())


def _members_by_seed(conn, table):
    systems = collections.defaultdict(list)
    for row in conn.execute(
            f"SELECT t.seed AS seed, m.member_no AS member_no, m.kind AS kind, m.mass_solar AS mass, m.seed AS member_seed"
            f" FROM {table} t JOIN scatter_members m ON m.source = ? AND m.source_id = t.id"
            " ORDER BY t.seed, m.member_no", (scatter_systems.SOURCES[table],)).fetchall():
        systems[row["seed"]].append((row["kind"], row["mass"], row["member_seed"]))
    return dict(systems)


def check_bookkeeping(conn, table):
    """The stored members and nodes of every row of `table` against the row's own columns; returns how many rows have a system."""
    source = scatter_systems.SOURCES[table]
    rows = conn.execute(f"SELECT id, member_count, system_mass_solar FROM {table} WHERE member_count IS NOT NULL").fetchall()
    assert rows, f"no row of {table} has members"
    members = {row["source_id"]: (row["n"], row["mass"]) for row in conn.execute(
        "SELECT source_id, COUNT(*) AS n, SUM(mass_solar) AS mass FROM scatter_members WHERE source = ? GROUP BY source_id", (source,)).fetchall()}
    nodes = {row["source_id"]: row["n"] for row in conn.execute(
        "SELECT source_id, COUNT(*) AS n FROM scatter_nodes WHERE source = ? GROUP BY source_id", (source,)).fetchall()}
    for row in rows:
        count, mass = members[row["id"]]
        assert count == row["member_count"] >= 2, (table, row["id"])
        assert mass == pytest.approx(row["system_mass_solar"], rel=1e-9), (table, row["id"])
        assert nodes[row["id"]] == count - 1, (table, row["id"])
    ids = {row["id"] for row in rows}
    assert set(members) == ids and set(nodes) == ids, f"members or nodes of {table} without a row with members"
    kinds = {row["kind"] for row in conn.execute("SELECT DISTINCT kind FROM scatter_members WHERE source = ?", (source,)).fetchall()}
    assert kinds <= set(multistar.MEMBER_KINDS), kinds - set(multistar.MEMBER_KINDS)
    return len(rows)


def test_the_scatter_tables_hold_no_body_columns(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        for table in ("scatter_members", "scatter_nodes"):
            columns = {row["Field"] for row in conn.execute(f"SHOW COLUMNS FROM {table}").fetchall()}
            assert not {name for name in columns if any(word in name for word in ("planet", "moon", "belt", "orbit_slot"))}, (table, columns)
    finally:
        conn.close()


def test_the_phenomenon_scatters_stored_members_are_booked_and_stored_the_same_on_a_second_plan(mysql_config):
    phenomenon_galaxy(mysql_config)
    run_plan.scatter_phenomena(phenomenon_args(mysql_config))
    conn = store.get_connection(mysql_config)
    try:
        systems = check_bookkeeping(conn, "phenomenon_scatter")
        assert systems > 20
        first = _members_by_seed(conn, "phenomenon_scatter")
        kinds = {row["kind"] for row in conn.execute("SELECT DISTINCT kind FROM phenomenon_scatter WHERE member_count IS NOT NULL").fetchall()}
        assert kinds <= {"black-hole", "neutron-star", "planetary-nebula"}
        before = _snapshot(conn)
        clock = store.orbit_clock(conn, 1.0e6)
        positions = {r["id"]: (r["position_x_mpc"], r["position_y_mpc"]) for r in conn.execute(
            "SELECT id, position_x_mpc, position_y_mpc FROM phenomenon_scatter").fetchall()}
        with conn:
            store.advance_scatter_points(conn, clock)
        moved = {r["id"]: (r["position_x_mpc"], r["position_y_mpc"]) for r in conn.execute(
            "SELECT id, position_x_mpc, position_y_mpc FROM phenomenon_scatter").fetchall()}
        assert moved != positions, "the nightly update moved no barycentre"
        assert _snapshot(conn) == before, "the nightly update wrote a member or a node row"
    finally:
        conn.close()
    run_plan.scatter_phenomena(phenomenon_args(mysql_config))
    conn = store.get_connection(mysql_config)
    try:
        assert _members_by_seed(conn, "phenomenon_scatter") == first
        check_bookkeeping(conn, "phenomenon_scatter")
    finally:
        conn.close()


def test_a_built_system_equals_its_stored_members_and_the_child_rows_are_gone(mysql_config):
    star_galaxy(mysql_config)
    run_plan.scatter_bright_stars(star_args(mysql_config, "--mass", "--workers", "1"))
    conn = store.get_connection(mysql_config)
    try:
        assert check_bookkeeping(conn, "bright_stars") > 20
        first = _members_by_seed(conn, "bright_stars")
        sectors = conn.execute(
            "SELECT b.shell_index, b.band_index, b.wedge_index, COUNT(*) AS n FROM bright_stars b WHERE b.member_count > 1"
            " GROUP BY b.shell_index, b.band_index, b.wedge_index ORDER BY n DESC, b.shell_index, b.band_index, b.wedge_index LIMIT ?",
            (BUILT_SECTORS,)).fetchall()
        assert sectors
        epoch = store.get_orbit_epoch_unix(conn)
        plans = {}
        for sector in sectors:
            address = (sector["shell_index"], sector["band_index"], sector["wedge_index"])
            rows = store.bright_stars_for_sector(conn, *address)
            trees = scatter_systems.read_stored_trees(conn, "bright_stars", [r["id"] for r in rows])
            plans[address] = {row["id"]: stored_trees.plan(trees.get(row["id"]), epoch) for row in rows}
    finally:
        conn.close()
    built = 0
    for address, planned in plans.items():
        args = run_galaxy._default_generation_args(config=mysql_config)
        args.num_systems = 0
        args.planets = False     # the stars of a system are what is compared; its planets would take most of the run
        run_galaxy.generate_and_save_sector_at(args, address, sector_position_pc(*address, EDGE_PC), EDGE_PC)
        conn = store.get_connection(mysql_config)
        try:
            for row_id, plan in planned.items():
                system_id = conn.execute("SELECT star_system_id FROM bright_stars WHERE id = ?", (row_id,)).fetchone()["star_system_id"]
                assert system_id is not None, (address, row_id)
                stars = conn.execute("SELECT role, mass_kg FROM stars WHERE star_system_id = ? ORDER BY mass_kg DESC", (system_id,)).fetchall()
                nodes = conn.execute("SELECT * FROM orbit_nodes WHERE star_system_id = ? ORDER BY node_no", (system_id,)).fetchall()
                left = scatter_systems.read_stored_trees(conn, "bright_stars", [row_id])
                if plan.members:
                    lit = plan.lit_members
                    assert sorted(s["role"] for s in stars) == (
                        ["single"] if len(lit) == 1 else sorted(["primary", "secondary"] + ["member"] * (len(lit) - 2)))
                    assert [s["mass_kg"] / SOLAR for s in stars] == pytest.approx(sorted((m.mass_solar for m in lit), reverse=True), rel=1e-6)
                    assert orbit_trees.read_dark_members(conn, system_id) == plan.dark_members
                    assert len(nodes) == len(plan.nodes)
                    assert not left, "the child rows of a built row are still stored"
                    built += 1
                elif plan.complete:
                    assert [s["role"] for s in stars] == ["single"] and not nodes and not left
                else:
                    assert [s["role"] for s in stars] == ["single"] and not nodes
                    assert left, "the rows of a tree that was not built stay for the build that can"
        finally:
            conn.close()
    assert built > 0, "no row with members was built"
    conn = store.get_connection(mysql_config)
    try:
        remaining = _members_by_seed(conn, "bright_stars")
        assert set(remaining) <= set(first), "a build stored members of a row that had none"
        # every row of a sector that was built has no member or node left over unless its tree was not built
        orphans = _scalar(conn, "SELECT COUNT(*) AS n FROM scatter_members m WHERE m.source = 1 AND NOT EXISTS"
                                " (SELECT 1 FROM bright_stars b WHERE b.id = m.source_id)")
        assert orphans == 0
        assert _scalar(conn, "SELECT COUNT(*) AS n FROM scatter_nodes m WHERE m.source = 1 AND NOT EXISTS"
                             " (SELECT 1 FROM bright_stars b WHERE b.id = m.source_id)") == 0
    finally:
        conn.close()
