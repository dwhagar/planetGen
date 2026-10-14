# tests/test_nav_chart.py

"""
NAV.48: the plan for charting the uncharted sectors that block a course
(`queryDb.nav_chart_plan`) and its bypass test (`corridor.UnknownEdges`,
`nav_graph.shortest_path(blocked=...)`).
"""

import math

import pytest

from planetgen import tuning
from planetgen.db import query, store
from planetgen.db.corridor import BypassBudgetExceeded, UnknownEdges
from planetgen.galaxy.density import build_galaxy_shape
from planetgen.galaxy.geometry import neighbor_addresses, sector_position_pc
from planetgen.galaxy.nav_graph import shortest_path
from planetgen.galaxy.sector import SpaceSector
from planetgen.physics.units import ly_to_pc
from tests.test_route_edge_cases import EDGE_LY, RING, _cluster, _make_system

EDGE_PC = ly_to_pc(tuning.DEFAULT_SECTOR_EDGE_LY)
SHAPE = build_galaxy_shape(
    disk_scale_length_pc=40.0, disk_scale_height_pc=12.0, bulge_scale_radius_pc=10.0,
    bulge_amplitude=2.0, arm_count=2, pitch_angle_rad=math.radians(15), arm_amplitude=0.4,
)
LAYERS = [(-1, RING + 20), (0, RING + 20), (1, RING + 20)]
"""list: The test galaxy's outline: three layers reaching a little past the ring the sectors sit in."""


@pytest.fixture
def galaxy(mysql_config):
    store.save_galaxy_shape(SHAPE, edge_pc=EDGE_PC, outer_ring_index=RING + 20,
                            expected_system_count_at_density_1=2000.0, config=mysql_config)
    store.replace_galaxy_layers(LAYERS, config=mysql_config)
    return mysql_config


def _save_at(config, name, address, local_positions):
    """Saves one sector at `address` with systems at the given sector-local positions; returns their ids."""
    ring, layer, slot = address
    sector = SpaceSector(name, edge_ly=EDGE_LY)
    for position in local_positions:
        system, cfg = _make_system()
        sector.add_system(system, position=position, system_config=cfg)
    center = sector_position_pc(ring, layer, slot, EDGE_PC)
    sector_id = store.save_sector(sector, config=config, galaxy_position={
        "center_x_pc": center[0], "center_y_pc": center[1], "center_z_pc": center[2],
        "galactic_radius_pc": math.hypot(center[0], center[1]),
        "vertices_pc": {"inner": [], "outer": []},
        "ring_index": ring, "layer_index": layer, "ring_slot_index": slot,
    })
    conn = store.get_connection(config)
    try:
        return [row["id"] for row in conn.execute(
            "SELECT id FROM star_systems WHERE sector_id = ? ORDER BY id", (sector_id,)).fetchall()]
    finally:
        conn.close()


def _plan(config, west, east, **kwargs):
    conn = store.get_connection(config)
    try:
        return query.nav_chart_plan(conn, f"system:{west}", f"system:{east}", **kwargs)
    finally:
        conn.close()


def _ends(config, west_slot=0, east_slot=5, east_layer=0):
    west = _save_at(config, "West", (RING, 0, west_slot), _cluster((0.0, 0.0, 0.0)))
    east = _save_at(config, "East", (RING, east_layer, east_slot), _cluster((0.0, 0.0, 0.0)))
    return west[0], east[0]


def test_the_plan_lists_the_unfilled_cells_the_unknown_hop_crosses(galaxy):
    west, east = _ends(galaxy)
    plan = _plan(galaxy, west, east)

    assert plan["unknown_hops"] == 1
    assert plan["cells"], "slots 1 to 4 between the two sectors were never generated"
    assert (RING, 0, 3) in plan["cells"] and (RING, 0, 0) not in plan["cells"] and (RING, 0, 5) not in plan["cells"]
    assert len(plan["cells"]) == len(set(plan["cells"]))
    assert plan["outside_galaxy"] == 0
    assert plan["route_distance_ly"] > 0


def test_charting_the_cells_leaves_nothing_to_chart(galaxy):
    west, east = _ends(galaxy)
    for index, address in enumerate(_plan(galaxy, west, east)["cells"]):
        _save_at(galaxy, f"Chart {index}", address, _cluster((0.0, 0.0, 0.0), count=7))

    plan = _plan(galaxy, west, east)
    assert plan["unknown_hops"] == 0 and plan["cells"] == [] and plan["bypass"] is None


def test_cells_past_the_galaxys_edge_are_counted_but_never_offered(galaxy):
    west, east = _ends(galaxy, east_layer=2)    # layer 2 is outside the three-layer outline
    plan = _plan(galaxy, west, east)

    assert plan["unknown_hops"] == 1
    assert plan["outside_galaxy"] > 0
    assert all(layer <= 1 for _ring, layer, _slot in plan["cells"])


def test_the_border_adds_the_uncharted_cells_next_to_the_line(galaxy):
    west, east = _ends(galaxy)
    line = _plan(galaxy, west, east)["cells"]
    bordered = _plan(galaxy, west, east, border=True)["cells"]

    assert bordered[:len(line)] == line
    assert len(bordered) > len(line) and len(bordered) == len(set(bordered))
    extra = set(bordered) - set(line)
    assert all(any(near in line for near in neighbor_addresses(*cell)) for cell in extra)
    assert all(layer in (-1, 0, 1) for _ring, layer, _slot in bordered)
    assert (RING, 0, 0) not in bordered and (RING, 0, 5) not in bordered   # the sectors that exist


def test_the_bypass_test_says_no_route_through_known_space_exists(galaxy):
    west, east = _ends(galaxy)
    bypass = _plan(galaxy, west, east)["bypass"]

    assert bypass == {"found": False, "distance_ly": None, "checked": True}


def test_unknown_edges_answers_per_hop_and_keeps_to_its_budget(galaxy):
    _ends(galaxy)
    conn = store.get_connection(galaxy)
    try:
        def ly(address, offset=(0.0, 0.0, 0.0)):
            return tuple(c / ly_to_pc(1.0) + o for c, o in zip(sector_position_pc(*address, EDGE_PC), offset))

        edges = UnknownEdges(conn, 2)
        assert edges(ly((RING, 0, 0)), ly((RING, 0, 5))) is True         # through slots 1 to 4
        assert edges(ly((RING, 0, 5)), ly((RING, 0, 0))) is True         # the same hop, asked the other way
        assert edges(ly((RING, 0, 0)), ly((RING, 0, 0), (1.0, 0.0, 0.0))) is False   # inside one cell
        assert edges(ly((RING, 0, 0)), ly((RING, 0, 0), (1.0, 0.0, 0.0))) is False   # remembered, not counted again
        with pytest.raises(BypassBudgetExceeded):
            edges(ly((RING, 0, 0)), ly((RING, 0, 5), (0.0, 1.0, 0.0)))
    finally:
        conn.close()


def test_the_path_search_skips_blocked_edges():
    graph = {"a": {"b": 1.0, "c": 5.0}, "b": {"a": 1.0, "c": 1.0}, "c": {"a": 5.0, "b": 1.0}}
    assert shortest_path(graph, "a", "c") == (["a", "b", "c"], 2.0)
    assert shortest_path(graph, "a", "c", blocked=lambda x, y: {x, y} == {"a", "b"}) == (["a", "c"], 5.0)
    assert shortest_path(graph, "a", "c", blocked=lambda x, y: "c" in (x, y)) is None


# ---------------------------------------------------------------------------
# The command line, the Generate page and the job
# ---------------------------------------------------------------------------

def _galaxy_args(config, argv):
    from planetgen.cli import generate as generate_cli

    _parser, parsers = generate_cli.build_parser()
    parser = parsers["galaxy"]
    args = parser.parse_args([*argv, "--mysql-host", config.host, "--mysql-port", str(config.port),
                              "--mysql-user", config.user, "--mysql-password", config.password,
                              "--mysql-database", config.database])
    generate_cli.validate_shared_generation_args(args, parser)
    generate_cli.validate_galaxy_args(args, parser)
    return args


def test_the_command_line_takes_a_course_and_refuses_what_cannot_go_with_it(mysql_config):
    args = _galaxy_args(mysql_config, ["--course", "system:1", "system:2", "--course-border"])
    assert args.course == ["system:1", "system:2"] and args.course_border

    for argv in (["--course", "system:1", "system:2", "--ring", "3"],
                 ["--course", "system:1", "system:2", "--limit", "5"],
                 ["--course-border"]):
        with pytest.raises(SystemExit):
            _galaxy_args(mysql_config, argv)


def test_the_generate_page_builds_the_course_job():
    from planetgen.web import generate_page

    argv, description = generate_page.galaxy_argv({"mode": "course", "course_from": "fe81000a2b-0000005-000",
                                                   "course_to": "system:FE81000A2B-0000006-000"})
    assert argv == ["--course", "system:FE81000A2B-0000005-000", "system:FE81000A2B-0000006-000"]
    assert "uncharted sectors" in description

    argv, description = generate_page.galaxy_argv({
        "mode": "course", "course_from": "FE81000A2B-0000005-000", "course_to": "FE81000A2B-0000006-000",
        "course_border": "on", "course_confirm": "on"})
    assert argv[-2:] == ["--course-border", "--yes"] and "border" in description

    argv, _description = generate_page.galaxy_argv({"mode": "course", "course_from": "planet:fe81000a2b-0000005-001",
                                                    "course_to": "nebula:FE81000A2B-0000006-000"})
    assert argv == ["--course", "planet:FE81000A2B-0000005-001", "nebula:FE81000A2B-0000006-000"]

    for form in ({"mode": "course", "course_from": "FE81000A2B-0000005-000"},
                 {"mode": "course", "course_from": "sector:FE81000A2B", "course_to": "FE81000A2B-0000006-000"}):
        with pytest.raises(generate_page.FormError):
            generate_page.galaxy_argv(form)


def test_charting_a_course_generates_its_uncharted_sectors(galaxy):
    from planetgen.api import ids
    from planetgen.generation import run_common, run_galaxy

    west, east = _ends(galaxy)
    conn = store.get_connection(galaxy)
    try:
        refs = [ids.printed(conn, "system", system_id) for system_id in (west, east)]
    finally:
        conn.close()
    plan = run_galaxy.course_cells(galaxy, *refs)
    assert plan["cells"] and plan["unknown_hops"] == 1

    args = _galaxy_args(galaxy, ["--course", *refs, "--workers", "1"])
    args.num_systems = 2
    with run_common._generation_progress(disable=True) as progress:
        run_galaxy.run_course(args, EDGE_PC, progress)

    conn = store.get_connection(galaxy)
    try:
        made = store.get_occupied_addresses(conn, {cell[0] for cell in plan["cells"]})
    finally:
        conn.close()
    assert set(plan["cells"]) <= set(made)
    assert run_galaxy.course_cells(galaxy, *refs)["cells"] == []   # nothing left on this course to chart


def test_a_course_past_the_confirmation_size_needs_yes(galaxy, monkeypatch):
    from planetgen.api import ids
    from planetgen.generation import run_common, run_galaxy

    west, east = _ends(galaxy)
    conn = store.get_connection(galaxy)
    try:
        refs = [ids.printed(conn, "system", system_id) for system_id in (west, east)]
    finally:
        conn.close()
    monkeypatch.setattr(tuning, "NAV_CHART_CONFIRM_SECTORS", 1)
    args = _galaxy_args(galaxy, ["--course", *refs])
    with run_common._generation_progress(disable=True) as progress, pytest.raises(SystemExit):
        run_galaxy.run_course(args, EDGE_PC, progress)
    conn = store.get_connection(galaxy)
    try:
        assert conn.execute("SELECT COUNT(*) AS n FROM sectors").fetchone()["n"] == 2   # nothing was generated
    finally:
        conn.close()
