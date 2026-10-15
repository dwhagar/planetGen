# tests/test_nav_stones.py

"""
NAV.47: a jump through unknown space stops at the scattered bright stars,
black holes, neutron stars, quasars and hypervelocity stars of the
uncharted cells (`planetgen.db.stones`, `queryDb.nav_between`).
"""

import pytest

from planetgen import tuning
from planetgen.db import query, stones, store
from planetgen.galaxy.geometry import sector_address_at
from planetgen.physics.units import ly_to_pc
from tests.publicids import pid
from tests.test_api import client  # noqa: F401  (the fixture)
from tests.test_nav_chart import EDGE_PC, _ends, galaxy  # noqa: F401  (the fixture)

STONES_ALONG_THE_LINE = 20


def _route(config, west, east, **kwargs):
    conn = store.get_connection(config)
    try:
        return query.nav_between(conn, west, east, **kwargs)
    finally:
        conn.close()


def _line(config, west, east):
    """The galaxy-frame light-year positions of the two ends of the route's one jump through unknown space."""
    route = _route(config, west, east)["route"]
    gap = [hop for hop in route["hops"] if hop["unknown_space"]]
    assert len(gap) == 1, route["path"]
    return route["positions"][gap[0]["from"]], route["positions"][gap[0]["to"]]


def _mpc(point_ly):
    return tuple(round(ly_to_pc(c) * 1000) for c in point_ly)


def _point(a, b, fraction, off=0.0):
    return tuple(x + (y - x) * fraction + (off if axis == 2 else 0.0) for axis, (x, y) in enumerate(zip(a, b)))


def _insert_star(config, point_ly, system_id=None):
    x, y, z = _mpc(point_ly)
    ring, layer, slot = sector_address_at(tuple(ly_to_pc(c) for c in point_ly), EDGE_PC)
    conn = store.get_connection(config)
    try:
        with conn:
            conn.execute(
                "INSERT INTO bright_stars (ring_index, layer_index, ring_slot_index, position_x_mpc, position_y_mpc,"
                " position_z_mpc, population, star_type, yerkes_class, mass_kg, radius_km, temperature_k,"
                " luminosity_w, age_gy, initial_mass_sol, seed, star_system_id)"
                " VALUES (?, ?, ?, ?, ?, ?, 'young', 'B2V', 'V', 1e31, 1e6, 20000, 1e30, 0.1, 5.0, 1, ?)",
                (ring, layer, slot, x, y, z, system_id))
    finally:
        conn.close()


def _insert_scatter(config, kind, point_ly, built=False):
    x, y, z = _mpc(point_ly)
    ring, layer, slot = sector_address_at(tuple(ly_to_pc(c) for c in point_ly), EDGE_PC)
    conn = store.get_connection(config)
    try:
        with conn:
            conn.execute(
                "INSERT INTO phenomenon_scatter (ring_index, layer_index, ring_slot_index, kind, position_x_mpc,"
                " position_y_mpc, position_z_mpc, seed, built_at)"
                f" VALUES (?, ?, ?, ?, ?, ?, ?, 1, {'CURRENT_TIMESTAMP' if built else 'NULL'})",
                (ring, layer, slot, kind, x, y, z))
    finally:
        conn.close()


def _scatter_the_line(config, a, b, kinds=("black-hole", "bright-star", "neutron-star", "quasar")):
    for index in range(STONES_ALONG_THE_LINE):
        point = _point(a, b, (index + 1) / (STONES_ALONG_THE_LINE + 1))
        kind = kinds[index % len(kinds)]
        if kind == "bright-star":
            _insert_star(config, point)
        else:
            _insert_scatter(config, kind, point)


def test_a_jump_through_unknown_space_stops_at_the_scattered_objects(galaxy):  # noqa: F811
    west, east = _ends(galaxy)
    a, b = _line(galaxy, west, east)
    before = _route(galaxy, west, east)["route"]

    _scatter_the_line(galaxy, a, b)
    route = _route(galaxy, west, east)["route"]

    middle = [node for node in route["path"] if stones.is_stone(node)]
    assert middle and len(route["path"]) == len(before["path"]) + len(middle)
    assert [node for node in route["path"] if not stones.is_stone(node)] == before["path"]
    assert route["distance_ly"] == pytest.approx(before["distance_ly"], rel=0.05)
    assert sum(hop["unknown_space"] for hop in route["hops"]) == len(middle) + 1
    named = {place["node"]: place for place in route["stop_places"] if place.get("kind")}
    assert set(named) == set(middle)
    assert {place["kind"] for place in named.values()} <= {"black-hole", "bright-star", "neutron-star", "quasar"}
    assert all("uncharted sector" in place["name"] for place in named.values())
    assert route["stops"] == len(before["path"]) - 2 + len(middle)


def test_only_the_kinds_boss_listed_and_only_uncharted_ones_are_stops(galaxy):  # noqa: F811
    west, east = _ends(galaxy)
    a, b = _line(galaxy, west, east)
    _scatter_the_line(galaxy, a, b, kinds=("planetary-nebula",))     # not a stop
    for index in range(STONES_ALONG_THE_LINE):
        point = _point(a, b, (index + 1) / (STONES_ALONG_THE_LINE + 1))
        _insert_scatter(galaxy, "black-hole", point, built=True)     # built into a real object: not a stop
        _insert_star(galaxy, point, system_id=west)                  # claimed by a generated system: not a stop
    _scatter_the_line(galaxy, a, b, kinds=("supernova-remnant",))    # not a stop

    assert not any(stones.is_stone(node) for node in _route(galaxy, west, east)["route"]["path"])


def test_a_chain_longer_than_the_jump_allows_is_not_used(galaxy, monkeypatch):  # noqa: F811
    west, east = _ends(galaxy)
    a, b = _line(galaxy, west, east)
    _insert_scatter(galaxy, "neutron-star", _point(a, b, 0.5, off=20.0))     # off the line, inside the corridor

    assert any(stones.is_stone(node) for node in _route(galaxy, west, east)["route"]["path"])   # a detour is allowed
    monkeypatch.setattr(tuning, "NAV_UNKNOWN_JUMP_PENALTY", 0.0)
    assert not any(stones.is_stone(node) for node in _route(galaxy, west, east)["route"]["path"])  # none is


def test_stops_in_unknown_space_are_not_stayed_at(galaxy):  # noqa: F811
    west, east = _ends(galaxy)
    a, b = _line(galaxy, west, east)
    _scatter_the_line(galaxy, a, b)
    plain = _route(galaxy, west, east, stay_minutes=0.0)["route"]
    staying = _route(galaxy, west, east, stay_minutes=600.0)["route"]
    systems = sum(1 for node in plain["path"][1:-1] if not stones.is_stone(node))
    assert plain["stops"] > systems > 0
    extra = staying["warp_times"][0].years - plain["warp_times"][0].years
    assert extra == pytest.approx(600.0 * systems / (365.25 * 24 * 60), rel=0.01)


def test_the_chart_plan_still_names_the_cells_of_the_direct_jump(galaxy):  # noqa: F811
    west, east = _ends(galaxy)
    a, b = _line(galaxy, west, east)
    conn = store.get_connection(galaxy)
    try:
        before = query.nav_chart_plan(conn, f"system:{west}", f"system:{east}")
    finally:
        conn.close()
    _scatter_the_line(galaxy, a, b)
    conn = store.get_connection(galaxy)
    try:
        after = query.nav_chart_plan(conn, f"system:{west}", f"system:{east}")
    finally:
        conn.close()
    assert after["cells"] == before["cells"] and after["unknown_hops"] == before["unknown_hops"] == 1


def test_a_hypervelocity_star_is_where_it_is_at_the_orbit_epoch(galaxy):  # noqa: F811
    west, east = _ends(galaxy)
    a, b = _line(galaxy, west, east)
    on_the_line = _point(a, b, 0.5)
    seconds = 1.0e5 * 365.25 * 24 * 3600          # it left its stored position 100,000 years before the epoch
    speed_kms = 1000.0
    behind_ly = speed_kms * seconds / (3.0856775814913673e13) * 3.26156   # km -> parsecs -> light-years
    conn = store.get_connection(galaxy)
    try:
        epoch = float(conn.execute("SELECT UNIX_TIMESTAMP(NOW(6)) AS now").fetchone()["now"])
        stored = (on_the_line[0], on_the_line[1], on_the_line[2] - behind_ly)
        x, y, z = _mpc(stored)
        ring, layer, slot = sector_address_at(tuple(ly_to_pc(c) for c in on_the_line), EDGE_PC)
        with conn:
            conn.execute(
                "INSERT INTO phenomenon_scatter (ring_index, layer_index, ring_slot_index, kind, position_x_mpc,"
                " position_y_mpc, position_z_mpc, velocity_x_kms, velocity_y_kms, velocity_z_kms, seed, epoch_unix)"
                " VALUES (?, ?, ?, 'hypervelocity-star', ?, ?, ?, 0, 0, ?, 1, ?)",
                (ring, layer, slot, x, y, z, speed_kms, epoch - seconds))
    finally:
        conn.close()

    route = _route(galaxy, west, east)["route"]
    kinds = [place["kind"] for place in route["stop_places"] if place.get("kind")]
    assert kinds == ["hypervelocity-star"]
    stop = next(node for node in route["path"] if stones.is_stone(node))
    assert route["positions"][stop] == pytest.approx(on_the_line, abs=0.01)


def test_the_api_carries_the_stops_in_unknown_space(galaxy, client):  # noqa: F811
    west, east = _ends(galaxy)
    a, b = _line(galaxy, west, east)
    _scatter_the_line(galaxy, a, b)
    body = client.get(f"/api/nav?from={pid('system', west)}&to={pid('system', east)}").get_json()
    route = body["route"]
    stops = [place for place in route["stop_places"] if place.get("kind")]
    assert stops and all(place["node"].startswith("stone:") and place["node"] in route["positions"] for place in stops)
    assert all(hop["unknown_space"] for hop in route["hops"] if hop["from"].startswith("stone:"))
