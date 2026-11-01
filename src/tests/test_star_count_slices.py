# tests/test_star_count_slices.py

"""UX.239, the later slices: the stars filter and the star count in `db.near` (Nearby), `GET /api/near?stars=`, `planetgen.cli.query systems|near --stars` and the
Nearby page's menu. The Systems list, its `stars` filter and sort are in `test_star_count_systems.py`."""

import sys

import pytest

from planetgen.cli import query as query_cli
from planetgen.db import near, query
from planetgen.galaxy.spherical import sector_position_pc
from planetgen.util import star_count
from tests.bughunt_support import mysql_argv
from tests.test_near import ADDRESSES, EDGE_PC, _conn, _everything, client, galaxy  # noqa: F401


def _counts(conn):
    """`{system id: stars it holds}` read straight from the `stars` table."""
    return {row["id"]: row["n"] for row in conn.execute(
        "SELECT ss.id, (SELECT COUNT(*) FROM stars st WHERE st.star_system_id = ss.id) AS n FROM star_systems ss").fetchall()}


# --- Nearby -------------------------------------------------------------------------------------------------------------------

def test_a_stars_filter_keeps_the_systems_that_hold_that_many_and_what_is_in_them(galaxy):
    conn = _conn(galaxy)
    try:
        counts = _counts(conn)
        place = near.place_from_point(sector_position_pc(*ADDRESSES[0], EDGE_PC))
        _page, everything = _everything(conn, place, 20.0)
        all_systems = {row["id"] for row in everything if row["kind"] == "system"}
        assert all_systems
        for stars in sorted(set(counts[i] for i in all_systems)):
            page, rows = _everything(conn, place, 20.0, stars=(stars, stars))
            kept = {row["id"] for row in rows if row["kind"] == "system"}
            assert kept == {i for i in all_systems if counts[i] == stars}, stars
            assert all(row["star_count"] == stars for row in rows if row["kind"] == "system")
            assert all(row["system_id"] in kept for row in rows), "what is left is in a kept system"
            assert page["total"] == len(rows)
    finally:
        conn.close()


def test_or_more_and_exactly_split_the_systems_and_a_count_nothing_has_finds_none(galaxy):
    conn = _conn(galaxy)
    try:
        place = near.place_from_point(sector_position_pc(*ADDRESSES[0], EDGE_PC))
        _page, everything = _everything(conn, place, 20.0, kinds=["system"])
        one = _everything(conn, place, 20.0, kinds=["system"], stars=(1, 1))[1]
        more = _everything(conn, place, 20.0, kinds=["system"], stars=(2, None))[1]
        assert {r["ref"] for r in one} | {r["ref"] for r in more} == {r["ref"] for r in everything}
        assert not {r["ref"] for r in one} & {r["ref"] for r in more}
        none = _everything(conn, place, 20.0, stars=(99, None))[1]
        assert none == []
    finally:
        conn.close()


def test_a_filter_leaves_out_what_is_not_in_a_system(galaxy):
    conn = _conn(galaxy)
    try:
        place = near.place_from_point(sector_position_pc(*ADDRESSES[0], EDGE_PC))
        _page, everything = _everything(conn, place, 20.0)
        assert any(row["kind"] == "facility" and row["system_id"] is None for row in everything), "the fixture has a standalone facility"
        _page, filtered = _everything(conn, place, 20.0, stars=(1, None))
        assert not [row for row in filtered if row["system_id"] is None]
    finally:
        conn.close()


def test_a_system_row_carries_its_star_count_even_without_the_filter(galaxy):
    conn = _conn(galaxy)
    try:
        counts = _counts(conn)
        place = near.place_from_point(sector_position_pc(*ADDRESSES[0], EDGE_PC))
        _page, rows = _everything(conn, place, 20.0, kinds=["system"])
        assert rows and all(row["star_count"] == counts[row["id"]] for row in rows)
    finally:
        conn.close()


def test_the_nearest_few_follow_the_filter(galaxy):
    conn = _conn(galaxy)
    try:
        place = near.place_from_point(sector_position_pc(*ADDRESSES[0], EDGE_PC))
        result = near.nearest(conn, place, 3, kinds=["system"], stars=(1, None))
        assert result["rows"] and all(row["star_count"] >= 1 for row in result["rows"])
        assert near.nearest(conn, place, 3, kinds=["system"], stars=(99, None))["rows"] == []
    finally:
        conn.close()


def test_api_near_takes_the_stars_filter_and_refuses_a_bad_one(client, galaxy):  # noqa: F811
    point = ",".join(f"{v:.3f}" for v in sector_position_pc(*ADDRESSES[0], EDGE_PC))
    every = client.get(f"/api/near?point={point}&distance=20&kinds=system&limit=200").get_json()
    kept = client.get(f"/api/near?point={point}&distance=20&kinds=system&limit=200&stars=1%2B").get_json()
    assert kept["total"] == every["total"] and all(row["star_count"] >= 1 for row in kept["rows"])
    none = client.get(f"/api/near?point={point}&distance=20&stars=quintenary").get_json()
    assert none["total"] == 0
    bad = client.get(f"/api/near?point={point}&distance=20&stars=lots")
    assert bad.status_code == 400 and "stars must be" in bad.get_json()["error"]


def test_the_nearby_page_has_the_stars_menu_and_keeps_the_choice(client, galaxy):  # noqa: F811
    conn = _conn(galaxy)
    try:
        system_id = next(iter(_counts(conn)))
    finally:
        conn.close()
    from tests.publicids import pid
    place = f"system:{pid('system', system_id)}"
    asking = client.get(f"/nearby?place={place}").get_data(as_text=True)
    assert 'name="stars"' in asking and all(label in asking for _value, label in star_count.STAR_COUNT_CHOICES)
    listed = client.get(f"/nearby?place={place}&distance=20&kinds=system&stars=1%2B").get_data(as_text=True)
    assert '<option value="1+"' not in listed and '<option value="5+" selected' not in listed
    assert "<th scope=\"col\">Stars</th>" in listed
    ignored = client.get(f"/nearby?place={place}&distance=20&kinds=system&stars=lots")
    assert ignored.status_code == 200, "a hand-edited stars value is ignored, not an error"


# --- the command line ---------------------------------------------------------------------------------------------------------

@pytest.fixture
def cli(galaxy, monkeypatch, capsys):  # noqa: F811
    def run(*argv):
        monkeypatch.setattr(sys, "argv", ["query", *mysql_argv(galaxy), *argv])
        query_cli.main()
        return capsys.readouterr().out
    return run


def test_the_command_line_lists_systems_by_stars_and_says_how_many_each_holds(cli, galaxy):  # noqa: F811
    conn = _conn(galaxy)
    try:
        counts = sorted(set(_counts(conn).values()))
    finally:
        conn.close()
    for stars in counts:
        out = cli("systems", "--stars", str(stars), "--limit", "0")
        label = star_count.star_count_label((stars, stars)).lower()
        lines = [line for line in out.splitlines() if line.startswith("[")]
        assert lines and all(f"({label}," in line for line in lines), out
    assert "No matching systems." in cli("systems", "--stars", "quintenary")


def test_the_command_line_near_takes_the_filter(cli, galaxy):  # noqa: F811
    point = ",".join(f"{v:.3f}" for v in sector_position_pc(*ADDRESSES[0], EDGE_PC))
    out = cli("near", point, "--distance", "20", "--stars", "99+")
    assert "Nothing within" in out
    out = cli("near", point, "--distance", "20", "--stars", "1+", "--kinds", "system")
    assert "(single)" in out or "(binary)" in out or "(trinary)" in out


def test_the_command_line_refuses_a_bad_stars_value_and_names_every_choice_in_its_help(cli, capsys):
    with pytest.raises(SystemExit) as refused:
        cli("systems", "--stars", "lots")
    assert refused.value.code == 2 and "stars must be" in capsys.readouterr().err
    with pytest.raises(SystemExit):
        cli("systems", "--help")
    help_text = capsys.readouterr().out
    for name, _count in star_count.STAR_COUNT_NAMES:
        assert name in help_text
