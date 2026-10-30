# tests/test_query_cli_units.py

"""UX.107: the query command's distances and radii go through the unit ladders, not their own `:.2f` formats."""

import sys

import pytest

from planetgen.cli import query


class _Connection:
    def close(self):
        pass


@pytest.fixture
def cli(monkeypatch, capsys):
    monkeypatch.setattr(query, "open_readonly", lambda config: _Connection())
    monkeypatch.setattr(query, "_printed", lambda conn, kind, row_id: f"{kind}-{row_id}")

    def run(*argv):
        monkeypatch.setattr(sys, "argv", ["query", *argv])
        query.main()
        return capsys.readouterr().out
    return run


def test_a_sector_lists_its_edge_on_the_distance_ladder(cli, monkeypatch):
    monkeypatch.setattr(query, "list_sectors", lambda conn: [{"id": 1, "name": "A", "edge_ly": 12.0, "system_count": 3}])
    assert cli("sectors").strip() == "[sector-1] A (edge 3.68 pc (12 ly), 3 systems)"


def test_a_planet_and_a_moon_list_their_radii_on_the_distance_ladder(cli, monkeypatch):
    row = {"id": 2, "name": "P", "planet_class": "M", "radius_km": 6371.2, "system_name": "S", "planet_name": "Q"}
    monkeypatch.setattr(query, "list_planets", lambda conn, **filters: [row])
    monkeypatch.setattr(query, "list_moons", lambda conn, **filters: [row])
    assert cli("planets").strip() == "[planet-2] P (Class M, 6,371 km) -- S"
    assert cli("moons").strip() == "[moon-2] P (Class M, 6,371 km) -- S / Q"


def test_a_listing_asks_for_one_row_past_its_limit_and_prints_the_after_to_carry_on_with(cli, monkeypatch):
    """PERF.207: `planets` and `moons` show `--limit` rows (100 unless told), and when another row follows print the
    `--after` that carries on from the last shown; `--limit 0` shows every row and `--after` is passed on as a row id."""
    asked = []

    def rows(conn, **filters):
        asked.append(filters)
        return [{"id": n, "name": f"P{n}", "planet_class": "M", "radius_km": 6371.2, "system_name": "S", "planet_name": "Q"}
                for n in range(1, 5)]

    monkeypatch.setattr(query, "list_planets", rows)
    monkeypatch.setattr(query, "list_moons", rows)
    monkeypatch.setattr(query, "_row", lambda conn, kind, printed: None if printed is None else int(printed.split("-")[1]))
    out = cli("planets", "--limit", "3")
    assert asked[-1]["limit"] == 4 and asked[-1]["after"] is None
    assert out.count("[planet-") == 3 and "--after planet-3" in out
    out = cli("moons", "--limit", "5", "--after", "moon-9")
    assert asked[-1]["limit"] == 6 and asked[-1]["after"] == 9
    assert out.count("[moon-") == 4 and "--after" not in out
    out = cli("planets", "--limit", "0")
    assert asked[-1]["limit"] is None and out.count("[planet-") == 4 and "--after" not in out
    cli("planets")
    assert asked[-1]["limit"] == query.DEFAULT_PAGE + 1
    with pytest.raises(SystemExit):
        cli("planets", "--limit", "-1")
