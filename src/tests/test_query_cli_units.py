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
