# tests/test_population_slot.py

"""
MAP.318: population is a metric slot. A world with a population ranks by it, before every world without one, and the worlds without
one rank by habitability score (so the inhabited slice is the habitable worlds until head counts are stored, POP.30); the
charted-regions list carries `population` on a region only when its people are counted.
"""

import pytest

from planetgen.db import population_slot, region_stats
from planetgen.galaxy import drill
from tests.test_charted_regions import client, _stored  # noqa: F401 -- fixtures and helper
from tests.test_region_wire_names import galaxy, keyed  # noqa: F401 -- fixtures


def _names(worlds):
    return [world["name"] for world in population_slot.rank(worlds)]


def test_a_populated_fixture_ranks_by_population():
    worlds = [{"name": "Small", "population": 2.0e6, "habitability": 9.0},
              {"name": "Big", "population": 8.0e9, "habitability": 4.0},
              {"name": "Mid", "population": 3.0e8, "habitability": 6.0}]
    assert _names(worlds) == ["Big", "Mid", "Small"]


def test_an_unpopulated_habitable_fixture_ranks_by_habitability_score():
    worlds = [{"name": "Poor", "habitability": 2.5}, {"name": "Best", "habitability": 8.0}, {"name": "Fair", "habitability": 5.0},
              {"name": "Unscored"}]
    assert _names(worlds) == ["Best", "Fair", "Poor", "Unscored"]


def test_a_populated_world_outranks_a_better_habitable_world_without_people():
    worlds = [{"name": "Empty paradise", "habitability": 9.9}, {"name": "Colony", "population": 1200, "habitability": 3.0}]
    assert _names(worlds) == ["Colony", "Empty paradise"]


def test_no_zero_missing_or_nan_value_is_a_population():
    for value in (None, 0, -5, float("nan")):
        assert population_slot.population_of({"population": value}) is None
    assert population_slot.population_of({}) is None
    assert population_slot.population_of({"population": 7}) == 7.0
    worlds = [{"name": "B", "population": 0, "habitability": 5.0}, {"name": "A", "population": float("nan"), "habitability": 6.0}]
    assert _names(worlds) == ["A", "B"]


def test_ties_do_not_depend_on_the_order_the_worlds_came_in():
    worlds = [{"name": "Beta", "id": 2, "population": 10}, {"name": "Alpha", "id": 3, "population": 10},
              {"name": "Alpha", "id": 1, "population": 10}]
    forward = [(w["name"], w["id"]) for w in population_slot.rank(worlds)]
    backward = [(w["name"], w["id"]) for w in population_slot.rank(list(reversed(worlds)))]
    assert forward == backward == [("Alpha", 1), ("Alpha", 3), ("Beta", 2)]


def test_no_region_has_a_population_while_nothing_stores_head_counts(galaxy, client, keyed):
    conn_regions = client.get("/api/galaxy/charted-regions").get_json()["regions"]
    assert conn_regions and all("population" not in entry for entry in conn_regions)


def test_a_region_whose_people_are_counted_carries_its_population(galaxy, client, keyed, monkeypatch):
    stored = _stored(galaxy)
    cell = stored[0]["block"]
    monkeypatch.setattr(population_slot, "region_totals",
                        lambda conn: {(cell.shell, cell.band, cell.wedge): 4.2e9, (-1, -1, -1): 5.0})
    listed = client.get("/api/galaxy/charted-regions").get_json()["regions"]
    with_people = [entry for entry in listed if "population" in entry]
    assert [entry["key"] for entry in with_people] == [drill.format_drill_key(cell)]
    assert with_people[0]["population"] == 4.2e9
    assert [entry["key"] for entry in listed] == [drill.format_drill_key(r["block"]) for r in stored]       # the order is the share's still
