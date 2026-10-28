"""
DOC.82: the region-name audit (`tests.audit.region_names`): around several points of the galaxy and under several naming
keys no region shares a name, a family word or an ordinal with a neighbour, and the check does find a clash when one is made.
"""

import pytest

from planetgen.galaxy import named_regions
from tests.audit import region_names

KEYS = ("DEADBEEF", "00000001", "C0FFEE42")
POINTS = ((0.0, 0.0, 0.0), (8000.0, 0.0, 0.0), (-3000.0, 6000.0, 400.0), (100.0, -100.0, 25000.0))


@pytest.mark.parametrize("key", KEYS)
@pytest.mark.parametrize("point", POINTS)
def test_regions_around_a_point_stay_unique_among_neighbours(point, key):
    cells = region_names.cells_around(point, radius_cells=1)
    assert len(cells) >= 8
    assert region_names.findings(cells, key) == []


def test_a_cell_has_neighbours_and_never_itself():
    cell = region_names.cell_at((8000.0, 0.0, 0.0))
    around = region_names.neighbour_cells(cell)
    assert 6 <= len(around) <= 26 and cell not in around


def test_the_check_finds_two_neighbours_given_one_name(monkeypatch):
    cells = region_names.cells_around((8000.0, 0.0, 0.0), radius_cells=1)
    monkeypatch.setattr(named_regions, "region_name", lambda block, key, edge_pc=4: "Same Name Region")
    lines = region_names.findings(cells, "DEADBEEF")
    assert lines and all(line.startswith("same name") for line in lines)


def test_the_check_finds_a_family_word_shared_by_neighbouring_families(monkeypatch):
    cells = region_names.cells_around((8000.0, 0.0, 0.0), radius_cells=1)
    monkeypatch.setattr(named_regions, "region_name",
                        lambda block, key, edge_pc=4: f"Same {named_regions.ordinal_in_family(block)} Region")
    lines = region_names.findings(cells, "DEADBEEF")
    assert any(line.startswith("same family word") for line in lines)
