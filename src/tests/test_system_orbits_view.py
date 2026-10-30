# tests/test_system_orbits_view.py

"""MAP.213: the rows of a built system's "Orbits of the stars" panel (`web/system_orbits_view.py`), from the detail's `orbits` and
the `parent_node_no` and `slot` of its stars."""

import pytest

from planetgen.web.system_orbits_view import orbit_rows

AU_KM = 149_597_870.7


def _node(node_no, parent, slot, kind, mass, au, ecc, period, interaction="detached"):
    return {"node_no": node_no, "parent_node_no": parent, "slot": slot, "kind": kind, "mass_solar": mass,
            "separation_km": au * AU_KM, "eccentricity": ecc, "inclination_deg": 0.0, "period_years": period,
            "interaction": interaction, "flags": 0}


def _star(name, node, slot):
    return {"name": name, "parent_node_no": node, "slot": slot}


TRIPLE = {
    "stars": [_star("Alpha A", 0, 0), _star("Alpha B", 0, 1), _star("Alpha C", 1, 1)],
    "orbits": [_node(0, 1, 0, "close", 1.5, 0.3, 0.0, 0.1, "contact"), _node(1, None, 0, "wide", 2.0, 400.0, 0.3, 4000.0)],
}


def test_a_triple_lists_each_orbit_with_what_it_joins():
    rows = orbit_rows(TRIPLE)
    assert [row[0] for row in rows] == ["1", "2"]
    assert rows[0][1] == "Alpha A and Alpha B" and rows[0][2] == "Close pair"
    assert rows[1][1] == "orbit 1 and Alpha C" and rows[1][2] == "Wide pair"


def test_a_row_reads_the_orbit_in_the_pages_units():
    inner, outer = orbit_rows(TRIPLE)
    assert inner[6] == "0.00" and outer[6] == "0.30"
    assert inner[7] == "Contact" and outer[7] == "Detached"
    assert all(cell for row in (inner, outer) for cell in row)


def test_a_place_nothing_hangs_in_is_left_out_of_the_joins():
    detail = {"stars": [_star("Solo A", 0, 0)], "orbits": [_node(0, None, 0, "wide", 1.0, 30.0, 0.1, 90.0)]}
    assert orbit_rows(detail)[0][1] == "Solo A"


@pytest.mark.parametrize("detail", [{"stars": [], "orbits": []}, {"stars": [_star("Lone", None, None)]}, {}])
def test_a_single_star_or_a_system_with_no_stored_tree_lists_no_orbits(detail):
    assert orbit_rows(detail) == []
