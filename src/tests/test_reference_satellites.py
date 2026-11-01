# tests/test_reference_satellites.py

"""
GEN.366: the reference table of the Milky Way's satellites is checked against [R2972]: position and distance from
Table 2, half-light radius (pc, semi-major axis) from Table 3, stellar mass from Table 4 (a stellar mass-to-light ratio of 1) and mean [Fe/H] from Table 5.
The numbers below are transcribed from the paper's tables, not from the table under test.
"""

import pytest

from tests.reference_satellites import MILKY_WAY_SATELLITES

# name: (l, b, D kpc, M* solar, r_h pc or None, [Fe/H])
PAPER = {
    "Large Magellanic Cloud": (280.5, -32.9, 51.0, 1.5e9, None, -0.5),
    "Small Magellanic Cloud": (302.8, -44.3, 64.0, 4.6e8, None, -1.00),
    "Sagittarius dwarf": (5.6, -14.2, 26.0, 2.1e7, 2587.0, -0.4),
    "Fornax dwarf": (237.1, -65.7, 147.0, 2.0e7, 710.0, -0.99),
    "Sculptor dwarf": (287.5, -83.2, 86.0, 2.3e6, 283.0, -1.68),
    "Carina dwarf": (260.1, -22.2, 105.0, 3.8e5, 250.0, -1.72),
    "Draco dwarf": (86.4, 34.7, 76.0, 2.9e5, 221.0, -1.93),
    "Ursa Minor dwarf": (105.0, 44.8, 76.0, 2.9e5, 181.0, -2.13),
    "Sextans dwarf": (243.5, 42.3, 86.0, 4.4e5, 695.0, -1.93),
    "Leo I": (226.0, 49.1, 254.0, 5.5e6, 251.0, -1.43),
    "Leo II": (220.2, 67.2, 233.0, 7.4e5, 176.0, -1.62),
    "Canes Venatici I": (74.3, 79.8, 218.0, 2.3e5, 564.0, -1.98),
}


def test_every_row_is_a_row_of_the_paper_and_the_paper_has_no_other_row_here():
    assert {row.name for row in MILKY_WAY_SATELLITES} == set(PAPER)


@pytest.mark.parametrize("row", MILKY_WAY_SATELLITES, ids=lambda row: row.name)
def test_the_row_matches_the_paper(row):
    longitude, latitude, distance, mass, radius, fe_h = PAPER[row.name]
    assert row.longitude_deg == pytest.approx(longitude, abs=0.05)
    assert row.latitude_deg == pytest.approx(latitude, abs=0.05)
    assert row.distance_kpc == pytest.approx(distance, abs=0.05)
    assert row.stellar_mass_sol == pytest.approx(mass, rel=0.01)
    if radius is not None:                       # the paper gives no half-light radius for the Magellanic Clouds
        assert row.half_light_radius_pc == pytest.approx(radius, abs=0.5)
    assert row.fe_h == pytest.approx(fe_h, abs=0.005)


def test_the_table_holds_the_twelve_satellites_above_a_hundred_thousand_solar_masses_and_not_canis_major():
    """TEST.507: [R2972] lists 27 satellites within 300 kpc, 12 above 1e5 solar masses (Canis Major, 4.9e7, a disputed disk substructure, is not one of the table's rows)."""
    assert len(MILKY_WAY_SATELLITES) == 12
    assert all(row.stellar_mass_sol > 1.0e5 for row in MILKY_WAY_SATELLITES)
    assert "Canis Major" not in {row.name for row in MILKY_WAY_SATELLITES}
