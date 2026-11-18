"""
GEN.1572: the ice-capped flag. A body is ice-capped when ice (sea ice and land ice together) covers more than a quarter of its surface: the larger of the cover its hydrosphere stored and the ice share its own
albedo law (`climate_branches`, CC16) puts on a world of its mean surface temperature. It is derived, shown on the sheet, and changes no class letter.
"""

import pytest

from planetgen.describe import sheet
from planetgen.physics import climate_branches, hydrosphere


def test_the_edge_falls_where_the_ice_ramp_reaches_a_quarter():
    law = climate_branches.IceAlbedoLaw()
    edge = law.free_above_k - (law.free_above_k - law.covered_below_k) * hydrosphere.ICE_CAPPED_MIN_COVER ** (1.0 / law.exponent)
    assert edge == pytest.approx(264.05, abs=0.05)
    assert hydrosphere.ice_capped("surface ocean", 0.0, edge - 0.5)
    assert not hydrosphere.ice_capped("surface ocean", 0.0, edge + 0.5)


@pytest.mark.parametrize("state, cover, temperature, expected", [
    ("surface ocean", 0.0, 288.0, False),        # an Earth analog
    ("surface ocean", 0.0, 250.0, True),         # a cold sea: the law puts about 40 percent of the world under ice
    ("surface ocean", 0.0, 310.0, False),
    ("brine seep", 0.0, 230.0, True),
    ("ice-covered ocean", 1.0, 110.0, True),
    ("ice", 0.6, 120.0, True),
    ("ice-covered ocean", 0.9, 285.0, True),     # a stored cover counts whatever the temperature says
    ("ice-covered ocean", 0.2, 285.0, False),
    ("dry", 0.0, 150.0, False),                  # no surface water, no ice to speak of
    ("vapor", 0.0, 150.0, False),
    (None, None, 150.0, False),                  # a giant has no hydrosphere
    ("surface ocean", 0.0, None, False),
])
def test_the_flag_follows_the_stored_cover_and_the_temperature(state, cover, temperature, expected):
    assert hydrosphere.ice_capped(state, cover, temperature) is expected


def test_the_share_is_the_larger_of_the_stored_cover_and_the_law():
    law = climate_branches.IceAlbedoLaw()
    assert hydrosphere.ice_cover_share("surface ocean", 0.0, 250.0) == pytest.approx(law.ice_fraction(250.0))
    assert hydrosphere.ice_cover_share("ice", 0.9, 250.0) == pytest.approx(0.9)
    assert hydrosphere.ice_cover_share("dry", 0.9, 100.0) == 0.0


class _Row(dict):
    """A stored row as the sheet's providers see it."""


def test_the_sheet_row_appears_only_for_an_ice_capped_body():
    capped = list(sheet._ice_capped_provider(None, "planet", 1, _Row(hydrosphere="surface ocean", ice_cover_fraction=0.0, surface_temperature_k=250.0)))
    assert [(group, row["key"], row["text"]) for group, row in capped] == [("water", sheet.ICE_CAPPED_KEY, "Yes")]
    assert capped[0][1]["raw"] == pytest.approx(climate_branches.IceAlbedoLaw().ice_fraction(250.0))
    assert not list(sheet._ice_capped_provider(None, "planet", 1, _Row(hydrosphere="surface ocean", ice_cover_fraction=0.0, surface_temperature_k=288.0)))
    assert not list(sheet._ice_capped_provider(None, "moon", 1, _Row(hydrosphere=None, ice_cover_fraction=None, surface_temperature_k=100.0)))
    assert not list(sheet._ice_capped_provider(None, "planet", 1, _Row(hydrosphere="ice")))      # a row without the columns gives none


def test_the_provider_serves_planets_and_moons_only():
    kinds = {kind for kind in ("planet", "moon", "star", "asteroid_belt") for applies, provide in sheet.ROW_PROVIDERS if provide is sheet._ice_capped_provider and applies(kind)}
    assert kinds == {"planet", "moon"}
