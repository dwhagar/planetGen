"""GEN.462: ground albedo by band and space weathering (`planetgen/physics/ground_albedo.py`)."""

import pytest

from planetgen.physics import ground_albedo as ga
from planetgen.util import draw


def test_a_pure_ground_type_is_its_table_row():
    albedo = ga.band_albedo({"basalt_powder": 1.0})
    assert albedo.values == pytest.approx((0.207, 0.214, 0.224))
    assert ga.row_albedo("calcite").values == pytest.approx((0.911, 0.941, 0.946))
    assert ga.row_albedo("kerogen").visible_mean == pytest.approx((0.013 + 0.016 + 0.021) / 3.0)


def test_the_cover_is_the_area_weighted_mean_and_scales_its_shares():
    mixed = ga.band_albedo({"basalt_powder": 3.0, "water_ice": 1.0})
    assert mixed.b2 == pytest.approx(0.75 * 0.214 + 0.25 * 0.832)
    assert mixed.values == pytest.approx(ga.band_albedo({"basalt_powder": 0.75, "water_ice": 0.25}).values)


def test_a_frost_cover_brightens_and_weathering_darkens():
    bare = ga.band_albedo({"basalt_powder": 1.0})
    frosted = ga.band_albedo({"basalt_powder": 0.8, "water_ice": 0.2})
    assert all(f > b for f, b in zip(frosted.values, bare.values))
    old = ga.weather(bare, 4.0, 1.5, 0.5)
    assert all(o < b for o, b in zip(old.values, bare.values))
    assert ga.weather(bare, 0.0, 1.5, 0.5).values == pytest.approx(bare.values)
    assert ga.weathering_factor(1e6, 1.5, 0.5) == pytest.approx(0.5)


def test_an_ocean_comes_in_as_an_explicit_cover_and_a_missing_albedo_is_an_error():
    wet = ga.band_albedo({"basalt_powder": 0.3, "ocean": 0.7}, extra={"ocean": (0.07, 0.06, 0.05)})
    assert wet.b1 == pytest.approx(0.3 * 0.207 + 0.7 * 0.07)
    with pytest.raises(ga.MissingAlbedo, match="co2_ice"):
        ga.band_albedo({"co2_ice": 1.0})
    with pytest.raises(ValueError):
        ga.band_albedo({})
    with pytest.raises(ValueError):
        ga.band_albedo({"basalt_powder": -1.0, "water_ice": 2.0})
    with pytest.raises(ValueError):
        ga.weathering_factor(-1.0, 1.0, 0.5)


def test_the_weathering_draw_stays_in_its_ranges_and_repeats():
    with draw.bound(1):
        values = [ga.draw_weathering() for _ in range(100)]
    assert all(1.0 <= tau <= 2.0 and 0.3 <= d <= 0.7 for tau, d in values)
    with draw.bound(1):
        assert [ga.draw_weathering() for _ in range(100)] == values
