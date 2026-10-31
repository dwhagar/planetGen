"""
Ozone column against oxygen and host star (GEN.413, `planetgen.physics.ozone`): the Kozakis et al. 2022 fit against the
values the research read from its Table 3, the host mix, the spread and the abiotic branch.
"""

import pytest

from planetgen.physics import ozone
from planetgen.util import draw


@pytest.mark.parametrize("teff, pal, du", [
    (5800.0, 1.0, 192.0), (5800.0, 1.5, 183.0), (5800.0, 0.25, 210.0), (5800.0, 0.1, 199.0), (5800.0, 0.01, 108.0),
    (6000.0, 1.0, 263.0), (6000.0, 0.01, 169.0), (5000.0, 1.0, 160.0), (5000.0, 0.01, 63.0),
    (4500.0, 1.0, 179.0), (4500.0, 0.01, 38.0), (3000.0, 1.0, 162.0), (3000.0, 0.01, 17.0)])
def test_the_fit_reproduces_the_values_read_from_kozakis_table_3(teff, pal, du):
    assert ozone.median_column_du(teff, pal) == pytest.approx(du, abs=1.0)


def test_the_suns_maximum_is_the_papers_5_64e18_at_a_quarter_of_the_oxygen():
    best = max(ozone.median_column_du(5800.0, pal / 100.0) for pal in range(1, 151))
    assert best * ozone.DOBSON_CM2 == pytest.approx(5.64e18, rel=0.01)
    assert ozone.median_column_du(5800.0, 0.25) == pytest.approx(best, rel=0.01)


def test_hot_hosts_peak_at_low_oxygen_and_cool_ones_rise_to_the_highest():
    assert ozone.median_column_du(6000.0, 0.25) > ozone.median_column_du(6000.0, 1.5)
    assert ozone.median_column_du(3000.0, 1.5) > ozone.median_column_du(3000.0, 1.0) > ozone.median_column_du(3000.0, 0.25)


def test_between_hosts_the_column_mixes_in_temperature_and_cooler_than_the_coolest_keeps_its_row():
    low, high = ozone.median_column_du(5000.0, 0.5), ozone.median_column_du(5800.0, 0.5)
    assert ozone.median_column_du(5400.0, 0.5) == pytest.approx((low + high) / 2)
    assert ozone.median_column_du(2500.0, 0.5) == ozone.median_column_du(3000.0, 0.5)
    assert ozone.median_column_du(5772.0, 1.0) == pytest.approx(191.0, abs=1.5)  # the Sun


def test_hotter_than_g0v_scales_to_the_rugheimer_f_star_columns_and_holds_beyond():
    assert ozone.median_column_du(6240.0, 1.0) * ozone.DOBSON_CM2 == pytest.approx(8.62e18, rel=1e-6)
    assert ozone.median_column_du(7000.0, 1.0) * ozone.DOBSON_CM2 == pytest.approx(1.15e19, rel=1e-6)
    assert ozone.median_column_du(9000.0, 1.0) == ozone.median_column_du(7000.0, 1.0)
    assert ozone.median_column_du(6000.0, 1.0) < ozone.median_column_du(6100.0, 1.0) < ozone.median_column_du(6240.0, 1.0)
    assert ozone.median_column_du(6500.0, 0.25) > ozone.median_column_du(6500.0, 1.5)  # the G0V shape in oxygen


def test_the_fit_is_held_at_its_ends_and_never_negative():
    assert ozone.median_column_du(5800.0, 5.0) == ozone.median_column_du(5800.0, ozone.PAL_MAX)
    assert ozone.median_column_du(5800.0, 1e-9) == ozone.median_column_du(5800.0, ozone.PAL_MIN)
    for teff in (3000.0, 4500.0, 5000.0, 5800.0, 6000.0, 6500.0):
        assert all(ozone.median_column_du(teff, pal) >= 0.0 for pal in (1e-9, 1e-4, 0.001, 0.5, 1.0, 1.5))


def test_oxygen_level_follows_the_column_of_oxygen():
    assert ozone.o2_level_pal(21200.0, 1.0) == pytest.approx(1.0)
    assert ozone.o2_level_pal(21200.0, 2.0) == pytest.approx(0.5)
    assert ozone.o2_level_pal(0.0, 1.0) == 0.0 and ozone.o2_level_pal(None, 1.0) == 0.0 and ozone.o2_level_pal(5000.0, 0.0) == 0.0


def test_an_earth_world_draws_within_the_sampled_models_including_the_measured_columns():
    """GEN.748: at the present oxygen the spread is the models' own, 279 DU (Cooke) to 327 DU (Grenfell), and holds the record's 300 DU; the
    earlier 163 DU that no model gives is out."""
    median = ozone.median_column_du(5772.0, 1.0)
    low, high = (median * f for f in ozone.structure_range(1.0))
    assert low == pytest.approx(279.0, abs=3.0) and high == pytest.approx(327.0, abs=3.0) and low < 300.0 < high
    with draw.bound(draw.Stream("ozone-earth")):
        values = [ozone.draw_column_du(5772.0, 21200.0, 1.0) for _ in range(400)]
    assert low - 1e-9 <= min(values) and max(values) <= high + 1e-9
    assert any(279.0 <= v <= 300.0 for v in values) and any(300.0 <= v <= 327.0 for v in values)
    assert min(values) > 190.0


def test_the_spread_is_the_widest_where_no_model_narrows_it_and_the_models_ends_between():
    assert ozone.structure_range(1.0) == (1.45, 1.70) == ozone.structure_range(1.5)
    assert ozone.structure_range(0.1) == (0.85, 1.34)
    assert ozone.structure_range(0.01) == (0.85, 1.71) == ozone.structure_range(0.0002)
    previous = None
    for pal in (0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0):
        low, high = ozone.structure_range(pal)
        assert 0.85 <= low <= high <= 1.71
        if previous:
            assert low >= previous[0] - 1e-12     # the lowest ratio never falls as the oxygen rises
        previous = (low, high)
    # nothing is wider than the earlier range at any level, and the whole range is only where the models stop
    assert all(ozone.structure_range(pal)[1] - ozone.structure_range(pal)[0] <= 1.71 - 0.85 + 1e-12 for pal in (0.01, 0.1, 0.5, 1.0))


def test_the_columns_across_the_five_hosts_follow_the_paper_at_every_temperature_of_the_sweep():
    """The median across the host temperatures of the fit (Kozakis Table 1: 3000, 4500, 5000, 5800, 6000 K) at the present oxygen, and the draw
    stays inside the models' range around each."""
    for teff, du in ((3000.0, 162.0), (4500.0, 179.0), (5000.0, 160.0), (5800.0, 192.0), (6000.0, 263.0)):
        median = ozone.median_column_du(teff, 1.0)
        assert median == pytest.approx(du, abs=1.0)
        with draw.bound(draw.Stream("sweep-%d" % teff)):
            values = [ozone.draw_column_du(teff, 21200.0, 1.0) for _ in range(60)]
        assert all(1.45 * median - 1e-9 <= v <= 1.70 * median + 1e-9 for v in values)


def test_no_oxygen_has_no_ozone_and_draws_nothing():
    with draw.bound(draw.Stream("no-o2")):
        before = draw.random()
    with draw.bound(draw.Stream("no-o2")):
        assert ozone.draw_column_du(5772.0, 0.0, 1.0) == 0.0
        assert ozone.draw_column_du(None, 21200.0, 1.0) == 0.0
        assert draw.random() == before


def test_abiotic_oxygen_gives_a_thousandth_of_the_full_column_or_the_f_share():
    full = ozone.median_column_du(5772.0, 1.0)
    low, high = ozone.abiotic_column_range_du(5772.0)
    assert low == pytest.approx(1e-4 * full) and high == pytest.approx(1e-3 * full)
    f_low, f_high = ozone.abiotic_column_range_du(7000.0)
    assert f_low == f_high == pytest.approx(0.15 * ozone.median_column_du(7000.0, 1.0))
    with draw.bound(draw.Stream("abiotic")):
        values = [ozone.draw_column_du(5772.0, 0.5, 1.0) for _ in range(50)]  # 0.5 Pa of oxygen
    assert all(low <= v <= high for v in values) and max(values) > 2 * min(values)
    with draw.bound(draw.Stream("abiotic")):
        assert ozone.draw_column_du(7000.0, 0.5, 1.0) == pytest.approx(f_low)


def test_the_draw_is_reproducible_and_one_draw():
    with draw.bound(draw.Stream("repeat")):
        first = ozone.draw_column_du(5772.0, 21200.0, 1.0)
        after = draw.random()
    with draw.bound(draw.Stream("repeat")):
        assert ozone.draw_column_du(5772.0, 21200.0, 1.0) == first
        assert draw.random() == after
    with draw.bound(draw.Stream("repeat")):
        draw.uniform(0, 1)
        assert draw.random() == after


