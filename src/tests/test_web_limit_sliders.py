# tests/test_web_limit_sliders.py

"""UX.116: the Generate page's limit sliders (web/limit_sliders.py), without a page."""

import math

import pytest

from planetgen import tuning
from planetgen.web import limit_sliders
from planetgen.web.limit_sliders import SLIDERS, LimitSlider, ladder, plain, steps


def test_a_ladder_steps_one_two_five_from_end_to_end():
    assert ladder(100.0, 1000.0) == (100.0, 200.0, 500.0, 1000.0)
    assert ladder(50.0, 1000.0) == (50.0, 100.0, 200.0, 500.0, 1000.0)
    assert ladder(1.0, 10.0) == (1.0, 2.0, 5.0, 10.0)
    wide = ladder(1.0, 1.0e8)
    assert wide[0] == 1.0 and wide[-1] == 1.0e8 and list(wide) == sorted(set(wide))


def test_steps_are_even_and_exact():
    assert steps(8.0, 24.0, 0.2)[0] == 8.0 and steps(8.0, 24.0, 0.2)[-1] == 24.0
    assert len(steps(8.0, 24.0, 0.2)) == 81
    assert 15.4 in steps(8.0, 24.0, 0.2)


def test_plain_writes_a_number_without_an_exponent():
    assert plain(4_000_000.0) == "4000000" and plain(15.4) == "15.4" and plain(9900.000000001) == "9900"


@pytest.mark.parametrize("slider", list(SLIDERS.values()), ids=list(SLIDERS))
def test_every_slider_has_ascending_stops_and_a_position_for_any_value(slider):
    assert list(slider.stops) == sorted(set(slider.stops))
    assert slider.position(slider.low) == 0 and slider.position(slider.high) == len(slider.stops) - 1
    assert slider.position(slider.low / 10) == 0 and slider.position(slider.high * 10) == len(slider.stops) - 1
    for index, stop in enumerate(slider.stops):
        assert slider.position(stop) == index
    if slider.default is not None:
        assert slider.low <= slider.default <= slider.high


def test_a_log_slider_picks_the_nearer_stop_by_ratio_and_a_linear_one_by_difference():
    log = LimitSlider("x", "X", "log", (100.0, 1000.0))
    assert log.position(300.0) == 0 and log.position(320.0) == 1       # the geometric middle is about 316
    linear = LimitSlider("y", "Y", "linear", (100.0, 1000.0))
    assert linear.position(540.0) == 0 and linear.position(560.0) == 1


def test_the_defaults_come_from_the_tunings():
    assert SLIDERS["phenomenon_min_mass"].default == tuning.PHENOMENON_MIN_MASS_SOLAR
    assert SLIDERS["bright_min_luminosity"].default == tuning.BRIGHT_STAR_MIN_LUMINOSITY_SOL
    assert SLIDERS["halo_min_luminosity"].default == tuning.HALO_MIN_LUMINOSITY_SOL
    assert SLIDERS["cell_star_ceiling"].default == tuning.CELL_STAR_CEILING
    assert SLIDERS["fill_luminosity_cap"].default == tuning.FILL_LUMINOSITY_CAP_SOL
    assert SLIDERS["core_min_luminosity"].default is None and SLIDERS["core_mass_limit"].default is None


def test_the_ranges_match_what_the_settings_accept():
    halo = SLIDERS["halo_min_luminosity"]
    assert (halo.low, halo.high) == tuning.HALO_MIN_LUMINOSITY_RANGE_SOL
    for name in ("fill_luminosity_cap", "core_min_luminosity"):
        assert SLIDERS[name].low == tuning.WD_LUMINOSITY_RANGE_SOL[1]
    assert SLIDERS["bright_min_luminosity"].low == tuning.BRIGHT_STAR_FLOOR_MIN_SOL
    assert SLIDERS["bright_min_luminosity"].high == tuning.BRIGHT_STAR_FLOOR_MAX_SOL
    assert SLIDERS["phenomenon_min_mass"].snap and not SLIDERS["bright_min_luminosity"].snap


def test_a_long_log_ladder_marks_only_its_decades():
    ceiling = SLIDERS["cell_star_ceiling"]
    marks = [ceiling.stops[index] for index in ceiling.tick_positions]
    assert marks == [10.0 ** exponent for exponent in range(0, 9)]
    assert len(SLIDERS["phenomenon_min_mass"].tick_positions) == 7
