"""GEN.1491: a comet's activity reads the host's luminosity (`planetgen/physics/comet_activity.py`, `Comet`)."""

import types

import pytest

from planetgen import tuning
from planetgen.generation.comet import Comet, _activity_chance
from planetgen.generation.config import SystemConfig
from planetgen.physics import comet_activity, constants
from planetgen.util import draw

SUN = constants.SOLAR_LUMINOSITY


@pytest.fixture(autouse=True)
def _keep_draw_state():
    state = draw.getstate()
    yield
    draw.setstate(state)


def _star(luminosity):
    # No abundances and no companions: the comet reads only the luminosity of it.
    return types.SimpleNamespace(luminosity=luminosity, snow_line_au=2.7)


def _active_share(luminosity, trials=400):
    star = _star(luminosity)
    active = 0
    for seed in range(trials):
        draw.set_run_seed(seed)
        comet = Comet(SystemConfig(), 1.0, orbit_type="elliptical", star=star)
        active += comet.is_active
    return active / trials


def test_the_sun_and_no_star_give_the_distance_they_were_given():
    for au in (0.1, 1.0, 3.0, 9.0):
        assert comet_activity.solar_equivalent_au(au, SUN) == pytest.approx(au)
        assert comet_activity.solar_equivalent_au(au, None) == au
        assert comet_activity.solar_equivalent_au(au, 0.0) == au


def test_a_bright_star_lights_a_comet_from_further_out_in_proportion_to_the_square_root_of_its_luminosity():
    assert comet_activity.solar_equivalent_au(30.0, 100.0 * SUN) == pytest.approx(3.0)
    assert comet_activity.solar_equivalent_au(0.3, 0.01 * SUN) == pytest.approx(3.0)


def test_a_comet_at_the_activity_distance_has_the_floor_chance_whatever_the_star():
    for luminosity in (0.01 * SUN, SUN, 400.0 * SUN):
        at = tuning.COMET_ACTIVITY_PERIHELION_THRESHOLD_AU * (luminosity / SUN) ** 0.5
        floor = _activity_chance(comet_activity.solar_equivalent_au(at, luminosity))
        assert floor == pytest.approx(tuning.COMET_ACTIVITY_MIN_CHANCE)


def test_the_comet_reads_its_hosts_luminosity():
    seen = []
    real = comet_activity.solar_equivalent_au
    comet_activity.solar_equivalent_au = lambda lit, luminosity: seen.append((lit, luminosity)) or real(lit, luminosity)
    try:
        draw.set_run_seed(3)
        comet = Comet(SystemConfig(), 1.0, orbit_type="elliptical", star=_star(5.0 * SUN))
    finally:
        comet_activity.solar_equivalent_au = real
    assert seen == [(pytest.approx(comet.perihelion_distance_au), 5.0 * SUN)]


def test_a_bright_host_makes_more_of_its_comets_active_and_a_faint_one_fewer():
    faint, sun, bright = (_active_share(factor * SUN) for factor in (0.01, 1.0, 100.0))
    assert faint < sun < bright
