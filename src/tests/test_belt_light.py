# tests/test_belt_light.py

"""
GEN.1444: belts and comets read the radiation field (`physics/radiation_field.py`, GEN.1441) for the place they are at: a belt takes the light averaged over its ring for the
temperature its bodies' makeup goes by, the taxa of its members and the side of the snow line it is on (its size law); a comet takes the light at its perihelion for its activity.
A belt or a comet of a lone star must give the numbers of before the change to the last digit.
"""

import math

import pytest

from planetgen.generation import belt as belt_generation
from planetgen.generation import comet as comet_generation
from planetgen.generation.belt import AsteroidBelt
from planetgen.generation.comet import Comet
from planetgen.generation.system import StarSystem
from planetgen.physics import belt_ledger, condensate_shares, radiation_field as rf
from planetgen.util import draw
from tests.test_radiation_field import star as plain_star
from tests.test_systems import make_config


@pytest.fixture()
def system():
    return StarSystem(system_config=make_config("G2V", PLANETS=False))


@pytest.fixture()
def host(system):
    """The system's star, given back with no companions however a test left it."""
    star = system.star
    rf.attach(star, ())
    yield star
    rf.attach(star, ())


def companion_of(host, axis_au, eccentricity=0.0, luminosity_w=None):
    rf.attach(host, [rf.Companion(plain_star(luminosity_w or host.luminosity), 1, axis_au, eccentricity)])
    return host


def brute_ratio(host, lower, upper, annuli=2000):
    """The ring's mean total light over its mean host light, area-weighted, on a fine grid."""
    own = total = 0.0
    width = (upper - lower) / annuli
    for index in range(annuli):
        distance = lower + (index + 0.5) * width
        own += distance * rf.starlight_w_m2(rf.Companionless(host), distance)
        total += distance * rf.starlight_w_m2(host, distance)
    return total / own


def test_the_ring_of_a_lone_star_has_a_light_ratio_of_exactly_one(host):
    assert rf.ring_light_ratio(host, 2.2, 3.4) == 1.0
    assert rf.lit_distance_au(host, 2.5) == 2.5
    assert belt_generation.lit_luminosity_w(host, 2.2, 3.4) is host.luminosity


def test_the_ring_ratio_is_the_mean_light_over_the_ring(host):
    companion_of(host, 8.0, 0.2)
    assert rf.ring_light_ratio(host, 2.2, 3.4) == pytest.approx(brute_ratio(host, 2.2, 3.4), rel=2e-3)
    inner, outer = rf.light_ratio(host, 2.2), rf.light_ratio(host, 3.4)
    assert inner < rf.ring_light_ratio(host, 2.2, 3.4) < outer, "the companion's share grows outward, so the ring's mean lies between its edges"
    assert rf.ring_light_ratio(host, 2.5, 2.5) == rf.light_ratio(host, 2.5)
    assert rf.ring_light_ratio(host, 2.5, 2.5001) == pytest.approx(rf.light_ratio(host, 2.5), rel=1e-4)


def test_a_lit_distance_is_the_distance_of_a_lone_star_that_lights_as_much(host):
    companion_of(host, 6.0)
    lit = rf.lit_distance_au(host, 2.0)
    assert lit < 2.0
    assert rf.starlight_w_m2(rf.Companionless(host), lit) == pytest.approx(rf.starlight_w_m2(host, 2.0), rel=1e-9)


def test_a_belt_of_a_lone_star_is_unchanged(system, host):
    with draw.bound(5):
        belt = AsteroidBelt(system.system_config, host, 2.5, 2.2, 3.4)
    assert belt.makeup.detail["temperature_k"] == condensate_shares.disk_temperature_k(2.5, host.luminosity)
    assert belt_generation.size_law(2.5, host.snow_line_au) == (belt_ledger.ICY_LAW if 2.5 >= host.snow_line_au else belt_ledger.ROCKY_LAW)


def test_a_belt_close_to_a_companion_is_as_warm_as_the_ring_is_lit(system, host):
    with draw.bound(5):
        lone = AsteroidBelt(system.system_config, host, 2.5, 2.2, 3.4)
    companion_of(host, 7.0)
    ratio = rf.ring_light_ratio(host, 2.2, 3.4)
    assert ratio > 1.1
    with draw.bound(5):
        warm = AsteroidBelt(system.system_config, host, 2.5, 2.2, 3.4)
    assert warm.makeup.detail["temperature_k"] == pytest.approx(lone.makeup.detail["temperature_k"] * ratio ** 0.25, rel=1e-9)
    assert warm.makeup.detail["temperature_k"] > lone.makeup.detail["temperature_k"]


def test_a_belt_far_from_its_companion_is_barely_moved(system, host):
    with draw.bound(5):
        lone = AsteroidBelt(system.system_config, host, 2.5, 2.2, 3.4)
    companion_of(host, 2000.0)
    with draw.bound(5):
        far = AsteroidBelt(system.system_config, host, 2.5, 2.2, 3.4)
    assert far.makeup.detail["temperature_k"] == pytest.approx(lone.makeup.detail["temperature_k"], rel=1e-6)
    assert far.mass_kg == lone.mass_kg


def test_a_belt_just_past_the_snow_line_is_rocky_when_a_companion_lights_it_like_one_inside(system, host):
    snow = host.snow_line_au
    distance = 1.05 * snow
    assert belt_generation.size_law(distance, snow) is belt_ledger.ICY_LAW
    companion_of(host, 2.5 * snow)
    ratio = rf.ring_light_ratio(host, 0.95 * distance, 1.05 * distance)
    assert distance / math.sqrt(ratio) < snow
    assert belt_generation.size_law(distance, snow, ratio) is belt_ledger.ROCKY_LAW
    with draw.bound(5):
        belt = AsteroidBelt(system.system_config, host, distance, 0.95 * distance, 1.05 * distance)
    assert belt.largest_body_km == belt_generation._present_state(host.mass, snow, host.age, belt_generation.bulk_draw.metal_ratio_to_solar(host), distance,
                                                                  0.95 * distance, 1.05 * distance, None, ratio)[0]


def test_the_members_of_a_belt_are_drawn_under_the_same_light_as_its_makeup(system, host):
    companion_of(host, 7.0)
    expected = host.luminosity * rf.ring_light_ratio(host, 2.2, 3.4)
    assert belt_generation.lit_luminosity_w(host, 2.2, 3.4) == expected
    with draw.bound(5):
        belt = AsteroidBelt(system.system_config, host, 2.5, 2.2, 3.4)
    assert belt.makeup.detail["temperature_k"] == pytest.approx(condensate_shares.disk_temperature_k(2.5, expected), rel=1e-9)


@pytest.fixture()
def seen_activity(monkeypatch):
    """The perihelion distance each comet's activity chance was asked for."""
    asked = []
    real = comet_generation._activity_chance

    def watch(perihelion_distance_au):
        asked.append(perihelion_distance_au)
        return real(perihelion_distance_au)

    monkeypatch.setattr(comet_generation, "_activity_chance", watch)
    return asked


def test_a_comet_of_a_lone_star_is_active_by_its_perihelion(system, host, seen_activity):
    with draw.bound(8):
        comet = Comet(system.system_config, 1.0, orbit_type="elliptical", star=host)
    assert seen_activity == [comet.perihelion_distance_au]


def test_a_comet_is_active_by_the_light_at_its_perihelion(system, host, seen_activity):
    with draw.bound(8):
        lone = Comet(system.system_config, 1.0, orbit_type="elliptical", star=host)
    companion_of(host, 4.0)
    with draw.bound(8):
        near = Comet(system.system_config, 1.0, orbit_type="elliptical", star=host)
    assert near.perihelion_distance_au == lone.perihelion_distance_au
    assert seen_activity[1] == pytest.approx(rf.lit_distance_au(host, near.perihelion_distance_au), rel=1e-12)
    assert seen_activity[1] < seen_activity[0]
    assert comet_generation._activity_chance(seen_activity[1]) >= comet_generation._activity_chance(seen_activity[0])

