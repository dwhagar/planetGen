# tests/test_trojan_swarm.py

"""GEN.1150 slice 0: the two regions of a trojan swarm and the libration overlay (physics/trojan_swarm.py), pure and unwired."""

import math

import pytest

from planetgen.physics import trojan_swarm as ts
from planetgen.physics import trojans
from planetgen.util import draw

SUN, JUPITER = 1.0, 9.546e-4
SWARM = 1.0e-9       # in solar masses


def regions(seed, trojan_mass=SWARM):
    return ts.draw_regions(SUN, JUPITER, trojan_mass, draw.Stream(f"swarm:{seed}"))


def test_a_swarm_gets_a_leading_and_a_trailing_region_that_share_its_mass():
    lead, trail = regions(1)
    assert (lead.point, trail.point) == ("L4", "L5")
    assert lead.mass + trail.mass == pytest.approx(SWARM)
    assert trojans.LEADING_SHARE_RANGE[0] <= lead.mass / SWARM <= trojans.LEADING_SHARE_RANGE[1]


def test_every_region_is_a_tadpole_inside_the_orbit_limits():
    for seed in range(300):
        for region in regions(seed):
            assert 0.0 <= region.libration_amplitude_deg < trojans.TADPOLE_LIMIT_DEG
            assert 0.0 <= region.libration_phase_deg < 360.0
            assert trojans.orbit_allowed(region.eccentricity, region.inclination_deg)
            assert trojans.motion(region.libration_amplitude_deg, trojans.mass_parameter(SUN, JUPITER)) == "tadpole"


def test_the_draw_is_repeatable_and_takes_exactly_nine_draws():
    assert regions(5) == regions(5)

    class Counting:
        def __init__(self):
            self.count, self.stream = 0, draw.Stream("count")

        def random(self):
            self.count += 1
            return self.stream.random()

    counter = Counting()
    ts.draw_regions(SUN, JUPITER, SWARM, counter)
    assert counter.count == 9
    counter = Counting()
    ts.draw_regions(SUN, JUPITER, 5 * SWARM, counter)
    assert counter.count == 9          # the mass of the swarm does not change how much is drawn


def test_a_pair_that_cannot_hold_the_swarm_or_a_swarm_with_no_mass_gives_none():
    assert ts.draw_regions(SUN, 0.05, SWARM, draw.Stream("a")) == ()                 # a 5 percent companion is past the limit
    assert ts.draw_regions(SUN, JUPITER, 0.0, draw.Stream("a")) == ()
    assert ts.draw_regions(SUN, JUPITER, 0.1, draw.Stream("a")) == ()                # a swarm that would itself push the pair past it


def test_the_truncated_rayleigh_stays_under_its_limit_and_follows_the_scale():
    assert ts.truncated_rayleigh(0.0, 10.0, 24.0) == 0.0
    assert ts.truncated_rayleigh(1.0 - 1e-12, 10.0, 24.0) == pytest.approx(24.0, rel=1e-4)
    values = sorted(ts.truncated_rayleigh(i / 1000.0, 10.0, 24.0) for i in range(1000))
    assert values[500] == pytest.approx(10.0 * math.sqrt(-2.0 * math.log(1.0 - 0.5 * (1.0 - math.exp(-0.5 * 2.4 ** 2)))), rel=1e-2)


def test_the_overlay_swings_about_sixty_degrees_with_the_libration_period_of_the_pair():
    mu = trojans.mass_parameter(SUN, JUPITER)
    region = ts.Region("L4", SWARM, 15.0, 0.0, 0.05, 10.0)
    period = trojans.libration_period(11.86, mu)
    assert period == pytest.approx(148.0, rel=0.05)                                  # Jupiter's: about 150 years
    assert ts.longitude_offset_deg(region, 0.0, 11.86, mu) == pytest.approx(60.0)
    assert ts.longitude_offset_deg(region, period / 4.0, 11.86, mu) == pytest.approx(75.0)
    assert ts.longitude_offset_deg(region, period, 11.86, mu) == pytest.approx(60.0, abs=1e-9)
    trailing = ts.Region("L5", SWARM, 15.0, 90.0, 0.05, 10.0)
    assert ts.longitude_offset_deg(trailing, 0.0, 11.86, mu) == pytest.approx(-60.0 + 15.0)
    swing = [ts.longitude_offset_deg(region, year, 11.86, mu) for year in range(0, 400, 3)]
    assert min(swing) >= 60.0 - 15.0 - 1e-9 and max(swing) <= 60.0 + 15.0 + 1e-9


def test_the_longitude_is_the_secondarys_plus_the_offset_and_a_pair_that_holds_nothing_has_none():
    mu = trojans.mass_parameter(SUN, JUPITER)
    region = ts.Region("L5", SWARM, 0.0, 0.0, 0.0, 0.0)
    assert ts.longitude_deg(region, 100.0, 7.0, 11.86, mu) == pytest.approx(40.0)
    assert ts.longitude_deg(region, 10.0, 7.0, 11.86, mu) == pytest.approx(310.0)
    with pytest.raises(trojans.TrojanError):
        ts.longitude_offset_deg(region, 0.0, 1.0, 0.05)
