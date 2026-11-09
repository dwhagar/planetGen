# tests/test_radiation_field.py

"""GEN.1441: the radiation field of a system (`planetgen.physics.radiation_field`): every radiating object's part of the light on a body, by position."""

import math
from types import SimpleNamespace

import numpy as np
import pytest

from planetgen.physics import constants, planet_zones, radiation_field as rf

SUN_W = constants.SOLAR_LUMINOSITY


def star(luminosity_w=SUN_W, **extra):
    return SimpleNamespace(luminosity=luminosity_w, temperature=5772.0, mass=1.989e30, radius=6.957e8, **extra)


def with_companion(host, companion, axis_au, eccentricity=0.0):
    rf.attach(host, [rf.Companion(companion, 1, axis_au, eccentricity)])
    return host


def old_light_w_m2(luminosity_w, distance_au):
    """The expression `planets.calculate_atmospheric_conditions` has always used for the host's light."""
    km = distance_au * constants.AU_TO_KM
    return (luminosity_w / (4 * math.pi * km ** 2)) / 1e6


def test_a_lone_star_gives_the_host_light_to_the_last_digit():
    host = star()
    assert rf.companions_of(host) == ()
    for distance in (0.3, 1.0, 1.52, 38.0):
        assert rf.starlight_w_m2(host, distance) == old_light_w_m2(SUN_W, distance)
    assert rf.companion_share(host, 1.0) == 0.0


def test_an_equal_companion_at_five_au_adds_about_a_tenth_to_a_planet_at_one_and_a_half():
    host = with_companion(star(), star(), 5.0)
    share = rf.companion_share(host, 1.5)
    assert share == pytest.approx(1.0 / (5.0 ** 2 - 1.5 ** 2) * 1.5 ** 2, rel=1e-6)
    assert 0.09 <= share <= 0.11, "the point-mass sum is 9 percent, the average over the planet's orbit 9.9"


def test_a_companion_at_a_hundred_au_moves_the_temperature_by_under_a_tenth_of_a_percent():
    host = with_companion(star(), star(), 100.0)
    ratio = 1.0 + rf.companion_share(host, 1.0)
    assert ratio ** 0.25 - 1.0 < 1.0e-3


def test_an_eccentric_companion_averages_to_more_light_than_a_circular_one_of_the_same_axis():
    circular = rf.companion_share(with_companion(star(), star(), 6.0, 0.0), 1.0)
    eccentric = rf.companion_share(with_companion(star(), star(), 6.0, 0.5), 1.0)
    assert eccentric > circular
    point_mass = 1.0 / (6.0 ** 2 * math.sqrt(1.0 - 0.5 ** 2)) * 1.0 ** 2
    assert eccentric == pytest.approx(point_mass, rel=0.10), "close to the 1 / (a^2 sqrt(1 - e^2)) of `planet_zones.flux_at` where the planet is far inside"


def test_the_average_matches_a_brute_force_average_over_both_orbits():
    axis, e, d = 7.0, 0.3, 2.0
    grid = 400
    mean = (np.arange(grid) + 0.5) * 2 * math.pi / grid
    anomaly = mean.copy()
    for _ in range(40):
        anomaly = anomaly - (anomaly - e * np.sin(anomaly) - mean) / (1 - e * np.cos(anomaly))
    radius = axis * (1 - e * np.cos(anomaly))
    true = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(anomaly / 2), np.sqrt(1 - e) * np.cos(anomaly / 2))
    phase = np.linspace(0, 2 * math.pi, 361)[:-1]
    total = 0.0
    for r, nu in zip(radius, true):
        sep2 = (r * np.cos(nu) - d * np.cos(phase)) ** 2 + (r * np.sin(nu) - d * np.sin(phase)) ** 2
        total += np.mean(1.0 / sep2)
    assert rf.mean_inverse_square(d, axis, e) == pytest.approx(total / grid, rel=2e-3)


def test_the_field_has_each_band_as_a_total_and_per_source():
    host = star(l_xuv_w=1.0e-3 * SUN_W, xuv_fluence_j=1.0e13, flare_n33_per_yr=0.5)
    companion = star(0.5 * SUN_W, l_xuv_w=1.0e-3 * 0.5 * SUN_W, xuv_fluence_j=2.0e13, flare_n33_per_yr=20.0)
    rf.attach(host, [rf.Companion(companion, 1, 4.0, 0.0)])
    field = rf.field_at(host, 1.0)
    assert set(field.bolometric_w_m2.by_source) == {None, 1}
    assert field.bolometric_w_m2.total == pytest.approx(sum(field.bolometric_w_m2.by_source.values()))
    assert field.bolometric_w_m2.by_source[None] == pytest.approx(old_light_w_m2(SUN_W, 1.0), rel=1e-12)
    assert field.flare_irradiation_index.by_source[None] == pytest.approx(0.5)
    assert field.flare_irradiation_index.by_source[1] == pytest.approx(20.0 / (16.0 - 1.0), rel=1e-3), "a flaring companion adds more flares than the quiet host"
    assert field.flare_irradiation_index.total > field.flare_irradiation_index.by_source[None]


def test_a_source_with_no_value_for_a_band_gives_none_of_it():
    host = star()
    companion = star(l_xuv_w=1.0e22)
    rf.attach(host, [rf.Companion(companion, 1, 5.0, 0.0)])
    field = rf.field_at(host, 1.0)
    assert set(field.xuv_flux_earth.by_source) == {1}
    assert field.exposure_index.total == 0.0


def test_a_body_with_no_luminosity_host_gets_no_light():
    assert rf.starlight_w_m2(star(None), 1.0) == 0.0 and rf.starlight_w_m2(star(), 0.0) == 0.0


def _tree():
    """Three stars: a pair (members 0, 1 at node 1) and a third (member 2) round the pair, at node 0."""
    nodes = [SimpleNamespace(node_no=0, parent_node_no=None, slot=0, kind="orbit", separation_km=60.0 * constants.AU_TO_KM, eccentricity=0.2),
             SimpleNamespace(node_no=1, parent_node_no=0, slot=0, kind="orbit", separation_km=0.5 * constants.AU_TO_KM, eccentricity=0.0)]
    members = [SimpleNamespace(member_no=0, parent_node_no=1, slot=0, mass_solar=1.0),
               SimpleNamespace(member_no=1, parent_node_no=1, slot=1, mass_solar=0.5),
               SimpleNamespace(member_no=2, parent_node_no=0, slot=1, mass_solar=0.8)]
    return nodes, members


def test_the_orbits_a_host_sees_are_those_of_the_nodes_it_meets_the_others_at():
    nodes, members = _tree()
    zone = planet_zones.Zone(planet_zones.S_TYPE, 0, 0.05, 0.1)
    found = {number: (axis, e) for number, axis, e in planet_zones.companion_orbits(zone, nodes, members)}
    assert found[1] == (pytest.approx(0.5), 0.0) and found[2] == (pytest.approx(60.0), 0.2)
    assert planet_zones.host_members(zone, nodes, members) == [0]
    pair = planet_zones.Zone(planet_zones.P_TYPE, 1, 1.0, 5.0)
    assert [number for number, *_ in planet_zones.companion_orbits(pair, nodes, members)] == [2]
    assert planet_zones.host_members(pair, nodes, members) == [0, 1]
