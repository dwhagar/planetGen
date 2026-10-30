# tests/test_irregular_orbits.py

"""
GEN.1078: the orbits of irregular moons (`planetgen.physics.irregular_orbits`): the retrograde to prograde split from stability alone, the observed medians of a / r_H, no body in the Kozai band or inside the
outermost regular moon, and families.
"""

import math

import pytest

from planetgen.physics import irregular_orbits as orbits
from planetgen.physics import irregular_zone as zone
from planetgen.util import draw

SUN = 1.989e30
JUPITER = (5.203, 1.898e27, 1.8827e6 / 1.495979e8)
SATURN = (9.537, 5.683e26, 3.5613e6 / 1.495979e8)
URANUS = (19.19, 8.681e25, 5.8352e5 / 1.495979e8)
NEPTUNE = (30.07, 1.024e26, 0.0)
THE_FOUR = (JUPITER, SATURN, URANUS, NEPTUNE)


def _zone(planet):
    distance, mass, moon = planet
    return zone.zone(distance, mass, SUN, moon)


def _median(values):
    ordered = sorted(values)
    return ordered[len(ordered) // 2]


def test_the_median_before_survival_runs_through_the_four_giants_and_is_held_at_the_ends():
    for distance, median in orbits.MEDIAN_BEFORE_SURVIVAL:
        assert orbits.median_before_survival(distance) == pytest.approx(median)
    assert orbits.median_before_survival(1.0) == 1.1 and orbits.median_before_survival(100.0) == 0.193
    assert 0.172 < orbits.median_before_survival(25.0) < 0.332


def test_stability_alone_leaves_four_and_a_half_retrograde_for_each_prograde_at_a_jupiter():
    draw.set_run_seed(3)
    z = _zone(JUPITER)
    population = orbits.draw_irregulars(z, JUPITER[0], 40000)
    retrograde = sum(1 for o in population.orbits if o.retrograde)
    prograde = len(population.orbits) - retrograde
    assert retrograde / prograde == pytest.approx(4.5, rel=0.15)       # the known irregulars: 88 to 19
    assert population.lost["unstable"] > 0 and "kozai" not in population.lost


def test_the_survivors_have_the_observed_medians_of_a_over_the_hill_radius():
    draw.set_run_seed(4)
    for planet, observed in zip(THE_FOUR, (0.44, 0.29, 0.17, 0.19)):
        z = _zone(planet)
        population = orbits.draw_irregulars(z, planet[0], 20000)
        median = _median([o.semi_major_axis_au / z.hill_radius_au for o in population.orbits])
        assert median == pytest.approx(observed, rel=0.15)


def test_no_survivor_is_in_the_kozai_band_or_inside_the_outermost_regular_moon_or_beyond_its_limit():
    draw.set_run_seed(5)
    for planet in THE_FOUR:
        z = _zone(planet)
        for o in orbits.draw_irregulars(z, planet[0], 3000).orbits:
            assert not zone.in_kozai_band(o.inclination_deg)
            assert o.semi_major_axis_au * (1.0 - o.eccentricity) > z.inner_edge_au
            assert o.semi_major_axis_au <= zone.stable_limit_au(z, o.inclination_deg)
            assert 0.0 <= o.eccentricity <= orbits.ECCENTRICITY_MAX


def test_the_draw_is_deterministic_from_the_units_stream():
    z = _zone(JUPITER)
    draw.set_run_seed(8)
    first = orbits.draw_irregulars(z, JUPITER[0], 50)
    draw.set_run_seed(8)
    assert orbits.draw_irregulars(z, JUPITER[0], 50) == first


def test_the_gauss_changes_of_a_tangential_kick_on_a_circular_orbit_are_the_textbook_ones():
    parent = orbits.Orbit(0.1, 0.0, 30.0, 0.0, 0.0, 0.0)
    da, de, di = orbits.gauss_changes(parent, 1.898e27, 0.0, 100.0, 0.0)
    n = math.sqrt(6.674e-11 * 1.898e27 / (0.1 * 1.495979e11) ** 3)
    assert da / 0.1 == pytest.approx(2.0 * 100.0 / (n * 0.1 * 1.495979e11), rel=1e-3)
    assert de == pytest.approx(2.0 * 100.0 / (n * 0.1 * 1.495979e11), rel=1e-3)        # cos f + cos E = 2 at the pericentre of a circle
    assert di == 0.0


def test_a_family_clusters_round_its_parent_and_loses_the_fragments_that_no_longer_stay():
    draw.set_run_seed(9)
    z = _zone(JUPITER)
    parent = orbits.Orbit(0.15, 0.3, 150.0, 40.0, 100.0, 10.0)
    assert zone.verdict(z, parent.semi_major_axis_au, parent.eccentricity, parent.inclination_deg) is None
    family = orbits.family(z, parent, JUPITER[1], 200)
    assert len(family.orbits) + sum(family.lost.values()) == 200 and len(family.orbits) > 100
    spread = [abs(o.semi_major_axis_au / parent.semi_major_axis_au - 1.0) for o in family.orbits]
    assert max(spread) < 0.4 and sum(spread) / len(spread) < 0.1      # 300 m/s on a 2.4 km/s orbit moves the axis by up to a quarter
    assert all(abs(o.inclination_deg - parent.inclination_deg) < 10.0 for o in family.orbits)
