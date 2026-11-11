# tests/test_trojan_neighbours.py

"""GEN.1455: another giant or a companion star takes part of a giant's trojan swarm (`trojans.neighbour_stable_share`, `scatter.scatter_shares(trojan_survival=...)`, `belt_birth._fates`)."""

import math
import types

import pytest

from planetgen.physics import constants, scatter, trojans

SUN = 1.0
JUPITER = 9.546e-4
SATURN = 2.858e-4


def test_a_lone_giant_keeps_its_whole_swarm():
    assert trojans.surviving_share(SUN, JUPITER, 5.2, []) == 1.0
    assert trojans.neighbour_stable_share(SUN, JUPITER, 5.2, 0.0, 9.5) == 1.0


def test_a_neighbour_inside_the_hill_stable_separation_takes_the_swarm():
    hill = ((JUPITER + JUPITER) / 3.0) ** (1.0 / 3.0) * 0.5 * (5.2 + 5.8)
    assert trojans.neighbour_stable_share(SUN, JUPITER, 5.2, JUPITER, 5.8) == 0.0
    near = 5.2 + trojans.HILL_SEPARATION_K * hill * 0.99
    assert trojans.neighbour_stable_share(SUN, JUPITER, 5.2, JUPITER, near) < 0.15


def test_a_second_giant_at_twice_the_distance_cuts_the_swarm_by_the_expected_amount():
    share = trojans.neighbour_stable_share(SUN, JUPITER, 5.2, JUPITER, 10.4)
    hill = ((2 * JUPITER) / 3.0) ** (1.0 / 3.0) * 0.5 * (5.2 + 10.4)
    x = 5.2 / (trojans.HILL_SEPARATION_K * hill)
    assert share == pytest.approx(1.0 - math.exp(-trojans.STABLE_RAMP * (x - 1.0)))
    assert 0.5 < share < 0.95, "a loss of tens of percent, as the audit guessed"


def test_the_share_grows_with_the_separation_and_a_distant_neighbour_leaves_it_whole():
    shares = [trojans.neighbour_stable_share(SUN, JUPITER, 5.2, SATURN, au) for au in (6.0, 7.0, 9.5, 20.0, 200.0)]
    assert shares == sorted(shares) and shares[-1] > 0.99


def test_every_neighbour_takes_its_part():
    one = trojans.neighbour_stable_share(SUN, JUPITER, 5.2, SATURN, 9.5)
    two = trojans.neighbour_stable_share(SUN, JUPITER, 5.2, SATURN, 2.0)
    assert trojans.surviving_share(SUN, JUPITER, 5.2, [(SATURN, 9.5), (SATURN, 2.0)]) == pytest.approx(one * two)


def test_the_trojans_the_swarm_loses_are_ejected_and_the_shares_still_sum_to_one():
    args = (1.0, 317.8, 30.0, 11.2, 5.0, 15.0, 50.0)
    whole = scatter.scatter_shares(*args)
    half = scatter.scatter_shares(*args, trojan_survival=0.5)
    assert half["trojans"] == pytest.approx(0.5 * whole["trojans"])
    assert half["ejected"] == pytest.approx(whole["ejected"] + 0.5 * whole["trojans"])
    assert sum(half.values()) == pytest.approx(1.0) and sum(whole.values()) == pytest.approx(1.0)
    assert scatter.scatter_shares(*args, trojan_survival=1.0) == whole, "a lone giant gives the shares of before to the last digit"


def test_a_companion_star_limit_takes_the_swarm_of_a_giant_near_it():
    assert trojans.star_limit_stable_share(30.0, None) == 1.0
    assert trojans.star_limit_stable_share(30.0, 50.0) == 0.0, "a giant at 30 AU under a limit of 50 AU is past half of it"
    share = trojans.star_limit_stable_share(30.0, 150.0)
    assert share == pytest.approx(1.0 - math.exp(-trojans.STABLE_RAMP * (150.0 / (trojans.STAR_LIMIT_MARGIN * 30.0) - 1.0))) and 0.5 < share < 1.0


def test_the_disk_builder_hands_the_survival_to_the_scatter():
    from planetgen.generation import belt_birth
    from planetgen.physics import belt_origin

    giant = types.SimpleNamespace(body_type='g', distance=30.0, mass=317.8 * constants.EARTH_MASS_TO_KG, radius=11.2 * 6371.0, formation_distance_au=25.0)
    inner = types.SimpleNamespace(body_type='g', distance=12.0, mass=100.0 * constants.EARTH_MASS_TO_KG, radius=9.0 * 6371.0, formation_distance_au=12.0)
    origin = belt_origin.BeltOrigin("outer-disk", (), 15.0, 50.0, 2.7, 1.0, 1.0, 0.0, 0.0, False, (), "")
    sun = constants.SOLAR_MASS_TO_KG
    lone = belt_birth._fates(types.SimpleNamespace(a_crit_au=None), sun, [giant], origin, 50.0, 0.0)[0]
    assert lone["trojans"] > 0.0
    with_giant = belt_birth._fates(types.SimpleNamespace(a_crit_au=None), sun, [inner, giant], origin, 50.0, 0.0)[0]
    assert with_giant["trojans"] == pytest.approx(lone["trojans"] * trojans.neighbour_stable_share(sun, giant.mass, 30.0, inner.mass, 12.0))
    pair = belt_birth._fates(types.SimpleNamespace(a_crit_au=120.0), sun, [giant], origin, 50.0, 0.0)[0]
    assert pair["trojans"] == pytest.approx(lone["trojans"] * trojans.star_limit_stable_share(30.0, 120.0))
