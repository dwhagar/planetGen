# tests/test_resonant_chains.py

"""
GEN.1043: resonant chains and their forced eccentricity (`planetgen.physics.resonant_chains`): the published Laplace-coefficient table, TRAPPIST-1, Kepler-223 and the
Galilean moons as chains, the Solar System as none, Io's and Europa's forced eccentricity, the fall of the eccentricity with the detuning, the width and crossing bounds,
and the hand-off to `tidal_heating`.
"""

import math

import pytest

from planetgen.physics import constants, tidal_heating
from planetgen.physics import resonant_chains as rc

EARTH_MASS = 5.972e24
SUN = constants.SOLAR_MASS_TO_KG
JUPITER = 1.8982e27

TRAPPIST_PERIODS = [1.510826, 2.421937, 4.049219, 6.101013, 9.207540, 12.352446, 18.772866]
TRAPPIST_MASSES = [m * EARTH_MASS for m in (1.374, 1.308, 0.388, 0.692, 1.039, 1.321, 0.326)]
TRAPPIST_STAR = 0.0898 * SUN
KEPLER_223 = [7.3845, 9.8456, 14.7887, 19.7257]
GALILEAN_PERIODS = [1.769138, 3.551181, 7.154553]
GALILEAN_MASSES = [8.9319e22, 4.7998e22, 1.4819e23]
SOLAR_SYSTEM = [87.97, 224.7, 365.25, 686.98, 4332.6, 10759.0, 30688.0, 60182.0]


@pytest.mark.parametrize("p, inner, outer", [(2, -1.1905, 0.4284), (3, -2.0252, 2.4840), (4, -2.8404, 3.2833)])
def test_the_laplace_coefficients_are_the_published_table_at_exact_resonance(p, inner, outer):
    alpha = ((p - 1) / p) ** (2.0 / 3.0)
    assert rc.inner_coefficient(p, alpha) == pytest.approx(inner, abs=2e-4)
    assert rc.outer_coefficient(p, alpha) == pytest.approx(outer, abs=2e-4)


def test_a_ratio_is_taken_to_the_nearest_target_inside_five_percent():
    assert rc.classify_pair(1.0, 1.5)[:2] == (3, 2) and rc.classify_pair(1.0, 1.5)[2] == pytest.approx(0.0)
    assert rc.classify_pair(1.0, 2.0)[:2] == (2, 1)
    assert rc.classify_pair(1.0, 4.0 / 3.0 * 1.04)[:2] == (4, 3)
    assert rc.classify_pair(1.0, 1.6)[:2] == (5, 3)                         # 8:5 is 4 percent under 5:3
    assert rc.classify_pair(1.0, 1.5 * 0.94) is None and rc.classify_pair(1.0, 2.5) is None
    assert rc.classify_pair(1.0, 1.0) is None


def test_trappist_1_is_one_chain_of_seven_planets():
    chains = rc.find_chains(TRAPPIST_PERIODS)
    assert len(chains) == 1 and chains[0].members == tuple(range(7))
    assert [(p.p, p.q) for p in chains[0].pairs] == [(5, 3), (5, 3), (3, 2), (3, 2), (4, 3), (3, 2)]
    assert all(abs(p.delta) < 0.04 for p in chains[0].pairs)


def test_the_galilean_moons_and_kepler_223_are_chains_and_the_solar_system_is_none():
    laplace = rc.find_chains(GALILEAN_PERIODS)
    assert len(laplace) == 1 and laplace[0].members == (0, 1, 2) and all((p.p, p.q) == (2, 1) for p in laplace[0].pairs)
    assert [p.delta for p in laplace[0].pairs] == pytest.approx([0.0037, 0.0074], abs=0.0005)
    kepler = rc.find_chains(KEPLER_223)
    assert len(kepler) == 1 and kepler[0].members == (0, 1, 2, 3) and [(p.p, p.q) for p in kepler[0].pairs] == [(4, 3), (3, 2), (4, 3)]
    assert rc.find_chains(SOLAR_SYSTEM) == ()


def test_two_resonant_planets_are_a_pair_and_not_a_chain_and_a_gap_splits_a_chain():
    assert len(rc.find_pairs([1.0, 1.5])) == 1 and rc.find_chains([1.0, 1.5]) == ()
    two = rc.find_chains([1.0, 1.5, 2.25, 5.0, 7.5, 11.25])
    assert [c.members for c in two] == [(0, 1, 2), (3, 4, 5)]


def test_the_chain_indexes_the_input_whatever_its_order():
    shuffled = [GALILEAN_PERIODS[2], GALILEAN_PERIODS[0], GALILEAN_PERIODS[1]]
    chain = rc.find_chains(shuffled)[0]
    assert chain.members == (1, 2, 0)


def test_io_is_forced_by_europa_to_about_its_observed_eccentricity_and_europa_by_io_and_ganymede():
    e = rc.forced_eccentricities(GALILEAN_PERIODS, GALILEAN_MASSES, JUPITER)
    # observed forced eccentricities 0.0041, 0.0094, 0.0013 (Yoder and Peale 1981); Io's precession, not in the formula, cuts its detuning's effect
    assert e[0] == pytest.approx(0.0041, rel=0.35)
    assert e[1] == pytest.approx(0.0094, rel=0.25)
    assert e[2] < e[0] < e[1] and 0.0003 < e[2] < 0.003


def test_trappist_1s_forced_eccentricities_are_the_design_notes_thousandths():
    e = rc.forced_eccentricities(TRAPPIST_PERIODS, TRAPPIST_MASSES, TRAPPIST_STAR)
    assert e[0] == 0.0 and e[1] == 0.0                         # the 5:3 pairs are second order and force nothing at this order
    assert all(0.001 < v < 0.01 for v in e[2:])


def test_a_planet_in_no_chain_has_a_circular_orbit():
    e = rc.forced_eccentricities([1.0, 1.5, 2.25, 40.0], [EARTH_MASS] * 4, SUN)
    assert e[3] == 0.0 and all(v > 0.0 for v in e[:3])


def test_the_eccentricity_falls_as_the_chain_moves_off_exact_resonance():
    def first(delta):
        periods = [1.0, 1.5 * (1.0 + delta), 2.25 * (1.0 + delta) ** 2]
        return rc.forced_eccentricities(periods, [3.0 * EARTH_MASS] * 3, 0.5 * SUN)[0]

    values = [first(d) for d in (0.01, 0.02, 0.04)]
    assert values == sorted(values, reverse=True) and values[0] == pytest.approx(2.0 * values[1], rel=0.1)


def test_exact_resonance_does_not_blow_up_it_is_held_at_the_width_of_the_resonance():
    exact = rc.forced_eccentricities([1.0, 1.5, 2.25], [3.0 * EARTH_MASS] * 3, 0.5 * SUN)
    mu = 3.0 * EARTH_MASS / (0.5 * SUN)
    assert all(math.isfinite(v) and v < 0.1 for v in exact)
    near = rc.forced_eccentricities([1.0, 1.5 * (1.0 + mu ** (2.0 / 3.0)), 2.25 * (1.0 + mu ** (2.0 / 3.0)) ** 2], [3.0 * EARTH_MASS] * 3, 0.5 * SUN)
    assert exact[0] == pytest.approx(near[0], rel=0.15)


def test_heavy_planets_are_held_to_orbits_that_do_not_cross():
    masses = [1.0e-3 * SUN, 1.0e-3 * SUN, 1.0e-3 * SUN]
    periods = [1.0, 1.5, 2.25]
    e = rc.forced_eccentricities(periods, masses, SUN)
    a = [rc.semi_major_axis_m(p, SUN, m) for p, m in zip(periods, masses)]
    assert all(a[i] * (1.0 + e[i]) < a[i + 1] * (1.0 - e[i + 1]) for i in range(2))
    limit = rc.CROSSING_SHARE * (a[1] - a[0]) / (a[1] + a[0])
    assert e[0] <= limit + 1e-12


def test_the_forced_orbits_of_trappist_1_do_not_cross_and_keep_the_planets_inside_a_hill_separation():
    e = rc.forced_eccentricities(TRAPPIST_PERIODS, TRAPPIST_MASSES, TRAPPIST_STAR)
    a = [rc.semi_major_axis_m(p, TRAPPIST_STAR, m) for p, m in zip(TRAPPIST_PERIODS, TRAPPIST_MASSES)]
    assert all(a[i] * (1.0 + e[i]) < a[i + 1] * (1.0 - e[i + 1]) for i in range(6))


def test_the_forced_eccentricity_feeds_the_tide_as_it_does_for_a_moon():
    # Io's observed forced eccentricity through the melt-limited power gives the 1e14 W the model was fitted to
    io = rc.forcing(JUPITER, GALILEAN_PERIODS[0], GALILEAN_MASSES[0], 0.0041)
    assert tidal_heating.molten_power_w(io, 1.8216e6) == pytest.approx(1.0e14, rel=0.6)
    # TRAPPIST-1 d: a solid rocky planet, tide from the forced eccentricity, a modest fraction of Io's flux
    e = rc.forced_eccentricities(TRAPPIST_PERIODS, TRAPPIST_MASSES, TRAPPIST_STAR)
    d = rc.forcing(TRAPPIST_STAR, TRAPPIST_PERIODS[2], TRAPPIST_MASSES[2], e[2])
    power = tidal_heating.solid_power_w(d, TRAPPIST_MASSES[2], 0.788 * 6.371e6)
    flux = power / (4.0 * math.pi * (0.788 * 6.371e6) ** 2)
    assert 0.0 < flux < 1.0
    circular = rc.forcing(TRAPPIST_STAR, TRAPPIST_PERIODS[2], TRAPPIST_MASSES[2], 0.0)
    assert tidal_heating.solid_power_w(circular, TRAPPIST_MASSES[2], 0.788 * 6.371e6) == 0.0


def test_inputs_the_model_cannot_use_are_refused():
    with pytest.raises(rc.ResonantChainError):
        rc.find_pairs([1.0, -2.0])
    with pytest.raises(rc.ResonantChainError):
        rc.forced_eccentricities([1.0, 1.5, 2.25], [EARTH_MASS] * 2, SUN)
    with pytest.raises(rc.ResonantChainError):
        rc.forced_eccentricities([1.0, 1.5, 2.25], [EARTH_MASS, 0.0, EARTH_MASS], SUN)
    with pytest.raises(rc.ResonantChainError):
        rc.forced_eccentricities([1.0, 1.5, 2.25], [EARTH_MASS] * 3, 0.0)
