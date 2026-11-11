# tests/test_moon_resonances.py

"""GEN.1449 part 2: the capture of moon pairs in a first-order resonance on the way out and the eccentricity it forces (`physics/moon_resonances`), against Jupiter's Galilean moons."""

import math

import pytest

from planetgen.physics import moon_resonances, resonant_chains

JUPITER_KG = 1.898e27
GALILEAN_KG = [8.93e22, 4.80e22, 1.482e23]          # Io, Europa, Ganymede
GALILEAN_DAYS = [1.769, 3.551, 7.155]


def _ratio(distances, i):
    return (distances[i + 1] / distances[i]) ** 1.5


def test_a_pair_that_swept_two_to_one_is_held_there_and_keeps_its_angular_momentum():
    masses, born, ended = [8.9e22, 4.8e22], [3.0e8, 9.0e8], [5.0e8, 5.3e8]
    held, made = moon_resonances.capture(masses, born, ended)
    assert [(c.p, c.q) for c in made] == [(2, 1)]
    assert _ratio(held, 0) == pytest.approx(2.0 * (1.0 + moon_resonances.CAPTURE_OFFSET), rel=1e-12)
    assert sum(m * math.sqrt(a) for m, a in zip(masses, held)) == pytest.approx(sum(m * math.sqrt(a) for m, a in zip(masses, ended)), rel=1e-12)
    assert held[0] >= born[0] and held[1] > ended[1]                  # the inner moon is held back a little, the outer pushed out


def test_the_first_resonance_met_on_the_way_in_is_the_one_that_captures():
    """A pair born at a ratio of 1.45 meets 4:3 and not 3:2 (the ratio never was above 3:2); one born at 3 meets 2:1 first."""
    assert moon_resonances.swept_resonance(1.45, 1.1) == (4, 3)
    assert moon_resonances.swept_resonance(3.0, 1.1) == (2, 1)
    assert moon_resonances.swept_resonance(1.2, 1.1) is None             # born inside 4:3
    assert moon_resonances.swept_resonance(1.45, 1.4) is None            # never reached it
    assert moon_resonances.swept_resonance(1.1, 3.0) is None             # a pair that moved apart is not captured


def test_a_chain_of_three_moons_is_held_in_both_resonances_with_the_momentum_conserved():
    masses, born, ended = GALILEAN_KG, [2.0e8, 6.0e8, 1.2e9], [4.0e8, 4.6e8, 5.5e8]
    held, made = moon_resonances.capture(masses, born, ended)
    assert [(c.p, c.q) for c in made] == [(2, 1), (2, 1)]
    assert _ratio(held, 0) == pytest.approx(2.004, rel=1e-9) and _ratio(held, 1) == pytest.approx(2.004, rel=1e-9)
    assert sum(m * math.sqrt(a) for m, a in zip(masses, held)) == pytest.approx(sum(m * math.sqrt(a) for m, a in zip(masses, ended)), rel=1e-12)


def test_a_heavy_outer_moon_that_would_hold_the_inner_one_inside_its_birth_orbit_is_not_captured():
    masses, born, ended = [1.0e20, 1.0e24], [3.0e8, 9.0e8], [3.01e8, 3.1e8]
    held, made = moon_resonances.capture(masses, born, ended)
    assert made == () and held == ended


def test_no_swept_pair_leaves_the_distances_as_they_were():
    held, made = moon_resonances.capture([1.0e22, 2.0e22, 3.0e22], [1.0e8, 2.0e8, 4.0e8], [1.0e8, 2.1e8, 4.2e8])
    assert made == () and held == [1.0e8, 2.1e8, 4.2e8]


def test_the_galilean_moons_are_forced_to_the_observed_order_of_eccentricity():
    forced = resonant_chains.pair_forced_eccentricities(GALILEAN_DAYS, GALILEAN_KG, JUPITER_KG)
    assert forced[0] == pytest.approx(0.0041, rel=0.4)           # Io (observed forced eccentricity 0.0041)
    assert forced[1] == pytest.approx(0.0094, rel=0.4)           # Europa (0.0094)
    assert 0.0 < forced[2] < forced[1]                            # Ganymede, the least forced
    chain = resonant_chains.forced_eccentricities(GALILEAN_DAYS, GALILEAN_KG, JUPITER_KG)
    assert chain == pytest.approx(forced, rel=1e-12)             # a chain of three gives the same by the chain and the pair route


def test_two_moons_are_enough_for_a_pair_and_a_chain_needs_three():
    days, kg = GALILEAN_DAYS[:2], GALILEAN_KG[:2]
    assert resonant_chains.forced_eccentricities(days, kg, JUPITER_KG) == (0.0, 0.0)
    assert all(e > 0.0 for e in resonant_chains.pair_forced_eccentricities(days, kg, JUPITER_KG))


def test_second_order_pairs_and_distant_pairs_force_nothing():
    assert resonant_chains.pair_forced_eccentricities([3.0, 5.0], [1.0e22, 1.0e22], JUPITER_KG) == (0.0, 0.0)      # 5:3
    assert resonant_chains.pair_forced_eccentricities([3.0, 8.0], [1.0e22, 1.0e22], JUPITER_KG) == (0.0, 0.0)
