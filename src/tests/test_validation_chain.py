"""GEN.1521: a chain of four or more planets is held to an effective spacing, not only the pairwise floor."""

import types

import pytest

from planetgen import tuning
from planetgen.generation import validation
from planetgen.physics import clearance_ledger as ledger

MSUN_KG = 1.98847e30

# Two systems of Lane 04's N-body spot-check (GEN.1468) that scattered within 2e4 years though every pair passed the screen: (a AU, e, m Msun), host Msun
SEED_282 = (0.44670248057178824, [(0.042543, 0.112, 7.589939893381247e-07), (0.054283, 0.0001, 1.8903618049226502e-06),
                                  (0.062347, 0.0007, 1.5886330444662975e-06), (0.073356, 0.0027, 3.7699195416536864e-06)])
SEED_264 = (0.1537945953236674, [(0.036569, 0.0001, 3.840722160054173e-07), (0.043187, 0.0366, 5.436476860842643e-07), (0.053526, 0.0037, 1.1605208716864998e-06),
                                 (0.068737, 0.0795, 5.966364074775377e-07), (0.088321, 0.0013, 2.01842853536718e-06), (0.125652, 0.0711, 3.3155636182705254e-06)])


def system(system_of):
    host, rows = system_of
    star = types.SimpleNamespace(mass=host * MSUN_KG)
    return [types.SimpleNamespace(body_type="t", distance=a, eccentricity=e, mass=m * MSUN_KG, star=star) for a, e, m in rows]


def gaps_of(planets):
    central = planets[0].star.mass
    return [ledger.edge_gap_hill_radii(a.distance, a.eccentricity, a.mass, b.distance, b.eccentricity, b.mass, central) for a, b in zip(planets, planets[1:])]


def test_the_effective_spacing_is_the_mean_less_half_the_spread():
    assert ledger.effective_spacing([12.0, 12.0, 12.0]) == pytest.approx(12.0)
    assert ledger.effective_spacing([10.0, 12.0, 14.0]) == pytest.approx(12.0 - 0.5 * (8.0 / 3.0) ** 0.5)


@pytest.mark.parametrize("data", [SEED_282, SEED_264], ids=["seed 282", "seed 264"])
def test_the_chains_that_scattered_in_the_nbody_check_are_refused_though_every_pair_passes_the_floor(data):
    planets = system(data)
    assert min(gaps_of(planets)) >= tuning.MUTUAL_HILL_RADII_SEPARATION - 0.2            # the pairwise screen is met
    problems = validation.check_orbit_spacing(planets)
    assert problems, "a chain at the floor is a problem"


def test_the_nearest_distance_for_a_chain_holds_the_effective_spacing_and_no_nearer_one_does():
    planets = system(SEED_282)
    body = planets[-1]
    inner = [(p.distance, p.eccentricity, p.mass) for p in planets[:-1]]
    start = body.distance * 0.9
    found = ledger.chain_distance_au(inner, body.mass, body.eccentricity, body.star.mass, start)
    assert found > start

    def k_eff(distance):
        gaps = gaps_of(planets[:-1]) + [ledger.edge_gap_hill_radii(*inner[-1], distance, body.eccentricity, body.mass, body.star.mass)]
        return ledger.effective_spacing(gaps)

    assert k_eff(found) == pytest.approx(tuning.CHAIN_K_EFF, abs=1e-6)
    assert k_eff(found * 0.999) < tuning.CHAIN_K_EFF


def test_fewer_than_four_planets_a_belt_in_the_chain_or_the_test_switched_off_leave_the_pairwise_rule(monkeypatch):
    planets = system(SEED_282)
    body = planets[-1]
    inner = [(p.distance, p.eccentricity, p.mass) for p in planets[:-1]]
    assert ledger.chain_distance_au(inner[1:], body.mass, body.eccentricity, body.star.mass, 0.05) == 0.05      # two inner planets: a chain of three
    belt = types.SimpleNamespace(body_type="a", lower_limit=0.05, upper_limit=0.06, distance=0.055, eccentricity=0.0, mass=1e20, star=body.star)
    with_belt = planets[:2] + [belt] + planets[2:3]
    pairwise = validation.min_distance_after_au(body, with_belt)
    assert pairwise == validation._chain_needed(body, with_belt, pairwise)
    monkeypatch.setattr(tuning, "CHAIN_K_EFF", 0.0)
    assert ledger.chain_distance_au(inner, body.mass, body.eccentricity, body.star.mass, 0.05) == 0.05


def test_a_chain_spaced_as_the_rule_asks_validates_and_a_loose_one_is_untouched():
    planets = system(SEED_282)
    inner = planets[:-1]
    body = planets[-1]
    body.distance = validation.min_distance_after_au(body, inner)
    assert validation.check_orbit_spacing(inner + [body]) == []
    loose = [types.SimpleNamespace(body_type="t", distance=0.1 * 2.5 ** i, eccentricity=0.0, mass=1e24, star=body.star) for i in range(5)]
    assert validation.check_orbit_spacing(loose) == []
