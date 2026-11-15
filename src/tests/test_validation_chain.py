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


# The window of M3V seed 1317 as it stood after the settle pass: gaps of 10.7 and 11.2 mutual Hill radii and a 276 Earth-mass giant last, so the spacing a planet beyond can reach saturates at 12.0 (a last gap
# tops out near 22 mutual Hill radii) and only a distance of 18.5 AU holds it. A tenth of a percent of eccentricity either way turned this between "no distance closes it" and 18.5 AU.
SEED_1317_WINDOW = [(0.17736489065683342, 0.00010530538382627255, 2.1589670278845407e-05 * MSUN_KG), (0.2663213215446025, 0.03867301011741685, 3.1198164823724346e-05 * MSUN_KG),
                    (0.780980104277751, 0.011812512097469878, 0.000828067362132774 * MSUN_KG)]


def test_a_window_only_a_far_distance_closes_is_left_to_the_pairwise_rule():
    start = 1.9469395113748464
    found = ledger.chain_distance_au(SEED_1317_WINDOW, 3.0742534385850345e-06 * MSUN_KG, 0.010332050185493899, 0.4787564990065853 * MSUN_KG, start)
    assert found == start, "a push of more than CHAIN_MAX_PUSH times the pairwise distance is not the chain rule"


def test_the_push_is_the_nearest_distance_when_a_near_one_closes_the_window(monkeypatch):
    planets = system(SEED_282)
    body = planets[-1]
    inner = [(p.distance, p.eccentricity, p.mass) for p in planets[:-1]]
    start = body.distance * 0.9
    found = ledger.chain_distance_au(inner, body.mass, body.eccentricity, body.star.mass, start)
    monkeypatch.setattr(ledger, "CHAIN_SCAN_STEP", 1.001)
    assert ledger.chain_distance_au(inner, body.mass, body.eccentricity, body.star.mass, start) == pytest.approx(found, rel=1e-12)


@pytest.mark.parametrize("seed", [1273, 1317])
def test_m3v_systems_that_failed_the_final_check_after_the_settle_pass_validate(seed):
    """Seed 1273: an eccentricity the oddities forced after the spacing pass closed a belt's gap by 0.13 percent; seed 1317: the knife edge above."""
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from planetgen.util import draw

    config = SystemConfig()
    config.STAR_TYPE, config.PLANETS, config.COMETS, config.BINARY_SYSTEM = "M3V", True, False, False
    with draw.bound(seed):
        system = StarSystem(system_config=config)
    assert validation.check_star_system(system) == []


def test_a_push_that_would_open_a_gap_the_belt_rule_reads_as_a_belt_is_left_to_the_pairwise_rule(monkeypatch):
    """GEN.1521, CH5: the effective spacing (mean less half the spread) rewards one wide gap among narrow ones, so a window of two gaps of 10 mutual Hill radii pushed the next planet to a gap of 30.5 and
    windows beside a wide gap pushed it to 43 to 60, which the belt rule (BG1) reads as a belt. A push whose first distance opens a gap of `CHAIN_MAX_GAP_K` or more leaves the pairwise rule."""
    from planetgen.physics import belt_origin
    assert ledger.CHAIN_MAX_GAP_K == belt_origin.GAP_SPACING_K
    sun = 1.989e30
    mass = 3e-6 * sun
    share = (2 * mass / (3 * sun)) ** (1 / 3)

    def beyond(inner_au, gap):
        return inner_au * (1 + 0.5 * gap * share) / (1 - 0.5 * gap * share)

    first = 1.0
    second = beyond(first, 10.0)
    third = beyond(second, 10.0)
    window = [(first, 0.0, mass), (second, 0.0, mass), (third, 0.0, mass)]
    start = beyond(third, 10.0)
    found = ledger.chain_distance_au(window, mass, 0.0, sun, start)
    opened = ledger.edge_gap_hill_radii(third, 0.0, mass, found, 0.0, mass, sun)
    assert found > start and opened == pytest.approx(30.5, abs=0.1)
    monkeypatch.setattr(ledger, "CHAIN_MAX_GAP_K", 25.0)
    assert ledger.chain_distance_au(window, mass, 0.0, sun, start) == start
