# tests/test_belt_populations.py

"""
GEN.1092: each stored member of an outer disk belongs to one of the scatter step's populations, and its orbit follows from it (`physics/belt_populations.py`, `belt_roster.draw_roster`).
"""

import math

import pytest

from planetgen.physics import belt_populations as pops
from planetgen.physics import belt_roster, constants, scatter
from planetgen.util import draw

STAR = constants.SOLAR_MASS_TO_KG
GIANT_AU = 30.0
LOWER, UPPER = 36.0, 150.0


def _shares(**kwargs):
    shares = dict.fromkeys(scatter.POPULATIONS, 0.0)
    shares.update(kwargs)
    return shares


def test_the_weights_are_the_in_disk_shares_and_a_belt_without_shares_is_all_cold():
    weights = dict(pops.member_weights(_shares(cold_classical=0.4, resonant=0.1, ejected=0.3, implanted_inner_belt=0.05, scattered_disk=0.15)))
    assert set(weights) == set(pops.IN_DISK) and weights["cold_classical"] == 0.4 and weights["resonant"] == 0.1
    assert pops.member_weights(None) == (("cold_classical", 1.0),)
    assert pops.member_weights(_shares(ejected=1.0)) == (("cold_classical", 1.0),)


def test_the_uniform_number_picks_by_the_weights():
    weights = (("a", 0.2), ("b", 0.3), ("c", 0.5))
    assert [pops.pick(weights, u) for u in (0.0, 0.19, 0.2, 0.49, 0.5, 0.999999)] == ["a", "a", "b", "b", "c", "c"]


def test_a_resonant_body_sits_at_the_giants_resonance_and_cold_falls_back_when_the_belt_misses_it():
    tag, a, scale, fixed = pops.orbit("resonant", 80.0, LOWER, UPPER, GIANT_AU, 0.0, 0.5)
    assert tag == "resonant" and a in {scatter.resonance_au(GIANT_AU, p, q) for p, q in scatter.RESONANCES} and scale is None
    assert pops.RESONANT_E_RANGE[0] <= fixed <= pops.RESONANT_E_RANGE[1]
    assert pops.orbit("resonant", 80.0, 60.0, 70.0, GIANT_AU, 0.0, 0.5)[0] == "cold_classical"


def test_scattered_crossing_and_detached_bodies_have_the_perihelion_their_population_gives():
    for u in (0.05, 0.5, 0.95):
        tag, a, scale, e = pops.orbit("scattered_disk", 90.0, LOWER, UPPER, GIANT_AU, 0.5, u)
        assert tag == "scattered_disk" and GIANT_AU <= a * (1.0 - e) <= LOWER + 1e-9
        tag, a, scale, e = pops.orbit("centaurs", 90.0, LOWER, UPPER, GIANT_AU, u, u)
        assert LOWER <= a <= 1.5 * LOWER + 1e-9 and a * (1.0 - e) < GIANT_AU        # crosses the giant's orbit
        tag, a, scale, e = pops.orbit("detached", 90.0, LOWER, UPPER, GIANT_AU, u, u)
        assert tag == "detached" and a * (1.0 - e) >= LOWER - 1e-9 and a <= UPPER and e <= pops.E_MAX_BOUND


def test_a_hot_body_is_stirred_like_the_belt_and_a_cold_one_less():
    assert pops.orbit("hot_classical", 80.0, LOWER, UPPER, GIANT_AU, 0.5, 0.5)[2] == 1.0
    assert pops.orbit("cold_classical", 80.0, LOWER, UPPER, GIANT_AU, 0.5, 0.5)[2] == pops.COLD_STIRRING < 1.0
    assert pops.orbit("scattered_disk", 80.0, LOWER, UPPER, None, 0.5, 0.5)[0] == "cold_classical"      # no giant: all cold


def test_an_unknown_population_is_refused():
    with pytest.raises(ValueError):
        pops.orbit("trojans", 80.0, LOWER, UPPER, GIANT_AU, 0.5, 0.5)


def _roster(populations, count=200, giant=GIANT_AU):
    return belt_roster.draw_roster(count, LOWER, UPPER, STAR, 2.7, 1.0, draw.Stream("pops"), populations, giant)


def test_the_roster_tags_the_bodies_in_the_shares_and_keeps_them_inside_the_belt():
    shares = _shares(cold_classical=0.4, hot_classical=0.1, resonant=0.2, scattered_disk=0.1, centaurs=0.1, detached=0.1)
    entries = _roster(pops.member_weights(shares))
    counts = {tag: sum(entry.tag == tag for entry in entries) for tag in pops.IN_DISK}
    assert all(counts[tag] > 0 for tag in pops.IN_DISK)
    assert counts["cold_classical"] > counts["hot_classical"]
    for entry in entries:
        assert LOWER <= entry.placement.distance_au <= UPPER and 0.0 <= entry.eccentricity <= pops.E_MAX_BOUND


def test_the_tags_do_not_move_the_draws_the_other_numbers_come_from():
    plain = _roster(None)
    tagged = _roster(pops.member_weights(_shares(cold_classical=1.0)))
    assert [entry.tag for entry in plain] == [None] * len(plain) and all(entry.tag == "cold_classical" for entry in tagged)
    assert [entry.node_deg for entry in plain] == [entry.node_deg for entry in tagged]
    assert [entry.placement.u_taxon for entry in plain] == [entry.placement.u_taxon for entry in tagged]
    assert _roster(pops.member_weights(_shares(cold_classical=0.5, resonant=0.5)), 50) == _roster(pops.member_weights(_shares(cold_classical=0.5, resonant=0.5)), 200)[:50]


def test_cold_members_are_stirred_less_than_the_belt():
    cold = _roster(pops.member_weights(_shares(cold_classical=1.0)), 300)
    hot = _roster(pops.member_weights(_shares(hot_classical=1.0)), 300)
    assert sum(entry.eccentricity for entry in hot) > 2.0 * sum(entry.eccentricity for entry in cold)


def _outer_disk_with_a_giant():
    from planetgen.generation.system import StarSystem
    from tests.test_belt_birth import _belts, _config
    config = _config()
    for seed in range(1, 60):
        draw.set_run_seed(seed)
        system = StarSystem(system_config=config)
        for belt in _belts(system):
            if belt.origin == "outer-disk" and belt.scatter_giant_au is not None and belt.members:
                return belt
    pytest.fail("no outer disk swept by a giant in 59 systems")


def test_a_generated_outer_disk_tags_its_members_and_the_tag_is_the_orbit_it_has():
    belt = _outer_disk_with_a_giant()
    tags = {member["scatter_tag"] for member in belt.members}
    assert tags <= set(pops.IN_DISK) and None not in tags
    for member in belt.members:
        assert belt.lower_limit <= member["distance_au"] <= belt.upper_limit
        if member["scatter_tag"] == "resonant":
            assert any(member["distance_au"] == pytest.approx(scatter.resonance_au(belt.scatter_giant_au, p, q)) for p, q in scatter.RESONANCES)
        if member["scatter_tag"] in ("scattered_disk", "centaurs", "detached"):
            assert member["eccentricity"] > 0.0


def test_a_belt_that_is_not_an_outer_disk_has_untagged_members():
    from planetgen.generation.system import StarSystem
    from tests.test_belt_birth import _belts, _config
    config = _config()
    for seed in range(1, 40):
        draw.set_run_seed(seed)
        for belt in _belts(StarSystem(system_config=config)):
            if belt.origin not in (None, "outer-disk") and belt.members:
                assert {member["scatter_tag"] for member in belt.members} == {None}
                return
    pytest.fail("no inner belt with members in 39 systems")
