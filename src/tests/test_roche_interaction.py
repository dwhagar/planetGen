# tests/test_roche_interaction.py

"""
GEN.269: Roche geometry (`planetgen/physics/roche.py`) and the interaction
classifier of the member tree (`planetgen/generation/multistar.py`), with the
audit of docs/design/systems/multi-star-systems.md section 15 (13, interaction
consistency).
"""

import math

import pytest

from planetgen.generation import multistar, star_labels
from planetgen.physics import roche
from planetgen.util import draw


def test_the_roche_lobe_matches_eggleton_for_known_ratios():
    assert roche.roche_lobe_fraction(1.0, 1.0) == pytest.approx(0.3789, abs=1e-3)
    assert roche.roche_lobe_fraction(1.0, 10.0) == pytest.approx(0.2068, abs=1e-3)
    assert roche.roche_lobe_fraction(10.0, 1.0) == pytest.approx(0.578, abs=2e-3)
    assert roche.roche_lobe_radius_au(2.0, 1.0, 1.0) == pytest.approx(2.0 * 0.3789, abs=2e-3)


def test_the_filling_factor_and_its_inverse_agree():
    separation = roche.separation_for_filling(0.01, 1.0, 1.0, 0.5)
    assert roche.filling_factor(0.01, separation, 1.0, 0.5) == pytest.approx(1.0)
    assert roche.filling_factor(0.01, 2 * separation, 1.0, 0.5) == pytest.approx(0.5)


def _leaf(kind, mass, radius_au):
    return multistar._Leaf(kind, None, mass, mass, radius_au)


def _classify(a, b, axis, e=0.0):
    builder = multistar._Builder([a, b], 5e9, draw.Stream(1))
    return builder.interaction(0, 1, axis, e)


def test_the_classifier_names_each_kind_of_pair_from_its_radii():
    star = _leaf("main-sequence", 1.0, 0.0046)
    near = roche.separation_for_filling(0.0046, 1.0, 1.0, 1.0)
    assert _classify(star, star, near * 3)[0] == "detached"
    assert _classify(star, star, near)[0] == "contact"
    giant = _leaf("giant", 1.0, 0.05)
    assert _classify(giant, star, roche.separation_for_filling(0.05, 1.0, 1.0, 1.0))[0] == "semi-detached"
    assert _classify(giant, star, roche.separation_for_filling(0.05, 1.0, 1.0, 1.0))[1] == {1: multistar.MEMBER_MASS_GAINER}
    hole = _leaf("black-hole", 9.0, 1e-7)
    assert _classify(giant, hole, roche.separation_for_filling(0.05, 1.0, 1.0, 9.0))[0] == "accreting"
    assert _classify(star, hole, 0.2)[0] == "detached"
    # A 0.6 solar mass white dwarf has a past: a 2 solar mass star that swelled to 1.5 AU and lost 1.4 solar masses (GEN.1012 reads the
    # history, so a dwarf with no larger past than its present is never forced).
    dwarf = multistar._Leaf("white-dwarf", None, 0.6, 2.0, 4e-5, radius_max_au=1.5, mass_max_solar=2.0)
    assert _classify(dwarf, star, 0.05)[0] == "post-common-envelope"
    assert _classify(dwarf, star, 20.0)[0] == "detached"
    # The same dwarf next to a black hole at 0.01 AU had its giant swallow that orbit; a dwarf with no larger past is only a close, dark pair.
    assert _classify(dwarf, hole, 0.01)[0] == "post-common-envelope"
    assert _classify(_leaf("white-dwarf", 0.6, 4e-5), hole, 0.01)[0] == "detached"


def test_eccentricity_brings_the_closest_approach_into_the_lobe():
    star = _leaf("main-sequence", 1.0, 0.0046)
    near = roche.separation_for_filling(0.0046, 1.0, 1.0, 1.0)
    assert _classify(star, star, near * 1.5, 0.0)[0] == "detached"
    assert _classify(star, star, near * 1.5, 0.5)[0] == "contact"


def test_every_pair_the_generator_draws_is_consistent_with_its_radii():
    seen = {}
    for seed in range(1500):
        for placed, age in ((multistar.Placed("star", 1.0, 1.0), 4.6e9), (multistar.Placed("star", 20.0), 1e6),
                            (multistar.Placed("white-dwarf", 0.6, 1.5), 5e9),
                            (multistar.Placed("neutron-star", 1.4, 10.0, recycled=True), 3e9)):
            tree = multistar.draw_member_tree(placed, age, draw.Stream(seed))
            for node in tree.nodes:
                assert node.interaction in multistar.INTERACTIONS
                members = [m for m in tree.members if m.parent_node_no == node.node_no]
                seen[node.interaction] = seen.get(node.interaction, 0) + 1
                if len(members) != 2:
                    assert node.interaction == "detached"
                    continue
                a, b = sorted(members, key=lambda m: m.slot)
                peri = node.separation_au * (1 - node.eccentricity)
                fills = [roche.filling_factor(x.radius_au, peri, x.mass_solar, y.mass_solar)
                         for x, y in ((a, b), (b, a))]
                assert max(fills) <= multistar.MAX_FILL_FACTOR + 1e-9
                if node.interaction == "semi-detached":
                    assert max(fills) >= 1.0 and min(fills) < 1.0
                    assert all(m.kind not in multistar.ACCRETOR_KINDS for m in members)
                if node.interaction == "contact":
                    assert min(fills) >= 1.0
                if node.interaction == "accreting":
                    assert sum(m.kind in multistar.ACCRETOR_KINDS for m in members) == 1
                    assert max(fills) >= 1.0
                if node.interaction == "post-common-envelope":
                    # GEN.1012: any member that shrank (a remnant, or a young star that contracted) can force it, not only a white dwarf.
                    assert any(m.kind in ("white-dwarf", "neutron-star", "black-hole", "pre-main-sequence") for m in members)
    assert {"detached", "post-common-envelope"} <= set(seen)


def test_an_ultra_wide_system_is_flagged_fragile_and_only_that_one():
    placed = multistar.Placed("star", 1.0, 1.0)
    wide = narrow = 0
    for seed in range(3000):
        tree = multistar.draw_member_tree(placed, 4.6e9, draw.Stream(seed), n_stars=2)
        flagged = [bool(m.flags & multistar.MEMBER_FRAGILE) for m in tree.members]
        if tree.nodes[0].separation_au > multistar.FRAGILE_AXIS_AU if tree.nodes else False:
            assert all(flagged)
            wide += 1
        else:
            assert not any(flagged)
            narrow += 1
    assert wide > 0 and narrow > 0


def test_a_merger_product_has_a_rejuvenated_age_and_the_system_age_stays():
    tree = multistar.draw_member_tree(multistar.Placed("star", 1.0, 1.0), 4.6e9, draw.Stream(2), n_stars=2)
    low, high = multistar.REJUVENATION_FRACTION_RANGE
    for seed in range(200):
        assert low * 4.6e9 <= multistar.rejuvenated_age_yr(tree.age_yr, draw.Stream(seed)) <= high * 4.6e9
    assert tree.age_yr == 4.6e9 and tree.members[0].effective_age_yr is None


def test_the_new_member_kinds_exist_and_have_labels():
    for kind in ("stripped-helium", "thorne-zytkow", "wolf-rayet"):
        assert kind in multistar.MEMBER_KINDS
        assert star_labels.star_label("O5", kind)[0]
    assert star_labels.star_label("O5", "wolf-rayet") == ("Wolf-Rayet star", "Wolf-Rayet stars")
    assert star_labels.star_label("G2", "V") == ("G-type", "G-type")


def test_the_flag_bits_are_distinct_powers_of_two():
    values = list(multistar.MEMBER_FLAGS.values())
    assert len(set(values)) == len(values)
    assert all(v & (v - 1) == 0 for v in values)
    assert math.prod(values) > 0
