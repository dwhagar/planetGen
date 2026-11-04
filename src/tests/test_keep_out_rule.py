# tests/test_keep_out_rule.py

"""
GEN.1010: one keep-out rule at every depth of a star tree (`multistar.pair_floor_au`, `pair_clears`, `node_extent_au`), the same
rule for the tree as drawn, the build (`built_tree_clears`) and the independent verifier (`physics/clearance.py`).
"""

import dataclasses

import pytest

from planetgen.generation import multistar
from planetgen.generation.multistar import Side, node_extent_au, pair_clears, pair_floor_au, single_side
from planetgen.physics import clearance, constants, stability
from planetgen.util import draw

AU_KM = constants.AU_TO_KM
SOLAR_KG = constants.SOLAR_MASS_TO_KG


def _tree(count, seed):
    placed = multistar.Placed("star", 1.0 + 0.03 * (seed % 20), fe_h=0.0)
    return multistar.draw_member_tree(placed, (1.0 + 0.4 * (seed % 15)) * 1.0e9, draw.Stream(seed), n_stars=count)


def _rows(tree):
    stars = [{"id": m.member_no, "mass_kg": m.mass_solar * SOLAR_KG, "radius_km": m.radius_au * AU_KM,
              "parent_node_no": m.parent_node_no, "slot": m.slot} for m in tree.members]
    return {"stars": stars, "nodes": [dataclasses.asdict(node) for node in tree.nodes]}


def test_two_singles_are_held_by_the_roche_lobe_and_other_pairs_by_the_extent_rule():
    star, partner = single_side(1.0, 0.01), single_side(0.5, 0.004)
    floor = pair_floor_au(star, partner)
    assert pair_clears(floor * 1.001, 0.0, star, partner) and not pair_clears(floor * 0.9, 0.0, star, partner)
    node = Side(1.5, 0.3)
    outer = single_side(0.8, 0.005)
    assert pair_floor_au(node, outer) == pytest.approx(multistar.EXTENT_AXIS_FACTOR * (0.3 + 0.005))
    assert pair_clears(1.0, 0.0, node, outer)
    assert not pair_clears(1.0, 0.5, node, outer)           # the pericentre is 0.5 AU, under 2 x 0.305


def test_the_extent_of_a_node_is_the_barycentric_one_the_verifier_uses():
    low, high = single_side(1.0, 0.01), single_side(0.25, 0.004)
    extent = node_extent_au(0.2, 0.3, low, high)
    rows = {"stars": [{"id": 0, "mass_kg": SOLAR_KG, "radius_km": 0.01 * AU_KM, "parent_node_no": 1, "slot": 0},
                      {"id": 1, "mass_kg": 0.25 * SOLAR_KG, "radius_km": 0.004 * AU_KM, "parent_node_no": 1, "slot": 1}],
            "nodes": [{"node_no": 1, "parent_node_no": None, "slot": 0, "mass_solar": 1.25, "separation_km": 0.2 * AU_KM,
                       "eccentricity": 0.3}]}
    tree = clearance._Tree(rows["stars"], rows["nodes"])
    assert extent == pytest.approx(tree.extent_au(("node", rows["nodes"][0])), rel=1.0e-12)


@pytest.mark.parametrize("count", [2, 3, 4, 5])
def test_a_drawn_tree_clears_by_the_build_rule_and_the_verifier(count):
    checked = 0
    for seed in range(40):
        tree = _tree(count, seed)
        if tree.member_count < 2:
            continue
        stars = [(m.parent_node_no, m.radius_au, m.mass_solar) for m in tree.members]
        assert multistar.built_tree_clears(tree.nodes, stars), (count, seed)
        bad = [item for item in clearance.audit(_rows(tree)) if item.hard and item.kind in ("overlap", "roche", "stability")]
        assert not bad, (count, seed, bad[:2])
        checked += 1
    assert checked >= 20


def test_the_stability_margin_holds_the_inner_orbit_in():
    bare = stability.max_inner_axis(1.0, 0.8, 0.6, 100.0, 0.0, 0.3, 0.4)
    margin = stability.max_inner_axis(1.0, 0.8, 0.6, 100.0, 0.0, 0.3, 0.4, safety=multistar.STABILITY_MARGIN)
    assert margin == pytest.approx(bare / multistar.STABILITY_MARGIN, rel=0.05)
    assert stability.is_triple_stable(1.0, 0.8, 0.6, margin * 0.99, 100.0, 0.0, 0.3, 0.4, safety=multistar.STABILITY_MARGIN)
    assert not stability.is_triple_stable(1.0, 0.8, 0.6, bare * 0.99, 100.0, 0.0, 0.3, 0.4, safety=multistar.STABILITY_MARGIN)
