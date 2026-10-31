# tests/test_truncated_sampling.py

"""
GEN.1011: a pair's period and eccentricity are drawn once from their laws cut to the room the pair has
(`multistar._Builder.place`), so no draw is repeated because it overlaps a partner or breaks the stability of the triple above it.
`draw.counting` counts the redraws that are left, by kind (`redraw_kinds.KINDS`).
"""

import math

import pytest

from planetgen.generation import multistar
from planetgen.generation.multistar import Side, pair_clears, single_side
from planetgen.util import draw


def _tree(count, seed):
    placed = multistar.Placed("star", 1.0 + 0.03 * (seed % 20), fe_h=0.0)
    return multistar.draw_member_tree(placed, (1.0 + 0.4 * (seed % 15)) * 1.0e9, draw.Stream(seed), n_stars=count)


def _builder(seed=1):
    return multistar._Builder([], 1.0e9, draw.Stream(seed))


def test_the_thermal_eccentricity_cut_at_a_limit_is_the_conditional_law():
    builder = _builder()
    draws = [builder.eccentricity("star", 400.0, 1.0, True, 0.3) for _ in range(4000)]
    assert max(draws) <= 0.3
    # the law is 2 e / cap^2 on [0, cap]: its mean is 2 cap / 3
    assert sum(draws) / len(draws) == pytest.approx(0.2, abs=0.01)
    assert builder.eccentricity("star", 2.0, 0.05, True, 0.3) == 0.0       # a pair under a week stays circular


def test_the_room_limit_follows_the_pericentre_the_budget_and_the_stability():
    builder = _builder()
    assert builder._eccentricity_limit(1.0, 0.4, math.inf, None) == pytest.approx(0.6)
    assert builder._eccentricity_limit(1.0, 0.4, 1.3, None) == pytest.approx(0.3)           # apocentre at most 1.3 AU
    assert builder._eccentricity_limit(1.0, 1.5, math.inf, None) < 0.0                      # the pericentre floor is past the axis
    assert builder._eccentricity_limit(1.0, 0.4, math.inf, lambda e: e <= 0.25) == pytest.approx(0.25, abs=1.0e-6)
    assert builder._eccentricity_limit(1.0, 0.4, math.inf, lambda e: False) == -1.0


def test_the_apocentre_budget_is_the_extent_a_node_may_have():
    builder = _builder()
    low, high = single_side(1.0, 0.01), single_side(0.25, 0.004)
    limit = builder._apocentre_limit_au(0.5, low, high)
    assert multistar.node_extent_au(limit, 0.0, low, high) == pytest.approx(0.5)
    assert builder._apocentre_limit_au(math.inf, low, high) == math.inf


@pytest.mark.parametrize("count", [2, 3, 4, 5, 6])
def test_no_draw_is_repeated_for_overlap_or_stability(count):
    with draw.counting() as redraws:
        trees = [_tree(count, seed) for seed in range(60)]
    assert all(redraws.get(kind, 0) == 0 for kind in ("pair.overlap", "pair.unstable_triple", "tree.restart_overlap")), redraws
    assert all(tree.member_count == count for tree in trees)
    for tree in trees:
        stars = [(m.parent_node_no, m.radius_au, m.mass_solar) for m in tree.members]
        assert multistar.built_tree_clears(tree.nodes, stars)


def test_the_redraws_left_are_the_named_kinds():
    with draw.counting() as redraws:
        for seed in range(40):
            _tree(4, seed)
    assert set(redraws) <= {"pair.no_allowed_period", "pair.axis_out_of_range", "tree.restart_pair_failed", "tree.member_dropped",
                            "pair.kozai", "pair.interaction_gate", "pair.gravitational_waves", "pair.brown_dwarf_desert"}, redraws


def test_a_pair_with_no_room_is_reported_not_forced():
    builder = multistar._Builder([multistar._Leaf("main-sequence", "V", 1.0, 1.0, 0.1), multistar._Leaf("main-sequence", "V", 1.0, 1.0, 0.1)],
                                 1.0e9, draw.Stream(3), max_axis_au=0.01)
    with draw.counting() as redraws:
        with pytest.raises(multistar._Failed):
            builder.place((0, 1), None, 0, None)
    assert redraws.get("pair.no_allowed_period", 0) + redraws.get("pair.axis_out_of_range", 0) == multistar.PAIR_RETRIES, redraws
    assert not builder.nodes
