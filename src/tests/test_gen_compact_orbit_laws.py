# tests/test_gen_compact_orbit_laws.py

"""
GEN.266: the orbit law of each pairing with a compact object, as docs/design/systems/multi-star-systems.md 4.6 gives them.
The laws live in the member-tree generator (`planetgen/generation/multistar.py`, GEN.261); the Peters-time and
wide-companion filters and the V404 Cygni tertiary are tested in `test_multistar.py`. This pins the rest: a black
hole or neutron star with a star (log-uniform in period from 0.3 to 1e4 days, circular under 7 days, Roche floor),
a tight compact pair (log-uniform in axis from 0.003 to 100 AU, circular under 0.1 AU) and the relativistic quench
of the Kozai screen.
"""

import math
from collections import Counter

import pytest

from planetgen.generation import multistar
from planetgen.physics import roche, stability
from planetgen.util import draw

OLD = 8.0e9
YOUNG = 1.0e7


def _force(monkeypatch, state):
    monkeypatch.setattr(multistar, "companion_state", lambda placed, rng=draw: state)


def _pairs(placed, age, seeds, **kwargs):
    """The single node of every two-member tree drawn for `placed`, with its tree."""
    for seed in range(seeds):
        tree = multistar.draw_member_tree(placed, age, draw.Stream(seed), **kwargs)
        if tree.member_count == 2:
            yield tree, tree.nodes[0]


def _flat(values, decades, tolerance):
    """Whether the values fall evenly across the decades `decades` (a list of integer powers of ten)."""
    counts = Counter(math.floor(math.log10(v)) for v in values)
    mean = sum(counts[d] for d in decades) / len(decades)
    return all(abs(counts[d] - mean) <= tolerance * mean for d in decades), counts


@pytest.mark.parametrize("placed", [multistar.Placed("black-hole", 10.0, 45.0), multistar.Placed("neutron-star", 1.4, 10.0)])
def test_a_compact_object_with_a_star_draws_its_period_log_uniformly_from_0_3_to_1e4_days(monkeypatch, placed):
    _force(monkeypatch, "star")
    periods = [node.period_years * 365.25 for _tree, node in _pairs(placed, OLD, 2500, n_stars=1)]
    assert len(periods) > 1000
    assert min(periods) >= multistar.INTERACTED_PERIOD_DAYS[0] * (1 - 1e-9)
    assert max(periods) <= multistar.INTERACTED_PERIOD_DAYS[1] * (1 + 1e-9)
    flat, counts = _flat(periods, [0, 1, 2, 3], 0.2)
    assert flat, counts


def test_a_compact_object_with_a_star_is_circular_under_seven_days_and_eccentric_above(monkeypatch):
    _force(monkeypatch, "star")
    nodes = [node for _tree, node in _pairs(multistar.Placed("black-hole", 10.0, 45.0), OLD, 2000, n_stars=1)]
    short = [n for n in nodes if n.period_years * 365.25 < multistar.CIRCULAR_BELOW_DAYS]
    long = [n for n in nodes if n.period_years * 365.25 >= multistar.CIRCULAR_BELOW_DAYS]
    assert short and long
    assert all(n.eccentricity == 0.0 for n in short)
    assert all(0.0 < n.eccentricity < 1.0 for n in long)
    for n in long:
        assert n.eccentricity <= 1.0 - (n.period_years * 365.25 / 2.0) ** (-2.0 / 3.0) + 1e-12


def test_a_compact_object_with_a_star_is_marked_interacted_and_keeps_clear_of_the_stars_roche_lobe(monkeypatch):
    _force(monkeypatch, "star")
    checked = 0
    for tree, node in _pairs(multistar.Placed("black-hole", 10.0, 45.0), OLD, 800, n_stars=1):
        assert node.kind == "compact" and node.flags & multistar.NODE_INTERACTED
        star = next(m for m in tree.members if m.kind != "black-hole")
        hole = next(m for m in tree.members if m.kind == "black-hole")
        peri = node.separation_au * (1.0 - node.eccentricity)
        assert roche.filling_factor(star.radius_au, peri, star.mass_solar, hole.mass_solar) <= multistar.MAX_FILL_FACTOR
        checked += 1
    assert checked > 300


def test_two_stars_are_not_marked_interacted_and_use_the_star_law(monkeypatch):
    _force(monkeypatch, None)
    placed = multistar.Placed("star", 1.0)
    nodes = [node for _tree, node in _pairs(placed, 5e9, 1500, n_stars=2)]
    assert len(nodes) > 500
    assert all(n.kind in ("close", "wide") and not n.flags & multistar.NODE_INTERACTED for n in nodes)


def test_a_tight_compact_pair_draws_its_axis_log_uniformly_from_0_003_to_100_au(monkeypatch):
    _force(monkeypatch, "neutron-star")
    nodes = [node for _tree, node in _pairs(multistar.Placed("neutron-star", 1.4, 10.0), YOUNG, 3000)]
    assert len(nodes) > 1500
    axes = [n.separation_au for n in nodes]
    low, high = multistar.COMPACT_PAIR_AXIS_AU
    assert min(axes) >= low * (1 - 1e-9) and max(axes) <= high * (1 + 1e-9)
    flat, counts = _flat(axes, [-2, -1, 0, 1], 0.2)
    assert flat, counts
    assert all(n.kind == "compact" and n.flags & multistar.NODE_INTERACTED for n in nodes)


def test_a_tight_compact_pair_is_circular_under_a_tenth_of_an_au_and_eccentric_beyond(monkeypatch):
    _force(monkeypatch, "neutron-star")
    nodes = [node for _tree, node in _pairs(multistar.Placed("neutron-star", 1.4, 10.0), YOUNG, 3000)]
    inside = [n for n in nodes if n.separation_au < multistar.COMPACT_PAIR_CIRCULAR_BELOW_AU]
    beyond = [n for n in nodes if n.separation_au >= multistar.COMPACT_PAIR_CIRCULAR_BELOW_AU]
    assert inside and beyond
    assert all(n.eccentricity == 0.0 for n in inside)
    assert sum(n.eccentricity > 0.0 for n in beyond) > 0.99 * len(beyond)


def test_a_compact_pair_that_would_merge_inside_the_system_age_is_not_drawn(monkeypatch):
    _force(monkeypatch, "black-hole")
    placed = multistar.Placed("black-hole", 15.0, 45.0)
    for tree, node in _pairs(placed, OLD, 400):
        if all(m.kind == "black-hole" for m in tree.members):
            assert stability.peters_merger_time_yr(tree.members[0].mass_solar, tree.members[1].mass_solar,
                                                   node.separation_au, node.eccentricity) > OLD


@pytest.mark.parametrize("a_in, a_out", [(0.003, 1.0), (0.01, 20.0), (0.05, 20.0), (0.05, 300.0)])
def test_the_relativistic_precession_quenches_the_kozai_cycle_of_a_tight_compact_pair(a_in, a_out):
    verdict = stability.kozai_verdict(1.4, 1.4, 1.0, a_in, a_out, 0.3, 0.3, math.radians(60.0), 1.0e10, 1.0e-8)
    assert verdict == "safe"


def test_a_wider_pair_in_the_kozai_window_is_not_quenched():
    verdict = stability.kozai_verdict(1.4, 1.4, 1.0, 0.3, 20.0, 0.3, 0.3, math.radians(60.0), 1.0e10, 1.0e-8)
    assert verdict == "kl_active"
