"""
Bound compact pairs promoted to the scatter, whatever their mass (GEN.263, docs/design/systems/multi-star-systems.md 4.5): a neutron star with a neutron star or a black hole
with a black hole or a neutron star is a row of the scatter although the pair is under the mass cut, its tree drawn with a compact companion, and the sector's fill leaves those
objects out so the objects of a class add up to the whole population.
"""

import collections

import pytest

from planetgen import tuning
from planetgen.generation import multistar, phenomenon_scatter as scatter, scatter_companions
from planetgen.galaxy import remnant_density

NS_PAIR = ("neutron-star", "pair")
BH_PAIR = ("black-hole", "stellar-pair")
CUT = tuning.PHENOMENON_MIN_MASS_SOLAR


def test_the_pair_classes_come_before_the_magnetars_and_carry_the_base_class():
    classes = list(scatter.SCATTER_CLASSES)
    assert classes.index(NS_PAIR) < classes.index(BH_PAIR) < classes.index(("neutron-star", scatter.MAGNETAR_SUBTYPE))
    assert scatter.base_subtype(*NS_PAIR) is None and scatter.base_subtype(*BH_PAIR) == "stellar"
    assert scatter.base_subtype("black-hole", "intermediate") == "intermediate"
    assert remnant_density.follows_metals(*NS_PAIR) and remnant_density.follows_metals(*BH_PAIR)


def test_a_pair_row_is_drawn_below_the_cut_where_its_class_is_not_scattered():
    assert scatter.mass_range(*BH_PAIR, CUT, True) == (scatter.mass_law("black-hole", "stellar")[0][0], CUT)
    assert scatter.mass_range(*BH_PAIR, CUT, True) == scatter.mass_range("black-hole", "stellar", CUT, False)
    low, high = scatter.mass_law(*NS_PAIR)[0]
    assert scatter.mass_range(*NS_PAIR, CUT, True) == (low, high)      # a cut over the heaviest neutron star leaves them all below it
    assert scatter.share_above(*BH_PAIR, CUT) == pytest.approx(1.0 - scatter.share_above("black-hole", "stellar", CUT))
    assert scatter.share_above(*NS_PAIR, CUT) == pytest.approx(1.0)
    assert scatter.share_above(*BH_PAIR, 1.0) == 0.0      # a cut under the lightest black hole scatters them all, so there is no pair below it


def test_the_pair_rates_are_the_class_rate_times_the_below_cut_share_times_the_companion_odds():
    rates = {(kind, subtype): rate for kind, subtype, rate in scatter._scattered_rates(CUT)}
    share = scatter.NEUTRON_STAR_RECYCLED_SHARE
    expected_ns = ((1.0 - share) * multistar.MEMBER_STATE_BY_PLACED["neutron-star"]["neutron-star"][0]
                   + share * multistar.MEMBER_STATE_BY_PLACED["recycled-pulsar"]["neutron-star"][0])
    assert scatter.pair_odds(*NS_PAIR) == pytest.approx(expected_ns)
    assert scatter.pair_odds(*BH_PAIR) == pytest.approx(0.002 + 0.0005)
    assert rates[NS_PAIR] == pytest.approx(scatter.class_rate_per_star("neutron-star", None) * scatter.pair_odds(*NS_PAIR))
    below = 1.0 - scatter.share_above("black-hole", "stellar", CUT)
    assert rates[BH_PAIR] == pytest.approx(scatter.class_rate_per_star("black-hole", "stellar") * below * 0.0025)
    # About 1.7e5 neutron star pairs of 9.5e8 neutron stars (multi-star-systems.md 4.5)
    assert rates[NS_PAIR] / scatter.class_rate_per_star("neutron-star", None) == pytest.approx(1.8e-4, rel=0.2)


def test_the_fill_draws_the_objects_of_the_pairs_no_more():
    for kind, subtype in (("black-hole", "stellar"), ("neutron-star", None)):
        assert 0.99 < scatter.fill_share(kind, subtype) < 1.0
    assert scatter.fill_share("black-hole", "stellar") == pytest.approx(1.0 - (2 * 0.002 + 0.0005))
    assert scatter.fill_share("black-hole", "intermediate") == 1.0 and scatter.fill_share("planetary-nebula", None) == 1.0
    draws, _rng = scatter.below_cut_draws((9, 0, 0), (8200.0, 0.0, 0.0), _toy(), 1.0e8, CUT, 3)
    assert not [kind for kind, subtype, _range, _count in draws if scatter.is_pair(kind, subtype)]


def _toy():
    from tests.toy_shape import toy_shape
    import math
    return toy_shape(disk_scale_length_pc=40.0, disk_scale_height_pc=12.0, bulge_scale_radius_pc=10.0, bulge_amplitude=2.0, arm_count=2,
                     pitch_angle_rad=math.radians(15), arm_amplitude=0.4)


def test_every_pair_row_has_a_compact_companion_of_the_right_kind():
    states = {NS_PAIR: collections.Counter(), BH_PAIR: collections.Counter()}
    recycled = 0
    for seed in range(1, 6001):
        for pair in (NS_PAIR, BH_PAIR):
            placed_kind, was_recycled, state = scatter_companions.phenomenon_companion(*pair, seed)
            assert placed_kind == scatter.PLACED_KINDS[(pair[0], scatter.base_subtype(*pair))]
            states[pair][state] += 1
            recycled += was_recycled and pair == NS_PAIR
    assert set(states[NS_PAIR]) == {"neutron-star"}
    assert set(states[BH_PAIR]) <= {"black-hole", "neutron-star"}
    assert states[BH_PAIR]["black-hole"] / sum(states[BH_PAIR].values()) == pytest.approx(0.8, abs=0.03)
    assert scatter.pair_recycled_share(*NS_PAIR) == pytest.approx(
        scatter.NEUTRON_STAR_RECYCLED_SHARE * 2.0e-3 / scatter.pair_odds(*NS_PAIR))
    assert scatter.pair_recycled_share(*BH_PAIR) == 0.0


def test_a_pair_row_stores_a_tree_with_two_compact_members():
    columns = scatter.PHENOMENON_SCATTER_COLUMNS
    for kind, subtype in (NS_PAIR, BH_PAIR):
        low, high = scatter.mass_range(kind, subtype, CUT, True)
        row = [None] * len(columns)
        row[columns.index("kind")], row[columns.index("subtype")] = kind, subtype
        row[columns.index("mass_solar")] = min(high, max(low, 1.4 if kind == "neutron-star" else 8.0))
        trees = 0
        for seed in range(40, 70):
            row[columns.index("seed")] = seed
            tree = scatter_companions.phenomenon_tree(tuple(row))
            if tree is None:
                continue      # no room for it in the draw (the filters of multi-star-systems.md 3.3)
            trees += 1
            compact = [member for member in tree.members if member.kind in ("neutron-star", "black-hole")]
            assert len(compact) >= 2
        assert trees >= 10
