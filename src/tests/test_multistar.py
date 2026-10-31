# tests/test_multistar.py

"""
GEN.261: the member-tree generator (`planetgen/generation/multistar.py`).
The odds are those of docs/design/systems/multi-star-systems.md section 4 and
docs/design/systems/multistar-and-compact-systems.md 3.3; the stability, Peters and
wide-companion rules are section 3.3 and 9 of the first note.
"""

import math
from collections import Counter

import pytest

from planetgen import tuning
from planetgen.generation import close_binary_trend, mass_ratio, multistar
from planetgen.generation.phenomena import hosts
from planetgen.physics import brown_dwarf, mass_boundaries, ms_lifetime, stability, stellar_evolution as evolution
from planetgen.util import draw

OLD = 8.0e9


def _tree(placed, age, seed, **kwargs):
    return multistar.draw_member_tree(placed, age, draw.Stream(seed), **kwargs)


def _children(tree, node_no):
    """The members and nodes that are children of the node, by slot: `{slot: (is_node, record)}`."""
    out = {}
    for member in tree.members:
        if member.parent_node_no == node_no:
            out[member.slot] = (False, member)
    for node in tree.nodes:
        if node.parent_node_no == node_no:
            out[node.slot] = (True, node)
    return out


def _mass(child):
    return child[1].mass_solar


def _check_structure(tree, placed):
    assert tree.members[0].kind == (placed.kind if placed.kind != "intermediate-mass-black-hole" else "black-hole") \
        or placed.kind == "star"
    assert [member.member_no for member in tree.members] == list(range(tree.member_count))
    assert len(tree.nodes) == tree.member_count - 1
    if not tree.nodes:
        assert tree.member_count == 1 and tree.members[0].parent_node_no is None
        return
    roots = [node for node in tree.nodes if node.parent_node_no is None]
    assert [node.node_no for node in roots] == [1]
    for node in tree.nodes:
        kids = _children(tree, node.node_no)
        assert sorted(kids) == [0, 1]
        assert node.mass_solar == pytest.approx(_mass(kids[0]) + _mass(kids[1]))
        assert _mass(kids[0]) >= _mass(kids[1])
        assert node.separation_km == pytest.approx(node.separation_au * 149_597_870.7)
        assert 0.0 <= node.eccentricity < 1.0
        assert node.period_years == pytest.approx(stability.orbital_period_yr(node.separation_au, node.mass_solar))
    assert tree.system_mass_solar == pytest.approx(tree.nodes[0].mass_solar)


def _check_stable(tree):
    for node in tree.nodes:
        kids = _children(tree, node.node_no)
        for child_slot in (0, 1):
            is_node, child = kids[child_slot]
            sibling = kids[1 - child_slot]
            if is_node:
                inner = _children(tree, child.node_no)
                m1, m2 = _mass(inner[0]), _mass(inner[1])
                mutual = stability.mutual_inclination(
                    math.radians(child.inclination_deg), math.radians(child.ascending_node_deg),
                    math.radians(node.inclination_deg), math.radians(node.ascending_node_deg))
                assert stability.is_triple_stable(m1, m2, _mass(sibling), child.separation_au, node.separation_au,
                                                  child.eccentricity, node.eccentricity, mutual)
        assert node.separation_au * (1.0 - node.eccentricity) > 0.0


def test_a_single_object_has_one_member_and_no_orbit():
    tree = _tree(multistar.Placed("star", 1.0), 4.6e9, 1, n_stars=1)
    assert tree.member_count == 1 and tree.nodes == () and tree.dropped == 0
    assert tree.members[0].kind == "main-sequence" and tree.members[0].parent_node_no is None


def test_a_forced_single_object_has_no_companion_whatever_it_became():
    """GEN.719: n_stars=1 (no binaries, one star) is one member for a star that is now a white dwarf or a neutron star,
    and for a placed white dwarf, neutron star or black hole, where the companion table used to add one anyway."""
    placed = [(multistar.Placed("star", 1.0), 14.0e9), (multistar.Placed("star", 1.0), 4.6e9),
              (multistar.Placed("white-dwarf", 0.6, 1.5), 5.0e9), (multistar.Placed("neutron-star", 1.4, 15.0), 5.0e9),
              (multistar.Placed("black-hole", 10.0, 40.0), 5.0e9), (multistar.Placed("star", 0.4), 13.0e9)]
    for object_, age in placed:
        for seed in range(150):
            tree = _tree(object_, age, seed, n_stars=1)
            assert tree.member_count == 1 and tree.nodes == (), (object_, age, seed)


def test_the_same_stream_gives_the_same_tree_and_another_gives_another():
    placed = multistar.Placed("star", 20.0)
    assert _tree(placed, 1.0e6, 5) == _tree(placed, 1.0e6, 5)
    assert _tree(placed, 1.0e6, 5, n_stars=6) != _tree(placed, 1.0e6, 6, n_stars=6)


@pytest.mark.parametrize("placed, age", [
    (multistar.Placed("star", 20.0), 1.0e6),
    (multistar.Placed("star", 1.0), 4.6e9),
    (multistar.Placed("star", 0.3), 5.0e9),
    (multistar.Placed("black-hole", 15.0, 45.0), OLD),
    (multistar.Placed("neutron-star", 1.4, 10.0, recycled=True), 3.0e9),
    (multistar.Placed("white-dwarf", 0.6, 1.5), 5.0e9),
])
def test_trees_are_well_formed_and_every_level_is_stable(placed, age):
    seen = Counter()
    for seed in range(300):
        tree = _tree(placed, age, seed)
        _check_structure(tree, placed)
        _check_stable(tree)
        seen[tree.member_count] += 1
    assert len(seen) >= 1


def test_every_member_kind_is_known_and_masses_are_positive():
    for seed in range(300):
        tree = _tree(multistar.Placed("star", 20.0), 1.0e6, seed)
        for member in tree.members:
            assert member.kind in multistar.MEMBER_KINDS
            assert member.mass_solar > 0 and member.radius_au > 0


def test_the_number_of_stars_follows_the_class_table(monkeypatch):
    monkeypatch.setattr(multistar, "draw_merger_product", lambda *args: None)  # a merger has one star fewer (GEN.270)
    placed = multistar.Placed("star", 20.0)
    counts = Counter(_tree(placed, 1.0e6, seed).star_count for seed in range(4000))
    table = hosts.STAR_COUNT_PROBABILITIES["O"]
    for n in (1, 2, 3):
        expected = table[n - 1]
        sigma = math.sqrt(expected * (1 - expected) / 4000)
        assert abs(counts[n] / 4000 - expected) < 5 * sigma + 0.01, n
    assert max(counts) <= 12


def test_the_star_count_can_be_forced_capped_and_scaled():
    placed = multistar.Placed("star", 20.0)
    assert {_tree(placed, 1.0e6, seed, n_stars=4).star_count for seed in range(40)} <= {4, 3}
    assert max(_tree(placed, 1.0e6, seed, max_stars=3).member_count for seed in range(300)) <= 3
    flat = Counter(_tree(placed, 1.0e6, seed, prevalence=0.0).star_count for seed in range(600))
    assert max(flat) <= 4
    rich = Counter(_tree(placed, 1.0e6, seed, prevalence=5.0).star_count for seed in range(600))
    base = Counter(_tree(placed, 1.0e6, seed).star_count for seed in range(600))
    assert sum(n for k, n in rich.items() if k >= 5) > sum(n for k, n in base.items() if k >= 5)


def test_the_tail_reaches_eight_stars_and_stops_at_the_cap():
    placed = multistar.Placed("star", 20.0)
    seen = Counter(_tree(placed, 1.0e6, seed, prevalence=20.0).member_count for seed in range(1500))
    assert max(seen) >= 8 and max(seen) <= 12


def test_members_take_their_state_from_the_age():
    # GEN.242: a star of 8 solar masses or more has no contracting stage; the lighter ones are protostars.
    young = _tree(multistar.Placed("star", 6.0), 2.0e5, 3, n_stars=3)
    assert {member.kind for member in young.members} == {"protostar"}
    massive = _tree(multistar.Placed("star", 20.0), 2.0e5, 3, n_stars=3)
    assert {m.kind for m in massive.members if m.initial_mass_solar >= 8.0} <= {"main-sequence"}
    assert {m.kind for m in massive.members if m.initial_mass_solar < 8.0} <= {"protostar"}
    mid = {m.kind for seed in range(80) for m in _tree(multistar.Placed("star", 20.0), 1.0e6, seed, n_stars=4).members}
    assert "pre-main-sequence" in mid and "main-sequence" in mid
    old = {m.kind for seed in range(80) for m in _tree(multistar.Placed("star", 1.0), 12.0e9, seed, n_stars=3).members}
    assert "pre-main-sequence" not in old and "protostar" not in old and "giant" in old | {"giant"}


def test_a_collapsed_placed_star_is_refused():
    with pytest.raises(ValueError):
        _tree(multistar.Placed("star", 25.0), 2.0 * ms_lifetime.lifetime_gy(25.0) * 1.0e9, 1)
    with pytest.raises(ValueError):
        multistar.placed_class(multistar.Placed("planet", 1.0))


def test_brown_dwarfs_keep_clear_of_the_stars_they_orbit():
    inside, total = 0, 0
    for seed in range(4000):
        tree = _tree(multistar.Placed("star", 0.5), 3.0e9, seed)
        for member in tree.members[1:]:
            if member.kind == "brown-dwarf" and tree.member_count == 2:     # a pair of the star and the dwarf
                total += 1
                inside += tree.nodes[0].separation_au < multistar.BROWN_DWARF_DESERT_AU
    assert total > 20 and inside / total < 0.1


@pytest.mark.parametrize("birth", [0.12, 0.25, 0.5, 1.0])
def test_the_stars_of_a_tree_are_the_tables_count_and_a_brown_dwarf_is_one_more_member(birth):
    """GEN.899: a companion under the hydrogen-burning limit used to use up one of the table's places (an M2V primary kept
    13% of its 27%). Every star of a tree is a star; the brown dwarfs come at `tuning.BROWN_DWARF_COMPANION_CHANCE`."""
    draws = 3000
    table = multistar.star_count_weights(birth)
    stars, dwarfs = Counter(), 0
    for seed in range(draws):
        tree = _tree(multistar.Placed("star", birth), 3.0e9, seed)
        kinds = Counter(member.kind == "brown-dwarf" for member in tree.members[1:])
        dwarfs += kinds[True]
        stars[1 + kinds[False]] += 1
        assert all(member.initial_mass_solar >= mass_boundaries.sol("hydrogen_burning_limit")
                   for member in tree.members[1:] if member.kind != "brown-dwarf")
    share = 1.0 - table[0] / sum(table)
    assert abs((1.0 - stars[1] / draws) - share) < 4.0 * math.sqrt(share * (1.0 - share) / draws) + 0.005
    chance = tuning.BROWN_DWARF_COMPANION_CHANCE
    assert abs(dwarfs / draws - chance) < 4.0 * math.sqrt(chance * (1.0 - chance) / draws) + 0.005


def test_a_forced_count_is_of_stars_and_adds_no_brown_dwarf():
    for seed in range(300):
        tree = _tree(multistar.Placed("star", 0.2), 3.0e9, seed, n_stars=2)
        assert tree.member_count == 2 and tree.members[1].kind != "brown-dwarf"


def test_the_companion_table_is_consistent():
    for name, table in multistar.MEMBER_STATE_BY_PLACED.items():
        assert sum(default for default, _, _ in table.values()) <= 1.0, name
        for state, (default, low, high) in table.items():
            assert state in multistar.MEMBER_KINDS or state in multistar.COMPANION_STATES, (name, state)
            assert 0.0 <= low <= default <= high <= 1.0, (name, state)
    black_hole = multistar.MEMBER_STATE_BY_PLACED["black-hole"]
    assert sum(default for default, _, _ in black_hole.values()) == pytest.approx(0.0125)


@pytest.mark.parametrize("state, expected", [("star", 0.010), ("black-hole", 0.002), ("neutron-star", 0.0005)])
def test_a_black_hole_has_a_companion_at_the_table_odds(state, expected):
    draws = 300_000
    rng = draw.Stream(11)
    placed = multistar.Placed("black-hole", 15.0)
    hits = sum(multistar.companion_state(placed, rng) == state for _ in range(draws))
    sigma = math.sqrt(draws * expected * (1 - expected))
    assert abs(hits - draws * expected) < 4 * sigma


def test_a_neutron_star_pair_odds_and_recycled_pulsars_mostly_have_white_dwarfs():
    rng = draw.Stream(3)
    recycled = multistar.Placed("neutron-star", 1.4, recycled=True)
    states = Counter(multistar.companion_state(recycled, rng) for _ in range(20_000))
    assert states[None] / 20_000 == pytest.approx(0.22, abs=0.02)
    assert states["white-dwarf"] > states["star"] > states["brown-dwarf"]
    young = Counter(multistar.companion_state(multistar.Placed("neutron-star", 1.4), rng) for _ in range(20_000))
    assert young[None] / 20_000 == pytest.approx(0.98, abs=0.01)


def test_a_living_massive_star_gets_a_compact_companion_only_in_the_last_third_of_its_life():
    placed = multistar.Placed("star", 20.0)
    life_yr = ms_lifetime.lifetime_gy(20.0) * evolution.remnant_age_fraction(20.0) * 1.0e9
    late = [_tree(placed, 0.75 * life_yr, seed) for seed in range(6000)]
    early = [_tree(placed, 0.15 * life_yr, seed) for seed in range(2000)]
    compact = sum(tree.compact_count > 0 for tree in late)
    assert 0.5 * 0.009 * 6000 < compact < 2.0 * 0.009 * 6000
    assert sum(tree.compact_count > 0 for tree in early) == 0
    for tree in late:
        if tree.compact_count:
            assert tree.compact_count == 1


def test_a_giant_can_have_a_white_dwarf_companion():
    placed = multistar.Placed("star", 3.0)
    lifetime = ms_lifetime.lifetime_gy(3.0) * 1.0e9 * 1.2
    kinds = Counter(m.kind for seed in range(6000) for m in _tree(placed, lifetime * 0.9, seed).members[1:])
    assert kinds["white-dwarf"] > 0


def _force(monkeypatch, state):
    monkeypatch.setattr(multistar, "companion_state", lambda placed, rng=draw: state)


def test_compact_pairs_never_have_merged_by_gravitational_waves(monkeypatch):
    _force(monkeypatch, "neutron-star")
    placed = multistar.Placed("neutron-star", 1.4, 10.0)
    pairs = 0
    for seed in range(600):
        tree = _tree(placed, 10.0e9, seed)
        if tree.member_count >= 2:
            pairs += 1
            node = tree.nodes[0]
            assert stability.peters_merger_time_yr(
                tree.members[0].mass_solar, tree.members[1].mass_solar, node.separation_au, node.eccentricity) > 10.0e9
    assert pairs > 300


def test_a_young_compact_pair_may_be_as_tight_as_the_double_pulsar(monkeypatch):
    _force(monkeypatch, "neutron-star")
    placed = multistar.Placed("neutron-star", 1.4, 10.0)
    old_min = min(_tree(placed, 10.0e9, s).nodes[0].separation_au for s in range(1500))
    young_min = min(_tree(placed, 1.0e7, s).nodes[0].separation_au for s in range(1500))
    assert young_min < 0.0059 < old_min


def test_a_wide_companion_of_a_hole_that_did_not_collapse_directly_is_refused():
    light = multistar._Leaf("black-hole", None, 14.0, 25.0, 1e-9)
    direct = multistar._Leaf("black-hole", None, 14.0, 45.0, 1e-9, multistar.MEMBER_DIRECT_COLLAPSE)
    star = multistar._Leaf("main-sequence", "V", 0.9, 0.9, 0.005)
    for leaf, cap in ((light, multistar.WIDE_COMPACT_AXIS_AU), (direct, math.inf)):
        builder = multistar._Builder([leaf, star, star], 5e9, draw.Stream(1))
        assert builder.wide_cap(((0, 1), 2)) == cap
        assert builder.wide_cap(2) == math.inf


def test_black_hole_trees_carry_wide_tertiaries_only_around_direct_collapse_holes(monkeypatch):
    _force(monkeypatch, "star")
    wide = {}
    for label, birth in (("light", 30.0), ("direct", 40.0)):
        placed = multistar.Placed("black-hole", 12.0, birth)
        axes = []
        for seed in range(1500):
            tree = _tree(placed, OLD, seed)
            axes.extend(node.separation_au for node in tree.nodes)
        wide[label] = max(axes)
    assert wide["light"] <= multistar.WIDE_COMPACT_AXIS_AU
    assert wide["direct"] > 3 * multistar.WIDE_COMPACT_AXIS_AU


def test_a_black_hole_binary_tertiary_is_stable_like_v404_cygni(monkeypatch):
    _force(monkeypatch, "star")
    placed = multistar.Placed("black-hole", 9.0, 45.0)
    found = 0
    for seed in range(1500):
        tree = _tree(placed, OLD, seed)
        if tree.member_count == 3:
            found += 1
            _check_structure(tree, placed)
            _check_stable(tree)
    assert found > 100


def test_compact_members_stay_within_three(monkeypatch):
    _force(monkeypatch, "black-hole")
    placed = multistar.Placed("black-hole", 15.0, 45.0)
    for seed in range(500):
        assert _tree(placed, OLD, seed).compact_count <= multistar.MAX_COMPACT_MEMBERS


def test_scattered_objects_carry_no_planets_or_belts():
    kinds = set()
    for placed, age in ((multistar.Placed("star", 20.0), 1.0e6), (multistar.Placed("black-hole", 15.0, 45.0), OLD)):
        for seed in range(200):
            kinds.update(member.kind for member in _tree(placed, age, seed).members)
    assert kinds <= set(multistar.MEMBER_KINDS)


def test_the_b_star_rows_reproduce_moe_and_di_stefanos_table_13_fractions():
    # GEN.449: single, binary and triple-or-more fractions of the mid-B and early-B bins.
    for name, row in (("mid B", "B-mid"), ("early B", "B-early")):
        bin_ = tuning.MULTIPLICITY_BY_PRIMARY_MASS_BIN[name]
        table = hosts.STAR_COUNT_PROBABILITIES[row]
        assert sum(table) == pytest.approx(1.0, abs=1e-4)
        assert table[0] == pytest.approx(bin_["single"], abs=1e-4)
        assert table[1] == pytest.approx(bin_["binary"], abs=1e-4)
        assert sum(table[2:]) == pytest.approx(bin_["multiple"], abs=1e-4)


def test_the_table_13_bins_hold_their_own_totals_and_rise_with_mass():
    bins = tuning.MULTIPLICITY_BY_PRIMARY_MASS_BIN.values()
    for bin_ in bins:
        assert bin_["single"] + bin_["binary"] + bin_["multiple"] == pytest.approx(1.0, abs=0.011)
    ordered = sorted(bins, key=lambda b: b["masses"][0])
    assert [b["single"] for b in ordered] == sorted((b["single"] for b in ordered), reverse=True)
    assert [b["multiplicity"] for b in ordered] == sorted(b["multiplicity"] for b in ordered)


@pytest.mark.parametrize("mass, letter, key", [(0.5, "K", "K"), (3.0, "B", "B"), (4.99, "B", "B"), (5.0, "B", "B-mid"),
                                                (8.9, "B", "B-mid"), (9.0, "B", "B-early"), (15.9, "B", "B-early"),
                                                (16.0, "O", "O"), (1.5, "A", "A")])
def test_the_companion_class_splits_the_b_stars_by_birth_mass(mass, letter, key):
    assert multistar.companion_class(mass, letter) == key
    assert multistar.companion_class(mass) == (key if letter != "B" else multistar.companion_class(mass, "B"))


def test_the_chance_of_a_companion_rises_through_the_b_sub_classes_to_o():
    chain = [1.0 - hosts.STAR_COUNT_PROBABILITIES[k][0] / sum(hosts.STAR_COUNT_PROBABILITIES[k])
             for k in ("A", "B", "B-mid", "B-early", "O")]
    assert chain == sorted(chain) and chain[2] == pytest.approx(0.76) and chain[3] == pytest.approx(0.84)


def test_a_primary_near_the_hydrogen_burning_limit_never_draws_a_companion_under_the_deuterium_limit():
    """GEN.596: `birth * U(0.1, 1)` reached 0.0085 solar masses for a 0.085 primary, below the lowest cooling track."""
    floor = brown_dwarf.mass_range_solar()[0]
    for n in range(3000):
        rng = draw.Stream(f"f{n}".encode())
        assert _tree_draw_ok(multistar.Placed("star", 0.085), rng, n)
    for birth in (0.015, floor, 0.02, 0.085, 0.124, 0.5):
        rng = draw.Stream(596)
        masses = [multistar.companion_mass(birth, rng) for _ in range(200)]
        assert all(min(birth, floor) <= mass <= birth for mass in masses)
        assert all(mass >= 0.1 * birth - 1e-12 for mass in masses)


def _tree_draw_ok(placed, rng, n):
    tree = multistar.draw_member_tree(placed, 5.0e9, rng)
    return tree.member_count >= 1 and all(member.mass_solar >= 0.0 for member in tree.members)


def test_a_companion_is_the_primarys_mass_times_a_ratio_of_the_mass_ratio_law():
    """GEN.951: the draw takes one random number and gives it to `mass_ratio.draw`; a heavier primary keeps the whole range."""
    rng, again = draw.Stream(1), draw.Stream(1)
    assert multistar.companion_mass(1.0, rng) == 1.0 * mass_ratio.draw(1.0, again.random())
    assert multistar.companion_ratio_range(1.0) == tuning.BINARY_MASS_RATIO_RANGE


# GEN.705: the [Fe/H] factor on the close companions reaches the tree draw.

def test_the_tree_draw_leaves_the_table_alone_for_a_star_at_the_solar_composition():
    assert multistar.close_trend(1.0, None) == (None, 1.0)
    assert multistar.close_trend(1.0, 0.0) == (None, 1.0)


@pytest.mark.parametrize("mass", [0.3, 1.0, 5.0, 20.0])
@pytest.mark.parametrize("fe_h", [-2.0, -0.5, 0.5])
def test_close_trend_moves_the_odds_of_a_companion_inside_10_au_by_the_factor_and_leaves_the_rest(mass, fe_h):
    table = multistar.star_count_weights(mass)
    weights, inside_weight = multistar.close_trend(mass, fe_h)
    chance, moved = 1.0 - table[0] / sum(table), 1.0 - weights[0] / sum(weights)
    inside = multistar._inside_share(mass)
    outside = chance * (1.0 - inside)
    inside_after = close_binary_trend.scaled_close_chance(chance * inside, fe_h)
    assert moved == pytest.approx(min(inside_after + outside, 1.0), abs=1e-8)      # the outside part gives way past 1
    # the counts of 2 to 7 stars keep their proportions, and the share of the companions inside moves with the chance
    assert [w / sum(weights[1:]) for w in weights[1:]] == pytest.approx([w / sum(table[1:]) for w in table[1:]])
    odds = lambda p: p / (1.0 - p)
    assert odds(inside_after / moved) / odds(inside) == pytest.approx(inside_weight, rel=1e-9)
    assert (inside_weight > 1.0) == (fe_h < 0.0)


@pytest.mark.parametrize("share", [0.1, 0.5, 0.9])
@pytest.mark.parametrize("weight", [0.4, 1.0, 2.5])
def test_a_tilted_draw_puts_the_inside_share_at_the_weighted_odds(share, weight):
    draws = [(i + 0.5) / 4000 for i in range(4000)]
    inside = sum(multistar._tilted(u, share, weight) < share for u in draws) / len(draws)
    odds = share / (1.0 - share) * weight
    assert inside == pytest.approx(odds / (1.0 + odds), abs=1e-3)
    assert multistar._tilted(0.0, share, weight) == 0.0 and multistar._tilted(1.0, share, weight) == pytest.approx(1.0)


def test_metal_poor_stars_have_more_companions_inside_10_au_than_metal_rich_ones():
    def rates(fe_h, n=250):
        multiple = close = 0
        for seed in range(n):
            tree = _tree(multistar.Placed("star", 1.0, 1.0, fe_h=fe_h), 4.0e9, seed)
            if tree.member_count > 1:
                multiple += 1
                close += tree.nodes[0].separation_au < multistar.CLOSE_COMPANION_AXIS_AU
        return multiple / n, close / n
    poor, rich = rates(-2.0), rates(0.5)
    assert poor[0] > rich[0] + 0.1      # the chance of a companion (0.60 and 0.41 at the table's 0.47)
    assert poor[1] > rich[1] + 0.1      # and of one inside 10 AU (0.29 and 0.08)


def test_the_lightest_companion_follows_the_places_own_hydrogen_burning_limit():
    """GEN.706: at low [Fe/H] the limit between a star and a brown dwarf is above 0.08 (GEN.936), and a companion drawn between the two
    would be a brown dwarf using up a place of the `star-count` table."""
    from planetgen.physics import stellar_abundances
    assert multistar.companion_floor_solar(None) == multistar.companion_floor_solar(0.3) == mass_boundaries.sol("hydrogen_burning_limit")
    assert multistar.companion_floor_solar(-1.0) == pytest.approx(stellar_abundances.hydrogen_burning_limit_for_fe_h(-1.0))
    assert multistar.companion_floor_solar(-1.0) > 0.09
    for seed in range(300):
        tree = _tree(multistar.Placed("star", 0.3, 0.3, fe_h=-1.0), 3.0e9, seed, n_stars=2)
        assert tree.members[1].kind != "brown-dwarf" and tree.members[1].initial_mass_solar >= multistar.companion_floor_solar(-1.0)


@pytest.mark.parametrize("birth,age_yr", [(1.3, 3.0e9), (2.0, 5.0e8), (5.0, 5.0e7), (15.0, 3.0e6), (40.0, 1.0e6)])
def test_the_stars_of_a_tree_follow_the_table_for_the_hotter_classes_too(birth, age_yr):
    """GEN.706: the share of systems with a companion star is the `star-count` table's for every class, to four standard errors of
    the draws (a brown dwarf is a member of its own, not a star of the count)."""
    draws = 1500
    share = 1.0 - multistar.star_count_weights(birth)[0] / sum(multistar.star_count_weights(birth))
    pairs = sum(sum(member.kind != "brown-dwarf" for member in _tree(multistar.Placed("star", birth), age_yr, seed).members) > 1
                for seed in range(draws))
    assert abs(pairs / draws - share) < 4.0 * math.sqrt(share * (1.0 - share) / draws) + 0.02     # +0.02: the compact companions
