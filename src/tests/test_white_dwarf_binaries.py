# tests/test_white_dwarf_binaries.py

"""
GEN.271: the white-dwarf binary types of the member tree (`planetgen/generation/multistar.py`): close double white dwarfs,
double white dwarf merger products, post-common-envelope pairs with hot subdwarfs, cataclysmic variables, AM CVn binaries,
and barium and symbiotic giants. The audit is that of docs/design/multi-star-systems.md section 15, number 10 (type
counts inside their windows). No database or Redis: the draws are seeded.
"""

import math

import pytest

from planetgen import tuning
from planetgen.generation import multistar
from planetgen.generation.phenomena import system_types
from planetgen.physics import stability
from planetgen.util import draw

SIGMAS = 4.0
AGE_YR = 8.0e9
WHITE_DWARFS = 120_000
DAYS_PER_YEAR = 365.25


def _within(count, total, expected_share):
    """Whether `count` of `total` is within `SIGMAS` binomial sigmas of the share."""
    sigma = math.sqrt(total * expected_share * (1.0 - expected_share))
    return abs(count - total * expected_share) <= SIGMAS * sigma


def _hours(node):
    return node.period_years * DAYS_PER_YEAR * 24.0


@pytest.fixture(scope="module")
def white_dwarf_trees():
    """The trees of `WHITE_DWARFS` old white dwarfs of 0.6 solar masses born from 2 (the ones with a companion)."""
    rng = draw.Stream("gen271-white-dwarfs")
    placed = multistar.Placed("white-dwarf", 0.6, 2.0)
    trees = [multistar.draw_member_tree(placed, AGE_YR, rng) for _ in range(WHITE_DWARFS)]
    return [tree for tree in trees if tree.nodes], trees[0].members[0].flags


def _pairs(trees, interaction=None):
    return [tree for tree in trees if tree.nodes and (interaction is None or tree.nodes[0].interaction == interaction)]


# --- the type draws, with the uniforms forced ---

def test_a_white_dwarf_companion_is_a_close_double_white_dwarf_when_the_state_is_white_dwarf():
    assert multistar.draw_white_dwarf_binary("white-dwarf", "white-dwarf", AGE_YR, 0.99) == "double-white-dwarf"
    assert multistar.draw_white_dwarf_binary("white-dwarf", None, AGE_YR, 0.0) is None


def test_the_star_companion_types_follow_the_rows_rates_over_the_odds_of_a_star():
    rates = system_types.SYSTEM_TYPE_RATES
    star_odds = multistar.MEMBER_STATE_BY_PLACED["white-dwarf"]["star"][0]
    am_cvn = rates["am-cvn"].default / star_odds
    cv = rates["cataclysmic-variable"].default / star_odds
    assert multistar.draw_white_dwarf_binary("white-dwarf", "star", AGE_YR, am_cvn - 1e-9) == "am-cvn"
    assert multistar.draw_white_dwarf_binary("white-dwarf", "star", AGE_YR, am_cvn + 1e-9) == "cataclysmic-variable"
    assert multistar.draw_white_dwarf_binary("white-dwarf", "star", AGE_YR, am_cvn + cv + 1e-9) == "post-common-envelope"
    assert multistar.draw_white_dwarf_binary("white-dwarf", "star", AGE_YR, 0.99) is None


def test_a_post_common_envelope_pair_is_the_knob_share_of_white_dwarf_companions():
    table = multistar.MEMBER_STATE_BY_PLACED["white-dwarf"]
    companions = sum(default for default, _low, _high in table.values())
    chance = system_types.TUNING_KNOBS["post_ce_fraction"] * companions / table["star"][0]
    rates = system_types.SYSTEM_TYPE_RATES
    start = (rates["am-cvn"].default + rates["cataclysmic-variable"].default) / table["star"][0]
    assert multistar.draw_white_dwarf_binary("white-dwarf", "star", AGE_YR, start + chance - 1e-9) == (
        "post-common-envelope")
    assert multistar.draw_white_dwarf_binary("white-dwarf", "star", AGE_YR, start + chance + 1e-9) is None


def test_no_young_system_has_a_cataclysmic_variable_a_post_common_envelope_pair_or_an_am_cvn():
    assert multistar.draw_white_dwarf_binary("white-dwarf", "star", 1.0e7, 0.0) is None
    assert multistar.draw_white_dwarf_binary("white-dwarf", "star", 1.0e7, 0.1) is None
    # 2 Gy is intermediate: a cataclysmic variable, but no AM CVn (the old population's)
    assert multistar.draw_white_dwarf_binary("white-dwarf", "star", 2.0e9, 0.0) == "cataclysmic-variable"


def test_a_hot_subdwarf_is_the_stand_in_share_of_post_common_envelope_pairs():
    share = tuning.HOT_SUBDWARF_POST_CE_SHARE
    assert multistar.draw_stripped_helium(share - 1e-9)
    assert not multistar.draw_stripped_helium(share + 1e-9)


def test_an_old_white_dwarf_is_a_double_white_dwarf_merger_product_at_the_knob_share_and_never_when_young():
    share = system_types.TUNING_KNOBS["merger_product_fraction"]["white-dwarf"]
    assert multistar.draw_white_dwarf_merger(AGE_YR, share - 1e-9)
    assert not multistar.draw_white_dwarf_merger(AGE_YR, share + 1e-9)
    assert not multistar.draw_white_dwarf_merger(1.0e7, 0.0)


def test_a_g_or_k_giant_is_symbiotic_or_barium_at_the_rows_shares():
    entry = system_types.SYSTEM_TYPE_RATES["symbiotic-barium"]
    barium = entry.by_class["barium and CH (G and K giants)"]
    assert multistar.draw_giant_companion("giant", 1.2, entry.default - 1e-12) == "symbiotic"
    assert multistar.draw_giant_companion("giant", 1.2, entry.default + barium - 1e-9) == "barium"
    assert multistar.draw_giant_companion("giant", 1.2, entry.default + barium + 1e-9) is None
    assert multistar.draw_giant_companion("main-sequence", 1.2, 0.0) is None
    assert multistar.draw_giant_companion("giant", tuning.BARIUM_GIANT_MAX_INITIAL_MASS_SOL + 0.1, 0.0) is None


# --- the trees of white dwarfs: the audit counts ---

def test_a_white_dwarf_has_a_companion_at_thirty_percent(white_dwarf_trees):
    trees, _ = white_dwarf_trees
    assert _within(len(trees), WHITE_DWARFS, 0.30)


def test_close_double_white_dwarfs_are_9_5_percent_of_white_dwarfs_inside_4_au_and_alive_by_the_peters_rule(
        white_dwarf_trees):
    trees, _ = white_dwarf_trees
    pairs = [tree for tree in trees if all(member.kind == "white-dwarf" for member in tree.members)]
    assert _within(len(pairs), WHITE_DWARFS, system_types.SYSTEM_TYPE_RATES["close-double-white-dwarf"].default)
    low, high = tuning.DOUBLE_WHITE_DWARF_AXIS_AU
    for tree in pairs:
        node = tree.nodes[0]
        assert low <= node.separation_au <= high
        assert stability.survives_gravitational_waves(tree.members[0].mass_solar, tree.members[1].mass_solar,
                                                      node.separation_au, node.eccentricity, AGE_YR)


def test_double_white_dwarf_axes_are_log_uniform_up_to_the_peters_cut(white_dwarf_trees):
    trees, _ = white_dwarf_trees
    axes = sorted(tree.nodes[0].separation_au for tree in trees
                  if all(member.kind == "white-dwarf" for member in tree.members))
    low, high = tuning.DOUBLE_WHITE_DWARF_AXIS_AU
    assert axes[0] >= low and axes[-1] <= high
    # The Peters rule removes the pairs closer than the cut (the circular pair of the typical mass that has just merged by the
    # age); above it the axes are log-uniform, so the share above 0.1 AU is the decades above 0.1 over the decades above the cut.
    cut = 0.001
    while not stability.survives_gravitational_waves(0.6, 0.6, cut, 0.0, AGE_YR):
        cut *= 1.01
    expected = math.log(high / 0.1) / math.log(high / cut)
    wide = sum(1 for axis in axes if axis >= 0.1)
    assert _within(wide, len(axes), expected) or abs(wide / len(axes) - expected) < 0.03


def test_double_white_dwarf_merger_products_are_ten_percent_of_old_white_dwarfs_and_flagged_and_rejuvenated():
    rng = draw.Stream("gen271-mergers")
    placed = multistar.Placed("white-dwarf", 0.9, 4.0)
    flagged = 0
    for _ in range(20_000):
        member = multistar.draw_member_tree(placed, AGE_YR, rng).members[0]
        if member.flags & multistar.MEMBER_MERGER_PRODUCT:
            flagged += 1
            assert 0.0 < member.effective_age_yr < AGE_YR
        else:
            assert member.effective_age_yr is None
    assert _within(flagged, 20_000, system_types.TUNING_KNOBS["merger_product_fraction"]["white-dwarf"])


def test_cataclysmic_variables_are_6_5_in_10000_white_dwarfs_with_periods_1_3_to_10_hours_and_accreting(
        white_dwarf_trees):
    trees, _ = white_dwarf_trees
    low, high = tuning.CATACLYSMIC_PERIOD_HOURS
    variables = [tree for tree in trees if tree.nodes[0].interaction == "accreting"
                 and tree.members[1].kind == "main-sequence"]
    assert _within(len(variables), WHITE_DWARFS, system_types.SYSTEM_TYPE_RATES["cataclysmic-variable"].default)
    for tree in variables:
        assert low <= _hours(tree.nodes[0]) <= high
        assert tree.nodes[0].eccentricity == 0.0
        assert any(member.kind == "white-dwarf" for member in tree.members)
    magnetic = sum(1 for tree in variables
                   if next(m for m in tree.members if m.kind == "white-dwarf").flags & multistar.MEMBER_MAGNETIC)
    assert _within(magnetic, len(variables), tuning.CATACLYSMIC_MAGNETIC_SHARE)


def test_am_cvn_binaries_are_7_5_in_100000_white_dwarfs_with_a_helium_donor_and_periods_5_to_65_minutes(
        white_dwarf_trees):
    trees, _ = white_dwarf_trees
    low, high = tuning.AM_CVN_PERIOD_MINUTES
    binaries = [tree for tree in trees if tree.nodes[0].interaction == "accreting"
                and any(member.kind == "stripped-helium" for member in tree.members)]
    assert _within(len(binaries), WHITE_DWARFS, system_types.SYSTEM_TYPE_RATES["am-cvn"].default)
    for tree in binaries:
        assert low <= _hours(tree.nodes[0]) * 60.0 <= high


def test_a_white_dwarf_pair_is_accreting_only_as_a_cataclysmic_variable_or_an_am_cvn(white_dwarf_trees):
    trees, _ = white_dwarf_trees
    accreting = _pairs(trees, "accreting")
    assert accreting
    lowest = tuning.AM_CVN_PERIOD_MINUTES[0] / 60.0
    for tree in accreting:
        hours = _hours(tree.nodes[0])
        kinds = {member.kind for member in tree.members}
        assert kinds in ({"white-dwarf", "main-sequence"}, {"white-dwarf", "stripped-helium"})
        assert lowest <= hours <= tuning.CATACLYSMIC_PERIOD_HOURS[1]


def test_a_cataclysmic_variable_has_the_donor_inside_its_lobe_to_the_fill_factor(white_dwarf_trees):
    trees, _ = white_dwarf_trees
    for tree in _pairs(trees, "accreting")[:50]:
        node = tree.nodes[0]
        donor = next(member for member in tree.members if member.kind != "white-dwarf")
        accretor = next(member for member in tree.members if member.kind == "white-dwarf")
        fill = multistar.roche.filling_factor(donor.radius_au, node.separation_au, donor.mass_solar,
                                              accretor.mass_solar)
        assert fill == pytest.approx(tuning.ACCRETING_FILL_FACTOR, rel=1e-6)


# --- post-common-envelope pairs ---

def test_a_post_common_envelope_pair_has_its_axis_in_the_window_and_the_interaction(monkeypatch):
    monkeypatch.setattr(multistar, "draw_white_dwarf_binary", lambda *args: "post-common-envelope")
    rng = draw.Stream("gen271-post-ce")
    placed = multistar.Placed("white-dwarf", 0.7, 3.0)
    low, high = tuning.POST_COMMON_ENVELOPE_AXIS_AU
    seen = kinds = 0
    stripped = 0
    for _ in range(3_000):
        tree = multistar.draw_member_tree(placed, AGE_YR, rng, companion="star")
        if not tree.nodes:
            continue
        seen += 1
        node = tree.nodes[0]
        assert low <= node.separation_au <= high
        assert node.interaction in ("post-common-envelope", "detached")  # a donor that fills its lobe is a CV, not here
        assert node.interaction != "accreting"
        stripped += any(member.kind == "stripped-helium" for member in tree.members)
        kinds += 1
    assert seen > 2_500
    assert _within(stripped, kinds, tuning.HOT_SUBDWARF_POST_CE_SHARE)


# --- barium and symbiotic giants ---

@pytest.mark.parametrize("binary,flag,days", [
    ("barium", multistar.MEMBER_BARIUM, tuning.BARIUM_PERIOD_DAYS),
    ("symbiotic", multistar.MEMBER_SYMBIOTIC, tuning.SYMBIOTIC_PERIOD_DAYS)])
def test_a_barium_or_symbiotic_giant_has_a_white_dwarf_companion_on_the_types_orbit(monkeypatch, binary, flag, days):
    monkeypatch.setattr(multistar, "draw_giant_companion", lambda *args: binary)
    rng = draw.Stream(f"gen271-{binary}")
    mass = 1.2
    age_yr = multistar.evolution.main_sequence_lifetime_gy(mass) * 1.0e9 * 1.08
    found = 0
    for _ in range(300):
        tree = multistar.draw_member_tree(multistar.Placed("star", mass), age_yr, rng)
        giant = tree.members[0]
        if giant.kind != "giant" or not giant.flags & flag:
            continue
        found += 1
        dwarf = next(member for member in tree.members if member.kind == "white-dwarf")
        assert dwarf.initial_mass_solar > giant.initial_mass_solar
        node = next(n for n in tree.nodes if n.node_no == dwarf.parent_node_no)
        assert days[0] * 0.99 <= node.period_years * DAYS_PER_YEAR <= days[1] * 1.01
    assert found > 150


def test_a_giant_tree_is_barium_at_about_one_percent_of_g_and_k_giants():
    rng = draw.Stream("gen271-barium-share")
    giants = barium = 0
    for _ in range(12_000):
        mass = 1.0 + 0.5 * rng.random()
        age_yr = multistar.evolution.main_sequence_lifetime_gy(mass) * 1.0e9 * 1.08
        giant = multistar.draw_member_tree(multistar.Placed("star", mass), age_yr, rng).members[0]
        if giant.kind == "giant":
            giants += 1
            barium += bool(giant.flags & multistar.MEMBER_BARIUM)
    assert giants > 11_000
    # a white dwarf that cannot have formed yet, or a tree the stability screens dropped, takes a few of the draws
    assert 0.6 * 0.01 * giants <= barium <= 1.2 * 0.01 * giants


def test_the_new_flags_are_distinct_bits_of_the_member_flags_smallint():
    assert multistar.MEMBER_MAGNETIC | multistar.MEMBER_BARIUM | multistar.MEMBER_SYMBIOTIC < 2 ** 15
    assert {"magnetic", "barium", "symbiotic"} <= set(multistar.MEMBER_FLAGS)
