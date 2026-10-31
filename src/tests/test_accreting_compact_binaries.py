# tests/test_accreting_compact_binaries.py

"""
GEN.272: the accreting compact binaries of the member tree (`planetgen/generation/multistar.py`, `planetgen/physics/accretion.py`):
high-mass X-ray binaries (supergiant and Be/X-ray), low-mass X-ray binaries, the microquasar flag, the spider millisecond
pulsars, the dormant neutron stars and black holes of the Gaia kind, and the compact triples. The audit is that of
docs/design/systems/multi-star-systems.md section 15, number 10 (type counts inside their windows). No database or Redis: the
draws are seeded.
"""

import math

import pytest

from planetgen import tuning
from planetgen.generation import multistar
from planetgen.generation.phenomena import system_types
from planetgen.physics import accretion, ms_lifetime
from planetgen.physics import stellar_evolution as evolution
from planetgen.util import draw

SIGMAS = 4.0
OLD_AGE_YR = 8.0e9
DAYS_PER_YEAR = 365.25
MASSIVE_STARS = 30_000
RATES = system_types.SYSTEM_TYPE_RATES
KNOBS = system_types.TUNING_KNOBS


def _within(count, total, share):
    """Whether `count` of `total` is within `SIGMAS` binomial sigmas of `share`."""
    sigma = math.sqrt(total * share * (1.0 - share))
    return abs(count - total * share) <= SIGMAS * sigma


def _days(node):
    return node.period_years * DAYS_PER_YEAR


def _accreting_node(tree):
    return next(node for node in tree.nodes if node.interaction == "accreting")


@pytest.fixture(scope="module")
def massive_star_trees():
    """The trees of `MASSIVE_STARS` 20 solar mass stars three quarters of the way to their collapse (the window in which a
    compact companion exists): `(trees with a compact member, all trees)`."""
    rng = draw.Stream("gen272-massive")
    placed = multistar.Placed("star", 20.0)
    age_yr = 0.75 * ms_lifetime.lifetime_gy(20.0) * evolution.remnant_age_fraction(20.0) * 1.0e9
    trees = [multistar.draw_member_tree(placed, age_yr, rng) for _ in range(MASSIVE_STARS)]
    return [tree for tree in trees if tree.compact_count], trees, age_yr


@pytest.fixture(scope="module")
def low_mass_trees():
    """`(neutron star trees, black hole trees)` of 30,000 low-mass X-ray binary pairs each (the draw of the pair forced), on the
    old population."""
    rng = draw.Stream("gen272-low-mass")
    found = []
    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(multistar, "draw_low_mass_xray_pair", lambda placed, age_yr, unit: True)
        for kind, mass in (("neutron-star", 1.4), ("black-hole", 9.0)):
            placed = multistar.Placed(kind, mass, 20.0)
            found.append([multistar.draw_member_tree(placed, OLD_AGE_YR, rng, companion="star") for _ in range(30_000)])
    return found


# --- the draws, with the uniforms forced ---

def test_a_high_mass_compact_companion_accretes_at_the_knob_and_is_a_be_xray_binary_at_the_rows_share():
    knob = KNOBS["accreting_fraction"]["high-mass-xray-binary"]
    be_share = RATES["high-mass-xray-binary"].by_class["be_xray_share"]
    young = 1.0e7
    assert multistar.draw_accreting_compact("high-mass-xray-binary", knob + 1e-9, 0.0, young) is None
    assert multistar.draw_accreting_compact("high-mass-xray-binary", knob - 1e-9, be_share - 1e-9, young) == "disc"
    assert multistar.draw_accreting_compact("high-mass-xray-binary", knob - 1e-9, be_share + 1e-9, young) == "wind"


def test_no_old_population_system_has_a_be_xray_binary():
    assert multistar.draw_accreting_compact("high-mass-xray-binary", 0.0, 0.0, OLD_AGE_YR) == "wind"


def test_a_low_mass_compact_companion_accretes_at_the_knob_and_not_in_the_young_population():
    knob = KNOBS["accreting_fraction"]["low-mass-xray-binary"]
    assert multistar.draw_accreting_compact("low-mass-xray-binary", knob - 1e-9, age_yr=OLD_AGE_YR) == "overflow"
    assert multistar.draw_accreting_compact("low-mass-xray-binary", knob + 1e-9, age_yr=OLD_AGE_YR) is None
    assert multistar.draw_accreting_compact("low-mass-xray-binary", 0.0, age_yr=1.0e7) is None


def test_a_microquasar_is_the_stand_in_share_of_accreting_compact_objects():
    share = tuning.MICROQUASAR_SHARE_OF_XRAY_BINARIES
    assert multistar.draw_accreting_compact("microquasar", share - 1e-9) == "jets"
    assert multistar.draw_accreting_compact("microquasar", share + 1e-9) is None
    with pytest.raises(ValueError):
        multistar.draw_accreting_compact("quasar", 0.0)


def test_a_low_mass_x_ray_binary_pair_is_the_rows_rate_over_the_star_odds_of_a_neutron_star_or_a_black_hole():
    row, table = RATES["low-mass-xray-binary"], multistar.MEMBER_STATE_BY_PLACED
    neutron_star, black_hole = multistar.Placed("neutron-star", 1.4), multistar.Placed("black-hole", 9.0)
    for placed, rate, class_name in ((neutron_star, row.default, "neutron-star"), (black_hole, row.by_class["black hole"],
                                                                                     "black-hole")):
        share = rate / table[class_name]["star"][0]
        assert multistar.draw_low_mass_xray_pair(placed, OLD_AGE_YR, share - 1e-9)
        assert not multistar.draw_low_mass_xray_pair(placed, OLD_AGE_YR, share + 1e-9)
        assert not multistar.draw_low_mass_xray_pair(placed, 5.0e7, 0.0)
    recycled = multistar.Placed("neutron-star", 1.5, recycled=True)
    assert not multistar.draw_low_mass_xray_pair(recycled, OLD_AGE_YR, 0.0)
    assert not multistar.draw_low_mass_xray_pair(multistar.Placed("white-dwarf", 0.6), OLD_AGE_YR, 0.0)


def test_the_black_holes_that_accrete_from_a_star_are_the_rows_rate_times_the_knob():
    rng = draw.Stream("gen272-bh-star")
    placed = multistar.Placed("black-hole", 9.0, 25.0)
    trees = 120_000
    accreting = sum(bool(multistar.draw_member_tree(placed, OLD_AGE_YR, rng, companion="star").x_ray_luminosities())
                    for _ in range(trees))
    expected = RATES["low-mass-xray-binary"].by_class["black hole"] * KNOBS["accreting_fraction"]["low-mass-xray-binary"] / (
        multistar.MEMBER_STATE_BY_PLACED["black-hole"]["star"][0])
    # the donors that fit no period leave the pair out, so the count sits a little under the expected share
    assert 0.7 * expected * trees <= accreting <= expected * trees + SIGMAS * math.sqrt(expected * trees)


def test_the_recycled_pulsar_classes_are_the_rows_shares():
    table = multistar.MEMBER_STATE_BY_PLACED
    pulsar, row = table["recycled-pulsar"], RATES["millisecond-pulsar-companion"]
    companions = sum(default for default, _low, _high in pulsar.values())
    assert companions == pytest.approx(row.default, abs=0.001)
    for state, key in (("black-widow", "black-widow"), ("redback", "redback"), ("star", "other")):
        assert pulsar[state][0] == pytest.approx(row.default * row.by_class[key], rel=1e-3)
    assert 0.80 * row.default < pulsar["white-dwarf"][0] < 0.86 * row.default


def test_the_flag_bits_fit_their_columns_and_do_not_collide():
    bits = list(multistar.MEMBER_FLAGS.values())
    assert len(set(bits)) == len(bits) and max(bits) < 2 ** 16
    nodes = [multistar.NODE_INTERACTED, multistar.NODE_KL_ACTIVE, multistar.NODE_WIND_FED, multistar.NODE_DISC_FED]
    assert len(set(nodes)) == len(nodes) and max(nodes) < 2 ** 8


# --- the X-ray luminosity ---

def test_the_eddington_limit_is_the_textbook_value():
    assert accretion.eddington_luminosity_w(1.0) == pytest.approx(1.26e31, rel=0.01)
    assert accretion.eddington_luminosity_w(10.0) == pytest.approx(10.0 * accretion.eddington_luminosity_w(1.0))


def test_a_neutron_stars_efficiency_is_its_surface_potential_and_a_black_holes_is_the_stand_in():
    radius_au = 12.0e3 / 1.495978707e11
    assert accretion.efficiency("neutron-star", 1.4, radius_au) == pytest.approx(0.17, abs=0.01)
    assert accretion.efficiency("black-hole", 10.0, 1.0e-9) == tuning.ACCRETION_EFFICIENCY_BLACK_HOLE


def test_a_flow_far_past_the_limit_radiates_exactly_the_limit():
    power = accretion.x_ray_luminosity_w("overflow", "neutron-star", 1.4, 8.0e-8, "giant", 20.0, 0.5, 0.05, 0.0)
    assert power == accretion.eddington_luminosity_w(1.4)


def test_a_wind_fed_pair_is_brighter_the_closer_it_is_and_an_unknown_mode_is_refused():
    def power(axis_au):
        return accretion.x_ray_luminosity_w("wind", "neutron-star", 1.4, 8.0e-8, "giant", 25.0, 0.1, axis_au, 0.0)
    assert 0.0 < power(0.4) < power(0.2) <= accretion.eddington_luminosity_w(1.4)
    with pytest.raises(ValueError):
        accretion.x_ray_luminosity_w("magic", "neutron-star", 1.4, 8.0e-8, "giant", 25.0, 0.1, 0.2, 0.0)


# --- the high-mass X-ray binaries ---

def test_the_compact_companion_of_a_massive_star_is_the_rows_share_and_a_share_of_those_accrete(massive_star_trees):
    compact, trees, _age = massive_star_trees
    window_share = RATES["high-mass-xray-binary"].by_class["compact_companion"] / multistar.REMNANT_WINDOW_FRACTION
    assert _within(len(compact), len(trees), window_share)
    accreting = [tree for tree in compact if tree.x_ray_luminosities()]
    assert _within(len(accreting), len(compact), KNOBS["accreting_fraction"]["high-mass-xray-binary"])


def test_sixty_percent_of_accreting_massive_binaries_are_be_xray_binaries_with_a_be_donor(massive_star_trees):
    compact, _trees, age_yr = massive_star_trees
    accreting = [tree for tree in compact if tree.x_ray_luminosities()]
    be_xray = [tree for tree in accreting if tree.be_xray_nodes()]
    assert _within(len(be_xray), len(accreting), RATES["high-mass-xray-binary"].by_class["be_xray_share"])
    assert age_yr / 1.0e9 < tuning.YOUNG_POPULATION_MAX_AGE_GY
    for tree in be_xray:
        node = tree.nodes[tree.be_xray_nodes()[0] - 1]
        donor = next(m for m in tree.members if m.parent_node_no == node.node_no and m.kind not in multistar.COMPACT_KINDS)
        assert donor.flags & multistar.MEMBER_BE
        assert _days(node) >= tuning.BE_XRAY_PERIOD_DAYS[0] * 0.999


def test_an_accreting_pair_shows_an_x_ray_luminosity_at_or_under_the_limit_of_its_accretor(massive_star_trees):
    compact, _trees, _age = massive_star_trees
    seen = 0
    for tree in compact:
        for node_no, power in tree.x_ray_luminosities().items():
            node = tree.nodes[node_no - 1]
            accretor = next(m for m in tree.members if m.parent_node_no == node_no and m.kind in multistar.COMPACT_KINDS)
            assert node.interaction == "accreting"
            assert 0.0 < power <= accretion.eddington_luminosity_w(accretor.mass_solar)
            seen += 1
    assert seen > 50


def test_a_wind_fed_pair_has_a_supergiant_period_and_the_wind_bit(massive_star_trees):
    compact, _trees, _age = massive_star_trees
    winds = [(tree, node) for tree in compact for node in tree.nodes
             if node.interaction == "accreting" and node.flags & multistar.NODE_WIND_FED]
    assert len(winds) > 20
    for _tree, node in winds:
        assert _days(node) >= tuning.SUPERGIANT_XRAY_PERIOD_DAYS[0] * 0.999


# --- the low-mass X-ray binaries and the microquasars ---

def test_one_percent_of_low_mass_companions_of_a_neutron_star_or_a_black_hole_accrete(low_mass_trees):
    for trees in low_mass_trees:
        accreting = [tree for tree in trees if tree.x_ray_luminosities()]
        # the retries that find no lobe-filling donor leave the pair out, so the count sits a little under the knob
        knob = KNOBS["accreting_fraction"]["low-mass-xray-binary"]
        assert 0.8 * knob * len(trees) <= len(accreting) <= knob * len(trees) + SIGMAS * math.sqrt(knob * len(trees))


def test_a_low_mass_x_ray_binary_has_a_main_sequence_donor_on_a_period_of_hours(low_mass_trees):
    for trees in low_mass_trees:
        for tree in trees:
            for node_no, power in tree.x_ray_luminosities().items():
                node = tree.nodes[node_no - 1]
                donor = next(m for m in tree.members if m.parent_node_no == node_no and m.kind not in multistar.COMPACT_KINDS)
                assert donor.kind == "main-sequence" and donor.mass_solar <= tuning.LOW_MASS_XRAY_DONOR_MASS_RANGE_SOL[1]
                low, high = tuning.LOW_MASS_XRAY_PERIOD_DAYS
                assert low * 0.99 <= _days(node) <= high * 1.01
                assert not node.flags & (multistar.NODE_WIND_FED | multistar.NODE_DISC_FED)


def test_microquasars_are_the_share_of_accreting_neutron_stars_and_black_holes_and_only_those(massive_star_trees, low_mass_trees):
    accreting, flagged = 0, 0
    trees = [tree for tree in massive_star_trees[0]] + [tree for group in low_mass_trees for tree in group]
    for tree in trees:
        jets = any(m.flags & multistar.MEMBER_MICROQUASAR for m in tree.members)
        if tree.x_ray_luminosities():
            accreting += 1
            flagged += jets
        else:
            assert not jets
    assert accreting > 400
    assert _within(flagged, accreting, tuning.MICROQUASAR_SHARE_OF_XRAY_BINARIES)


def test_a_white_dwarf_accretor_has_no_x_ray_luminosity():
    rng = draw.Stream("gen272-cv")
    placed = multistar.Placed("white-dwarf", 0.6, 2.0)
    seen = 0
    for _ in range(40_000):
        tree = multistar.draw_member_tree(placed, OLD_AGE_YR, rng)
        if tree.nodes and tree.nodes[0].interaction == "accreting":
            assert tree.x_ray_luminosities() == {}
            seen += 1
    assert seen > 5


# --- the millisecond pulsars and their spiders ---

@pytest.mark.parametrize("state,low,high,label", [
    ("black-widow", *tuning.BLACK_WIDOW_MASS_RANGE_SOL, "black-widow"),
    ("redback", *tuning.REDBACK_MASS_RANGE_SOL, "redback")])
def test_a_spider_companion_has_the_classes_mass_and_a_period_of_hours_and_the_pulsar_ablates(state, low, high, label):
    rng = draw.Stream(f"gen272-{state}")
    placed = multistar.Placed("neutron-star", 1.5, 8.0, recycled=True)
    seen = 0
    for _ in range(600):
        tree = multistar.draw_member_tree(placed, OLD_AGE_YR, rng, companion=state)
        for node_no, spider in tree.spider_nodes().items():
            node = tree.nodes[node_no - 1]
            pair = [m for m in tree.members if m.parent_node_no == node_no]
            companion = next(m for m in pair if m.kind != "neutron-star")
            pulsar = next(m for m in pair if m.kind == "neutron-star")
            assert spider == label and pulsar.flags & multistar.MEMBER_ABLATING
            assert low * 0.99 <= companion.initial_mass_solar <= high * 1.01
            assert not 0.07 < companion.mass_solar < 0.1
            assert tuning.SPIDER_PERIOD_DAYS[0] * 0.99 <= _days(node) <= tuning.SPIDER_PERIOD_DAYS[1] * 1.01
            # A pair of hours cannot have held the pulsar's giant progenitor, whose largest radius is far past the whole orbit, so the history
            # rule (GEN.1012, `pair_history`) makes it post-common-envelope: the usual channel to a low-mass X-ray binary and its spider.
            assert node.interaction == "post-common-envelope"
            seen += 1
    assert seen > 500


def test_a_white_dwarf_companion_does_not_ablate_the_pulsar():
    rng = draw.Stream("gen272-msp-wd")
    placed = multistar.Placed("neutron-star", 1.5, 8.0, recycled=True)
    for _ in range(200):
        tree = multistar.draw_member_tree(placed, OLD_AGE_YR, rng, companion="white-dwarf")
        assert tree.spider_nodes() == {}
        assert not tree.members[0].flags & multistar.MEMBER_ABLATING


# --- dormant companions, compact triples and magnetars ---

def test_a_black_hole_with_a_star_can_be_dormant_like_the_gaia_systems():
    rng = draw.Stream("gen272-dormant")
    placed = multistar.Placed("black-hole", 9.0, 25.0)
    gaia = 0
    for _ in range(4000):
        tree = multistar.draw_member_tree(placed, OLD_AGE_YR, rng, companion="star")
        for node_no in tree.dormant_nodes():
            node = tree.nodes[node_no - 1]
            star = next(m for m in tree.members if m.parent_node_no == node_no and m.kind not in multistar.COMPACT_KINDS)
            gaia += 100.0 <= _days(node) <= 5000.0 and 0.5 <= star.mass_solar <= 1.5
    assert gaia > 100


@pytest.mark.parametrize("kind,mass", [("neutron-star", 1.4), ("black-hole", 9.0)])
def test_thirty_percent_of_compact_binaries_have_a_third_member_and_never_more_than_three_compact_members(kind, mass):
    rng = draw.Stream(f"gen272-triples-{kind}")
    placed = multistar.Placed(kind, mass, 20.0)
    trees = [multistar.draw_member_tree(placed, OLD_AGE_YR, rng, companion="star") for _ in range(6000)]
    pairs = [tree for tree in trees if tree.nodes]
    triples = sum(tree.member_count >= 3 for tree in pairs)
    assert _within(triples, len(pairs), RATES["compact-triple"].default)
    assert max(tree.compact_count for tree in pairs) <= multistar.MAX_COMPACT_MEMBERS


def test_magnetar_binaries_are_zero_by_default():
    assert RATES["magnetar-binary"].default == 0.0
