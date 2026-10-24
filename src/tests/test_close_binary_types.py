# tests/test_close_binary_types.py

"""
GEN.270: the close-binary types of the member tree (`planetgen/generation/multistar.py`): contact binaries and their
tertiaries, Algols, merger products, Be stars with a stripped companion, R Coronae Borealis, FK Comae and
Thorne-Zytkow stars, and the derived eclipsing and spectroscopic flags. The audits are those of
docs/design/multi-star-systems.md section 15: 10 (type counts inside their windows) and 11 (no type where its gate
forbids it). No database or Redis: the draws are seeded.
"""

import math

import pytest

from planetgen import tuning
from planetgen.generation import multistar
from planetgen.generation.phenomena import system_types
from planetgen.util import draw

SIGMAS = 4.0


def _within(count, total, expected_share):
    """Whether `count` of `total` is within `SIGMAS` binomial sigmas of the share."""
    sigma = math.sqrt(total * expected_share * (1.0 - expected_share))
    return abs(count - total * expected_share) <= SIGMAS * sigma


def _trees(key, count, mass_range, age_gy, **kwargs):
    rng = draw.Stream(key)
    for _ in range(count):
        mass = multistar._alive_star_mass(*mass_range, age_gy * 1.0e9, rng)
        yield multistar.draw_member_tree(multistar.Placed("star", mass), age_gy * 1.0e9, rng, **kwargs)


# --- the single-star type draws, with the uniforms forced ---

def test_a_massive_young_star_is_a_merger_product_at_the_knob_share():
    share = system_types.TUNING_KNOBS["merger_product_fraction"]["massive"]
    assert multistar.draw_merger_product(15.0, 5.0e6, share - 1e-9) == "massive"
    assert multistar.draw_merger_product(15.0, 5.0e6, share + 1e-9) is None
    assert multistar.draw_merger_product(15.0, 2.0e9, 0.0) is None


def test_a_field_star_is_a_blue_straggler_at_the_knob_share_and_never_when_young():
    share = system_types.TUNING_KNOBS["merger_product_fraction"]["field"]
    assert multistar.draw_merger_product(1.0, 6.0e9, share - 1e-9) == "field"
    assert multistar.draw_merger_product(1.0, 6.0e9, share + 1e-9) is None
    assert multistar.draw_merger_product(1.0, 5.0e6, 0.0) is None


def test_a_young_b_star_is_a_be_star_at_the_row_share_with_a_stripped_companion_share_of_those():
    entry = system_types.SYSTEM_TYPE_RATES["be-star"]
    assert multistar.draw_be_star(5.0, 5.0e7, entry.default - 1e-9, 0.0) == (True, True)
    assert multistar.draw_be_star(5.0, 5.0e7, entry.default - 1e-9,
                                  entry.by_class["stripped_companion"] + 1e-9) == (True, False)
    assert multistar.draw_be_star(5.0, 5.0e7, entry.default + 1e-9, 0.0) == (False, False)
    assert multistar.draw_be_star(1.0, 5.0e7, 0.0, 0.0) == (False, False)
    assert multistar.draw_be_star(5.0, 2.0e9, 0.0, 0.0) == (False, False)


def test_only_a_young_massive_star_gains_mass():
    share = system_types.SYSTEM_TYPE_RATES["massive-merger-product"].by_class["mass_gainer"]
    assert multistar.draw_mass_gainer(15.0, 5.0e6, share - 1e-9)
    assert not multistar.draw_mass_gainer(15.0, 5.0e6, share + 1e-9)
    assert not multistar.draw_mass_gainer(1.0, 5.0e6, 0.0)
    assert not multistar.draw_mass_gainer(15.0, 2.0e9, 0.0)


def test_an_old_giant_is_an_exotic_star_at_the_counts_of_the_catalogue():
    share = tuning.OLD_GIANT_SYSTEM_SHARE
    rcb = system_types.fraction_of_systems("rcb-fk-comae") / share
    fk = system_types.SYSTEM_TYPE_RATES["rcb-fk-comae"].by_class["FK Comae"] / system_types.GALAXY_SYSTEMS / share
    tz = system_types.fraction_of_systems("thorne-zytkow") / share
    assert multistar.draw_exotic_star("giant", 8.0e9, rcb - 1e-12) == "r-crb"
    assert multistar.draw_exotic_star("giant", 8.0e9, rcb + 1e-12) == "fk-comae"
    assert multistar.draw_exotic_star("giant", 8.0e9, rcb + fk + 1e-12) == "thorne-zytkow"
    assert multistar.draw_exotic_star("giant", 8.0e9, rcb + fk + tz + 1e-12) is None
    assert multistar.draw_exotic_star("main-sequence", 8.0e9, 0.0) is None
    assert multistar.draw_exotic_star("giant", 1.0e9, 0.0) is None


# --- the roster, with a type forced ---

def _roster(monkeypatch, key, placed=None, age_yr=5.0e6, **forced):
    for name, value in forced.items():
        monkeypatch.setattr(multistar, name, value)
    placed = placed or multistar.Placed("star", 15.0)
    return multistar._roster(placed, age_yr, draw.Stream(key), multistar.MAX_MEMBERS, None, 1.0)


def test_a_massive_merger_product_has_one_star_fewer_and_a_rejuvenated_age(monkeypatch):
    plain_leaves, _ = _roster(monkeypatch, "merger", draw_merger_product=lambda *a: None)
    leaves, _ = _roster(monkeypatch, "merger", draw_merger_product=lambda *a: "massive")
    assert len(leaves) == max(1, len(plain_leaves) - 1)
    assert leaves[0].flags & multistar.MEMBER_MERGER_PRODUCT
    low, high = multistar.REJUVENATION_FRACTION_RANGE
    assert low * 5.0e6 <= leaves[0].effective_age_yr <= high * 5.0e6
    assert plain_leaves[0].effective_age_yr is None


def test_a_be_star_with_a_stripped_companion_has_one_and_is_a_mass_gainer(monkeypatch):
    leaves, units = _roster(monkeypatch, "be", Placed_B(), draw_be_star=lambda *a: (True, True))
    assert leaves[0].flags & multistar.MEMBER_BE
    assert leaves[0].flags & multistar.MEMBER_MASS_GAINER
    assert leaves[1].kind == "stripped-helium"
    low, high = tuning.BE_STRIPPED_COMPANION_MASS_RANGE_SOL
    assert low <= leaves[1].mass <= high
    assert leaves[1].initial > leaves[0].initial


def Placed_B():
    return multistar.Placed("star", 5.0)


def test_a_be_star_without_a_stripped_companion_is_flagged_alone(monkeypatch):
    leaves, _ = _roster(monkeypatch, "be", Placed_B(), draw_be_star=lambda *a: (True, False))
    assert leaves[0].flags & multistar.MEMBER_BE
    assert all(leaf.kind != "stripped-helium" for leaf in leaves)


def _giant_placed():
    for mass in (1.0 + 0.02 * step for step in range(40)):
        if multistar._star_leaf(mass, 9.0e9, draw.Stream("probe"), lighten=False).kind == "giant":
            return multistar.Placed("star", mass)
    pytest.fail("no giant among 1.0 to 1.8 solar masses at 9 Gy")


@pytest.mark.parametrize("exotic,flag", [("r-crb", multistar.MEMBER_R_CRB), ("fk-comae", multistar.MEMBER_FK_COMAE)])
def test_an_rcb_or_fk_comae_star_is_a_flagged_merger_product(monkeypatch, exotic, flag):
    monkeypatch.setattr(multistar, "draw_merger_product", lambda *a: None)
    monkeypatch.setattr(multistar, "draw_exotic_star", lambda *a: exotic)
    placed = _giant_placed()
    leaves, _ = multistar._roster(placed, 9.0e9, draw.Stream("exotic"), multistar.MAX_MEMBERS, None, 1.0)
    assert leaves[0].flags & flag
    assert leaves[0].flags & multistar.MEMBER_MERGER_PRODUCT


def test_a_thorne_zytkow_object_is_a_single_theoretical_member(monkeypatch):
    monkeypatch.setattr(multistar, "draw_exotic_star", lambda *a: "thorne-zytkow")
    tree = multistar.draw_member_tree(_giant_placed(), 9.0e9, draw.Stream("tz"))
    assert tree.member_count == 1 and not tree.nodes
    member = tree.members[0]
    assert member.kind == "thorne-zytkow"
    assert member.flags & multistar.MEMBER_THEORETICAL


# --- derived flags ---

def test_a_pair_eclipses_when_its_inclination_is_inside_the_grazing_angle():
    assert multistar.pair_eclipses(1.0, 90.0, 0.005, 0.005)
    grazing = math.degrees(math.asin(0.01))
    assert multistar.pair_eclipses(1.0, 90.0 - grazing + 0.01, 0.005, 0.005)
    assert not multistar.pair_eclipses(1.0, 90.0 - grazing - 0.01, 0.005, 0.005)
    assert not multistar.pair_eclipses(1.0, 0.0, 0.005, 0.005)
    assert multistar.pair_eclipses(0.01, 10.0, 0.006, 0.006)


def test_eclipsing_and_spectroscopic_nodes_are_derived_from_the_tree():
    seen_eclipsing = seen_spectroscopic = 0
    for tree in _trees("derived", 3000, (0.8, 1.4), 2.0):
        by_no = {node.node_no: node for node in tree.nodes}
        for node_no in tree.spectroscopic_nodes():
            assert by_no[node_no].separation_au < tuning.SPECTROSCOPIC_MAX_AXIS_AU
            seen_spectroscopic += 1
        for node_no in tree.eclipsing_nodes():
            node = by_no[node_no]
            pair = [member for member in tree.members if member.parent_node_no == node_no]
            assert len(pair) == 2
            assert multistar.pair_eclipses(node.separation_au, node.inclination_deg, pair[0].radius_au,
                                           pair[1].radius_au)
            seen_eclipsing += 1
    assert seen_spectroscopic and seen_eclipsing


# --- audit 11: no type where its gate forbids it ---

def test_no_algol_is_drawn_among_old_stars_and_no_solar_contact_binary_among_the_youngest():
    for tree in _trees("gate-old", 4000, (0.8, 1.2), 5.0):
        assert all(node.interaction != "semi-detached" for node in tree.nodes)
    for tree in _trees("gate-young", 4000, (0.8, 1.2), 0.05):
        assert all(node.interaction != "contact" for node in tree.nodes)


def test_a_massive_over_contact_binary_is_allowed_among_the_youngest():
    assert not multistar._gate_forbids("contact", multistar.MASSIVE_PRIMARY_MIN_MASS_SOL + 1.0, 1.0e7)
    assert multistar._gate_forbids("contact", 1.0, 1.0e7)
    assert not multistar._gate_forbids("contact", 1.0, 5.0e9)
    assert not multistar._gate_forbids("semi-detached", 3.0, 1.0e9)
    assert multistar._gate_forbids("semi-detached", 1.0, 5.0e9)


# --- audit 10: counts inside their windows ---

def test_the_contact_share_of_solar_type_systems_is_inside_the_catalogue_window():
    entry = system_types.SYSTEM_TYPE_RATES["contact-binary"]
    total = 40000
    contact = sum(any(node.interaction == "contact" for node in tree.nodes)
                  for tree in _trees("contact-share", total, (0.8, 1.2), 5.0))
    assert entry.low * 0.9 <= contact / total <= entry.high


def test_the_algol_share_of_disk_stars_is_inside_the_catalogue_window():
    entry = system_types.SYSTEM_TYPE_RATES["algol"]
    total = 60000
    rng = draw.Stream("algol-share")
    algol = 0
    for _ in range(total):
        age = rng.uniform(0.05, 10.0) * 1.0e9  # a disk of constant star formation
        mass = multistar._alive_star_mass(0.1, 60.0, age, rng)
        tree = multistar.draw_member_tree(multistar.Placed("star", mass), age, rng)
        algol += any(node.interaction == "semi-detached" for node in tree.nodes)
    assert entry.low * 0.7 * total <= algol <= entry.high * total


def test_an_algol_or_contact_pair_holds_no_star_that_is_not_fusing():
    for tree in _trees("fusing", 6000, (0.08, 3.0), 0.5):
        for node in tree.nodes:
            if node.interaction in ("contact", "semi-detached"):
                assert all(member.kind not in multistar.UNFUSING_KINDS for member in tree.members
                           if member.parent_node_no == node.node_no)


def test_an_algol_donor_that_was_the_heavier_at_birth_is_flagged_mass_reversed():
    reversed_count = checked = 0
    for tree in _trees("reversed", 150000, (0.8, 8.0), 1.0):
        for node in tree.nodes:
            if node.interaction != "semi-detached":
                continue
            pair = [member for member in tree.members if member.parent_node_no == node.node_no]
            checked += 1
            flagged = [member for member in pair if member.flags & multistar.MEMBER_MASS_REVERSED]
            assert len(flagged) <= 1
            if flagged:
                reversed_count += 1
                other = pair[0] if flagged[0] is pair[1] else pair[1]
                assert flagged[0].initial_mass_solar > other.initial_mass_solar
                assert other.flags & multistar.MEMBER_MASS_GAINER
    assert checked and reversed_count


def test_a_contact_binary_has_a_tertiary_at_about_the_catalogue_share():
    share = system_types.SYSTEM_TYPE_RATES["contact-binary"].by_class["tertiary"]
    contact = triples = 0
    for tree in _trees("tertiary", 120000, (0.8, 1.2), 5.0):
        if any(node.interaction == "contact" for node in tree.nodes):
            contact += 1
            triples += tree.star_count >= 3
    assert contact > 200
    assert abs(triples / contact - share) < 0.1


def test_a_forced_star_count_keeps_the_stars_of_a_contact_pair():
    builder = multistar._Builder([], 5.0e9, draw.Stream("alone"))
    builder.nodes = [{"interaction": "contact"}]
    leaves = [multistar._Leaf("main-sequence", "V", 1.0, 1.0, 0.005)] * 3
    assert multistar._contact_pair_alone(leaves, builder, 3, draw.Stream("alone")) is None
    assert multistar._contact_pair_alone(leaves[:2], builder, None, draw.Stream("alone")) is None


def test_the_be_and_stripped_companion_shares_of_young_b_stars_are_at_the_row():
    entry = system_types.SYSTEM_TYPE_RATES["be-star"]
    total = 4000
    be = stripped = 0
    for tree in _trees("be-share", total, (3.0, 12.0), 0.03):
        star = tree.members[0]
        if star.flags & multistar.MEMBER_BE:
            be += 1
            stripped += any(member.kind == "stripped-helium" for member in tree.members)
    assert _within(be, total, entry.default)
    assert _within(stripped, be, entry.by_class["stripped_companion"])
