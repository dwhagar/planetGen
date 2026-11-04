# tests/test_captured_composition.py

"""GEN.1080: what a captured body is made of and how a family differs (`planetgen.physics.captured_composition`)."""

import pytest

from planetgen.physics import belt_member, captured_composition as cc, compaction, constants, irregular_zone, moon_composition
from planetgen.physics import small_body_composition as sbc

LUMINOSITY = 3.828e26
AL26 = 5.0e-5
SUN = constants.SOLAR_MASS_TO_KG
JUPITER = constants.JUPITER_MASS_TO_KG


def _taxa(distance_au, count=60):
    return {belt_member.pick_taxon((k + 0.5) / count, belt_member.shares(distance_au, LUMINOSITY, AL26)) for k in range(count)}


def test_a_captured_body_takes_the_kernels_taxon_at_its_source_distance():
    near = cc.compose("encounter", 1e17, 2.5, LUMINOSITY, AL26, 0.05, 0.5, 0.5)
    far = cc.compose("encounter", 1e17, 25.0, LUMINOSITY, AL26, 0.05, 0.5, 0.5)
    assert near.member.taxon == belt_member.member(0.05, 0.5, 1e17, 2.5, LUMINOSITY, AL26).taxon
    assert far.taxon in ("D", "P") and near.taxon not in ("D",)
    assert far.makeup.water > 0.0 and near.makeup.water < far.makeup.water
    assert _taxa(2.5) >= {"C", "S"} and _taxa(25.0) <= {"D", "P", "C", "B", "S", "E", "L", "K"}


def test_capture_removes_the_ice_share_the_moon_rule_says_and_nothing_else():
    low, high = moon_composition.CAPTURE_ICE_LOSS
    assert cc.ice_loss_share(0.0) == low and cc.ice_loss_share(1.0) == high
    kept = cc.compose("encounter", 1e17, 25.0, LUMINOSITY, AL26, 0.1, 0.5, 0.0)
    lost = cc.compose("encounter", 1e17, 25.0, LUMINOSITY, AL26, 0.1, 0.5, 1.0)
    assert kept.ice_loss == 0.0 and kept.makeup.water == pytest.approx(kept.member.makeup.water)
    assert lost.ice_loss == high
    assert lost.makeup.water < kept.makeup.water
    water = kept.member.makeup.water
    assert lost.makeup.water == pytest.approx(water * (1.0 - high) / (1.0 - water * high))
    assert lost.makeup.rock / lost.makeup.water > kept.makeup.rock / kept.makeup.water
    assert lost.makeup.metal + lost.makeup.rock + lost.makeup.water + lost.makeup.ice == pytest.approx(1.0)


def test_a_dry_body_loses_no_ice_it_did_not_have():
    body = cc.compose("gas", 1e17, 1.0, LUMINOSITY, AL26, 0.99, 0.5, 1.0)         # a rocky taxon in the dry zone
    assert body.member.makeup.water == 0.0 and body.makeup.water == 0.0


def test_density_and_radius_follow_the_loss_and_the_filling():
    body = cc.compose("encounter", 1e17, 25.0, LUMINOSITY, AL26, 0.1, 0.5, 1.0)
    assert body.density_kg_m3 == pytest.approx(body.makeup.grain_density_kg_m3 * body.member.filling)
    assert body.radius_km == pytest.approx((3.0 * 1e17 / (4.0 * 3.141592653589793 * body.density_kg_m3)) ** (1.0 / 3.0) / 1.0e3)
    heavier = cc.compose("encounter", 1e20, 25.0, LUMINOSITY, AL26, 0.1, 0.5, 1.0)
    assert heavier.density_kg_m3 >= body.density_kg_m3 or heavier.member.filling >= body.member.filling      # the taper closes pores for a heavy body
    with pytest.raises(compaction.OutOfRange):
        cc.compose("exchange", 1e23, 25.0, LUMINOSITY, AL26, 0.1, 0.5, 0.5)


def test_the_window_names_where_a_body_came_from():
    gas = cc.source_annulus("gas", 5.2, JUPITER, SUN, 20.0, 35.0)
    hill = irregular_zone.hill_radius_au(5.2, JUPITER, SUN)
    assert gas[0] == pytest.approx(5.2 - cc.FEEDING_ZONE_HILL_RADII * hill) and gas[1] == pytest.approx(5.2 + cc.FEEDING_ZONE_HILL_RADII * hill)
    assert cc.source_annulus("encounter", 5.2, JUPITER, SUN, 20.0, 35.0) == (20.0, 35.0)
    assert cc.source_annulus("exchange", 5.2, JUPITER, SUN, 20.0, 35.0) == (20.0, 35.0)
    # a planet so close its feeding zone would reach the star is held at a tenth of its distance
    assert cc.source_annulus("gas", 0.05, JUPITER, SUN, 20.0, 35.0)[0] > 0.0
    with pytest.raises(cc.CapturedCompositionError):
        cc.source_annulus("temporary", 5.2, JUPITER, SUN, 20.0, 35.0)
    with pytest.raises(cc.CapturedCompositionError):
        cc.source_annulus("encounter", 5.2, JUPITER, SUN, 35.0, 20.0)


def test_the_source_distance_follows_the_mass_of_the_solids():
    assert cc.draw_source_distance(20.0, 35.0, 0.0) == pytest.approx(20.0)
    assert cc.draw_source_distance(20.0, 35.0, 1.0) == pytest.approx(35.0)
    median = cc.draw_source_distance(20.0, 35.0, 0.5)
    assert median < 0.5 * (20.0 + 35.0)                             # the solids go as a^-1.5: the mass up to a goes as a^0.5, so more of it is near the inner edge
    k = 0.5
    assert (median ** k - 20.0 ** k) / (35.0 ** k - 20.0 ** k) == pytest.approx(0.5)
    with pytest.raises(cc.CapturedCompositionError):
        cc.draw_source_distance(20.0, 35.0, 1.5)


def test_a_family_shares_the_parents_capture_and_has_the_parents_taxa():
    us_d = [0.05, 0.35, 0.65, 0.95]
    us_t = [0.05, 0.35, 0.65, 0.95]
    parent = cc.draw_parent("gas", 2.0, 3.4, us_d, us_t, 0.4)
    assert len(parent.lumps) == cc.PARENT_LUMPS
    assert [l.distance_au for l in parent.lumps] == sorted(l.distance_au for l in parent.lumps)
    lump_taxa = {belt_member.member(l.u_taxon, 0.5, 1e15, l.distance_au, LUMINOSITY, AL26).taxon for l in parent.lumps}
    members = [cc.fragment(parent, u, 1e15, LUMINOSITY, AL26, 0.5) for u in (0.1, 0.3, 0.6, 0.9)]
    assert {m.taxon for m in members} == lump_taxa
    assert len(lump_taxa) >= 2                                      # the members differ
    assert {m.ice_loss for m in members} == {cc.ice_loss_share(0.4)}  # one capture, one loss
    assert all(m.window == "gas" for m in members)
    # a fragment is exactly the lump's composition: the same body composed alone
    lump = cc.fragment_lump(parent, 0.6)
    alone = cc.compose("gas", 1e15, lump.distance_au, LUMINOSITY, AL26, lump.u_taxon, 0.5, 0.4)
    assert members[2] == alone


def test_a_family_from_a_one_sided_kernel_is_nearly_uniform():
    parent = cc.draw_parent("encounter", 20.0, 35.0, [0.1, 0.4, 0.6, 0.9], [0.1, 0.3, 0.5, 0.7], 0.2)
    taxa = {cc.fragment(parent, u, 1e15, LUMINOSITY, AL26, 0.5).taxon for u in (0.1, 0.3, 0.6, 0.9)}
    assert taxa <= {"D", "P"}


def test_bad_inputs_are_refused():
    with pytest.raises(cc.CapturedCompositionError):
        cc.draw_parent("gas", 2.0, 3.0, [0.1], [0.1], 0.5)
    with pytest.raises(cc.CapturedCompositionError):
        cc.compose("temporary", 1e17, 25.0, LUMINOSITY, AL26, 0.1, 0.5, 0.5)
    with pytest.raises(cc.CapturedCompositionError):
        cc.ice_loss_share(-0.1)
    with pytest.raises(cc.CapturedCompositionError):
        cc.fragment_lump(cc.draw_parent("gas", 2.0, 3.0, [0.1] * 4, [0.1] * 4, 0.5), 1.5)
