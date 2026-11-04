# tests/test_belt_member.py

"""GEN.1067: the composition of one stored belt member (`planetgen.physics.belt_member`): the taxon gradient across a belt, the undifferentiated case, the density from the compaction taper."""

import random
from collections import Counter

import pytest

from planetgen.physics import belt_member as bm, compaction, small_body_composition as sbc, taxon_kernel

SUN_W = 3.828e26
AL26_SOLAR = 5.2e-5                       # the Solar System's, over the melting ratio of about 1.1e-5
AL26_POOR = 5.0e-6                        # under the melting ratio: no body melts


def _counts(distance_au, ratio, n=3000, mass=1e20):
    rng = random.Random(11)
    return Counter(bm.member(rng.random(), 0.5, mass, distance_au, SUN_W, ratio).taxon for _ in range(n))


def test_the_shares_are_the_kernel_at_the_members_own_distance():
    assert bm.shares(2.5, SUN_W, AL26_SOLAR) == taxon_kernel.kernel_at_distance(2.5, SUN_W, AL26_SOLAR)
    assert sum(bm.shares(3.1, SUN_W, AL26_SOLAR).values()) == pytest.approx(1.0)


def test_the_picks_follow_the_weights_and_the_order_of_the_taxa():
    weights = {"C": 0.5, "S": 0.5}
    assert bm.pick_taxon(0.0, weights) == "C" and bm.pick_taxon(0.49, weights) == "C"
    assert bm.pick_taxon(0.51, weights) == "S" and bm.pick_taxon(1.0, weights) == "S"
    with pytest.raises(bm.BeltMemberError):
        bm.pick_taxon(1.5, weights)
    with pytest.raises(bm.BeltMemberError):
        bm.pick_taxon(0.5, {})


def test_the_belt_is_drier_in_and_icier_out():
    inner, middle, outer = (_counts(a, AL26_SOLAR) for a in (1.8, 2.7, 3.6))
    dry = ("S", "A", "M", "V", "K", "L", "E")
    share = lambda c, group: sum(c[t] for t in group) / sum(c.values())
    assert share(inner, dry) > share(middle, dry) > share(outer, dry)
    assert share(outer, ("P", "D")) > share(middle, ("P", "D")) > share(inner, ("P", "D"))
    assert share(middle, ("C", "B")) > share(inner, ("C", "B"))                    # the hydrated zone sits in the middle
    assert share(inner, dry) > 0.6 and share(outer, ("P", "D")) > 0.6


def test_a_star_born_without_enough_aluminium_26_makes_no_differentiated_members():
    assert taxon_kernel.melting_ratio() > AL26_POOR
    counts = _counts(2.0, AL26_POOR)
    assert not any(counts[t] for t in taxon_kernel.IGNEOUS_SPLIT)
    assert counts["S"] > 0
    rich = _counts(2.0, AL26_SOLAR)
    assert sum(rich[t] for t in taxon_kernel.IGNEOUS_SPLIT) / sum(rich.values()) > 0.2


def test_the_differentiated_flag_is_the_melted_taxa():
    assert [bm.is_differentiated(t) for t in ("M", "A", "V", "S", "C", "P")] == [True, True, True, False, False, False]
    assert bm.member(0.0, 0.5, 1e20, 2.0, SUN_W, AL26_SOLAR).differentiated == bm.is_differentiated(bm.member(0.0, 0.5, 1e20, 2.0, SUN_W, AL26_SOLAR).taxon)


def test_the_density_closes_with_mass_through_the_compaction_taper():
    densities = [bm.member(0.4, 0.5, m, 2.2, SUN_W, AL26_SOLAR).density_kg_m3 for m in (2e19, 1e20, 5e20, 9.9e20)]
    assert densities == sorted(densities)
    top = bm.member(0.4, 0.5, 9.9e20, 2.2, SUN_W, AL26_SOLAR)
    assert top.filling > 0.97 and top.density_kg_m3 < top.makeup.grain_density_kg_m3 * 1.0001


def test_a_stored_member_takes_the_taper_and_a_light_one_the_family_porosity():
    stored = bm.member(0.4, 0.0, 1e19, 2.2, SUN_W, AL26_SOLAR)
    stored_other_u = bm.member(0.4, 1.0, 1e19, 2.2, SUN_W, AL26_SOLAR)
    assert stored.density_kg_m3 == stored_other_u.density_kg_m3                    # the porosity number is not read above the taper's start
    light_low = bm.member(0.4, 0.0, 1e17, 2.2, SUN_W, AL26_SOLAR)
    light_high = bm.member(0.4, 1.0, 1e17, 2.2, SUN_W, AL26_SOLAR)
    low, high = sbc.POROSITY_RANGES[sbc.TAXON_DATA[light_low.taxon].family]
    assert light_low.filling == pytest.approx(1.0 - low) and light_high.filling == pytest.approx(1.0 - high)
    assert light_low.density_kg_m3 > light_high.density_kg_m3


def test_the_taxon_sets_the_grain_density_so_metal_is_denser_than_rock():
    metal = bm.member(0.0, 0.5, 1e20, 2.0, SUN_W, AL26_SOLAR)
    assert metal.taxon in taxon_kernel.TAXA
    by_taxon = {}
    rng = random.Random(5)
    for _ in range(400):
        m = bm.member(rng.random(), 0.5, 1e20, 2.2, SUN_W, AL26_SOLAR)
        by_taxon[m.taxon] = m
    assert by_taxon["M"].density_kg_m3 > by_taxon["S"].density_kg_m3 > by_taxon["C"].density_kg_m3


def test_the_radius_is_that_of_the_sphere_of_the_mass_and_density():
    m = bm.member(0.3, 0.5, 4e20, 2.7, SUN_W, AL26_SOLAR)
    assert 4.0 / 3.0 * 3.14159265 * (m.radius_km * 1e3) ** 3 * m.density_kg_m3 == pytest.approx(4e20, rel=1e-4)
    assert 300.0 < m.radius_km < 500.0


def test_a_member_over_the_taper_is_for_the_structure_solver():
    with pytest.raises(compaction.OutOfRange):
        bm.member(0.3, 0.5, 5e21, 2.7, SUN_W, AL26_SOLAR)


def test_bad_inputs_are_refused():
    with pytest.raises(bm.BeltMemberError):
        bm.member(0.3, 0.5, 0.0, 2.7, SUN_W, AL26_SOLAR)
    with pytest.raises(bm.BeltMemberError):
        bm.member(0.3, 1.5, 1e17, 2.7, SUN_W, AL26_SOLAR)


def test_a_typical_body_has_the_family_porosity_middle_below_the_taper_and_the_tapers_filling_above():
    light = bm.typical_density_kg_m3("C", 1.0e15)
    grain = bm.sbc.makeup_from_components("belt", bm.sbc.TAXON_DATA["C"].fractions, 0.0).grain_density_kg_m3
    low, high = bm.sbc.POROSITY_RANGES[bm.sbc.TAXON_DATA["C"].family]
    assert light == pytest.approx(grain * (1.0 - 0.5 * (low + high)))
    assert bm.typical_density_kg_m3("C", 1.0e20) > light
    with pytest.raises(bm.BeltMemberError):
        bm.typical_density_kg_m3("Q", 1.0e15)
    with pytest.raises(bm.compaction.OutOfRange):
        bm.typical_density_kg_m3("C", 1.0e23)
