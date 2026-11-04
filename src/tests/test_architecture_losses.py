# tests/test_architecture_losses.py

"""
GEN.1119 (RP1): every loss of the planet architecture is classed ejected, collided or bound by the Safronov number of the encounter
(`scatter.loss_fates`, `planet_architecture.losses`), so no body is simply deleted and the ejected ones are a record the rogue builder reads
(docs/design/research/rogue-planet-review.md).
"""

import math
import random

import pytest

from planetgen.physics import constants, planet_architecture as pa, scatter
from planetgen.util import draw

SUN = constants.SOLAR_MASS_TO_KG
EARTH = constants.EARTH_MASS_TO_KG
JUPITER_R = constants.JUPITER_RADIUS_KM / constants.EARTH_RADIUS_KM


def seed(distance, mass_earth, core_earth, runaway, migrated=False, formed=None):
    return pa.Seed(distance, mass_earth * EARTH, core_earth * EARTH, runaway, 1, migrated, distance if formed is None else formed)


def lost_of(**kwargs):
    base = dict(chain_mass_kg=0.0, merge_debris_kg=0.0, merges=0, capped_kg=0.0, capped=0, small=(), swept=(), thrown=(), snow_line_au=2.7, star_mass_kg=SUN)
    base.update(kwargs)
    return pa.Lost(**base)


def architecture_of(planets, lost):
    return pa.Architecture(tuple(planets), 1.0, 12.0, len(planets), lost)


# --- the rule -----------------------------------------------------------------------------------------------------------------


def test_the_shares_of_a_scattered_body_sum_to_one_and_ejection_rises_with_the_safronov_number():
    last = -1.0
    for theta in (0.0, 0.05, 0.3, 1.0, 3.0, 10.0, 100.0):
        ejected, collided, bound = scatter.loss_fates(theta)
        assert ejected + collided + bound == pytest.approx(1.0)
        assert ejected > last
        last = ejected
    assert scatter.loss_fates(0.0)[0] == 0.0


def test_a_giant_throws_out_what_it_scatters_and_a_rocky_planet_does_not():
    jupiter = scatter.safronov_number(318.0, JUPITER_R, 5.2, 1.0)
    earth = scatter.safronov_number(1.0, 1.0, 1.0, 1.0)
    assert jupiter == pytest.approx(10.0, rel=0.3)         # Jupiter's Safronov number is near 10 [R3867]
    assert earth == pytest.approx(0.07, rel=0.3)
    assert scatter.loss_fates(jupiter)[0] > 0.95
    assert scatter.loss_fates(earth)[0] < 0.01


def test_the_reach_of_a_giant_is_its_chaotic_zone_and_the_whole_path_it_moved_along():
    still = scatter.reach_band_au(5.0, 5.0, 318.0, 1.0)
    moved = scatter.reach_band_au(1.0, 5.0, 318.0, 1.0)
    width = scatter.zone_half_width_au(5.0, 318.0, 1.0)
    assert still == pytest.approx((5.0 - width, 5.0 + width))
    assert moved[1] == pytest.approx(5.0 + scatter.zone_half_width_au(1.0, 318.0, 1.0))
    assert moved[0] < 1.0 < 5.0 < moved[1] + 1e-9 + width


def test_the_hill_radius_of_a_body_is_the_textbook_one():
    assert scatter.hill_radius_au(1.0, 1.0, 1.0) == pytest.approx(0.01, rel=0.02)         # Earth's Hill sphere is 0.01 AU


# --- the ledger of one architecture -------------------------------------------------------------------------------------------


def test_a_planet_a_crowded_giant_throws_out_is_an_ejection_record_with_its_encounter():
    heavy = seed(1.57, 300.0, 30.0, True)
    light = seed(1.59, 80.0, 20.0, True)
    arch = architecture_of([heavy], lost_of(thrown=((light, heavy),)))
    (loss,) = pa.losses(arch)
    assert loss.kind == "thrown"
    assert loss.mass_kg == pytest.approx(80.0 * EARTH) and loss.core_kg == pytest.approx(20.0 * EARTH)
    assert loss.perturber_au == 1.57 and loss.perturber_mass_kg == pytest.approx(300.0 * EARTH)
    assert loss.theta == pytest.approx(scatter.safronov_number(300.0, JUPITER_R, 1.59, 1.0))
    assert loss.ejected > 0.7 and loss.ejected + loss.collided + loss.bound == pytest.approx(1.0)
    assert loss.hill_radius_au == pytest.approx(scatter.hill_radius_au(80.0, 1.59, 1.0))
    assert loss.time_yr == scatter.DISC_LIFETIME_YR and loss.ice_rich is False


def test_a_planet_a_hot_jupiter_sweeps_up_collides_instead_of_leaving_because_the_giant_is_close_to_the_star():
    hot = seed(0.05, 300.0, 30.0, True, migrated=True, formed=4.0)
    inside = seed(0.03, 2.0, 2.0, False)
    arch = architecture_of([hot], lost_of(swept=((inside, hot),)))
    (loss,) = pa.losses(arch)
    assert loss.kind == "swept" and loss.theta < 0.1
    assert loss.ejected < 0.01 and loss.collided > 0.3


def test_chain_bodies_a_giant_reaches_are_classed_and_the_rest_stay_bound_as_one_record():
    giant = seed(5.0, 300.0, 30.0, True)
    width = scatter.zone_half_width_au(5.0, 300.0, 1.0)
    reached = (5.0 + 0.5 * width, 0.05 * EARTH)
    far = [(40.0 + i, 0.02 * EARTH) for i in range(7)]
    arch = architecture_of([giant], lost_of(small=tuple([reached] + far)))
    result = pa.losses(arch)
    assert [(l.kind, l.count) for l in result] == [("sub_minimum", 1), ("sub_minimum", 7)]
    assert result[0].ejected > 0.9 and result[0].perturber_au == 5.0
    assert (result[1].ejected, result[1].collided, result[1].bound) == (0.0, 0.0, 1.0)
    assert result[1].mass_kg == pytest.approx(7 * 0.02 * EARTH)


def test_the_merges_are_a_collision_that_keeps_the_lost_mass():
    arch = architecture_of([], lost_of(merge_debris_kg=3.0 * EARTH, merges=5))
    (loss,) = pa.losses(arch)
    assert (loss.kind, loss.count, loss.collided, loss.ejected) == ("merge", 5, 1.0, 0.0)
    assert loss.mass_kg == pytest.approx(3.0 * EARTH)


def test_the_mass_a_planets_core_cap_cuts_off_is_a_collision_that_keeps_the_mass():
    arch = architecture_of([], lost_of(capped_kg=2.0 * EARTH, capped=3))
    (loss,) = pa.losses(arch)
    assert (loss.kind, loss.count, loss.collided, loss.ejected) == ("capped", 3, 1.0, 0.0)


def test_the_ejection_yield_counts_expected_bodies_above_a_mass_floor():
    heavy = seed(1.57, 300.0, 30.0, True)
    light = seed(1.59, 80.0, 20.0, True)
    arch = architecture_of([heavy], lost_of(thrown=((light, heavy),), small=((1.6, 0.05 * EARTH),)))
    count_all, mass_all = pa.ejection_yield(arch)
    count_planets, mass_planets = pa.ejection_yield(arch, 0.1)
    assert count_all > count_planets > 0.7
    assert mass_planets == pytest.approx(80.0 * pa.losses(arch)[0].ejected)
    assert mass_all > mass_planets


# --- built architectures ------------------------------------------------------------------------------------------------------


def built(star_mass, seeds_wanted, fe_h=0.0):
    out = []
    for i in range(seeds_wanted):
        with draw.bound(5000 + i):
            out.append(pa.build(star_mass * SUN, 2.7 * star_mass ** 2, 10.0 ** fe_h, fe_h, 100.0))
    return out


def test_the_planets_and_the_losses_add_up_to_the_whole_chain():
    for arch in built(1.0, 40) + built(0.4, 40, -0.2):
        total = sum(p.core_kg for p in arch.planets) + sum(l.core_kg for l in pa.losses(arch))
        assert total == pytest.approx(arch.lost.chain_mass_kg, rel=1e-9)
        for loss in pa.losses(arch):
            assert loss.ejected + loss.collided + loss.bound == pytest.approx(1.0)


def test_a_star_without_a_giant_ejects_next_to_nothing_and_one_with_giants_ejects_giants():
    stars = built(1.0, 150, 0.2)
    giant_free = [a for a in stars if not any(p.runaway for p in a.planets) and not a.lost.thrown and not a.lost.swept]
    rich = [a for a in stars if a.lost.thrown or a.lost.swept]
    assert giant_free and rich
    free_mass = sum(pa.ejection_yield(a)[1] for a in giant_free) / len(giant_free)
    assert free_mass < 1e-3
    rich_mass = sum(pa.ejection_yield(a)[1] for a in rich) / len(rich)
    assert rich_mass > 10.0 * max(free_mass, 1e-6)


def test_the_yield_per_star_is_reported_and_is_far_below_the_flat_rate():
    stars = built(1.0, 300)
    planets = sum(pa.ejection_yield(a, 0.1)[0] for a in stars) / len(stars)
    mass = sum(pa.ejection_yield(a, 0.1)[1] for a in stars) / len(stars)
    # The architecture alone throws out well under one planet per star (the review report measures about 0.2): the rest of the free-floating
    # population is not made by it (docs/design/research/rogue-planet-review.md).
    assert 0.0 < planets < 1.0
    assert mass > 0.0
