# tests/test_rogue_radiogenic.py

"""GEN.1134: the long-lived heat of a thrown-out planet from the star it left (`planetgen.physics.rogue_radiogenic`), identical to a bound planet's for the same star."""

import random
from types import SimpleNamespace

import pytest

from planetgen.physics import constants, interior_heat, radiogenic, rogue_radiogenic as rr, stellar_composition

KEY = stellar_composition.abundance_key(mass=1.0, distance_ly=26000.0)


def _body(fe_h, age):
    star = SimpleNamespace(metallicity_fe_h=fe_h, age=age, mass=1.0, galactic_center_dist_ly=26000.0)
    return SimpleNamespace(body_type="r", star=star, system_config=None)


@pytest.mark.parametrize("fe_h, age", [(0.0, 4.6), (-0.5, 10.0), (0.3, 1.0)])
def test_a_rogue_has_the_content_a_bound_planet_of_the_same_star_has(fe_h, age):
    body = _body(fe_h, age)
    rng_bound, rng_rogue = random.Random(9), random.Random(9)
    bound = radiogenic.draw_content(radiogenic.star_offset_dex(body.star), rng=rng_bound, th_u_factor=radiogenic.star_th_u_factor(body))
    rogue = rr.content(fe_h, age, KEY, rng=rng_rogue)
    assert rogue == bound


def test_the_host_offset_is_the_stars_and_the_same_for_every_planet_it_made():
    assert rr.host_offset_dex(0.0, 4.6, KEY) == rr.host_offset_dex(0.0, 4.6, KEY)
    assert rr.host_offset_dex(0.0, 4.6, KEY) == radiogenic.host_offset_dex(0.0, 4.6, KEY)
    assert abs(rr.host_offset_dex(0.0, 4.6, KEY)) <= 0.5
    other = stellar_composition.abundance_key(mass=0.9, distance_ly=20000.0)
    assert rr.host_offset_dex(0.0, 4.6, other) != rr.host_offset_dex(0.0, 4.6, KEY)


def test_an_older_stars_rogue_has_less_thorium_over_uranium_enrichment():
    assert rr.th_u_factor(1.0) >= 1.0
    assert rr.th_u_factor(10.0) <= rr.th_u_factor(1.0)


def test_the_spread_over_stars_is_inside_the_sourced_limit():
    offsets = [rr.host_offset_dex(fe, age, stellar_composition.abundance_key(mass=0.5 + 0.01 * i, distance_ly=24000.0 + 100.0 * i))
               for i, (fe, age) in enumerate([(f, a) for f in (-0.4, 0.0, 0.3) for a in (1.0, 5.0, 10.0)] * 4)]
    assert max(offsets) - min(offsets) > 0.05
    assert all(abs(o) <= 0.5 for o in offsets)                                     # RADIOGENIC_OFFSET_LIMIT_DEX


def test_the_rogues_heat_is_the_interior_models_for_its_content():
    content = rr.content(0.0, 4.6, KEY, rng=random.Random(1))
    rock = radiogenic.Rock(0.675, 0.325, 0.0, radiogenic.BULK_SILICATE_EARTH_AL_MASS_FRACTION)
    interior = interior_heat.draw_interior(random.Random(2))
    heat = rr.heat(content, interior, constants.EARTH_MASS_TO_KG, 6371.0, 4.5, rock)
    assert heat == interior_heat.heat_flow(content, interior, constants.EARTH_MASS_TO_KG, 6371.0, 4.5, rock)
    assert 0.03 < heat.flux_w_m2 < 0.2                                             # an Earth twin: the Earth's 0.087 W/m2 within the model's spread
    older = rr.heat(content, interior, constants.EARTH_MASS_TO_KG, 6371.0, 10.0, rock)
    assert older.flux_w_m2 < heat.flux_w_m2


def test_the_activity_cut_sits_between_mars_and_the_earth():
    assert rr.is_active(0.087) and not rr.is_active(0.02)
    assert rr.is_active(rr.ACTIVE_HEAT_FLUX_W_M2) and not rr.is_active(0.039)


def test_bad_inputs_are_refused():
    with pytest.raises(rr.RogueRadiogenicError):
        rr.host_offset_dex(0.0, 0.0, KEY)
    with pytest.raises(rr.RogueRadiogenicError):
        rr.th_u_factor(-1.0)
