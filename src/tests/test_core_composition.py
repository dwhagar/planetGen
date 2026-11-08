"""
The make-up of a metal core (GEN.1364, `planetgen.physics.core_composition`): the split of the light elements by the redox of the rock and the
pressure at the base of the mantle, the liquidus the light elements move and the compositional energy the freezing inner core releases.
"""

import pytest

from planetgen.physics import constants, core_composition as cc, core_heat, interior_heat as ih, radiogenic
from planetgen.util import draw

EARTH_MASS = constants.EARTH_MASS_TO_KG
MID = ih.Interior(1.0e12, 250.0)


def test_the_redox_of_earths_and_mars_rock_is_where_their_cores_formed():
    assert -2.8 < cc.EARTH_DELTA_IW < -1.8        # Earth's core formed near IW-2
    assert -2.0 < cc.MARS_DELTA_IW < -1.0         # Mars' near IW-1 to IW-2
    assert cc.MARS_DELTA_IW > cc.EARTH_DELTA_IW
    assert cc.delta_iw(0.001) < -4.5 and cc.delta_iw(0.0) == -8.0     # a rock with next to no FeO is far below the buffer


def test_the_anchor_shares_give_back_the_cores_they_were_made_from():
    earth = cc.mean_shares(cc.EARTH_DELTA_IW, 136.0)
    for name in cc.ELEMENTS:
        assert earth[name] == pytest.approx(getattr(cc.EARTH, name) / cc.EARTH.light)
    mars = cc.mean_shares(cc.MARS_DELTA_IW, 20.0)
    assert mars["si"] == 0.0 and mars["c"] == 0.0
    assert mars["s"] == pytest.approx(6.6 / 12.7, abs=0.01) and mars["o"] == pytest.approx(5.2 / 12.7, abs=0.01)


def test_a_reduced_body_puts_silicon_in_its_core_and_an_oxidized_one_sulfur_and_oxygen():
    reduced, oxidized = cc.mean_shares(-4.0, 60.0), cc.mean_shares(-1.6, 60.0)
    assert reduced["si"] > 0.5 > oxidized["si"]
    assert oxidized["s"] + oxidized["o"] > reduced["s"] + reduced["o"]
    assert all(sum(cc.mean_shares(redox, pressure).values()) == pytest.approx(1.0) for redox in (-6.0, -3.0, -1.0) for pressure in (5.0, 80.0, 400.0))


def test_pressure_pushes_oxygen_and_silicon_into_the_metal():
    low, high = cc.mean_shares(-3.0, 30.0), cc.mean_shares(-3.0, 300.0)
    assert high["o"] > low["o"] and high["si"] > low["si"]
    assert high["s"] < low["s"]


def test_the_split_adds_up_to_the_light_total_is_the_same_for_the_same_body_and_moves_no_other_draw():
    with draw.bound("core-split"):
        before = draw.random()
    with draw.bound("core-split"):
        parts = cc.split(0.09, -2.5, 100.0, "body-a")
        after = draw.random()
    assert sum(parts.values()) == pytest.approx(0.09)
    assert before == after
    assert parts == cc.split(0.09, -2.5, 100.0, "body-a")
    assert parts != cc.split(0.09, -2.5, 100.0, "body-b")
    assert set(parts) == set(cc.ELEMENTS) and all(value >= 0.0 for value in parts.values())


def test_the_scatter_of_the_split_spans_the_models_of_earths_core():
    shares = [cc.split(1.0, cc.EARTH_DELTA_IW, 136.0, f"body-{n}") for n in range(400)]
    silicon = sorted(part["si"] for part in shares)
    assert 0.15 < silicon[20] < 0.35 and 0.6 < silicon[380] < 0.85         # the 15 models run Si from 1.9 to 7.4 wt percent of a 10 percent light total
    assert sum(part["si"] for part in shares) / 400 == pytest.approx(cc.EARTH.si / cc.EARTH.light, abs=0.05)


def test_earths_mix_is_the_reference_and_leaves_the_liquidus_and_the_energy_as_they_were():
    assert cc.liquidus_factor(cc.EARTH) == pytest.approx(1.0)
    assert cc.composition_energy_factor(cc.EARTH) == pytest.approx(1.0)
    assert cc.EARTH.fe + cc.EARTH.ni + cc.EARTH.light == pytest.approx(1.0)
    assert cc.EARTH.fe == pytest.approx(82.8 / 98.157, rel=0.001)


def test_a_core_with_more_sulfur_freezes_later_and_one_with_less_light_elements_sooner():
    sulfur = cc.Mix(ni=0.06, s=0.13, si=0.0, o=0.0, c=0.0, h=0.0)
    clean = cc.Mix(ni=0.06, s=0.0, si=0.0, o=0.0, c=0.0, h=0.0)
    assert cc.liquidus_factor(sulfur) < 1.0 < cc.liquidus_factor(clean)
    assert cc.liquidus_factor(sulfur) > cc.LIQUIDUS_FACTOR_LIMITS[0]


def test_the_fe_fes_eutectic_depression_is_within_a_third_of_the_one_at_one_bar():
    # 31 wt percent S against pure Fe: the eutectic is at 1,261 K against 1,811 K, a fall of 550 K; the model scales its slope with the temperature
    sulfur_share = 0.31
    mix = cc.Mix(ni=0.0, s=sulfur_share, si=0.0, o=0.0, c=0.0, h=0.0)
    slope_at_bar = cc.LIQUIDUS_DEPRESSION_K * core_heat.LIQUIDUS_ZERO_K / cc.ICB_TEMPERATURE_K
    fall = slope_at_bar * mix.mole_fraction_left
    assert fall == pytest.approx(550.0, rel=0.35)


def test_the_light_elements_the_solid_leaves_in_the_liquid_drive_the_compositional_energy():
    sulfur = cc.Mix(ni=0.06, s=0.13, si=0.0, o=0.0, c=0.0, h=0.0)
    silicon = cc.Mix(ni=0.06, s=0.0, si=0.13, o=0.0, c=0.0, h=0.0)
    assert cc.composition_energy_factor(sulfur) > 1.0 > cc.composition_energy_factor(silicon)
    assert cc.composition_energy_factor(cc.Mix(ni=0.06, s=0.0, si=0.0, o=0.0, c=0.0, h=0.0)) == 0.0


def _stepper(mix, mass_earth=1.0, radius_km=6371.0, metal=0.325):
    dry_m = radius_km * 1000.0
    return core_heat.Stepper(MID, mass_earth * EARTH_MASS * metal, mass_earth * EARTH_MASS * (1.0 - metal), dry_m * metal ** 0.5, dry_m, mix)


def test_the_stepper_with_earths_mix_is_the_stepper_it_was():
    earth = _stepper(cc.EARTH)
    assert earth.liquidus_zero_k == pytest.approx(core_heat.LIQUIDUS_ZERO_K)
    assert earth.liquidus_k(core_heat.LIQUIDUS_ANCHOR[0]) == pytest.approx(core_heat.LIQUIDUS_ANCHOR[1])
    assert earth.gravitational_j_kg == pytest.approx(core_heat.GRAVITATIONAL_ENERGY_J_KG)
    assert _stepper(cc.MARS).gravitational_j_kg == pytest.approx(core_heat.GRAVITATIONAL_ENERGY_J_KG * cc.composition_energy_factor(cc.MARS))


def test_a_sulfur_rich_core_lowers_the_liquidus_and_the_inner_core_is_smaller_at_the_same_temperature():
    earth, sulfur = _stepper(cc.EARTH), _stepper(cc.MARS)
    assert sulfur.liquidus_k(200.0) < earth.liquidus_k(200.0)
    assert sulfur.inner_radius_m(4500.0) <= earth.inner_radius_m(4500.0)


def test_the_pressure_at_the_base_of_the_mantle_agrees_with_the_stepper_and_with_a_body_with_no_metal():
    earth = _stepper(cc.EARTH)
    rock = radiogenic.Rock(0.675, 0.325, 0.0, radiogenic.BULK_SILICATE_EARTH_AL_MASS_FRACTION)
    assert ih.cmb_pressure_gpa(EARTH_MASS, 6371.0, rock) == pytest.approx(earth.cmb_gpa)
    assert ih.cmb_pressure_gpa(EARTH_MASS, 6371.0, radiogenic.Rock(1.0, 0.0, 0.0, 0.0)) == 0.0


def test_a_body_that_stores_no_split_takes_earths_mix_and_one_that_stores_it_takes_its_own():
    class Body:
        pass
    assert cc.mix_of(Body()) is cc.EARTH
    body = Body()
    body.core_ni_fraction = 0.07
    for field in cc.BODY_FIELDS:
        setattr(body, field, 0.01)
    mix = cc.mix_of(body)
    assert mix.ni == 0.07 and mix.s == 0.01 and mix.fe == pytest.approx(1.0 - 0.07 - 0.05)
