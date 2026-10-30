# tests/test_field_composition.py

"""What a collision field is made of (GEN.1102, `physics/field_composition.py`): the outermost layers leave first, the cores stay in the remnant, ice only from water-rich parents, and the
mass is accounted for."""

import math

import pytest

from planetgen.physics import constants, field_composition as fc, small_body_composition

EARTH = constants.EARTH_MASS_TO_KG
DRY = fc.Parent(EARTH, 1.0 / 3.0, 0.0)
WET = fc.Parent(EARTH, 1.0 / 3.0, 0.05)


def test_the_layers_of_a_parent_add_up_to_its_mass():
    layers = fc.layer_masses(WET)
    assert layers.ice == pytest.approx(0.05 * EARTH) and layers.core == pytest.approx(EARTH / 3.0)
    assert layers.ice + layers.mantle + layers.core == pytest.approx(EARTH)
    for bad in (fc.Parent(0.0, 0.3, 0.0), fc.Parent(EARTH, 0.7, 0.4), fc.Parent(EARTH, -0.1, 0.0)):
        with pytest.raises(fc.FieldCompositionError):
            fc.layer_masses(bad)


def test_a_parent_gives_up_its_surface_first():
    lost = fc.escaping_layers(WET, 0.04)
    assert lost.ice == pytest.approx(0.04 * EARTH) and lost.mantle == 0.0 and lost.core == 0.0     # all ice while the shell lasts
    lost = fc.escaping_layers(WET, 0.13)
    assert lost.ice == pytest.approx(0.05 * EARTH) and lost.mantle == pytest.approx(0.08 * EARTH) and lost.core == 0.0
    lost = fc.escaping_layers(DRY, 0.8)                                                          # past the mantle the core goes
    assert lost.mantle == pytest.approx(EARTH * 2 / 3) and lost.core == pytest.approx(EARTH * (0.8 - 2 / 3))
    assert sum(fc.escaping_layers(DRY, 1.0)) == pytest.approx(EARTH)
    assert fc.escaping_layers(DRY, 0.0) == (0.0, 0.0, 0.0)
    with pytest.raises(fc.FieldCompositionError):
        fc.escaping_layers(DRY, 1.5)


def test_a_core_leak_takes_iron_before_the_mantle_is_gone():
    lost = fc.escaping_layers(DRY, 0.1, core_leak_share=0.2)
    assert lost.core == pytest.approx(0.02 * EARTH) and lost.mantle == pytest.approx(0.08 * EARTH)


def test_the_field_of_a_dry_pair_is_silicate_with_no_iron():
    layers = fc.field_layers(DRY, DRY, 0.13)
    assert layers.core == 0.0 and layers.ice == 0.0
    assert layers.mantle == pytest.approx(0.26 * EARTH)
    components = fc.field_components(layers)
    assert components == pytest.approx({"forsterite": 0.5, "enstatite": 0.5})


def test_a_water_rich_pair_makes_a_field_of_ice_shards():
    layers = fc.field_layers(WET, WET, 0.13)
    assert layers.ice / sum(layers) == pytest.approx(0.05 / 0.13)                              # the whole shell leaves first
    components = fc.field_components(layers)
    assert sum(components.values()) == pytest.approx(1.0) and components["water_ice"] == pytest.approx(0.05 / 0.13)
    assert "kamacite" not in components


def test_the_remnant_keeps_the_cores_by_the_average_of_two_models():
    big = fc.remnant_core_fraction(DRY, DRY, 0.235)                  # a remnant of 0.47 Earth can hold 0.67 Earth of core only in part
    assert big == pytest.approx(0.5 * (1.0 + 1.0 / 3.0))
    whole = fc.remnant_core_fraction(DRY, DRY, 1.0)
    assert whole == pytest.approx(1.0 / 3.0)                         # all the mass: the pair's own core share
    assert 1.0 / 3.0 <= fc.remnant_core_fraction(DRY, DRY, 0.1) <= 1.0
    with pytest.raises(fc.FieldCompositionError):
        fc.remnant_core_fraction(DRY, DRY, 0.0)


def test_two_earths_make_a_field_only_in_the_rubble_window():
    slow = fc.collision_field(DRY, DRY, 10e3)
    assert slow.regime in ("merger", "partial_erosion") and slow.field is None and slow.field_mass_kg == 0.0
    assert slow.remnant.mass_kg + slow.bound_debris_kg == pytest.approx(2 * EARTH)
    window = fc.collision_field(DRY, DRY, 25e3)
    assert window.regime == "catastrophic" and window.field is not None
    assert 0.0 < window.field_mass_kg < 0.2 * EARTH
    fast = fc.collision_field(DRY, DRY, 60e3)
    assert fast.regime == "super_catastrophic" and fast.field_mass_kg / (2 * EARTH) == pytest.approx(0.133, abs=0.015)     # 11 to 13.6 percent of the mass
    assert fast.remnant.mass_kg < 0.05 * EARTH


def test_a_collision_accounts_for_all_the_mass():
    for speed in (25e3, 28e3, 35e3, 60e3):
        o = fc.collision_field(WET, DRY, speed)
        assert o.remnant.mass_kg + o.field_mass_kg + o.bound_debris_kg == pytest.approx(2 * EARTH)
        assert o.bound_debris_kg > 0.0 and sum(o.field) == pytest.approx(o.field_mass_kg)
        assert o.field.core == 0.0                                                          # the cores stay out of the pieces that fly free


def test_a_graze_breaks_less_than_a_head_on_impact():
    head_on = fc.collision_field(DRY, DRY, 28e3, impact_parameter=0.0)
    graze = fc.collision_field(DRY, DRY, 28e3, impact_parameter=0.6)
    assert graze.largest_fraction > head_on.largest_fraction and graze.field_mass_kg <= head_on.field_mass_kg


def test_the_makeup_is_the_shards_in_rows_of_the_materials_table():
    o = fc.collision_field(WET, WET, 35e3)
    makeup = fc.field_makeup(o.field, 0.2)
    assert isinstance(makeup, small_body_composition.SmallBodyMakeup) and makeup.kind == "asteroid-field"
    assert math.fsum(makeup.components.values()) == pytest.approx(1.0)
    assert makeup.water > 0.1 and makeup.metal == 0.0 and makeup.rock > 0.1
    assert makeup.porosity == 0.2 and makeup.bulk_density_kg_m3 == pytest.approx(makeup.grain_density_kg_m3 * 0.8)
    assert makeup.detail["ice_share"] == pytest.approx(makeup.water)
    assert fc.field_makeup(fc.field_layers(DRY, DRY, 0.13), 0.0).water == 0.0
    with pytest.raises(fc.FieldCompositionError):
        fc.field_makeup(o.field, 1.0)
    with pytest.raises(fc.FieldCompositionError):
        fc.field_makeup(fc.Layers(0.0, 0.0, 0.0), 0.2)


def test_bad_arguments_are_refused():
    with pytest.raises(fc.FieldCompositionError):
        fc.collision_field(DRY, DRY, 0.0)
    with pytest.raises(fc.FieldCompositionError):
        fc.collision_field(DRY, DRY, 25e3, impact_parameter=1.5)
