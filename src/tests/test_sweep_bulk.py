# tests/test_sweep_bulk.py

"""
TEST.210: the bulk composition sweep. The bodies of the shared sweep (`tests/audit/sweep.py`: every class and zone round three host
stars, and 24 generated systems) obey the rules of `tests/audit/sweep_bulk.py`: the shares make one; the details are shares; the density
and gravity follow from the mass and radius; the radius is the structure model's; the body is at least as dense as its materials at zero
pressure; the phrase names the materials the shares hold. The last test breaks each rule in turn and checks that the sweep reports it.

The first half tests each rule on synthetic bodies, where the answer is known.
"""

import types

import pytest

from planetgen.physics import bulk_composition, composition_text, constants
from tests.audit import sweep
from tests.audit import sweep_bulk as rules


def _rule(name):
    return next(rule for rule in rules.RULES if rule.name == name)


def _body(**changes):
    """A synthetic Earth-like planet: the shares, the details, and a radius, mass, density and gravity that agree."""
    mass, radius = 5.972e24, 6371.0
    volume = 4.0 / 3.0 * 3.141592653589793 * (radius * 1000.0) ** 3
    fields = dict(
        metal_mass_fraction=0.32, rock_mass_fraction=0.68, water_mass_fraction=0.0, ice_mass_fraction=0.0, envelope_mass_fraction=0.0,
        rock_feo_fraction=0.08, rock_al2o3_fraction=0.045, rock_mg_si_ratio=1.2, envelope_helium_fraction=None, envelope_z_fraction=None,
        porosity=None, carbon_rock=False, composition_origin="drawn", mass=mass, radius=radius, density=mass / volume / 1000.0,
        gravity=constants.G * mass / (radius * 1000.0) ** 2 / constants.EARTH_GRAVITY)
    fields.update(changes)
    body = types.SimpleNamespace(**fields)
    try:
        body.composition = composition_text.text(body)
    except bulk_composition.CompositionError:    # shares made wrong on purpose
        body.composition = "mostly silicate rock"
    return body


def _entry(body=None, moon=False):
    return sweep.Entry(body or _body(), object() if moon else None, None, "x")


def _found(name, body, moon=False):
    return _rule(name).check(_entry(body, moon))


# --- The rules on bodies made for them -----------------------------------------------------------------------------------


def test_an_earth_like_body_obeys_every_rule_but_the_radius_of_a_made_up_planet():
    """The synthetic planet is Earth's size and mass but not the structure model's Earth (it is a toy); every other rule holds."""
    body = _body()
    for rule in rules.RULES:
        if rule.name != "radius-from-structure":
            assert rule.check(_entry(body)) == [], rule.name


def test_shares_that_do_not_add_up_or_are_negative_are_reported():
    assert _found("fractions-sum", _body(metal_mass_fraction=0.5))
    assert _found("fractions-sum", _body(metal_mass_fraction=-0.1, rock_mass_fraction=1.1))
    assert _found("fractions-sum", _body(metal_mass_fraction=None))


def test_a_detail_that_is_not_a_share_or_an_envelope_detail_without_an_envelope_is_reported():
    assert _found("detail-ranges", _body(rock_feo_fraction=1.5))
    assert _found("detail-ranges", _body(envelope_helium_fraction=0.25, envelope_z_fraction=0.02))
    giant = dict(metal_mass_fraction=0.05, rock_mass_fraction=0.05, envelope_mass_fraction=0.9, rock_feo_fraction=None)
    assert _found("detail-ranges", _body(**giant)) != []                  # an envelope with no helium and no heavy elements
    assert _found("detail-ranges", _body(**giant, envelope_helium_fraction=0.25, envelope_z_fraction=0.05)) == []
    assert _found("detail-ranges", _body(**giant, envelope_helium_fraction=0.8, envelope_z_fraction=0.4))
    assert _found("detail-ranges", _body(rock_feo_fraction=None))         # rock with no FeO


def test_a_porosity_outside_zero_to_one_is_reported():
    assert _found("porosity-range", _body(porosity=1.2)) and _found("porosity-range", _body(porosity=-0.1))
    assert _found("porosity-range", _body(porosity=0.3)) == [] and _found("porosity-range", _body()) == []


def test_density_and_gravity_must_follow_the_mass_and_radius():
    assert _found("density-from-mass-radius", _body(density=6.0)) and _found("gravity-from-mass-radius", _body(gravity=2.0))
    assert _found("density-from-mass-radius", _body()) == [] and _found("gravity-from-mass-radius", _body()) == []


def test_the_radius_follows_the_structure_to_a_millionth_for_a_planet_and_to_four_percent_for_a_moon():
    from planetgen.physics import bulk_structure

    body = _body()
    bulk = bulk_composition.bulk_fractions(body)
    x_fe = bulk_structure.x_fe_from_feo_mass_fraction(bulk.rock_feo)
    radius = bulk_structure.structure(bulk, body.mass, x_fe=x_fe).radius_km
    assert _found("radius-from-structure", _body(radius=radius)) == []
    assert _found("radius-from-structure", _body(radius=radius * 1.001))
    assert _found("radius-from-structure", _body(radius=radius * 1.03), moon=True) == []
    assert _found("radius-from-structure", _body(radius=radius * 1.05), moon=True)
    assert _found("radius-from-structure", _body(radius=radius * 1.5, envelope_mass_fraction=0.0, metal_mass_fraction=0.32))


def test_a_body_with_an_envelope_is_left_out_of_the_radius_and_grain_rules():
    giant = _body(metal_mass_fraction=0.05, rock_mass_fraction=0.05, envelope_mass_fraction=0.9, radius=70000.0, density=0.1,
                  envelope_helium_fraction=0.25, envelope_z_fraction=0.05, rock_feo_fraction=None)
    assert _found("radius-from-structure", giant) == [] and _found("denser-than-grains", giant) == []
    assert rules._break_radius(_entry(giant)) is None and rules._break_denser_than_grains(_entry(giant)) is None


def test_a_body_lighter_than_its_materials_at_zero_pressure_is_reported():
    assert _found("denser-than-grains", _body(density=1.0))
    assert _found("denser-than-grains", _body(density=20.0)) == []
    assert _found("denser-than-grains", _body(porosity=0.5, metal_mass_fraction=0.0, rock_mass_fraction=1.0, density=1.7)) == []
    assert _found("denser-than-grains", _body(porosity=0.5, metal_mass_fraction=0.0, rock_mass_fraction=1.0, density=3.0 * 0.95)) == []
    assert _found("denser-than-grains", _body(porosity=0.5, metal_mass_fraction=0.0, rock_mass_fraction=1.0, density=1.0))


def test_the_phrase_must_name_the_materials_the_shares_hold_and_no_others():
    body = _body(metal_mass_fraction=0.0, rock_mass_fraction=0.5, water_mass_fraction=0.5, rock_feo_fraction=0.08)
    body.composition = composition_text.text(body)
    assert _found("text-materials", body) == []
    body.composition = "mostly silicate rock"
    assert _found("text-materials", body)
    body.composition = "mostly silicate rock and water, with a small share of iron-rich metal"
    assert any("iron-rich metal" in problem for problem in _found("text-materials", body))
    loose = _body(porosity=0.2)
    assert loose.composition.endswith("loosely packed") and _found("text-materials", loose) == []
    loose.composition = loose.composition.replace(", loosely packed", "")
    assert _found("text-materials", loose)


def test_every_rule_has_a_name_a_check_and_a_breaker_and_the_names_are_unique():
    assert len({rule.name for rule in rules.RULES}) == len(rules.RULES) >= 8
    assert all(callable(rule.check) and callable(rule.breaker) for rule in rules.RULES)


# --- The sweep ---------------------------------------------------------------------------------------------------------------


@pytest.fixture(scope="module")
def everything():
    return sweep.everything()


@pytest.mark.slow
def test_every_body_of_the_sweep_obeys_every_bulk_rule(everything):
    found = sweep.run_rules(rules.RULES, everything)
    failures = {name: violations[:5] for name, violations in found.items() if violations}
    assert not failures, "\n".join(f"{name}: {len(found[name])} bodies, first {violations}" for name, violations in failures.items())


@pytest.mark.slow
def test_the_sweep_holds_rocky_icy_porous_moon_and_enveloped_bodies_so_every_rule_is_exercised(everything):
    bodies = [entry.body for entry in everything]
    assert any(body.envelope_mass_fraction > 0.0 for body in bodies)
    assert any(getattr(body, "porosity", None) for body in bodies)
    assert any(entry.is_moon and entry.body.radius for entry in everything)
    assert any(body.water_mass_fraction > 0.0 for body in bodies) and any(body.ice_mass_fraction > 0.0 for body in bodies)


@pytest.mark.slow
def test_each_bulk_rule_reports_a_body_with_exactly_that_rule_broken(everything):
    result = sweep.prove_rules_can_fail(rules.RULES, everything)
    assert result == {rule.name: True for rule in rules.RULES}, result
