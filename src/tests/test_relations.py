# tests/test_relations.py

"""
The named laws (GEN.586): `planetgen.physics.relations` declares each law once, and the kernel (plain SI floats) and the
Quantity form come from that one declaration. These tests hold every law against a known value, that the two forms
agree bit for bit, that the dimension test finds a deliberately wrong formula, and that the declaration refuses
mistakes.
"""

import math

import pytest

from planetgen.physics import constants
from planetgen.physics import relations
from planetgen.physics import unit_registry as units
from planetgen.physics.quantity import Quantity
from planetgen.physics.relations import DimensionError, Law, Probe

_SUN = Quantity("mass", 1.0, "Msun")
_EARTH = Quantity("mass", 1.0, "Mearth")
_EARTH_RADIUS = Quantity("length", 1.0, "Rearth")
_AU = Quantity("length", 1.0, "au")


def _wrong(formula, output, **inputs):
    return Law("wrong", tuple(inputs.items()), output, formula)


def test_every_law_has_a_kernel_a_quantity_form_and_its_dimension():
    assert len(relations.LAWS) >= 15
    for name, declared in relations.LAWS.items():
        assert getattr(relations, name) is declared.formula
        assert callable(getattr(relations.quantities, name))
        declared.check_dimensions()
    assert relations.check_dimensions() == []


def test_a_wrong_formula_is_found_by_the_dimension_test():
    missing_square = _wrong(lambda mass, radius, *, c=relations.SI: c.G * mass / radius, "acceleration",
                            mass="mass", radius="length")
    with pytest.raises(DimensionError, match="acceleration"):
        missing_square.check_dimensions()
    missing_constant = _wrong(lambda mass, radius, *, c=relations.SI: mass / radius ** 2, "acceleration",
                              mass="mass", radius="length")
    with pytest.raises(DimensionError):
        missing_constant.check_dimensions()
    wrong_root = _wrong(lambda mass, radius, *, c=relations.SI: (c.G * mass / radius) ** (1 / 3), "speed",
                        mass="mass", radius="length")
    with pytest.raises(DimensionError):
        wrong_root.check_dimensions()
    unlike_sum = _wrong(lambda mass, radius, *, c=relations.SI: mass + radius, "mass", mass="mass", radius="length")
    with pytest.raises(DimensionError, match="cannot add"):
        unlike_sum.check_dimensions()
    right = _wrong(lambda mass, radius, *, c=relations.SI: c.G * mass / radius ** 2, "acceleration",
                   mass="mass", radius="length")
    right.check_dimensions()


def test_a_wrong_law_is_reported_by_the_gate_function(monkeypatch):
    wrong = _wrong(lambda mass, *, c=relations.SI: mass, "length", mass="mass")
    monkeypatch.setitem(relations.LAWS, "wrong", wrong)
    assert relations.check_dimensions() == ["wrong"]


def test_probes_follow_the_arithmetic():
    length, mass = Probe(units.dimension(length=1)), Probe(units.dimension(mass=1))
    assert (length * mass).vector == Probe(units.dimension(length=1, mass=1)).vector
    assert (length / length).vector == Probe().vector
    assert (1 / length).vector == Probe(units.dimension(length=-1)).vector
    assert relations.sqrt(length * length).vector == length.vector
    assert (length ** (1 / 3)).vector[0] == pytest.approx(1 / 3)
    assert relations.log(length / length).vector == Probe().vector
    with pytest.raises(DimensionError):
        relations.log(length)
    with pytest.raises(DimensionError):
        length + mass
    with pytest.raises(DimensionError):
        length + 1.0
    assert (Probe() + 1.0).vector == Probe().vector
    with pytest.raises(TypeError):
        length * "x"
    assert relations.sqrt(4.0) == 2.0 and relations.log(1.0) == 0.0


@pytest.mark.parametrize("name, args, expected, rel", [
    ("schwarzschild_radius", (constants.SOLAR_MASS_TO_KG,), 2953.25, 1e-4),
    ("surface_gravity", (constants.EARTH_MASS_TO_KG, 6.371e6), 9.82, 5e-3),
    ("escape_speed", (constants.EARTH_MASS_TO_KG, 6.371e6), 11186.0, 1e-3),
    ("circular_speed", (constants.EARTH_MASS_TO_KG, 6.371e6), 7905.0, 1e-3),
    ("orbital_period", (constants.AU_M, constants.SOLAR_MASS_TO_KG), constants.SECONDS_PER_YEAR, 1e-4),
    ("mean_motion", (constants.AU_M, constants.SOLAR_MASS_TO_KG), 2 * math.pi / constants.SECONDS_PER_YEAR, 1e-4),
    ("semi_major_axis", (constants.SECONDS_PER_YEAR, constants.SOLAR_MASS_TO_KG), constants.AU_M, 1e-4),
    ("rotation_rate", (86164.1,), 7.2921e-5, 1e-4),
    ("rotation_period", (7.2921e-5,), 86164.1, 1e-4),
    ("hill_radius", (constants.AU_M, constants.EARTH_MASS_TO_KG, constants.SOLAR_MASS_TO_KG), 1.4967e9, 2e-3),
    ("roche_limit", (6.371e6, 5514.0, 3344.0), 2.44 * 6.371e6 * (5514.0 / 3344.0) ** (1 / 3), 1e-12),
    ("sphere_volume", (6.371e6,), 1.0833e21, 1e-3),
    ("mean_density", (constants.EARTH_MASS_TO_KG, 6.371e6), 5514.0, 1e-3),
    ("luminosity", (constants.SOLAR_RADIUS_M, 5772.0), 3.828e26, 1e-2),
    ("effective_temperature", (constants.SOLAR_LUMINOSITY, constants.SOLAR_RADIUS_M), 5772.0, 1e-2),
    ("inverse_square_flux", (constants.SOLAR_LUMINOSITY, constants.AU_M), 1361.0, 1e-2),
    ("equilibrium_temperature", (constants.SOLAR_LUMINOSITY, constants.AU_M, 0.3), 254.6, 1e-3),
    ("mass_flow", (constants.SOLAR_MASS_TO_KG * 2e-14 * constants.SECONDS_PER_YEAR, constants.SECONDS_PER_YEAR),
     constants.SOLAR_MASS_TO_KG * 2e-14, 1e-12),
])
def test_a_law_gives_the_textbook_value(name, args, expected, rel):
    assert getattr(relations, name)(*args) == pytest.approx(expected, rel=rel)


def test_the_laws_that_invert_one_another_do():
    radius = 6.371e6
    mass = constants.EARTH_MASS_TO_KG
    assert relations.sphere_radius(mass, relations.mean_density(mass, radius)) == pytest.approx(radius, rel=1e-12)
    temperature = 5772.0
    power = relations.luminosity(radius, temperature)
    assert relations.effective_temperature(power, radius) == pytest.approx(temperature, rel=1e-12)
    assert relations.rotation_period(relations.rotation_rate(86400.0)) == pytest.approx(86400.0, rel=1e-12)
    period = relations.orbital_period(constants.AU_M, constants.SOLAR_MASS_TO_KG)
    assert relations.semi_major_axis(period, constants.SOLAR_MASS_TO_KG) == pytest.approx(constants.AU_M, rel=1e-12)
    assert relations.mean_motion(constants.AU_M, constants.SOLAR_MASS_TO_KG) * period == pytest.approx(2 * math.pi)
    assert relations.mean_density(mass, radius) * relations.sphere_volume(radius) == pytest.approx(mass, rel=1e-12)


def test_vis_viva_on_a_circle_is_the_circular_speed_and_at_infinity_the_escape_speed():
    mass, radius = constants.SOLAR_MASS_TO_KG, constants.AU_M
    assert relations.vis_viva_speed(radius, radius, mass) == pytest.approx(relations.circular_speed(mass, radius),
                                                                           rel=1e-12)
    assert relations.vis_viva_speed(radius, math.inf, mass) == pytest.approx(relations.escape_speed(mass, radius),
                                                                              rel=1e-12)


def test_a_breakup_rate_is_the_circular_speed_over_the_radius():
    mass, radius = constants.EARTH_MASS_TO_KG, 6.371e6
    assert relations.breakup_rate(mass, radius) == pytest.approx(relations.circular_speed(mass, radius) / radius,
                                                                  rel=1e-12)


def test_the_quantity_form_and_the_kernel_agree_bit_for_bit():
    shown = {"mass": (_SUN, _EARTH, Quantity("mass", 3.2e5, "g")), "length": (_AU, _EARTH_RADIUS, Quantity("length", 7.0, "km")),
             "duration": (Quantity("duration", 3.0, "d"), Quantity("duration", 0.7, "yr")),
             "temperature": (Quantity("temperature", 300.0, "K"), Quantity("temperature", 10.0, "degC")),
             "power": (Quantity("power", 1.0, "Lsun"),), "density": (Quantity("density", 5.5, "g_cm3"),),
             "angular_rate": (Quantity("angular_rate", 2.0, "turn_d"),)}
    for name, declared in relations.LAWS.items():
        for index in range(2):
            given = []
            for argument, kind in declared.inputs:
                if kind is None:
                    given.append(0.3)
                else:
                    options = shown[kind]
                    given.append(options[index % len(options)])
            bases = [value if isinstance(value, float) else value.base() for value in given]
            try:
                kernel = getattr(relations, name)(*bases)
            except (ValueError, ZeroDivisionError):
                continue
            answer = getattr(relations.quantities, name)(*given)
            assert answer.kind == declared.output
            assert answer.unit == units.KINDS[declared.output].base
            assert answer.value == kernel, name


def test_the_quantity_form_answers_in_every_unit_of_its_kind():
    gravity = relations.quantities.surface_gravity(_EARTH, _EARTH_RADIUS)
    assert gravity.kind == "acceleration" and gravity.unit == "m_s2"
    assert gravity.g0 == pytest.approx(1.0014, rel=1e-4)
    assert gravity.ft_s2 == pytest.approx(gravity.m_s2 / 0.3048)
    named = relations.quantities.surface_gravity(mass=_EARTH, radius=_EARTH_RADIUS)
    assert named == gravity
    temperature = relations.quantities.equilibrium_temperature(Quantity("power", 1.0, "Lsun"), _AU, albedo=0.3)
    assert temperature.degC == pytest.approx(temperature.K - 273.15)


def test_the_quantity_form_refuses_the_wrong_kind_and_the_wrong_arguments():
    with pytest.raises(ValueError, match="radius must be a Quantity of kind length"):
        relations.quantities.surface_gravity(_EARTH, _EARTH)
    with pytest.raises(ValueError, match="mass must be a Quantity"):
        relations.quantities.surface_gravity(6.0e24, _EARTH_RADIUS)
    with pytest.raises(TypeError, match="takes mass, radius"):
        relations.quantities.surface_gravity(_EARTH)
    with pytest.raises(TypeError, match="takes mass, radius"):
        relations.quantities.surface_gravity(_EARTH, _EARTH_RADIUS, _EARTH_RADIUS)
    with pytest.raises(TypeError, match="unexpected or repeated"):
        relations.quantities.surface_gravity(_EARTH, radius=_EARTH_RADIUS, mass=_EARTH)
    with pytest.raises(TypeError, match="unexpected or repeated"):
        relations.quantities.surface_gravity(_EARTH, size=_EARTH_RADIUS)


def test_a_declaration_names_each_argument_once_and_every_kind_exists():
    with pytest.raises(ValueError, match="are named for"):
        @relations.law("length", mass="mass")
        def stray(radius, *, c=relations.SI):
            return radius
    with pytest.raises(ValueError, match="declared twice"):
        @relations.law("length", mass="mass")
        def schwarzschild_radius(mass, *, c=relations.SI):
            return mass
    with pytest.raises(ValueError, match="unknown kind"):
        @relations.law("beauty", mass="mass")
        def vain(mass, *, c=relations.SI):
            return mass
    assert "stray" not in relations.LAWS and "vain" not in relations.LAWS


def test_the_constants_come_from_the_registry_table():
    for name, (value, vector) in units.PHYSICAL_CONSTANTS.items():
        assert getattr(relations.SI, name) == value
        assert getattr(relations.PROBES, name).vector == Probe(vector).vector
    assert relations.SI.turn == 2 * math.pi and relations.SI.radian == 1.0


def test_an_orbit_in_au_years_and_solar_masses_follows_from_the_same_law():
    period_s = relations.orbital_period(constants.AU_M, constants.SOLAR_MASS_TO_KG)
    in_years = units.in_system(period_s, units.dimension(time=1), "au-yr-msun")
    assert in_years == pytest.approx(1.0, rel=1e-4)
    gauss = units.constant_in("G", "au-yr-msun")
    assert 2 * math.pi / math.sqrt(gauss) == pytest.approx(in_years, rel=1e-12)
