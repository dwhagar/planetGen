# tests/test_unit_registry.py

"""
The unit registry (GEN.568): `planetgen.physics.unit_registry` holds every unit once, with its kind, factor, offset,
symbol, ASCII key and name. These tests hold the table itself (every factor against its definition written out in
full, compound units against their parts, kinds against their dimensions), the conversions (round trips, offsets,
other kinds refused), the bounds, the declared relations, and that every rung of every display ladder
(`planetgen.util.ladders`) is a unit of the registry with the same factor, so the two cannot drift apart.
"""

import math

import pytest
from hypothesis import given
from hypothesis import strategies as st

from planetgen.physics import constants
from planetgen.physics import unit_registry as units
from planetgen.util import ladders
from tests import fuzz_support  # noqa: F401  (loads the hypothesis profile)

_JULIAN_YEAR_S = 365.25 * 86400
_AU_M = 149_597_870_700
_PARSEC_M = _AU_M * 648_000 / math.pi
_LIGHT_YEAR_M = 299_792_458 * _JULIAN_YEAR_S

# Every factor and offset stated here from its definition, not read from `constants`.
DEFINITIONS = {
    "m": 1.0, "km": 1e3, "Mm": 1e6, "Gm": 1e9, "mm": 1e-3, "um": 1e-6, "nm": 1e-9,
    "au": _AU_M, "ly": _LIGHT_YEAR_M, "pc": _PARSEC_M, "mpc": _PARSEC_M * 1e-3, "cpc": _PARSEC_M * 1e-2,
    "kpc": _PARSEC_M * 1e3, "Mpc": _PARSEC_M * 1e6, "Gpc": _PARSEC_M * 1e9,
    "ft": 0.3048, "mi": 5280 * 0.3048,
    "kg": 1.0, "g": 1e-3, "lb": 0.45359237,
    "s": 1.0, "min": 60.0, "h": 3600.0, "d": 86400.0, "yr": _JULIAN_YEAR_S, "kyr": _JULIAN_YEAR_S * 1e3,
    "Myr": _JULIAN_YEAR_S * 1e6, "Gyr": _JULIAN_YEAR_S * 1e9, "ms": 1e-3, "us": 1e-6, "ns": 1e-9,
    "m_s": 1.0, "km_s": 1e3, "km_h": 1e3 / 3600, "mph": 5280 * 0.3048 / 3600, "c": 299_792_458.0,
    "au_d": _AU_M / 86400, "au_yr": _AU_M / _JULIAN_YEAR_S,
    "K": 1.0, "degC": 1.0, "degF": 5 / 9,
    "rad": 1.0, "deg": math.pi / 180, "arcmin": math.pi / 10800, "arcsec": math.pi / 648000, "turn": 2 * math.pi,
    "W": 1.0, "kW": 1e3, "erg_s": 1e-7, "erg": 1e-7, "J": 1.0, "eV": 1.602176634e-19,
    "Pa": 1.0, "kPa": 1e3, "bar": 1e5, "atm": 101325.0, "psi": 0.45359237 * 9.80665 / 0.0254 ** 2,
    "ld": 299_792_458 * 86_400.0, "g0": 9.80665, "m_s2": 1.0, "ft_s2": 0.3048,
    "T": 1.0, "G": 1e-4, "W_m2": 1.0, "kg_m3": 1.0, "g_cm3": 1e3, "m-3": 1.0, "cm-3": 1e6,
    "m3_s2": 1.0, "km3_s2": 1e9,
    "m2": 1.0, "km2": 1e6, "cm2": 1e-4, "m3": 1.0, "km3": 1e9, "cm3": 1e-6, "L": 1e-3,
}
OFFSETS = {"degC": 273.15, "degF": 459.67 * 5 / 9}


def test_every_defined_factor_is_the_registry_factor():
    for key, factor in DEFINITIONS.items():
        assert units.unit(key).factor == pytest.approx(factor, rel=1e-12), key
    for key, offset in OFFSETS.items():
        assert units.unit(key).offset == pytest.approx(offset, rel=1e-12), key


def test_the_astronomical_masses_and_radii_read_the_constants():
    assert units.unit("Msun").factor == constants.SOLAR_MASS_TO_KG
    assert units.unit("Mearth").factor == constants.EARTH_MASS_TO_KG
    assert units.unit("Mjup").factor == constants.JUPITER_MASS_TO_KG
    assert units.unit("Mmoon").factor == constants.MOON_MASS_KG
    assert units.unit("Lsun").factor == constants.SOLAR_LUMINOSITY
    assert units.unit("Rsun").factor == constants.SOLAR_RADIUS_M
    assert units.unit("Rearth").factor == pytest.approx(6.371e6)
    assert units.unit("Rjup").factor == pytest.approx(7.1492e7)


def test_the_registry_has_a_unit_for_each_prefix_and_no_key_twice():
    assert units.unit("uW").factor == 1e-6
    assert units.unit("uW").symbol == "µW"
    assert units.unit("GW").factor == 1e9
    assert len(units.UNITS) == len({row.key for row in units.UNITS.values()})
    with pytest.raises(ValueError, match="declared twice"):
        units._unit("km", "length", "km", "kilometers", 1e3)


def test_a_key_is_ascii_and_a_symbol_is_unique_within_its_kind():
    for row in units.UNITS.values():
        assert row.key.isascii() and " " not in row.key, row.key
    for name in units.KINDS:
        symbols = [row.symbol for row in units.units_of(name)]
        assert len(symbols) == len(set(symbols)), name


def test_every_kind_has_a_base_unit_of_factor_one():
    for name, kind in units.KINDS.items():
        base = units.unit(kind.base)
        assert base.kind == name
        assert (base.factor, base.offset) == (1.0, 0.0), name


def test_a_compound_unit_has_the_dimension_and_factor_of_its_parts():
    compound = [key for key, row in units.UNITS.items() if row.parts]
    assert {"km_s", "au_d", "pc_Myr", "rad_s", "kg_m2_s", "Msun_pc3", "km3_s2"} <= set(compound)
    for key in compound:
        assert units.parts_dimension(key) == units.unit_dimension(key), key
        assert units.parts_factor(key) == pytest.approx(units.unit(key).factor, rel=1e-12), key


def test_a_kilometer_per_second_is_a_kilometer_over_a_second():
    assert units.convert(1.0, "km_s", "m_s") == units.convert(1.0, "km", "m") / units.convert(1.0, "s", "s")
    assert units.convert(1.0, "pc_Myr", "km_s") == pytest.approx(0.97779, rel=1e-4)


def test_dimensions_of_the_kinds_that_derive_from_others():
    d = units.dimension
    assert units.KINDS["speed"].dimension == units.multiply_dimensions(d(length=1), d(time=1), -1)
    assert units.KINDS["power"].dimension == units.multiply_dimensions(units.KINDS["energy"].dimension, d(time=1), -1)
    assert units.KINDS["flux"].dimension == units.multiply_dimensions(units.KINDS["power"].dimension, d(length=2), -1)
    assert units.KINDS["pressure"].dimension == units.multiply_dimensions(units.KINDS["energy"].dimension,
                                                                         d(length=3), -1)
    assert units.KINDS["density"].dimension == units.multiply_dimensions(d(mass=1), d(length=3), -1)
    assert units.KINDS["angular_momentum"].dimension == d(length=2, mass=1, time=-1)
    assert units.KINDS["angular_rate"].dimension == d(angle=1, time=-1)
    assert units.KINDS["gravitational_parameter"].dimension == d(length=3, time=-2)


@pytest.mark.parametrize("value, source, target, expected", [
    (100.0, "degC", "K", 373.15), (0.0, "degC", "degF", 32.0), (212.0, "degF", "degC", 100.0),
    (-40.0, "degC", "degF", -40.0), (1.0, "au", "km", 149_597_870.7), (1.0, "Msun", "Mjup", 1047.57),
    (1.0, "yr", "d", 365.25), (180.0, "deg", "rad", math.pi), (1.0, "atm", "kPa", 101.325),
    (1.0, "turn", "deg", 360.0), (1.0, "ly", "pc", 0.306601), (365.25, "ld", "ly", 1.0),
])
def test_conversions(value, source, target, expected):
    assert units.convert(value, source, target) == pytest.approx(expected, rel=1e-5)


def test_a_difference_is_not_a_temperature():
    assert units.convert(10.0, "dK", "ddegC") == 10.0
    assert units.convert(9.0, "ddegF", "dK") == pytest.approx(5.0)
    with pytest.raises(ValueError, match="cannot convert"):
        units.convert(10.0, "dK", "K")


def test_a_unit_of_another_kind_is_refused():
    with pytest.raises(ValueError, match="cannot convert"):
        units.convert(1.0, "km", "kg")
    with pytest.raises(ValueError, match="unknown unit"):
        units.convert(1.0, "furlong", "km")
    with pytest.raises(ValueError, match="unknown kind"):
        units.units_of("beauty")


@given(st.floats(min_value=-1e12, max_value=1e12, allow_nan=False),
       st.sampled_from(sorted(units.UNITS)))
def test_a_value_round_trips_through_the_base_unit(value, key):
    assert units.from_base(units.to_base(value, key), key) == pytest.approx(value, rel=1e-9, abs=1e-9)


@given(st.floats(min_value=-1e9, max_value=1e9, allow_nan=False),
       st.sampled_from(sorted(units.UNITS)), st.data())
def test_a_conversion_and_its_reverse_agree(value, key, data):
    other = data.draw(st.sampled_from([row.key for row in units.units_of(units.unit(key).kind)]))
    forward = units.convert(value, key, other)
    assert units.convert(forward, other, key) == pytest.approx(value, rel=1e-9, abs=1e-6)


def test_the_bounds_of_a_kind():
    assert units.check(5.0, "Msun") == 5.0
    assert units.check(0.0, "km") == 0.0
    assert units.check(-40.0, "degC") == -40.0
    for bad in (-1.0, math.nan, math.inf):
        with pytest.raises(ValueError):
            units.check(bad, "kg")
    with pytest.raises(ValueError, match="at least 0"):
        units.check(-300.0, "degC")
    with pytest.raises(ValueError, match="below"):
        units.check(1.0, "c")
    assert units.check(0.5, "c") == 0.5
    assert units.check(-5.0, "ddegC") == -5.0
    assert units.check(-5.0, "unix_s") == -5.0


def test_the_relations_are_inverse_pairs_in_base_units():
    for name, relation in units.RELATIONS.items():
        for sample in (0.3, 1.0, 7.5e4):
            there = units.apply_relation(name, sample, relation.first_role)
            assert units.apply_relation(name, there, relation.second_role) == pytest.approx(sample, rel=1e-12), name
    mass = 1.989e30
    assert units.apply_relation("mass-mu", mass, "mass") == pytest.approx(1.327e20, rel=1e-3)
    assert units.apply_relation("radius-diameter", 3.0, "radius") == 6.0
    assert units.apply_relation("period-frequency", 0.5, "period") == 2.0
    assert units.apply_relation("rotation_period-angular_rate", 86400.0, "rotation_period") == pytest.approx(7.2722e-5,
                                                                                                          rel=1e-4)
    with pytest.raises(ValueError, match="not 'colour'"):
        units.apply_relation("mass-mu", 1.0, "colour")
    with pytest.raises(ValueError, match="unknown relation"):
        units.apply_relation("mass-charm", 1.0, "mass")


def test_every_relation_names_kinds_that_exist():
    for relation in units.RELATIONS.values():
        assert relation.first_kind in units.KINDS and relation.second_kind in units.KINDS


_LADDER_KINDS = {"gravity": "acceleration", "time": "duration", "number density": "number_density",
                 "magnetic field": "magnetic_field"}


@pytest.mark.parametrize("ladder", ladders.LADDERS, ids=lambda ladder: ladder.quantity)
def test_every_ladder_rung_is_a_unit_of_the_registry(ladder):
    kind = _LADDER_KINDS.get(ladder.quantity, ladder.quantity)
    rows = units.units_of(kind)
    for rung in ladder.rungs:
        assert any(row.factor == pytest.approx(rung.factor, rel=1e-12) and row.offset == pytest.approx(rung.offset)
                   for row in rows), f"{ladder.quantity}: {rung.symbol} ({rung.factor})"


def test_a_power_of_ten_ratio_is_one_correctly_rounded_step():
    assert units.convert(1.5, "km", "m") == 1500.0
    assert units.convert(1500.0, "m", "km") == 1.5
    assert units.convert(2.5, "um", "nm") == 2500.0
    assert units.scale(7.0, 1e3, 1e-6) == 7e9
    assert units.scale(7.0, 1.0, 1e3) == 7.0 / 1e3
    for value in (0.1, 0.7, 3.3, 123.456, 9.99e-3):
        assert units.convert(value, "km", "um") == value * 1e9
        assert units.convert(value, "um", "km") == value / 1e9
        assert units.convert(value, "g_cm3", "kg_m3") == value * 1e3


def test_a_ratio_that_is_not_a_power_of_ten_takes_the_two_factors():
    assert units.scale(2.0, 149_597_870_700.0, 1.0) == 2.0 * 149_597_870_700.0 / 1.0
    assert units.scale(1.0, 3.0, 1.0) == 3.0
    assert units.scale(5.0, 1.0, 0.0 + 4.0) == 1.25


# --- Dimensions and unit systems (GEN.585) --------------------------------------------------------------------------

def test_every_kind_has_a_five_integer_dimension():
    for name, kind in units.KINDS.items():
        assert len(kind.dimension) == len(units.DIMENSION_BASES), name
        assert all(isinstance(exponent, int) for exponent in kind.dimension), name


def test_the_derived_kinds_the_code_uses_are_listed_with_their_dimensions():
    d = units.dimension
    wanted = {
        "density": d(length=-3, mass=1), "acceleration": d(length=1, time=-2), "energy": d(length=2, mass=1, time=-2),
        "force": d(length=1, mass=1, time=-2), "mass_flow": d(mass=1, time=-1), "flux": d(mass=1, time=-3),
        "number_density": d(length=-3), "angular_momentum": d(length=2, mass=1, time=-1),
        "angular_rate": d(angle=1, time=-1), "power": d(length=2, mass=1, time=-3),
    }
    for name, vector in wanted.items():
        assert units.KINDS[name].dimension == vector, name
    # A force is a mass times an acceleration, a power a force times a speed.
    assert units.KINDS["force"].dimension == units.multiply_dimensions(units.KINDS["mass"].dimension,
                                                                      units.KINDS["acceleration"].dimension)
    assert units.KINDS["power"].dimension == units.multiply_dimensions(units.KINDS["force"].dimension,
                                                                      units.KINDS["speed"].dimension)
    assert units.convert(1.0, "dyn", "N") == 1e-5
    assert units.convert(1.0, "Msun_yr", "kg_s") == pytest.approx(constants.SOLAR_MASS_TO_KG / _JULIAN_YEAR_S)


def test_a_unit_system_is_four_registry_units_of_the_right_kinds():
    si, orbits = units.UNIT_SYSTEMS["si"], units.UNIT_SYSTEMS["au-yr-msun"]
    assert si.base_factors() == (1.0, 1.0, 1.0, 1.0)
    assert orbits.base_factors() == (units.unit("au").factor, units.unit("Msun").factor, units.unit("yr").factor, 1.0)
    for bad in (("kg", "kg", "s"), ("m", "m", "s"), ("m", "kg", "m"), ("m", "kg", "s", "degC")):
        with pytest.raises(ValueError, match="not a"):
            units._system("bad", *bad)
    assert "bad" not in units.UNIT_SYSTEMS


def test_a_constant_in_si_is_the_constant_itself():
    for name, (value, _vector) in units.PHYSICAL_CONSTANTS.items():
        assert units.constant_in(name) == value, name
        assert units.constant_in(name, "si") == value, name


def test_the_constants_in_au_year_and_solar_mass_follow_from_their_dimensions():
    au, msun, year = _AU_M, constants.SOLAR_MASS_TO_KG, _JULIAN_YEAR_S
    assert units.constant_in("G", "au-yr-msun") == pytest.approx(constants.G * msun * year ** 2 / au ** 3, rel=1e-14)
    assert units.constant_in("c", "au-yr-msun") == pytest.approx(299_792_458.0 * year / au, rel=1e-14)
    assert units.constant_in("sigma", "au-yr-msun") == pytest.approx(
        constants.STEFAN_BOLTZMANN_CONSTANT * year ** 3 / msun, rel=1e-14)
    # The typed 4 pi squared of `kepler.AU3_PER_YR2_PER_SOLAR_MASS` is G for the Sun in the Gaussian year, not in the Julian
    # year with this solar mass: the two differ by 3.8e-5 (its 39.4784 against the computed 39.4769), the difference this
    # note states.
    typed = 4 * math.pi ** 2
    difference = units.constant_in("G", "au-yr-msun") / typed - 1.0
    assert difference == pytest.approx(-3.78e-5, abs=1e-6)


def test_a_value_in_a_system_agrees_with_its_parts():
    d = units.dimension
    mu_sun_si = constants.G * constants.SOLAR_MASS_TO_KG
    assert units.in_system(mu_sun_si, d(length=3, time=-2), "au-yr-msun") == pytest.approx(
        units.constant_in("G", "au-yr-msun"), rel=1e-14)
    assert units.system_factor(d(), "au-yr-msun") == 1.0
    assert units.system_factor(d(angle=1), "au-yr-msun") == 1.0
    assert units.system_factor(units.KINDS["speed"].dimension, "au-yr-msun") == pytest.approx(
        units.unit("au_yr").factor, rel=1e-14)
    with pytest.raises(ValueError, match="unknown constant"):
        units.constant_in("alpha")
