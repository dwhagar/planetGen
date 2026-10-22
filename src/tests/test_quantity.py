# tests/test_quantity.py

"""
The Quantity value class (GEN.569): `planetgen.physics.quantity.Quantity` keeps a number as given and answers it in every
other unit of its kind. These tests hold that every form is in step after any `set`, that a form is the truth scaled and
never one form from another, that a power-of-ten conversion is exact, the declared relations, the refusals (a unit of
another kind, a value outside the kind's bounds, an unknown name) and the dict, copy and pickle round trips.
"""

import copy
import math
import pickle

import pytest
from hypothesis import given
from hypothesis import strategies as st

from planetgen.physics import unit_registry as units
from planetgen.physics.quantity import Quantity
from tests import fuzz_support  # noqa: F401  (loads the hypothesis profile)

_VALUES = st.floats(min_value=1e-6, max_value=1e9, allow_nan=False)
_FLOAT_KINDS = [name for name, kind in units.KINDS.items() if kind.bounds.low is not None and name != "speed"]


def test_the_truth_is_kept_as_given():
    power = Quantity("power", 3.828e26, "W")
    assert (power.kind, power.value, power.unit) == ("power", 3.828e26, "W")
    assert power.W == 3.828e26
    assert power.Lsun == pytest.approx(1.0, rel=1e-15)
    assert power.kW == 3.828e23
    assert power.as_("erg_s") == pytest.approx(3.828e33, rel=1e-15)
    assert power.to_dict() == {"kind": "power", "value": 3.828e26, "unit": "W"}


def test_a_form_of_the_units_own_kind_is_the_truth_itself():
    quantity = Quantity("length", 0.1 + 0.2, "km")
    assert quantity.km == 0.1 + 0.2
    assert quantity.value == 0.1 + 0.2


def test_every_unit_key_of_a_kind_is_an_attribute_where_it_is_a_python_name():
    quantity = Quantity("mass", 2.0, "Msun")
    for row in units.units_of("mass"):
        assert getattr(quantity, row.key) == quantity.as_(row.key)
    number_density = Quantity("number_density", 3.0, "cm-3")
    assert number_density.as_("m-3") == 3e6


@given(st.sampled_from(_FLOAT_KINDS), _VALUES, st.data())
def test_every_form_follows_a_set(kind, value, data):
    keys = [row.key for row in units.units_of(kind) if row.offset == 0.0]
    first, second, third = (data.draw(st.sampled_from(keys)) for _ in range(3))
    quantity = Quantity(kind, value, first)
    before = quantity.as_(second)
    assert before == units.convert(value, first, second)
    quantity.set(value * 2.0, third)
    assert (quantity.value, quantity.unit) == (value * 2.0, third)
    for key in keys:
        assert quantity.as_(key) == units.convert(value * 2.0, third, key)
        assert getattr(quantity, key) == quantity.as_(key) if key.isidentifier() else True


def test_set_keeps_the_unit_when_none_is_given():
    quantity = Quantity("speed", 10.0, "km_s")
    quantity.set(20.0)
    assert (quantity.value, quantity.unit) == (20.0, "km_s")
    assert quantity.m_s == 20_000.0


def test_a_form_is_the_truth_scaled_not_one_form_from_another():
    quantity = Quantity("length", 1.0, "au")
    km = quantity.km
    assert km == pytest.approx(149_597_870.7, rel=1e-15)
    # Reading the km form first and then pc is the same as reading pc alone.
    assert quantity.pc == Quantity("length", 1.0, "au").pc


def test_a_power_of_ten_conversion_is_exact():
    assert Quantity("length", 1.5, "km").m == 1500.0
    assert Quantity("length", 1500.0, "m").km == 1.5
    assert Quantity("mass", 7.0, "kg").g == 7000.0
    for value in (0.1, 0.7, 3.3, 123.456):
        assert Quantity("length", value, "km").um == value * 1e9


def test_a_temperature_keeps_its_offset():
    quantity = Quantity("temperature", 100.0, "degC")
    assert quantity.K == pytest.approx(373.15)
    assert quantity.degF == pytest.approx(212.0)
    quantity.set(32.0, "degF")
    assert quantity.degC == pytest.approx(0.0, abs=1e-12)
    assert Quantity("delta_temperature", 10.0, "dK").ddegC == 10.0


def test_a_unit_of_another_kind_is_refused():
    with pytest.raises(ValueError, match="measures mass, not length"):
        Quantity("length", 1.0, "kg")
    quantity = Quantity("length", 1.0, "km")
    with pytest.raises(ValueError, match="cannot convert"):
        quantity.as_("kg")
    with pytest.raises(ValueError, match="measures mass, not length"):
        quantity.set(1.0, "Msun")
    assert (quantity.value, quantity.unit) == (1.0, "km")
    with pytest.raises(ValueError, match="unknown kind"):
        Quantity("beauty", 1.0, "km")
    with pytest.raises(ValueError, match="unknown unit"):
        Quantity("length", 1.0, "furlong")


def test_a_value_outside_the_kinds_bounds_is_refused():
    for bad in (-1.0, math.nan, math.inf):
        with pytest.raises(ValueError):
            Quantity("mass", bad, "kg")
    with pytest.raises(ValueError, match="below"):
        Quantity("speed", 1.0, "c")
    with pytest.raises(ValueError, match="at least 0"):
        Quantity("temperature", -300.0, "degC")
    quantity = Quantity("mass", 1.0, "kg")
    with pytest.raises(ValueError):
        quantity.set(-5.0)
    assert quantity.value == 1.0
    assert Quantity("instant", -5.0, "unix_s").value == -5.0


def test_an_unknown_attribute_is_an_attribute_error():
    quantity = Quantity("mass", 1.0, "kg")
    with pytest.raises(AttributeError):
        quantity.furlong
    with pytest.raises(AttributeError):
        quantity._nonsense
    assert not hasattr(quantity, "Lsun")
    assert not hasattr(quantity, "beauty")


def test_the_declared_relations():
    mass = Quantity("mass", 1.0, "Msun")
    mu = mass.related("mass-mu")
    assert (mu.kind, mu.unit) == ("gravitational_parameter", "m3_s2")
    assert mu.value == pytest.approx(1.32712e20, rel=1e-4)
    assert mu.related("mass-mu").kg == pytest.approx(mass.kg, rel=1e-12)
    radius = Quantity("length", 3.0, "km")
    assert radius.related("radius-diameter", "radius").m == 6000.0
    assert Quantity("length", 6.0, "km").related("radius-diameter", "diameter").m == 3000.0
    assert Quantity("duration", 0.5, "s").related("period-frequency").Hz == 2.0
    assert Quantity("frequency", 4.0, "Hz").related("period-frequency").s == 0.25
    rate = Quantity("duration", 1.0, "d").related("rotation_period-angular_rate")
    assert rate.rad_s == pytest.approx(2.0 * math.pi / 86400.0, rel=1e-12)


def test_a_relation_refuses_what_does_not_fit_it():
    with pytest.raises(ValueError, match="needs a role"):
        Quantity("length", 1.0, "km").related("radius-diameter")
    with pytest.raises(ValueError, match="needs a role: this mass could be neither side"):
        Quantity("mass", 1.0, "kg").related("period-frequency")
    with pytest.raises(ValueError, match="not 'colour'"):
        Quantity("mass", 1.0, "kg").related("mass-mu", "colour")
    with pytest.raises(ValueError, match="is a gravitational_parameter, not a mass"):
        Quantity("mass", 1.0, "kg").related("mass-mu", "mu")
    with pytest.raises(ValueError, match="unknown relation"):
        Quantity("mass", 1.0, "kg").related("mass-charm")


def test_the_dict_shape_round_trips():
    quantity = Quantity("pressure", 3.5, "atm")
    assert Quantity.from_dict(quantity.to_dict()) == quantity
    for bad in ({}, {"kind": "pressure"}, None, {"kind": "pressure", "value": 1.0, "unit": "kg"}):
        with pytest.raises(ValueError):
            Quantity.from_dict(bad)


def test_equality_across_units_and_no_hash():
    assert Quantity("length", 1500.0, "m") == Quantity("length", 1.5, "km")
    assert Quantity("length", 1.0, "au") != Quantity("length", 1.0, "ly")
    assert Quantity("length", 1.0, "m") != Quantity("mass", 1.0, "kg")
    assert Quantity("length", 1.0, "m") != 1.0
    with pytest.raises(TypeError):
        hash(Quantity("length", 1.0, "m"))


def test_text_and_slots():
    quantity = Quantity("speed", 12.5, "km_s")
    assert str(quantity) == "12.5 km/s"
    assert repr(quantity) == "Quantity('speed', 12.5, 'km_s')"
    assert not hasattr(quantity, "__dict__")
    with pytest.raises(AttributeError):
        quantity.colour = "red"


def test_copy_and_pickle_keep_the_truth_and_drop_the_forms():
    quantity = Quantity("length", 2.0, "pc")
    quantity.ly
    for clone in (copy.copy(quantity), copy.deepcopy(quantity), pickle.loads(pickle.dumps(quantity))):
        assert clone == quantity and clone is not quantity
        assert (clone.value, clone.unit) == (2.0, "pc")
        assert clone._forms is None
