# tests/test_field_units.py

"""
Field declarations (GEN.570): every persisted attribute of a mass object class declares the unit the instance holds it
in (`MassObject.FIELD_UNITS`), the persisted list is that table's names, and `quantity(name)` builds the `Quantity`.
These tests fail when a persisted field has no declaration, when a declaration names a unit the registry does not hold,
when a name's unit suffix and its declaration disagree (Boss's card Q4), and when a declared unit measures something
other than the column the attribute is stored in (DB.76's column registry).
"""

import pytest

from planetgen.db import column_units, store
from planetgen.generation.binary import BinaryStarProxy
from planetgen.generation.config import SystemConfig
from planetgen.generation.mass_object import PLAIN, MassObject
from planetgen.generation.phenomena.compact_remnant import BlackHole, NeutronStar
from planetgen.generation.star import Star
from planetgen.generation.stellar_object import StellarObject
from planetgen.physics import unit_registry as units
from planetgen.physics.quantity import Quantity

STELLAR = (Star, BlackHole, NeutronStar, BinaryStarProxy)

# The kind table each class is stored in, with the `Field` list that writes and loads it and how many unit-valued
# attributes the comparison must reach (a remnant's own row holds only its own columns).
STORED = ((Star, "stars", store.STAR_FIELDS, 25), (NeutronStar, "neutron_stars", store.NEUTRON_STAR_FIELDS, 4),
          (BlackHole, "black_holes", store.BLACK_HOLE_FIELDS, 4))


class _Row(dict):
    """A stored row that answers 1.0 for any column, to learn which keys a `Field` loads."""

    def __missing__(self, key):
        return 1.0


@pytest.mark.parametrize("cls", STELLAR)
def test_the_persisted_list_is_the_declared_table(cls):
    assert cls.SERIALIZABLE_FIELDS == list(cls.FIELD_UNITS)
    assert len(set(cls.SERIALIZABLE_FIELDS)) == len(cls.SERIALIZABLE_FIELDS)


@pytest.mark.parametrize("cls", STELLAR)
def test_every_declaration_is_a_registry_unit_or_a_plain_word(cls):
    for name, unit in cls.FIELD_UNITS.items():
        assert unit in units.UNITS or unit in PLAIN, f"{cls.__name__}.{name}: {unit!r}"
    assert not set(PLAIN) & set(units.UNITS)


@pytest.mark.parametrize("cls", STELLAR)
def test_a_unit_in_a_name_is_the_unit_declared(cls):
    for name, unit in cls.FIELD_UNITS.items():
        stated = column_units.unit_matching_suffix(name)
        if stated is not None:
            assert unit == stated, f"{cls.__name__}.{name} is declared {unit!r}, its name says {stated!r}"


@pytest.mark.parametrize("cls, table, fields, minimum", STORED)
def test_a_declared_unit_measures_what_its_column_measures(cls, table, fields, minimum):
    checked = 0
    for field in fields:
        if field.load is None:
            continue
        for key in field.load(_Row()):
            unit = cls.FIELD_UNITS.get(key)
            if unit is None or unit in PLAIN:
                continue
            columns = [key] if key in field.columns else list(field.columns)
            kinds = {column_units.kind_of(table, column) for column in columns if (table, column) in column_units.COLUMNS}
            if kinds:
                assert units.unit(unit).kind in kinds, f"{cls.__name__}.{key}: {unit!r} is not {sorted(kinds)}"
                checked += 1
    assert checked >= minimum


def test_the_shared_tables_are_one_declaration():
    for cls in (BlackHole, NeutronStar):
        for name, unit in StellarObject.FIELD_UNITS.items():
            assert cls.FIELD_UNITS[name] == unit
    from planetgen.physics import activity, pre_main_sequence, spin
    for table in (spin.SPIN_UNITS, activity.STAR_ACTIVITY_UNITS, pre_main_sequence.DISK_UNITS):
        for name, unit in table.items():
            assert StellarObject.FIELD_UNITS[name] == unit
    assert spin.SPIN_FIELDS == tuple(spin.SPIN_UNITS)


def test_a_quantity_is_the_attribute_in_its_native_unit():
    star = Star(SystemConfig())
    star.mass = 1.98841e30
    star.radius = 695_700.0
    star.age = 4.6
    star.habitable_zone = (0.95, 1.7)
    mass = star.quantity("mass")
    assert isinstance(mass, Quantity) and (mass.kind, mass.unit, mass.value) == ("mass", "kg", 1.98841e30)
    assert mass.Msun == pytest.approx(1.0, rel=1e-4)
    assert star.quantity("radius").Rsun == pytest.approx(1.0, rel=1e-4)
    assert star.quantity("age").yr == pytest.approx(4.6e9, rel=1e-9)
    inner, outer = star.quantity("habitable_zone")
    assert (inner.kind, inner.unit, outer.value) == ("length", "au", 1.7)
    assert inner.km == pytest.approx(0.95 * 149_597_870.7)


def test_a_quantity_refuses_what_is_not_one():
    star = Star(SystemConfig())
    star.mass = 2e30
    star.metallicity_fe_h = -0.2
    star.a_crit_au = None
    with pytest.raises(ValueError, match="holds no quantity"):
        star.quantity("name")
    with pytest.raises(ValueError, match="holds no quantity"):
        star.quantity("metallicity_fe_h")
    with pytest.raises(ValueError, match="no value yet"):
        star.quantity("a_crit_au")
    with pytest.raises(KeyError):
        star.quantity("not_a_field")
    with pytest.raises(KeyError):
        Star.field_unit("not_a_field")
    assert Star.field_unit("mass") == "kg"


def test_a_remnant_and_a_binary_proxy_answer_in_their_own_units():
    hole = BlackHole(SystemConfig())
    hole.mass_solar = 10.0
    hole.event_horizon_radius_km = 29.5
    assert hole.quantity("mass_solar").kg == pytest.approx(1.98841e31, rel=1e-4)
    assert hole.quantity("event_horizon_radius_km").m == pytest.approx(29_500.0)
    neutron_star = NeutronStar(SystemConfig())
    neutron_star.spin_period_ms = 1.4
    assert neutron_star.quantity("spin_period_ms").s == pytest.approx(1.4e-3)
    assert neutron_star.quantity("magnetic_field_gauss").T == pytest.approx(neutron_star.magnetic_field_gauss * 1e-4)
    assert BinaryStarProxy.field_unit("_effective_mass") == "kg"
    assert BinaryStarProxy.field_unit("binary_mutual_orbital_speed_kms") == "km_s"


def test_the_base_class_declares_nothing_and_a_family_that_has_not_declared_yet_is_empty():
    assert MassObject.FIELD_UNITS == {} and MassObject.SERIALIZABLE_FIELDS == []


def _mass_object_classes(base=MassObject):
    for subclass in base.__subclasses__():
        yield subclass
        yield from _mass_object_classes(subclass)


# Classes whose fields wait for their own declaration slice of GEN.570 (the planetary and small-body families); a class
# leaves this list when it declares, and the test below fails if one declares and stays on it.
WAITING = {"Planet", "RoguePlanet", "Comet", "InterstellarComet", "AsteroidField"}


def test_a_class_with_persisted_fields_declares_them_or_waits_in_the_list():
    import planetgen.generation.comet  # noqa: F401
    import planetgen.generation.phenomena.asteroid_field  # noqa: F401
    import planetgen.generation.phenomena.rogue  # noqa: F401
    import planetgen.generation.planet  # noqa: F401
    seen = set()
    for cls in _mass_object_classes():
        if not cls.SERIALIZABLE_FIELDS:
            continue
        seen.add(cls.__name__)
        if cls.__name__ in WAITING:
            assert not cls.FIELD_UNITS, f"{cls.__name__} declares its fields: take it off WAITING"
        else:
            assert cls.SERIALIZABLE_FIELDS == list(cls.FIELD_UNITS), f"{cls.__name__} persists fields it does not declare"
    assert WAITING <= seen
