# tests/test_unit_converters.py

"""
`unit_registry.converter` (GEN.571): a conversion the sweep moved from a hand-written factor to the registry gives the
same float bit for bit as the expression it replaced, on random inputs of every size a kernel meets. Each row is the
unit pair a swept module names and the old expression.
"""

import random

import pytest

from planetgen.physics import constants
from planetgen.physics import unit_registry as units

SECONDS_PER_YEAR = constants.SECONDS_PER_YEAR

# (from, to, the expression the sweep replaced)
REPLACED = [
    ("au", "m", lambda x: x * constants.AU_TO_M),
    ("m", "au", lambda x: x / constants.AU_M),
    ("au", "km", lambda x: x * constants.AU_TO_KM),
    ("km", "au", lambda x: x / constants.AU_TO_KM),
    ("ly", "m", lambda x: x * constants.LY_TO_M),
    ("m", "pc", lambda x: x / constants.PARSEC_M),
    ("yr", "s", lambda x: x * SECONDS_PER_YEAR),
    ("s", "yr", lambda x: x / SECONDS_PER_YEAR),
    ("Gyr", "yr", lambda x: x * 1.0e9),
    ("Gyr", "s", lambda x: x * (1e9 * SECONDS_PER_YEAR)),
    ("ms", "s", lambda x: x / 1000.0),
    ("km_s", "m_s", lambda x: x * 1000.0),
    ("km", "m", lambda x: x * 1000.0),
    ("m", "km", lambda x: x / 1000.0),
    ("pc", "m", lambda x: x * constants.PARSEC_M),
    ("yr", "d", lambda x: 365.25 * x),
    ("d", "yr", lambda x: x / 365.25),
    ("kpc", "m", lambda x: x * (constants.PARSEC_M * 1000.0)),
    ("km", "pc", lambda x: x / (constants.PARSEC_M / 1000.0)),
    ("ly", "au", lambda x: x * constants.LY_TO_AU),
    ("yr", "Gyr", lambda x: x / 1e9),
    ("au_yr", "km_s", lambda x: x * (constants.AU_TO_KM / SECONDS_PER_YEAR)),
]


def _inputs():
    generator = random.Random(571)
    values = [0.0, 1.0, 0.5, 1e-300, 1e300]
    values += [generator.uniform(0.0, 10.0) for _ in range(2000)]
    values += [10.0 ** generator.uniform(-12.0, 18.0) * generator.random() for _ in range(5000)]
    return values


@pytest.mark.parametrize("source,target,old", REPLACED, ids=[f"{s}->{t}" for s, t, _ in REPLACED])
def test_a_swept_conversion_is_bit_identical_to_the_expression_it_replaced(source, target, old):
    convert = units.converter(source, target)
    for value in _inputs():
        assert convert(value) == old(value), (value, source, target)


def test_a_converter_agrees_with_convert_to_rounding():
    for source, target, _old in REPLACED:
        convert = units.converter(source, target)
        assert convert(3.7) == pytest.approx(units.convert(3.7, source, target), rel=1e-14)


def test_units_of_another_kind_are_refused_when_the_converter_is_made():
    with pytest.raises(ValueError, match="cannot convert"):
        units.converter("km", "kg")


def test_temperature_scales_keep_their_offsets():
    to_c = units.converter("K", "degC")
    assert to_c(273.15) == pytest.approx(0.0, abs=1e-12)
    assert units.converter("m", "m")(5.5) == 5.5
