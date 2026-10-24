# tests/test_describe_quantities.py

"""
GEN.348: the describer and the page spans hold a value in the unit it is stored in, as a `Quantity`, and no call site
multiplies by a unit factor. The Document still carries `si` in the ladder's base unit, so each helper's number is checked
against the registry's own factor, and a bad stored value gives a dash on a page and never raises.
"""

import math

import pytest

from planetgen.describe import build, model
from planetgen.describe.quantity import display_value, quantity_text
from planetgen.physics import constants
from planetgen.physics.quantity import Quantity
from planetgen.web.lib import fmt


def _qty(item):
    return item["qty"]


@pytest.mark.parametrize("helper, value, kind, si", [
    (build._au, 2.5, "length", 2.5 * constants.AU_M),
    (build._km, 384_400.0, "length", 384_400.0 * 1000.0),
    (build._radius_km, 6371.0, "radius", 6371.0 * 1000.0),
    (build._years, 11.86, "time", 11.86 * constants.SECONDS_PER_YEAR),
    (build._gy, 4.6, "time", 4.6 * 1e9 * constants.SECONDS_PER_YEAR),
    (build._kms, 29.78, "speed", 29.78 * 1000.0),
    (build._gravity, 0.38, "gravity", 0.38 * constants.EARTH_GRAVITY),
])
def test_a_helper_holds_the_stored_unit_and_writes_the_base_unit(helper, value, kind, si):
    item = _qty(helper(value))
    assert item["kind"] == kind and display_value(item) == pytest.approx(si, rel=1e-14)
    assert _qty(helper(None)) == {"kind": kind, "value": None, "unit": model.DISPLAY_UNITS[kind]}


def test_a_mass_and_a_luminosity_carry_their_ratio_to_the_sun():
    kilograms = 2.0 * constants.SOLAR_MASS_TO_KG
    assert _qty(build._mass(kilograms, solar=True))["sol"] == pytest.approx(2.0, rel=1e-14)
    assert "sol" not in _qty(build._mass(kilograms))
    held = _qty(build._mass_solar(0.5))
    assert display_value(held) == pytest.approx(0.5 * constants.SOLAR_MASS_TO_KG, rel=1e-14) and held["sol"] == 0.5
    assert _qty(build._luminosity(constants.SOLAR_LUMINOSITY))["sol"] == pytest.approx(1.0, rel=1e-14)


def test_an_offset_is_a_length_in_au():
    assert build._offset(3.0, 4.0, 0.0) == 5.0
    assert display_value(_qty(build._au(build._offset(3.0, 4.0, 0.0)))) == pytest.approx(5.0 * constants.AU_M, rel=1e-14)


def test_a_quantity_item_holds_the_value_in_the_unit_it_was_stored_in():
    assert _qty(build._au(2.5)) == {"kind": "length", "value": 2.5, "unit": "au"}
    assert _qty(build._gy(4.6)) == {"kind": "time", "value": 4.6, "unit": "Gyr"}
    assert _qty(build._radius_km(6371.0)) == {"kind": "radius", "value": 6371.0, "unit": "km"}
    assert _qty(model.qty("ratio", 0.5)) == {"kind": "ratio", "value": 0.5, "unit": "ratio"}
    with pytest.raises(ValueError):
        model.qty("length", 1.0, "furlong")


def test_a_measured_item_keeps_the_unit_the_quantity_was_held_in_and_the_text_converts_it():
    held_in_megametres = model.measured("length", Quantity("length", 384.4, "Mm"))
    assert held_in_megametres == model.qty("length", 384.4, "Mm")
    assert display_value(_qty(held_in_megametres)) == pytest.approx(384_400_000.0, rel=1e-14)
    assert quantity_text(_qty(held_in_megametres)) == quantity_text(_qty(model.qty("length", 384_400_000.0)))
    assert model.measured("length", None) == model.qty("length", None)
    with pytest.raises(ValueError):
        model.measured("not a kind", Quantity("length", 1.0, "m"))


def test_a_quasar_broad_line_region_is_held_in_light_days():
    item = _qty(build._held("ld", 30.0))
    assert display_value(item) == pytest.approx(30.0 * 86_400.0 * constants.SPEED_OF_LIGHT_M_S, rel=1e-14)
    assert quantity_text(item) == "2.52 cpc (0.0821 ly)"


def test_a_body_radius_is_written_in_kilometers_from_meters():
    assert quantity_text({"kind": "radius", "value": 6_371_000.0, "unit": "m"}) == "6.37 × 10³ km"


@pytest.mark.parametrize("value", [None, float("nan"), float("inf"), -1.0, -5e10, 10 ** 400])
def test_a_page_span_of_a_value_its_kind_refuses_is_a_dash(value):
    assert fmt.measured_html("km", value) == "\u2013"
    assert fmt.measured_html("Msun", value) == "\u2013"


def test_a_speed_at_or_above_light_is_a_dash_not_a_traceback():
    assert fmt.measured_html("km_s", 299_792.458) == "\u2013"
    assert fmt.measured_html("km_s", 100.0) != "\u2013"


def test_a_page_span_is_the_same_whatever_unit_the_value_was_held_in():
    assert fmt.quantity_html(Quantity("length", 1.5e8, "km")) == fmt.measured_html("au", 1.5e8 / 149_597_870.7)
    span = str(fmt.measured_html("Gyr", 4.6))
    assert span.startswith('<span class="qty" data-q="time" data-si="') and span.endswith(">4.6 Gyr</span>")
    assert math.isclose(float(span.split('data-si="')[1].split('"')[0]), 4.6e9 * constants.SECONDS_PER_YEAR, rel_tol=1e-14)
