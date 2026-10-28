# tests/test_body_colour.py

"""
DB.96: the atmosphere colour a planet or moon stores as one sRGB integer (`physics/body_colour.py`).
"""

from types import SimpleNamespace

import pytest

from planetgen.physics import aerosols, atmosphere, body_colour


def _body(**kwargs):
    base = dict(body_type="t", gravity=1.0, atmospheric_pressure=101325.0, planet_class="M", ozone_column_du=300.0,
                star=SimpleNamespace(temperature=5772.0),
                tau_haze=None, tau_dust=None, tau_cloud=None, cloud_cover=None, aerosol_draw=None)
    base.update(kwargs)
    body = SimpleNamespace(**base)
    setattr(body, atmosphere.GAS_FIELD, {"N2": 78084.0, "O2": 20946.0, "Ar": 934.0, "CO2": 41.0})
    return body


def test_the_integer_is_the_three_bytes_in_red_green_blue_order():
    assert body_colour.rgb_integer((1, 2, 3)) == 0x010203
    assert body_colour.rgb8_of(0xABCDEF) == (0xAB, 0xCD, 0xEF)
    assert body_colour.rgb8_of(body_colour.rgb_integer((0, 255, 7))) == (0, 255, 7)


def test_an_earth_like_sky_is_blue_and_is_stored_as_the_integer():
    body = _body()
    body_colour.generate(body)
    red, green, blue = body_colour.rgb8_of(body.atmosphere_colour)
    assert blue > green > red


def test_a_body_with_no_air_or_no_star_temperature_stores_none():
    for body in (_body(atmospheric_pressure=0.0), _body(star=SimpleNamespace(temperature=None)), _body(star=None)):
        setattr(body, atmosphere.GAS_FIELD, {} if body.atmospheric_pressure == 0.0 else getattr(body, atmosphere.GAS_FIELD))
        body_colour.generate(body)
        assert body.atmosphere_colour is None


def test_a_thick_cloud_deck_whitens_the_sky():
    clear = _body()
    cloudy = _body(tau_haze=0.0, tau_dust=0.0, tau_cloud=40.0, cloud_cover=1.0, aerosol_draw=0.5)
    body_colour.generate(clear)
    body_colour.generate(cloudy)
    saturation = [max(c) - min(c) for c in (body_colour.rgb8_of(clear.atmosphere_colour), body_colour.rgb8_of(cloudy.atmosphere_colour))]
    assert saturation[1] < saturation[0]


def test_the_stored_aerosols_become_layers_with_their_own_depths():
    haze, dust, cloud = aerosols.colour_layers(aerosols.Aerosols(0.3, 0.5, 2.0, 0.7, 0.4), body_colour._CENTRES)
    assert (haze.tau_550, dust.tau_550, cloud.tau_550) == (0.3, 0.5, 2.0)
    assert aerosols.colour_layers(aerosols.Aerosols(0.0, 0.0, 0.0, 0.0, 0.0), body_colour._CENTRES) == (None, None, None)


def test_a_giant_takes_the_colour_of_its_air_at_the_reference_level():
    giant = _body(body_type="g", gravity=2.5, planet_class="J", atmospheric_pressure=None)
    setattr(giant, atmosphere.GAS_FIELD, {"H2": 8.0e4, "He": 2.0e4, "CH4": 2.0e3})
    body_colour.generate(giant)
    assert giant.atmosphere_colour is None or 0 <= giant.atmosphere_colour <= 0xFFFFFF
