"""MAP.86: a generated sector's Galaxy Map color (stellarObjects.sectorLook)."""

import colorsys

import pytest

from stellarObjects import sectorLook


def _hls(rgb):
    return colorsys.rgb_to_hls(*rgb)


def test_blackbody_colors_run_from_red_to_blue_white():
    red = sectorLook.blackbody_rgb(3000)
    sun = sectorLook.blackbody_rgb(5800)
    blue = sectorLook.blackbody_rgb(30000)
    assert red[0] == 1.0 and red[2] < 0.5
    assert sun[0] == 1.0 and sun[2] > 0.8
    assert blue[2] == 1.0 and blue[0] < blue[1] < blue[2]
    assert sectorLook.blackbody_rgb(None) == sun


def test_fill_share_is_log_scaled_from_empty_to_the_densest():
    assert sectorLook.fill_share(0, 500) == 0.0
    assert sectorLook.fill_share(None, 500) == 0.0
    assert sectorLook.fill_share(500, 500) == 1.0
    assert sectorLook.fill_share(900, 500) == 1.0
    assert 0.3 < sectorLook.fill_share(10, 500) < 0.5
    assert sectorLook.fill_share(3, None) == 1.0


def test_hue_follows_the_stars_saturation_the_density_and_lightness_the_luminosity():
    cool = _hls(sectorLook.sector_color(3200, 1.0, 0.5))
    hot = _hls(sectorLook.sector_color(20000, 1.0, 0.5))
    assert cool[0] < 0.15 and 0.5 < hot[0] < 0.75  # orange against blue
    sparse = _hls(sectorLook.sector_color(5800, 1.0, 0.0))
    dense = _hls(sectorLook.sector_color(5800, 1.0, 1.0))
    assert sparse[2] == pytest.approx(sectorLook.EMPTY_SATURATION, abs=0.01)
    assert dense[2] == pytest.approx(1.0, abs=0.01)
    assert sparse[1] == pytest.approx(sectorLook.LIGHTNESS_AT_SUN, abs=0.01)
    dim = _hls(sectorLook.sector_color(5800, 0.01, 0.5))
    bright = _hls(sectorLook.sector_color(5800, 100.0, 0.5))
    assert dim[1] < sparse[1] < bright[1]
    assert _hls(sectorLook.sector_color(5800, 1e-9, 0.5))[1] == pytest.approx(sectorLook.LIGHTNESS_RANGE[0], abs=0.01)
    assert _hls(sectorLook.sector_color(5800, 1e9, 0.5))[1] == pytest.approx(sectorLook.LIGHTNESS_RANGE[1], abs=0.01)


def test_a_sector_without_stars_has_no_color():
    assert sectorLook.sector_color(None, None, 0.0) is None
    assert sectorLook.sector_color(5800, None, 0.2) is None
