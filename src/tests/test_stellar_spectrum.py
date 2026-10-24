"""
The star's light at a planet, band by band (GEN.443, `planetgen.physics.stellar_spectrum`): the blackbody shares the research
computed, the open ultraviolet bracket and the grid sums.
"""

import math
import random

import numpy as np
import pytest

from planetgen.physics import constants, stellar_spectrum as ss, uv_cross_sections

SUN_K = 5772.0
SUN_AT_1AU_W_M2 = 1361.17


def test_the_sun_at_one_au_is_1361_watts_per_square_metre():
    assert ss.bolometric_flux_w_m2(constants.SOLAR_LUMINOSITY, 1.0) == pytest.approx(SUN_AT_1AU_W_M2, rel=1e-4)
    assert ss.bolometric_flux_w_m2(constants.SOLAR_LUMINOSITY, 2.0) == pytest.approx(SUN_AT_1AU_W_M2 / 4.0, rel=1e-4)


def test_a_blackbody_is_all_of_its_power():
    assert ss.blackbody_fraction(0.0, math.inf, SUN_K) == pytest.approx(1.0, abs=1e-5)
    assert ss.blackbody_fraction(0.0, 1e9, 300.0) == pytest.approx(1.0, abs=1e-5)
    assert ss.blackbody_fraction(400.0, 700.0, SUN_K) + ss.blackbody_fraction(700.0, math.inf, SUN_K) \
        + ss.blackbody_fraction(0.0, 400.0, SUN_K) == pytest.approx(1.0, abs=1e-5)


def test_the_blackbody_ultraviolet_of_the_sun_is_what_the_research_computed():
    # radiation-model.md section 12: 5772 K blackbody, W/m2 at 1 AU (the measured Sun, ATLAS-1: 6.85, 17.44, 86.00)
    got = [SUN_AT_1AU_W_M2 * ss.blackbody_fraction(low, high, SUN_K) for low, high in ((200, 280), (280, 315), (315, 400))]
    measured = (6.85, 17.44, 86.00)
    assert [g / m for g, m in zip(got, measured)] == pytest.approx([3.9, 1.6, 1.27], rel=0.03)


def test_the_four_micron_edge_leaks_under_one_percent_for_sun_like_and_cool_earth_blackbodies():
    # radiation-model.md R6 evidence: 0.98 percent of a 5772 K blackbody lies beyond 4 micron, 0.14 percent of 288 K below it
    assert ss.blackbody_fraction(4000.0, math.inf, SUN_K) == pytest.approx(0.0098, rel=0.01)
    assert ss.blackbody_fraction(0.0, 4000.0, 288.0) == pytest.approx(0.0014, rel=0.05)
    assert ss.blackbody_fraction(4000.0, math.inf, 3000.0) == pytest.approx(0.055, rel=0.02)


def test_a_blackbody_needs_a_temperature():
    with pytest.raises(ValueError):
        ss.blackbody_fraction(100.0, 200.0, 0.0)


def test_the_ultraviolet_shares_follow_one_family_through_the_suns_measured_values():
    sun = ss.uv_shares(0.256)
    assert sun["uv_c"] == pytest.approx(0.256)
    assert sun["uv_b"] == pytest.approx(0.625, rel=0.01)
    assert sun["uv_a"] == pytest.approx(0.787, rel=0.01)
    assert ss.uv_shares(1.0) == {"uv_c": 1.0, "uv_b": 1.0, "uv_a": 1.0}
    shares = ss.uv_shares(0.1)
    assert shares["uv_c"] < shares["uv_b"] < shares["uv_a"] < 1.0


def test_the_ultraviolet_share_is_drawn_log_uniform_in_the_research_bracket():
    rng = random.Random(3)
    draws = [ss.draw_uv_c_share(rng) for _ in range(4000)]
    assert min(draws) >= 0.1 and max(draws) <= 1.0
    assert np.median(draws) == pytest.approx(math.sqrt(0.1), rel=0.05)


def test_the_cells_of_a_spectrum_add_up_to_the_bolometric_flux_below_the_grid_top():
    edges = uv_cross_sections.standard_edges()
    cells = ss.irradiance_w_m2(edges, SUN_K, SUN_AT_1AU_W_M2, uv_c_share=1.0)
    inside = ss.blackbody_fraction(edges[0], edges[-1], SUN_K)
    assert cells.sum() == pytest.approx(SUN_AT_1AU_W_M2 * inside, rel=1e-4)
    cut = ss.irradiance_w_m2(edges, SUN_K, SUN_AT_1AU_W_M2, uv_c_share=0.256)
    assert ss.band_flux_w_m2(edges, cut, "uv_c") == pytest.approx(0.256 * ss.band_flux_w_m2(edges, cells, "uv_c"), rel=1e-9)
    assert ss.band_flux_w_m2(edges, cut, "visible") == pytest.approx(ss.band_flux_w_m2(edges, cells, "visible"), rel=1e-9)
