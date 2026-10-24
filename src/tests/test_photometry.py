# tests/test_photometry.py

"""GEN.416: a star's colour from its temperature, against the table of docs/design/composition/atmosphere-composition-parts/02-gas-optics-uv-visible.md section 2.9."""

import math

import pytest

from planetgen.physics import photometry

# T K: (x, y, 8-bit sRGB) of a blackbody, CIE 1931 2-degree tables (section 2.9).
TABLE = {
    3000: (0.4369, 0.4041, (255, 184, 109)),
    3500: (0.4053, 0.3907, (255, 199, 139)),
    4000: (0.3804, 0.3767, (255, 211, 165)),
    4500: (0.3608, 0.3635, (255, 222, 188)),
    5000: (0.3451, 0.3516, (255, 230, 208)),
    5772: (0.3266, 0.3358, (255, 241, 234)),
    6504: (0.3135, 0.3236, (255, 249, 254)),
    7000: (0.3064, 0.3165, (244, 242, 255)),
    8000: (0.2952, 0.3048, (227, 231, 255)),
    10000: (0.2806, 0.2883, (205, 217, 255)),
    15000: (0.2637, 0.2673, (181, 201, 255)),
    20000: (0.2565, 0.2576, (171, 193, 255)),
    30000: (0.2501, 0.2488, (162, 187, 255)),
}

# T K: (dx, dy) BT-Settl minus blackbody, section 2.9.
BT_SETTL = {2300: (-0.0100, -0.0222), 2500: (-0.0495, -0.0154), 3000: (-0.0208, -0.0079), 3500: (0.0026, -0.0039), 4000: (0.0180, -0.0023)}


def test_c2_is_the_second_radiation_constant():
    assert photometry.SECOND_RADIATION_CONSTANT_M_K == pytest.approx(1.438776877e-2, rel=1e-9)


@pytest.mark.parametrize("temperature", sorted(TABLE))
def test_a_blackbody_matches_the_cie_table(temperature):
    x, y, rgb8 = TABLE[temperature]
    colour = photometry.blackbody_colour(temperature)
    assert colour.x == pytest.approx(x, abs=0.001)
    assert colour.y == pytest.approx(y, abs=0.001)
    assert all(abs(a - b) <= 2 for a, b in zip(colour.rgb8, rgb8))


def test_the_sun_is_a_warm_white():
    colour = photometry.blackbody_colour(5772)
    assert colour.hex == "#fff1ea"
    assert max(colour.linear) == 1.0


def test_a_hotter_star_is_bluer_and_a_cooler_one_redder():
    colours = [photometry.star_colour(t) for t in (2500, 3500, 5000, 6500, 10000, 30000)]
    assert [c.x for c in colours] == sorted((c.x for c in colours), reverse=True)
    reds = [c.linear[0] / c.linear[2] for c in colours]
    assert reds == sorted(reds, reverse=True)


def test_the_brightest_channel_is_1_and_none_is_negative():
    for temperature in (1500, 2000, 3000, 5772, 9000, 45000):
        linear = photometry.blackbody_colour(temperature).linear
        assert max(linear) == 1.0 and min(linear) >= 0.0


def test_the_hex_is_the_8_bit_colour():
    colour = photometry.star_colour(4200)
    assert colour.hex == "#{:02x}{:02x}{:02x}".format(*colour.rgb8)


@pytest.mark.parametrize("temperature", sorted(BT_SETTL))
def test_a_cool_star_is_the_blackbody_plus_the_bt_settl_offset(temperature):
    dx, dy = BT_SETTL[temperature]
    bb = photometry.blackbody_colour(temperature)
    star = photometry.star_colour(temperature)
    assert star.x - bb.x == pytest.approx(dx, abs=2e-4)
    assert star.y - bb.y == pytest.approx(dy, abs=2e-4)


def test_a_bt_settl_star_matches_its_8_bit_table_row():
    # 3000 K: BT-Settl sRGB 255, 194, 129; 4000 K: 255, 194, 151 (section 2.9).
    assert all(abs(a - b) <= 2 for a, b in zip(photometry.star_colour(3000).rgb8, (255, 194, 129)))
    assert all(abs(a - b) <= 2 for a, b in zip(photometry.star_colour(4000).rgb8, (255, 194, 151)))


def test_the_colour_does_not_step_at_4500_k():
    below = photometry.star_colour(4499.0)
    above = photometry.star_colour(4500.0)
    assert above == photometry.blackbody_colour(4500.0)
    assert abs(below.x - above.x) < 1e-4 and abs(below.y - above.y) < 1e-4


def test_below_the_grid_the_offset_stays_at_the_last_row():
    assert photometry.cool_offset(1800.0)[:2] == (-0.0100, -0.0222)


def test_the_spread_moves_a_cool_star_by_its_half_width_and_a_hot_one_not_at_all():
    mean = photometry.star_colour(3000)
    high = photometry.star_colour(3000, spread_x=1.0, spread_y=1.0)
    assert high.x - mean.x == pytest.approx(0.02, abs=1e-3)
    assert high.y - mean.y == pytest.approx(0.01, abs=1e-3)
    assert photometry.star_colour(6000, spread_x=1.0, spread_y=-1.0) == photometry.star_colour(6000)
    half = photometry.star_colour(4250, spread_x=1.0)
    assert half.x - photometry.star_colour(4250).x == pytest.approx(0.01, abs=1e-3)


def test_a_spectrum_on_its_own_grid_gives_the_same_colour():
    grid = [380.0 + 5.0 * i for i in range(81)]
    values = [photometry.planck_relative(lam, 5772.0) for lam in grid]
    coarse = photometry.colour_from_xyz(photometry.spectrum_to_xyz(values, grid))
    fine = photometry.blackbody_colour(5772.0)
    assert coarse.x == pytest.approx(fine.x, abs=2e-4) and coarse.y == pytest.approx(fine.y, abs=2e-4)


def test_an_equal_energy_spectrum_is_the_equal_energy_white():
    xyz = photometry.spectrum_to_xyz(lambda lam: 1.0)
    x, y = photometry.chromaticity(xyz)
    assert x == pytest.approx(1 / 3, abs=0.01) and y == pytest.approx(1 / 3, abs=0.01)


def test_the_cie_curves_peak_where_the_tables_do():
    assert photometry.cie_xyz_bar(555)[1] == pytest.approx(1.0, abs=0.01)
    assert photometry.cie_xyz_bar(600)[0] == pytest.approx(1.06, abs=0.03)
    assert photometry.cie_xyz_bar(445)[2] == pytest.approx(1.78, abs=0.05)


def test_the_srgb_curve_is_the_standard_one():
    assert photometry.srgb_encode(0.0) == 0.0 and photometry.srgb_encode(1.0) == pytest.approx(1.0)
    assert photometry.srgb_encode(0.002) == pytest.approx(12.92 * 0.002)
    assert photometry.srgb_encode(0.5) == pytest.approx(0.7354, abs=1e-4)


def test_the_same_temperature_gives_the_same_colour_to_the_last_bit():
    assert photometry.star_colour(3777.7, 0.3, -0.2) == photometry.star_colour(3777.7, 0.3, -0.2)


@pytest.mark.parametrize("bad", [0.0, -5.0, math.nan, math.inf])
def test_a_temperature_must_be_positive_and_finite(bad):
    with pytest.raises(ValueError):
        photometry.star_colour(bad)


def test_a_spectrum_with_no_visible_light_raises():
    with pytest.raises(ValueError):
        photometry.spectrum_to_xyz(lambda lam: 0.0)
    with pytest.raises(ValueError):
        photometry.spectrum_to_xyz([1.0, 2.0])
