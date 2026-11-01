"""
GEN.418: the atmosphere colour of a world with a ground (`planetgen.physics.atmosphere_colour`), against the worked example of
docs/design/composition/atmosphere-colour.md section 14.1 and the targets of section 15.

Solar System rows (the Earth-like test air, the thin carbon dioxide and the dense carbon dioxide cases) live here only; the module
reads none of them. The ozone's Chappuis band is a table in the repository, so no test needs a download.
"""

import math

import numpy as np
import pytest

from planetgen.physics import air_column, atmosphere_colour as ac, photometry

EARTH_AIR_PA = {"N2": 0.78084 * 101325.0, "O2": 0.209476 * 101325.0, "Ar": 0.00934 * 101325.0, "CO2": 4e-4 * 101325.0}
SUN_K = 5772.0
DUST = ac.Aerosol(0.9, ((440, 0.772), (460, 0.788), (500, 0.831), (600, 0.945), (700, 0.971), (800, 0.968)), 0.673)


def _earth(scale=1.0):
    return air_column.from_partials({gas: p * scale for gas, p in EARTH_AIR_PA.items()}, 9.80665)


def _thin_co2():
    return air_column.from_partials({"CO2": 600.0}, 3.7)


def _dense_co2():
    return air_column.from_partials({"CO2": 9.2e6}, 8.87)


def _channels(hex_colour):
    return [int(hex_colour[i:i + 2], 16) for i in (1, 3, 5)]


def test_earth_like_air_under_a_sun_like_star_is_the_blue_of_the_worked_example():
    result = ac.atmosphere_colour(_earth(), temperature_k=SUN_K)
    assert result.hex == "#99bfff"
    assert (result.colour.x, result.colour.y) == pytest.approx((0.2452, 0.2521), abs=0.001)
    assert result.word == "pale blue"
    assert result.tau_rayleigh_550 == pytest.approx(0.09721, rel=0.02)         # [R0575: Table 4]
    assert result.regime == "thin"
    assert not result.aerosol_dominated and result.aerosol_shift == 0.0
    assert result.label == "pale blue"


def test_the_same_air_under_a_red_dwarf_is_cream_not_blue():
    result = ac.atmosphere_colour(_earth(), temperature_k=3300.0)
    assert (result.colour.x, result.colour.y) == pytest.approx((0.3243, 0.3441), abs=0.006)   # a blackbody there; BT-Settl is 0.002 to 0.003 off
    assert result.word == "cream"
    assert _channels(result.hex)[0] == 255 and min(_channels(result.hex)) >= 0xe8


def test_thin_carbon_dioxide_is_a_blue_hue_at_a_thirtieth_of_the_brightness():
    thin = ac.atmosphere_colour(_thin_co2(), temperature_k=SUN_K)
    earth = ac.atmosphere_colour(_earth(), temperature_k=SUN_K)
    assert thin.hex == "#8bb1ff"
    assert thin.tau_550 == pytest.approx(0.0026, abs=0.0002)
    assert thin.luminance_ratio / earth.luminance_ratio == pytest.approx(0.03, abs=0.01)


def test_dust_turns_the_thin_blue_sky_warm():
    clear = ac.atmosphere_colour(_thin_co2(), temperature_k=SUN_K)
    dusty = ac.atmosphere_colour(_thin_co2(), temperature_k=SUN_K, dust=DUST)
    assert dusty.colour.x > clear.colour.x + 0.08
    assert dusty.word in ("peach", "cream", "orange")
    assert dusty.aerosol_dominated and dusty.label.endswith("(aerosol-dominated)")
    assert dusty.luminance_ratio > 20.0 * clear.luminance_ratio


def test_dense_carbon_dioxide_is_orange_never_the_clipped_red_of_single_scatter():
    result = ac.atmosphere_colour(_dense_co2(), temperature_k=SUN_K)
    assert result.regime == "thick"
    assert result.tau_550 == pytest.approx(16.5, rel=0.03)
    assert (result.colour.x, result.colour.y) == pytest.approx((0.4259, 0.3972), abs=0.002)
    assert result.hex == "#ffbb79"
    assert result.word == "orange"
    assert result.sun_luminance_ratio < 1e-4          # the direct beam is dead; the disc is a faint deep red


def test_a_megapascal_of_nitrogen_is_white_cyan_with_an_orange_sun():
    air = air_column.from_partials({"N2": 1e6}, 9.8)
    result = ac.atmosphere_colour(air, temperature_k=SUN_K)
    r, g, b = _channels(result.hex)
    assert g == 255 and b >= 0xf0 and r >= 0xd8       # white tinted toward cyan
    assert result.word == "near white"
    sun_r, sun_g, sun_b = _channels(result.sun.hex)
    assert sun_r == 255 and sun_b < sun_g < sun_r      # orange
    assert 0.1 < result.sun_luminance_ratio < 0.5


def test_the_colour_does_not_step_between_thin_and_thick_air():
    previous = None
    for scale in np.geomspace(0.3, 40.0, 400):
        result = ac.atmosphere_colour(_earth(scale), temperature_k=SUN_K, ground_albedo=0.2)
        channels = _channels(result.hex)
        if previous is not None:
            assert max(abs(a - b) for a, b in zip(channels, previous[0])) <= 3
            assert abs(math.log(result.luminance_ratio / previous[1])) < 0.03
        previous = (channels, result.luminance_ratio)


def test_the_regime_follows_the_depth():
    regimes = [ac.atmosphere_colour(_earth(scale), temperature_k=SUN_K).regime for scale in (0.5, 2.5, 20.0)]
    assert regimes == ["thin", "blend", "thick"]


def test_a_brighter_ground_brightens_the_sky_and_a_red_ground_warms_it():
    black = ac.atmosphere_colour(_earth(8.0), temperature_k=SUN_K)
    white = ac.atmosphere_colour(_earth(8.0), temperature_k=SUN_K, ground_albedo=0.6)
    red = ac.atmosphere_colour(_earth(8.0), temperature_k=SUN_K, ground_albedo=np.where(ac.WAVELENGTHS_NM > 580.0, 0.6, 0.05))
    assert white.luminance_ratio > 1.2 * black.luminance_ratio
    assert red.colour.x > white.colour.x + 0.01


def test_cloud_whitens_and_brightens_by_its_cover():
    cloud = ac.Aerosol(10.0, 1.0, 0.85)
    clear = ac.atmosphere_colour(_earth(), temperature_k=SUN_K)
    none = ac.atmosphere_colour(_earth(), temperature_k=SUN_K, cloud=cloud, cloud_cover=0.0)
    half = ac.atmosphere_colour(_earth(), temperature_k=SUN_K, cloud=cloud, cloud_cover=0.5)
    full = ac.atmosphere_colour(_earth(), temperature_k=SUN_K, cloud=cloud, cloud_cover=1.0)
    assert none.hex == clear.hex and none.luminance_ratio == clear.luminance_ratio
    assert clear.luminance_ratio < half.luminance_ratio < full.luminance_ratio
    assert full.word in ("white", "near white", "cream") and full.colour.x > clear.colour.x
    assert full.aerosol_dominated


def test_ozone_takes_a_few_percent_of_the_light_and_shifts_the_hue_little():
    ozone = ac.atmosphere_colour(_earth(), temperature_k=SUN_K, ozone_du=300.0)
    bare = ac.atmosphere_colour(_earth(), temperature_k=SUN_K)
    drop = 1.0 - ozone.luminance_ratio / bare.luminance_ratio
    assert 0.0 < drop < 0.15
    assert abs(ozone.colour.x - bare.colour.x) < 0.02 and abs(ozone.colour.y - bare.colour.y) < 0.02


def test_ozone_matters_more_when_the_sun_is_low_and_the_table_is_the_chappuis_band():
    table = ac.chappuis_cross_section_cm2()
    assert ac.WAVELENGTHS_NM[int(np.argmax(table))] == pytest.approx(600.0, abs=10.0)
    assert table.max() == pytest.approx(5.1e-21, rel=0.03)
    low = ac.atmosphere_colour(_earth(), temperature_k=SUN_K, mu0=0.2)
    low_ozone = ac.atmosphere_colour(_earth(), temperature_k=SUN_K, mu0=0.2, ozone_du=300.0)
    assert 0.05 < 1.0 - low_ozone.luminance_ratio / low.luminance_ratio < 0.2


def test_the_chappuis_band_is_weaker_in_the_cold_at_its_blue_edge_and_unchanged_at_its_peak():
    """GEN.748: the ratio of [R0599] at 220 K to 298 K is 0.64 to 0.8 between 410 and 450 nm, within 2 percent of 1 near
    the 600 nm peak, and the table at 298 K is the Brion-team one."""
    cold, warm, default = ac.chappuis_cross_section_cm2(220.0), ac.chappuis_cross_section_cm2(298.0), ac.chappuis_cross_section_cm2()
    nm = ac.WAVELENGTHS_NM
    blue = (nm >= 410.0) & (nm <= 440.0)
    assert np.all(cold[blue] / warm[blue] < 0.85) and np.all(cold[blue] / warm[blue] > 0.55)
    peak = (nm >= 580.0) & (nm <= 620.0)
    assert np.all(np.abs(cold[peak] / warm[peak] - 1.0) < 0.03)
    assert np.allclose(default, ac.chappuis_cross_section_cm2(220.0))
    assert np.allclose(ac.chappuis_cross_section_cm2(150.0), cold) and np.allclose(ac.chappuis_cross_section_cm2(400.0), warm)
    middle = ac.chappuis_cross_section_cm2(230.0)
    assert np.allclose(middle, 0.5 * (cold + ac.chappuis_cross_section_cm2(240.0)))
    assert warm.max() == pytest.approx(5.1e-21, rel=0.03)


def test_a_colder_ozone_layer_lets_a_little_more_blue_through_the_air_at_a_low_sun():
    warm = ac.atmosphere_colour(_earth(), temperature_k=SUN_K, mu0=0.2, ozone_du=300.0, ozone_temperature_k=298.0)
    cold = ac.atmosphere_colour(_earth(), temperature_k=SUN_K, mu0=0.2, ozone_du=300.0, ozone_temperature_k=220.0)
    assert cold.luminance_ratio >= warm.luminance_ratio
    assert abs(cold.colour.x - warm.colour.x) < 0.01 and abs(cold.colour.y - warm.colour.y) < 0.01


def test_a_cool_star_has_the_bt_settl_offset_and_no_step_at_4500_kelvin():
    for temperature in (2500.0, 3000.0, 3500.0, 4000.0):
        spectrum = ac.star_spectrum(temperature)
        reference = photometry.star_colour(temperature)
        assert photometry.chromaticity(ac._xyz(spectrum)) == pytest.approx((reference.x, reference.y), abs=0.001)
    below = ac.atmosphere_colour(_earth(), temperature_k=4499.0)
    above = ac.atmosphere_colour(_earth(), temperature_k=4501.0)
    assert max(abs(a - b) for a, b in zip(_channels(below.hex), _channels(above.hex))) <= 2


def test_the_seeded_spread_of_a_cool_star_moves_its_air_and_not_a_hot_one():
    low = ac.atmosphere_colour(_earth(), temperature_k=3000.0, spread_x=-1.0, spread_y=-1.0)
    high = ac.atmosphere_colour(_earth(), temperature_k=3000.0, spread_x=1.0, spread_y=1.0)
    assert high.colour.x - low.colour.x > 0.03
    hot_low = ac.atmosphere_colour(_earth(), temperature_k=SUN_K, spread_x=-1.0)
    hot_high = ac.atmosphere_colour(_earth(), temperature_k=SUN_K, spread_x=1.0)
    assert hot_low.hex == hot_high.hex


def test_the_sun_at_the_zenith_is_the_limit_of_the_sun_near_it():
    overhead = ac.atmosphere_colour(_earth(), temperature_k=SUN_K, mu0=1.0)
    near = ac.atmosphere_colour(_earth(), temperature_k=SUN_K, mu0=1.0 - 1e-5)
    assert overhead.hex == near.hex
    assert overhead.luminance_ratio == pytest.approx(near.luminance_ratio, rel=1e-3)


def test_every_air_gives_a_finite_colour_and_repeats():
    airs = [{"N2": 1.0}, {"CO2": 1.0}, {"H2": 0.86, "He": 0.14}, {"CH4": 0.05, "N2": 0.95}, {"H2O": 1.0}]
    for mix in airs:
        for pressure in (1e-1, 1e3, 1e5, 1e7, 1e9):
            column = air_column.from_partials({gas: share * pressure for gas, share in mix.items()}, 10.0)
            for temperature in (2400.0, 5772.0, 12000.0):
                for mu0 in (0.2, 0.7):
                    for aerosol in (None, ac.Aerosol(3.0, 0.9, 0.5, 1.0)):
                        result = ac.atmosphere_colour(column, temperature_k=temperature, mu0=mu0, haze=aerosol, ground_albedo=0.3)
                        assert all(math.isfinite(v) for v in (result.colour.x, result.colour.y, result.luminance_ratio, result.tau_550))
                        assert len(result.hex) == 7 and 0.0 <= result.luminance_ratio <= 1.0 and result.word
    first = ac.atmosphere_colour(_earth(), temperature_k=SUN_K, haze=ac.Aerosol(0.5, 0.9, 0.5, 1.5))
    second = ac.atmosphere_colour(_earth(), temperature_k=SUN_K, haze=ac.Aerosol(0.5, 0.9, 0.5, 1.5))
    assert first == second


def test_wrong_arguments_are_refused():
    with pytest.raises(ValueError):
        ac.atmosphere_colour(_earth(), temperature_k=SUN_K, mu0=0.1)
    with pytest.raises(ValueError):
        ac.Aerosol(-1.0)
    with pytest.raises(ValueError):
        ac.Aerosol(1.0, omega=1.5)
    with pytest.raises(ValueError):
        ac.Aerosol(1.0, g=1.0)


def test_the_words():
    assert ac.colour_word((255, 255, 255)) == "white"
    assert ac.colour_word((153, 191, 255)) == "pale blue"
    assert ac.colour_word((255, 251, 236)) == "cream"
    assert ac.colour_word((240, 255, 255)) == "near white"
    assert ac.colour_word((255, 187, 121)) == "orange"
    assert ac.colour_word((255, 224, 192)) == "peach"
    assert ac.colour_word((120, 140, 255)) == "blue"
    assert ac.colour_word((255, 20, 20)) == "deep red"


@pytest.mark.parametrize("pressure_pa", [1e11, 1e12, 1e14])
def test_air_that_lets_no_light_down_still_has_a_colour(pressure_pa):
    """GEN.844: a dense steam world under a hot star had no visible light at the ground and raised; it takes the colour of
    the light the air sends back out."""
    column = air_column.from_partials({"H2O": 0.95 * pressure_pa, "H2S": 0.03 * pressure_pa, "CH4": 0.02 * pressure_pa}, 9.8)
    result = ac.atmosphere_colour(column, temperature_k=6700.0, ground_albedo=ac.STAND_IN_GROUND_ALBEDO)
    assert len(result.colour.hex) == 7 and result.colour.hex != "#000000"


def test_air_that_returns_no_light_either_is_black_on_the_stars_chromaticity():
    zero = np.zeros_like(ac.WAVELENGTHS_NM)
    star = ac.star_spectrum(5772.0)
    sky = ac._sky_colour(star, zero, zero)
    assert sky.rgb8 == (0, 0, 0) and sky.hex == "#000000"
    assert (sky.x, sky.y) == pytest.approx(photometry.chromaticity(ac._xyz(star)))


def test_an_air_that_leaves_no_light_at_the_zenith_has_no_colour_to_show(monkeypatch):
    """A steam world of 200 kilobars (GEN.646 makes them) can take all the light: the model then raises, and a body shows no colour."""
    import types

    def dark(*args, **kwargs):
        raise ValueError("the spectrum has no light between 380 and 780 nm")

    monkeypatch.setattr(ac, "atmosphere_colour", dark)
    monkeypatch.setattr(air_column, "column_of", lambda body: _earth())
    body = types.SimpleNamespace(body_type="t", star=types.SimpleNamespace(temperature=SUN_K), ozone_column_du=0.0)
    assert ac.body_atmosphere_colour(body) is None
