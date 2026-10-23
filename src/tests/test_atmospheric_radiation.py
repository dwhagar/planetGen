"""
Atmospheric radiation band by band (GEN.443, `planetgen.physics.atmospheric_radiation`): the research's surface-UV table
(radiation-model.md section 7 and part 02 2.7), the action spectra, the layers and the limits.

The numbers of the research table need the real ozone cross sections (`datasets/ozone_uv.py`, downloaded by the setup step); those
tests skip on a machine that has not fetched them. Everything else runs on the fixture archive the tests build offline.
"""

import numpy as np
import pytest

from planetgen.datasets import fetch, ozone_uv
from planetgen.physics import air_column, atmospheric_radiation as ar, constants, stellar_spectrum as ss, uv_cross_sections as uvx

SUN_FLUX = ss.bolometric_flux_w_m2(constants.SOLAR_LUMINOSITY, 1.0)
SUN_K = 5772.0
SUN_UV_C_SHARE = 0.256        # the Sun's measured UV-C over its blackbody's (a target, radiation-model.md section 12)
EARTH_AIR_PA = {"N2": 0.78084 * 101325.0, "O2": 0.20946 * 101325.0, "Ar": 0.00934 * 101325.0, "CO2": 4e-4 * 101325.0}

_real = fetch.locate(ozone_uv.DATASET)
needs_real_ozone = pytest.mark.skipif(_real.missing, reason="the ozone cross sections are not downloaded (python3 -m planetgen.cli.datasets fetch)")


def _earth_air():
    return air_column.from_partials(EARTH_AIR_PA, 9.80665)


def _sun(column, **kwargs):
    options = dict(temperature_k=SUN_K, flux_w_m2=SUN_FLUX, uv_c_share=SUN_UV_C_SHARE, ground_albedo=0.0)
    options.update(kwargs)
    return ar.surface_radiation(column, **options)


def _direct(spectrum):
    c = spectrum.centres_nm
    return {"uvi": ar.UVI_PER_W_M2 * float(np.sum(spectrum.direct_w_m2 * ar.erythemal_weight(c))),
            "dna": float(np.sum(spectrum.direct_w_m2 * ar.setlow_weight(c))),
            "uv_b": float(np.sum(spectrum.direct_w_m2[(c >= 280) & (c < 315)]))}


# ---- the research table: overhead Sun, no aerosol, no cloud, no ground (part 02 2.7, ATLAS-1 spectrum) ----

RESEARCH_DIRECT = {      # DU: (UVI direct, Setlow-weighted direct W/m2, UV-B 280 to 315 direct W/m2)
    0: (133.6, 1.23, 5.24), 100: (23.4, 0.0339, 2.47), 200: (11.5, 9.7e-3, 1.70), 300: (7.27, 4.4e-3, 1.27), 500: (4.11, 1.6e-3, 0.776),
}
RESEARCH_GLOBAL_UVI = {100: 48.9, 200: 22.8, 300: 13.9, 500: 7.4}      # the research's crude global (plus or minus 50 percent)


@needs_real_ozone
@pytest.mark.parametrize("du", sorted(RESEARCH_DIRECT))
def test_the_direct_uv_follows_the_research_table_on_a_star_with_the_suns_ultraviolet_shares(du):
    uvi, dna, uv_b = RESEARCH_DIRECT[du]
    got = _direct(_sun(_earth_air(), ozone_du=du))
    assert got["uvi"] == pytest.approx(uvi, rel=0.12)
    assert got["dna"] == pytest.approx(dna, rel=0.15)
    assert got["uv_b"] == pytest.approx(uv_b, rel=0.15)


@needs_real_ozone
@pytest.mark.parametrize("du", sorted(RESEARCH_GLOBAL_UVI))
def test_the_global_uv_index_is_inside_the_crude_global_estimate(du):
    assert _sun(_earth_air(), ozone_du=du).uv_index() == pytest.approx(RESEARCH_GLOBAL_UVI[du], rel=0.5)


@needs_real_ozone
def test_the_uv_index_falls_as_the_ozone_column_to_about_minus_one_power():
    # the research: direct UVI falls as DU^-1.08 from 100 to 500 DU
    low, high = (_direct(_sun(_earth_air(), ozone_du=du))["uvi"] for du in (100, 500))
    exponent = -np.log(high / low) / np.log(5.0)
    assert exponent == pytest.approx(1.08, abs=0.15)


@needs_real_ozone
def test_a_planet_with_oxygen_and_no_ozone_has_lethal_uv_c():
    # part 02 2.7: the Setlow dose with no ozone is about 280 times that of 300 DU
    ratio = _direct(_sun(_earth_air(), ozone_du=0.0))["dna"] / _direct(_sun(_earth_air(), ozone_du=300.0))["dna"]
    assert ratio == pytest.approx(280.0, rel=0.25)
    assert _sun(_earth_air(), ozone_du=0.0).band_ground_w_m2("uv_c") > 0.1
    assert _sun(_earth_air(), ozone_du=300.0).band_ground_w_m2("uv_c") < 1e-9


@needs_real_ozone
def test_a_star_with_no_atmosphere_between_has_the_uv_index_of_its_light_alone():
    sun = _sun(_earth_air(), ozone_du=300.0)
    assert sun.top_uv_index() == pytest.approx(666.0, rel=0.2)      # ATLAS-1, part 02 2.7
    assert sun.top_dna_dose_w_m2() == pytest.approx(8.66, rel=0.15)


# ---- the action spectra ----

def test_the_erythemal_action_spectrum_is_the_cie_one():
    nm = np.array([249.0, 260.0, 298.0, 300.0, 305.0, 310.0, 320.0, 328.0, 330.0, 350.0, 380.0, 401.0])
    expected = [0.0, 1.0, 1.0, 0.649, 0.220, 0.0745, 0.0086, 10 ** (0.094 * (298 - 328.0)), 1.41e-3, 7.1e-4, 10 ** (0.015 * (140 - 380.0)), 0.0]      # part 02 2.6 prints 3.1e-4 at 380 nm; its own formula gives 2.5e-4
    assert ar.erythemal_weight(nm) == pytest.approx(expected, rel=0.01)


def test_the_setlow_spectrum_is_the_tabulated_one_per_quantum_and_times_lambda_over_300_per_energy():
    assert ar.setlow_weight(np.array([300.0]))[0] == pytest.approx(0.0333)
    assert ar.setlow_weight(np.array([260.0]))[0] == pytest.approx(1.12 * 260.0 / 300.0)
    assert ar.setlow_weight(np.array([200.0]))[0] == pytest.approx(1.07 * 200.0 / 300.0)         # held below the tabulated range
    assert ar.setlow_weight(np.array([450.0]))[0] == 0.0
    nm = np.linspace(256.0, 360.0, 200)
    assert np.all(np.diff(np.log(ar.setlow_weight(nm)[nm > 270.0])) < 0.0)


# ---- the layers ----

def test_the_ozone_layer_sits_at_the_chapman_peak_for_the_columns_oxygen():
    earth = _earth_air()
    assert ar.ozone_peak_pressure_pa(earth) == pytest.approx(6300.0 / (0.20946 * 101325.0 / 21200.0), rel=1e-6)
    tenth = air_column.from_partials({"N2": 0.8 * 101325.0, "O2": 0.02095 * 101325.0}, 9.80665)
    assert ar.ozone_peak_pressure_pa(tenth) == pytest.approx(6300.0 / (tenth.partials_pa["O2"] / 21200.0), rel=1e-6)
    assert ar.ozone_peak_pressure_pa(tenth) == pytest.approx(63000.0, rel=0.01)       # part 04 3.3: 630 hPa at 0.1 PAL
    no_oxygen = air_column.from_partials({"N2": 1e5}, 9.8)
    assert ar.ozone_peak_pressure_pa(no_oxygen) == pytest.approx(0.9e5)
    thin = air_column.from_partials({"CO2": 600.0, "O2": 0.1}, 3.7)
    assert ar.ozone_peak_pressure_pa(thin) <= 0.9 * 600.1


def test_the_sun_angle_must_be_inside_the_range_the_two_stream_is_good_for():
    for mu0 in (0.1, 1.5):
        with pytest.raises(ValueError):
            _sun(_earth_air(), mu0=mu0)


def test_a_dense_carbon_dioxide_world_lets_no_ultraviolet_through():
    with fetch.disregard_missing_data():
        venus = air_column.from_partials({"CO2": 0.965 * 9.2e6, "N2": 0.035 * 9.2e6, "SO2": 1.5e-4 * 9.2e6, "H2S": 1e-4 * 9.2e6}, 8.87)
        spectrum = ar.surface_radiation(venus, temperature_k=SUN_K, flux_w_m2=SUN_FLUX / 0.723 ** 2, ground_albedo=0.1)
    assert spectrum.uv_index() < 1e-6
    assert spectrum.band_transmission("uv_c") == 0.0
    assert spectrum.band_transmission("visible") < 0.5


def test_the_bands_of_a_clear_sky_are_ordered_by_how_much_air_scatters_and_absorbs():
    with fetch.disregard_missing_data():
        fractions = _sun(_earth_air(), ozone_du=300.0).band_fractions()
    assert fractions["uv_c"] < fractions["uv_b"] < fractions["uv_a"] < fractions["visible"] < fractions["near_infrared"] <= 1.0
    assert fractions["visible"] == pytest.approx(0.92, abs=0.06)           # 0.097 optical depth at 550 nm


def test_the_bands_sum_to_the_whole_ground_flux_inside_the_grid():
    with fetch.disregard_missing_data():
        spectrum = _sun(_earth_air(), ozone_du=300.0)
    parts = sum(spectrum.band_ground_w_m2(band) for band in ("uv_c", "uv_b", "uv_a", "visible", "near_infrared"))
    assert parts == pytest.approx(spectrum.ground_w_m2[spectrum.centres_nm >= 100.0].sum(), rel=1e-9)
    assert spectrum.band_top_w_m2("visible") == pytest.approx(SUN_FLUX * ss.blackbody_fraction(400.0, 700.0, SUN_K), rel=1e-3)


def test_cloud_cuts_the_light_in_proportion_to_its_cover_and_depth():
    with fetch.disregard_missing_data():
        clear = _sun(_earth_air(), ozone_du=300.0)
        edges = clear.edges_nm
        cloud = ar.grey_layer(20.0, ar.CLOUD_OMEGA, ar.CLOUD_ASYMMETRY, uvx.centres(edges))
        overcast = _sun(_earth_air(), ozone_du=300.0, low_aerosol=cloud, cloud_cover=1.0)
        half = _sun(_earth_air(), ozone_du=300.0, low_aerosol=cloud, cloud_cover=0.5)
    assert overcast.uv_index() < 0.5 * clear.uv_index()
    assert half.uv_index() == pytest.approx(0.5 * (clear.uv_index() + overcast.uv_index()), rel=1e-6)
    assert overcast.band_ground_w_m2("visible") < 0.5 * clear.band_ground_w_m2("visible")
    assert overcast.direct_w_m2.sum() < 0.01 * clear.direct_w_m2.sum()


def test_bright_ground_raises_the_diffuse_light_and_a_dark_absorbing_haze_lowers_it():
    with fetch.disregard_missing_data():
        dark = _sun(_earth_air(), ozone_du=300.0, ground_albedo=0.0)
        snow = _sun(_earth_air(), ozone_du=300.0, ground_albedo=0.8)
        edges = dark.edges_nm
        haze = ar.grey_layer(0.5, 0.6, 0.6, uvx.centres(edges), exponent=1.5)
        hazy = _sun(_earth_air(), ozone_du=300.0, ground_albedo=0.0, haze=haze)
    assert snow.band_ground_w_m2("uv_a") > dark.band_ground_w_m2("uv_a")
    assert hazy.band_ground_w_m2("visible") < dark.band_ground_w_m2("visible")
    assert hazy.band_ground_w_m2("uv_a") < dark.band_ground_w_m2("uv_a")


def test_a_lower_sun_means_less_light_on_a_horizontal_surface_and_a_longer_path():
    with fetch.disregard_missing_data():
        overhead = _sun(_earth_air(), ozone_du=300.0, mu0=1.0)
        low = _sun(_earth_air(), ozone_du=300.0, mu0=0.4)
    assert low.band_ground_w_m2("visible") < 0.5 * overhead.band_ground_w_m2("visible")
    assert low.band_transmission("uv_b") < overhead.band_transmission("uv_b")


def test_without_the_download_the_kernel_still_runs_on_the_gas_bins_when_told_to():
    with fetch.disregard_missing_data():
        spectrum = _sun(_earth_air(), ozone_du=300.0)
    assert spectrum.band_ground_w_m2("uv_b") >= 0.0 and spectrum.uv_index() >= 0.0
