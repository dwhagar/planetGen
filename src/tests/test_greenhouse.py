"""
Greenhouse closure (GEN.443, `planetgen.physics.greenhouse`): the Robinson and Catling grey atmosphere against the
targets of docs/design/composition/radiation-model.md section 13. Earth, Venus and the class multipliers are test
targets, never inputs of the code.
"""

import pytest

from planetgen.physics import greenhouse as gh

EARTH = dict(absorbed=240.0, gamma=1.4, alpha=0.6, n=2)
VENUS = dict(absorbed=160.0, gamma=1.3, alpha=0.8, n=1)


def _implied(temperature, world):
    return gh.tau0_implied(temperature, world["absorbed"], world["gamma"], world["alpha"], world["n"])


def test_earth_needs_the_published_optical_depth():
    tau0 = _implied(288.0, EARTH)
    assert tau0 == pytest.approx(1.82, rel=0.03)


def test_earth_greenhouse_factor_and_boundary():
    effective = gh.effective_temperature_k(EARTH["absorbed"])
    assert gh.greenhouse_factor(288.0, effective) == pytest.approx(0.625, rel=0.01)
    result = gh.check_stored_multiplier(288.0, EARTH["absorbed"], EARTH["gamma"], EARTH["alpha"], EARTH["n"])
    assert result["tau_rc"] == pytest.approx(0.0863, rel=0.02)


def test_venus_like_world_needs_about_four_hundred():
    effective = gh.effective_temperature_k(VENUS["absorbed"])
    assert gh.greenhouse_factor(730.0, effective) == pytest.approx(99.6, rel=0.01)
    assert _implied(730.0, VENUS) == pytest.approx(400.0, rel=0.15)


def test_the_closure_and_its_inverse_agree():
    for tau0 in (0.01, 0.3, 1.8, 40.0, 360.0):
        temperature, _ = gh.surface_temperature_k(tau0, EARTH["absorbed"], EARTH["gamma"], EARTH["alpha"], EARTH["n"])
        assert _implied(temperature, EARTH) == pytest.approx(tau0, rel=1e-6)


def test_more_opacity_is_warmer():
    temperatures = [gh.surface_temperature_k(t, EARTH["absorbed"], EARTH["gamma"], EARTH["alpha"], EARTH["n"])[0]
                    for t in (0.1, 0.5, 1.8, 10.0, 100.0)]
    assert temperatures == sorted(temperatures)
    assert temperatures[0] > gh.effective_temperature_k(EARTH["absorbed"])


def test_no_opacity_is_the_effective_temperature():
    effective = gh.effective_temperature_k(EARTH["absorbed"])
    assert gh.surface_temperature_k(0.0, EARTH["absorbed"], EARTH["gamma"], EARTH["alpha"], EARTH["n"]) == (effective, None)


def test_an_atmosphere_thinner_than_the_search_floor_is_purely_radiative():
    temperature, boundary = gh.surface_temperature_k(5e-7, EARTH["absorbed"], EARTH["gamma"], EARTH["alpha"], EARTH["n"])
    assert boundary is None
    assert temperature ** 4 == pytest.approx(
        EARTH["absorbed"] / gh.constants.STEFAN_BOLTZMANN_CONSTANT * (1 + gh.DIFFUSIVITY * 5e-7 / 2), rel=1e-9)


def test_a_carbon_dioxide_doubling_warms_earth_by_a_degree():
    """A forcing of 3.7 W/m2 needs tau0 raised by about 5 percent and warms the state by about 1.1 K (section 4)."""
    tau0 = _implied(288.0, EARTH)
    warmed = gh.surface_temperature_k(tau0 * 1.05, EARTH["absorbed"], EARTH["gamma"], EARTH["alpha"], EARTH["n"])[0]
    assert warmed - 288.0 == pytest.approx(1.12, rel=0.2)


def test_the_check_reports_what_a_stored_multiplier_needs():
    result = gh.check_stored_multiplier(288.0, EARTH["absorbed"], EARTH["gamma"], EARTH["alpha"], EARTH["n"],
                                        stored_multiplier=1.168)
    assert result["stored_surface_temperature_k"] == pytest.approx(309.5, abs=0.2)
    assert result["stored_tau0"] == pytest.approx(4.3, rel=0.02)


def test_a_surface_no_warmer_than_the_effective_temperature_is_refused():
    effective = gh.effective_temperature_k(EARTH["absorbed"])
    with pytest.raises(ValueError):
        gh.tau0_implied(effective, EARTH["absorbed"], EARTH["gamma"], EARTH["alpha"], EARTH["n"])


def test_a_temperature_beyond_the_search_range_is_refused():
    with pytest.raises(ValueError):
        gh.tau0_implied(1e6, EARTH["absorbed"], EARTH["gamma"], EARTH["alpha"], EARTH["n"])


def test_absorbed_flux_of_the_sun_at_one_au():
    assert gh.absorbed_flux_w_m2(1361.17, 0.30) == pytest.approx(238.2, rel=0.001)
