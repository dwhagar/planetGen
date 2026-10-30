# src/tests/test_atmosphere_retention.py

"""GEN.355: the retention library answers 'when is an atmosphere possible' on the 13 bodies the design table was written from."""

import pytest

from planetgen.physics import atmosphere_retention as ar

# name: (escape velocity km/s, insolation in Earth's, ratio at K = 3.85e-4 from the design note, verdict)
BODIES = {
    "Earth": (11.19, 1.0, 0.17, "dense"),
    "Venus": (10.36, 1.91, 0.43, "dense"),
    "Titan": (2.64, 0.011, 0.59, "dense"),
    "Triton": (1.45, 0.0011, 0.65, "dense"),
    "Pluto": (1.21, 0.00064, 0.78, "moderate"),
    "Ganymede": (2.74, 0.037, 1.7, "tenuous"),
    "Mars": (5.03, 0.431, 1.77, "tenuous"),
    "Io": (2.56, 0.037, 2.2, "tenuous"),
    "Callisto": (2.44, 0.037, 2.7, "trace"),
    "Europa": (2.03, 0.037, 5.7, "trace"),
    "Mercury": (4.25, 6.67, 53.0, "none"),
    "Moon": (2.38, 1.0, 82.0, "none"),
    "Ceres": (0.51, 0.13, 4990.0, "none"),
}


@pytest.mark.parametrize("name", list(BODIES))
def test_the_shoreline_ratio_and_verdict_of_the_thirteen_bodies(name):
    velocity, insolation, ratio, verdict = BODIES[name]
    found = ar.shoreline_ratio(velocity, insolation)
    assert found == pytest.approx(ratio, rel=0.05)
    assert ar.verdict_of(found) == verdict


def test_the_window_of_the_constant_holds_the_bodies_that_bound_it():
    low, high = ar.SHORELINE_K_WINDOW
    assert low < ar.SHORELINE_K < high
    assert ar.shoreline_ratio(*BODIES["Titan"][:2], k=2.26e-4) == pytest.approx(1.0, abs=0.02)      # Titan holds air: K at least this
    assert ar.shoreline_ratio(*BODIES["Ganymede"][:2], k=6.56e-4) == pytest.approx(1.0, abs=0.02)   # Ganymede holds none: K at most this


def test_the_pressure_cap_steps_down_with_the_ratio():
    assert ar.pressure_cap_pa(0.17) is None
    caps = [ar.pressure_cap_pa(r) for r in (0.8, 1.7, 5.7, 53.0)]
    assert caps == [1.0e4, 1.0e3, 1.0, 0.0]
    assert caps == sorted(caps, reverse=True)


def test_the_escape_velocity_of_earth_and_a_super_earth():
    assert ar.escape_velocity_km_s(1.0, 1.0) == pytest.approx(11.186)
    assert ar.escape_velocity_km_s(8.0, 1.8) == pytest.approx(11.186 * (8.0 / 1.8) ** 0.5)


def test_earth_holds_everything_but_hydrogen_and_helium():
    found = ar.retention(1.0, 1.0, 1.0)
    assert "H2" not in found.held and "He" not in found.held
    assert {"N2", "O2", "CO2", "H2O", "Ar"} <= set(found.held)
    assert found.verdict == "dense" and found.pressure_cap_pa is None


def test_a_body_of_the_moons_size_holds_none_of_the_gases():
    found = ar.retention(0.0123, 0.2724, 1.0)
    assert found.verdict == "none" and found.pressure_cap_pa == 0.0
    assert found.held == ()


def test_the_lightest_gas_kept_falls_with_the_escape_velocity_and_rises_with_the_exobase():
    assert ar.minimum_molar_mass_kg_mol(11.2, 1000.0) < ar.minimum_molar_mass_kg_mol(5.0, 1000.0)
    assert ar.minimum_molar_mass_kg_mol(5.0, 500.0) < ar.minimum_molar_mass_kg_mol(5.0, 1000.0)
    # the design table: 5 km/s at 1000 K keeps a molar mass of 20 and over; 10 km/s at 1000 K 5 and over
    assert ar.minimum_molar_mass_kg_mol(5.0, 1000.0) * 1000.0 == pytest.approx(20.0, rel=0.02)
    assert ar.minimum_molar_mass_kg_mol(10.0, 1000.0) * 1000.0 == pytest.approx(5.0, rel=0.02)


def test_the_cold_side_of_carbon_dioxide_nitrogen_and_water():
    assert ar.cold_side_k("CO2", 600.0) == pytest.approx(148.0, abs=6.0)       # the frost point at Mars's poles
    assert ar.cold_side_k("N2", 1.0) == pytest.approx(33.0, abs=6.0)           # Pluto's nitrogen is about 37 K
    assert ar.cold_side_k("H2O", 1.0e5) == pytest.approx(373.0, abs=3.0)
    assert ar.cold_side_k("N2", 1.0e5) < ar.cold_side_k("CO2", 1.0e5) < ar.cold_side_k("H2O", 1.0e5)


def test_the_energy_limited_loss_is_small_for_a_dense_body_and_large_for_a_light_core_near_its_star():
    sun_cumulative = 3.1e16
    assert ar.energy_limited_loss_fraction(1.0, 1.0, sun_cumulative) < 0.01
    assert ar.energy_limited_loss_fraction(2.0, 1.2, 100.0 * sun_cumulative) > ar.energy_limited_loss_fraction(2.0, 1.2, sun_cumulative)
    assert ar.energy_limited_loss_fraction(2.0, 1.2, sun_cumulative) > ar.energy_limited_loss_fraction(8.0, 1.8, sun_cumulative)


def test_an_active_star_moves_the_shoreline_by_its_cumulative_xuv():
    quiet = ar.retention(1.0, 1.0, 0.4)
    active = ar.retention(1.0, 1.0, 0.4, xuv_index=40.0, xuv_flux_earth=20.0)
    assert active.ratio == pytest.approx(40.0 * quiet.ratio)
    assert active.exobase_k > quiet.exobase_k
