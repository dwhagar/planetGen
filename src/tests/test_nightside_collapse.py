# tests/test_nightside_collapse.py

"""GEN.1044: whether the air of a tidally locked planet freezes onto its night side (`physics/nightside_collapse.py`). The numbers checked are the ones [R4767] prints; Mars is a test point, never an input."""

import pytest

from planetgen.physics import nightside_collapse as nc

EARTH_FLUX, EARTH_G = 1366.0, 9.8
MARS_FLUX_W_M2, MARS_G = 590.0 * 0.43, 3.71   # Mars's flux in the units of the fit and its gravity


def test_the_fit_gives_the_paper_s_critical_pressure_at_earth_s_flux():
    pressure, extrapolated = nc.fit_critical_pressure_pa(EARTH_FLUX)
    assert pressure == pytest.approx(0.14e5, rel=0.05)   # the abstract: 0.14 bar
    assert not extrapolated


def test_a_heavier_planet_collapses_at_a_higher_pressure_and_a_brighter_star_at_a_lower_one():
    assert nc.fit_critical_pressure_pa(EARTH_FLUX, 10.0)[0] > nc.fit_critical_pressure_pa(EARTH_FLUX, 1.0)[0]
    assert nc.fit_critical_pressure_pa(EARTH_FLUX, 3.0)[0] == pytest.approx(
        (nc.fit_critical_pressure_pa(EARTH_FLUX, 1.0)[0] * nc.fit_critical_pressure_pa(EARTH_FLUX, 10.0)[0]) ** 0.5, rel=0.3)
    assert nc.fit_critical_pressure_pa(3.0 * EARTH_FLUX)[0] < nc.fit_critical_pressure_pa(EARTH_FLUX)[0] < nc.fit_critical_pressure_pa(0.4 * EARTH_FLUX)[0]


def test_outside_what_the_paper_ran_the_fit_is_held_at_its_edge_and_marked():
    low, high = nc.FIT_FLUX_RANGE_W_M2
    assert nc.fit_critical_pressure_pa(0.5 * low) == (nc.fit_critical_pressure_pa(low)[0], True)
    assert nc.fit_critical_pressure_pa(2.0 * high)[1] is True
    assert nc.fit_critical_pressure_pa(EARTH_FLUX, 100.0)[1] is True and nc.fit_critical_pressure_pa(EARTH_FLUX, 0.1)[1] is True


def test_a_brighter_surface_collapses_as_a_dimmer_star_would():
    assert nc.fit_critical_pressure_pa(EARTH_FLUX, albedo=0.5)[0] == pytest.approx(nc.fit_critical_pressure_pa(EARTH_FLUX * 0.5 / 0.8)[0])
    assert nc.fit_critical_pressure_pa(EARTH_FLUX, albedo=0.5)[0] > nc.fit_critical_pressure_pa(EARTH_FLUX, albedo=0.2)[0]


def test_the_night_side_of_a_tenth_of_a_bar_is_the_paper_s_125_k():
    assert nc.night_side_temperature_k(1.0e4, EARTH_G, EARTH_FLUX, 0.2) == pytest.approx(125.2, rel=0.02)


def test_the_night_side_warms_with_the_air_and_stops_at_the_equilibrium_temperature():
    values = [nc.night_side_temperature_k(p, EARTH_G, EARTH_FLUX) for p in (1.0e3, 1.0e4, 5.0e4, 2.0e5, 1.0e7)]
    assert values == sorted(values)
    equilibrium = (0.8 * EARTH_FLUX / 4.0 / 5.670374419e-8) ** 0.25
    assert values[-1] == pytest.approx(equilibrium, rel=1e-3) and values[-2] == values[-1]
    assert nc.night_side_temperature_k(0.0, EARTH_G, EARTH_FLUX) == 0.0


def test_an_earth_like_bar_of_co2_keeps_all_its_gas_and_a_thin_one_loses_it():
    assert nc.collapse("CO2", 1.0e5, EARTH_G, EARTH_FLUX).retained_fraction == 1.0
    thin = nc.collapse("CO2", 5.0e3, EARTH_G, EARTH_FLUX)
    assert thin.collapsed and thin.retained_fraction < 0.05
    assert nc.collapse("CO2", 1.0e4, EARTH_G, EARTH_FLUX).collapsed       # under the paper's 0.1 bar
    assert nc.collapse("CO2", 1.0e4, EARTH_G, EARTH_FLUX).critical_pressure_pa == pytest.approx(0.14e5, rel=0.05)


def test_the_retained_share_rises_without_a_jump_to_all_of_it_at_the_critical_pressure():
    critical = nc.fit_critical_pressure_pa(EARTH_FLUX)[0]
    shares = [nc.collapse("CO2", critical * f, EARTH_G, EARTH_FLUX).retained_fraction for f in (0.05, 0.3, 0.7, 0.95, 0.999, 1.0, 3.0)]
    assert shares == sorted(shares)
    assert shares[-3] > 0.99 and shares[-1] == 1.0 and shares[0] < 1e-3


def test_a_mars_like_air_on_a_locked_world_is_lost_to_the_night_side():
    result = nc.collapse("CO2", 600.0, MARS_G, MARS_FLUX_W_M2, mass_earth=0.107)
    assert result.collapsed and result.retained_fraction < 1e-6
    assert result.extrapolated      # Mars's flux is under the lowest the paper ran


def test_n2_holds_at_earth_s_flux_and_freezes_on_a_much_colder_world():
    assert nc.collapse("N2", 1.0e3, EARTH_G, EARTH_FLUX).retained_fraction == 1.0
    cold = nc.collapse("N2", 1.0e3, EARTH_G, 30.0)
    assert cold.collapsed and cold.source == "radiative-balance"
    assert nc.collapse("N2", 1.0e5, EARTH_G, 30.0).retained_fraction == 1.0


def test_no_air_is_nothing_lost_and_bad_inputs_are_refused():
    assert nc.collapse("CO2", 0.0, EARTH_G, EARTH_FLUX).retained_fraction == 1.0
    with pytest.raises(ValueError):
        nc.collapse("CO2", -1.0, EARTH_G, EARTH_FLUX)
    with pytest.raises(ValueError):
        nc.fit_critical_pressure_pa(0.0)
    with pytest.raises(ValueError):
        nc.fit_critical_pressure_pa(EARTH_FLUX, albedo=1.0)
