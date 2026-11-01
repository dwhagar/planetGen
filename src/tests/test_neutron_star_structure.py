# tests/test_neutron_star_structure.py

"""
GEN.968: the Tolman-Oppenheimer-Volkoff integration of `physics.neutron_star_structure`, on equations of state built here,
so that it needs no downloaded file. The 21 real tables against the radii NICER measured are in
`test_neutron_star_structure_real.py`.
"""

import math

import numpy as np
import pytest

from planetgen.datasets.neutron_star_matter import EquationOfState
from planetgen.physics import neutron_star_structure as ns

_BARYON_MASS_MEV = 939.0


def _polytrope(k_km2, energy_low_km2=1e-12, energy_high_km2=2e-2, rows=400):
    """The table of P = K e^2 (a Newtonian n = 1 polytrope when P is small against e), rows spaced evenly in log e."""
    energy = np.geomspace(energy_low_km2, energy_high_km2, rows)
    pressure = k_km2 * energy ** 2
    energy_mev = energy * ns.MEV_FM3_PER_KM2
    return EquationOfState(0, energy_mev / _BARYON_MASS_MEV, energy_mev, pressure * ns.MEV_FM3_PER_KM2)


def test_a_solar_mass_is_one_point_four_seven_seven_kilometres_of_gravitational_length():
    assert ns.SOLAR_MASS_KM == pytest.approx(1.4766, abs=2e-4)


def test_one_mev_per_cubic_femtometre_is_thirteen_millionths_of_an_inverse_square_kilometre():
    assert 1.0 / ns.MEV_FM3_PER_KM2 == pytest.approx(1.3238e-6, rel=1e-4)


def test_a_light_polytropic_star_has_the_radius_of_the_newtonian_n_equals_one_polytrope():
    """The Lane-Emden solution of n = 1 has R = sqrt(pi K / 2 G), whatever the mass; at a central energy density of 1e-5
    km^-2 (a compactness of 0.2%) general relativity changes it by a few parts in a thousand."""
    k = 100.0
    eos = ns._Eos(_polytrope(k))
    star = ns.structure(eos, k * 1e-5 ** 2 * ns.MEV_FM3_PER_KM2)
    assert star.radius_km == pytest.approx(math.sqrt(math.pi * k / 2.0), rel=0.005)
    assert star.compactness < 0.005


def test_a_heavier_polytropic_star_is_more_compact_and_stays_outside_the_buchdahl_limit():
    eos = ns._Eos(_polytrope(100.0))
    light = ns.structure(eos, 100.0 * 1e-4 ** 2 * ns.MEV_FM3_PER_KM2)
    heavy = ns.structure(eos, 100.0 * 1e-3 ** 2 * ns.MEV_FM3_PER_KM2)
    assert heavy.compactness > light.compactness
    assert heavy.mass_solar > light.mass_solar
    assert heavy.compactness < 4.0 / 9.0


def test_the_polytrope_has_a_maximum_mass_and_a_curve_that_rises_to_it():
    eos = ns._Eos(_polytrope(100.0))
    masses = [ns.structure(eos, 100.0 * e ** 2 * ns.MEV_FM3_PER_KM2).mass_solar for e in (1e-4, 3e-4, 1e-3, 3e-3, 1e-2)]
    assert masses[0] < masses[1] < masses[2]
    assert max(masses) < 3.0 * math.sqrt(100.0) / ns.SOLAR_MASS_KM


def test_a_central_pressure_outside_the_table_is_refused():
    eos = ns._Eos(_polytrope(100.0))
    with pytest.raises(ValueError):
        ns.structure(eos, 1e9)


def test_an_unknown_equation_of_state_number_is_refused():
    with pytest.raises(ValueError):
        ns.structure(22, 100.0)


def test_a_table_whose_pressure_falls_is_refused():
    table = _polytrope(100.0)
    reversed_table = EquationOfState(0, table.baryon_density_fm3, table.energy_density_mev_fm3, table.pressure_mev_fm3[::-1])
    with pytest.raises(ValueError):
        ns._Eos(reversed_table)


# --- the committed table generation reads (no download needed) ---------------------------------------------------------------------

INDICES = range(1, 22)


def test_the_table_holds_all_21_equations_of_state_with_curves_that_rise():
    for i in INDICES:
        model = ns._table()[i]
        assert (np.diff(model["mass_solar"]) > 0).all() and 0.1 < model["mass_solar"][0] < 0.7, i
        assert model["mass_solar"][-1] == pytest.approx(model["maximum_mass_solar"], abs=2e-3), i
        assert (model["mass_solar"] * ns.SOLAR_MASS_KM / model["radius_km"] < 4.0 / 9.0).all(), i
        assert (model["crust_thickness_km"] > 0).all() and (model["crust_thickness_km"] < model["radius_km"]).all(), i


def test_nearly_every_equation_of_state_holds_up_a_two_solar_mass_pulsar():
    """PSR J0740+6620 weighs 2.08 +- 0.07 solar masses [R1404]; the compilation is built to hold it."""
    assert len([i for i in INDICES if ns.maximum_mass_solar(i) >= 2.0]) >= 19
    assert all(ns.maximum_mass_solar(i) < 3.0 for i in INDICES)
    assert len(ns.supporting(2.08)) == 19 and ns.supporting(2.9) == ()


def test_the_radius_of_a_one_point_four_solar_mass_star_spans_the_nicer_band():
    """[R3075]: R(1.4) = 12.45 +- 0.65 km. The 21 tables run from the soft ones near 11.8 km to the stiff ones near 13.8."""
    low, median, high = ns.radius_band_km(1.4)
    assert 11.5 < low < 12.3 and 13.3 < high < 14.2
    assert 12.45 - 0.65 < median < 12.45 + 1.0


def test_the_radii_agree_with_the_two_nicer_pulsars():
    """J0030+0451, 1.34 solar masses: 12.71 +1.14 -1.19 km [R3766]. J0740+6620, 2.08: 13.7 +2.6 -1.5 km [R3075]."""
    low, median, high = ns.radius_band_km(1.34)
    assert 12.71 - 1.19 < median < 12.71 + 1.14
    low, median, high = ns.radius_band_km(2.08)
    assert 13.7 - 1.5 < median < 13.7 + 2.6
    assert 12.35 - 0.75 < median < 12.35 + 0.75    # the combined band of [R3075] at 2.08 solar masses


def test_a_one_point_four_solar_mass_star_has_a_crust_about_a_kilometre_thick():
    for i in INDICES:
        assert 0.4 < ns.crust_thickness_km(i, 1.4) < 2.0, i


def test_a_mass_outside_an_equations_curve_has_no_radius():
    assert ns.radius_km(6, 2.5) is None and ns.radius_km(1, 0.1) is None and ns.crust_thickness_km(6, 2.5) is None
    with pytest.raises(ValueError):
        ns.radius_km(22, 1.4)
