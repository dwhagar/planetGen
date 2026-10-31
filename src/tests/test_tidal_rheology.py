# tests/test_tidal_rheology.py

"""
GEN.954: the tidal Love number and quality factor from temperature, melt and tide frequency (`planetgen.physics.tidal_rheology`), checked on
the Moon, Earth, Io and Europa and against the constant-Q rheology of `tidal_heating` that GEN.955 will replace.
"""

import ast
import inspect
import math

import pytest

from planetgen.physics import tidal_heating, tidal_rheology

ETA0 = 3.6e9
"""float: A viscosity prefactor in the middle of the drawn range (9.3e7 to 1.9e11 Pa s, geometric middle 1.3e9): with EA the viscosity of the old 1e12 and 250 kJ/mol at 1,600 K."""
EA = 3.25e5
"""float: And an activation energy in the middle of the drawn range (260 to 390 kJ/mol)."""
SOLIDUS, LIQUIDUS = 1700.0, 2200.0
"""float: A solidus and a liquidus 500 K apart (`tidal_heating.MELT_RANGE_K`) in the middle of the box's range; the wiring passes the body's own."""

MOON = (7.342e22, 1.7374e6)
EARTH = (5.972e24, 6.371e6)
IO = (8.932e22, 1.8216e6)
EUROPA = (4.80e22, 1.5608e6)
MONTH_RATE = 2.0 * math.pi / (27.32 * 86400.0)
"""float: The tide the Earth raises in the Moon, rad/s (one orbit)."""
SEMIDIURNAL_RATE = 2.0 * math.pi / (12.42 * 3600.0)
IO_RATE = 2.0 * math.pi / (1.769 * 86400.0)
EUROPA_RATE = 2.0 * math.pi / (3.551 * 86400.0)


def _rock(body, frequency, temperature_k, **keywords):
    return tidal_rheology.rock_response(*body, frequency, temperature_k, ETA0, EA, SOLIDUS, LIQUIDUS, **keywords)


def test_a_cold_stiff_mantle_has_the_elastic_love_number_and_a_large_q():
    """Far below the solidus the mantle does not relax in a tide: k2 is the rigid sphere's `tidal_heating.rigid_k2` and Q is in the hundreds."""
    cold = _rock(MOON, MONTH_RATE, 1000.0)
    assert cold.k2 == pytest.approx(tidal_heating.rigid_k2(*MOON), rel=0.01)
    assert cold.q > 200.0


def test_the_moon_comes_out_with_the_lunar_laser_ranging_k2_and_q():
    """The Moon: k2 0.0242 and Q 37.5 at the monthly tide (Williams and Boggs 2015, recalled); a mantle at 1,450 K gives them to a factor of 2."""
    moon = _rock(MOON, MONTH_RATE, 1450.0)
    assert moon.k2 == pytest.approx(0.0242, rel=0.15)
    assert 20.0 < moon.q < 75.0


def test_the_solid_earth_has_the_semidiurnal_tide_of_a_q_of_a_few_hundred():
    """The solid Earth's body tide has Q near 280 at the semi-diurnal frequency; the homogeneous sphere's k2 of 0.5
    is above the real 0.30 (the core and the layers are not here), which is why the module keeps to the form of `rigid_k2`."""
    earth = _rock(EARTH, SEMIDIURNAL_RATE, 1500.0)
    assert 100.0 < earth.q < 800.0
    assert earth.k2 == pytest.approx(tidal_heating.rigid_k2(*EARTH), rel=0.1)


def test_a_very_cold_mantle_keeps_the_loss_of_cold_rock():
    """Andrade `Q` grows without limit as a mantle cools; the ceiling keeps a cold small planet from being frictionless."""
    assert _rock(MOON, MONTH_RATE, 200.0).q == tidal_rheology.COLD_Q_CEILING
    assert _rock(MOON, MONTH_RATE, 200.0, maxwell=True).q > 1.0e6


def test_q_falls_as_the_mantle_warms_toward_the_solidus():
    """The signature of a real rheology and the thing constant Q cannot show: the same body dissipates more as it warms."""
    qs = [_rock(MOON, MONTH_RATE, t).q for t in (1000.0, 1200.0, 1400.0, 1600.0, 1700.0)]
    assert qs == sorted(qs, reverse=True)
    assert qs[0] / qs[-1] > 10.0


def test_a_faster_tide_is_flexed_less_dissipatively_when_the_mantle_is_warm():
    """Andrade creep: at fixed temperature a higher frequency is further from the relaxation time, so Q rises as the tide quickens."""
    slow = _rock(MOON, MONTH_RATE / 100.0, 1500.0).q
    fast = _rock(MOON, MONTH_RATE * 100.0, 1500.0).q
    assert fast > slow


def test_partial_melt_gives_io_the_k2_over_q_its_orbital_acceleration_needs():
    """Io at 30 percent melt: Lainey et al. 2009 find k2/Q = 0.015 +- 0.003 from its orbital acceleration (recalled); the rheology gives a few
    hundredths with no number fitted to Io (the fit `MOLTEN_K2_OVER_Q` of 0.015 is gone)."""
    io = _rock(IO, IO_RATE, SOLIDUS + 0.3 * (LIQUIDUS - SOLIDUS))
    assert 0.003 < io.k2_over_q < 0.1
    assert io.k2_over_q > 30.0 * _rock(IO, IO_RATE, 1500.0).k2_over_q


def test_dissipation_peaks_in_the_partial_melt_and_falls_at_both_ends():
    """k2/Q (the tide's power) rises from the solid mantle into the mush and falls again in a liquid, where the body no longer dissipates."""
    ratios = {t: _rock(IO, IO_RATE, t).k2_over_q for t in (1500.0, 1700.0, 1850.0, 2000.0, 2150.0)}
    assert ratios[1850.0] > ratios[1700.0] > ratios[1500.0]
    assert ratios[2150.0] < 0.05 * ratios[1850.0]


def test_a_magma_ocean_tends_to_the_fluid_love_number_and_no_dissipation():
    liquid = _rock(IO, IO_RATE, 2200.0)
    assert liquid.k2 == pytest.approx(1.5, rel=0.01)
    assert liquid.k2_over_q <= 1.5 / tidal_rheology.COLD_Q_CEILING


def test_the_maxwell_form_reproduces_henning_peak_and_the_andrade_form_fixes_its_cold_end():
    """Maxwell (Henning et al. 2009) dissipates only near its relaxation time; at 1,000 K its Q is beyond any measured rock's."""
    maxwell = _rock(MOON, MONTH_RATE, 1000.0, maxwell=True)
    andrade = _rock(MOON, MONTH_RATE, 1000.0)
    assert maxwell.q > 1.0e4 > andrade.q
    peak = max(_rock(MOON, MONTH_RATE, t, maxwell=True).k2_over_q for t in range(1500, 1800, 20))
    assert peak > 20.0 * maxwell.k2_over_q


def test_a_warm_ice_sphere_is_far_more_dissipative_than_a_cold_one():
    """A homogeneous ice sphere of Europa's size near its melting point is a hundred times more dissipative than one at 140 K, which is what the ice shell
    (`icy_tide`, GEN.960) integrates over its temperatures."""
    cold = tidal_rheology.ice_response(*EUROPA, EUROPA_RATE, 140.0)
    warm = tidal_rheology.ice_response(*EUROPA, EUROPA_RATE, 260.0)
    assert warm.k2_over_q > 100.0 * cold.k2_over_q


def test_k2_over_q_of_the_response_matches_its_parts():
    response = _rock(MOON, MONTH_RATE, 1450.0)
    assert response.k2_over_q == pytest.approx(response.k2 / response.q)
    assert tidal_rheology.Response(0.3, math.inf).k2_over_q == 0.0


def test_melt_fraction_runs_from_the_solidus_to_the_liquidus():
    assert tidal_rheology.melt_fraction(SOLIDUS - 100.0, SOLIDUS, LIQUIDUS) == 0.0
    assert tidal_rheology.melt_fraction(0.5 * (SOLIDUS + LIQUIDUS), SOLIDUS, LIQUIDUS) == pytest.approx(0.5)
    assert tidal_rheology.melt_fraction(LIQUIDUS + 100.0, SOLIDUS, LIQUIDUS) == 1.0
    with pytest.raises(ValueError):
        tidal_rheology.melt_fraction(1800.0, 2000.0, 1900.0)


def test_viscosity_is_arrhenius_then_weakened_by_melt_then_a_liquid_at_full_melt():
    solid = tidal_rheology.viscosity_pa_s(1600.0, ETA0, EA)
    assert solid == pytest.approx(ETA0 * math.exp(EA / (tidal_rheology.GAS_CONSTANT * 1600.0)))
    assert tidal_rheology.viscosity_pa_s(1500.0, ETA0, EA) > solid > tidal_rheology.viscosity_pa_s(1700.0, ETA0, EA)
    assert tidal_rheology.viscosity_pa_s(1700.0, ETA0, EA, 0.2) == pytest.approx(
        tidal_rheology.viscosity_pa_s(1700.0, ETA0, EA) * math.exp(-tidal_rheology.MELT_VISCOSITY_FACTOR * 0.2))
    critical = tidal_rheology.CRITICAL_MELT_FRACTION
    below = tidal_rheology.viscosity_pa_s(1900.0, ETA0, EA, critical - 1e-9)
    above = tidal_rheology.viscosity_pa_s(1900.0, ETA0, EA, critical + 1e-9)
    assert above == pytest.approx(below, rel=1e-6)
    assert tidal_rheology.viscosity_pa_s(1900.0, ETA0, EA, 1.0) == pytest.approx(tidal_rheology.LIQUID_VISCOSITY_PA_S)
    melts = [tidal_rheology.viscosity_pa_s(1900.0, ETA0, EA, m / 20.0) for m in range(21)]
    assert melts == sorted(melts, reverse=True)


def test_rigidity_falls_with_melt_to_a_floor():
    assert tidal_rheology.rigidity_pa(0.0) == tidal_rheology.ROCK.rigidity_pa == tidal_heating.ROCK_RIGIDITY_PA
    assert tidal_rheology.rigidity_pa(0.2) == pytest.approx(0.25 * tidal_rheology.ROCK.rigidity_pa)
    assert tidal_rheology.rigidity_pa(0.9) == pytest.approx(tidal_rheology.MELT_RIGIDITY_FLOOR * tidal_rheology.ROCK.rigidity_pa)


def test_ice_viscosity_is_pinned_at_the_melting_point():
    eta0 = tidal_rheology.arrhenius_prefactor(tidal_rheology.ICE_MELTING_VISCOSITY_PA_S, tidal_rheology.ICE_MELTING_K,
                                               tidal_rheology.ICE_ACTIVATION_J_MOL)
    at_melting = tidal_rheology.viscosity_pa_s(tidal_rheology.ICE_MELTING_K, eta0, tidal_rheology.ICE_ACTIVATION_J_MOL)
    assert at_melting == pytest.approx(tidal_rheology.ICE_MELTING_VISCOSITY_PA_S)


def test_the_tidal_frequency_of_a_locked_and_of_a_spinning_body():
    assert tidal_rheology.tidal_frequency_rad_s(1.0e-5) == 1.0e-5
    assert tidal_rheology.tidal_frequency_rad_s(1.0e-5, 7.0e-5) == pytest.approx(1.2e-4)
    assert tidal_rheology.tidal_frequency_rad_s(1.0e-5, 1.0e-5) == 0.0


def test_a_zero_frequency_is_refused():
    with pytest.raises(ValueError):
        _rock(MOON, 0.0, 1500.0)


def test_the_module_imports_none_of_the_users_it_will_have():
    """`tidal_heating`, `interior_heat` and `orbit_elements` will import this module (GEN.955), so it must not import any of them."""
    tree = ast.parse(inspect.getsource(tidal_rheology))
    imported = {alias.name for node in ast.walk(tree) if isinstance(node, ast.ImportFrom) for alias in node.names}
    assert imported.isdisjoint({"tidal_heating", "interior_heat", "orbit_elements", "tidal_spin"})


def test_a_mantle_far_colder_than_any_body_does_not_overflow_and_is_elastic_with_the_cold_ceiling():
    """A mantle at 50 K with an activation energy of 300 kJ/mol has an Arrhenius exponent over 700: its response is the elastic one."""
    frozen = tidal_rheology.rock_response(*MOON, MONTH_RATE, 50.0, ETA0, 3.0e5, SOLIDUS, LIQUIDUS)
    assert frozen.k2 == pytest.approx(tidal_heating.rigid_k2(*MOON), rel=0.01)
    assert frozen.q == tidal_rheology.COLD_Q_CEILING
