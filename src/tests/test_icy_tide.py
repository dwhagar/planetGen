# tests/test_icy_tide.py

"""
GEN.960: the tide in an ice shell (`planetgen.physics.icy_tide`): the heat flux and shell in equilibrium with their own tide, and the effective `k2/Q` the despinning and
the circularisation read.
"""

import math
import types

import pytest

from planetgen.physics import constants as pc
from planetgen.physics import hydrosphere, icy_tide, interior_heat, tidal_rheology

JUPITER, SATURN = pc.JUPITER_MASS_TO_KG, 5.683e26
EUROPA = (4.8e22, 1560.8e3, JUPITER, 6.711e8, 0.0094, 7.0e-4, 100.0)
ENCELADUS = (1.08e20, 252.1e3, SATURN, 2.3802e8, 0.0047, 3.0e-4, 75.0)
CALLISTO = (1.08e23, 2410.0e3, JUPITER, 1.883e9, 0.0074, 3.0e-4, 130.0)


def test_the_loss_modulus_peaks_where_the_ice_relaxes_at_the_tides_period():
    """Im(mu) is small in cold ice (it does not flow) and small in ice near its melting point (it flows too freely), with a maximum between."""
    n = 2.0e-5
    values = {t: icy_tide.loss_modulus_pa(t, n) for t in (100.0, 150.0, 200.0, 230.0, 250.0, 265.0, 272.0)}
    peak = max(values, key=values.get)
    assert 150.0 < peak < 272.0
    assert values[100.0] < 0.1 * values[peak] and values[272.0] < values[peak]


def test_the_shell_is_in_equilibrium_with_its_own_tide():
    """F = F_own + P / (4 pi R^2), with the tide's power from the strain and the shell the flux makes."""
    flux, power = icy_tide.flux_from_tide(*EUROPA)
    mass, radius, host, a, e, own, top = EUROPA
    assert flux == pytest.approx(own + power / (4.0 * math.pi * radius ** 2), rel=1e-9)
    n = math.sqrt(pc.G * host / a ** 3)
    strain = 3.0 * e * (1.0 + icy_tide.love_k2(mass, radius, n)) * (host / mass) * (radius / a) ** 3
    area = 4.0 * math.pi * radius ** 2
    direct = area * 0.5 * n * strain ** 2 * (hydrosphere.ICE_CONDUCTIVITY_A_W_M / flux) * icy_tide.dissipation_integral_pa(top, hydrosphere.WATER_MELTING_K, n)
    assert power == pytest.approx(direct, rel=1e-6)


def test_europa_enceladus_and_callisto_come_out_the_size_the_shells_are():
    """Nothing is fitted: Europa 1e11 to 1e12 W and a shell of 10 to 30 km, Enceladus within ten of its 16 GW, Callisto's tide a tenth of Europa's or less."""
    europa = icy_tide.flux_from_tide(*EUROPA)
    thickness = hydrosphere.ice_shell_thickness_km(europa[0], EUROPA[-1])
    assert 1.0e11 < europa[1] < 1.0e12 and 10.0 < thickness < 40.0
    assert 1.6e9 < icy_tide.flux_from_tide(*ENCELADUS)[1] < 1.6e11
    assert icy_tide.flux_from_tide(*CALLISTO)[1] < 0.2 * europa[1]


def test_a_stronger_tide_makes_a_thinner_shell_and_more_heat_and_no_eccentricity_none():
    mass, radius, host, a, e, own, top = EUROPA
    weak = icy_tide.flux_from_tide(mass, radius, host, a, e / 4.0, own, top)
    strong = icy_tide.flux_from_tide(mass, radius, host, a, e, own, top)
    assert strong[1] > weak[1] > 0.0 and strong[0] > weak[0]
    assert icy_tide.flux_from_tide(mass, radius, host, a, 0.0, own, top) == (own, 0.0)
    assert icy_tide.flux_from_tide(mass, radius, host, a, 1.0e-203, own, top)[1] == 0.0


def test_a_surface_at_the_melting_point_has_no_shell_to_heat():
    mass, radius, host, a, e, own, _ = EUROPA
    assert icy_tide.dissipation_integral_pa(hydrosphere.WATER_MELTING_K, hydrosphere.WATER_MELTING_K, 1.0e-5) == 0.0
    assert icy_tide.flux_from_tide(mass, radius, host, a, e, own, 280.0) == (own, 0.0)
    assert icy_tide.k2_over_q(mass, radius, 1.0e-5, 1.0e-3, 280.0) == 0.0


def test_the_effective_k2_over_q_carries_the_same_power_through_the_peale_law():
    """The despinning and the circularisation read `k2/Q`; with it the Peale formula returns the shell's power."""
    from planetgen.physics import tidal_heating
    for case in (EUROPA, ENCELADUS, CALLISTO):
        mass, radius, host, a, e, own, top = case
        flux, power = icy_tide.flux_from_tide(*case)
        n = math.sqrt(pc.G * host / a ** 3)
        q_eff = icy_tide.k2_over_q(mass, radius, n, flux, top)
        assert tidal_heating.peale_power_w(q_eff, tidal_heating.Forcing(host, a, e), radius) == pytest.approx(power, rel=1e-3)


def test_a_thin_hot_shell_cannot_dissipate_more_than_a_fluid_body_at_a_q_of_one():
    mass, radius, host, a, e, own, top = EUROPA
    n = 2.0e-5
    assert icy_tide.k2_over_q(mass, radius, n, icy_tide.MIN_FLUX_W_M2, top) == icy_tide.MAX_K2_OVER_Q
    assert icy_tide.k2_over_q(mass, radius, n, 1.0, top) < icy_tide.MAX_K2_OVER_Q


def test_an_icy_body_reads_its_shell_through_the_interior_heat_helpers():
    """`tidal_response` of an icy body is the homogeneous ice sphere's k2 and the Q its shell gives; the shell's flux is the stored heat plus the tide."""
    body = types.SimpleNamespace(mass=4.8e22, radius=1560.8, body_type="t", water_mass_fraction=0.5, surface_temperature=100.0,
                                 internal_heat_flux_w_m2=7.0e-4, tidal_flux_w_m2=0.0174)
    assert interior_heat.icy_shell_top_k(body) == 100.0
    assert interior_heat.icy_shell_flux_w_m2(body) == pytest.approx(0.0181)
    n = 2.0e-5
    k2, q = interior_heat.tidal_response(body, n)
    assert k2 == pytest.approx(icy_tide.love_k2(4.8e22, 1560.8e3, n))
    assert k2 / q == pytest.approx(icy_tide.k2_over_q(4.8e22, 1560.8e3, n, 0.0181, 100.0))
    del body.surface_temperature
    assert interior_heat.icy_shell_top_k(body) == icy_tide.TOP_FALLBACK_K
