"""
The two-stream layers and their addition (GEN.443, `planetgen.physics.two_stream`): the closed form the research checked
against a Monte Carlo run, energy conservation, the exact equations solved another way, and the ground bounce.
"""

import numpy as np
import pytest
from scipy.integrate import solve_bvp
from scipy.special import expn

from planetgen.physics import two_stream as ts

TAU = np.array([0.1, 0.5, 1.0, 3.0, 10.0])


def _eddington(tau, mu0):
    """Part 02 2.8: the transmission of a conservative isotropic layer (direct plus diffuse)."""
    return ((2.0 / 3.0 + mu0) + (2.0 / 3.0 - mu0) * np.exp(-tau / mu0)) / (4.0 / 3.0 + tau)


@pytest.mark.parametrize("mu0", [1.0, 0.5, 0.2])
def test_a_conservative_isotropic_layer_is_the_research_closed_form(mu0):
    layer = ts.scattering_layer(TAU, np.ones(5), np.zeros(5), mu0)
    total = layer.beam_transmission + layer.beam_diffuse
    assert total == pytest.approx(_eddington(TAU, mu0), rel=1e-3)
    assert layer.beam_reflectance + total == pytest.approx(np.ones(5), abs=1e-3)


def test_a_layer_of_no_depth_changes_nothing():
    layer = ts.scattering_layer(np.array([1e-9]), np.array([0.9]), np.array([0.5]), 0.6)
    assert layer.beam_reflectance[0] == pytest.approx(0.0, abs=1e-6)
    assert layer.beam_transmission[0] + layer.beam_diffuse[0] == pytest.approx(1.0, abs=1e-6)


@pytest.mark.parametrize("omega, g", [(0.9, 0.0), (0.999, 0.85), (0.5, 0.7), (0.1, 0.0), (0.0, 0.0)])
def test_energy_is_conserved_and_absorption_is_not_negative(omega, g):
    layer = ts.scattering_layer(TAU, np.full(5, omega), np.full(5, g), 0.7)
    assert np.all(layer.absorbed_beam() > -1e-9)
    assert np.all(layer.beam_reflectance >= 0.0)
    assert np.all(layer.reflectance_above + layer.transmission_diffuse <= 1.0 + 1e-6)


def test_an_absorber_in_a_scattering_layer_lowers_the_transmission_and_the_reflection():
    clear = ts.scattering_layer(TAU, np.full(5, 0.999999), np.zeros(5), 1.0)
    dirty = ts.scattering_layer(TAU, np.full(5, 0.7), np.zeros(5), 1.0)
    assert np.all(dirty.beam_transmission + dirty.beam_diffuse < clear.beam_transmission + clear.beam_diffuse)
    assert np.all(dirty.beam_reflectance < clear.beam_reflectance)


def test_forward_scattering_clouds_reflect_less_than_isotropic_ones_and_follow_the_monte_carlo_table():
    isotropic = ts.scattering_layer(TAU, np.ones(5), np.zeros(5), 0.5)
    cloud = ts.scattering_layer(TAU, np.ones(5), np.full(5, 0.85), 0.5)
    assert np.all(cloud.beam_reflectance < isotropic.beam_reflectance)
    # part 06 3.10 (scalar Monte Carlo, g 0.85, omega 1, mu0 0.5): 0.462 at tau 5 and 0.870 at tau 50
    thick = ts.scattering_layer(np.array([5.0, 50.0]), np.ones(2), np.full(2, 0.85), 0.5)
    assert thick.beam_reflectance == pytest.approx([0.462, 0.870], rel=0.10)


def test_thick_layers_do_not_overflow():
    layer = ts.scattering_layer(np.array([1e3, 1e5]), np.array([0.3, 0.999]), np.array([0.0, 0.85]), 0.5)
    assert np.all(np.isfinite(layer.beam_reflectance)) and np.all(np.isfinite(layer.transmission_diffuse))
    assert layer.beam_reflectance[0] < 0.1 and layer.beam_diffuse[0] < 1e-12


def test_the_layer_matches_the_two_stream_equations_solved_numerically():
    tau_total, omega, g, mu0 = 2.0, 0.9, 0.5, 0.6
    f = g * g
    tau_s = (1 - omega * f) * tau_total
    omega_s = omega * (1 - f) / (1 - omega * f)
    g_s = g / (1 + g)
    g1, g2 = (7 - omega_s * (4 + 3 * g_s)) / 4, -(1 - omega_s * (4 - 3 * g_s)) / 4
    g3 = (2 - 3 * g_s * mu0) / 4
    g4 = 1 - g3

    def rhs(t, y):
        beam = np.exp(-t / mu0) / mu0
        return np.vstack([g1 * y[0] - g2 * y[1] - omega_s * g3 * beam, g2 * y[0] - g1 * y[1] + omega_s * g4 * beam])

    solution = solve_bvp(rhs, lambda ya, yb: np.array([yb[0], ya[1]]), np.linspace(0, tau_s, 50), np.zeros((2, 50)), tol=1e-10)
    assert solution.success
    layer = ts.scattering_layer(np.array([tau_total]), np.array([omega]), np.array([g]), mu0)
    assert layer.beam_reflectance[0] == pytest.approx(solution.sol(0)[0], rel=1e-6)
    assert layer.beam_diffuse[0] == pytest.approx(solution.sol(tau_s)[1], rel=1e-6)


def test_the_absorbing_layer_takes_the_beam_by_beer_lambert_and_diffuse_light_exactly():
    layer = ts.absorbing_layer(np.array([0.0, 1.0, 5.0]), 0.5)
    assert layer.beam_transmission == pytest.approx(np.exp(-np.array([0.0, 1.0, 5.0]) / 0.5))
    assert layer.transmission_diffuse == pytest.approx(2.0 * expn(3, np.array([0.0, 1.0, 5.0])))
    assert layer.transmission_diffuse[0] == pytest.approx(1.0)
    assert layer.transmission_diffuse[1] == pytest.approx(0.2194, abs=1e-3)


def _vacuum(wavelengths, mu0):
    """A stack of no layers: everything passes."""
    one, zero = np.ones(len(wavelengths)), np.zeros(len(wavelengths))
    return ts.Stack(zero, one, zero, zero, zero, one, mu0)


def test_nothing_changes_when_a_vacuum_is_added():
    layer = ts.scattering_layer(TAU, np.full(5, 0.9), np.full(5, 0.3), 0.8)
    both = ts.add(_vacuum(TAU, 0.8), layer)
    assert both.beam_reflectance == pytest.approx(layer.beam_reflectance)
    assert both.beam_diffuse == pytest.approx(layer.beam_diffuse)
    assert both.transmission_diffuse == pytest.approx(layer.transmission_diffuse)
    other = ts.add(layer, _vacuum(TAU, 0.8))
    assert other.beam_reflectance == pytest.approx(layer.beam_reflectance)


def test_two_halves_of_a_layer_add_up_to_the_whole_within_the_diffuse_quadrature():
    whole = ts.scattering_layer(np.array([1.0, 4.0]), np.array([0.95, 0.99]), np.array([0.0, 0.0]), 0.7)
    half = ts.scattering_layer(np.array([0.5, 2.0]), np.array([0.95, 0.99]), np.array([0.0, 0.0]), 0.7)
    joined = ts.add(half, half)
    assert joined.beam_reflectance == pytest.approx(whole.beam_reflectance, rel=0.03)
    assert joined.beam_transmission + joined.beam_diffuse == pytest.approx(whole.beam_transmission + whole.beam_diffuse, rel=0.03)
    assert joined.absorbed_beam() == pytest.approx(whole.absorbed_beam(), abs=0.01)


def test_adding_in_either_order_with_mirrored_layers_conserves_energy():
    a = ts.scattering_layer(np.array([0.3]), np.array([0.8]), np.array([0.0]), 0.5)
    b = ts.absorbing_layer(np.array([0.7]), 0.5)
    both = ts.add(b, a)
    assert both.absorbed_beam()[0] > 0.0
    assert both.beam_reflectance[0] + both.beam_transmission[0] + both.beam_diffuse[0] < 1.0


def test_a_bright_ground_under_a_conservative_layer_follows_the_research_formula():
    # part 02 2.8: R_tot = R + T A (1 - s) / (1 - A s), s the layer's spherical albedo
    mu0 = 0.5
    tau = np.array([0.1, 0.5, 1.0, 2.0])
    layer = ts.scattering_layer(tau, np.ones(4), np.zeros(4), mu0)
    for albedo in (0.0, 0.3, 0.8):
        arriving, outgoing = ts.ground_flux(layer, np.full(4, albedo))
        s = layer.reflectance_above
        transmission = layer.beam_transmission + layer.beam_diffuse
        expected = layer.beam_reflectance + transmission * albedo * (1.0 - s) / (1.0 - albedo * s)
        assert outgoing == pytest.approx(expected, rel=2e-3)
        assert arriving + outgoing - albedo * arriving == pytest.approx(np.ones(4), abs=2e-3)       # nothing is absorbed


def test_the_ground_flux_of_a_clear_sky_is_the_beam_and_the_albedo_does_not_matter():
    arriving, outgoing = ts.ground_flux(_vacuum(TAU, 0.5), np.full(5, 0.4))
    assert arriving == pytest.approx(np.ones(5))
    assert outgoing == pytest.approx(np.full(5, 0.4))


def test_stacks_made_for_different_sun_angles_are_refused():
    with pytest.raises(ValueError):
        ts.add(ts.absorbing_layer(np.array([1.0]), 0.5), ts.absorbing_layer(np.array([1.0]), 0.7))
