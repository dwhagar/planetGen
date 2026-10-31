# planetgen/physics/two_stream.py

"""
Two-stream radiative transfer for a plane-parallel atmosphere of homogeneous layers (GEN.443;
docs/design/radiation-model.md section 3, docs/design/atmosphere-composition-parts/02-gas-optics-uv-visible.md 2.8).

A `Stack` is the response of a set of layers to sunlight at every wavelength of a grid at once (every field is an array
over the grid, except the incidence cosine `mu0`, a number): the reflectance and transmittance for a beam from the
star at cosine `mu0`, and for diffuse light from above and from below. Layers are added top to bottom with `add`, and
`ground_flux` closes the stack with a Lambertian ground. All fluxes are fractions of the flux the star puts on a
horizontal surface at the top of the atmosphere (`mu0` times the flux through a surface facing the star).

* **A scattering layer** (`scattering_layer`) is solved in the delta-Eddington two-stream approximation (Joseph,
  Wiscombe and Weinman 1976; coefficients of Meador and Weaver 1980): optical depth, single-scattering albedo and
  asymmetry factor per wavelength. With no absorption, no asymmetry and no ground it gives the research's closed form
  `T = [(2/3 + mu0) + (2/3 - mu0) exp(-tau/mu0)] / (4/3 + tau)`, which the research checked against a Monte Carlo run
  within 1.3 percent for `mu0` of 0.5 and 1 (part 02 2.8). The equations are solved exactly for the layer (not
  iterated), with the exponentials scaled so a thick or strongly absorbing layer cannot overflow.
* **An absorbing layer** (`absorbing_layer`), such as the ozone, takes the beam through `exp(-tau/mu0)` and diffuse light
  through `2 E3(tau)` (the exact value of the isotropic field; the factor-1.66 shortcut is 3 percent off at `tau` of 1).
* **The ground** (`ground`) reflects a fraction `albedo` of what reaches it, in all directions.

The diffuse response of a scattering layer is the beam response integrated over incidence (`2 mu dmu`, Gauss-Legendre).

Limits: the approximation is good to a few percent for `mu0` above 0.3 and under-reflects at low sun (3.5 percent at
`mu0` of 0.2, research 2.8); a layer is one homogeneous slab, so a gas that is high in the air and one that is low are
separate layers; there is no polarisation and no spherical geometry (use `mu0` of at least 0.2).
"""

from dataclasses import dataclass

import numpy as np
from scipy.special import expn

MAX_ALBEDO_OF_SINGLE_SCATTERING = 1.0 - 1e-6
"""float: The conservative limit; a layer with less absorption than this is solved as if it had this much, so the equations stay
non-degenerate (the error is of the order of the absorption removed, under one part in a million per optical depth)."""

_NODES, _WEIGHTS = np.polynomial.legendre.leggauss(8)
_MU = 0.5 * (_NODES + 1.0)
_MU_WEIGHT = 0.5 * _WEIGHTS * 2.0 * _MU
"""Gauss-Legendre nodes on (0, 1) and the weights of `integral 2 mu f(mu) dmu`."""


@dataclass(frozen=True)
class Stack:
    """
    A set of layers' response at every grid wavelength.

    Attributes:
        beam_reflectance (ndarray): Fraction of the star's horizontal flux reflected back out (all of it diffuse).
        beam_transmission (ndarray): Fraction passing through unscattered, `exp(-tau/mu0)`.
        beam_diffuse (ndarray): Fraction passing through scattered (diffuse light leaving the bottom).
        reflectance_above (ndarray): Reflectance for diffuse light arriving from above.
        reflectance_below (ndarray): Reflectance for diffuse light arriving from below.
        transmission_diffuse (ndarray): Transmittance for diffuse light (the same in both directions).
        mu0 (float): The cosine of the angle of the star from the zenith.
    """

    beam_reflectance: np.ndarray
    beam_transmission: np.ndarray
    beam_diffuse: np.ndarray
    reflectance_above: np.ndarray
    reflectance_below: np.ndarray
    transmission_diffuse: np.ndarray
    mu0: float

    def absorbed_beam(self):
        """The fraction of the beam the stack absorbs: `1 - reflected - transmitted`."""
        return 1.0 - self.beam_reflectance - self.beam_transmission - self.beam_diffuse


def _solve_beam(tau, omega, g, mu0):
    """
    Reflectance and the two transmissions of a homogeneous delta-Eddington layer for a beam of cosine `mu0`.

    Returns:
        tuple: (reflectance, direct transmission, diffuse transmission), each an array like `tau`.
    """
    tau = np.asarray(tau, dtype=float)
    omega = np.minimum(np.asarray(omega, dtype=float), MAX_ALBEDO_OF_SINGLE_SCATTERING)
    g = np.asarray(g, dtype=float)
    f = g * g                                           # delta scaling: the forward peak counts as unscattered
    tau_s = (1.0 - omega * f) * tau
    omega_s = omega * (1.0 - f) / (1.0 - omega * f)
    g_s = g / (1.0 + g)
    g1 = (7.0 - omega_s * (4.0 + 3.0 * g_s)) / 4.0
    g2 = -(1.0 - omega_s * (4.0 - 3.0 * g_s)) / 4.0
    g3 = (2.0 - 3.0 * g_s * mu0) / 4.0
    g4 = 1.0 - g3
    k = np.sqrt(np.maximum(g1 * g1 - g2 * g2, 1e-30))
    direct = np.exp(-tau_s / mu0)
    # particular solution c exp(-tau/mu0) of dy/dtau = A y + b exp(-tau/mu0), A = [[g1, -g2], [g2, -g1]]
    source_up, source_down = -omega_s / mu0 * g3, omega_s / mu0 * g4
    mu_safe = np.where(np.abs(1.0 - k * mu0) < 1e-7, mu0 * (1.0 + 1e-6), mu0)     # the resonance k mu0 = 1 has no particular solution of this form
    shift = 1.0 / mu_safe
    # (-shift I - A) c = b
    a11, a12, a21, a22 = -shift - g1, g2, -g2, -shift + g1
    det = a11 * a22 - a12 * a21
    c_up = (source_up * a22 - a12 * source_down) / det
    c_down = (a11 * source_down - a21 * source_up) / det
    # homogeneous solutions, scaled: v+ exp(-k (tau_s - tau)) and v- exp(-k tau); v+ = [g2, g1 - k], v- = [g2, g1 + k]
    scale = np.exp(-k * tau_s)
    vp_up, vp_down = g2, g1 - k
    vm_up, vm_down = g2, g1 + k
    # no downward diffuse at the top: c_down + a vp_down scale + b vm_down = 0
    # no upward diffuse at the bottom: c_up direct + a vp_up + b vm_up scale = 0
    m11, m12 = vp_down * scale, vm_down
    m21, m22 = vp_up, vm_up * scale
    r1, r2 = -c_down, -c_up * direct
    d = m11 * m22 - m12 * m21
    a = (r1 * m22 - m12 * r2) / d
    b = (m11 * r2 - m21 * r1) / d
    reflectance = c_up + a * vp_up * scale + b * vm_up
    diffuse = c_down * direct + a * vp_down + b * vm_down * scale
    return np.clip(reflectance, 0.0, 1.0), direct, np.clip(diffuse, 0.0, 1.0)


def scattering_layer(tau, omega, g, mu0):
    """
    A homogeneous scattering and absorbing layer.

    Args:
        tau (ndarray): Optical depth (scattering plus absorption) at each wavelength.
        omega (ndarray): Single-scattering albedo (the share of the extinction that scatters), 0 to 1.
        g (ndarray): Asymmetry factor of the scattering (Rayleigh 0, cloud droplets 0.85).
        mu0 (float): Cosine of the star's zenith angle, 0.2 to 1.
    """
    tau = np.asarray(tau, dtype=float)
    reflectance, direct, diffuse = _solve_beam(tau, omega, g, mu0)
    reflectance_diffuse = np.zeros_like(tau)
    transmission_diffuse = np.zeros_like(tau)
    for mu, weight in zip(_MU, _MU_WEIGHT):
        r, t_dir, t_dif = _solve_beam(tau, omega, g, mu)
        reflectance_diffuse += weight * r
        transmission_diffuse += weight * (t_dir + t_dif)
    return Stack(reflectance, direct, diffuse, reflectance_diffuse, reflectance_diffuse, transmission_diffuse, mu0)


def absorbing_layer(tau, mu0):
    """A layer that only absorbs (no scattering): the ozone."""
    tau = np.asarray(tau, dtype=float)
    zero = np.zeros_like(tau)
    return Stack(zero, np.exp(-tau / mu0), zero, zero, zero, 2.0 * expn(3, np.maximum(tau, 0.0)), mu0)


def ground(albedo, mu0):
    """A Lambertian ground of reflectance `albedo` (array or number) under every layer."""
    albedo = np.asarray(albedo, dtype=float)
    zero = np.zeros_like(albedo)
    return Stack(albedo, zero, zero, albedo, albedo, zero, mu0)


def add(top, bottom):
    """
    The stack that has `top` over `bottom`, with the light bouncing between them.

    Both must have been made for the same `mu0`. The result's `reflectance_below` is that of the combination seen from the
    bottom, so layers can be added one after the other from the top.
    """
    if top.mu0 != bottom.mu0:
        raise ValueError("stacks made for different sun angles cannot be added")
    denominator = 1.0 - top.reflectance_below * bottom.reflectance_above
    down = (top.beam_diffuse + top.reflectance_below * bottom.beam_reflectance * top.beam_transmission) / denominator
    up = bottom.beam_reflectance * top.beam_transmission + bottom.reflectance_above * down
    return Stack(
        beam_reflectance=top.beam_reflectance + top.transmission_diffuse * up,
        beam_transmission=top.beam_transmission * bottom.beam_transmission,
        beam_diffuse=top.beam_transmission * bottom.beam_diffuse + down * bottom.transmission_diffuse,
        reflectance_above=top.reflectance_above + top.transmission_diffuse ** 2 * bottom.reflectance_above / denominator,
        reflectance_below=bottom.reflectance_below + bottom.transmission_diffuse ** 2 * top.reflectance_below / denominator,
        transmission_diffuse=top.transmission_diffuse * bottom.transmission_diffuse / denominator,
        mu0=top.mu0,
    )


def stack_layers(layers):
    """The stack of `layers` (top first) added one after the other."""
    stack = layers[0]
    for layer in layers[1:]:
        stack = add(stack, layer)
    return stack


def ground_flux(atmosphere, albedo):
    """
    The flux that reaches the ground under `atmosphere` over a Lambertian ground of `albedo`, as a fraction of the star's
    horizontal flux at the top (direct plus diffuse, multiple bounces between ground and air included), and the
    fraction of the star's flux the ground then reflects back out of the atmosphere.

    Returns:
        tuple: (flux at the ground, flux leaving the top of the atmosphere).
    """
    albedo = np.asarray(albedo, dtype=float)
    denominator = 1.0 - atmosphere.reflectance_below * albedo
    down = (atmosphere.beam_diffuse + atmosphere.reflectance_below * albedo * atmosphere.beam_transmission) / denominator
    arriving = atmosphere.beam_transmission + down
    up = albedo * arriving
    outgoing = atmosphere.beam_reflectance + atmosphere.transmission_diffuse * up
    return arriving, outgoing
