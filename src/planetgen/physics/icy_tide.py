# planetgen/physics/icy_tide.py

"""
The tide in an ice shell (GEN.960)
==================================

An icy body's tide heats its ice shell, not its whole volume (the ocean beneath moves with the tide and dissipates almost nothing). This replaces the one
`k2/Q` of 3e-4 that was tuned to Europa's 1.4e11 W (`tidal_heating.ICY_K2_OVER_Q`) with the shell's own rheology (`tidal_rheology.ICE`, Andrade, the viscosity
pinned at its melting point), after a standard treatment (recalled): the heat is made where the ice is warm enough that its relaxation time
is near the tide's period, a thin layer just above the shell's base, so it depends on the shell's thickness and the shell's thickness depends on the heat.

**The shell.** From the surface at `top_k` to the base at the melting point `base_k` the ice conducts a heat flux `F` (`k = A / T`, as `hydrosphere.ice_shell_thickness_km`),
so `T(z) = top exp(F z / A)` and `dz = A dT / (F T)`: the thickness is `D = (A / F) ln(base / top)`. The flux is the body's own (`F_own`, radiogenic from its rock) plus
the tide's: `F = F_own + P / (4 pi R^2)`.

**The dissipation.** The shell is strained by the tide through the whole of its thickness, by `epsilon` of the tide's amplitude over the radius, and the ice loses
`(1/2) omega epsilon^2 Im(mu(T, omega))` per volume, with `mu` the complex rigidity of `tidal_rheology.compliance` at the local temperature. For an eccentricity tide at
the orbital rate `n` the strain is `epsilon = 3 e (1 + k2) (M_h / M) (R / a)^3` (the amplitude of the tide's change over an orbit, `3 e h2` times the equilibrium bulge, with
`h2 = 1 + k2` and `k2` of a homogeneous ice sphere, `tidal_rheology.ice_response`), and for the semi-diurnal tide of a spinning body (`dissipation`) the same bulge with the
factor `sqrt(9/7)` that makes the two agree for the same `k2/Q` (the 21/2 of [R3464] against the 3 of the despinning torque). The power is

    P = 4 pi R^2 (1/2) omega epsilon^2 (A / F) I(omega),        I(omega) = integral of Im(mu) dT / T from top to base

and `I` does not depend on `F` (the profile only sets where in the shell each temperature is), so `F = F_own + c / F` with `c = (1/2) omega epsilon^2 A I`:
`F = (F_own + sqrt(F_own^2 + 4 c)) / 2`, one root, the shell in equilibrium with its own tide. Nothing is tuned: Europa comes out at 5e11 W and a 31 km shell, Enceladus
at 1e11 W (6 km), Callisto at 3e10 W; the observed figures are Europa 1e11 to 1e12 W and a shell of 10 to 30 km, Enceladus 1.6e10 W. Stand-ins are IC1 to IC5
of `docs/design/research/gap-filling-plan.md` section 6.23.
"""

import math

from scipy import integrate

from planetgen.physics import constants, hydrosphere, tidal_rheology

STRAIN_FACTOR = 1.0
"""float: The strain in the shell as a multiple of the tide's amplitude over the body's radius (the shell follows the ocean under it, which follows the equilibrium tide). Source tier
T6, tweak level L1 (an order-of-one choice; Europa's strain of 1e-5 to 1e-4 is the range it gives)."""
ENERGY_RATIO = 9.0 / 7.0
"""float: The ratio of the strain squared of the semi-diurnal tide of a spinning body to that of the eccentricity tide for the same `k2/Q`, from the two standard laws
(the 21/2 of the power and the 3 of the torque, `spin_torques.GRAVITATIONAL_TORQUE_FACTOR`): `P_e = 10.5 K n e^2` and `P_s = N sigma / 2 = 1.5 K sigma`. Source tier T5, tweak level L0."""
LOVE_TEMPERATURE_K = 250.0
"""float: The temperature at which the ice sphere's `k2` is taken for `h2` of the tide (a body with an ocean under its shell is flexed as a fluid one is; the homogeneous sphere's
value varies from 0.25 to 0.5 over the shell's temperatures). Source tier T6, tweak level L1."""
TOP_FALLBACK_K = 100.0
"""float: The surface temperature of an icy body whose climate is not yet known (an icy moon's orbit is circularised before its air and climate are made). Source tier T6, tweak
level L1 (the cold end of the Galilean moons, 90 to 130 K)."""
MAX_K2_OVER_Q = 1.5
"""float: The most `k2/Q` a body can have: the fluid Love number over a quality factor of 1. A very thin hot shell would otherwise return more. Source tier T5, tweak level L0."""
OWN_FLUX_FALLBACK_W_M2 = 1.0e-3
"""float: The heat flux of an icy body whose rock's heat is not known (the radiogenic flux of the icy moons is 1e-3 to 1e-2 W/m2). Source tier T6, tweak level L1."""
MIN_FLUX_W_M2 = 1.0e-6
"""float: The least heat flux the shell is given, W/m2 (`hydrosphere` uses the same floor)."""


def loss_modulus_pa(temperature_k, frequency_rad_s):
    """`Im(mu)` of ice at `temperature_k` flexed at `frequency_rad_s`, Pa: the part of the complex rigidity out of phase with the strain."""
    eta0 = tidal_rheology.arrhenius_prefactor(tidal_rheology.ICE_MELTING_VISCOSITY_PA_S, tidal_rheology.ICE_MELTING_K, tidal_rheology.ICE_ACTIVATION_J_MOL)
    viscosity = tidal_rheology.viscosity_pa_s(temperature_k, eta0, tidal_rheology.ICE_ACTIVATION_J_MOL)
    compliance = tidal_rheology.compliance(frequency_rad_s, tidal_rheology.ICE.rigidity_pa, viscosity, tidal_rheology.ICE)
    return (1.0 / compliance).imag


def dissipation_integral_pa(top_k, base_k, frequency_rad_s):
    """`I(omega)`: the integral of `Im(mu) dT / T` from `top_k` to `base_k`, Pa; 0 when the surface is at or above the melting point (no shell)."""
    if top_k >= base_k:
        return 0.0
    value, _ = integrate.quad(lambda ln_t: loss_modulus_pa(math.exp(ln_t), frequency_rad_s), math.log(top_k), math.log(base_k), epsabs=0.0, epsrel=1.0e-6, limit=200)
    return max(value, 0.0)


def love_k2(mass_kg, radius_m, frequency_rad_s):
    """`k2` of the homogeneous ice sphere at `LOVE_TEMPERATURE_K`, which gives the tide's `h2 = 1 + k2`."""
    return tidal_rheology.ice_response(mass_kg, radius_m, frequency_rad_s, LOVE_TEMPERATURE_K).k2


def flux_from_tide(mass_kg, radius_m, host_mass_kg, semi_major_axis_m, eccentricity, own_flux_w_m2, top_k, base_k=hydrosphere.WATER_MELTING_K):
    """
    `(F, P)`: the heat flux the shell conducts, W/m2, and the power its tide makes, W, for a body on an orbit of `eccentricity` about a host of `host_mass_kg`, with `own_flux_w_m2` of
    its own heat coming up into the shell. The tide's flux is `F - own_flux_w_m2`.
    """
    own = max(own_flux_w_m2, MIN_FLUX_W_M2)
    if eccentricity <= 0.0 or top_k >= base_k:
        return own, 0.0
    rate = math.sqrt(constants.G * host_mass_kg / semi_major_axis_m ** 3)
    h2 = 1.0 + love_k2(mass_kg, radius_m, rate)
    strain = STRAIN_FACTOR * 3.0 * eccentricity * h2 * (host_mass_kg / mass_kg) * (radius_m / semi_major_axis_m) ** 3
    c = 0.5 * rate * strain ** 2 * hydrosphere.ICE_CONDUCTIVITY_A_W_M * dissipation_integral_pa(top_k, base_k, rate)
    tide = 2.0 * c / (own + math.sqrt(own ** 2 + 4.0 * c))      # F - F_own without the cancellation
    return own + tide, 4.0 * math.pi * radius_m ** 2 * tide


def k2_over_q(mass_kg, radius_m, frequency_rad_s, flux_w_m2, top_k, base_k=hydrosphere.WATER_MELTING_K):
    """
    The effective `k2/Q` of the shell for a tide of `frequency_rad_s` (the semi-diurnal frequency of a spinning body, or the orbital rate of an eccentricity tide) with the shell
    conducting `flux_w_m2`: the dissipation the strain above makes, over what a body of that `k2/Q` would make in the torque law (`spin_torques.Torques.gravitational`),
    `(4 pi / 3) (9/7) (1 + k2)^2 R^3 A I / (F G M^2)`, never more than `MAX_K2_OVER_Q` (IC5). 0 when there is no shell.
    """
    integral = dissipation_integral_pa(top_k, base_k, frequency_rad_s)
    if integral <= 0.0:
        return 0.0
    h2 = 1.0 + love_k2(mass_kg, radius_m, frequency_rad_s)
    flux = max(flux_w_m2, MIN_FLUX_W_M2)
    shell = (4.0 * math.pi / 3.0 * ENERGY_RATIO * STRAIN_FACTOR ** 2 * h2 ** 2 * radius_m ** 3 * hydrosphere.ICE_CONDUCTIVITY_A_W_M * integral
             / (flux * constants.G * mass_kg ** 2))
    return min(shell, MAX_K2_OVER_Q)
