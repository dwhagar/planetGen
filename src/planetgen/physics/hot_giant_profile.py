# planetgen/physics/hot_giant_profile.py

"""
The temperature-pressure profile of a hot giant's atmosphere (GEN.1577, built by Lane 07 for Lane 14)
======================================================================================================

An analytic, semi-grey, irradiated atmosphere in radiative equilibrium: Guillot 2010 [R1715] equation 29, the profile the cloud deck of `hot_cloud_deck` is placed on.
The design note is docs/design/composition/hot-giant-profile.md; this is the module of it.

    T^4 = (3 T_int^4 / 4) (2/3 + tau)
        + (3 f T_irr^4 / 4) [ 2/3 + 1/(gamma sqrt 3) + (gamma/sqrt 3 - 1/(gamma sqrt 3)) exp(-gamma tau sqrt 3) ]

with `tau = kappa_th P / g` (a constant thermal opacity `kappa_th`, so the column mass is `P / g`), `gamma = kappa_v / kappa_th` the visible-to-thermal opacity ratio,
`T_int` the temperature of the planet's own heat flux (`(F_int / sigma)^(1/4)`) and `f T_irr^4` the absorbed starlight in temperature form. `T_irr` is taken from the
equilibrium temperature `T_eq` of a planet that shares the light round the whole globe (`climate_draws.equilibrium_temperature_k`, which the caller passes with no internal flux):
`T_irr^4 = 4 T_eq^4`, so `f = 1/4` gives back `T_eq` and `f = 1/2` is the dayside.

What it gives: a profile that rises with depth to a deep isotherm at `T_irr`-and-`gamma` set `((3/4) f T_irr^4 (2/3 + 1/(gamma sqrt 3)))^(1/4)` (for `T_int` = 0), a cool skin at the
top (`skin_k`), and, when `gamma` is large (strong visible absorbers, TiO and VO), a temperature that falls with depth: the inversion. Pressure runs over the whole range the deck
searches (1e-6 to 1e3 bar). Above about 10 bar real collision-induced and electron opacity make `kappa` grow with pressure, so the constant opacity underestimates the deep
temperature [R1715]; the deck's lowest clouds (chromium, manganese sulfide) can sit there, which is why the deep end is a stand-in (PF3) and the caller may floor the result with the
interior's adiabat. Pure functions: no draw, nothing stored, no call into the generator.

Use (the caller composes the inputs from what is already built):

    t_eq = climate_draws.equilibrium_temperature_k(albedo, starlight_w_m2)            # no internal flux: the planet's own heat enters as T_int
    profile = profile_function(t_eq, giant_draw.own_flux_w_m2(planet), gravity_m_s2)  # K at a pressure in bar
    hot_cloud_deck.deck(host, profile)
"""

import math

from planetgen.physics import constants

KAPPA_THERMAL_M2_KG = 1.0e-3
"""float: The constant thermal (infrared) opacity, m2/kg (Guillot's 1e-2 cm2/g for HD 209458 b). Stand-in PF1, tier T2, tweak L0."""

GAMMA_AT_2000_K = 0.6
"""float: The visible-to-thermal opacity ratio at an irradiation temperature of 2,000 K; `gamma` follows `GAMMA_AT_2000_K * sqrt(T_irr / 2000 K)` (Guillot's
`kappa_v = 6e-3 sqrt(T_irr / 2000 K)` cm2/g over `kappa_th` = 1e-2 cm2/g). Stand-in PF2, tier T2, tweak L0."""

HEAT_REDISTRIBUTION = 0.25
"""float: `f`, the share of the starlight spread over the planet: 1 for the substellar point, 1/2 for the dayside, 1/4 for the whole globe (Guillot's `f`)."""

SQRT3 = math.sqrt(3.0)


def irradiation_temperature_k(equilibrium_k):
    """`T_irr` for a planet whose whole-globe equilibrium temperature is `equilibrium_k`: `T_irr^4 = 4 T_eq^4`, so that `f = 1/4` returns `T_eq`."""
    return 4.0 ** 0.25 * equilibrium_k


def internal_temperature_k(internal_flux_w_m2):
    """`T_int` from the planet's own heat flux, `(F_int / sigma)^(1/4)`, K; 0 for no flux."""
    if internal_flux_w_m2 < 0.0:
        raise ValueError("the internal flux cannot be negative")
    return (internal_flux_w_m2 / constants.STEFAN_BOLTZMANN_CONSTANT) ** 0.25


def opacity_ratio(irradiation_k):
    """`gamma` for an irradiation temperature: `GAMMA_AT_2000_K * sqrt(T_irr / 2000 K)` (Guillot's rule for the visible opacity)."""
    return GAMMA_AT_2000_K * math.sqrt(irradiation_k / 2000.0)


def temperature_k(pressure_bar, equilibrium_k, internal_flux_w_m2, gravity_m_s2, gamma=None, redistribution=HEAT_REDISTRIBUTION, kappa_thermal=KAPPA_THERMAL_M2_KG):
    """
    The temperature, K, at `pressure_bar` (a positive pressure; `tau = kappa_th P / g`).

    Args:
        equilibrium_k: The whole-globe equilibrium temperature `T_eq`, K, with no internal flux.
        internal_flux_w_m2: The planet's own heat flux, W/m2 (`giant_draw.own_flux_w_m2`).
        gravity_m_s2: The gravity at the level, m/s2 (taken constant).
        gamma: The visible-to-thermal opacity ratio; `opacity_ratio` of the irradiation temperature when None.
        redistribution: `f`; 1/4 whole globe, 1/2 dayside, 1 substellar.
        kappa_thermal: The thermal opacity, m2/kg.

    Raises:
        ValueError: A non-positive pressure, gravity, opacity or `gamma`, or a negative temperature or flux.
    """
    if not (pressure_bar > 0.0 and gravity_m_s2 > 0.0 and kappa_thermal > 0.0 and equilibrium_k >= 0.0 and redistribution > 0.0):
        raise ValueError("the pressure, gravity, opacity and redistribution must be positive and the equilibrium temperature not negative")
    t_irr = irradiation_temperature_k(equilibrium_k)
    if gamma is None:
        gamma = opacity_ratio(t_irr)
    if not gamma > 0.0:
        raise ValueError("gamma must be positive")
    tau = kappa_thermal * pressure_bar * 1.0e5 / gravity_m_s2
    t_int4 = internal_temperature_k(internal_flux_w_m2) ** 4
    light = 0.75 * redistribution * t_irr ** 4 * (2.0 / 3.0 + 1.0 / (gamma * SQRT3) + (gamma / SQRT3 - 1.0 / (gamma * SQRT3)) * math.exp(-gamma * tau * SQRT3))
    return (0.75 * t_int4 * (2.0 / 3.0 + tau) + light) ** 0.25


def skin_k(equilibrium_k, internal_flux_w_m2, gamma=None, redistribution=HEAT_REDISTRIBUTION):
    """The temperature at the top (`tau` = 0), K: `T^4 = T_int^4 / 2 + 3 f T_irr^4 / 4 (2/3 + gamma/sqrt 3)` of equation 29 at `tau` = 0, with the `1/(gamma sqrt 3)` terms cancelling."""
    t_irr = irradiation_temperature_k(equilibrium_k)
    if gamma is None:
        gamma = opacity_ratio(t_irr)
    return (0.5 * internal_temperature_k(internal_flux_w_m2) ** 4 + 0.75 * redistribution * t_irr ** 4 * (2.0 / 3.0 + gamma / SQRT3)) ** 0.25


def profile_function(equilibrium_k, internal_flux_w_m2, gravity_m_s2, gamma=None, redistribution=HEAT_REDISTRIBUTION, kappa_thermal=KAPPA_THERMAL_M2_KG):
    """A function of the pressure in bar giving K, the form `hot_cloud_deck.cloud_base` and `deck` take as a profile."""
    return lambda pressure_bar: temperature_k(pressure_bar, equilibrium_k, internal_flux_w_m2, gravity_m_s2, gamma, redistribution, kappa_thermal)
