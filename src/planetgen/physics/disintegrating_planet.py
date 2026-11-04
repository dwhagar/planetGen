# planetgen/physics/disintegrating_planet.py

"""
A disintegrating rocky planet: its wind, its lifetime and when a body can be one (GEN.1041, W14 of docs/design/systems/weird-worlds.md, registry 6.73)
=======================================================================================================================================================

Kepler-1520 b (KIC 12557548 b), K2-22 b and KOI-2700 b dim by a changing fraction with a dust tail. [R3487]
show the dayside of a rocky planet under about 0.1 Earth mass at over 2000 K lifts a wind of its own rock vapour, and the gravity of such a small body is
too weak to hold it: the mass loss goes as the vapour density at the surface times the exponential of the escape energy over the thermal energy.

The wind is the isothermal Parker wind of the vapour (eq. 27 of the paper, a solid angle of one steradian):

    c^2 = gamma k T / mu                 mu = 30 m_H and gamma = 1.3 (eq. 8), the mean mass and index of the silicate vapour
    P_vapor = exp(-m L / (k T) + b)      eq. 7 with the olivine constants of 2.1
    rho_0 = mu P_vapor / (k T)           eq. 9
    r_c = G M / (2 c^2)                  the sonic radius
    rho_c = rho_0 exp(3/2 - G M / (R c^2))
    Mdot = Omega r_c^2 rho_c c

The paper's full model, with its radiative cooling of the wind, gives a rate about 40 times lower than the isothermal one (3.2), taken here as the factor
`NON_ISOTHERMAL_FACTOR`. The dayside temperature is `T_star sqrt(R_star / a)` (eq. 6, no optical depth). The planet loses mass for `DUTY_CYCLE` of
the time (3.3). The lifetime integrates `dM / (f Mdot(M))` from nothing to the mass: the lighter the body, the faster it goes.

The tail (the dust the wind carries) and the transit depth are not here: the paper gives no length or depth, only the dust loading as a ratio.

Use:

    teq = dayside_temperature_k(star_teff_k, star_radius_m, distance_m)
    rate = mass_loss_rate_kg_s(mass_earth, teq)            # kg/s
    life = lifetime_gy(mass_earth, teq)                    # Gyr
    can_disintegrate(mass_earth, teq, period_d)            # the gate of the oddity
"""

import math

from planetgen.physics import bulk_structure, constants, unit_registry as units

_kg_to_mearth = units.converter("kg", "Mearth")
_mearth_to_kg = units.converter("Mearth", "kg")

OLIVINE_MOLECULAR_MASS_H = 169.0
"""float: The mass of a molecule of olivine vapour in hydrogen masses, `m` of eq. 7 (row DP1)."""

OLIVINE_SUBLIMATION_ERG_G = 3.21e10
"""float: The latent heat of sublimation of olivine, erg/g, `L_sub` of eq. 7 (row DP1)."""

OLIVINE_VAPOUR_CONSTANT_DYNE_CM2 = 6.72e14
"""float: `exp(b)` of eq. 7 for olivine, dyne/cm^2 (row DP1)."""

VAPOUR_MEAN_MASS_H = 30.0
"""float: The mean mass of the vapour in hydrogen masses, `mu` of eq. 8 (row DP2)."""

ADIABATIC_INDEX = 1.3
"""float: The index of the vapour, `gamma` of eq. 8, in the sound speed (row DP2). With 1 the planet that survives 10 Gyr is half as heavy as the paper's."""

NON_ISOTHERMAL_FACTOR = 1.0 / 40.0
"""float: The full model's maximum rate over the isothermal one (3.2: "about 40 times lower"), row DP3."""

SOLID_ANGLE_SR = 1.0
"""float: The solid angle of the wind, `Omega` of eq. 27 (row DP3)."""

DUTY_CYCLE = 0.5
"""float: The share of the time the planet loses mass at its maximum rate (3.3, row DP4)."""

DIVIDING_MASS_EARTH = 0.11
"""float: The mass over which a planet at 2100 K or so survives 10 Gyr in the isothermal wind (3.3, row DP5); the paper's abstract says about 0.1."""

MIN_DAYSIDE_TEMPERATURE_K = 2000.0
"""float: The dayside temperature over which the surface vaporises and drives the wind (the abstract, row DP5)."""

MAX_PERIOD_D = 1.0
"""float: The orbital period under which a disintegrating planet orbits (the design note's gate, row DP5)."""

CORE_MASS_FRACTION = 0.32
"""float: The core fraction of the Earth-like body whose radius the wind is made on (`bulk_structure.small_body_radius_earth`, row DP2)."""

GONE_BELOW = 1.0e-2
"""float: The share of its mass under which a planet is gone (the lower end of the lifetime integral; the last hundredth goes in a thousandth of the time)."""

_GY_S = 1.0e9 * constants.SECONDS_PER_YEAR
_H = constants.HYDROGEN_ATOM_MASS_KG
_K = constants.BOLTZMANN
_G = constants.G


def dayside_temperature_k(star_teff_k, star_radius_m, distance_m):
    """The temperature of the dayside surface, K: `T_star sqrt(R_star / a)` (eq. 6 with no optical depth)."""
    return star_teff_k * math.sqrt(star_radius_m / distance_m)


def vapour_pressure_pa(temperature_k):
    """The olivine vapour pressure at `temperature_k`, Pa (eq. 7)."""
    energy_j = OLIVINE_MOLECULAR_MASS_H * _H * OLIVINE_SUBLIMATION_ERG_G * 1.0e-4      # erg/g to J/kg is 1e-4
    return OLIVINE_VAPOUR_CONSTANT_DYNE_CM2 * 0.1 * math.exp(-energy_j / (_K * temperature_k))


def radius_m(mass_earth):
    """The radius of the body, m: an Earth-like rock of its mass."""
    return bulk_structure.small_body_radius_earth(mass_earth, CORE_MASS_FRACTION) * constants.EARTH_RADIUS_KM * 1.0e3


def mass_loss_rate_kg_s(mass_earth, temperature_k, factor=NON_ISOTHERMAL_FACTOR):
    """The mass the wind carries away while it blows, kg/s: the isothermal Parker wind of the vapour times `factor` (the full model's, or 1 for the
    isothermal wind of the paper's first estimate)."""
    mass_kg = _mearth_to_kg(mass_earth)
    mu = VAPOUR_MEAN_MASS_H * _H
    sound_squared = ADIABATIC_INDEX * _K * temperature_k / mu
    surface_density = mu * vapour_pressure_pa(temperature_k) / (_K * temperature_k)
    sonic_m = _G * mass_kg / (2.0 * sound_squared)
    exponent = 1.5 - _G * mass_kg / (radius_m(mass_earth) * sound_squared)
    density = surface_density * math.exp(max(exponent, -700.0))
    return factor * SOLID_ANGLE_SR * sonic_m ** 2 * density * math.sqrt(sound_squared)


def mass_loss_earth_per_gy(mass_earth, temperature_k, factor=NON_ISOTHERMAL_FACTOR):
    """`mass_loss_rate_kg_s` in Earth masses per billion years."""
    return _kg_to_mearth(mass_loss_rate_kg_s(mass_earth, temperature_k, factor) * _GY_S)


def lifetime_gy(mass_earth, temperature_k, duty_cycle=DUTY_CYCLE, steps=400, factor=NON_ISOTHERMAL_FACTOR):
    """The time the planet lasts, Gyr: `dM / (f Mdot(M))` summed (trapezoid on log M) from a hundredth of its mass to its mass. The lighter body
    loses mass faster, so what is under a hundredth of the mass adds a negligible time."""
    low = math.log(mass_earth * 1.0e-2)
    step = (math.log(mass_earth) - low) / steps
    total, previous = 0.0, None
    for i in range(steps + 1):
        mass = math.exp(low + i * step)
        value = mass / (duty_cycle * mass_loss_earth_per_gy(mass, temperature_k, factor))      # d(ln M) weighted
        if previous is not None:
            total += 0.5 * (previous + value) * step
        previous = value
    return total


def remaining_mass_earth(mass_earth, temperature_k, age_gy, duty_cycle=DUTY_CYCLE):
    """The mass left after `age_gy` of loss, Earth masses: 0 when the planet is gone. The mass falls at its own rate, which rises as it falls."""
    mass, elapsed = mass_earth, 0.0
    while mass > mass_earth * GONE_BELOW and elapsed < age_gy:
        rate = duty_cycle * mass_loss_earth_per_gy(mass, temperature_k)
        dt = min(age_gy - elapsed, 0.02 * mass / rate)
        mass -= rate * dt
        elapsed += dt
    return mass if mass > mass_earth * GONE_BELOW else 0.0


def can_disintegrate(mass_earth, temperature_k, period_d):
    """Whether a body can be a disintegrating planet: rock under the dividing mass, a dayside over 2000 K and an orbit under a day."""
    return mass_earth <= DIVIDING_MASS_EARTH and temperature_k >= MIN_DAYSIDE_TEMPERATURE_K and period_d <= MAX_PERIOD_D
