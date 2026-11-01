# planetgen/physics/envelope_escape.py

"""
How much of a hydrogen-helium envelope a close-in planet loses to its star's XUV (GEN.822, docs/design/composition/envelope-escape-wiring.md)
=========================================================================================================================================

The share of its starting envelope that a planet of 1 to 30 Earth masses has lost by its age is the number rows e1 (the planet is a bare
core) and e2 (it keeps an envelope) of the atmosphere recipes weigh by. It comes from the hydro-based mass-loss rate of Kubyshkina et al.
2018 (`thermosphere.hydro_based_mass_loss_g_s`, GEN.738), integrated over the star's XUV history (`activity.coronal_xuv_at`) while the envelope
shrinks and the planet contracts with it (the radius is the law of [R2741], `bulk_structure.envelope_radius_earth`): as the
envelope thins, the radius falls, the escape parameter rises and the rate drops steeply, so a planet either loses its envelope early or
settles.

The step is the smaller of a fraction of the age and a few percent of the envelope, so a fast loss is followed and not jumped over; the loss
stops when the envelope falls below `BARE_BELOW` of the planet's mass (a bare core).

The loss starts at `ONSET_GY` (1 Gyr), not at the first ten million years. Taken from 10 Myr, the rate of Kubyshkina et al. (35 to 60 times the
energy-limited rate of 0.15 at a Lambda of 15) strips every core below 15 Earth masses inside 30 days whatever the starting envelope, and no
valley comes out (docs/design/composition/envelope-escape-wiring.md sections 1 and 5). The loss before `ONSET_GY`, the boil-off of the hot, inflated envelope and
the first, saturated hundred million years, is part of the starting envelope `f_env0`: it is the envelope that is left when the model begins.
"""

import math

from planetgen.physics import activity, bulk_structure, constants, thermosphere

ONSET_GY = 1.0
"""float: The age at which the integration of the hydro-based rate starts, Gyr; what is lost earlier is in `f_env0`. Tuned to the radius valley
(sweeps of 0.01 to 1 Gyr: from 10 Myr the rate strips every core under 15 Earth masses inside 30 days, and the valley comes out at 2.0
Earth radii or more with a slope of -0.16 or steeper at any weaker or later start but 1 Gyr), tweak level 3."""

BARE_BELOW = 0.002
"""float: An envelope under this fraction of the planet's mass is gone (a bare core)."""

MAX_ENVELOPE = 0.5
"""float: A starting envelope fraction is held at this (the law of Lopez and Fortney covers a few percent to tens of percent)."""

TEQ_OVER_INSOLATION_FOURTH = 255.0
"""float: `T_eq / S^(1/4)`, K, for a Bond albedo of 0.3 (the Earth)."""

AGE_STEP = 0.08
"""float: The largest step of the integration as a fraction of the age (numerical)."""

ENVELOPE_STEP = 0.03
"""float: The largest share of the envelope one step may remove (numerical)."""

MAX_STEPS = 4000
"""int: A bound on the steps of one integration (numerical; a planet that has not settled by then counts as stripped)."""

_GY_S = 1e9 * constants.SECONDS_PER_YEAR


def equilibrium_temperature_k(insolation):
    """The equilibrium temperature of a planet at `insolation` Earth units, K, for an albedo of 0.3."""
    return TEQ_OVER_INSOLATION_FOURTH * insolation ** 0.25


def envelope_left_earth(core_mass_earth, f_env0, insolation, host_mass_solar, luminosity_solar, age_gy, teq_k=None):
    """
    The envelope the planet still has at `age_gy`, Earth masses (0 when it has fallen under `BARE_BELOW` of the planet's mass).

    Args:
        core_mass_earth (float): The mass of the core, Earth masses (the rows pass the planet's mass, as before).
        f_env0 (float): The starting envelope as a fraction of the planet's mass.
        insolation (float): The bolometric flux at the planet, Earth units.
        host_mass_solar, luminosity_solar (float): The star's, for the XUV it puts out and the distance of the planet.
        age_gy (float): The age.
        teq_k (float, optional): The equilibrium temperature; from the insolation when not given. It is held to the 300 to 2,000 K of the
            fit of Kubyshkina et al. (a colder planet loses little in any case).
    """
    f_env0 = min(f_env0, MAX_ENVELOPE)
    start = f_env0 * core_mass_earth / (1.0 - f_env0)
    if age_gy <= ONSET_GY or f_env0 <= 0.0:
        return start
    teq = min(max(teq_k if teq_k else equilibrium_temperature_k(insolation), thermosphere.HYDRO_BASED_TEQ_RANGE_K[0]),
              thermosphere.HYDRO_BASED_TEQ_RANGE_K[1])
    luminosity_w = luminosity_solar * constants.SOLAR_LUMINOSITY
    distance_au = math.sqrt(luminosity_solar / insolation)
    area = 4.0 * math.pi * (distance_au * constants.AU_TO_M) ** 2
    envelope = start
    t = ONSET_GY * _GY_S
    end = age_gy * _GY_S
    earth_kg = constants.EARTH_MASS_TO_KG
    for _ in range(MAX_STEPS):
        age = t / _GY_S
        flux_erg_cm2_s = activity.coronal_xuv_at(host_mass_solar, age, luminosity_w) / area * 1e3
        mass = core_mass_earth + envelope
        radius = bulk_structure.envelope_radius_earth(mass, envelope / mass, insolation, age)
        rate = thermosphere.hydro_based_mass_loss_g_s(mass, radius, teq, distance_au, flux_erg_cm2_s) * 1e-3 / earth_kg
        step = min(AGE_STEP * t, end - t)
        if rate > 0.0:
            step = min(step, ENVELOPE_STEP * envelope / rate)
        envelope -= rate * step
        t += step
        if envelope <= BARE_BELOW * (core_mass_earth + envelope):
            return 0.0
        if t >= end:
            return envelope
    return 0.0


def lost_fraction(core_mass_earth, f_env0, insolation, host_mass_solar, luminosity_solar, age_gy, teq_k=None):
    """The share of its starting envelope that the planet has lost by `age_gy`, 0 to 1 (arguments as `envelope_left_earth`); 1 for a bare core."""
    f_env0 = min(f_env0, MAX_ENVELOPE)
    start = f_env0 * core_mass_earth / (1.0 - f_env0)
    if start <= 0.0:
        return 0.0
    return 1.0 - envelope_left_earth(core_mass_earth, f_env0, insolation, host_mass_solar, luminosity_solar, age_gy, teq_k) / start
