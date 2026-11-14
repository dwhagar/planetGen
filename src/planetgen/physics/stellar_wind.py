# planetgen/physics/stellar_wind.py

"""
A star's wind (the tiered mass-loss rate and terminal speed), shared by the heliosphere (`Star.calculate_heliosphere`) and by the light of a black hole that a companion's wind feeds
(`dark_light`, GEN.1501): one model, so the star that sets the size of its bubble is the one that sets how much a black hole beside it takes in.
"""

import math

from planetgen.physics import constants
from planetgen.physics import unit_registry as units

_km_to_m = units.converter("km", "m")
_kg_to_msun = units.converter("kg", "Msun")
_msun_to_kg = units.converter("Msun", "kg")


def wind(mass, luminosity, radius_km, star_type, yerkes_class):
    """
    `(mass_loss_rate_kg_s, wind_velocity_m_s)` of a star of `mass` (kg), `luminosity` (W) and `radius_km`, by its type string (`star_type`, "G2V ...") and Yerkes class: evolved stars by the
    empirical fit of `constants.NDJ_*`, hot O and B dwarfs by a luminosity power law, and cool dwarfs and white dwarfs scaled from the Sun's coronal wind.
    """
    # --- 1. Get Fundamental Stellar Properties (from arguments) ---
    radius_m = _km_to_m(radius_km)
    lum_sol = luminosity / constants.SOLAR_LUMINOSITY
    mass_sol = _kg_to_msun(mass)
    radius_sol = radius_m / constants.SOLAR_RADIUS_M
    spectral_class = star_type[0].upper() if star_type else None

    # --- 2. Calculate Mass-Loss Rate (M-dot) and Wind Velocity (v_inf) ---
    # A single formula for mass loss is insufficient. We use a tiered system based on
    # the star's Yerkes luminosity class (evolutionary stage) and spectral type.
    escape_velocity = math.sqrt(constants.ESCAPE_VELOCITY_CONSTANT * constants.G * mass / radius_m)


    # TIER 1: Evolved stars -- giants (III), subgiants (IV), bright giants
    # (II), supergiants (IB/IAB/IA/IA+), and hypergiants (0). These were
    # previously three separately hand-fit tiers (a hypergiant power-law,
    # plus Reimers' Law with two different eta values for "giants" vs
    # "supergiants") stitched together at Yerkes-class boundaries -- the
    # exact structure that let the old hypergiant constant drift ~9
    # orders of magnitude out of calibration unnoticed. An
    # empirical mass-loss fit spanning this
    # whole regime, so one formula and one continuous curve now covers
    # all of it (see constants.py for the correction applied
    # above its known high-luminosity overestimation threshold).
    if yerkes_class in ["0", "IA+", "IA", "IAB", "IB", "II", "III", "IV"]:
        mass_loss_rate_smyr = (constants.NDJ_MASS_LOSS_COEFFICIENT
                               * (lum_sol**constants.NDJ_LUMINOSITY_EXPONENT)
                               * (mass_sol**constants.NDJ_MASS_EXPONENT)
                               * (radius_sol**constants.NDJ_RADIUS_EXPONENT))
        if lum_sol > constants.NDJ_HIGH_LUMINOSITY_THRESHOLD_LSUN:
            mass_loss_rate_smyr *= constants.NDJ_HIGH_LUMINOSITY_CORRECTION_FACTOR

        if yerkes_class == "0":
            # Hypergiant winds are an extreme multiple of escape velocity.
            wind_velocity = constants.HYPERGIANT_WIND_VELOCITY_FACTOR * escape_velocity
        else:
            # Wind velocity is a smaller fraction of escape velocity for these cooler giants/supergiants.
            wind_velocity = constants.GIANT_WIND_VELOCITY_FACTOR * escape_velocity

    # TIER 2a: Hot main-sequence dwarfs (O and B, Class V/VI)
    # Unlike cool dwarfs, O/B stars have powerful radiation-driven winds (line-driven,
    # via the same physical mechanism as hypergiants/supergiants, just far less extreme):
    # mass loss scales steeply with luminosity rather than the coronal/magnetic-activity
    # scaling used for cool dwarfs below. Without this tier, a hot O-type star would be
    # (incorrectly) modeled with a weaker wind than a cool, dim M dwarf.
    elif spectral_class in ("O", "B") and yerkes_class in ["V", "VI"]:
        mass_loss_rate_smyr = constants.OB_DWARF_MASS_LOSS_RATE_FACTOR * (lum_sol**constants.OB_DWARF_MASS_LOSS_RATE_EXPONENT)
        wind_velocity = constants.OB_DWARF_WIND_VELOCITY_FACTOR * escape_velocity

    # TIER 2b: Cool main-sequence dwarfs (A through M, Class V/VI) and White Dwarfs (VII)
    # For sun-like and cooler stars (plus stellar remnants), mass loss is very low. We
    # scale directly from the Sun's known properties for a stable, physically grounded result.
    else:
        # For cool main-sequence stars, mass loss is driven by magnetic activity, not
        # radiation pressure. The Reimers' law scaling (L*R/M) is inaccurate here.
        # We use a more modern, empirically-derived formula from [R4756],
        # which relates mass-loss rate to surface X-ray flux.
        # M-dot ~ R^2 * F_x^1.34, and F_x ~ L^-0.5
        # This results in a scaling relationship of M-dot ~ R^2 * L^-0.5
        # We normalize this to the Sun's known mass-loss rate.
        # The mass term is added to account for gravitational binding.
        scaling_factor = (radius_sol**constants.RADIUS_SOL_EXPONENT) * (lum_sol**constants.LUMINOSITY_SOL_EXPONENT) * (mass_sol**constants.MASS_SOL_EXPONENT)
        mass_loss_rate_smyr = constants.SUN_MASS_LOSS_RATE_SOLAR_MASS_PER_YEAR * scaling_factor # 2e-14 is the Sun's M-dot in M_sol/yr

        # Scale wind velocity based on the star's escape velocity relative to the Sun's.
        wind_velocity = constants.SOLAR_WIND_VELOCITY * (escape_velocity / constants.SOLAR_ESCAPE_VELOCITY)

    # --- 3. Convert and Calculate Final Radius ---

    # Convert M-dot from (solar masses/year) to (kg/s)
    mass_loss_rate_kgs = _msun_to_kg(mass_loss_rate_smyr) / constants.SECONDS_PER_YEAR
    return mass_loss_rate_kgs, wind_velocity
