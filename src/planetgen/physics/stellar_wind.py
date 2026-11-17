# planetgen/physics/stellar_wind.py

"""
A star's wind (the tiered mass-loss rate and terminal speed), shared by the heliosphere (`Star.calculate_heliosphere`) and by the light of a black hole that a companion's wind feeds
(`dark_light`, GEN.1501): one model, so the star that sets the size of its bubble is the one that sets how much a black hole beside it takes in.
"""

import math

from planetgen import tuning
from planetgen.physics import constants
from planetgen.physics import unit_registry as units

WHITE_DWARF_CLASSES = ("VII", "D")
"""tuple: The Yerkes classes of a white dwarf (`D` is an alias for `VII`, `stellar_evolution.YERKES_CLASS_NAMES`)."""

_km_to_m = units.converter("km", "m")
_kg_to_msun = units.converter("kg", "Msun")
_msun_to_kg = units.converter("Msun", "kg")


def coronal_mass_loss_solar(radius_sol, lum_sol, log_lx_lbol, spectral_class):
    """
    The mass loss, in the Sun's, of a cool dwarf's coronal wind (GEN.1625). An A to K star's is its surface area times its X-ray flux over the Sun's to
    `tuning.CORONAL_WIND_X_RAY_EXPONENT` (1.29 [R5303], 1.34 [R5304]), the flux no more than `tuning.CORONAL_WIND_X_RAY_FLUX_CAP_ERG_CM2_S` (the dividing line at
    which the winds of more active stars turn weak [R5303]). An M dwarf's is `tuning.M_DWARF_WIND_PER_AREA_SOLAR` times its surface area: the measured winds
    of M dwarfs do not follow their X-ray flux (two orders of magnitude of scatter at a flux, [R5302]).

    Raises:
        ValueError: An A to K star without `log_lx_lbol`.
    """
    area = radius_sol ** 2
    if spectral_class == "M":
        return tuning.M_DWARF_WIND_PER_AREA_SOLAR * area
    if log_lx_lbol is None:
        raise ValueError("a cool dwarf's coronal wind needs its X-ray luminosity (log_lx_lbol)")
    flux_over_solar = 10 ** log_lx_lbol * lum_sol * constants.SOLAR_LUMINOSITY / constants.SOLAR_X_RAY_LUMINOSITY_W / area
    cap = (tuning.CORONAL_WIND_X_RAY_FLUX_CAP_ERG_CM2_S * constants.ERG_CM2_S_TO_W_M2
           * 4.0 * math.pi * constants.SOLAR_RADIUS_M ** 2 / constants.SOLAR_X_RAY_LUMINOSITY_W)
    return area * min(flux_over_solar, cap) ** tuning.CORONAL_WIND_X_RAY_EXPONENT


def wind(mass, luminosity, radius_km, star_type, yerkes_class, log_lx_lbol=None):
    """
    `(mass_loss_rate_kg_s, wind_velocity_m_s)` of a star of `mass` (kg), `luminosity` (W) and `radius_km`, by its type string (`star_type`, "G2V ...") and Yerkes class: evolved stars by the
    empirical fit of `constants.NDJ_*`, hot O and B dwarfs by a luminosity power law, cool dwarfs by their coronal wind (`coronal_mass_loss_solar`, from the star's X-ray luminosity
    `log_lx_lbol`, log10 of its share of the luminosity, which a star that has a corona must give) and white dwarfs none.
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

    # TIER 2b: a white dwarf (VII) has no corona and no wind to speak of (GEN.1625: it used to get the cool dwarfs' scaled wind).
    elif yerkes_class in WHITE_DWARF_CLASSES:
        mass_loss_rate_smyr = 0.0
        wind_velocity = 0.0

    # TIER 2c: cool main-sequence dwarfs (A through M, Class V/VI) and the stars before the main sequence: a coronal wind,
    # from the X-ray flux at the surface for A to K stars and a measured mass loss per area for M dwarfs (`coronal_mass_loss_solar`).
    else:
        mass_loss_rate_smyr = constants.SUN_MASS_LOSS_RATE_SOLAR_MASS_PER_YEAR * coronal_mass_loss_solar(
            radius_sol, lum_sol, log_lx_lbol, spectral_class)
        # The wind speed of the lower main sequence is near the surface escape speed [R5302], here the Sun's 400 km/s scaled by the escape speed.
        wind_velocity = constants.SOLAR_WIND_VELOCITY * (escape_velocity / constants.SOLAR_ESCAPE_VELOCITY)

    # --- 3. Convert and Calculate Final Radius ---

    # Convert M-dot from (solar masses/year) to (kg/s)
    mass_loss_rate_kgs = _msun_to_kg(mass_loss_rate_smyr) / constants.SECONDS_PER_YEAR
    return mass_loss_rate_kgs, wind_velocity
