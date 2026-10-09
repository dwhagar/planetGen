# planetgen/physics/rogue_surface.py

"""
Rogue Planet Surface Conditions
===============================

A rogue planet has no star, so its surface is set by its own heat alone
(Boss's research, 2026-10-01; Stevenson 1999, Burrows et al. 2001,
Guillot 2005, Abbot & Switzer 2011, Turcotte & Schubert 2014). This
module follows that research's four steps:

1. Internal heat flux: radioactive decay plus leftover formation heat
   for a rocky rogue (plus tides if it kept a moon), Kelvin-Helmholtz
   cooling for a giant or brown dwarf.
2. Effective temperature: sigma T_eff^4 = F_int + sigma T_CMB^4.
3. Atmosphere: a hydrogen envelope stays only if the planet can hold it
   (Jeans escape parameter); any other air freezes out as surface frost.
4. Surface: down the envelope's adiabat to its base, through a
   conducting ice shell to a buried ocean, or the bare frozen surface at
   T_eff. A giant has no surface; its conditions are given at 1 bar.

`rogue_surface_conditions` returns plain numbers, so the same code fills
a new `RoguePlanet` and back-fills stored rows (`_db._migrate_v47_to_v48`).
"""

import math

from planetgen.physics import constants
from planetgen import tuning
from planetgen.util.random import log_uniform
from planetgen.util import draw

SURFACE_REGIMES = (
    "bare-rock", "frozen-atmosphere", "ice-shell-ocean", "ice-world",
    "hydrogen-envelope", "gas-giant", "brown-dwarf",
)
"""tuple: Every `surface_regime`:
- `bare-rock`: a small, dry world with no air left, at T_eff.
- `frozen-atmosphere`: its air has frozen onto the ground as frost; a
  trace of nitrogen sublimates from it.
- `ice-shell-ocean`: a water world whose ice crust insulates a liquid
  ocean beneath (Abbot & Switzer's "Steppenwolf").
- `ice-world`: a water world too cold, or with too little water, for any
  of it to stay liquid: ice down to the rock.
- `hydrogen-envelope`: a thick hydrogen envelope, opaque from collision-
  induced absorption, traps the heat so the ground is warm (Stevenson).
- `gas-giant` / `brown-dwarf`: no surface; conditions at 1 bar."""

SURFACE_REGIME_LABELS = {
    "bare-rock": "Bare frozen rock",
    "frozen-atmosphere": "Frozen atmosphere (surface frost)",
    "ice-shell-ocean": "Ice shell over a liquid ocean",
    "ice-world": "Frozen solid (ice to the rock)",
    "hydrogen-envelope": "Thick hydrogen envelope",
    "gas-giant": "No solid surface (1 bar level)",
    "brown-dwarf": "No solid surface (1 bar level)",
}
"""dict: A regime's display label."""

WATER_MELTING_K = 273.15
WATER_CRITICAL_K = 647.1
WATER_BOILING_K = 373.15
WATER_VAPORIZATION_J_MOL = 40.7e3
NITROGEN_TRIPLE_POINT_K = 63.15
NITROGEN_TRIPLE_POINT_PA = 12.5e3
NITROGEN_SUBLIMATION_J_MOL = 6.9e3
GAS_CONSTANT_J_MOL_K = 8.314
BAR_TO_PA = 1e5
HYDROGEN_MOLECULE_KG = 2 * constants.HYDROGEN_ATOM_MASS_KG
LIQUID_WATER_DENSITY_KG_M3 = 1000.0



def _sigma():
    return constants.STEFAN_BOLTZMANN_CONSTANT


def radiogenic_flux_w_m2(mass_kg, radius_km, age_gy, abundance=1.0):
    """
    Radioactive heat reaching the surface, W/m^2: Boss's
    F = (M_mantle / 4 pi R^2) sum C_i H_i e^(-lambda_i t), with Earth's
    mantle concentrations (`ROGUE_RADIOGENIC_ISOTOPES`) as they stood at
    `age_gy`, times `abundance`.
    """
    pc = tuning
    mantle_kg = mass_kg * pc.ROGUE_MANTLE_MASS_FRACTION
    area_m2 = 4 * math.pi * (radius_km * 1000.0) ** 2
    heat_w_kg = 0.0
    for _name, now_kg_kg, heat_w_per_kg, half_life_gy in pc.ROGUE_RADIOGENIC_ISOTOPES:
        decay = math.log(2) / half_life_gy
        heat_w_kg += now_kg_kg * heat_w_per_kg * math.exp(decay * (pc.ROGUE_EARTH_REFERENCE_AGE_GY - age_gy))
    return abundance * mantle_kg * heat_w_kg / area_m2


def _giant_luminosity_w(mass_kg, age_gy):
    """A giant's or brown dwarf's internal luminosity, W: a power law in
    mass through Neptune and Jupiter below 1 Mjup, from Jupiter to Burrows
    & Liebert's brown dwarf law at 13 Mjup, then that law; all fading as
    t^-1.3 (`ROGUE_GIANT_COOLING_TIME_EXPONENT`)."""
    pc = tuning
    phys = constants
    jupiter_l = pc.ROGUE_JUPITER_INTERNAL_FLUX_W_M2 * 4 * math.pi * (phys.JUPITER_RADIUS_KM * 1000.0) ** 2
    neptune_mass_kg, neptune_radius_m = 1.024e26, 2.4764e7
    neptune_l = _sigma() * pc.ROGUE_NEPTUNE_INTERNAL_TEMPERATURE_K ** 4 * 4 * math.pi * neptune_radius_m ** 2

    def brown_dwarf_l(m_kg, t_gy):
        m_solar = m_kg / phys.SOLAR_MASS_TO_KG
        return (pc.ROGUE_BROWN_DWARF_LUMINOSITY_SOL * phys.SOLAR_LUMINOSITY
                * t_gy ** pc.ROGUE_GIANT_COOLING_TIME_EXPONENT * (m_solar / 0.05) ** pc.ROGUE_BROWN_DWARF_MASS_EXPONENT)

    m_jup = mass_kg / phys.JUPITER_MASS_TO_KG
    age_factor = (age_gy / pc.ROGUE_EARTH_REFERENCE_AGE_GY) ** pc.ROGUE_GIANT_COOLING_TIME_EXPONENT
    if m_jup >= 13.0:
        return brown_dwarf_l(mass_kg, age_gy)
    if m_jup >= 1.0:
        at_13 = brown_dwarf_l(13.0 * phys.JUPITER_MASS_TO_KG, pc.ROGUE_EARTH_REFERENCE_AGE_GY)
        exponent = math.log(at_13 / jupiter_l) / math.log(13.0)
        return jupiter_l * m_jup ** exponent * age_factor
    exponent = math.log(jupiter_l / neptune_l) / math.log(phys.JUPITER_MASS_TO_KG / neptune_mass_kg)
    return jupiter_l * m_jup ** exponent * age_factor


def giant_internal_flux_w_m2(mass_kg, radius_km, age_gy):
    """A giant's or brown dwarf's internal heat flux at its photosphere,
    W/m^2, capped at `ROGUE_MAX_EFFECTIVE_TEMPERATURE_K`'s."""
    flux = _giant_luminosity_w(mass_kg, age_gy) / (4 * math.pi * (radius_km * 1000.0) ** 2)
    return min(flux, _sigma() * tuning.ROGUE_MAX_EFFECTIVE_TEMPERATURE_K ** 4)


def effective_temperature_k(flux_w_m2):
    """Step 2: T_eff = ((F_int + sigma T_CMB^4) / sigma)^(1/4)."""
    cmb = constants.COSMIC_BACKGROUND_TEMPERATURE_K
    return ((flux_w_m2 + _sigma() * cmb ** 4) / _sigma()) ** 0.25


def surface_gravity_m_s2(mass_kg, radius_km):
    return constants.G * mass_kg / (radius_km * 1000.0) ** 2


def escape_parameter(mass_kg, radius_km, molecule_kg, temperature_k):
    """Jeans escape parameter lambda = G M m / (k T R)."""
    return (constants.G * mass_kg * molecule_kg
            / (constants.BOLTZMANN * temperature_k * radius_km * 1000.0))


def adiabat_temperature_k(top_temperature_k, top_pressure_pa, pressure_pa):
    """Down a dry H2/He adiabat: T = T_top (P / P_top)^((gamma - 1) / gamma)."""
    gamma = constants.ADIABATIC_INDEX_H2_HE
    return top_temperature_k * (pressure_pa / top_pressure_pa) ** ((gamma - 1) / gamma)


def giant_photosphere_pressure_pa(gravity_m_s2):
    """Where tau = 2/3: P = (2/3) g / kappa_R, with one Rosseland opacity
    for every giant fixed by Jupiter (`ROGUE_GIANT_PHOTOSPHERE_PRESSURE_BAR`)."""
    jupiter_g = surface_gravity_m_s2(constants.JUPITER_MASS_TO_KG, constants.JUPITER_RADIUS_KM)
    kappa = (2 / 3) * jupiter_g / (tuning.ROGUE_GIANT_PHOTOSPHERE_PRESSURE_BAR * BAR_TO_PA)
    return (2 / 3) * gravity_m_s2 / kappa


def giant_temperature_at_k(t_eff_k, gravity_m_s2, pressure_pa):
    """A giant's temperature at `pressure_pa`: the gray radiative profile
    T^4 = (3/4) T_eff^4 (tau + 2/3) above the photosphere, the adiabat
    below it."""
    photosphere_pa = giant_photosphere_pressure_pa(gravity_m_s2)
    if pressure_pa <= photosphere_pa:
        tau = (2 / 3) * pressure_pa / photosphere_pa
        return t_eff_k * (0.75 * (tau + 2 / 3)) ** 0.25
    return adiabat_temperature_k(t_eff_k, photosphere_pa, pressure_pa)


def nitrogen_vapor_pressure_pa(temperature_k):
    """Nitrogen frost's vapor pressure (Clausius-Clapeyron from N2's triple
    point): ~1-6 Pa at 37-40 K, Pluto's air; nothing at 20 K."""
    if temperature_k >= NITROGEN_TRIPLE_POINT_K:
        return NITROGEN_TRIPLE_POINT_PA
    return NITROGEN_TRIPLE_POINT_PA * math.exp(
        -NITROGEN_SUBLIMATION_J_MOL / GAS_CONSTANT_J_MOL_K * (1 / temperature_k - 1 / NITROGEN_TRIPLE_POINT_K))


def water_boiling_k(pressure_pa):
    """Water's boiling point at `pressure_pa` (Clausius-Clapeyron), capped at
    its critical point: above that there is no liquid, only a supercritical
    fluid."""
    if pressure_pa <= 0:
        return WATER_MELTING_K
    inverse = 1 / WATER_BOILING_K - GAS_CONSTANT_J_MOL_K * math.log(pressure_pa / 101325.0) / WATER_VAPORIZATION_J_MOL
    return WATER_CRITICAL_K if inverse <= 0 else min(1 / inverse, WATER_CRITICAL_K)


def ice_shell_thickness_km(flux_w_m2, top_temperature_k, base_temperature_k=WATER_MELTING_K):
    """Step 4's conducting lid: with k = A / T, Fourier's law integrates to
    D = (A / F) ln(T_base / T_top). 0 when the top is already at melting."""
    if top_temperature_k >= base_temperature_k:
        return 0.0
    return (tuning.ROGUE_ICE_CONDUCTIVITY_A_W_M / flux_w_m2
            * math.log(base_temperature_k / top_temperature_k) / 1000.0)


def water_layer_depth_km(mass_kg, radius_km, water_mass_fraction):
    """How deep `water_mass_fraction` of the planet's mass would lie spread
    over its surface, at liquid water's density (ignores compression)."""
    area_m2 = 4 * math.pi * (radius_km * 1000.0) ** 2
    return water_mass_fraction * mass_kg / (LIQUID_WATER_DENSITY_KG_M3 * area_m2) / 1000.0


def rogue_surface_conditions(mass_kg, radius_km, planet_type, mass_bin, has_moons, rng=draw):
    """
    A rogue planet's surface conditions (the module docstring's four
    steps). `rng` supplies every draw (age, abundance, tides, envelope,
    water), so a seeded `draw.Stream` gives the same answer each time.

    Args:
        mass_kg (float), radius_km (float): The body.
        planet_type (str): `'t'` or `'g'`.
        mass_bin (str): Its `ROGUE_PLANET_MASS_BIN_CHOICES` bin.
        has_moons (bool): Whether it kept a moon (tidal heating, rocky only).
        rng: Anything with `uniform` and `random` (default: `random`).

    Returns:
        dict: `age_gy`, `internal_heat_flux_w_m2`, `effective_temperature_k`,
            `surface_regime` (`SURFACE_REGIMES`), `surface_temperature_k`
            (a giant's at 1 bar), `surface_pressure_pa` (1 bar for a giant;
            0 for none), `ice_shell_thickness_km` and `ocean_depth_km`
            (`None` when there is no ice or no ocean), `has_liquid_water`
            and `has_internal_heat` (still geologically active, or a giant).
    """
    pc = tuning
    age_gy = rng.uniform(*pc.ROGUE_PLANET_AGE_RANGE_GY)
    gravity = surface_gravity_m_s2(mass_kg, radius_km)
    result = {
        "age_gy": age_gy, "ice_shell_thickness_km": None, "ocean_depth_km": None, "has_liquid_water": False,
    }

    if planet_type == 'g':
        flux = giant_internal_flux_w_m2(mass_kg, radius_km, age_gy)
        t_eff = effective_temperature_k(flux)
        result.update({
            "internal_heat_flux_w_m2": flux,
            "effective_temperature_k": t_eff,
            "surface_regime": "brown-dwarf" if mass_bin == "brown-dwarf" else "gas-giant",
            "surface_temperature_k": giant_temperature_at_k(t_eff, gravity, BAR_TO_PA),
            "surface_pressure_pa": BAR_TO_PA,
            "has_internal_heat": True,
        })
        return result

    # Step 1: radioactive heat, leftover heat in Earth's proportion, tides.
    abundance = log_uniform(*pc.ROGUE_RADIOGENIC_ABUNDANCE_RANGE, rng=rng)
    flux = radiogenic_flux_w_m2(mass_kg, radius_km, age_gy, abundance) / pc.ROGUE_UREY_RATIO
    if has_moons:
        flux += log_uniform(*pc.ROGUE_TIDAL_FLUX_RANGE_W_M2, rng=rng)
    # Step 2.
    t_eff = effective_temperature_k(flux)

    # Step 3: a hydrogen envelope, if it formed one and can hold it.
    envelope_pa = None
    if rng.random() < pc.ROGUE_HYDROGEN_ENVELOPE_CHANCE.get(mass_bin, 0.0):
        envelope_pa = log_uniform(*pc.ROGUE_HYDROGEN_ENVELOPE_PRESSURE_BAR[mass_bin], rng=rng) * BAR_TO_PA
        if escape_parameter(mass_kg, radius_km, HYDROGEN_MOLECULE_KG, t_eff) < pc.ROGUE_MIN_ESCAPE_PARAMETER:
            envelope_pa = None
    water_fraction = None
    if rng.random() < pc.ROGUE_WATER_RICH_CHANCE.get(mass_bin, 0.0):
        water_fraction = log_uniform(*pc.ROGUE_WATER_MASS_FRACTION_RANGE, rng=rng)
    had_air = mass_kg / constants.EARTH_MASS_TO_KG >= pc.ROGUE_FROZEN_ATMOSPHERE_MIN_MASS_EARTH

    # Step 4: down the envelope's adiabat, or the bare surface at T_eff.
    if envelope_pa is not None:
        photosphere_pa = pc.ROGUE_ENVELOPE_PHOTOSPHERE_PRESSURE_BAR * BAR_TO_PA
        surface_k = adiabat_temperature_k(t_eff, photosphere_pa, max(envelope_pa, photosphere_pa))
        regime, pressure_pa = "hydrogen-envelope", envelope_pa
    else:
        surface_k = t_eff
        pressure_pa = nitrogen_vapor_pressure_pa(surface_k) if had_air else 0.0
        regime = "frozen-atmosphere" if had_air else "bare-rock"

    if water_fraction is not None:
        depth_km = water_layer_depth_km(mass_kg, radius_km, water_fraction)
        if surface_k < WATER_MELTING_K:
            shell_km = ice_shell_thickness_km(flux, surface_k)
            if shell_km >= depth_km:
                result["ice_shell_thickness_km"] = depth_km
            else:
                result.update({"ice_shell_thickness_km": shell_km, "ocean_depth_km": depth_km - shell_km,
                               "has_liquid_water": True})
        elif surface_k < water_boiling_k(pressure_pa):
            result.update({"ocean_depth_km": depth_km, "has_liquid_water": True})
        if regime != "hydrogen-envelope":
            regime = "ice-shell-ocean" if result["has_liquid_water"] else "ice-world"

    result.update({
        "internal_heat_flux_w_m2": flux,
        "effective_temperature_k": t_eff,
        "surface_regime": regime,
        "surface_temperature_k": max(surface_k, constants.COSMIC_BACKGROUND_TEMPERATURE_K),
        "surface_pressure_pa": pressure_pa,
        "has_internal_heat": flux >= pc.ROGUE_ACTIVE_HEAT_FLUX_W_M2,
    })
    return result


ROGUE_SURFACE_FIELDS = (
    "age_gy", "internal_heat_flux_w_m2", "effective_temperature_k", "surface_regime",
    "surface_temperature_k", "surface_pressure_pa", "ice_shell_thickness_km", "ocean_depth_km",
    "has_liquid_water",
)
"""tuple: The `rogue_surface_conditions` keys stored on `RoguePlanet` and
in `rogue_planets` (schema v48); `has_internal_heat` is stored already."""
