# planetgen/physics/rogue_surface.py

"""
The surface conditions of a free-floating brown dwarf
=====================================================

A brown dwarf has no star and no solid surface, so what is stored for it is its own heat (Kelvin-Helmholtz cooling, `giant_internal_flux_w_m2`), the effective
temperature that heat gives and the temperature at the 1 bar level (`giant_temperature_at_k`). The rogue planets that systems throw out are built by the
planet chain with no star (GEN.1111, `generation/phenomena/starless_host.py`), not here; this module's older rocky branch, which drew an envelope, water and tides at a flat
rate, is gone.
"""

import math

from planetgen.physics import brown_dwarf, constants, mass_boundaries, relations, substellar_ceiling, unit_registry as units
from planetgen import tuning
from planetgen.util import draw

_kg_to_mjup = units.converter("kg", "Mjup")
_kg_to_msun = units.converter("kg", "Msun")

SURFACE_REGIMES = ("brown-dwarf",)
"""tuple: Every `surface_regime` this module gives: `brown-dwarf`, no surface, conditions at 1 bar."""

SURFACE_REGIME_LABELS = {"brown-dwarf": "No solid surface (1 bar level)"}
"""dict: A regime's display label."""

BAR_TO_PA = 1e5



def _sigma():
    return constants.STEFAN_BOLTZMANN_CONSTANT


def _giant_luminosity_w(mass_kg, age_gy):
    """A giant's or brown dwarf's internal luminosity, W: a power law in
    mass through Neptune and Jupiter below 1 Mjup, from Jupiter to the brown
    dwarf cooling track (`physics/brown_dwarf.py`) at 13 Mjup, then that
    track; below 13 Mjup fading as t^-1.3 (`ROGUE_GIANT_COOLING_TIME_EXPONENT`)."""
    pc = tuning
    phys = constants
    deuterium_jup = mass_boundaries.jupiter("deuterium_burning_limit")
    jupiter_l = pc.ROGUE_JUPITER_INTERNAL_FLUX_W_M2 * 4 * math.pi * (phys.JUPITER_RADIUS_KM * 1000.0) ** 2
    neptune_mass_kg, neptune_radius_m = 1.024e26, 2.4764e7
    neptune_l = _sigma() * pc.ROGUE_NEPTUNE_INTERNAL_TEMPERATURE_K ** 4 * 4 * math.pi * neptune_radius_m ** 2

    def brown_dwarf_l(m_kg, t_gy):
        return brown_dwarf.luminosity_sol(_kg_to_msun(m_kg), t_gy) * phys.SOLAR_LUMINOSITY

    m_jup = _kg_to_mjup(mass_kg)
    age_factor = (age_gy / pc.ROGUE_EARTH_REFERENCE_AGE_GY) ** pc.ROGUE_GIANT_COOLING_TIME_EXPONENT
    if m_jup >= deuterium_jup:
        return brown_dwarf_l(mass_kg, age_gy)
    if m_jup >= 1.0:
        at_13 = brown_dwarf_l(mass_boundaries.kg("deuterium_burning_limit"), pc.ROGUE_EARTH_REFERENCE_AGE_GY)
        exponent = math.log(at_13 / jupiter_l) / math.log(deuterium_jup)
        return jupiter_l * m_jup ** exponent * age_factor
    exponent = math.log(jupiter_l / neptune_l) / math.log(phys.JUPITER_MASS_TO_KG / neptune_mass_kg)
    return jupiter_l * m_jup ** exponent * age_factor


def giant_internal_flux_w_m2(mass_kg, radius_km, age_gy):
    """A giant's or brown dwarf's internal heat flux at its photosphere,
    W/m^2, capped at the flux of the hottest a body of its mass and age can be (`substellar_ceiling.cap_k`)."""
    flux = _giant_luminosity_w(mass_kg, age_gy) / (4 * math.pi * (radius_km * 1000.0) ** 2)
    cap_k = substellar_ceiling.cap_k(_kg_to_msun(mass_kg), age_gy)
    return min(flux, _sigma() * cap_k ** 4)


def effective_temperature_k(flux_w_m2):
    """Step 2: T_eff = ((F_int + sigma T_CMB^4) / sigma)^(1/4)."""
    cmb = constants.COSMIC_BACKGROUND_TEMPERATURE_K
    return ((flux_w_m2 + _sigma() * cmb ** 4) / _sigma()) ** 0.25


def surface_gravity_m_s2(mass_kg, radius_km):
    return relations.surface_gravity(mass_kg, radius_km * 1000.0)


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




def rogue_surface_conditions(mass_kg, radius_km, age_gy=None, rng=draw):
    """
    A free-floating brown dwarf's conditions.

    Args:
        mass_kg (float), radius_km (float): The body.
        age_gy (float, optional): Its age; drawn from `ROGUE_PLANET_AGE_RANGE_GY` when not given. Its radius depends on it, so `RoguePlanet` draws it first.
        rng: Anything with `uniform` (default: `draw`).

    Returns:
        dict: `age_gy`, `internal_heat_flux_w_m2`, `effective_temperature_k`, `surface_regime` (`brown-dwarf`), `surface_temperature_k` (at 1 bar),
            `surface_pressure_pa` (1 bar), `ice_shell_thickness_km`, `ocean_depth_km` and `hp_ice_km` (`None`), `has_liquid_water` (False)
            and `has_internal_heat` (True).
    """
    if age_gy is None:
        age_gy = rng.uniform(*tuning.ROGUE_PLANET_AGE_RANGE_GY)
    gravity = surface_gravity_m_s2(mass_kg, radius_km)
    flux = giant_internal_flux_w_m2(mass_kg, radius_km, age_gy)
    t_eff = effective_temperature_k(flux)
    return {
        "age_gy": age_gy, "ice_shell_thickness_km": None, "ocean_depth_km": None, "hp_ice_km": None, "has_liquid_water": False,
        "internal_heat_flux_w_m2": flux,
        "effective_temperature_k": t_eff,
        "surface_regime": "brown-dwarf",
        "surface_temperature_k": giant_temperature_at_k(t_eff, gravity, BAR_TO_PA),
        "surface_pressure_pa": BAR_TO_PA,
        "has_internal_heat": True,
    }


ROGUE_SURFACE_UNITS = {
    "age_gy": "Gyr", "internal_heat_flux_w_m2": "W_m2", "effective_temperature_k": "K", "surface_regime": "text",
    "surface_temperature_k": "K", "surface_pressure_pa": "Pa", "ice_shell_thickness_km": "km",
    "ocean_depth_km": "km", "has_liquid_water": "flag", "hp_ice_km": "km",
}
ROGUE_SURFACE_FIELDS = tuple(ROGUE_SURFACE_UNITS)
"""tuple: The `rogue_surface_conditions` keys stored on `RoguePlanet` and
in `rogue_planets` (schema v48); `has_internal_heat` is stored already."""
