# planetgen/physics/rogue_surface.py

"""
Rogue Planet Surface Conditions
===============================

A rogue planet has no star, so its surface is set by its own heat alone
(Boss's research, 2026-10-01; Stevenson 1999, Burrows et al. 2001,
Guillot 2005, Abbot & Switzer 2011, Turcotte & Schubert 2014). This
module follows that research's four steps:

1. Internal heat flux: radioactive decay (`radiogenic`, the host-less
   draw: a rogue's parent star is unknown) plus the heat the body was born with, its mantle's cooling
   (`interior_heat`, GEN.699), for a rocky rogue (plus tides if it kept a moon), Kelvin-Helmholtz
   cooling for a giant or brown dwarf.
2. Effective temperature: sigma T_eff^4 = F_int + sigma T_CMB^4.
3. Atmosphere: a hydrogen envelope stays only if the planet can hold it
   (Jeans escape parameter); any other air freezes out as surface frost.
4. Surface: down the envelope's adiabat to its base, through a
   conducting ice shell to a buried ocean, or the bare frozen surface at
   T_eff. A giant has no surface; its conditions are given at 1 bar.

`rogue_surface_conditions` returns plain numbers, so the same code fills
a new `RoguePlanet` and back-fills stored rows (`store._migrate_v47_to_v48`).
"""

import math

from planetgen.physics import brown_dwarf, constants, gases, hydrosphere, interior_heat, mass_boundaries, radiogenic, relations, saturation, substellar_ceiling, water_layer, unit_registry as units
from planetgen import tuning
from planetgen.util.random import log_uniform
from planetgen.util import draw

_kg_to_mearth = units.converter("kg", "Mearth")
_kg_to_mjup = units.converter("kg", "Mjup")
_kg_to_msun = units.converter("kg", "Msun")

SURFACE_REGIMES = (
    "bare-rock", "frozen-atmosphere", "ice-shell-ocean", "ice-world",
    "hydrogen-envelope", "hycean", "gas-giant", "brown-dwarf",
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
- `hycean`: such an envelope over a liquid ocean shallow enough to touch
  rock (no high-pressure ice beneath; GEN.88).
- `gas-giant` / `brown-dwarf`: no surface; conditions at 1 bar."""

SURFACE_REGIME_LABELS = {
    "bare-rock": "Bare frozen rock",
    "frozen-atmosphere": "Frozen atmosphere (surface frost)",
    "ice-shell-ocean": "Ice shell over a liquid ocean",
    "ice-world": "Frozen solid (ice to the rock)",
    "hydrogen-envelope": "Thick hydrogen envelope",
    "hycean": "Ocean under a hydrogen envelope",
    "gas-giant": "No solid surface (1 bar level)",
    "brown-dwarf": "No solid surface (1 bar level)",
}
"""dict: A regime's display label."""

BAR_TO_PA = 1e5
HYDROGEN_MOLECULE_KG = 2 * constants.HYDROGEN_ATOM_MASS_KG



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


def frozen_nitrogen_pressure_pa(temperature_k):
    """Nitrogen frost's vapour pressure from the one saturation function (`physics/saturation.py`): ~1-6 Pa at 37-40 K, Pluto's
    air; nothing at 20 K. Above nitrogen's critical temperature there is no saturation and the critical pressure is the most the air
    can hold as a frozen-atmosphere regime."""
    pressure = saturation.pressure_pa("N2", temperature_k)
    if math.isinf(pressure):
        return saturation.pressure_pa("N2", gases.gas("N2").critical_temperature_k)
    return pressure


def rogue_surface_conditions(mass_kg, radius_km, planet_type, mass_bin, has_moons, rng=draw, age_gy=None):
    """
    A rogue planet's surface conditions (the module docstring's four
    steps). `rng` supplies every draw (age, abundance, tides, envelope,
    water), so a seeded `draw.Stream` gives the same answer each time.

    Args:
        mass_kg (float), radius_km (float): The body.
        planet_type (str): `'t'` or `'g'`.
        mass_bin (str): Its `ROGUE_PLANET_MASS_BIN_CHOICES` bin.
        has_moons (bool): Whether it kept a moon (tidal heating, rocky only).
        rng: Anything with `uniform`, `random` and `gauss` (default: `draw`).
        age_gy (float, optional): Its age; drawn from `ROGUE_PLANET_AGE_RANGE_GY` when not given. A brown dwarf's
            radius depends on it, so `RoguePlanet` draws it first.

    Returns:
        dict: `age_gy`, `internal_heat_flux_w_m2`, `effective_temperature_k`,
            `surface_regime` (`SURFACE_REGIMES`), `surface_temperature_k`
            (a giant's at 1 bar), `surface_pressure_pa` (1 bar for a giant;
            0 for none), `ice_shell_thickness_km` and `ocean_depth_km`
            (`None` when there is no ice or no ocean), `has_liquid_water`
            and `has_internal_heat` (still geologically active, or a giant).
    """
    pc = tuning
    if age_gy is None:
        age_gy = rng.uniform(*pc.ROGUE_PLANET_AGE_RANGE_GY)
    gravity = surface_gravity_m_s2(mass_kg, radius_km)
    result = {
        "age_gy": age_gy, "ice_shell_thickness_km": None, "ocean_depth_km": None, "hp_ice_km": None,
        "has_liquid_water": False,
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
    content = radiogenic.draw_content(None, rng=rng)
    rock = radiogenic.Rock(pc.RADIOGENIC_SILICATE_MASS_FRACTION, 1.0 - pc.RADIOGENIC_SILICATE_MASS_FRACTION, 0.0,
                           radiogenic.BULK_SILICATE_EARTH_AL_MASS_FRACTION)
    flux = interior_heat.heat_flow(content, interior_heat.draw_interior(rng), mass_kg, radius_km, age_gy, rock).flux_w_m2
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
    had_air = _kg_to_mearth(mass_kg) >= pc.ROGUE_FROZEN_ATMOSPHERE_MIN_MASS_EARTH

    # Step 4: down the envelope's adiabat, or the bare surface at T_eff.
    if envelope_pa is not None:
        photosphere_pa = pc.ROGUE_ENVELOPE_PHOTOSPHERE_PRESSURE_BAR * BAR_TO_PA
        surface_k = adiabat_temperature_k(t_eff, photosphere_pa, max(envelope_pa, photosphere_pa))
        regime, pressure_pa = "hydrogen-envelope", envelope_pa
    else:
        surface_k = t_eff
        pressure_pa = frozen_nitrogen_pressure_pa(surface_k) if had_air else 0.0
        regime = "frozen-atmosphere" if had_air else "bare-rock"

    if water_fraction is not None:
        water_kg = water_fraction * mass_kg
        column = None
        if surface_k < hydrosphere.WATER_MELTING_K:
            shell_km = hydrosphere.pressure_melted_shell_km(flux, surface_k, gravity)
            column = water_layer.column(mass_kg, radius_km, water_kg, shell_km=shell_km,
                                        bottom_temperature_k=hydrosphere.WATER_MELTING_K, pressure_pa=pressure_pa)
        elif surface_k < hydrosphere.water_boiling_k(pressure_pa):
            column = water_layer.column(mass_kg, radius_km, water_kg, bottom_temperature_k=surface_k, pressure_pa=pressure_pa)
        if column is not None:
            result.update({"ice_shell_thickness_km": column.ice_shell_km or None, "ocean_depth_km": column.liquid_km or None,
                           "hp_ice_km": column.hp_ice_km or None, "has_liquid_water": bool(column.liquid_km)})
        if regime == "hydrogen-envelope":
            if result["has_liquid_water"] and not result["hp_ice_km"]:
                regime = "hycean"
        else:
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


ROGUE_SURFACE_UNITS = {
    "age_gy": "Gyr", "internal_heat_flux_w_m2": "W_m2", "effective_temperature_k": "K", "surface_regime": "text",
    "surface_temperature_k": "K", "surface_pressure_pa": "Pa", "ice_shell_thickness_km": "km",
    "ocean_depth_km": "km", "has_liquid_water": "flag", "hp_ice_km": "km",
}
ROGUE_SURFACE_FIELDS = tuple(ROGUE_SURFACE_UNITS)
"""tuple: The `rogue_surface_conditions` keys stored on `RoguePlanet` and
in `rogue_planets` (schema v48); `has_internal_heat` is stored already."""
