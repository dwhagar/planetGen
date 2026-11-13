# planetgen/physics/air_final.py

"""
The air a body ends with (GEN.1657, layered fix J)
==================================================

`planets.calculate_atmospheric_conditions` used to solve the air (`surface_air.solve`) and then correct the solver's input and output in three places of its own: it chose the steam branch's
start temperature (`past_runaway_with_water`), replaced the air with an oddity's (`air_branches`, the rows of lever L7) and re-solved the albedo of a body left with no air
(`surface_air.airless_state`). This module takes those three into one step, `make`, which returns the air the body keeps as one `FinalAir`: the caller stores it and does not branch on
how it came out.

Nothing is drawn here that the old order did not draw, in the same order, so the output is the one it was (the golden is the proof).
"""

import dataclasses
import os
from typing import Optional

from planetgen import tuning
from planetgen.physics import air_branches, climate_draws, constants, oddities, orbits, surface_air
from planetgen.physics import unit_registry as units
from planetgen.util import log

CAP_TOLERANCE = 1e-9
"""float: How far above the stored temperature the air's uncapped warming must reach to count as held at the ionization ceiling (`surface_air.Air.ceiling_k`, GEN.1236, GEN.1580)."""

_m_to_km = units.converter("m", "km")


@dataclasses.dataclass(frozen=True)
class FinalAir:
    """
    The air a body keeps and the numbers the body stores from it: `air` (the solver's, or the oddity's, `None` when the body ended with no air), the surface temperature, the Bond albedo and
    greenhouse factor that make it, the pressure, the mean molar mass, the density, the scale height (km) and the partial pressures. A body with no air has pressure 0, no density
    or scale height and no partials.
    """

    air: object
    temperature_k: float
    albedo: float
    greenhouse_factor: float
    pressure_pa: float
    molar_density_kg_mol: Optional[float]
    density_kg_m3: Optional[float]
    scale_height_km: Optional[float]
    partials_pa: dict


def past_runaway_with_water(planet, starlight_w_m2, surface_w):
    """
    Whether the body holds surface water under starlight past its runaway greenhouse limit (GEN.1238): an effective flux (in the Sun's at the
    Earth, `starlight_w_m2` over `RUNAWAY_FLUX_REFERENCE_W_M2`) above `orbits.runaway_greenhouse_flux` of the star's temperature and the body's mass,
    the limit the water-loss history already follows (`runaway_history`, GEN.639). Past it no ocean can stand, but the air's fixed point has two
    answers there (a cool air with a few bar of vapor, whose greenhouse the table gives no limit to, and the steam air that holds the water),
    and which one it finds depends on where it starts, so a body past the limit starts hot and finds the steam one; a body below it, or with no
    water, starts from the airless temperature as before. No draw.
    """
    if not surface_w or planet.body_type == "g" or not getattr(planet.star, "temperature", None) or not planet.mass:
        return False
    effective = starlight_w_m2 / tuning.RUNAWAY_FLUX_REFERENCE_W_M2
    return effective > orbits.runaway_greenhouse_flux(planet.star.temperature, units.converter("kg", "Mearth")(planet.mass))


def make(planet, setting, no_air_temperature_k, starlight_w_m2, distance_au, albedo, giant, surface_w):
    """
    The air of `planet` as it ends: solved from the start temperature the body's water gives (`past_runaway_with_water`), replaced by the air of an oddity of lever L7 when a row names the body
    (`air_branches`), or the airless state when it keeps none (`surface_air.airless_state`, the ice-albedo law's own answer). Keeps on the planet, as before and not stored, the replay of the
    late accretion (`_air_replay`, GEN.1493), the retention cap of the first solve (`_retention_cap_pa`, GEN.1340) and the carbon sink of the air kept (`_co2_sink`, GEN.1309 A).

    Args:
        planet (Planet): The planet or moon.
        setting (surface_air.Setting): The solver's input (`surface_air.setting_of`).
        no_air_temperature_k (float): The temperature starlight alone gives it, which the solver starts from.
        starlight_w_m2 (float): The light it gets at the distance the air is made for.
        distance_au (float): That distance.
        albedo (float): The Bond albedo drawn before the air.
        giant (bool): Whether it has no surface (the air is solved at its 1 bar level).
        surface_w (float or None): The surface water drawn before the air.

    Returns:
        FinalAir: What the body keeps.
    """
    start_k = no_air_temperature_k
    if past_runaway_with_water(planet, starlight_w_m2, surface_w):
        start_k = max(start_k, tuning.RUNAWAY_START_TEMPERATURE_K)    # GEN.1238: the steam branch, not the cool one (see `past_runaway_with_water`)
    air = surface_air.solve(setting, start_k)
    if not planet.is_moon and not giant and air is not None and air.partials_pa:       # GEN.1493: the late accretion this air was made with and what it is small or large against, for `remake_air`
        planet._air_replay = (*planet._air_replay[:3], (setting, air.late_state, air.late, air.degassed))
    planet._retention_cap_pa = None if air is None else air.retention_cap_pa     # GEN.1340: kept, not stored; the life stage keeps its oxygen inside it
    air = _air_oddities(planet, setting, air, start_k, starlight_w_m2, distance_au, albedo, giant)
    planet._co2_sink = None if air is None else air.co2_sink                      # GEN.1309 A: kept, not stored; the life stage balances the CO2 again at the temperature the oxygen rise leaves
    planet._air_haze = None if air is None or setting.haze_draws is None else (setting.haze_draws, setting.star_temperature_k, setting.far_uv, setting.internal_flux_w_m2)    # GEN.1589: kept, not stored; the haze of the air kept, for the life stage
    if air is None or not air.partials_pa:
        temperature_k, albedo_kept = no_air_temperature_k, albedo
        if setting.albedo_law is not None:
            # GEN.972: with no air the ice the cold leaves is the only feedback; the temperature and the albedo that keep each other
            temperature_k, albedo_kept = surface_air.airless_state(setting.albedo_law, starlight_w_m2, no_air_temperature_k)
        return FinalAir(air=air, temperature_k=temperature_k, albedo=albedo_kept, greenhouse_factor=0.0, pressure_pa=0.0, molar_density_kg_mol=None, density_kg_m3=None,
                        scale_height_km=None, partials_pa={})
    if air.at_ceiling:
        _warn_at_ceiling(planet, air, distance_au)
    factor = air.greenhouse_factor
    albedo_kept = albedo if air.albedo is None else air.albedo      # a body with ice has the albedo of its last step (GEN.972)
    # GEN.1236: the fixed point holds the surface at the ionization ceiling (`ionization.ceiling_k`, GEN.1580; a guard no generated body should reach, the equilibrium of the gas does
    # not converge far above it), but the warming it keeps is the uncapped one, so a body at the cap kept a greenhouse factor that gave 2370 K
    # or 2600 K; its stored G is taken from the stored temperature too, so the three numbers agree.
    capped = air.warming is not None and air.warming.surface_temperature_k > air.temperature_k * (1.0 + CAP_TOLERANCE)
    if giant or capped or air.knee:        # (GEN.1344: the knee of an unsettled air is not the temperature its last step found)
        # GEN.412 reads the temperature as starlight alone times (1 + G)^(1/4); a giant also radiates its own heat, which the air's factor
        # leaves out of T_e, so the stored G is taken against the starlight alone (GEN.1030) and the pair reproduces the temperature
        factor = climate_draws.greenhouse_from_temperature(air.temperature_k, albedo_kept, starlight_w_m2)
        factor = air.greenhouse_factor if factor is None else max(factor, 0.0)
    return FinalAir(air=air, temperature_k=air.temperature_k, albedo=albedo_kept, greenhouse_factor=factor, pressure_pa=air.pressure_pa,
                    molar_density_kg_mol=air.state.mean_molar_mass_kg_mol, density_kg_m3=air.state.density_kg_m3, scale_height_km=_m_to_km(air.state.scale_height_m),
                    partials_pa=air.partials_pa)


class CeilingReached(RuntimeError):
    """A body's air asked for more than the ionization ceiling while `CEILING_IS_ERROR_VARIABLE` was set (the test suite sets it, GEN.1580)."""


CEILING_IS_ERROR_VARIABLE = "PLANETGEN_CEILING_IS_ERROR"
"""str: The environment variable that makes a body held at the ionization ceiling an error (`CeilingReached`) and not only a log line. Off in a run, set in the test suite
(`tests/conftest.py`, workers included), so every seeded census and fuzz run that makes a body at the ceiling fails and the body is explained (Boss, 2026-11-12, GEN.1580)."""


def _warn_at_ceiling(planet, air, distance_au):
    """
    GEN.1580: one normal-log line per body whose air asked for more than the ionization ceiling (`ionization.ceiling_k`) and was held at it: the body, the
    uncapped surface the greenhouse gave, the ceiling and the air. The ceiling is the temperature the main gas is half ionized, so a body here is a mistake of
    the model to explain, never a state to keep (the test suite sets `CEILING_IS_ERROR_VARIABLE`, so a seeded sweep that makes one fails).
    """
    if getattr(planet, "_ceiling_warned", False):
        return
    planet._ceiling_warned = True
    name = getattr(planet, "name", None) or "a body not yet named"
    composition = ", ".join(f"{formula} {pressure:.4g} Pa" for formula, pressure in sorted(air.partials_pa.items(), key=lambda item: -item[1])[:6])
    message = (f"Ionization ceiling: {name} ({planet.mass / constants.EARTH_MASS_TO_KG:.4g} Earth masses, {distance_au:.4g} au) asked for {air.warming.surface_temperature_k:.0f} K "
               f"and was held at {air.ceiling_k:.0f} K, the half-ionization temperature of its main gas ({composition}).")
    log.normal(message)
    if os.environ.get(CEILING_IS_ERROR_VARIABLE) == "1":
        raise CeilingReached(message)


def _air_oddities(planet, setting, air, start_k, starlight_w_m2, distance_au, albedo, giant):
    """
    The rows of the air lever (L7, GEN.1051: W19 the nightside collapse of a locked planet's air, W20 the other branch of a bistable climate) on the air just solved: they may replace it
    with the airless state (`None`) or with the air of the other branch (`physics/air_branches.py`). The rows are rolled again each time the air is made, from their own streams.
    """
    air_branches.clear(planet)
    if giant or planet.is_moon or air is None:
        return air
    planet._air_context = air_branches.Context(setting, air, start_k, starlight_w_m2, distance_au, albedo)
    oddities.drop(planet, ("L7",))
    oddities.apply(planet, levers=("L7",))
    choice = getattr(planet, air_branches.CHOICE_ATTRIBUTE, None)
    air_branches.clear(planet)
    if choice == "lost":
        return None
    return air if choice is None else choice
