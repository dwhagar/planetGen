# planetgen/physics/rogue_habitability.py

"""
The habitability score of a rogue planet.

A rogue planet has no star, so its score is `habitability_world.compute` read from its surface regime
(`rogue_surface.rogue_surface_conditions`), the way a star system's body is read from its air and water:

* **Surface.** The stored `surface_temperature_k` and `surface_pressure_pa`; the gas is the regime's: all hydrogen under a
  `hydrogen-envelope` (or a `hycean` ocean), all nitrogen under a `frozen-atmosphere`, none on `bare-rock`.
* **Water.** `has_liquid_water` and the ice lid give the state: `ice-covered ocean` when the liquid lies under ice, `hycean`
  when it lies under a hydrogen envelope, `ice` for ice with no liquid, else `dry`. The water mass fraction is worked back
  from the water column (liquid, lid and high-pressure ice) the body keeps. An ocean takes the middle of its class's pH and
  water activity (the class is `ice-sealed` under high-pressure ice, else `neutral`), so nothing is drawn.
* **Energy.** Its own heat only: 1 percent of `internal_heat_flux_w_m2` (no starlight).
* **Dose.** The galactic cosmic rays under the column of air with no magnetic field and no star (`radiation.gcr_dose_msv_yr`,
  the open-space multiplier of 1), plus the crust's own gamma dose for its age from the bulk silicate Earth's potassium, thorium and uranium (no felsic crust, and no
  radon without land under air). An ocean under ice sees the ground's dose only.
* **Life.** None: a rogue planet has no life draw, so its four life scores are 0 and its food is the lowest rung. The score
  says what the surface would allow a visitor, not what lives there.
* **Equipment.** The survival kit a visitor needs (`equipment.compute`); the regional temperatures are not kept.

A gas giant or brown dwarf has no surface and scores nothing (every field `None`).
"""

from planetgen.physics import atmosphere, constants, gases, habitability_world, hydrosphere, outside_hours, radiation, radiogenic, regional_temperature, rogue_surface, surface_composition

ROGUE_FIELDS = (*(name for name in habitability_world.SCORE_FIELDS
                  if name not in (*regional_temperature.FIELDS, "outside_hours", *outside_hours.FIELDS)), "energy_flux_w_m2")
"""tuple: What `generate` sets on a rogue planet, as stored in `rogue_planets`: the score of a star system's body without the
regional temperatures and the hours outside (a rogue planet has no star, so no latitude bands, day side or daylight) and its stored energy flux (there is
no star, so no flux variation)."""

ROGUE_UNITS = {
    **{name: "index" for name in ROGUE_FIELDS if name.startswith(("tier_", "life_", "eq_"))},
    **{name: "ratio" for name in ROGUE_FIELDS if name.startswith("phi4") or name in ("l_solv", "l_chem", "l_ener", "l_rad")},
    "energy_flux_w_m2": "W_m2", "shelter_column_g_cm2": "g_cm2",
}
"""dict: Each field of `ROGUE_FIELDS` and its native unit (GEN.570), as `RoguePlanet.FIELD_UNITS` carries it."""

HP_ICE_DENSITY_KG_M3 = 1310.0
"""float: The density of the high-pressure ice (ice VI), kg/m3, that turns its thickness into the water it holds."""

NEUTRAL_CONTENT = radiogenic.Content(*radiogenic.BULK_SILICATE_EARTH, al26_ratio=0.0)
"""radiogenic.Content: The bulk silicate Earth's potassium, thorium and uranium for the crust's dose: a rogue planet has no star and
so no stored content, and a mantle with a lid that does not recycle (no felsic crust)."""


def _gas(planet):
    """`{gas-table formula: Pa}` of the regime's gas."""
    pressure_pa = planet.surface_pressure_pa or 0.0
    if planet.surface_regime in ("hydrogen-envelope", "hycean"):
        return {gases.stored("h2").formula: pressure_pa}
    if planet.surface_regime == "frozen-atmosphere":
        return {gases.stored("n2").formula: pressure_pa}
    return {}


def _hydrosphere(planet):
    """`(state, water mass fraction)` of a rogue planet's water: the state of `hydrosphere.HYDROSPHERE_STATES`, and the
    fraction its water column makes of the mass (0 when it has none)."""
    shell_km = planet.ice_shell_thickness_km or 0.0
    liquid_km = planet.ocean_depth_km or 0.0
    hp_km = planet.hp_ice_km or 0.0
    column_km = (liquid_km + shell_km * hydrosphere.ICE_DENSITY_KG_M3 / hydrosphere.LIQUID_WATER_DENSITY_KG_M3
                 + hp_km * HP_ICE_DENSITY_KG_M3 / hydrosphere.LIQUID_WATER_DENSITY_KG_M3)
    fraction = column_km / hydrosphere.water_layer_depth_km(planet.mass_kg, planet.radius_km, 1.0) if column_km else 0.0
    if planet.has_liquid_water:
        frozen = (planet.surface_temperature_k or 0.0) < hydrosphere.WATER_MELTING_K
        return ("ice-covered ocean" if frozen else "hycean"), fraction
    return ("ice" if shell_km or hp_km else "dry"), fraction


def values(planet):
    """A rogue planet as the column-named dict `habitability_world.compute` reads."""
    state, water_fraction = _hydrosphere(planet)
    gravity_g = rogue_surface.surface_gravity_m_s2(planet.mass_kg, planet.radius_km) / constants.EARTH_GRAVITY
    pressure_pa = planet.surface_pressure_pa or 0.0
    column = radiation.column_g_cm2(pressure_pa, gravity_g)
    basalt = surface_composition.RESERVOIRS["basaltic"]
    rock = (basalt["k_ppm"], basalt["th_ppm"], basalt["u_ppm"])
    ground = radiation.ground_dose_msv_yr(NEUTRAL_CONTENT, planet.age_gy or radiation.EARTH_AGE_GY, rock)
    out = {
        "planet_class": planet.planet_class, "body_type": planet.planet_type,
        "atmospheric_pressure_pa": pressure_pa, "surface_temperature_k": planet.surface_temperature_k,
        "gravity_g": gravity_g, "radius_km": planet.radius_km, "mass_kg": planet.mass_kg,
        "hydrosphere": state, "water_mass_fraction": water_fraction or None,
        "surface_dose_msv_yr": radiation.gcr_dose_msv_yr(column) + ground, "dose_ground_msv_yr": ground,
        "energy_flux_w_m2": habitability_world.energy_flux_w_m2({"hydrosphere": state}, None, None,
                                                                planet.internal_heat_flux_w_m2),
        "native_life": "none", "flux_variation": 0.0,
    }
    out[atmosphere.GAS_FIELD] = _gas(planet)
    if planet.has_liquid_water:
        kind = "ice-sealed" if planet.hp_ice_km else "neutral"
        (ph_lo, ph_hi), (aw_lo, aw_hi) = hydrosphere.OCEAN_CHEMISTRY[kind]
        out.update(ocean_class=kind, ocean_ph=(ph_lo + ph_hi) / 2.0, water_activity=(aw_lo + aw_hi) / 2.0,
                   phosphorus="starved" if kind == "ice-sealed" else "limited")
    return out


def generate(planet):
    """Sets `ROGUE_FIELDS` on a rogue planet from its surface conditions (call after they are set); a gas giant or brown
    dwarf gets `None` in each."""
    if planet.planet_type == "g" or planet.surface_temperature_k is None:
        for name in ROGUE_FIELDS:
            setattr(planet, name, None)
        return
    given = values(planet)
    for name, value in habitability_world.compute(given).items():
        if name in ROGUE_FIELDS:
            setattr(planet, name, value)
    planet.energy_flux_w_m2 = given["energy_flux_w_m2"]
