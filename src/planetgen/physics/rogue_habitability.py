# planetgen/physics/rogue_habitability.py

"""
The habitability fields of a free-floating brown dwarf.

A brown dwarf has no surface and scores nothing: every field is `None`. (The rogue planets that systems throw out are built by the planet chain
with no star and scored like any planet, `generation/phenomena/starless_host.py`, GEN.1111; the surface regimes this module used to read for the
older flat-rate rogue are gone.)
"""

from planetgen.physics import habitability_world, outside_hours, regional_temperature, temperature_coverage

ROGUE_FIELDS = (*(name for name in habitability_world.SCORE_FIELDS
                  if name not in (*regional_temperature.FIELDS, "outside_hours", *outside_hours.FIELDS)), "energy_flux_w_m2")
"""tuple: What `generate` sets on a rogue planet, as stored in `rogue_planets`: the score of a star system's body without the
regional temperatures and the hours outside (a rogue planet has no star, so no latitude bands, day side or daylight) and its stored energy flux (there is
no star, so no flux variation)."""

ROGUE_UNITS = {
    **{name: "index" for name in ROGUE_FIELDS if name.startswith(("tier_", "life_", "eq_"))},
    **{name: "ratio" for name in ROGUE_FIELDS if name.startswith("phi4") or name in ("l_solv", "l_chem", "l_ener", "l_rad")},
    "energy_flux_w_m2": "W_m2", "shelter_column_g_cm2": "g_cm2",
    **{name: ("K" if name.endswith("_k") else "deg" if name.endswith("_deg") else "d" if name.endswith("_days") or name.endswith("_day") else "flag" if name == "temp_livable_belt_only" else "frac")
       for name in temperature_coverage.FIELDS},
}
"""dict: Each field of `ROGUE_FIELDS` and its native unit (GEN.570), as `RoguePlanet.FIELD_UNITS` carries it."""

def generate(planet):
    """Sets `ROGUE_FIELDS` on a brown dwarf: `None` in each."""
    for name in ROGUE_FIELDS:
        setattr(planet, name, None)
