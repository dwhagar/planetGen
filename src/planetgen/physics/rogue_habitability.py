# planetgen/physics/rogue_habitability.py

"""
The habitability fields of a free-floating brown dwarf.

A brown dwarf has no surface and scores nothing: every field is `None`. (The rogue planets that systems throw out are built by the planet chain
with no star and scored like any planet, `generation/phenomena/starless_host.py`, GEN.1111; the surface regimes this module used to read for the
older flat-rate rogue are gone.)
"""

from planetgen.db import column_units
from planetgen.physics import habitability_world, outside_hours, regional_temperature

ROGUE_FIELDS = (*(name for name in habitability_world.SCORE_FIELDS
                  if name not in (*regional_temperature.FIELDS, "outside_hours", *outside_hours.FIELDS)), "energy_flux_w_m2")
"""tuple: What `generate` sets on a rogue planet, as stored in `rogue_planets`: the score of a star system's body without the
regional temperatures and the hours outside (a rogue planet has no star, so no latitude bands, day side or daylight) and its stored energy flux (there is
no star, so no flux variation)."""

ROGUE_UNITS = {
    **{name: "index" for name in ROGUE_FIELDS if name.startswith(("tier_", "life_", "eq_"))},
    **{name: "ratio" for name in ROGUE_FIELDS if name.startswith("phi4") or name in ("l_solv", "l_chem", "l_ener", "l_rad")},
    **{name: "K" for name in column_units.COVERAGE_KELVIN}, **{name: "frac" for name in column_units.COVERAGE_SHARES},
    **{name: "deg" for name in column_units.COVERAGE_DEGREES}, **{name: "d" for name in column_units.COVERAGE_DAYS},
    **{name: "flag" for name in column_units.COVERAGE_FLAGS},       # GEN.1265: the temperature coverage columns a rogue row has (all None for a brown dwarf)
    "energy_flux_w_m2": "W_m2", "shelter_column_g_cm2": "g_cm2",
}
"""dict: Each field of `ROGUE_FIELDS` and its native unit (GEN.570), as `RoguePlanet.FIELD_UNITS` carries it."""

def generate(planet):
    """Sets `ROGUE_FIELDS` on a brown dwarf: `None` in each."""
    for name in ROGUE_FIELDS:
        setattr(planet, name, None)
