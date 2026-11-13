# planetgen/web/lib/planet_summary.py

"""
The short written summary of a planet or moon page (UX.220): a paragraph made from the stored data alone, so no two
sentences are canned text. It says what class the body is and the zone it is in, what its surface is like, how the PHI-6
rating and the Suitability score read and what safety equipment a person needs. The life of the world is its own paragraph
(`planet_layout.life_class_view`, UX.251), not part of this one.

The quantities arrive as HTML (`values`, a sheet key to the page's own rendering, so they read in the reader's unit
system); everything else is written here from the row of the system detail. Imports no Flask.
"""

from html import escape

from planetgen.describe import sheet as sheet_registry
from planetgen.physics import class_code, habitability, planet_type, roundness
from planetgen.physics.habitability_world import TIER_NAMES

def _class_sentence(name, row, values):
    code = values.get("planet_class") or escape(str(row.get("planet_class") or ""))
    if code and row.get("class_variant"):
        code = escape(class_code.format_code(row["planet_class"], row["class_variant"]))     # GEN.1220: the digit, when it is not 0
    giant = row.get("body_type") == "g"
    ratio = row.get("roundness_ratio")
    asteroid = not code and not giant and ratio is not None and not roundness.is_round(ratio)          # GEN.1160: under a roundness ratio of 1 a body has no class, it is an Asteroid
    kind = "asteroid" if asteroid else ("Class " + code + " " if code else "") + ("gas giant" if giant else "world")
    zone = sheet_registry.ZONE_NAMES.get(str(row.get("zone") or "").lower())
    where = f" in the {zone}" if zone else ""
    article = "an" if kind[:1] in "aeiouAEIOU" else "a"
    return f"{escape(name)} is {article} {kind}{where}.{_match_sentence(row, code)}"


def _match_sentence(row, code):
    """GEN.1207: how well the class fits, for example `It fits class M 92 percent and class O 41 percent.` (a body that fits no class says its nearest is shown)."""
    match = row.get("class_match")
    if match is None or not code:
        return ""
    if match <= 0.0:
        return f" No class fits it; class {code} is the nearest."
    text = f" It fits class {code} {match * 100:.0f} percent"
    runner = row.get("runner_up_class")
    if runner and row.get("runner_up_match"):
        text += f" and class {escape(str(runner))} {row['runner_up_match'] * 100:.0f} percent"
    return text + "."


def _surface_sentence(row, values):
    if row.get("body_type") == "g":
        return None
    parts = []
    ground = planet_type.has_ground_type(row.get("planet_type"))      # UX.297: no ground, no surface: the numbers are at the 1 bar reference level
    if values.get("surface_temperature_k"):
        parts.append(f"its surface is {values['surface_temperature_k']}" if ground else f"its temperature at the 1 bar reference level is {values['surface_temperature_k']}")
    pressure = row.get("atmospheric_pressure_pa")
    airless = pressure is not None and pressure < sheet_registry.VACUUM_LIMIT_PA
    if airless:
        parts.append("it has no atmosphere")
    elif values.get("atmospheric_pressure_pa"):
        parts.append(f"the air at the surface is {values['atmospheric_pressure_pa']}" if ground else f"the air at its reference level is {values['atmospheric_pressure_pa']}")
    if values.get("gravity_g"):
        parts.append(f"gravity is {values['gravity_g']}")
    if not parts:
        return None
    sentence = ", ".join(parts[:-1]) + (", and " if len(parts) > 2 else " and ") + parts[-1] if len(parts) > 1 else parts[0]
    return sentence[0].upper() + sentence[1:] + "."


def _day_sentence(values):
    """UX.205: how much of the time the star is up, from the stored shares. A world that turns once an orbit has a mean share and no
    solar day: it says so (day side, terminator, night side). A moon in its planet's shadow adds that night and the longest eclipse."""
    if not values.get("day_share_mean"):
        return None
    parts = []
    if values.get("day_share_equator"):
        parts.append(f"its star is up {values['day_share_equator']} of the time at the equator")
        if values.get("daylight_hours") and values.get("solar_day_hours"):
            parts.append(f"{values['daylight_hours']} of daylight in each solar day of {values['solar_day_hours']}")
        parts.append(f"{values['day_share_mean']} averaged over the whole surface")
        sentence = parts[0] + ", " + (parts[1] + ", and " if len(parts) > 2 else "") + parts[-1]
        return sentence[0].upper() + sentence[1:] + "."
    return ("It keeps one face to its star: the day side is always in daylight, the terminator ring sits on the line between day and "
            "night, and the night side is always dark.")


def _shadow_sentence(values, raw):
    """UX.205: a moon's extra night in its planet's shadow, with the longest eclipse."""
    if not (raw.get("parent_shadow_share") or 0) > 0 or not values.get("parent_shadow_share"):
        return None
    sentence = f"Its planet's shadow adds a night: it is in the shadow {values['parent_shadow_share']} of the time"
    if values.get("longest_eclipse_hours"):
        sentence += f", the longest eclipse lasting {values['longest_eclipse_hours']}"
    return sentence + "."


def _rating_sentence(row):
    if row.get("phi4") is None:
        return None
    overall = habitability.tier(row["phi4"])
    sentence = f"Its PHI-6 habitability rating is {overall}"
    sentence += ", so people can live here" if habitability.is_habitable(row["phi4"]) else ""
    suitability = row.get("tier_suitability")
    if suitability is not None:
        sentence += f", and its Suitability score is {TIER_NAMES[suitability]}"
    return sentence + "."


def _equipment_sentence(row):
    gear = row.get("equipment")
    if not gear:
        return None
    names = [escape(str(name)) for name in gear.get("required") or ()]
    if names:
        listed = ", ".join(names[:-1]) + (", and " if len(names) > 2 else " and ") + names[-1] if len(names) > 1 else names[0]
        return f"A person on the surface needs {listed}."
    return "A person on the surface needs no safety equipment."


def summary_sentences(name, row, values, raw=None):
    """The sentences of the page's summary paragraph, as HTML, in the order they are read. `raw` is the stored numbers of the sheet's
    rows (a sheet key to its number)."""
    sentences = [_class_sentence(name, row, values), _surface_sentence(row, values), _day_sentence(values),
                 _shadow_sentence(values, raw or {}), _rating_sentence(row),
                 _equipment_sentence(row)]
    return [sentence for sentence in sentences if sentence]
