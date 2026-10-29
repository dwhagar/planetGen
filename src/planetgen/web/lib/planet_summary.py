# planetgen/web/lib/planet_summary.py

"""
The short written summary of a planet or moon page (UX.220): a paragraph made from the stored data alone, so no two
sentences are canned text. It says what class the body is and the zone it is in, what its surface is like, how the PHI-6
rating and the Suitability score read, what safety equipment a person needs, and how much life it holds.

The quantities arrive as HTML (`values`, a sheet key to the page's own rendering, so they read in the reader's unit
system); everything else is written here from the row of the system detail. Imports no Flask.
"""

from html import escape

from planetgen.describe import sheet as sheet_registry
from planetgen.physics import habitability
from planetgen.physics.habitability_world import TIER_NAMES
from planetgen.web.lib.habitability_symbols import LIFE_MAX, LIFE_NAMES, LIFE_LABELS

PROSE_ROLES = {"life": "life", "timeline": "life"}
"""dict: The description's prose roles the Life paragraph of a page is made from."""


def _class_sentence(name, row, values):
    code = values.get("planet_class") or escape(str(row.get("planet_class") or ""))
    giant = row.get("body_type") == "g"
    kind = ("Class " + code + " " if code else "") + ("gas giant" if giant else "world")
    zone = sheet_registry.ZONE_NAMES.get(str(row.get("zone") or "").lower())
    where = f" in the {zone}" if zone else ""
    article = "an" if kind[:1] in "aeiouAEIOU" else "a"
    return f"{escape(name)} is {article} {kind}{where}."


def _surface_sentence(row, values):
    if row.get("body_type") == "g":
        return None
    parts = []
    if values.get("surface_temperature_k"):
        parts.append(f"its surface is {values['surface_temperature_k']}")
    pressure = row.get("atmospheric_pressure_pa")
    airless = pressure is not None and pressure < sheet_registry.VACUUM_LIMIT_PA
    if airless:
        parts.append("it has no atmosphere")
    elif values.get("atmospheric_pressure_pa"):
        parts.append(f"the air at the surface is {values['atmospheric_pressure_pa']}")
    if values.get("gravity_g"):
        parts.append(f"gravity is {values['gravity_g']}")
    if not parts:
        return None
    sentence = ", ".join(parts[:-1]) + (", and " if len(parts) > 2 else " and ") + parts[-1] if len(parts) > 1 else parts[0]
    return sentence[0].upper() + sentence[1:] + "."


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


def _life_sentence(row):
    scores = [row.get(f"life_{name}") for name in LIFE_NAMES]
    if any(score is None for score in scores):
        return None
    if not any(scores):
        return "No life has been found here."
    shown = ", ".join(f"{LIFE_LABELS[name].lower()} {score} of {LIFE_MAX}" for name, score in zip(LIFE_NAMES, scores) if score)
    return f"Its life scores are {shown}."


def summary_sentences(name, row, values):
    """The sentences of the page's summary paragraph, as HTML, in the order they are read."""
    sentences = [_class_sentence(name, row, values), _surface_sentence(row, values), _rating_sentence(row),
                 _equipment_sentence(row), _life_sentence(row)]
    return [sentence for sentence in sentences if sentence]
