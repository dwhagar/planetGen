# planetgen/web/lib/planet_layout.py

"""
The parts of a planet's or moon's page that the redesign (UX.215 to UX.222) takes out of the description and the property
sheet and sets in a new order: the written description's blocks sorted into the section each belongs to, the table of the
air by mole fraction (UX.219), the table of the body's make-up by mass (UX.219) and the civilization paragraph (UX.221).
Imports no Flask.
"""

import re
from html import escape

from planetgen import tuning
from planetgen.describe.html import inline_html
from planetgen.generation import evolution
from planetgen.physics import gases, life_class_text, life_classes
from planetgen.util.format import three_figures_or_tiny

PHI6_CAPTIONS = ("Radiation", "Survival equipment")
"""tuple: The captions of the description's facts tables that belong with the PHI-6 rating."""

ICE_KEYS = (("ice_co2_mass_fraction", "Carbon dioxide ice, part of the ice"), ("ice_ch4_mass_fraction", "Methane ice, part of the ice"),
            ("ice_nh3_mass_fraction", "Ammonia ice, part of the ice"), ("ice_n2_mass_fraction", "Nitrogen ice, part of the ice"),
            ("ice_so2_mass_fraction", "Sulfur dioxide ice, part of the ice"), ("ice_co_mass_fraction", "Carbon monoxide ice, part of the ice"),
            ("ice_hcn_mass_fraction", "Hydrogen cyanide ice, part of the ice"), ("ice_other_mass_fraction", "Other ices, part of the ice"))
"""tuple: The non-water ices of a body by species and how the table words each."""

MASS_KEYS = (("metal_mass_fraction", "Metal (iron and nickel core)"), ("rock_mass_fraction", "Rock (silicate mantle and crust)"),
             ("ice_mass_fraction", "Ice"), *ICE_KEYS, ("envelope_mass_fraction", "Gas envelope"))
"""tuple: The stored mass fractions of a body and what each is, in the order the table lists them; the ice's species (GEN.834) are the parts of the ice."""

MASS_SHARED_KEYS = frozenset(key for key, _label in MASS_KEYS)
"""frozenset: The sheet rows the mass table shows, which the Properties list leaves out."""


def split_blocks(blocks):
    """
    The description's blocks of a planet or moon sorted by where the page puts them: `phi6` (the ratings table and the
    radiation and survival equipment tables), `life` (what the life looks like) and `rest` (anything else the page
    still shows, such as the sensors' flavour line and a gas giant's depth). The facts table of the body, the class
    sentence, the atmosphere and composition summaries, and the life class, system age and milestone sentences (`life-class`, `timeline`
    and `milestone`, which the Life section says from the stored columns and the life readings, UX.251) are dropped.
    """
    parts = {"phi6": [], "life": [], "rest": []}
    for block in blocks:
        kind = block["type"]
        if kind == "ratings" or (kind == "facts" and block["caption"] in PHI6_CAPTIONS):
            parts["phi6"].append(block)
        elif kind == "prose" and block["role"] == "life":
            parts["life"].append(block)
        elif kind == "prose" and block["role"] in ("flavour", "depth"):
            parts["rest"].append(block)
    return parts


def _rows(sheet, group_key):
    return next((group["rows"] for group in sheet["groups"] if group["key"] == group_key), [])


def _number(row):
    if row.get("raw") is not None:
        return row["raw"]
    item = row["value"][0] if row["value"] else None
    return item["qty"]["value"] if isinstance(item, dict) and "qty" in item else None


def gas_table(sheet):
    """
    The air of a body as `[{"name", "fraction", "pressure"}]`, the most abundant gas first: the share of the molecules each
    gas makes (its mole fraction, which is also its share of the volume and of the pressure) as text, and its partial
    pressure as the page's own HTML. Empty for a body with no stored gases.
    """
    held = [(row["key"][2:], row) for row in _rows(sheet, "atmosphere") if row["key"].startswith("p_") and _number(row)]
    total = sum(_number(row) for _gas, row in held)
    table = []
    for formula, row in held:
        share = _number(row) / total * 100.0
        table.append({"name": gases.gas(formula).name.capitalize(),
                      "fraction": "less than 0.001 percent" if share < 0.001 else f"{three_figures_or_tiny(share)} percent",
                      "pressure": inline_html(row["value"])})
    return table


def mass_table(sheet):
    """The make-up of a body by mass as `[{"name", "share"}]`: each stored mass fraction above zero as a percent of its mass
    (the four parts add up to the whole body; the ice's species add up to the ice). Empty when none is stored."""
    rows = {row["key"]: row for row in _rows(sheet, "composition")}
    table = []
    for key, label in MASS_KEYS:
        row = rows.get(key)
        fraction = None if row is None else _number(row)
        if fraction:
            table.append({"name": label, "share": f"{three_figures_or_tiny(fraction * 100.0)} percent of its mass"})
    return table


def civilization_sentences(row, species):
    """The paragraph on the people of a world, as HTML sentences: its dominant species and whether a technological
    civilization or a colony stands on it. Empty for a world with neither, so the page leaves the section out."""
    sentences = []
    if species:
        sentences.append(f'Dominant species: <a href="{escape(species["url"], quote=True)}">{escape(species["name"])}</a>')
    if row is not None and row.get("inhabited"):
        sentences.append("A technological civilization or a colony lives here.")
    return sentences


LIFE_SCORE_MAX = 9
"""int: The top of the class scores (a class that reaches the line scores 1 to 9)."""


def _class_name(class_id):
    """The full name of a class, such as "Oxygenic light-eaters: visible pigment (Chlorophyll a)"."""
    return life_classes.BY_ID[class_id].name if class_id in life_classes.BY_ID else f"Class {class_id}"


def _short_name(class_id):
    """The name of a class before its variant, such as "oxygenic light-eaters"."""
    name = re.sub(r"\s*\([^)]*\)$", "", _class_name(class_id).partition(":")[0])
    return name[:1].lower() + name[1:]


GLANCE_OTHERS = 2
"""int: How many classes besides the best-suited one the At a glance paragraph names."""


def _label(class_id, class_url=None):
    """"Class 1a", a link to the class's help page when `class_url` (a function from a class id to an address) is given."""
    if class_url is None:
        return f"Class {class_id}"
    return f'<a class="help-link" href="{escape(class_url(class_id))}">Class {escape(class_id)}</a>'


def _class_phrase(class_id, class_url=None):
    return f"{_label(class_id, class_url)} ({escape(_short_name(class_id))})"


def _glance_classes(best, others, class_url=None):
    """The At a glance sentence naming the classes the conditions suit: "Conditions suit Class 1a (oxygenic light-eaters) best, then
    Class 2 (...) and Class 3a (...)." No scores; at most `GLANCE_OTHERS` others, then "and others"."""
    text = f"Conditions suit {_class_phrase(best, class_url)} best"
    named = [_class_phrase(class_id, class_url) for class_id in others[:GLANCE_OTHERS]]
    if len(others) > GLANCE_OTHERS:
        named.append("others")
    if named:
        text += ", then " + (", ".join(named[:-1]) + " and " if len(named) > 1 else "") + named[-1]
    return text + "."


def _glance(readings, living, best, best_score, others, class_url=None):
    """
    The Life paragraph of At a glance (UX.251): the sentence of each kind of life that is there ("Microbes: abundant, stable, secure."),
    or "Nothing lives here yet." when the conditions suit a class and none has appeared, then the classes the conditions suit. `None`
    for a world with no life that could not have any. Words only; the numbers are in the Life section.
    """
    parts = [reading["text"] for reading in readings if reading.get("abundance") not in (None, "Absent")]
    if not parts and best and best_score:
        parts.append("Nothing lives here yet." if not living else "Life here is too sparse to count.")
    if best and best_score:
        parts.append(_glance_classes(best, others, class_url))
    return " ".join(parts) or None


def _years_text(gyr):
    return None if gyr is None else f"about {gyr:.1f} billion years"


def _reading_rows(readings):
    """The per-kind rows of the Life section's table: `{"kind", "abundance", "trend", "security", "left"}`, each cell the plain word
    (empty for what a kind that is not there has no word for) and `left` the time the living class's conditions last, empty unless it is known."""
    rows = []
    for reading in readings:
        score = reading.get("score")
        shown = reading.get("abundance") or "Unknown"
        rows.append({"kind": reading["label"], "abundance": f"{shown} ({int(score)} of {LIFE_SCORE_MAX})" if score else shown,
                     "trend": reading.get("trend") or "", "security": reading.get("security") or "",
                     "left": _years_text(reading.get("years_left_gyr")) or ""})
    return rows


def _life_facts(row, living):
    """The stored facts of the world's life as `[(label, text)]` for the Life section: the stage, the animal grade, the pace and the ages."""
    facts = []
    stage = row.get("life_stage")
    grade = int(row.get("life_animal_grade") or 0) if stage == "multicellularity" else 0
    if living and stage:
        facts.append(("Stage", evolution.stage_name(stage, living, grade).lower()))
    if grade:
        facts.append(("Animal grade", life_class_text.GRADES[min(grade, 3) - 1]))
    pace = row.get("evolutionary_pace")
    if pace in tuning.EVOLUTIONARY_PACE_LABELS and (living or row.get("life_class_best")):
        facts.append(("Pace of evolution", tuning.EVOLUTIONARY_PACE_LABELS[pace]))
    if row.get("life_system_age_gyr") is not None and (living or row.get("life_class_best")):
        facts.append(("Age of the system", f"{row['life_system_age_gyr']:.2f} billion years"))
    if living and stage and row.get("life_milestone_age_gyr") is not None:
        facts.append(("Stage reached at", f"{row['life_milestone_age_gyr']:.2f} billion years"))
    return facts


def life_class_view(row, readings=(), class_url=None):
    """
    The life lines of a planet's or moon's page (UX.237, UX.251) from its system-detail row (API.37) and, for a world with life, the
    `readings` of `GET /api/planets/<id>/life-readings`: `glance` (the paragraph of At a glance), `paragraph` (the lead of the Life
    section: what lives there and the classes the conditions suit with their suitability), `readings` (the table of the three words for each
    kind of life), `facts` (the stage, grade, pace and ages as `[(label, text)]`), `supports`, `also`, `has` and `words` (what that class,
    or else the best-supported one, drinks, eats, breathes and gives off, as `{"class", "name", "rows": [(label, text)]}`). A line with
    nothing to say is `None` or empty; a body that supports no class and has no life gives `None` for the whole view. `glance` and
    `paragraph` are HTML: with `class_url` (a function from a class id to the address of its help page) each class name links to it.
    """
    if not row:
        return None
    scores = {class_id: int(score) for class_id, score in (row.get("life_class_scores") or {}).items() if score}
    best, best_score = row.get("life_class_best"), row.get("life_class_best_score")
    living = row.get("life_class_id")
    view = {"alongside": None, "supports": None, "also": [], "has": None, "words": None, "paragraph": None, "glance": None, "readings": [], "facts": []}
    if best and best_score:
        view["supports"] = f"Best suited to Class {best} life ({_short_name(best)}), suitability {best_score} of {LIFE_SCORE_MAX}."
    others = sorted(((class_id, score) for class_id, score in scores.items() if class_id != best), key=lambda pair: (-pair[1], pair[0]))
    view["also"] = [f"Class {class_id} ({_short_name(class_id)}, suitability {score} of {LIFE_SCORE_MAX})" for class_id, score in others]
    stage = row.get("life_stage")
    if living and stage:
        grade = int(row.get("life_animal_grade") or 0) if stage == "multicellularity" else 0
        text = f"Life here is Class {living}, at the stage of {evolution.stage_name(stage, living, grade).lower()}"
        view["has"] = text + (f" (animal grade {life_class_text.GRADES[min(grade, 3) - 1]})." if grade else ".")
    shown = living or best
    if shown and shown in life_classes.BY_ID:
        said = life_class_text.lines(shown)
        view["words"] = {"class": shown, "name": _class_name(shown), "lives": bool(living),
                         "rows": [(label, said[key]) for label, key in (("Drinks", "drinks"), ("Eats", "eats"),
                                                                       ("Breathes", "breathes"), ("Gives off", "gives_off"))]}
    view["glance"] = _glance(readings, living, best, best_score, [class_id for class_id, _score in others], class_url)
    view["readings"] = _reading_rows(readings)
    view["facts"] = _life_facts(row, living)
    view["alongside"] = _alongside(row, class_url)
    view["paragraph"] = _life_paragraph(view, best, best_score, living, row, others, class_url)
    return view if any((view["supports"], view["also"], view["has"], view["words"], view["glance"], view["alongside"])) else None


def _alongside(row, class_url=None):
    """The other classes that live on the body beside the dominant one (GEN.1209), as HTML: "Also living here: Class 3a (31 percent of the food production, at
    the stage of single-celled organisms), ...". `None` when only one class lives there."""
    others = [entry for entry in (row.get("life_classes") or [])[1:] if entry.get("life_stage")]
    if not others:
        return None
    parts = []
    for entry in others:
        class_id = entry["life_class_id"]
        grade = int(entry.get("animal_grade") or 0) if entry["life_stage"] == "multicellularity" else 0
        stage = escape(evolution.stage_name(entry["life_stage"], class_id, grade).lower())
        parts.append(f"{_label(class_id, class_url)} ({round(100.0 * entry['share'])} percent of the food production, at the stage of {stage})")
    return "Also living here: " + ", ".join(parts) + "."


def _life_paragraph(view, best, best_score, living, row, others, class_url=None):
    """
    The lead of the Life section (UX.237), as HTML: what lives there ("Nothing lives here yet." when the conditions suit a class and none
    has appeared), then the classes the conditions suit with their suitability. `None` when there is nothing to say. How common each kind
    of life is comes in the table under it (`_reading_rows`).
    """
    parts = []
    if view["has"]:
        stage = row.get("life_stage")
        grade = int(row.get("life_animal_grade") or 0) if stage == "multicellularity" else 0
        text = f"Life here is {_label(living, class_url)}, at the stage of {escape(evolution.stage_name(stage, living, grade).lower())}"
        parts.append(text + (f" (animal grade {life_class_text.GRADES[min(grade, 3) - 1]})." if grade else "."))
    elif best and best_score:
        parts.append("Nothing lives here yet.")
    if view.get("alongside"):
        parts.append(view["alongside"])
    if best and best_score:
        suited = f"Conditions here suit {_label(best, class_url)} life ({escape(_short_name(best))}) best, suitability {best_score} of {LIFE_SCORE_MAX}"
        if others:
            suited += "; they also suit " + ", ".join(
                f"{_label(class_id, class_url)} ({escape(_short_name(class_id))}, suitability {score} of {LIFE_SCORE_MAX})" for class_id, score in others)
        parts.append(suited + ".")
    return " ".join(parts) or None
