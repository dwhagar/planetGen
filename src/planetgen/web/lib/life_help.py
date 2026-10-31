# planetgen/web/lib/life_help.py

"""
The help pages about life (`/help/life`, DOC.84): the four kinds of life, one page for each scale a planet page uses (abundance,
trend, security, suitability), the stages of life and the animal grade, and one page for each of the nine life classes.

Every number on a page is read from the code that uses it (`physics.life_readings`, `physics.life_classes`, `physics.life_scores`,
`tuning`, the rule file `physics/data/life_classes.json`), so a page cannot disagree with the score it explains. Imports no Flask.

`page(topic, number, href)` returns the page as `{"title", "description", "intro": [text], "sections": [{"id", "heading", "paragraphs",
"table", "facts"}]}`, or `None` for a topic or class that does not exist. A table is `{"head": [text], "rows": [[cell]]}`; a cell is
text or `(text, url)`. `facts` is a list of `(label, value)`, a value being text or a list of phrases. `href(topic, number=None)` gives
the address of a page.
"""

from planetgen import tuning
from planetgen.physics import life_air, life_class_text, life_classes, life_readings
from planetgen.physics.life_readings import ABUNDANCE_BANDS
from planetgen.web.lib import help_registry
from planetgen.web.lib.help_registry import Child, Term

TOPICS = (
    ("kinds", "Kinds of life", "Microbes, flora, fauna and sentient life: what each word means and what its score is judged from."),
    ("abundance", "How much", "Absent, rare, occasional, frequent, abundant, dominant: how much of a kind of life a world holds."),
    ("trend", "Trend", "Growing, stable or declining: which way a kind of life is going."),
    ("security", "Security", "Secure, marginal, precarious or dying out: how long the conditions that keep the life alive will last."),
    ("suitability", "Suitability", "How well a world's conditions suit each life class, from 0 to 9."),
    ("stages", "Stages and animal grades", "The steps life takes on a world, from the first cells to a civilization."),
    ("classes", "Life classes", "The nine life classes: what each drinks, eats, breathes and gives off, and where it can live."),
)
"""tuple: `(slug, title, summary)` of each life help page, in the order the pages list each other."""

TIER_NAMES = {
    "T1": "Verified", "T2": "Confident unverified", "T3": "Marginal unverified", "T4": "Unverified",
    "T5": "Confident educated guess based on physics", "T6": "Marginal educated guess based on physics",
    "T7": "Educated guess based on inferred physics", "T8": "Educated guess based on reputable fringe or new physics",
    "T9": "Educated guess based on needed output",
}
"""dict: What each source tier of `docs/design/research/stand-in-values.md` is called."""

LEVEL_NAMES = {
    "L0": "As the source gives it", "L1": "Nudged inside the source's uncertainty", "L2": "Tuned to a published statistic",
    "L3": "Hand-tuned against many generated worlds", "L4": "Fitted to the outcome with no independent source",
}
"""dict: What each adjustment level is called."""

SOLVENT_WORDS = {
    "water": "water", "hydrocarbon": "liquid methane and ethane", "co2": "liquid or supercritical carbon dioxide",
    "sulfuric": "concentrated sulfuric acid", "ammonia_water": "water with ammonia in it",
}
"""dict: The liquid a class needs, by the rule file's `solvent.medium`."""

HYDROSPHERE_WORDS = {"surface ocean": "an open surface ocean", "hycean": "a hycean sea"}

SUPPLY_WORDS = {
    "h2_supply": "a supply of hydrogen from the rock", "radiolytic_h2_supply": "a supply of hydrogen freed by radiation",
    "donor_supply": "a supply of fuel to eat (hydrogen, iron or sulfur compounds)", "nitrate_supply": "a supply of nitrate to breathe",
    "perchlorate_supply": "a supply of perchlorate to breathe", "sulfur_supply": "a supply of reduced sulfur",
    "hydrogen_acetylene_supply": "a supply of hydrogen and acetylene from the haze",
}
"""dict: The fuel or oxidant a class needs delivered, by the rule file's `supply.quantity`."""


def _number(value):
    return f"{value:,.0f}" if float(value).is_integer() else f"{value:g}"


def _pressure(pa):
    return f"{_number(pa / 1000.0)} kPa" if pa >= 1000.0 else f"{_number(pa)} Pa"


def _billion_years(value):
    return f"{value:g} billion years"


def _table(head, rows):
    return {"head": head, "rows": rows}


def _section(key, heading, paragraphs=(), table=None, facts=None):
    return {"id": key, "heading": heading, "paragraphs": list(paragraphs), "table": table, "facts": facts}


# --- the scales -----------------------------------------------------------------------------------------------

def _kinds(href):
    score_max = tuning.LIFE_SCORE_MAX
    full_years = int(tuning.LIFE_SCORE_SENTIENT_FULL_GYR * 1000)
    intro = [
        f"A planet or moon page rates four kinds of life. Each rating is a whole number from 0 (none) to {score_max} (abundant).",
        "The page then gives three words for each kind: how much of it there is, which way it is going and how secure it is.",
    ]
    rows = [
        ["Microbes", "Single cells, mats and everything smaller than a plant or an animal.", "Any stage of life, from the first cells."],
        ["Flora", "Plant-like life: producers that make their own food, from algae and mats to land plants.",
         f"A class that eats light and has reached its producers. Mats and algae never score above {tuning.LIFE_SCORE_SIMPLE_FLORA_MAX}; "
         "a higher score needs life made of many cells."],
        ["Fauna", "Animals, of any size, from sponge-like feeders to large active hunters.",
         "Life made of many cells and an animal grade. Fauna never scores above the producers it feeds on."],
        ["Sentient life", "Species that build tools, settlements and a civilization.",
         f"The stage of a technological civilization. The score grows over {full_years:,} million years from that stage until it matches fauna."],
    ]
    sections = [
        _section("kinds", "The four kinds", table=_table(["Kind", "What it means", "What it needs"], rows)),
        _section("score", "How the number is made", [
            f"The number is {score_max} times how well the world's conditions suit that kind of life, rounded up. Below {tuning.HABITABILITY_LIFE_SCORE_MIN:.1%} "
            "it counts as none.",
            "The life class that lives on the world sets how well the conditions suit. The stage the life has reached sets which kinds can score at all.",
        ]),
        _section("next", "The words beside the number", [
            "How much, trend and security each have a page of their own. Suitability is a different score: it rates the life classes, not the kinds of life.",
        ]),
    ]
    return "Kinds of life", "What microbes, flora, fauna and sentient life mean on a planet page.", intro, sections


def _abundance(href):
    words = [word for _low, word in ABUNDANCE_BANDS]
    bands = []
    for index, (low, word) in enumerate(ABUNDANCE_BANDS):
        high = ABUNDANCE_BANDS[index + 1][0] - 1 if index + 1 < len(ABUNDANCE_BANDS) else tuning.LIFE_SCORE_MAX
        bands.append([word, str(low) if low == high else f"{low} to {high}"])
    meanings = {
        "Absent": "None of this kind of life is on the world.",
        "Rare": "A thin hold: a few places, or a trace.",
        "Occasional": "Present here and there, with gaps.",
        "Frequent": "Found in many places.",
        "Abundant": "Widespread and thriving.",
        "Dominant": "As much as the world's conditions allow.",
    }
    rows = [[word, score, meanings[word]] for word, score in bands]
    intro = [
        "How much is the first of the three words for each kind of life. It comes straight from the kind's score from 0 to "
        f"{tuning.LIFE_SCORE_MAX}.",
        "The words are the ecologists' abundance scale: dominant, abundant, frequent, occasional, rare.",
    ]
    sections = [
        _section("words", "The words", table=_table(["Word", "Score", "What it means"], rows)),
        _section("judged", "What it is judged from", [
            "The score is how well the world's conditions suit that kind of life, and which stage the life has reached. "
            "The page on kinds of life says how each score is made.",
            f"A world with no life class has no abundance to judge, so every kind reads {words[0].lower()}.",
        ]),
    ]
    return "How much", "What absent, rare, occasional, frequent, abundant and dominant mean.", intro, sections


def _trend(href):
    intro = ["Trend is the second of the three words for each kind of life. It says which way the life is going, as the star ages and brightens."]
    rows = [
        ["Growing", f"The conditions are getting better by more than {life_readings.TREND_BAND:.0%}, or the kind of life only just arrived."],
        ["Stable", "The conditions are about the same, and the kind of life has been there a while."],
        ["Declining", f"The conditions are getting worse by more than {life_readings.TREND_BAND:.0%}."],
    ]
    sections = [
        _section("words", "The words", table=_table(["Word", "What it means"], rows)),
        _section("judged", "What it is judged from", [
            f"The world's life class is scored again with its star {_billion_years(life_readings.TREND_SPAN_GYR)} older. The star is brighter, "
            "so the surface is warmer, and an ocean that boils away leaves the water-eating classes with nothing.",
            f"A change of more than {life_readings.TREND_BAND:.0%} of the score now is a trend. A kind whose stage came within "
            f"{_billion_years(life_readings.SPREAD_GYR)} of the system's age is still spreading, so it reads growing unless the conditions are declining.",
            "A world with no star to look at is judged by its timeline alone.",
        ]),
    ]
    return "Trend", "What growing, stable and declining mean.", intro, sections


def _security(href):
    intro = [
        "Security is the third of the three words for each kind of life. It says how long the conditions that keep the life alive will last "
        "as the star brightens.",
        "All the kinds of life on one world share its living class, so they share the same security.",
    ]
    rows = [
        ["Secure", f"The conditions last more than {_billion_years(life_readings.MARGINAL_GYR)}, and the class scores above {life_readings.MARGINAL_SCORE} now."],
        ["Marginal", f"The conditions last up to {_billion_years(life_readings.MARGINAL_GYR)}, or the class scores {life_readings.MARGINAL_SCORE} or less now."],
        ["Precarious", f"The conditions last {_billion_years(life_readings.PRECARIOUS_GYR)} or less."],
        ["Dying out", f"The class scores nothing now, or the conditions last {_billion_years(life_readings.DYING_GYR)} or less."],
    ]
    sections = [
        _section("words", "The words", table=_table(["Word", "What it means"], rows)),
        _section("judged", "What it is judged from", [
            f"The world's living class is scored again every {_billion_years(life_readings.STEP_GYR)} as the star brightens, up to "
            f"{_billion_years(life_readings.MARGINAL_GYR)} ahead. The first step where it scores nothing, or where the star leaves the main sequence, is the time left.",
            "A page shows the time left in the Life section when it is known and shorter than the span searched.",
        ]),
    ]
    return "Security", "What secure, marginal, precarious and dying out mean.", intro, sections


def _suitability(href):
    score_max = tuning.LIFE_SCORE_MAX
    intro = [
        f"Suitability rates how well a world's conditions suit a life class, from 0 to {score_max}. "
        "A planet page lists the classes that suit the world, best first.",
        "It is not the Suitability score of the habitability colors, which rates the world for a person.",
    ]
    sections = [
        _section("judged", "What it is judged from", [
            "The score looks at four things together: the liquid the class needs, the chemistry its life needs, the energy it can eat and the radiation it can stand. "
            "Each class also has a band of temperature and pressure, and a supply of fuel it must have.",
            f"The result is a number from 0 to 1. The score is {score_max} times that number, rounded up. Below {tuning.HABITABILITY_LIFE_SCORE_MIN:.1%} it is 0.",
        ]),
        _section("best", "Best suited and living there", [
            "The best-suited class is the one with the highest score. The class that lives on a world is the supported class that makes the most food, "
            "so it is not always the best-suited one.",
            "A world can suit a class and have no life of it yet: the star may be too young for it to have started.",
        ]),
        _section("classes", "The classes", ["Each of the nine life classes has its own page, listed on the page of life classes."]),
    ]
    return "Suitability", "How the 0 to 9 suitability of a life class is judged.", intro, sections


def _stages(href):
    rows = [[str(rung), life_class_text.DOCUMENT["rungs"][str(rung)]["name"],
             life_class_text.DOCUMENT["rungs"][str(rung)]["description"].replace("Class {n} life", "life of the class")] for rung in life_class_text.RUNGS]
    grades = [[key, entry["name"], entry["description"]] for key, entry in life_class_text.DOCUMENT["grades"].items()]
    pace = [[word, f"A world on the {key} schedule."] for key, word in tuning.EVOLUTIONARY_PACE_LABELS.items()]
    intro = [
        "Life on a world climbs a ladder of stages. A world's page names the highest stage its life has reached, using the words of its life class.",
        "The stages come from a schedule of five hard steps, drawn for each world. The stage reached is the latest step that came before the system's present age.",
    ]
    sections = [
        _section("stages", "The stages", table=_table(["Step", "Name", "What it means"], rows),
                 paragraphs=["Each life class has its own names for some steps. These are the names used when a class has none of its own."]),
        _section("grades", "Animal grades", table=_table(["Grade", "Name", "What it means"], grades),
                 paragraphs=["Animals are graded by how much energy and oxidant their world gives them. A world with little oxygen can reach grade A1 or A2, and one with plenty can reach A3."]),
        _section("pace", "Pace of evolution", table=_table(["Word", "Meaning"], pace),
                 paragraphs=["Each world draws a pace. It sets how soon each step arrives."]),
    ]
    return "Stages and animal grades", "The steps life takes on a world, and the three animal grades.", intro, sections


# --- the classes ----------------------------------------------------------------------------------------------

def family(number):
    """The rows of the life class `number` (its variants), in order; empty for a number that is not a class."""
    return [row for row in life_classes.CLASSES if row.number == number]


def family_name(number):
    """The name of a class family: the name before the colon, such as "Hydrogen eaters"."""
    rows = family(number)
    return rows[0].name.split(":")[0] if rows else None


def class_numbers():
    """The life class numbers with a page of their own, read from the rule file: `[(number, label, summary)]`, in order."""
    numbers = sorted({row.number for row in life_classes.CLASSES})
    return [(number, f"{family_label(number)}: {family_name(number)}",
             f"What the {family_name(number).lower()} drink, eat, breathe and give off, and where they can live.") for number in numbers]


def family_label(number):
    """How a table names a class: "Class 2" for a class with no letters, "Class 3 family" for one with variants."""
    rows = family(number)
    return f"Class {rows[0].id}" if len(rows) == 1 else f"Class {number} family"


def _variant_name(row):
    return row.name.split(":", 1)[1].strip() if ":" in row.name else row.name


def _gate_words(gate):
    quantity = gate["quantity"]
    if quantity == "hydrosphere":
        return "The surface has " + " or ".join(HYDROSPHERE_WORDS.get(v, v) for v in gate["in"])
    if quantity == "pressure_pa" and gate.get("min", 0.0) > 0.0:
        return f"Air pressure of at least {_pressure(gate['min'])}"
    if quantity == "co2_partial_pa":
        return f"Carbon dioxide in the air, at a partial pressure of at least {_pressure(gate['min'])}"
    if quantity == "band_light_vs_earth":
        return f"At least {gate['min']:.0%} of Earth's light in the band the class eats"
    if quantity == "share_o2" and "lt" in gate:
        return f"Free oxygen below {gate['lt']:.0%} of the air"
    if quantity == "share_o2":
        return f"Oxygen at least {gate['min']:.1%} of the air"
    if quantity == "acid_medium":
        return "A pool of liquid sulfuric acid at the ground" if gate["min"] >= 2.0 else "A cloud layer or pool of liquid sulfuric acid"
    if quantity == "share_co2":
        return f"Carbon dioxide at least {gate['min']:.0%} of the air"
    if quantity == "co2_phase":
        return "The carbon dioxide is " + " or ".join(gate["in"])
    if quantity == "ammonia_water_mass_fraction":
        return f"Ammonia at least {gate['min']:.0%} of the water's mass"
    if quantity == "silicon_rich_crust":
        return "A silicon-rich crust"
    return None


def _window_words(window):
    quantity, low, high = window["quantity"], window["low"], window["high"]
    if quantity == "habitat_temperature_k":
        return f"Where the life lives (it can be underground): {_number(low)} to {_number(high)} K"
    if quantity == "temperature_k":
        return f"Surface temperature {_number(low)} to {_number(high)} K"
    if quantity == "pressure_pa":
        return f"Air pressure {_pressure(low)} to {_pressure(high)}"
    return None


def needs(row):
    """What a class needs, as plain phrases: its gates, its temperature and pressure windows and the supply it must have."""
    found = [words for words in (_gate_words(gate) for gate in row.gates) if words]
    found += [words for words in (_window_words(window) for window in row.windows) if words]
    for supply in row.supplies:
        words = SUPPLY_WORDS.get(supply["quantity"])
        if words:
            found.append(words[0].upper() + words[1:])
    return found


def _energy_words(row):
    if row.energy["source"] == "light":
        low, high = row.energy["band_nm"]
        return f"Starlight between {_number(low)} and {_number(high)} nm"
    return "Chemical energy: a fuel and an oxidant meeting"


def _speculative(number):
    """`(tier, level)` of the class's air chemistry and whether it is speculative (level L3 or L4)."""
    entry = life_air.CLASSES[str(number)]
    return entry["tier"], entry["tweak"], entry["tweak"] in ("L3", "L4")


def _class(number, href):
    rows = family(number)
    if not rows:
        return None
    name = family_name(number)
    tier, level, speculative = _speculative(number)
    count = {1: "one form", 2: "two variants", 3: "three variants"}.get(len(rows), f"{len(rows)} variants")
    title = f"Class {rows[0].id}" if len(rows) == 1 else f"The Class {number} family"
    intro = [f"{title} is one way of making a living. It has {count}, each described below with what it drinks, eats, breathes and gives off."]
    sections = []
    for row in rows:
        lines = life_class_text.lines(row.id)
        facts = [("Class", row.id), ("Variant", _variant_name(row)) if len(rows) > 1 else None,
                 ("Drinks", lines["drinks"]), ("Eats", lines["eats"]), ("Breathes", lines["breathes"]), ("Gives off", lines["gives_off"]),
                 ("Liquid it needs", SOLVENT_WORDS.get(row.solvent["medium"], row.solvent["medium"])), ("Energy", _energy_words(row))]
        facts = [fact for fact in facts if fact]
        facts += [("Needs", needs(row) or ["Nothing beyond its liquid and its energy"]),
                  ("How far it can go", row.reaches),
                  ("Confidence it can exist", f"{row.tier}: {TIER_NAMES[row.tier]}")]
        sections.append(_section(f"class-{row.id}", f"Class {row.id}" + (f": {_variant_name(row)}" if len(rows) > 1 else ""), facts=facts))
    sections.append(_section("stages", "What the stages look like", table=_table(
        ["Step", "Name", "What it means"],
        [[str(rung), *_class_stage(number, rung)] for rung in life_class_text.RUNGS])))
    note = (f"Speculative chemistry: how this class changes the air is an educated guess. Its source tier is {tier} and its adjustment level is {level}."
            if speculative else f"How this class changes the air has source tier {tier} and adjustment level {level}.")
    sections.append(_section("sources", "How sure we are", [
        note, f"Tier {tier} means: {TIER_NAMES[tier].lower()}. Level {level} means: {LEVEL_NAMES[level].lower()}.",
        "The tiers run from T1, read in a source, to T9, chosen to give the outcome wanted. The levels run from L0, as the source gives it, to L4, fitted with no independent source.",
    ]))
    return f"{family_label(number)}: {name}", f"What the {name.lower()} drink, eat, breathe and give off, and where they can live.", intro, sections


def _class_stage(number, rung):
    entry = life_class_text.DOCUMENT["classes"][str(number)].get(str(rung)) or life_class_text.DOCUMENT["rungs"][str(rung)]
    return entry["name"], entry["description"].replace("Class {n} life", "life of this class")


def _classes(href):
    rows = []
    for number in sorted({row.number for row in life_classes.CLASSES}):
        members = family(number)
        said = life_class_text.lines(members[0].id)
        rows.append([(f"Class {number}" if len(members) == 1 else f"Classes {', '.join(m.id for m in members)}", href("classes", number)),
                     family_name(number), said["drinks"], said["eats"], members[0].reaches])
    intro = [
        "A life class is a way of making a living. The number is the way: what the life drinks, eats, breathes and gives off.",
        "A letter after the number is a variant of the same way, such as a different pigment or a different home.",
    ]
    sections = [
        _section("classes", "The nine life classes", table=_table(["Class", "Name", "Drinks", "Eats", "How far it can go"], rows)),
        _section("sources", "Speculative chemistry", [
            "Every class has a source tier, from T1 (read in a source) to T9 (chosen to give the outcome wanted). Each number that changes the air also has an adjustment level, from L0 to L4.",
            "A class whose air chemistry is at level L3 or L4 is marked speculative on the planet page, with its tier and level.",
        ], table=_table(["Class", "Air chemistry tier", "Adjustment level", "Speculative"],
                        [[family_label(number), f"{_speculative(number)[0]}: {TIER_NAMES[_speculative(number)[0]]}",
                          f"{_speculative(number)[1]}: {LEVEL_NAMES[_speculative(number)[1]]}", "Yes" if _speculative(number)[2] else "No"]
                         for number in sorted({row.number for row in life_classes.CLASSES})])),
    ]
    return "Life classes", "The nine life classes and what sets them apart.", intro, sections


BUILDERS = {"kinds": _kinds, "abundance": _abundance, "trend": _trend, "security": _security, "suitability": _suitability,
            "stages": _stages, "classes": _classes}


def page(topic, number=None, href=None):
    """
    The page of a life help topic (`TOPICS`) or, with `number`, of the life class of that number, as described at the top of this module;
    `None` for a topic or class that does not exist. `href(topic, number=None)` gives the address of a page (a plain path by default).
    """
    link = href or (lambda slug, number=None: f"/help/life/{slug}" + (f"/{number}" if number is not None else ""))
    if number is not None:
        built = _class(number, link) if topic == "classes" else None
    else:
        builder = BUILDERS.get(topic)
        built = builder(link) if builder else None
    if built is None:
        return None
    title, description, intro, sections = built
    return {"topic": topic, "number": number, "title": title, "description": description, "intro": intro, "sections": sections,
            "topics": [{"slug": slug, "title": name, "summary": summary, "url": link(slug)} for slug, name, summary in TOPICS]}


def help_children():
    """The life help pages (`TOPICS`, the first being the Life page itself) and one page for each life class number, as the pages under Life."""
    found = [Child("life_help_topic", {"topic": slug}, title, summary) for slug, title, summary in TOPICS if slug != "kinds"]
    found += [Child("life_help_class", {"number": number}, label, summary) for number, label, summary in class_numbers()]
    return found


def link_terms():
    """The words a description may link to a life help page: each life class id ("4c", "2") and each class family ("Class 3 family") links to
    the page of its number. A bare number ("2") is a term only after "Class "; an id with a letter ("4c") is one anywhere."""
    found = [Term(row.id, "life_help_class", {"number": row.number}, "" if row.id[-1].isalpha() else "Class ") for row in life_classes.CLASSES]
    found += [Term(family_label(number), "life_help_class", {"number": number}) for number, _label, _summary in class_numbers()
              if len(family(number)) > 1]
    return found


help_registry.register_category("worlds", "Worlds and life", "How a visitor's gear is judged on a world and how life on it is described.", order=20)
help_registry.register_topic("life_help_page", "Life", "worlds", "The kinds of life, the scales a planet page uses, the stages and the life classes.",
                             order=40, children=help_children)
help_registry.register_terms("life", link_terms)
