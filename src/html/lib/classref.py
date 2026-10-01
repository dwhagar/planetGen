# html/lib/classref.py

"""
The class reference catalog behind the `/classes` pages
(`web/class_pages.py`): every kind of class the generator assigns (star
spectral and luminosity classes, planet classes, nebula and supernova
remnant classes, asteroid field letters, black hole and rogue planet
mass classes, comet period classes), built from the generator's own
tables (`stellarObjects.program_constants`, `physical_constants`,
`stellarEvolution`, `nebulaData`, `cometData`) so the pages always show
the values the program uses. Nothing here repeats a number or a
description by hand.

`catalog()` builds the whole thing once per process (the site calls it
at startup, `web.init_app`) and returns the cached result after that:

    {type_slug: {"slug", "name", "label", "summary", "notes": [str],
                 "classes": {code: {"code", "name", "summary",
                                    "facts": [(label, text)]}}}}

A type's `name` is its page title ("Planet Classes") and `label` names
one of its classes ("Planet class" M). Every value is plain text; the
templates escape it.

The pages that show a class link it with `class_url_parts` (is this a
real class?) and, for a star, `star_type_classes` ("G2V Yellow Main
Sequence Star" -> ("G", "V")).
"""

import functools
import math

from stellarObjects import physical_constants as phys
from stellarObjects import program_constants as pc
from stellarObjects.cometData import PERIOD_CLASS_LABELS
from stellarObjects.nebulaData import NEBULA_CLASS_LETTERS, REMNANT_CLASS_LETTERS
from stellarObjects.starData import STAR_TYPE_PATTERN
from stellarObjects.stellarEvolution import YERKES_CLASS_NAMES

LUMINOSITY_CLASS_ALIASES = {"D": "VII"}
"""dict: Yerkes codes that are another spelling of a class
(`YERKES_CLASS_NAMES` marks "D" as an alias of "VII"); they get no page
of their own and link to the class they stand for."""

PLANET_TYPE_LABELS = {"t": "Terrestrial", "g": "Gas giant"}
"""dict: `PLANET_CLASSES[...]["type"]` -> label."""

PLANET_ZONE_LABELS = {"h": "Hot zone", "e": "Habitable zone (ecosphere)", "c": "Cold zone"}
"""dict: The `PLANET_CLASSES` zone flags, in order out from the star."""

ROGUE_MASS_CLASS_NAMES = {
    "terrestrial": "Terrestrial", "sub-neptune": "Sub-Neptune", "saturn": "Saturn-mass",
    "jupiter": "Jupiter-mass", "brown-dwarf": "Brown dwarf",
}
"""dict: Display name of each `rogue_planets.mass_bin`
(`program_constants.ROGUE_PLANET_MASS_BIN_CHOICES`)."""

_COMPACT_LABELS = {"neutron_star": "Neutron star", "black_hole": "Black hole", None: "None"}

_SUPERSCRIPT = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")


# ---------------------------------------------------------------------
# Number formatting
# ---------------------------------------------------------------------

def number(value):
    """
    A number as readable plain text: thousands separators, at most three
    significant figures after the point, and "a × 10ⁿ" for anything at or
    above ten thousand (UX.20: 5+ whole digits) or below a thousandth.

    >>> number(150), number(0.08), number(20_000), number(1e-5)
    ('150', '0.08', '2 × 10⁴', '1 × 10⁻⁵')
    """
    if value == 0:
        return "0"
    magnitude = abs(value)
    if magnitude >= 1e4 or magnitude < 1e-3:
        exponent = math.floor(math.log10(magnitude))
        mantissa = f"{value / 10 ** exponent:.2f}".rstrip("0").rstrip(".")
        if mantissa in ("10", "-10"):
            exponent += 1
            mantissa = mantissa[:-1]
        return f"{mantissa} × 10{str(exponent).translate(_SUPERSCRIPT)}"
    if magnitude >= 100:
        return f"{value:,.0f}" if value == round(value) else f"{value:,.1f}".rstrip("0").rstrip(".")
    return f"{value:.3g}"


def number_range(low, high, unit=""):
    """`"low to high unit"`, each end through `number`."""
    text = f"{number(low)} to {number(high)}"
    return f"{text} {unit}" if unit else text


def percent(share):
    """A share in 0-1 as a percentage, e.g. `"7.6%"` or `"0.0001%"`."""
    return f"{number(share * 100)}%"


def _sentence(text):
    """`text` with its first letter capitalized (the tables store
    mid-sentence phrases)."""
    return text[:1].upper() + text[1:] if text else text


def _label(key):
    """A table key as a label: `"jupiter_family"` -> `"Jupiter family"`."""
    return _sentence(key.replace("_", " "))


# ---------------------------------------------------------------------
# One builder per class type
# ---------------------------------------------------------------------

def _entry(code, name, summary, facts):
    return {"code": code, "name": name, "summary": summary, "facts": facts}


def _star_spectral():
    total = sum(pc.SPECTRAL_PROBABILITIES_NORMAL.values())
    classes = {}
    for letter, color in phys.SPECTRAL_CLASS_COLORS.items():
        low_k, high_k = phys.TEMP_RANGES[letter]
        share = pc.SPECTRAL_PROBABILITIES_NORMAL.get(letter, 0) / total
        classes[letter] = _entry(letter, f"{color} {letter}-type star",
                                 f"{color} stars with surfaces at {number_range(low_k, high_k, 'K')}.", [
            ("Color", color),
            ("Surface temperature", number_range(low_k, high_k, "K")),
            ("Mass (main sequence)", number_range(*phys.SPECTRAL_MASS_RANGES[letter], "solar masses")),
            ("Luminosity (main sequence)", number_range(*phys.SPECTRAL_LUMINOSITY_RANGES[letter], "times the Sun")),
            ("Share of ordinary stars", percent(share)),
            ("Subclasses", f"{letter}0 (hottest) to {letter}9 (coolest)"),
        ])
    return {
        "name": "Star Spectral Classes",
        "label": "Spectral class",
        "summary": "A star's color and surface temperature, hottest (O) to coolest (M).",
        "notes": [
            "A star's type reads like G2V: the spectral class letter, a subclass digit from 0 (hottest) to 9, "
            "and the luminosity class (see Star Luminosity Classes).",
            "Mass and luminosity here are for main-sequence (class V) stars; giants and white dwarfs of the "
            "same color differ.",
        ],
        "classes": classes,
    }


def _star_luminosity():
    classes = {}
    for code, name in YERKES_CLASS_NAMES.items():
        if code in LUMINOSITY_CLASS_ALIASES:
            continue
        facts = [
            ("Name", name),
            ("Luminosity", number_range(*phys.YERKES_LUMINOSITY_RANGES[code], "times the Sun")),
            ("Mass", number_range(*phys.YERKES_MASS_CONSTRAINTS[code], "solar masses")),
        ]
        aliases = [alias for alias, target in LUMINOSITY_CLASS_ALIASES.items() if target == code]
        if aliases:
            facts.append(("Also written", ", ".join(aliases)))
        classes[code] = _entry(code, name, f"Luminosity class {code}: {name.lower()} stars.", facts)
    return {
        "name": "Star Luminosity Classes",
        "label": "Luminosity class",
        "summary": "The Yerkes (MK) luminosity class: a star's size and brightness for its color.",
        "notes": [
            "The roman numeral at the end of a star's type (the V in G2V) says how big and bright the star is "
            "for its color, from hypergiants (0) down to white dwarfs (VII).",
        ],
        "classes": classes,
    }


def _planet():
    classes = {}
    for code, data in pc.PLANET_CLASSES.items():
        type_label = PLANET_TYPE_LABELS.get(data["type"], data["type"])
        habitable = code in pc.HABITABLE_PLANET_CLASSES
        zones = [label for key, label in PLANET_ZONE_LABELS.items() if data.get(key)]
        facts = [
            ("Description", _sentence(data["description"])),
            ("Type", type_label),
            ("Composition", _sentence(data["composition"])),
            ("Atmosphere", _sentence(data["atmosphere"]) if data.get("atmosphere") else "None"),
            ("Radius", number_range(*data["radius_range"], "km")),
            ("Forms in", ", ".join(zones) if zones else "None"),
            ("Habitable", "Yes" if habitable else "No"),
        ]
        if data.get("life_chemical"):
            facts.append(("Life's light-harvesting pigments", ", ".join(data["life_chemical"])))
        if code in pc.PLANET_CLASS_MAX_LIFE_STAGE:
            facts.append(("Most advanced life", _label(pc.PLANET_CLASS_MAX_LIFE_STAGE[code])))
        name = f"{type_label}, habitable" if habitable else type_label
        classes[code] = _entry(code, name, _sentence(data["description"]) + ".", facts)
    return {
        "name": "Planet Classes",
        "label": "Planet class",
        "summary": "One letter for what a planet or moon is like: its makeup, air and size.",
        "notes": [
            "Each planet and moon gets one class letter. A class forms only in the zones listed for it: "
            "the hot zone near the star, the habitable zone (ecosphere), or the cold zone beyond.",
        ],
        "classes": classes,
    }


def _cloud_facts(data):
    return [
        ("Density", number_range(*data["density_range_cm3"], "particles/cm³")),
        ("Gas temperature", number_range(*data["temperature_range_k"], "K")),
        ("Extinction", number_range(*data["extinction_range_av"], "magnitudes (visual)")),
    ]


def _nebula():
    classes = {}
    for code in NEBULA_CLASS_LETTERS:
        data = pc.NEBULA_CLASSES[code]
        facts = [
            ("Name", data["name"]),
            ("Family", _sentence(data["family"])),
            ("Contents", _sentence(data["species"])),
            ("Radius", number_range(*data["radius_range_ly"], "ly")),
            *_cloud_facts(data),
        ]
        if data.get("center"):
            facts.append(("At the center", _sentence(data["center"])))
        classes[code] = _entry(code, data["name"], f"{_sentence(data['family'])} nebula: {data['species']}.", facts)
    return {
        "name": "Nebula Classes",
        "label": "Nebula class",
        "summary": "Clouds of gas and dust, by family and contents.",
        "notes": [
            "Nebulae fall into five families: " + ", ".join(pc.NEBULA_FAMILIES) + ". "
            "The class letter says what is in the cloud. Ranges are what a new nebula of the class is drawn from.",
        ],
        "classes": classes,
    }


def _supernova_remnant():
    classes = {}
    for code in REMNANT_CLASS_LETTERS:
        data = pc.NEBULA_CLASSES[code]
        facts = [
            ("Name", data["name"]),
            ("Morphology", _sentence(data["morphology"])),
            ("Contents", _sentence(data["species"])),
            ("Age", number_range(*data["age_range_years"], "years")),
            ("Compact remnant allowed", ", ".join(_COMPACT_LABELS.get(kind, str(kind)) for kind in data["compact"])),
            *_cloud_facts(data),
        ]
        if data.get("center"):
            facts.append(("At the center", _sentence(data["center"])))
        summary = f"{_sentence(data['morphology'])} remnant: {data['species']}."
        classes[code] = _entry(code, data["name"], summary, facts)
    return {
        "name": "Supernova Remnant Classes",
        "label": "Remnant class",
        "summary": "The expanding wreckage of a supernova, by shape and age.",
        "notes": [
            "A remnant's shape is a shell, a plerion (filled by a pulsar wind) or a composite of both. "
            "Its radius follows from its age.",
        ],
        "classes": classes,
    }


def _asteroid_field():
    densities = {}  # letter -> (family, data, [density, ...])
    for family, data in pc.ASTEROID_FIELD_COMPOSITIONS.items():
        for density, letter in data["letters"].items():
            densities.setdefault(letter, (family, data, []))[2].append(density)
    classes = {}
    for letter in sorted(densities):
        family, data, density_list = densities[letter]
        components = data["components"]
        if components is None:
            components_text = "Any asteroid material: " + ", ".join(pc.ASTEROID_COMPONENTS)
        else:
            components_text = _sentence(", ".join(components))
        classes[letter] = _entry(letter, f"{_sentence(family)}, {' or '.join(density_list)}",
                                 _sentence(data["description"]) + ".", [
            ("Composition family", _sentence(family)),
            ("Density", _sentence(", ".join(density_list))),
            ("Description", _sentence(data["description"])),
            ("Made of", components_text),
        ])
    low_digit, high_digit = pc.ASTEROID_FIELD_SIZE_DIGIT_RANGE
    sizes = "; ".join(f"{digit}: {number_range(10 ** digit, 10 ** (digit + 1), 'AU')}"
                      for digit in range(low_digit, high_digit + 1))
    example = next(iter(pc.ASTEROID_FIELD_COMPOSITIONS.values()))["letters"]["dense"]
    return {
        "name": "Asteroid Field Classes",
        "label": "Asteroid field class",
        "summary": "A letter for what a free-floating field is made of and how dense it is, then a size digit.",
        "notes": [
            "An asteroid field's class is a letter and a digit. The letter (below) is its composition and "
            "density; the digit is its size.",
            f"Size digit d means a radius of 10^d to 10^(d+1) AU, from {low_digit} to {high_digit}: {sizes}. "
            f"So a full class reads like “{example}{min(low_digit + 2, high_digit)}”.",
        ],
        "classes": classes,
    }


def _black_hole():
    ranges = {
        "stellar": pc.BLACK_HOLE_MASS_RANGE_SOLAR,
        "intermediate": pc.BLACK_HOLE_INTERMEDIATE_MASS_RANGE_SOLAR,
        "supermassive": pc.BLACK_HOLE_SUPERMASSIVE_MASS_RANGE_SOLAR,
    }
    origins = {
        "stellar": f"The collapsed core of a massive star; {percent(1 - pc.BLACK_HOLE_INTERMEDIATE_MASS_CHANCE)} "
                   "of free-floating black holes",
        "intermediate": f"Rare: {percent(pc.BLACK_HOLE_INTERMEDIATE_MASS_CHANCE)} of free-floating black holes",
        "supermassive": "The quiet black hole at a galaxy's center (when it is not a quasar)",
    }
    classes = {}
    for code in pc.BLACK_HOLE_MASS_CLASSES:
        facts = []
        if code in ranges:
            facts.append(("Mass", number_range(*ranges[code], "solar masses")))
        if code in origins:
            facts.append(("Where found", origins[code]))
        facts.append(("Spin", number_range(*pc.BLACK_HOLE_SPIN_RANGE, "(dimensionless)")))
        summary = f"{_label(code)} black holes"
        if code in ranges:
            summary += f", {number_range(*ranges[code], 'solar masses')}"
        classes[code] = _entry(code, f"{_label(code)} black hole", summary + ".", facts)
    return {
        "name": "Black Hole Classes",
        "label": "Black hole class",
        "summary": "Black holes by mass.",
        "notes": ["A black hole's class is its mass range."],
        "classes": classes,
    }


def _rogue_planet():
    classes = {}
    for code, (low, high, rate) in pc.ROGUE_PLANET_MASS_BINS.items():
        name = ROGUE_MASS_CLASS_NAMES.get(code, _label(code))
        classes[code] = _entry(code, name, f"{name} planets with no star, {number_range(low, high, 'Earth masses')}.", [
            ("Mass", number_range(low, high, "Earth masses")),
            ("How common", f"About {number(rate)} per star in the galaxy"),
        ])
    low, high = pc.ROGUE_BROWN_DWARF_MASS_RANGE_JUPITER
    name = ROGUE_MASS_CLASS_NAMES["brown-dwarf"]
    summary = f"A failed star, {number_range(low, high, 'Jupiter masses')}."
    classes["brown-dwarf"] = _entry("brown-dwarf", name, summary, [
        ("Mass", number_range(low, high, "Jupiter masses")),
    ])
    return {
        "name": "Rogue Planet Classes",
        "label": "Rogue planet class",
        "summary": "Planets drifting between the stars, by mass.",
        "notes": [
            f"A rogue planet's mass class is drawn by how common each is (about "
            f"{number(sum(rate for _low, _high, rate in pc.ROGUE_PLANET_MASS_BINS.values()))} per star in all), "
            "then its mass within the class.",
        ],
        "classes": classes,
    }


def _comet():
    total = sum(data["weight"] for data in pc.COMET_PERIOD_CLASSES.values())
    classes = {}
    for code, data in pc.COMET_PERIOD_CLASSES.items():
        label = _sentence(PERIOD_CLASS_LABELS.get(code, _label(code)))
        classes[code] = _entry(code, f"{label} comet",
                               f"{label} comets, orbits of {number_range(*data['period_range_years'], 'years')}.", [
            ("Orbital period", number_range(*data["period_range_years"], "years")),
            ("Eccentricity", number_range(*data["eccentricity_range"])),
            ("Inclination", f"0 to {number(data['inclination_max_deg'])}°"),
            ("Share of periodic comets", percent(data["weight"] / total)),
        ])
    return {
        "name": "Comet Classes",
        "label": "Comet class",
        "summary": "A star's periodic comets, by the length of their orbits.",
        "notes": ["Only comets on closed (elliptical) orbits have a period class; single-apparition comets have none."],
        "classes": classes,
    }


_BUILDERS = (
    ("star-spectral", _star_spectral),
    ("star-luminosity", _star_luminosity),
    ("planet", _planet),
    ("nebula", _nebula),
    ("supernova-remnant", _supernova_remnant),
    ("asteroid-field", _asteroid_field),
    ("black-hole", _black_hole),
    ("rogue-planet", _rogue_planet),
    ("comet", _comet),
)


@functools.lru_cache(maxsize=None)
def catalog():
    """
    Every class type and its classes (see the module docstring), built
    once from the generator's tables and cached for the process. Treat
    the result as read-only.
    """
    types = {}
    for slug, builder in _BUILDERS:
        data = builder()
        data["slug"] = slug
        types[slug] = data
    return types


def class_type(type_slug):
    """One type's entry from `catalog()`, or `None`."""
    return catalog().get(type_slug)


def class_entry(type_slug, code):
    """One class's entry from `catalog()`, or `None`."""
    entry = class_type(type_slug)
    return entry["classes"].get(code) if entry else None


def class_url_parts(type_slug, code):
    """
    `(type_slug, code)` for a class page that exists, else `None`. A
    luminosity class alias ("D") resolves to its class ("VII"), and an
    asteroid field's full class ("C3") to its letter.
    """
    if code is None:
        return None
    code = str(code)
    if type_slug == "star-luminosity":
        code = LUMINOSITY_CLASS_ALIASES.get(code, code)
    elif type_slug == "asteroid-field":
        code = code.rstrip("0123456789")
    if class_entry(type_slug, code) is None:
        return None
    return type_slug, code


def star_type_classes(star_type):
    """
    A stored star type ("G2V Yellow Main Sequence Star", "B0IA", "O5VII",
    "M3D") -> `(spectral letter, luminosity class code)`, with "D" read
    as "VII" (both are white dwarfs). Parsed like a forced
    `SystemConfig.STAR_TYPE` (`starData.STAR_TYPE_PATTERN`), so a
    subclass digit is required. `None` for anything else.
    """
    words = (star_type or "").split()
    if not words:
        return None
    match = STAR_TYPE_PATTERN.fullmatch(words[0].upper())
    if not match:
        return None
    letter, _subclass, yerkes = match.groups()
    return letter, LUMINOSITY_CLASS_ALIASES.get(yerkes, yerkes)
