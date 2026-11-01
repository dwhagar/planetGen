# planetgen/web/lib/classref.py

"""
The class reference catalog behind the `/classes` pages
(`web/class_pages.py`): every kind of class the generator assigns (star
spectral and luminosity classes, planet classes, nebula and supernova
remnant classes, asteroid field letters, black hole and rogue planet
mass classes, comet period classes), built from the generator's own
tables (`planetgen.tuning`, `physical_constants`,
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
import logging
import math

from markupsafe import Markup, escape

from planetgen.physics import brown_dwarf_types, constants as phys, planet_classes, small_body_composition
from planetgen import tuning as pc
from planetgen.datasets import mamajek_sequence
from planetgen.datasets.fetch import DataMissing
from planetgen.generation import forced_star
from planetgen.generation.comet import PERIOD_CLASS_LABELS
from planetgen.generation.phenomena.nebula import NEBULA_CLASS_LETTERS, REMNANT_CLASS_LETTERS
from planetgen.generation.star import STAR_TYPE_PATTERN
from planetgen.physics.stellar_evolution import YERKES_CLASS_NAMES
from planetgen.util import chem
from planetgen.web.lib import help_registry

LUMINOSITY_CLASS_ALIASES = {"D": "VII"}
"""dict: Yerkes codes that are another spelling of a class
(`YERKES_CLASS_NAMES` marks "D" as an alias of "VII"); they get no page
of their own and link to the class they stand for."""

PLANET_TYPE_LABELS = {"t": "Terrestrial", "g": "Gas giant"}
"""dict: `PLANET_CLASSES[...]["type"]` -> label."""

PLANET_ZONE_LABELS = {"h": "Hot zone", "e": "Habitable zone (ecosphere)", "c": "Cold zone",
                      "r": "Interstellar space (rogue planet)"}
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

def _is_scientific(value):
    magnitude = abs(value)
    return value != 0 and (magnitude >= 1e4 or magnitude < 1e-3)


def number(value, scientific=None):
    """
    A number as readable plain text: thousands separators, at most three
    significant figures after the point, and "a × 10ⁿ" for anything at or
    above ten thousand (UX.20: 5+ whole digits) or below a thousandth.
    `scientific` forces the style (`number_range` gives both ends of a range
    the one style, so a range never reads "5,000 to 1 × 10⁴").

    >>> number(150), number(0.08), number(20_000), number(1e-5)
    ('150', '0.08', '2 × 10⁴', '1 × 10⁻⁵')
    """
    if value == 0:
        return "0"
    magnitude = abs(value)
    if _is_scientific(value) if scientific is None else scientific:
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
    """`"low to high unit"`, each end through `number`, both in scientific notation when either end needs it and neither is below 1."""
    scientific = (_is_scientific(low) or _is_scientific(high)) and abs(low) >= 1 and abs(high) >= 1
    text = f"{number(low, scientific or None)} to {number(high, scientific or None)}"
    return f"{text} {unit}" if unit else text


def percent(share):
    """A share in 0-1 as a percentage, e.g. `"7.6%"` or `"0.0001%"`."""
    return f"{number(share * 100)}%"


def _sentence(text):
    """`text` with its first letter capitalized (the tables store
    mid-sentence phrases)."""
    return text[:1].upper() + text[1:] if text else text


def _contents(text):
    """A cloud's contents (`"H+, e-, O2+, N+ and S+"`) as HTML, its formulas written out (UX.149)."""
    return Markup(chem.mark_html(text))


def _label(key):
    """A table key as a label: `"jupiter_family"` -> `"Jupiter family"`."""
    return _sentence(key.replace("_", " "))


# ---------------------------------------------------------------------
# One builder per class type
# ---------------------------------------------------------------------

def _entry(code, name, summary, facts):
    return {"code": code, "name": name, "summary": summary, "facts": facts}


def _star_spectral():
    classes = {}
    for letter, color in phys.SPECTRAL_CLASS_COLORS.items():
        low_k, high_k = mamajek_sequence.letter_range_k(letter)
        found = forced_star.letter_summary(letter)
        classes[letter] = _entry(letter, f"{color} {letter}-type star",
                                 f"{color} stars with surfaces at {number_range(low_k, high_k, 'K')}.", [
            ("Color", color),
            ("Surface temperature", number_range(low_k, high_k, "K")),
            ("Mass (main sequence)", number_range(*found.mass_sol, "solar masses")),
            ("Luminosity (main sequence)", number_range(*found.luminosity_sol, "times the Sun")),
            ("Share of ordinary stars", percent(found.share)),
            ("Subclasses", f"{letter}0 (hottest) to {letter}9 (coolest)"),
        ])
    return {
        "name": "Star Spectral Classes",
        "label": "Spectral class",
        "summary": "A star's color and surface temperature, hottest (O) to coolest (M).",
        "notes": [
            "A star's type reads like G2V. The letter is the spectral class, the digit is the subclass from 0 "
            "(hottest) to 9, and the roman numeral is the luminosity class (see Star Luminosity Classes).",
            "Mass and luminosity here are for main-sequence (class V) stars. Giants and white dwarfs of the "
            "same color differ.",
        ],
        "classes": classes,
    }


def _star_luminosity():
    classes = {}
    for code, name in YERKES_CLASS_NAMES.items():
        if code in LUMINOSITY_CLASS_ALIASES:
            continue
        found = forced_star.class_summary(code)
        facts = [
            ("Name", name),
            ("Luminosity", number_range(*found.luminosity_sol, "times the Sun")),
            ("Mass", number_range(*found.mass_sol, "solar masses")),
        ]
        if found.share is not None:
            facts.append(("Share of ordinary stars", percent(found.share)))
        aliases = [alias for alias, target in LUMINOSITY_CLASS_ALIASES.items() if target == code]
        if aliases:
            facts.append(("Also written", ", ".join(aliases)))
        classes[code] = _entry(code, name, f"Luminosity class {code}: {name.lower()} stars.", facts)
    return {
        "name": "Star Luminosity Classes",
        "label": "Luminosity class",
        "summary": "The Yerkes (MK) luminosity class: a star's size and brightness for its color.",
        "notes": [
            "The roman numeral at the end of a star's type, like the V in G2V, says how big and bright the star is "
            "for its color.",
            "It runs from hypergiants (0) down to white dwarfs (VII).",
        ],
        "classes": classes,
    }


def _planet():
    """The planet and moon classes from the rule rows (`physics/planet_classes.py`, GEN.645): a class is what the finished
    body's conditions say, so a page gives the rule, not a table of ranges the class was drawn from."""
    classes = {}
    ordered = planet_classes.CLASSES
    high_gravity = next((rule["min"] for c in ordered if c.letter == "V" for rule in c.when if rule["quantity"] == "gravity_g"), None)
    for cls in sorted(ordered, key=lambda c: c.letter):
        if cls.letter in classes:     # GEN.1142: T has a giant row and a terrestrial one; the page is the giant row's, with a line for the other
            classes[cls.letter]["facts"].append(("Also", f"{PLANET_TYPE_LABELS[cls.kind]} bodies that meet: " + "; ".join(planet_classes.rule_words(cls))))
            continue
        before = [c.letter for c in ordered if c.kind == cls.kind and ordered.index(c) < ordered.index(cls)]
        rule = planet_classes.rule_words(cls)
        facts = [
            ("Name", cls.name),
            ("Type", PLANET_TYPE_LABELS[cls.kind]),
            ("Rule", "; ".join(rule) if rule else "Every " + PLANET_TYPE_LABELS[cls.kind].lower()
                      + " that no class checked earlier takes"),
            ("Checked after", ", ".join(before) if before else "Nothing: checked first"),
            ("Core class", "Yes" if cls.core else "No"),
        ]
        if cls.letter == "M":     # GEN.1136: the class above takes over at the high-gravity line
            facts.append(("Gravity", f"Any gravity up to the {number(high_gravity)} g line of class V (a super-Earth can be a temperate world)"))
        classes[cls.letter] = _entry(cls.letter, cls.name, cls.sentence, facts)
    return {
        "name": "Planet Classes",
        "label": "Planet class",
        "summary": "One letter for what a planet or moon is like, worked out from its finished conditions.",
        "notes": [
            "A class is never drawn first: once a planet or moon is complete, its class is the first class in the check "
            "order whose rule its mass, gravity, temperature, air, water and starlight satisfy. All the rules must hold.",
            "Habitable is not a class. It is the overall habitability colour (Blue or Green); a world of any class can "
            "score it or not.",
            "Letters no class uses (" + ", ".join(planet_classes.RESERVE) + ") are kept free for new classes.",
            "A rogue planet, adrift with no star, gets a class too, from the classes that list "
            "interstellar space: with no starlight, no class with life qualifies.",
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
            ("Contents", _contents(_sentence(data["species"]))),
            ("Radius", number_range(*data["radius_range_ly"], "ly")),
            *_cloud_facts(data),
        ]
        if data.get("center"):
            facts.append(("At the center", _sentence(data["center"])))
        classes[code] = _entry(code, data["name"], escape(f"{_sentence(data['family'])} nebula: ") + _contents(data["species"]) + ".",
                               facts)
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
            ("Contents", _contents(_sentence(data["species"]))),
            ("Age", number_range(*data["age_range_years"], "years")),
            ("Compact remnant allowed", ", ".join(_COMPACT_LABELS.get(kind, str(kind)) for kind in data["compact"])),
            *_cloud_facts(data),
        ]
        if data.get("center"):
            facts.append(("At the center", _sentence(data["center"])))
        summary = escape(f"{_sentence(data['morphology'])} remnant: ") + _contents(data["species"]) + "."
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
        components_text = small_body_composition.field_family_text(family)
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
        elif code == "supermassive":
            facts.append(("Mass", f"{number(pc.BLACK_HOLE_SUPERMASSIVE_MIN_MASS_SOLAR)} solar masses or more "
                                  "(the galaxy's bulge or stellar mass sets it)"))
        if code in origins:
            facts.append(("Where found", origins[code]))
        facts.append(("Spin", number_range(*pc.BLACK_HOLE_SPIN_RANGE, "(dimensionless)")))
        summary = f"{_label(code)} black holes"
        if code in ranges:
            summary += f", {number_range(*ranges[code], 'solar masses')}"
        elif code == "supermassive":
            summary += f", {number(pc.BLACK_HOLE_SUPERMASSIVE_MIN_MASS_SOLAR)} solar masses or more"
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
        ("Type", f"L from {number(round(brown_dwarf_types.teff_of_type(0)))} K, T from "
                 f"{number(round(brown_dwarf_types.teff_of_type(10)))} K, Y from "
                 f"{number(round(brown_dwarf_types.teff_of_type(20)))} K down (hotter is a late-M dwarf)"),
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


PARABOLIC_COMET_CLASS = "parabolic"
"""str: The comet class code of a comet on an open (parabolic) orbit,
which has no period class (UX.29)."""


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
    # UX.29: a comet on an open (parabolic) orbit has no period class, so
    # it gets a class of its own and every comet row can link one.
    classes[PARABOLIC_COMET_CLASS] = _entry(
        PARABOLIC_COMET_CLASS, "Single-apparition comet",
        "Comets on open (parabolic) orbits: each passes its star once and never returns.", [
            ("Orbit", "Open (parabolic), eccentricity 1"),
            ("Inclination", f"0 to {number(pc.PARABOLIC_COMET_INCLINATION_MAX_DEG)}°"),
        ])
    return {
        "name": "Comet Classes",
        "label": "Comet class",
        "summary": "A star's comets: periodic ones by the length of their orbits, and single-apparition ones.",
        "notes": ["Comets on closed (elliptical) orbits have a period class; comets on open (parabolic) orbits "
                  "are single-apparition comets."],
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


def preload():
    """Builds `catalog()` now so no request does; `False` (the site still boots) when a dataset the star classes read
    has not been downloaded yet: `update.sh` swaps in the new code before its downloads finish, and a worker that
    started in between used to die at boot with every page a 500. The catalog is not cached on a failure, so the
    first /classes request after the download builds it."""
    try:
        catalog()
    except DataMissing as missing:
        logging.getLogger(__name__).warning("the class reference catalog was not built at start: %s", missing)
        return False
    return True


def alphabetical(rows, key="name"):
    """`rows` (dicts) in alphabetical order of `row[key]`, ignoring case: how
    the site lists the class categories (UX.77)."""
    return sorted(rows, key=lambda row: str(row[key]).casefold())


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
    `SystemConfig.STAR_TYPE` (`star.STAR_TYPE_PATTERN`), so a
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


def star_type_split(star_type):
    """
    A stored star type ("G2V Yellow Main Sequence Star") -> `(spectral text, luminosity text, rest)` as written: ("G2",
    "V", " Yellow Main Sequence Star"). The first two are what link to their class pages (UX.182), `rest` is the words
    after the code. `None` where `star_type_classes` has none.
    """
    if star_type_classes(star_type) is None:
        return None
    code, _, tail = star_type.strip().partition(" ")
    letter, subclass, yerkes = STAR_TYPE_PATTERN.fullmatch(code.upper()).groups()
    cut = len(code) - len(yerkes)
    return code[:cut], code[cut:], f" {tail}" if tail else ""


def help_children():
    """The class types, alphabetically, as the pages under Classes in the help index."""
    from planetgen.web.lib.help_registry import Child
    rows = alphabetical([{"name": entry["name"], "slug": slug, "summary": entry["summary"]} for slug, entry in catalog().items()])
    return [Child("class_type_page", {"type_slug": row["slug"]}, row["name"], row["summary"]) for row in rows]


help_registry.register_category("objects", "Objects and classes", "What the objects of the galaxy are, how they are sorted, and where their values come from.", order=10)
help_registry.register_topic("classes", "Classes", "objects", "Every class of star, planet, nebula and other object the survey assigns.",
                             order=10, children=help_children)
