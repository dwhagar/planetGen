# planetgen/physics/class_fit.py

"""
Planet and moon classes as percent matches (GEN.1205, GEN.1206; design: reports/class-m-and-trek-classes-2026-10-31.md sections 4 to 6)
========================================================================================================================================

`planet_classes.classify` walks the rows of `planet_classes.json` in order and returns the first whose conditions hold, so a row
quietly relies on the rows before it having said no. This module is the other way to name a class: every class is described by its
own positive conditions (`physics/data/planet_class_fit.json`), each condition is graded, and a body's class is the one it fits best.

* **Gates (kinds, not degrees):** giant versus terrestrial (`body_type`) and, for terrestrial classes, air versus no air (pressure
  at or above `gate_pressure_pa`, 10 Pa). A class is only scored against bodies that pass its gates.
* **Condition membership:** 1 inside the edge; across the shoulder beyond it a straight line down to 0; 0 past that. A `min` edge e
  with shoulder s is 1 at e or more, 0 at e - s or less. An `lt` edge e is 1 below e, 0 at e + s or more. `equals` is 0 or 1. A
  quantity the body does not store scores 0.
* **Domains (GEN.1221):** each quantity belongs to one domain of the class-by-domain matrix (`physics/class_domains.py`); `domain_scores` grades a body in each.
* **Class match:** the geometric mean of the memberships of the class's conditions, so one condition badly missed pulls the match
  down and a condition missed outright makes it 0.
* **Best class:** the highest match; a tie goes to the class with more conditions (the narrower description), then to the letter.

`classify` (GEN.1207) is the one that decides a stored class: the best class, its match, the runner-up and its match, and the digit of the class code.
A body that no class fits (every class has a condition missed outright) is a hole in the class set and is settled in two steps, both stand-ins (GEN.725, section 6.189):
the conditions of the gravity and the air-amount domains are set aside (`LIMIT_DOMAINS`) and the best class of what is left among the families (`LIMITED_FAMILIES`) is the letter, with the digit 9
(limited by another domain, `class_code`); a body that still fits no class takes the class it is nearest to (the highest geometric mean with every membership
lifted to `HOLE_FLOOR`, among the ordinary worlds `NEAREST_FAMILIES`), with a stored match of 0. `fit` and `best_class` are also used by the shadow scorer `src/tests/class_fit_cli.py` (GEN.1206).

Use:

    best_class(body)         # ("M", 0.97, "O", 0.41): letter, match, runner-up letter, runner-up match
    scores(body)             # {"M": 0.97, "O": 0.41, ...} for the classes whose gates the body passes
"""

import json
import math
import os
from collections import namedtuple
from collections.abc import Mapping

from planetgen.physics import class_domains, elements, planet_classes

DATA_FILE = os.path.join(elements.DATA_DIRECTORY, "planet_class_fit.json")


class FitError(ValueError):
    """A fit table that is refused, or a body that cannot be scored."""


class FitClass:
    """One class of the fit table: a letter, a kind, the air gate and the graded conditions."""

    __slots__ = ("letter", "name", "kind", "air", "when", "note")

    def __init__(self, row):
        self.letter = row["letter"]
        self.name = row["name"]
        self.kind = row["kind"]
        self.air = row["air"]
        self.when = tuple(row["when"])
        self.note = row.get("note", "")


def _load():
    with open(DATA_FILE, encoding="utf-8") as handle:
        document = json.load(handle)
    classes = tuple(FitClass(row) for row in document["classes"])
    shoulders = document["shoulders"]
    for cls in classes:
        if not cls.when:
            raise FitError(f"class {cls.letter} ({cls.kind}) has no conditions of its own")
        if (cls.kind == planet_classes.GIANT) != (cls.air is None):
            raise FitError(f"class {cls.letter} ({cls.kind}): the air gate is null for giants and true or false for terrestrial bodies")
        for condition in cls.when:
            keys = set(condition) - {"quantity", "min", "lt", "equals"}
            if keys or len(condition) != 2:
                raise FitError(f"class {cls.letter}: condition {condition!r} is not one quantity with one of min, lt or equals")
            if "equals" not in condition and condition["quantity"] not in shoulders:
                raise FitError(f"class {cls.letter}: {condition['quantity']!r} has no shoulder")
            try:
                class_domains.domain_of(condition["quantity"])
            except KeyError:
                raise FitError(f"class {cls.letter}: {condition['quantity']!r} belongs to no domain of the class-by-domain matrix") from None
    if len({(c.letter, c.kind) for c in classes}) != len(classes):
        raise FitError("a class letter is used twice for one kind")
    return document, classes


DOCUMENT, CLASSES = _load()
"""tuple: The `FitClass` rows; `DOCUMENT` the whole file."""

GATE_PRESSURE_PA = DOCUMENT["gate_pressure_pa"]
SHOULDERS = DOCUMENT["shoulders"]
LETTERS = tuple(sorted({c.letter for c in CLASSES}))


def shoulder_width(quantity, edge):
    """The width of the shoulder beyond `edge` of `quantity`, in the quantity's unit."""
    entry = SHOULDERS[quantity]
    if "absolute" in entry:
        return entry["absolute"]
    return max(entry["relative"] * abs(edge), entry.get("floor", 0.0))


def membership(condition, found):
    """The membership (0 to 1) of the named quantities `found` in one condition of a class."""
    value = found.get(condition["quantity"])
    if value is None or value != value:     # a missing or NaN quantity is no match: every comparison with NaN is false, which would otherwise read as a pass (GEN.1593)
        return 0.0
    if "equals" in condition:
        return 1.0 if value == condition["equals"] else 0.0
    if "min" in condition:
        edge = condition["min"]
        width = shoulder_width(condition["quantity"], edge)
        if value >= edge:
            return 1.0
        return 0.0 if width <= 0.0 or value <= edge - width else (value - (edge - width)) / width
    edge = condition["lt"]
    width = shoulder_width(condition["quantity"], edge)
    if value < edge:
        return 1.0
    return 0.0 if width <= 0.0 or value >= edge + width else 1.0 - (value - edge) / width


def match(cls, found, relax=(), floor=0.0):
    """The percent match of `cls` for the named quantities `found`, as a fraction (the geometric mean of its memberships). The conditions in the domains
    `relax` are left out (the class is judged on the rest); with a `floor` above 0 every membership is lifted to it, so a missed condition lowers the
    mean and does not zero it."""
    total, count = 0.0, 0
    for condition in cls.when:
        if relax and class_domains.domain_of(condition["quantity"]) in relax:
            continue
        m = max(membership(condition, found), floor)
        if m <= 0.0:
            return 0.0
        total += math.log(m)
        count += 1
    return math.exp(total / count) if count else 1.0


def passes_gates(cls, found):
    """Whether the body is of the class's kind and, for a terrestrial class, on the right side of the air gate."""
    if cls.kind != found["body_type"]:
        return False
    if cls.air is None:
        return True
    if found["pressure_pa"] != found["pressure_pa"]:        # a NaN pressure is on neither side of the gate
        return False
    return (found["pressure_pa"] >= GATE_PRESSURE_PA) == cls.air


def _found(body):
    if isinstance(body, Mapping) and "body_type" in body and "mass_earth" in body:
        return body
    values = body if isinstance(body, Mapping) else planet_classes.values_from_body(body)
    found = planet_classes.quantities(values)
    found.update(class_domains.derived(values))        # GEN.1221: the regional, day, night, wet-bulb, humidity, groundwater and phase quantities
    return found


def _ranked(found, relax=(), floor=0.0, only=None):
    scored = [(match(cls, found, relax, floor), cls) for cls in CLASSES if passes_gates(cls, found) and (only is None or cls.letter in only)]
    scored.sort(key=lambda pair: (-pair[0], -len(pair[1].when), pair[1].letter))
    return scored


def class_scores(body):
    """Every class whose gates the body passes, as a list of `(match, FitClass)` best first (ties: more conditions, then letter)."""
    return _ranked(_found(body))


def domain_scores(body, letter):
    """`{domain: score}` of a planet, moon or stored-column mapping for the class `letter` (a terrestrial class where both kinds use it), by domain (`class_domains`)."""
    found = _found(body)
    for cls in CLASSES:
        if cls.letter == letter and passes_gates(cls, found):
            return class_domains.domain_scores(cls, found, membership)
    raise FitError(f"no class {letter!r} that the body passes the gates of")


def scores(body):
    """`{letter: match}` for the classes the body passes the gates of."""
    return {cls.letter: m for m, cls in class_scores(body)}


def best_class(body):
    """
    `(letter, match, runner_up_letter, runner_up_match)` of a planet, moon or stored-column mapping. The runner-up is `None`, 0.0 when
    only one class passes the gates. When no class matches at all (best match 0) the letter is `None`: the body is a hole in the class
    set, and what the classifier does with it is GEN.1207's decision.
    """
    ranked = class_scores(body)
    if not ranked:
        raise FitError("no class passes the gates of the body")
    best_m, best = ranked[0]
    if best_m <= 0.0:
        return None, 0.0, None, 0.0
    if len(ranked) == 1:
        return best.letter, best_m, None, 0.0
    return best.letter, best_m, ranked[1][1].letter, ranked[1][0]


LIMIT_DOMAINS = ("gravity", "air_amount")
"""tuple: The domains a hole body is judged without (GEN.1207, registry HC1): a body whose class is found only with its gravity and its air amount set aside is limited by
another domain, digit 9 of the class code."""

HOLE_FLOOR = 0.01
"""float: The lift of every membership when a body that fits no class takes the class it is nearest to (GEN.1207, registry HC2)."""

NEAREST_FAMILIES = ("M", "R", "P", "O", "H", "F", "G")
"""tuple: The classes a body that fits none takes its nearest from (registry HC2): the ordinary worlds with air. A class with a narrow window (Y, Z, U, X) or one defined by an edge
(B, V, K, N, L, E) is never the nearest of a body that missed all of them."""

LIMITED_FAMILIES = ("M", "R", "P", "O", "H", "F", "G", "Y", "U", "X")
"""tuple: The classes a body can be limited in (`class_code.LETTERS`, less Z): a class defined by its gravity or its air amount (B, V, K, N, Z, L, E) has those conditions as its
own edges, so setting them aside would match it to everything."""

LIMITED_DIGIT = 9
"""int: The digit of the class code of a body whose class is found with its gravity and air amount set aside (`class_code.POSITION`, "limited by another domain")."""

Classification = namedtuple("Classification", ("letter", "kind", "match", "runner_up", "runner_up_match", "digit"))
"""namedtuple: The class of a body: its letter and kind, the match (0 when no class fits it), the runner-up letter (`None` when no other class matches at all) and its
match, and the digit of the class code (0, or 9 when the body is limited by another domain)."""


def classify(body):
    """
    The `Classification` of a planet, moon or stored-column mapping: the best class, or for a hole body (see the module note) the class its other domains give
    with digit 9, or the nearest class with a match of 0.

    Raises:
        FitError: No class passes the gates of the body.
        planet_classes.ClassError: The body stores a quantity that is NaN, or no mean temperature for a world with air.
    """
    found = _found(body)
    not_a_number = sorted(name for name, value in found.items() if isinstance(value, float) and value != value)
    if not_a_number:
        raise planet_classes.ClassError(f"the body stores a quantity that is not a number ({', '.join(not_a_number)}); a class is never fitted to it (GEN.1593)")
    if found["body_type"] == planet_classes.TERRESTRIAL and found["pressure_pa"] >= GATE_PRESSURE_PA and found["temperature_k"] is None:
        raise planet_classes.ClassError("the classes read the mean surface temperature of a world with air, which the body does not store")
    ranked = _ranked(found)
    if not ranked:
        raise FitError("no class passes the gates of the body")
    digit = 0
    if ranked[0][0] <= 0.0:
        relaxed = _ranked(found, relax=LIMIT_DOMAINS, only=LIMITED_FAMILIES)
        if relaxed and relaxed[0][0] > 0.0:
            ranked, digit = relaxed, LIMITED_DIGIT
        else:
            nearest = _ranked(found, floor=HOLE_FLOOR, only=NEAREST_FAMILIES) or _ranked(found, floor=HOLE_FLOOR)
            ranked = [(0.0, cls) for _m, cls in nearest]
    best_m, best = ranked[0]
    runner = ranked[1] if len(ranked) > 1 and ranked[1][0] > 0.0 else None       # a class that matches 0 percent is no runner-up
    return Classification(best.letter, best.kind, best_m, None if runner is None else runner[1].letter, 0.0 if runner is None else runner[0], digit)


def rule_words(cls, measure=None):
    """The conditions of a `FitClass` as a list of plain phrases (`planet_classes.condition_words`, which says what `measure` is)."""
    return [planet_classes.condition_words(condition, measure) for condition in cls.when]
