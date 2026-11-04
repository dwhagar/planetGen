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
* **Class match:** the geometric mean of the memberships of the class's conditions, so one condition badly missed pulls the match
  down and a condition missed outright makes it 0.
* **Best class:** the highest match; a tie goes to the class with more conditions (the narrower description), then to the letter.

Nothing here decides a stored letter yet: the first-match tree still does (GEN.1207 replaces it). `fit` and `best_class` are used
by the shadow scorer `src/tests/class_fit_cli.py` (GEN.1206).

Use:

    best_class(body)         # ("M", 0.97, "O", 0.41): letter, match, runner-up letter, runner-up match
    scores(body)             # {"M": 0.97, "O": 0.41, ...} for the classes whose gates the body passes
"""

import json
import math
import os
from collections.abc import Mapping

from planetgen.physics import elements, planet_classes

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
    if value is None:
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


def match(cls, found):
    """The percent match of `cls` for the named quantities `found`, as a fraction (the geometric mean of its memberships)."""
    total = 0.0
    for condition in cls.when:
        m = membership(condition, found)
        if m <= 0.0:
            return 0.0
        total += math.log(m)
    return math.exp(total / len(cls.when))


def passes_gates(cls, found):
    """Whether the body is of the class's kind and, for a terrestrial class, on the right side of the air gate."""
    if cls.kind != found["body_type"]:
        return False
    if cls.air is None:
        return True
    return (found["pressure_pa"] >= GATE_PRESSURE_PA) == cls.air


def _found(body):
    if isinstance(body, Mapping) and "body_type" in body and "mass_earth" in body:
        return body
    values = body if isinstance(body, Mapping) else planet_classes.values_from_body(body)
    return planet_classes.quantities(values)


def class_scores(body):
    """Every class whose gates the body passes, as a list of `(match, FitClass)` best first (ties: more conditions, then letter)."""
    found = _found(body)
    scored = [(match(cls, found), cls) for cls in CLASSES if passes_gates(cls, found)]
    scored.sort(key=lambda pair: (-pair[0], -len(pair[1].when), pair[1].letter))
    return scored


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
