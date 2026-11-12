# planetgen/physics/class_fit.py

"""
Planet and moon classes as percent matches (GEN.1205, GEN.1206; design: reports/class-m-and-trek-classes-2026-10-31.md sections 4 to 6)
========================================================================================================================================

`planet_classes.classify` walks the rows of `planet_classes.json` in order and returns the first whose conditions hold, so a row
quietly relies on the rows before it having said no. This module is the other way to name a class: every class is described by its
own positive conditions (`physics/data/planet_class_fit.json`), each condition is graded, and a body's class is the one it fits best.

* **Gates (kinds, not degrees):** giant versus terrestrial (`body_type`) and, for terrestrial classes, air versus no air (pressure
  at or above `gate_pressure_pa`, 10 Pa). A class is only scored against bodies that pass its gates.
* **Condition membership:** 1 inside the edge; across the shoulder beyond it a straight line down to 0; 0 past that. The shoulder of an edge is its width in the fudge table
  (`physics/data/fudge_table.json`, GEN.1313: how far the quantity could move if the stand-ins that make it were wrong by what the registry says), measured in the quantity's own
  unit, or as the natural log of the ratio for a pressure, a gravity or a mass; a body `z` widths past an edge has the membership `1 - z`. `equals` is 0 or 1. A
  quantity the body does not store scores 0.
* **Domains (GEN.1221):** each quantity belongs to one domain of the class-by-domain matrix (`physics/class_domains.py`); `domain_scores` grades a body in each.
* **Fudge (GEN.1314, plan section 4c):** we only fudge when we must and by as little as possible. A body that fits a class on every condition takes it and nothing is fudged. Otherwise a class
  needs some **values** moved to fit (the mean temperature with the regional, day and night temperatures counts as one value), each by some number of widths; the body takes the class that
  needs the fewest values moved, then the smallest total move, then the one with more conditions (the narrower description), then the letter. A class that needs a value moved past a cap is
  not eligible, and the caps come in tiers (`CAP_TIERS`): the smaller tier is tried first and the next opens only for a body the smaller cannot place.
* **Class match:** 1 for a full fit, and for a fudged class `1 - total / (values * largest cap)`, so it is above 0 for every eligible fit; 0 for a body no class is eligible for.
* **Best class:** the one the fudge rule takes.

`classify` (GEN.1207) is the one that decides a stored class: the best class, its match, the runner-up and its match, and the digit of the class code.
A body that no class is eligible for (every class needs a value moved past the largest cap) is a hole in the class set and is settled in two steps, both stand-ins (GEN.725, section 6.189):
the conditions of the gravity and the air-amount domains are set aside (`LIMIT_DOMAINS`) and the class of what is left among the families (`LIMITED_FAMILIES`) is the letter, with the digit 9
(limited by another domain, `class_code`); a body that still fits no class takes the class it is nearest to (the fewest values and the smallest move with no cap, among the ordinary worlds
`NEAREST_FAMILIES`), with a stored match of 0. The class a body takes keeps which values were fudged and by how many widths (`Classification.fudge_values`, `fudge_amount`). `fit` and `best_class` are also used by the shadow scorer `src/tests/class_fit_cli.py` (GEN.1206).

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
FUDGE_FILE = os.path.join(elements.DATA_DIRECTORY, "fudge_table.json")


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
    with open(FUDGE_FILE, encoding="utf-8") as handle:
        fudge = json.load(handle)
    classes = tuple(FitClass(row) for row in document["classes"])
    widths = {quantity: {"unit": entry["unit"], "edges": {row["edge"]: row["width"] for row in entry["edges"]}} for quantity, entry in fudge["quantities"].items()}
    for cls in classes:
        if not cls.when:
            raise FitError(f"class {cls.letter} ({cls.kind}) has no conditions of its own")
        if (cls.kind == planet_classes.GIANT) != (cls.air is None):
            raise FitError(f"class {cls.letter} ({cls.kind}): the air gate is null for giants and true or false for terrestrial bodies")
        for condition in cls.when:
            keys = set(condition) - {"quantity", "min", "lt", "equals"}
            if keys or len(condition) != 2:
                raise FitError(f"class {cls.letter}: condition {condition!r} is not one quantity with one of min, lt or equals")
            if "equals" not in condition:
                edge = condition.get("min", condition.get("lt"))
                if edge not in widths.get(condition["quantity"], {"edges": {}})["edges"]:
                    raise FitError(f"class {cls.letter}: the fudge table has no width for {condition['quantity']!r} at {edge}; run python scripts/build_fudge_table.py")
            try:
                class_domains.domain_of(condition["quantity"])
            except KeyError:
                raise FitError(f"class {cls.letter}: {condition['quantity']!r} belongs to no domain of the class-by-domain matrix") from None
    if len({(c.letter, c.kind) for c in classes}) != len(classes):
        raise FitError("a class letter is used twice for one kind")
    return document, classes, widths


DOCUMENT, CLASSES, WIDTHS = _load()
"""tuple: The `FitClass` rows; `DOCUMENT` the whole file; `WIDTHS` the fudge table's width at each edge, `{quantity: {"unit": "ln" or "plain", "edges": {edge: width}}}`."""

GATE_PRESSURE_PA = DOCUMENT["gate_pressure_pa"]
LETTERS = tuple(sorted({c.letter for c in CLASSES}))

THERMAL_QUANTITIES = ("temperature_k", "hottest_region_k", "coldest_region_k", "day_temperature_k", "night_temperature_k")
"""tuple: The quantities that move together when the mean surface temperature moves; they count as one value in the fudge rule."""
THERMAL_KEY = "temperature_k"
"""str: The name the thermal value goes by in `Classification.fudge_values`."""
CAP_TIERS = (0.25, 1.0)
"""tuple: The caps, in widths, a value may be moved by (registry CP1): the sweep of 1,624 bodies (reports/fudge-cap-sweep-2026-11-06.md) closed its near misses at 0.25 width and its remaining
fits by 1; the smaller tier is tried first, and the largest is also the scale of the match."""


def overshoot(condition, found):
    """
    How far past the edge of one class condition the quantities `found` are, in widths of the edge in the fudge table: 0 inside, infinity for a quantity the body does not store, an `equals`
    that differs, a width of 0 or a log of a quantity that is not above 0. A pressure, a gravity or a mass is measured as the natural log of the ratio to the edge.
    """
    quantity = condition["quantity"]
    value = found.get(quantity)
    if value is None or value != value:     # a missing or NaN quantity is no match: every comparison with NaN is false, which would otherwise read as a pass (GEN.1593)
        return math.inf
    if "equals" in condition:
        return 0.0 if value == condition["equals"] else math.inf
    below = "min" in condition
    edge = condition["min"] if below else condition["lt"]
    if (value >= edge) if below else (value < edge):
        return 0.0
    entry = WIDTHS[quantity]
    width = entry["edges"][edge]
    if width <= 0.0:
        return math.inf
    if entry["unit"] == "ln":
        if value <= 0.0 or edge <= 0.0:
            return math.inf
        return abs(math.log(edge / value)) / width
    return abs(edge - value) / width


def membership(condition, found):
    """The membership (0 to 1) of the named quantities `found` in one condition of a class: 1 inside the edge, `1 - z` for a body `z` widths past it (`overshoot`), 0 from one width."""
    return max(0.0, 1.0 - overshoot(condition, found))


Option = namedtuple("Option", ("match", "cls", "moves", "tier"))
"""namedtuple: One class a body could take: its `match`, the class, the values it would have to move (`{value: widths}`, the thermal ones as one under `THERMAL_KEY`) and the first tier of
`CAP_TIERS` that holds every move (`None` when a move is past the largest cap)."""


def moves_of(cls, found, relax=()):
    """`{value: widths}` the body would have to move to fit `cls`: the conditions in the domains `relax` are left out; the thermal quantities count as one value (the largest of their moves)."""
    moves = {}
    for condition in cls.when:
        if relax and class_domains.domain_of(condition["quantity"]) in relax:
            continue
        z = overshoot(condition, found)
        if z > 0.0:
            key = THERMAL_KEY if condition["quantity"] in THERMAL_QUANTITIES else condition["quantity"]
            moves[key] = max(moves.get(key, 0.0), z)
    return moves


def tier_of(moves):
    """The index of the first tier of `CAP_TIERS` that holds every move, or `None`."""
    largest = max(moves.values(), default=0.0)
    for index, cap in enumerate(CAP_TIERS):
        if largest <= cap and (index < len(CAP_TIERS) - 1 or largest < cap):        # the largest cap is the width at which the membership reaches 0: a move of it is no fit
            return index
    return None


def match_of(moves):
    """The match of a class that needs `moves`: 1 for none, `1 - total / (values * largest cap)` for moves within the caps, 0 past them."""
    if not moves:
        return 1.0
    if tier_of(moves) is None:
        return 0.0
    return 1.0 - sum(moves.values()) / (len(moves) * CAP_TIERS[-1])


def _order(option):
    moves = option.moves
    return (len(CAP_TIERS) if option.tier is None else option.tier, len(moves), round(sum(moves.values()), 9), -len(option.cls.when), option.cls.letter)


def match(cls, found, relax=()):
    """The match of `cls` for the named quantities `found` (`match_of` of the moves it needs); the conditions in the domains `relax` are left out."""
    return match_of(moves_of(cls, found, relax))


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


def _options(found, relax=(), only=None):
    """The `Option` of every class whose gates the body passes (among the letters `only`, if given), in no order."""
    options = []
    for cls in CLASSES:
        if passes_gates(cls, found) and (only is None or cls.letter in only):
            moves = moves_of(cls, found, relax)
            options.append(Option(match_of(moves), cls, moves, tier_of(moves)))
    return options


def _ranked(found, relax=(), only=None):
    """The options best first by the fudge rule: the smaller tier, then the fewest values moved, then the smallest total move, then more conditions, then the letter. Options past the
    largest cap come last, in the same order."""
    return sorted(_options(found, relax, only), key=_order)


def class_scores(body):
    """Every class whose gates the body passes, as a list of `(match, FitClass)` best first by the fudge rule."""
    return [(option.match, option.cls) for option in _ranked(_found(body))]


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
    if len(ranked) == 1 or ranked[1][0] <= 0.0:
        return best.letter, best_m, None, 0.0
    return best.letter, best_m, ranked[1][1].letter, ranked[1][0]


LIMIT_DOMAINS = ("gravity", "air_amount")
"""tuple: The domains a hole body is judged without (GEN.1207, registry HC1): a body whose class is found only with its gravity and its air amount set aside is limited by
another domain, digit 9 of the class code."""

NEAREST_FAMILIES = ("M", "R", "P", "O", "H", "F", "G")
"""tuple: The classes a body that fits none takes its nearest from (registry HC2): the ordinary worlds with air. A class with a narrow window (Y, Z, U, X) or one defined by an edge
(B, V, K, N, L, E) is never the nearest of a body that missed all of them."""

LIMITED_FAMILIES = ("M", "R", "P", "O", "H", "F", "G", "Y", "U", "X")
"""tuple: The classes a body can be limited in (`class_code.LETTERS`, less Z): a class defined by its gravity or its air amount (B, V, K, N, Z, L, E) has those conditions as its
own edges, so setting them aside would match it to everything."""

LIMITED_DIGIT = 9
"""int: The digit of the class code of a body whose class is found with its gravity and air amount set aside (`class_code.POSITION`, "limited by another domain")."""

Classification = namedtuple("Classification", ("letter", "kind", "match", "runner_up", "runner_up_match", "digit", "fudge_values", "fudge_amount"))
"""namedtuple: The class of a body: its letter and kind, the match (0 when no class fits it), the runner-up letter (`None` when no other class matches at all) and its
match, the digit of the class code (0, or 9 when the body is limited by another domain), the values the class needed moved (a sorted tuple of quantity names, the thermal ones as
`THERMAL_KEY`; empty for a full fit) and the total move in widths (0.0 for a full fit; `None` when a value the class reads is not stored)."""


def _nearest(found):
    """The class a hole body is nearest to: the fewest values moved and then the smallest move with no cap, among the ordinary worlds (or any class that passes the gates), as an `Option`."""
    for only in (NEAREST_FAMILIES, None):
        options = _options(found, only=only)
        if options:
            return sorted(options, key=lambda o: (len(o.moves), round(sum(o.moves.values()), 9), -len(o.cls.when), o.cls.letter))
    return []


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
        raise planet_classes.ClassError(f"the body stores a quantity that is not a number ({', '.join(not_a_number)}); a class is never fitted to it")
    if found["body_type"] == planet_classes.TERRESTRIAL and found["pressure_pa"] >= GATE_PRESSURE_PA and found["temperature_k"] is None:
        raise planet_classes.ClassError("the classes read the mean surface temperature of a world with air, which the body does not store")
    ranked = _ranked(found)
    if not ranked:
        raise FitError("no class passes the gates of the body")
    digit = 0
    if ranked[0].tier is None:
        relaxed = _ranked(found, relax=LIMIT_DOMAINS, only=LIMITED_FAMILIES)
        if relaxed and relaxed[0].tier is not None:
            ranked, digit = relaxed, LIMITED_DIGIT
        else:
            ranked = [option._replace(match=0.0, tier=None) for option in _nearest(found)]
    best = ranked[0]
    runner = ranked[1] if len(ranked) > 1 and ranked[1].match > 0.0 else None       # a class that matches 0 percent is no runner-up
    total = sum(best.moves.values())
    return Classification(best.cls.letter, best.cls.kind, best.match, None if runner is None else runner.cls.letter, 0.0 if runner is None else runner.match, digit,
                          tuple(sorted(best.moves)), None if math.isinf(total) else total)


def rule_words(cls):
    """The conditions of a `FitClass` as a list of plain phrases (`planet_classes.condition_words`)."""
    return [planet_classes.condition_words(condition) for condition in cls.when]
