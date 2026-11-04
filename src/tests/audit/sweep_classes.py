# tests/audit/sweep_classes.py

"""
The rule set of TEST.301 on the body sweep (`tests.audit.sweep`): the classes of planets and moons come from their conditions
(`physics/planet_classes.py`, GEN.645, GEN.647), and a request for a habitable world is answered by the conditions.

Rules (names as `RULES` lists them):

- `class-assigned`: every planet and moon has a class letter, the letter is one of the rule file's, and the letter's kind (terrestrial
  or giant) is the body's.
- `conditions-total`: `classify` returns a letter for the body (no `ClassError`), the same twice, the same from the body and from its
  stored columns (`values_from_body`), and it does not change the body.
- `class-from-conditions`: the stored class letter is the letter the conditions give, the first row of the body's kind whose conditions
  all hold (re-derived here from `quantities` and the rule rows, not by `classify`). Until the class leaves generation (GEN.646) the
  generator draws the letter first, so this rule is a known anomaly (`tests/data/class_sweep_expectations.json`).

`habitable_request_statuses` generates seeded systems that ask for a habitable world and returns how each answered
(`StarSystem.habitable_request_status`: `aim`, `fallback` or `unmet`); a request is met when the answer is `aim` or `fallback`.
"""

import copy

from planetgen.physics import planet_classes
from tests.audit.sweep import Rule, broken_copy, seed_for

REQUEST_SYSTEMS = 12
"""int: Systems generated with `+habitable_world`."""


def _conditions_letter(body):
    """The first row of the body's kind whose conditions hold, from the quantities of its stored values."""
    values = planet_classes.values_from_body(body)
    found = planet_classes.quantities(values)
    for cls in planet_classes.CLASSES:
        if cls.kind == found["body_type"] and planet_classes.holds(cls.when, found, owner=f"class {cls.letter}"):
            return cls.letter
    return None


def _asteroid(body):
    """Whether the conditions make the body an Asteroid (a terrestrial body under a roundness ratio of 1, GEN.1160): it is not a planet class and has no letter."""
    return planet_classes.classify_values(planet_classes.values_from_body(body)) is planet_classes.ASTEROID


def _check_assigned(entry):
    body = entry.body
    letter = getattr(body, "planet_class", None)
    if not letter:
        return [] if _asteroid(body) else ["no class"]
    if _asteroid(body):
        return [f"class {letter!r} on an Asteroid, which has no letter"]
    if letter not in planet_classes.BY_LETTER:
        return [f"class {letter!r} is not one of the rule file's"]
    # GEN.1142: a letter can have a row of each kind (class T is the giant "gas dwarf" and, for a rocky body whose air is mostly hydrogen, a
    # terrestrial one), so the letter fits the body when any of its rows is of the body's kind.
    kinds = sorted({cls.kind for cls in planet_classes.CLASSES if cls.letter == letter})
    if body.body_type not in kinds:
        return [f"class {letter!r} is a {'/'.join(kinds)!r} class and the body is {body.body_type!r}"]
    return []


def _break_assigned(entry):
    return broken_copy(entry.body, planet_class=None)


def _check_total(entry):
    body = entry.body
    before = copy.copy(vars(body))
    try:
        first = planet_classes.classify(body)
        again = planet_classes.classify(body)
        mapping = planet_classes.classify(planet_classes.values_from_body(body))
    except planet_classes.ClassError as error:
        return [f"cannot be classified: {error}"]
    problems = []
    if _asteroid(body):
        return [] if first is None and again is None and mapping is None else [f"an Asteroid gave {first!r}, then {again!r}, and {mapping!r} from the stored columns"]
    if first not in planet_classes.BY_LETTER:
        problems.append(f"classify gave {first!r}, not a letter of the rule file")
    if not first == again == mapping:
        problems.append(f"classify gave {first!r}, then {again!r}, and {mapping!r} from the stored columns")
    if vars(body) != before:
        problems.append("classify changed the body")
    return problems


def _break_total(entry):
    return broken_copy(entry.body, body_type="x")


def _check_from_conditions(entry):
    body = entry.body
    if _asteroid(body):
        return []                   # no letter to compare; `class-assigned` checks that it has none
    try:
        letter = _conditions_letter(body)
    except planet_classes.ClassError as error:
        return [f"cannot be classified: {error}"]
    stored = getattr(body, "planet_class", None)
    if letter != stored:
        return [f"stored class {stored!r}, conditions give {letter!r}"]
    return []


def _break_from_conditions(entry):
    try:
        letter = _conditions_letter(entry.body)
    except planet_classes.ClassError:
        return None
    other = "V" if letter != "V" else "A"
    return broken_copy(entry.body, planet_class=other)


RULES = (
    Rule("class-assigned", _check_assigned, _break_assigned),
    Rule("conditions-total", _check_total, _break_total),
    Rule("class-from-conditions", _check_from_conditions, _break_from_conditions),
)


def habitable_request_statuses(count=REQUEST_SYSTEMS):
    """`[status, ...]` for `count` seeded systems generated with `+habitable_world`."""
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from planetgen.util import draw

    statuses = []
    for index in range(count):
        draw.set_run_seed(seed_for("habitable-request", index))
        config = SystemConfig()
        config.HABITABLE_WORLD = True
        statuses.append(StarSystem(system_config=config).habitable_request_status())
    return statuses


def request_met_share(statuses):
    """The share of requests answered `aim` or `fallback`."""
    return sum(1 for status in statuses if status in ("aim", "fallback")) / len(statuses) if statuses else 0.0


def systems_with_a_class_from_conditions(count=2):
    """`True` where the same seeded systems, generated twice, hold the same classes (the class is determined by the draws)."""
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from planetgen.util import draw

    def classes(index):
        draw.set_run_seed(seed_for("class-repeat", index))
        system = StarSystem(system_config=SystemConfig())
        return [(body.planet_class, planet_classes.classify(body)) for planet in system.all_planets if planet.body_type != "a"
                for body in (planet, *planet.moons)]

    return all(classes(index) == classes(index) for index in range(count))
