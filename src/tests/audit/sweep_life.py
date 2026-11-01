# tests/audit/sweep_life.py

"""
The class-blind life test of TEST.452 (Boss, 2:46 PM EDT 2026-10-29: life depends only on the conditions, never on the planet class
letter; GEN.915): a world that can hold life can hold it whatever its class letter is.

Two checks, both on the shared body sweep (`tests.audit.sweep`):

- `scores-class-blind` (a `Rule`): every terrestrial planet and moon of the sweep is given each of the class letters in turn, its stored
  conditions unchanged, and the life model (`physics.life_classes.score_classes`) must answer the same scores for every letter. A failure
  names the body (its conditions), the letter and the life class whose score moved. Because the sweep's own bodies reach life-capable
  conditions on only a few letters, this is what proves that the conditions that host life on one class host it on every class.
- `class_references()`: the life code (`LIFE_MODULES`) is scanned for the places that read the class letter or the class tables. Every
  reference must be listed with its reason in `tests/data/class_blind_life_notes.json` (a function that only writes the class's description
  text, say); one that is not listed is a place where the letter can change the life, and fails the test. A listed reference that is gone is
  reported as stale.

The second layer (`draw_problems`) is on the generator itself, with the class store of epoch 120 (GEN.889): every planet and moon of the sweep that has not been through the life step (`fresh_planets`)
has the life step of the generator (`generation.life.apply_life_data`) run on a copy of it. A body the life model supports (a class scores above 0) must be given life (unless its star is younger than any life needs or the first step of its schedule comes after the star's age, GEN.997 and GEN.1001), one
that no class supports must not, and the class that lives there must be one the model supports. Then the same conditions are put under every
other class letter and the life step is run again on the same stream: the class, the stage and the schedule must be identical. The odds are
those of the model (a supported world always gets life), so a rate that is not the model's, or that depends on the letter, fails and names the
letter and the conditions.

What the sweep cannot see: draws that run before the life model and key on the letter (since the cut-over of epoch 121: the surface composition and the label of the aerosol stream; GEN.915). They change the conditions, not the way the conditions are judged, so they belong to
the cut-over; the census counts of capable bodies per class (`capable_by_class`) show what they leave reachable.
"""

import ast
import importlib
import json
import pathlib

from tests.audit.sweep import Rule, broken_copy

NOTES_PATH = pathlib.Path(__file__).resolve().parent.parent / "data" / "class_blind_life_notes.json"

LIFE_MODULES = ("planetgen.generation.life", "planetgen.generation.evolution", "planetgen.physics.life_classes", "planetgen.physics.life_scores",
                "planetgen.physics.life_supply", "planetgen.physics.life_complexity", "planetgen.physics.life_timeline",
                "planetgen.physics.ms_lifetime", "planetgen.physics.habitability", "planetgen.physics.habitability_world",
                "planetgen.physics.biogenic_oxygen")
"""tuple: The modules that judge whether, and what, life lives on a body. A module that does not exist yet (epoch 120 adds some) is skipped."""

CLASS_NAMES = ("planet_class", "PLANET_CLASSES", "HABITABLE_PLANET_CLASSES", "classify", "life_ceiling_by_class")
"""tuple: The names whose use in the life code is a read of the class letter or of the class tables."""


def letters():
    from planetgen import tuning

    return tuple(tuning.PLANET_CLASSES)


def _scores(body):
    """`({class id: score}, id of the best class, {habitability field: value})` the life model and the habitability index give `body`
    (the index is recomputed from the stored conditions, as `habitability_world.rescore` does after the life is drawn)."""
    from planetgen.physics import habitability_world, life_classes

    scores = life_classes.score_classes(body)
    best = life_classes.best(scores)
    index = {}
    if getattr(body, "surface_temperature", None) is not None:
        index = habitability_world.compute(habitability_world.values_from_body(body), habitability_world.compact_host_name(getattr(body, "star", None)))
    return ({class_id: result.score for class_id, result in scores.items()} if scores else {}), (best.id if best else None), index


def _check_scores_class_blind(entry):
    body = entry.body
    if getattr(body, "body_type", None) != "t" or getattr(body, "planet_class", None) is None:
        return []                   # an Asteroid (GEN.1160) has no class letter to swap
    reference = _scores(body)
    problems = []
    for letter in letters():
        if letter == body.planet_class:
            continue
        scores, best, index = _scores(broken_copy(body, planet_class=letter))
        moved_index = sorted(name for name in set(index) | set(reference[2]) if index.get(name) != reference[2].get(name))
        if moved_index:
            problems.append(f"as class {letter} the habitability index differs from that as class {body.planet_class} in {', '.join(moved_index)}")
        if best != reference[1]:
            problems.append(f"as class {letter} the best life class is {best!r}, as class {body.planet_class} it is {reference[1]!r}")
        moved = sorted(class_id for class_id in set(scores) | set(reference[0]) if scores.get(class_id) != reference[0].get(class_id))
        if moved:
            problems.append(f"as class {letter} the scores of life class {', '.join(moved)} differ from those as class {body.planet_class}")
    return problems


RULES = (Rule("scores-class-blind", _check_scores_class_blind, lambda entry: None,
              note="no breaker on a body: the class letter is not an input of the life model, so it cannot be broken from outside; "
                   "the test breaks the model itself"),)


def capable_by_class(entries):
    """`{class letter: (bodies, bodies the life model gives a class)}` over the terrestrial entries, by the letter each was drawn with."""
    out = {}
    for entry in entries:
        body = entry.body
        if getattr(body, "body_type", None) != "t" or body.planet_class is None:
            continue
        total, capable = out.get(body.planet_class, (0, 0))
        out[body.planet_class] = (total + 1, capable + (1 if _scores(body)[1] else 0))
    return out


def load_notes(path=NOTES_PATH):
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)


def class_references():
    """`{"module:function": [line, ...]}` for every use of a name of `CLASS_NAMES` in the life modules (`<module>` for module level)."""
    found = {}
    for name in LIFE_MODULES:
        try:
            module = importlib.import_module(name)
        except ModuleNotFoundError:
            continue
        path = pathlib.Path(module.__file__)
        tree = ast.parse(path.read_text(encoding="utf-8"))
        scopes = {}

        def walk(node, scope):
            for child in ast.iter_child_nodes(node):
                if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                    walk(child, f"{scope}.{child.name}" if scope != "<module>" else child.name)
                    continue
                label = None
                if isinstance(child, ast.Attribute) and child.attr in CLASS_NAMES:
                    label = child.attr
                elif isinstance(child, ast.Name) and child.id in CLASS_NAMES:
                    label = child.id
                elif isinstance(child, ast.Constant) and isinstance(child.value, str) and child.value in CLASS_NAMES:
                    label = child.value
                if label:
                    scopes.setdefault(scope, []).append(child.lineno)
                walk(child, scope)

        walk(tree, "<module>")
        for scope, lines in scopes.items():
            found[f"{name}:{scope}"] = sorted(set(lines))
    return found


def unlisted_references(found, notes):
    """The references with no reason in the notes."""
    return {place: lines for place, lines in found.items() if place not in notes["allowed_references"]}


def stale_references(found, notes):
    """The listed places that no longer read the class letter."""
    return sorted(place for place in notes["allowed_references"] if place not in found)


LIFE_STREAM_SEED = 4242
"""int: The run seed the life step is given before each run, so two runs on the same conditions use the same stream."""


def fresh_planets():
    """`[(label, planet), ...]`: the planets and moons of the shared sweep (`tests.audit.sweep.entries`: `PLANETS_PER_ZONE` planets in each zone
    round each host star and their moons), none of which has been through the life step (a bare `Planet` has none; `StarSystem` runs it). Since
    the class cut-over (GEN.646, epoch 121) the class letter is derived from the conditions and cannot be forced on a `Planet`, so the sampler
    takes the bodies as the generator builds them and the other letters are put on copies by `draw_problems`."""
    from tests.audit import sweep

    return [(entry.label, entry.body) for entry in sweep.entries()]


def _supported(planet):
    from planetgen.physics import life_classes

    scores = life_classes.score_classes(planet)
    return {class_id: result.score for class_id, result in (scores or {}).items() if result.score > 0}


def life_outcome(planet, letter=None):
    """What the life step of the generator gives a copy of `planet` (under `letter` when given): `(class id, stage, schedule ages)`."""
    return _run_life(planet, letter)[0]


def _run_life(planet, letter=None):
    """`(life_outcome, waiting)`: `waiting` is True when the body went without life only because its star is too young for any (GEN.997) or the
    first step of its schedule comes after the star's age (GEN.1001), a world the model supports but on which life has not begun yet."""
    import copy
    from unittest import mock

    from planetgen import tuning
    from planetgen.generation import evolution, life
    from planetgen.util import draw

    clone = copy.deepcopy(planet)
    assert clone is not planet and getattr(clone, "star", None) is not None
    if letter is not None:
        clone.planet_class = letter
    begun = []
    real = life._life_has_begun

    def recorded(body):
        begun.append(real(body))
        return begun[-1]

    draw.set_run_seed(LIFE_STREAM_SEED)
    with mock.patch.object(life, "_life_has_begun", recorded):
        life.apply_life_data(clone)
    outcome = clone.life_class_id, clone.life_stage, tuple(getattr(clone, field, None) for field in evolution.SCHEDULE_FIELDS)
    age = getattr(clone.star, "age", None)
    waiting = outcome[0] is None and ((age is not None and age < tuning.LIFE_MIN_STAR_AGE_GY) or begun == [False])
    return outcome, waiting


def draw_problems(planet, other_letters=None):
    """The problems of one planet: a supported world without life, life on an unsupported one, a class that lives where the model
    does not support it, and any difference in the outcome when the same conditions carry another class letter (`other_letters`, default all)."""
    if getattr(planet, "body_type", None) != "t" or planet.planet_class is None:
        return []
    supported = _supported(planet)
    letter = planet.planet_class
    outcome, waiting = _run_life(planet)
    problems = []
    if supported and outcome[0] is None and not waiting:
        problems.append(f"class {letter}: the model supports {', '.join(sorted(supported))} and the generator drew no life")
    if not supported and outcome[0] is not None:
        problems.append(f"class {letter}: no life class is supported and the generator drew {outcome[0]}")
    if outcome[0] is not None and outcome[0] not in supported:
        problems.append(f"class {letter}: the generator drew {outcome[0]}, which the model does not support ({', '.join(sorted(supported)) or 'none'})")
    if supported:
        for other in (other_letters if other_letters is not None else letters()):
            if other == letter:
                continue
            again = life_outcome(planet, other)
            if again != outcome:
                problems.append(f"the same conditions as class {other} give {again[0]!r}, stage {again[1]!r}, not {outcome[0]!r}, stage {outcome[1]!r} as class {letter}")
    return problems


def draw_rates(planets):
    """`{class letter: (terrestrial planets, supported by the model, given life by the generator, supported but waiting for the star's age)}`
    over `[(label, planet), ...]`."""
    out = {}
    for _, planet in planets:
        if getattr(planet, "body_type", None) != "t" or planet.planet_class is None:
            continue
        total, supported, drawn, waiting = out.get(planet.planet_class, (0, 0, 0, 0))
        has = bool(_supported(planet))
        outcome, young = _run_life(planet) if has else (None, False)
        out[planet.planet_class] = (total + 1, supported + has, drawn + (1 if has and outcome[0] is not None else 0), waiting + (1 if has and young else 0))
    return out
