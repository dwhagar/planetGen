# planetgen/physics/flavor_tree.py

"""
GEN.1430: the flavor text as a gated tree (`physics/data/flavor_tree.json`).

One tree per scope (`system`, `planet`, which also serves moons and rogue planets). Each node has a gate, a list of tests on stored quantities joined by "and"
(`quantity` with `min`, `lt`, `equals` or `not_equals`, the form `planet_classes.holds` reads for `planet_class_fit.json` and `oddities.json`), and a child inherits the gates
above it. A node's lines are offered only when every gate on the path passes, so a child never repeats what its parent guarantees. The walk gathers every line from every node
it reaches into one option list, takes out the ids already used in the system, and draws one by weight on the body's own sub-stream (`draw.Stream`, seeded from what is
final about the body, so the line does not depend on the order the system's bodies were made in); an empty list is no line, with no retry.

The gates read stored conditions only, never the class letter (the class is decided last, GEN.1397): `quantities_of` offers the quantity table of `oddities` and the four
flavor quantities of `QUANTITIES`. A node tagged `ancient_works` (or another family tag) is skipped when that tag is in the `off` set the galaxy plan hands the walk.

Line ids are permanent (a removed line stays with `retired` true, so a stored id still reads). The ids of the lines moved in from the old lists are the keys they were stored under
(`habitable:N`, `planet:N`, `orbital:N`, `system:N`), so no stored row changes.

Use:

    option = draw_line("planet", body, used_ids, off=())      # `Option(id, text, weight)` or None
    text_of("planet:3")                                       # the line a stored id names
"""

import json
import os
from dataclasses import dataclass

from planetgen.physics import elements, oddities, planet_classes, planet_type, surface_liquid_world
from planetgen.physics.elements import DataError
from planetgen.util import draw

DATA_FILE = os.path.join(elements.DATA_DIRECTORY, "flavor_tree.json")
"""str: The data file."""

SCOPES = ("system", "planet")
"""tuple: The scopes the file holds a tree for."""

QUANTITIES = ("has_ground", "life_present", "life_multicellular", "has_air")
"""tuple: The quantities this module adds to the table of `oddities.quantities_of`; each is 1 or 0."""

GATE_KEYS = {"quantity", "min", "lt", "equals", "not_equals"}


class FlavorError(DataError):
    """A flavor tree that is refused."""


@dataclass(frozen=True)
class Line:
    """One line of the tree."""
    id: str
    text: str
    weight: float
    retired: bool
    node: str
    tags: tuple


@dataclass(frozen=True)
class Option:
    """A line the walk offered, and the weight it carries after the multipliers on its path."""
    id: str
    text: str
    weight: float


def _check_gate(gate, owner):
    for condition in gate or []:
        if not isinstance(condition, dict) or "quantity" not in condition or set(condition) - GATE_KEYS or len(condition) < 2:
            raise FlavorError(f"{owner}: gate {condition!r} is not a quantity with min, lt, equals or not_equals")


def _read(node, scope, lines, nodes, seen_ids):
    node_id = node.get("id")
    if not isinstance(node_id, str) or not node_id or node_id in nodes:
        raise FlavorError(f"{scope}: node {node_id!r} needs an id of its own")
    _check_gate(node.get("gate"), node_id)
    multiplier = node.get("weight_mult", 1.0)
    if not isinstance(multiplier, (int, float)) or multiplier <= 0:
        raise FlavorError(f"{node_id}: weight_mult must be positive")
    tags = tuple(node.get("tags") or ())
    nodes[node_id] = node
    for raw in node.get("lines") or []:
        line_id, text, weight = raw.get("id"), raw.get("text"), raw.get("weight", 1.0)
        if not isinstance(line_id, str) or not line_id or line_id in seen_ids:
            raise FlavorError(f"{node_id}: line id {line_id!r} is missing or used twice")
        if not isinstance(text, str) or not text.strip() or not isinstance(weight, (int, float)) or weight <= 0:
            raise FlavorError(f"{line_id}: a line needs text and a positive weight")
        seen_ids.add(line_id)
        lines[line_id] = Line(line_id, text, float(weight), bool(raw.get("retired")), node_id, tags)
    for child in node.get("children") or []:
        _read(child, scope, lines, nodes, seen_ids)


def load(path=DATA_FILE):
    """Reads and checks the file: `(trees, lines)`, the tree of each scope as the file holds it and every line by id; a malformed file is refused (`FlavorError`)."""
    try:
        with open(path, encoding="utf-8") as handle:
            document = json.load(handle)
    except (OSError, ValueError) as error:
        raise FlavorError(f"cannot read the flavor tree {path}: {error}") from error
    trees, lines, nodes, seen = {}, {}, {}, set()
    for scope in SCOPES:
        tree = (document.get("scopes") or {}).get(scope)
        if not tree:
            raise FlavorError(f"the flavor tree has no {scope} scope")
        trees[scope] = tree
        _read(tree, scope, lines, nodes, seen)
    return trees, lines


_cache = {}


def _tree():
    if "tree" not in _cache:
        _cache["tree"] = load()
    return _cache["tree"]


def quantities_of(body):
    """The named quantities the gates read for `body`: the table of `oddities.quantities_of` and `QUANTITIES`. A quantity the body lacks is left out and a gate on it fails."""
    found = oddities.quantities_of(body)
    ground = planet_type.has_ground(body) and getattr(body, "surface_liquid", None) not in surface_liquid_world.LAVA_IDS
    found["has_ground"] = 1 if ground else 0
    life = getattr(body, "life_class", None) is not None
    found["life_present"] = 1 if life else 0
    found["life_multicellular"] = 1 if life and getattr(body, "life_stage", None) in ("multicellularity", "technological_civilization") else 0
    found["has_air"] = 1 if (getattr(body, "atmospheric_pressure", None) or 0.0) > 0.0 else 0
    return found


def walk(scope, body, off=()):
    """Every line the tree offers `body` in `scope`, as `Option`s in file order (retired lines are never offered)."""
    trees, lines = _tree()
    found = quantities_of(body) if body is not None else {}
    options = []

    def visit(node, multiplier):
        if set(node.get("tags") or ()) & set(off):
            return
        if not planet_classes.holds(node.get("gate") or [], found, node.get("id", ""), missing=None):
            return
        multiplier *= node.get("weight_mult", 1.0)
        for raw in node.get("lines") or []:
            line = lines[raw["id"]]
            if not line.retired:
                options.append(Option(line.id, line.text, line.weight * multiplier))
        for child in node.get("children") or []:
            visit(child, multiplier)

    visit(trees[scope], 1.0)
    return options


def seed_of(scope, body, star=None):
    """The seed of the body's own flavor stream: what is final about the body and its star, nothing about the order the system made its bodies in."""
    star = star if star is not None else getattr(body, "star", None)
    return (f"flavor:{scope}:{getattr(body, 'is_moon', False)}:{getattr(body, 'mass', None)!r}:{getattr(body, 'radius', None)!r}:"
            f"{getattr(body, 'distance', None)!r}:{getattr(star, 'age', None)!r}:{getattr(star, 'mass', None)!r}")


def draw_line(scope, body, used_ids=(), off=(), star=None):
    """
    The `Option` drawn for `body` in `scope`, or `None`: the walk's options less the ids in `used_ids`, one drawn by weight on the body's own stream (no retry on an empty list).
    """
    options = [option for option in walk(scope, body, off) if option.id not in set(used_ids)]
    if not options:
        return None
    with draw.bound(draw.Stream(seed_of(scope, body, star))):
        return draw.choices(options, weights=[option.weight for option in options], k=1)[0]


def text_of(line_id):
    """The text of the line a stored id names (retired lines still read), or `None` for no id."""
    if line_id is None:
        return None
    return _tree()[1][line_id].text
