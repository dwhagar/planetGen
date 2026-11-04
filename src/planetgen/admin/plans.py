# planetgen/admin/plans.py

"""
The change plan engine (ADM.58, the library under the 8.6 headline "edit everything", ADM.57;
docs/design/interface/edit-everything.md)
==============================================================================================

An admin change, small or large, is a **plan**: an ordered list of operations on stored objects. Nothing is written
by this module and it reads no database. It works on loaded objects (`World`: the `StarSystem`s, where they sit in
their sector, and the edit-version stamp of every object), so it is tested on its own; the store (ADM.59), the apply
job (ADM.60) and the screens come later and use what is here.

Operations (`OPERATIONS`), each aimed at an object by a `Ref` (system, star, planet, moon, belt):

* `Move`: a planet's or belt's distance from its star, a moon's distance from its planet, or a system's position in
  its sector.
* `SetAspect`: one aspect of an object (`ASPECTS`: the name, a single star's type; a body's class is derived from its conditions, so it is not set).
  The aspect catalogue (ADM.62) widens the list; `register_aspect` is where an aspect joins it.
* `Regenerate`: redraw a planet, moon or belt, with a keep set (`KEEPABLE`: the name, the orbit).
* `Add`: a planet, a moon or a belt. It gets a provisional id (`new-1`, ...) that later operations of the same plan
  may name.
* `Remove`: a planet (with its moons), a moon or a belt.

`simulate` applies the operations, in order, to copies of the systems with the edit functions of `admin/edits.py`,
which stabilise a system after every change (they move bodies and never remove one), so the effect is the whole
effect. Every draw comes from a stream seeded by the plan's seed and the operation's place in the plan, so the same
plan on the same stored state gives the same effect. `compute_effect` diffs the copies against the originals: one
`Change` per field that differs, `indirect` for an object no operation aimed at, with a plain sentence per
operation and the layout problems that are new (`Effect.warnings`). An operation that cannot be applied raises
`PlanError` naming why; a layout problem the stabilising could not remove is a warning for the approval screen.

A plan's state is one of `STATES`, moved only along `TRANSITIONS`; every move records who made it and in which role
(`Transition`). An admin's own plan goes draft, proposed, approved with no second approver (Boss, 2026-10-16, 1:49
PM EDT); `requested` is a plan someone below the approving role submitted (10.0). A proposed or approved plan is
**stale** when an object it touches has a different edit-version stamp than when it was proposed (`stale_keys`):
only an edit moves a stamp, never the nightly orbital update. After an orbital update a pending move shows its
**drift** (`Plan.drift`: the distance it was planned from and the one the object has now) and `Plan.rebase` plans
from the new values again; an approved plan that drifted is sent back to proposed, so the admin approves what is
shown.

A proposed plan also keeps the systems as they will be (`Plan.outcomes`, each in the `StarSystem.to_dict()` envelope
the remote upload uses, API.13), so the apply step writes exactly what the admin was shown.
"""

import copy
import dataclasses
import math
import secrets
import time
from dataclasses import dataclass, field
from typing import ClassVar

from planetgen import tuning
from planetgen.admin import edits
from planetgen.galaxy import sector as galaxy_sector
from planetgen.generation import system as system_module, validation
from planetgen.generation.belt import AsteroidBelt
from planetgen.generation.planet import Planet
from planetgen.names.bodies import moon_letters
from planetgen.physics import magnetism, planets as planetPhysics, radiation, unit_registry as units
from planetgen.util import detmath, draw

PLAN_VERSION = 1
"""int: The version of a plan's serialised form."""

ENVELOPE_VERSION = system_module.SERIALIZATION_SCHEMA_VERSION
"""int: The version of the system envelope (`StarSystem.to_dict()`) an operation and a plan outcome use."""

_km_to_au = units.converter("km", "au")

# --- states -------------------------------------------------------------------------------------------------------

DRAFT, REQUESTED, PROPOSED, APPROVED, APPLIED, UNDONE, STALE, REJECTED, WITHDRAWN = (
    "draft", "requested", "proposed", "approved", "applied", "undone", "stale", "rejected", "withdrawn")
STATES = (DRAFT, REQUESTED, PROPOSED, APPROVED, APPLIED, UNDONE, STALE, REJECTED, WITHDRAWN)

TRANSITIONS = {
    DRAFT: (REQUESTED, PROPOSED, WITHDRAWN),
    REQUESTED: (PROPOSED, REJECTED, WITHDRAWN),
    PROPOSED: (APPROVED, REJECTED, WITHDRAWN, STALE),
    APPROVED: (APPLIED, PROPOSED, WITHDRAWN, STALE),        # back to proposed when an orbital update moved what it plans
    APPLIED: (UNDONE,),
    STALE: (PROPOSED, WITHDRAWN),
    UNDONE: (),
    REJECTED: (),
    WITHDRAWN: (),
}
"""dict: Each state and the states it may go to."""

ROLE_USER, ROLE_MODERATOR, ROLE_ADMIN, ROLE_SYSTEM = "user", "moderator", "admin", "system"
ROLES = (ROLE_USER, ROLE_MODERATOR, ROLE_ADMIN, ROLE_SYSTEM)
APPROVER_ROLES = (ROLE_ADMIN,)
"""tuple: The roles that propose, approve and reject a plan. `system` is the engine itself (a stale mark)."""

# --- objects ------------------------------------------------------------------------------------------------------

SYSTEM, STAR, PLANET, MOON, BELT, COMET = "system", "star", "planet", "moon", "belt", "comet"
KINDS = (SYSTEM, STAR, PLANET, MOON, BELT, COMET)

KEEPABLE = ("name", "orbit")
"""tuple: The aspects a regeneration may be told to keep. The edit functions regenerate a body at its own orbit under
its own name, so those two are always kept for now; redrawing them comes with the rearranger (ADM.75)."""

DEFAULT_KEEP = ("name", "orbit")


class PlanError(ValueError):
    """An operation or a plan cannot be done; the message says why, in plain words."""


class PlanStateError(PlanError):
    """A plan was asked to move to a state it cannot reach from where it is, or by a role that may not."""


class PlanStale(PlanError):
    """The plan touches an object that was edited after the plan was proposed."""


@dataclass(frozen=True)
class Ref:
    """
    An object of a world: `system` (a system's id), `kind` and, for anything but the system itself, the object's row
    id (or its provisional `new-N` id while a plan adds it).
    """

    system: object
    kind: str = SYSTEM
    id: object = None

    def __post_init__(self):
        if self.kind not in KINDS:
            raise PlanError(f"unknown kind of object: {self.kind!r}")
        if (self.kind == SYSTEM) != (self.id is None):
            raise PlanError("a system has no id of its own; every other object needs one")

    @property
    def key(self):
        """str: One string for the object, the key of its edit-version stamp."""
        return str(self.system) if self.kind == SYSTEM else f"{self.system}/{self.kind}:{self.id}"

    def to_dict(self):
        return {"system": self.system, "kind": self.kind, "id": self.id}

    @classmethod
    def from_dict(cls, data):
        return cls(data["system"], data.get("kind", SYSTEM), data.get("id"))


@dataclass
class World:
    """
    What the engine reads: the loaded `systems` by id, where each sits (`placements`: id to `(x, y, z)` in light-years
    inside a sector cube of `edge_ly`), and the current edit-version `stamps` (object key to anything comparable: a
    journal row id, a version number).
    """

    systems: dict
    placements: dict = field(default_factory=dict)
    edge_ly: float = None
    stamps: dict = field(default_factory=dict)


# ---------------------------------------------------------------------
# Finding objects
# ---------------------------------------------------------------------

def _ident(obj, fallback):
    """The id a stored object goes by: its row id, else its place (`fallback`)."""
    row = getattr(obj, "db_id", None)
    return row if row is not None else fallback


def _host(system, owner):
    """The star whose planet list is `owner`."""
    for star, planets in validation.star_lists(system):
        if planets is owner:
            return star
    raise PlanError("that body is not in a star's planet list")


def _locate(sim, ref):
    """
    Finds `ref` in the simulated systems.

    Returns:
        tuple: `(system, obj, owner)`: `owner` is the list a planet or belt sits in, the planet of a moon, `None`
        for a system or a star.
    """
    system = sim.system(ref)
    if ref.kind == SYSTEM:
        return system, system, None
    if ref.kind == STAR:
        for index, star in enumerate(system.stars):
            if str(_ident(star, f"@{index}")) == str(ref.id):
                return system, star, None
        raise PlanError(f"no such star: {ref.id}")
    if ref.kind in (PLANET, MOON, BELT):
        try:
            body, owner = edits.find_body(system, ref.kind, ref.id)
        except edits.BodyNotFound as exc:
            raise PlanError(str(exc)) from exc
        return system, body, owner
    raise PlanError(f"{ref.kind} objects cannot be the target of an operation")


def _settle(system, summary):
    """Stabilises the changed system (`validation.stabilize_star_system`) and wraps up the result."""
    report = validation.stabilize_star_system(system)
    warnings = [f"{p.body}: {p.message}." for p in report.problems]
    return edits.EditResult(summary, report.moved, report.reclassified, report.removed, warnings)


def _letter_moons(planet):
    """Names a planet's moons `<planet name><letter>` in orbital order."""
    planet.moons.sort(key=lambda m: m.distance)
    for index, moon in enumerate(planet.moons):
        moon.name = f"{planet.name}{moon_letters(index)}"


def _number(value, label, positive=True):
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise PlanError(f"{label} must be a number, not {value!r}") from exc
    if not math.isfinite(number) or (positive and number <= 0):
        raise PlanError(f"{label} must be a positive number, not {value!r}")
    return number


# ---------------------------------------------------------------------
# Operations
# ---------------------------------------------------------------------

class Operation:
    """What every operation shares: its `OP` name, and the envelope it serialises in (`to_dict`, `from_dict`)."""

    OP: ClassVar[str] = ""

    def to_dict(self):
        data = {"op": self.OP, "schema_version": ENVELOPE_VERSION}
        for item in dataclasses.fields(self):
            value = getattr(self, item.name)
            if isinstance(value, Ref):
                value = value.to_dict()
            elif isinstance(value, (set, frozenset)):
                value = sorted(value)
            elif isinstance(value, tuple):
                value = list(value)
            data[item.name] = value
        return data

    @classmethod
    def from_dict(cls, data):
        """Rebuilds an operation of any kind from its envelope."""
        try:
            kind = OPERATIONS[data["op"]]
        except KeyError as exc:
            raise PlanError(f"unknown operation: {data.get('op')!r}") from exc
        if data.get("schema_version", ENVELOPE_VERSION) > ENVELOPE_VERSION:
            raise PlanError(f"operation envelope {data['schema_version']} is newer than this code understands "
                            f"(max {ENVELOPE_VERSION})")
        values = {}
        for item in dataclasses.fields(kind):
            if item.name not in data:
                continue
            value = data[item.name]
            if item.name == "target" or item.name == "owner":
                value = Ref.from_dict(value)
            elif item.name == "keep":
                value = tuple(value)
            elif item.name == "position" and value is not None:
                value = tuple(value)
            elif item.name == "before" and value is not None:
                value = {name: tuple(part) if isinstance(part, list) else part for name, part in value.items()}
            values[item.name] = value
        return kind(**values)


@dataclass
class Move(Operation):
    """
    Moves `target`: a planet or belt to `distance_au` from its star, a moon to `distance_au` from its planet, a
    system to `position` in its sector. `before` is what it was moved from, kept so the approval screen can show
    the drift an orbital update causes (`Plan.drift`).
    """

    OP: ClassVar[str] = "move"
    target: Ref
    distance_au: float = None
    position: tuple = None
    before: dict = None
    forced: bool = False


@dataclass
class SetAspect(Operation):
    """Sets one aspect of `target` to `value` (`ASPECTS`). `forced` takes a value the generator would not draw."""

    OP: ClassVar[str] = "set_aspect"
    target: Ref
    aspect: str = None
    value: object = None
    forced: bool = False


@dataclass
class Regenerate(Operation):
    """Redraws a planet, moon or belt, keeping what `keep` names (`KEEPABLE`)."""

    OP: ClassVar[str] = "regenerate"
    target: Ref
    keep: tuple = DEFAULT_KEEP
    forced: bool = False


@dataclass
class Add(Operation):
    """
    Adds a `planet` or a `belt` to the system `target`, or a `moon` to the planet `target`, at `distance_au`. A
    planet needs a `name`; a belt its `lower_au` and `upper_au`. A new body's class is derived from its conditions
    (epoch 121), not asked for. `new_id` is the provisional id the new object goes by in the plan.
    """

    OP: ClassVar[str] = "add"
    target: Ref
    what: str = PLANET
    distance_au: float = None
    name: str = None
    lower_au: float = None
    upper_au: float = None
    new_id: str = None
    forced: bool = False


@dataclass
class Remove(Operation):
    """Removes a planet (with its moons), a moon or a belt."""

    OP: ClassVar[str] = "remove"
    target: Ref
    forced: bool = False


OPERATIONS = {cls.OP: cls for cls in (Move, SetAspect, Regenerate, Add, Remove)}
"""dict: Each operation by its `OP` name."""


def _target_keys(operations):
    """The keys of the objects the operations are aimed at, and of the ones an `Add` makes."""
    keys = set()
    for op in operations:
        keys.add(op.target.key)
        if isinstance(op, Add) and op.new_id is not None:
            keys.add(Ref(op.target.system, op.what, op.new_id).key)
    return keys


# ---------------------------------------------------------------------
# Applying an operation to the copies
# ---------------------------------------------------------------------

class _Sim:
    """The systems of a world, copied the first time an operation touches one, and their placements."""

    def __init__(self, world):
        self.world = world
        self.systems = {}
        self.placements = dict(world.placements)

    def system(self, ref):
        if ref.system not in self.systems:
            if ref.system not in self.world.systems:
                raise PlanError(f"no such system: {ref.system}")
            self.systems[ref.system] = copy.deepcopy(self.world.systems[ref.system])
        return self.systems[ref.system]


def _outside_floor(system, star, body):
    floor = system_module.StarSystem._engulfment_radius_au(star)
    return floor and validation.inner_edge_au(body) < floor


def _move_system(sim, op):
    system = sim.system(op.target)
    if op.position is None or len(op.position) != 3:
        raise PlanError("moving a system needs a position of three coordinates")
    position = tuple(_number(value, "a coordinate", positive=False) for value in op.position)
    edge = sim.world.edge_ly
    if edge is not None and not all(0 <= value <= edge for value in position):
        raise PlanError(f"{position} is outside the sector (0 to {edge:g} ly on each axis)")
    for other_id, other_position in sim.placements.items():
        other = sim.world.systems.get(other_id)
        if other_id == op.target.system or other is None:
            continue
        gap = detmath.dist(position, other_position)
        needed = galaxy_sector.required_separation_ly(system, other)
        if gap < needed and not op.forced:
            raise PlanError(f"it would sit {gap:.3g} ly from system {other_id}, inside the {needed:.3g} ly their Hill "
                            "spheres need")
    sim.placements[op.target.system] = position
    return edits.EditResult(f"Moved system {op.target.system} to {position}.", [], [], [], [])


def _move_orbit(sim, op):
    system, body, owner = _locate(sim, op.target)
    distance = _number(op.distance_au, "the distance")
    if op.target.kind == MOON:
        low_km, high_km = planetPhysics.moon_orbit_bounds_km(owner)
        if not _km_to_au(low_km) <= distance <= _km_to_au(high_km) and not op.forced:
            raise PlanError(f"a moon of {owner.name} can sit {_km_to_au(low_km):.3g} to {_km_to_au(high_km):.3g} AU "
                            f"from it, not {distance:.3g}")
        body.distance = distance
        validation.refresh_moon_orbit(body, owner)
        owner.moons.sort(key=lambda moon: moon.distance)
        validation.refresh_planet_reflex(owner, force=True)
        return _settle(system, f"Moved {body.name} to {distance:.4g} AU from {owner.name}.")
    star = _host(system, owner)
    ceiling = validation.orbit_ceiling_au(star)
    if op.target.kind == BELT:
        shift = distance - body.distance
        body.distance, body.lower_limit, body.upper_limit = distance, body.lower_limit + shift, body.upper_limit + shift
    else:
        body.distance = distance
    if validation.reach_au(body) > ceiling * (1 + validation.RELATIVE_TOLERANCE):
        raise PlanError(f"{distance:.4g} AU is beyond the farthest stable orbit ({ceiling:.4g} AU)")
    if _outside_floor(system, star, body):
        raise PlanError(f"{distance:.4g} AU is inside the star")
    owner.sort(key=lambda item: item.distance)
    if op.target.kind == PLANET:
        validation.reconcile_moved_planet(body)
    label = body.name if getattr(body, "name", None) else "the asteroid belt"
    result = _settle(system, f"Moved {label} to {distance:.4g} AU.")
    if not math.isclose(body.distance, distance, rel_tol=1e-9):      # the spacing pass moves a body that crowds the one before it
        result = result._replace(summary=f"{result.summary[:-1]}, then to {body.distance:.4g} AU to keep the orbits spaced.")
    return result


def _apply_move(sim, op):
    if op.target.kind == SYSTEM:
        return _move_system(sim, op)
    if op.target.kind not in (PLANET, MOON, BELT):
        raise PlanError(f"a {op.target.kind} cannot be moved")
    return _move_orbit(sim, op)


@dataclass(frozen=True)
class Aspect:
    """One aspect an operation can set: the kinds it applies to and the function that sets it (`ASPECTS`)."""

    name: str
    kinds: tuple
    apply: object


def _set_name(sim, op, system, obj, owner):
    value = str(op.value or "").strip()
    if not value or len(value) > 120:
        raise PlanError("a name is 1 to 120 characters")
    old = getattr(obj, "name", None)
    obj.name = value
    if op.target.kind == PLANET:
        for index, moon in enumerate(sorted(obj.moons, key=lambda m: m.distance)):
            if moon.name == f"{old}{moon_letters(index)}":        # a moon the admin named by hand keeps its name
                moon.name = f"{value}{moon_letters(index)}"
    return edits.EditResult(f"Renamed {old} to {value}.", [], [], [], [])


def _set_star_type(sim, op, system, obj, owner):
    try:
        result, _star = edits.change_star(system, str(op.value))
    except ValueError as exc:
        raise PlanError(str(exc)) from exc
    return result


ASPECTS = {
    "name": Aspect("name", (SYSTEM, STAR, PLANET, MOON), _set_name),
    "star_type": Aspect("star_type", (STAR,), _set_star_type),
}
"""dict: The aspects `SetAspect` can set, by name. The aspect catalogue (ADM.62) grows it with `register_aspect`."""


def register_aspect(name, kinds, apply):
    """
    Adds an aspect `SetAspect` can set.

    Args:
        name (str): The aspect's name.
        kinds (tuple): The kinds of object it applies to.
        apply (callable): `apply(sim, op, system, obj, owner)`, which changes `obj` in the simulated `system` and
            returns an `edits.EditResult`, or raises `PlanError`.
    """
    ASPECTS[name] = Aspect(name, tuple(kinds), apply)


def _apply_set_aspect(sim, op):
    aspect = ASPECTS.get(op.aspect)
    if aspect is None:
        raise PlanError(f"unknown aspect {op.aspect!r}; the aspects are {', '.join(sorted(ASPECTS))}")
    if op.target.kind not in aspect.kinds:
        raise PlanError(f"the {aspect.name} of a {op.target.kind} cannot be set; it applies to "
                        f"{', '.join(aspect.kinds)}")
    system, obj, owner = _locate(sim, op.target)
    return aspect.apply(sim, op, system, obj, owner)


def _apply_regenerate(sim, op):
    if op.target.kind not in (PLANET, MOON, BELT):
        raise PlanError(f"a {op.target.kind} cannot be regenerated by a plan yet")
    keep = tuple(op.keep)
    unknown = [item for item in keep if item not in KEEPABLE]
    if unknown:
        raise PlanError(f"cannot keep {', '.join(unknown)}; the aspects to keep are {', '.join(KEEPABLE)}")
    if "name" not in keep or "orbit" not in keep:
        raise PlanError("a regeneration keeps the name and the orbit for now; redrawing them comes later"
                        "rearranger")
    system, body, owner = _locate(sim, op.target)
    return edits.regenerate_body(system, op.target.kind, body, owner)


def _apply_remove(sim, op):
    if op.target.kind not in (PLANET, MOON, BELT):
        raise PlanError(f"a {op.target.kind} cannot be removed by a plan yet")
    system, body, owner = _locate(sim, op.target)
    return edits.remove_body(system, body, owner)


def _add_planet(sim, op, system):
    if op.target.kind != SYSTEM:
        raise PlanError("a planet is added to a system")
    name = str(op.name or "").strip()
    if not name:
        raise PlanError("a new planet needs a name")
    distance = _number(op.distance_au, "the distance")
    star = system.star
    ceiling = validation.orbit_ceiling_au(star)
    if distance > ceiling * (1 + validation.RELATIVE_TOLERANCE):
        raise PlanError(f"{distance:.4g} AU is beyond the farthest stable orbit ({ceiling:.4g} AU)")
    planet = Planet(system.system_config, star, star.habitable_zone, distance)
    if planet.distance != distance:
        planet.distance = distance
        validation.reconcile_moved_planet(planet)
    planet.name, planet.db_id = name, op.new_id
    system.planets.append(planet)
    system.planets.sort(key=lambda item: item.distance)
    _letter_moons(planet)
    for moon in planet.moons:
        radiation.apply_parent_belts(moon, planet)  # GEN.983
    for body in [planet] + planet.moons:
        validation.reapply_life(body)
    magnetism.generate_stellar_induction(planet)  # GEN.981
    for moon in planet.moons:
        moon.db_id = f"{op.new_id}-{moon.name}"
        magnetism.generate_induced_field(moon, planet)
    return _settle(system, f"Added {name} at {distance:.4g} AU.")


def _add_moon(sim, op, system):
    if op.target.kind != PLANET:
        raise PlanError("a moon is added to a planet")
    _system, planet, _owner = _locate(sim, op.target)
    distance = _number(op.distance_au, "the distance")
    low_km, high_km = planetPhysics.moon_orbit_bounds_km(planet)
    if not _km_to_au(low_km) <= distance <= _km_to_au(high_km) and not op.forced:
        raise PlanError(f"a moon of {planet.name} can sit {_km_to_au(low_km):.3g} to {_km_to_au(high_km):.3g} AU "
                        f"from it, not {distance:.3g}")
    moon = edits.make_moon(planet, distance)
    if moon is None:
        raise PlanError(f"{planet.name} can hold no moon")
    moon.db_id = op.new_id
    planet.moons.append(moon)
    _letter_moons(planet)
    validation.reapply_life(moon)
    validation.refresh_planet_reflex(planet, force=True)
    return _settle(system, f"Added {moon.name} (class {moon.planet_class}) at {distance:.4g} AU from {planet.name}.")


def _add_belt(sim, op, system):
    if op.target.kind != SYSTEM:
        raise PlanError("a belt is added to a system")
    distance = _number(op.distance_au, "the distance")
    lower, upper = _number(op.lower_au, "the lower limit"), _number(op.upper_au, "the upper limit")
    if not lower < distance < upper:
        raise PlanError("a belt's distance lies between its lower and upper limits")
    star = system.star
    if upper > validation.orbit_ceiling_au(star) * (1 + validation.RELATIVE_TOLERANCE):
        raise PlanError(f"{upper:.4g} AU is beyond the farthest stable orbit")
    belt = AsteroidBelt(system.system_config, star, distance, lower, upper)
    belt.db_id = op.new_id
    system.planets.append(belt)
    system.planets.sort(key=lambda item: item.distance)
    return _settle(system, f"Added an asteroid belt at {distance:.4g} AU.")


def _apply_add(sim, op):
    if op.new_id is None:
        raise PlanError("an added object needs its provisional id; add the operation with Plan.add")
    system = sim.system(op.target)
    if op.what == PLANET:
        return _add_planet(sim, op, system)
    if op.what == MOON:
        return _add_moon(sim, op, system)
    if op.what == BELT:
        return _add_belt(sim, op, system)
    raise PlanError(f"cannot add a {op.what}; add a planet, a moon or a belt")


APPLY = {Move.OP: _apply_move, SetAspect.OP: _apply_set_aspect, Regenerate.OP: _apply_regenerate,
         Add.OP: _apply_add, Remove.OP: _apply_remove}
"""dict: The function that applies each operation to the simulated systems."""


def simulate(world, operations, seed):
    """
    Applies `operations`, in order, to copies of the systems of `world`.

    Args:
        world (World): What is stored.
        operations (list): `Operation`s.
        seed: The plan's seed; operation `i` draws from the stream `plan:<seed>:<i>`, so the same plan on the same
            stored state is the same plan.

    Returns:
        tuple: `(sim, results)`: the copies and an `edits.EditResult` per operation.

    Raises:
        PlanError: If an operation cannot be applied (the message says which, and why).
    """
    sim = _Sim(world)
    results = []
    for index, op in enumerate(operations):
        try:
            with draw.bound(draw.Stream(f"plan:{seed}:{index}")):
                results.append(APPLY[op.OP](sim, op))
        except PlanError as exc:
            raise PlanError(f"operation {index + 1} ({op.OP}): {exc}") from exc
    return sim, results


# ---------------------------------------------------------------------
# The effect diff
# ---------------------------------------------------------------------

@dataclass(frozen=True)
class Change:
    """
    One row of the effect: `field` of `ref` goes from `before` to `after`. An object that is added or removed is one
    change with `field` `"*"` and the object's summary on the side that has it. `indirect` is true for an object no
    operation was aimed at (a neighbour the stabilising moved, a moon that was renamed).
    """

    ref: Ref
    field: str
    before: object
    after: object
    indirect: bool = False

    def to_dict(self):
        return {"ref": self.ref.to_dict(), "field": self.field, "before": self.before, "after": self.after,
                "indirect": self.indirect}

    @classmethod
    def from_dict(cls, data):
        return cls(Ref.from_dict(data["ref"]), data["field"], data["before"], data["after"], data.get("indirect", False))


@dataclass
class Effect:
    """What a plan does: its `changes`, a sentence per operation (`summaries`), the layout problems it leaves that
    were not there before (`warnings`), and the keys of the objects it aims at (`direct`)."""

    changes: list = field(default_factory=list)
    summaries: list = field(default_factory=list)
    warnings: list = field(default_factory=list)
    direct: list = field(default_factory=list)

    @property
    def touched(self):
        """list: The keys of every object with a change, directly or indirectly."""
        return sorted({change.ref.key for change in self.changes})

    @property
    def indirect(self):
        """list: The keys of the objects touched that no operation was aimed at."""
        return sorted(set(self.touched) - set(self.direct))

    def to_dict(self):
        return {"changes": [change.to_dict() for change in self.changes], "summaries": list(self.summaries),
                "warnings": list(self.warnings), "direct": list(self.direct)}

    @classmethod
    def from_dict(cls, data):
        return cls([Change.from_dict(item) for item in data["changes"]], list(data["summaries"]),
                   list(data["warnings"]), list(data["direct"]))


def _flat(obj):
    data = obj.to_dict()
    data.pop("moons", None)
    return data


def _summary(obj):
    return {"name": getattr(obj, "name", None), "class": getattr(obj, "planet_class", None),
            "distance": getattr(obj, "distance", None)}


def _index(system_id, system):
    """`{key: (ref, obj)}` for the system and everything in it."""
    found = {}

    def put(kind, ident, obj):
        ref = Ref(system_id, kind, ident)
        found[ref.key] = (ref, obj)

    found[Ref(system_id).key] = (Ref(system_id), system)
    for index, star in enumerate(system.stars):
        put(STAR, _ident(star, f"@{index}"), star)
    for star_index, (_star, planets) in enumerate(validation.star_lists(system)):
        for index, body in enumerate(planets):
            if body.body_type == 'a':
                put(BELT, _ident(body, f"@{star_index}.{index}"), body)
                continue
            put(PLANET, _ident(body, f"@{star_index}.{index}"), body)
            for moon_index, moon in enumerate(body.moons):
                put(MOON, _ident(moon, f"@{star_index}.{index}.{moon_index}"), moon)
    for star_index, (_star, comets) in enumerate(validation.comet_lists(system)):
        for index, comet in enumerate(comets):
            put(COMET, _ident(comet, f"@{star_index}.{index}"), comet)
    return found


def _system_fields(system):
    """A system's own fields, without the objects the index holds apart."""
    data = system.to_dict()
    for name in ("star", "secondary_star", "wide_binary", "planets", "secondary_planets", "comets",
                 "secondary_comets", "member_stars", "member_planets", "system_config"):
        data.pop(name, None)
    return data


def compute_effect(world, sim, results, operations):
    """
    The effect of a simulated plan: `sim` and `results` are `simulate`'s, `operations` the plan's.

    Returns:
        Effect: One `Change` per field that differs between the stored systems and the copies, in the order of
        the plan's systems and objects.
    """
    direct = _target_keys(operations)
    changes = []
    for system_id, after_system in sim.systems.items():
        before = _index(system_id, world.systems[system_id])
        after = _index(system_id, after_system)
        for key in list(before) + [key for key in after if key not in before]:
            ref = (before.get(key) or after.get(key))[0]
            indirect = key not in direct
            if key not in after:
                changes.append(Change(ref, "*", _summary(before[key][1]), None, indirect))
            elif key not in before:
                changes.append(Change(ref, "*", None, _summary(after[key][1]), indirect))
            else:
                old, new = (_system_fields if ref.kind == SYSTEM else _flat)(before[key][1]), \
                    (_system_fields if ref.kind == SYSTEM else _flat)(after[key][1])
                for name in old.keys() | new.keys():
                    if old.get(name) != new.get(name):
                        changes.append(Change(ref, name, old.get(name), new.get(name), indirect))
    for system_id, position in sim.placements.items():
        if position != world.placements.get(system_id):
            ref = Ref(system_id)
            changes.append(Change(ref, "position", world.placements.get(system_id), position, ref.key not in direct))
    baseline = {system_id: {(p.body, p.message) for p in validation.check_star_system(world.systems[system_id])}
                for system_id in sim.systems}
    warnings = []
    for system_id, after_system in sim.systems.items():
        for problem in validation.check_star_system(after_system):
            if (problem.body, problem.message) not in baseline[system_id]:
                warnings.append(f"{problem.body}: {problem.message}.")
    return Effect(changes, [result.summary for result in results if result is not None], warnings, sorted(direct))


# ---------------------------------------------------------------------
# The plan
# ---------------------------------------------------------------------

@dataclass(frozen=True)
class Transition:
    """One move of a plan between states, with who made it and in which role."""

    from_state: str
    to_state: str
    who: str
    role: str
    at: float
    note: str = None

    def to_dict(self):
        return dataclasses.asdict(self)

    @classmethod
    def from_dict(cls, data):
        return cls(**data)


@dataclass(frozen=True)
class Drift:
    """A planned move whose start no longer is where the object is: `planned` was its `field`, `current` is now."""

    ref: Ref
    field: str
    planned: object
    current: object


def _current_orbit(world, ref):
    """What a `Move` of `ref` starts from now, as `{field: value}`."""
    if ref.kind == SYSTEM:
        return {"position": world.placements.get(ref.system)}
    if ref.kind not in (PLANET, MOON, BELT):
        raise PlanError(f"a {ref.kind} cannot be moved")
    system = world.systems.get(ref.system)
    if system is None:
        raise PlanError(f"no such system: {ref.system}")
    try:
        body, _owner = edits.find_body(system, ref.kind, ref.id)
    except edits.BodyNotFound as exc:
        raise PlanError(str(exc)) from exc
    return {"distance_au": body.distance}


def _same(a, b):
    if isinstance(a, (tuple, list)) and isinstance(b, (tuple, list)):
        return len(a) == len(b) and all(_same(x, y) for x, y in zip(a, b))
    if isinstance(a, float) and isinstance(b, float):
        return math.isclose(a, b, rel_tol=1e-12, abs_tol=0.0)
    return a == b


def new_plan_id():
    """A plan's id: sixteen hex digits."""
    return secrets.token_hex(8)


@dataclass
class Plan:
    """
    An ordered list of `operations`, its `state` and `history`, and, once proposed, its `effect`, the `outcomes`
    (each touched system as it will be, in the `StarSystem.to_dict()` envelope) and the `base_stamps` the stale
    check compares with. `seed` makes the effect repeatable.
    """

    id: str = field(default_factory=new_plan_id)
    title: str = ""
    seed: object = 0
    author: str = None
    author_role: str = None
    operations: list = field(default_factory=list)
    state: str = DRAFT
    history: list = field(default_factory=list)
    base_stamps: dict = field(default_factory=dict)
    effect: Effect = None
    outcomes: dict = field(default_factory=dict)

    # --- building ------------------------------------------------------------------------------------------

    def add(self, op):
        """
        Appends an operation (a draft or a requested plan only). An `Add` gets its provisional id here.

        Returns:
            Operation: `op`.
        """
        if self.state not in (DRAFT, REQUESTED):
            raise PlanStateError(f"a {self.state} plan cannot change; withdraw it and start another")
        if isinstance(op, Add) and op.new_id is None:
            taken = {item.new_id for item in self.operations if isinstance(item, Add)}
            number = len(taken) + 1
            while f"new-{number}" in taken:
                number += 1
            op.new_id = f"new-{number}"
        self.operations.append(op)
        return op

    # --- states --------------------------------------------------------------------------------------------

    def _move(self, to_state, who, role, now, note=None):
        if to_state not in TRANSITIONS[self.state]:
            raise PlanStateError(f"a {self.state} plan cannot become {to_state}")
        if role not in ROLES:
            raise PlanStateError(f"unknown role {role!r}")
        self.history.append(Transition(self.state, to_state, who, role, time.time() if now is None else now, note))
        self.state = to_state

    @staticmethod
    def _approver(role, what):
        if role not in APPROVER_ROLES:
            raise PlanStateError(f"only {' or '.join(APPROVER_ROLES)} can {what} a plan, not a {role}")

    def request(self, who, role, now=None):
        """Submits the plan as a request, by someone below the approving role (10.0): draft to requested."""
        if self.state != DRAFT:
            raise PlanStateError(f"a {self.state} plan cannot become requested")
        if not self.operations:
            raise PlanError("a plan needs at least one operation")
        self.author, self.author_role = self.author or who, self.author_role or role
        self._move(REQUESTED, who, role, now)

    def _evaluate(self, world):
        """Works the effect, the outcomes and the stamps out from `world` (raising `PlanError` when an operation
        cannot be applied, with the plan left as it was)."""
        if not self.operations:
            raise PlanError("a plan needs at least one operation")
        sim, results = simulate(world, self.operations, self.seed)
        effect = compute_effect(world, sim, results, self.operations)
        self.effect = effect
        self.outcomes = {system_id: system.to_dict() for system_id, system in sim.systems.items()}
        before = {key for system_id in sim.systems for key in _index(system_id, world.systems[system_id])}
        self.base_stamps = {key: world.stamps.get(key) for key in sorted(set(effect.touched) | set(effect.direct))
                            if key in before}

    def _capture_before(self, world):
        for op in self.operations:
            if isinstance(op, Move):
                op.before = _current_orbit(world, op.target)

    def propose(self, world, who, role, now=None):
        """
        Works the plan out against `world` and puts it up for approval: draft, requested or stale to proposed, by an
        approving role. Records what each move starts from and the stamps of the objects touched.

        Raises:
            PlanError: If an operation cannot be applied; the plan stays as it was.
        """
        if self.state not in (DRAFT, REQUESTED, STALE):
            raise PlanStateError(f"a {self.state} plan cannot become proposed")
        self._approver(role, "propose")
        self._capture_before(world)
        self._evaluate(world)
        self.author, self.author_role = self.author or who, self.author_role or role
        self._move(PROPOSED, who, role, now)
        return self.effect

    def stale_keys(self, world):
        """The keys of the objects the plan touches whose stamp is not what it was when the plan was proposed."""
        return [key for key, stamp in self.base_stamps.items() if world.stamps.get(key) != stamp]

    def check_stale(self, world, who="engine", now=None):
        """
        Marks a proposed or approved plan stale when an object it touches was edited since it was proposed.

        Returns:
            list: The keys that changed (empty: the plan stands).
        """
        if self.state not in (PROPOSED, APPROVED):
            return []
        keys = self.stale_keys(world)
        if keys:
            self._move(STALE, who, ROLE_SYSTEM, now, note="edited since the plan was proposed: " + ", ".join(keys))
        return keys

    def drift(self, world):
        """
        What an orbital update moved under the plan: for every `Move`, the value it was planned from against the
        object's value now.

        Returns:
            list: `Drift`s, empty when nothing moved.
        """
        found = []
        for op in self.operations:
            if not isinstance(op, Move) or not op.before:
                continue
            current = _current_orbit(world, op.target)
            for name, planned in op.before.items():
                if not _same(planned, current.get(name)):
                    found.append(Drift(op.target, name, planned, current.get(name)))
        return found

    def approve(self, world, who, role, now=None):
        """
        Confirms a proposed plan (the admin who proposed it may). Refused, and the plan marked stale, when an
        object it touches was edited since it was proposed; refused while an orbital update has moved what a move
        starts from (`rebase` first, then look at the effect again).

        Raises:
            PlanStale: If the plan is stale.
            PlanError: If the plan has drifted.
        """
        if self.state != PROPOSED:
            raise PlanStateError(f"a {self.state} plan cannot become approved")
        self._approver(role, "approve")
        keys = self.check_stale(world, who, now)
        if keys:
            raise PlanStale("the plan is stale: " + ", ".join(keys) + " changed after it was proposed")
        drifted = self.drift(world)
        if drifted:
            raise PlanError("an orbital update moved " + ", ".join(item.ref.key for item in drifted)
                            + " since the plan was proposed; rebase it and look at the effect again")
        self._move(APPROVED, who, role, now)

    def rebase(self, world, who, role, now=None):
        """
        Plans again from where the objects are now, after an orbital update moved them (the stamps do not change
        for that): every move starts from the current value and the effect is worked out afresh. An approved plan
        goes back to proposed, so what is approved is what is shown.

        Returns:
            list: The `Drift`s that were rebased.
        """
        if self.state not in (DRAFT, REQUESTED, PROPOSED, APPROVED):
            raise PlanStateError(f"a {self.state} plan cannot be rebased")
        self._approver(role, "rebase")
        drifted = self.drift(world)
        self._capture_before(world)
        self._evaluate(world)
        if self.state == APPROVED:
            self._move(PROPOSED, who, role, now, note="rebased after the orbital update")
        return drifted

    def reject(self, who, role, note=None, now=None):
        """Turns a requested or proposed plan down."""
        self._approver(role, "reject")
        self._move(REJECTED, who, role, now, note)

    def withdraw(self, who, role, note=None, now=None):
        """Takes the plan back: its author, or an approving role, until it is applied."""
        if role not in APPROVER_ROLES and who != self.author:
            raise PlanStateError("only the plan's author or an admin can withdraw it")
        self._move(WITHDRAWN, who, role, now, note)

    def mark_applied(self, who, role, now=None):
        """Records that the approved plan was applied (the apply job, ADM.60, calls it)."""
        self._approver(role, "apply")
        self._move(APPLIED, who, role, now)

    def mark_undone(self, who, role, now=None):
        """Records that the applied plan was undone."""
        self._approver(role, "undo")
        self._move(UNDONE, who, role, now)

    # --- the envelope --------------------------------------------------------------------------------------

    def to_dict(self):
        """The plan as JSON-ready data: operations in their envelope, the effect, the outcomes and the history."""
        return {
            "plan_version": PLAN_VERSION, "schema_version": ENVELOPE_VERSION, "id": self.id, "title": self.title,
            "seed": self.seed, "author": self.author, "author_role": self.author_role, "state": self.state,
            "operations": [op.to_dict() for op in self.operations],
            "history": [item.to_dict() for item in self.history], "base_stamps": dict(self.base_stamps),
            "effect": self.effect.to_dict() if self.effect is not None else None, "outcomes": self.outcomes,
        }

    @classmethod
    def from_dict(cls, data):
        """Rebuilds a plan from `to_dict`'s data."""
        if data.get("plan_version", PLAN_VERSION) > PLAN_VERSION:
            raise PlanError(f"plan version {data['plan_version']} is newer than this code understands (max {PLAN_VERSION})")
        if data["state"] not in STATES:
            raise PlanError(f"unknown plan state {data['state']!r}")
        return cls(id=data["id"], title=data.get("title", ""), seed=data.get("seed", 0), author=data.get("author"),
                   author_role=data.get("author_role"), operations=[Operation.from_dict(item) for item in data["operations"]],
                   state=data["state"], history=[Transition.from_dict(item) for item in data.get("history", [])],
                   base_stamps=dict(data.get("base_stamps", {})),
                   effect=Effect.from_dict(data["effect"]) if data.get("effect") is not None else None,
                   outcomes=dict(data.get("outcomes", {})))
