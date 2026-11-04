# planetgen/physics/system_bodies.py

"""
The ledger of a whole system: one list of every body with a mass, and the barycentre that includes the belts (GEN.1066; docs/design/systems/system-bodies-ledger.md, belt-evolution.md section 5)
================================================================================================================================================================================================

Planets, moons, belts and comets each live in their own typed table, so nothing today can answer "what is the total mass of this system" or "where is its barycentre with the belts in it". This module is the pure half of
the answer: a flat list of `Body` rows (`id`, `parent`, `kind`, `mass_kg`, `radius_km`, `density_kg_m3`, `a_au`, `e`, `i_deg`, `status`, `tier`, `position_au`), the checks that make it a tree, and the sums over it.
The list itself is a **view** over the typed tables (the design note decides view against table: a view, no new storage), so building the rows from the tables is the wiring, not this module.

1. **One tree, one rule.** The root is the central object (a star). Every other body names the body it orbits (`parent`), and its `position_au` is relative to that parent's centre, as the typed tables store it (the unified-system
   rule: stored relative, derived absolute, `unified-system-barycentre.md`). A body's place in the system frame (the root at the origin) is the sum of the offsets down its chain of parents (`absolute_positions_au`).
2. **Kinds and who may be whose parent.** `star` (root, or a companion star about another star), `planet` (about a star), `moon` (about a planet, or about a belt member that has a moon), `belt` (about a star; a *container*,
   its own mass is zero and its mass is the sum of its children), `belt_member` (an L1 body of a belt, one row each), `extended` (what is below the L1 cut: the L2 and L3 tiers of a belt as one extended mass, no
   individual orbit), `comet` (about a star), `ring` (about a planet; one extended mass). An `extended` body and a `ring` with no stored position sit at their parent's centre (a symmetric annulus has its centre of mass
   at its centre).
3. **Mass is conserved.** `subtree_mass_kg` of a belt is the sum of its members and its extended mass; `belt_mass_residual` is that against the belt's own mass (the budget `belt_ledger` spends), so a belt whose ledger
   does not add up is found by a number, not a look.
4. **The barycentre includes everything with a mass.** `barycentre_au` is `sum(m r) / sum(m)` over all bodies in the system frame; the central object sits at minus that, `central_offset_au`. `belt_wobble_au` is how much
   the belts move the central object away from where the planets alone put it, the test the design asked for. Every body is a point; an `extended` body is a point at its stored place.

Limits. Positions are the bodies' places at one instant (the caller gives them; `relative_position_au` turns elements and a mean anomaly into one), and the sums are exact for those points; a time series is a caller loop.
Close binaries (a barycentre as the origin of two stars) keep their own convention in `system_tree.py`; here a companion star is a body about the primary. Nothing calls this yet.

Use:

    rows = [Body(...), ...]                              # from the typed tables
    check(rows)                                          # one root, parents exist, no cycle, kind rules, belts hold no mass of their own
    b = barycentre_au(rows); offset = central_offset_au(rows)
    shift = belt_wobble_au(rows)                         # the star's displacement from the belts alone
"""

import math
from collections import namedtuple

from planetgen.physics.kepler import solve_eccentric_anomaly
from planetgen.util import detmath

Body = namedtuple(
    "Body",
    ("id", "parent", "kind", "mass_kg", "radius_km", "density_kg_m3", "a_au", "e", "i_deg", "status", "tier", "position_au"),
    defaults=(None,) * 8,
)
"""One row of the ledger. `parent` is the id of the body orbited (None for the root); `position_au` is `(x, y, z)` relative to the parent's centre, or None for a body that sits at its parent's centre."""

ROOT_KINDS = ("star",)
"""tuple: Kinds that can be the root."""

PARENT_KINDS = {
    "star": (None, "star"),
    "planet": ("star",),
    "moon": ("planet", "belt_member"),
    "belt": ("star",),
    "belt_member": ("belt",),
    "extended": ("belt",),
    "comet": ("star",),
    "ring": ("planet",),
}
"""dict: The kinds each kind may orbit (None is the root's parent)."""

KINDS = tuple(PARENT_KINDS)
"""tuple: Every kind of row."""

STATUSES = ("star", "planet", "dwarf planet", "minor body", "moon", "comet", "belt", "ring", "extended")
"""tuple: The status words a row may carry (belt-evolution.md section 3.2); None is allowed for a row whose status is not worked out."""

TIERS = ("typed", "L1", "L3")
"""tuple: Where a row comes from: a typed table (`typed`), an individual belt body (`L1`), or the extended remainder (`L3`)."""

ZERO = (0.0, 0.0, 0.0)


class SystemBodiesError(ValueError):
    """The rows are not a valid system tree."""


def _finite(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def check(bodies):
    """
    Check that `bodies` is one tree of valid rows and return the ids in an order where every parent comes before its children.

    Raises `SystemBodiesError` for: a duplicate id; a kind, status or tier not in the lists; a mass that is not finite and non-negative; a belt with a mass of its own; a position that is not three finite numbers;
    no root, or more than one; a parent that is missing or of a kind the row's kind may not orbit; a cycle.
    """
    by_id = {}
    for body in bodies:
        if body.id in by_id:
            raise SystemBodiesError(f"duplicate id {body.id!r}")
        by_id[body.id] = body
        if body.kind not in KINDS:
            raise SystemBodiesError(f"{body.id!r}: unknown kind {body.kind!r}")
        if body.status is not None and body.status not in STATUSES:
            raise SystemBodiesError(f"{body.id!r}: unknown status {body.status!r}")
        if body.tier is not None and body.tier not in TIERS:
            raise SystemBodiesError(f"{body.id!r}: unknown tier {body.tier!r}")
        if not _finite(body.mass_kg) or body.mass_kg < 0:
            raise SystemBodiesError(f"{body.id!r}: mass {body.mass_kg!r} must be finite and not negative")
        if body.kind == "belt" and body.mass_kg != 0:
            raise SystemBodiesError(f"{body.id!r}: a belt is a container, its mass is the sum of its children")
        if body.position_au is not None and (len(body.position_au) != 3 or not all(_finite(v) for v in body.position_au)):
            raise SystemBodiesError(f"{body.id!r}: position must be three finite numbers")
    roots = [body for body in by_id.values() if body.parent is None]
    if len(roots) != 1:
        raise SystemBodiesError(f"a system has exactly one root, found {len(roots)}")
    if roots[0].kind not in ROOT_KINDS:
        raise SystemBodiesError(f"{roots[0].id!r}: the root must be one of {ROOT_KINDS}, not {roots[0].kind!r}")
    for body in by_id.values():
        if body.parent is None:
            continue
        parent = by_id.get(body.parent)
        if parent is None:
            raise SystemBodiesError(f"{body.id!r}: parent {body.parent!r} is not in the ledger")
        if parent.kind not in PARENT_KINDS[body.kind]:
            raise SystemBodiesError(f"{body.id!r}: a {body.kind} cannot orbit a {parent.kind}")
    order, seen = [], set()
    children = {}
    for body in by_id.values():
        children.setdefault(body.parent, []).append(body.id)
    stack = [roots[0].id]
    while stack:
        identifier = stack.pop()
        seen.add(identifier)
        order.append(identifier)
        stack.extend(reversed(children.get(identifier, [])))
    if len(seen) != len(by_id):
        raise SystemBodiesError("a cycle: some bodies cannot be reached from the root")
    return order


def absolute_positions_au(bodies):
    """The place of every body in the system frame (the root's centre at the origin): `{id: (x, y, z)}` in AU, the sum of the offsets down the chain of parents."""
    order = check(bodies)
    by_id = {body.id: body for body in bodies}
    places = {}
    for identifier in order:
        body = by_id[identifier]
        offset = body.position_au if body.position_au is not None else ZERO
        base = places[body.parent] if body.parent is not None else ZERO
        places[identifier] = tuple(base[axis] + offset[axis] for axis in range(3))
    return places


def subtree_mass_kg(bodies, identifier):
    """The mass of a body and everything that orbits it, down the tree (a sub-system counts as the sum of its parts)."""
    order = check(bodies)
    by_id = {body.id: body for body in bodies}
    if identifier not in by_id:
        raise SystemBodiesError(f"{identifier!r} is not in the ledger")
    members = {identifier}
    for name in order:
        if by_id[name].parent in members:
            members.add(name)
    return math.fsum(by_id[name].mass_kg for name in members)


def total_mass_kg(bodies):
    """The mass of the whole system: every body with a mass, once."""
    check(bodies)
    return math.fsum(body.mass_kg for body in bodies)


def belt_mass_residual(bodies, belt_id, belt_mass_kg):
    """
    How far a belt's rows are from its mass: `(sum of its rows - belt_mass_kg) / belt_mass_kg`. Zero when the L1 members and the extended remainder add up to the belt (`belt_ledger` makes them), a positive
    number when the rows hold too much. Raises if `belt_id` is not a belt or `belt_mass_kg` is not positive.
    """
    by_id = {body.id: body for body in bodies}
    if belt_id not in by_id or by_id[belt_id].kind != "belt":
        raise SystemBodiesError(f"{belt_id!r} is not a belt in the ledger")
    if not _finite(belt_mass_kg) or belt_mass_kg <= 0:
        raise SystemBodiesError(f"belt mass {belt_mass_kg!r} must be finite and positive")
    return (subtree_mass_kg(bodies, belt_id) - belt_mass_kg) / belt_mass_kg


def barycentre_au(bodies, without_kinds=()):
    """
    The centre of mass of the system in the system frame (the root's centre at the origin), AU: `sum(m r) / sum(m)` over every body whose kind is not in `without_kinds`. A system whose bodies have no mass in this
    selection is a `SystemBodiesError`.
    """
    places = absolute_positions_au(bodies)
    chosen = [body for body in bodies if body.kind not in without_kinds]
    mass = math.fsum(body.mass_kg for body in chosen)
    if mass <= 0:
        raise SystemBodiesError("no mass to take a barycentre of")
    return tuple(math.fsum(body.mass_kg * places[body.id][axis] for body in chosen) / mass for axis in range(3))


def central_offset_au(bodies, without_kinds=()):
    """Where the root (the central object) sits relative to the barycentre, AU: minus the barycentre in the system frame."""
    return tuple(-value for value in barycentre_au(bodies, without_kinds))


def belt_wobble_au(bodies, kinds=("belt_member", "extended")):
    """
    How far the belts move the central object: the distance, AU, between where the barycentre puts the root with everything and where it puts it with the bodies of `kinds` (the belts' members and remainders) left out.
    Belts that are round and centred contribute nothing; a heavy member on one side does.
    """
    with_belts = central_offset_au(bodies)
    without = central_offset_au(bodies, kinds)
    return detmath.dist(with_belts, without)


def relative_position_au(a_au, e, i_deg, node_deg, argument_deg, mean_anomaly_deg):
    """
    Turn a Kepler orbit and a mean anomaly into the body's position relative to its parent, AU: the focus at the origin, the reference plane the x-y plane, the node and argument of periapsis and the inclination
    applied in the usual order. Raises for `a_au <= 0` or `e` outside `[0, 1)`.
    """
    if not _finite(a_au) or a_au <= 0:
        raise SystemBodiesError(f"semi-major axis {a_au!r} must be finite and positive")
    if not _finite(e) or not 0 <= e < 1:
        raise SystemBodiesError(f"eccentricity {e!r} must be in [0, 1)")
    anomaly = solve_eccentric_anomaly(math.radians(mean_anomaly_deg), e)
    x_orbit = a_au * (math.cos(anomaly) - e)
    y_orbit = a_au * math.sqrt(1.0 - e * e) * math.sin(anomaly)
    cos_w, sin_w = math.cos(math.radians(argument_deg)), math.sin(math.radians(argument_deg))
    cos_n, sin_n = math.cos(math.radians(node_deg)), math.sin(math.radians(node_deg))
    cos_i, sin_i = math.cos(math.radians(i_deg)), math.sin(math.radians(i_deg))
    x1, y1 = x_orbit * cos_w - y_orbit * sin_w, x_orbit * sin_w + y_orbit * cos_w
    return (
        x1 * cos_n - y1 * cos_i * sin_n,
        x1 * sin_n + y1 * cos_i * cos_n,
        y1 * sin_i,
    )
