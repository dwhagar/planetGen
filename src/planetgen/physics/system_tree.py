# planetgen/physics/system_tree.py

"""
Where the members of a multi-star system are at any time (GEN.265), computed
from the tree's stored orbits and never written back.

A scattered system is one point in space, its barycentre (the row's position,
velocity and epoch columns, `system_mass_solar` as its mass), and a tree of
orbits (`scatter_nodes`). Node `n` is the relative orbit of the child in slot 1
about the child in slot 0; a child is a member or another node. Its mean anomaly
at time `t` is `mean_anomaly_deg + 360 fmod(t - epoch_unix, period) / period` (mod 360), so
the nightly update, which only moves the barycentre, never has to write a node
or a member row.

A node's two children sit either side of the node's own centre in proportion to
the other's mass (slot 0 at `-r m1 / (m0 + m1)`, slot 1 at `+r m0 / (m0 + m1)`);
a member's offset from the barycentre is the sum of those along its ancestors.
The members of a `group` node (a young system that is not hierarchical) carry an
offset and a velocity of their own and fly straight from the node's epoch.

`orbit_nodes` of a built system have the same columns, so a fill that copies the
rows takes this same tree over.
"""

import math
from dataclasses import dataclass

from planetgen.physics import unit_registry as units
from planetgen.physics.kepler import comet_orbital_state, wrapped_turn_deg

_km_to_au = units.converter("km", "au")
_yr_to_s = units.converter("yr", "s")
_KM_PER_S_TO_AU_PER_S = _km_to_au(1.0)


@dataclass(frozen=True)
class OrbitNode:
    """One stored orbit, in the units of `scatter_nodes` and `orbit_nodes`: the relative orbit of the slot-1 child about
    the slot-0 child. `separation_km` is its semi-major axis; `interaction` and `flags` are the pair's state (GEN.269)."""
    node_no: int
    parent_node_no: int | None
    slot: int | None
    mass_solar: float
    separation_km: float
    eccentricity: float
    inclination_deg: float
    ascending_node_deg: float
    arg_periapsis_deg: float
    mean_anomaly_deg: float
    epoch_unix: float
    period_years: float
    kind: str = "close"
    interaction: str | None = None
    flags: int = 0


@dataclass(frozen=True)
class TreeMember:
    """One stored member: a child of node `parent_node_no` in `slot`. A member of a `group` node also has an offset
    from the node's centre (AU) and a velocity (km/s)."""
    member_no: int
    parent_node_no: int
    slot: int
    mass_solar: float
    offset_au: tuple | None = None
    velocity_kms: tuple | None = None


def node_relative_au(node, unix):
    """The slot-1 child's position relative to the slot-0 child at `unix`, in AU, in the system's frame."""
    if node.period_years > 0:
        anomaly = node.mean_anomaly_deg + wrapped_turn_deg(unix - node.epoch_unix, _yr_to_s(node.period_years))
    else:
        anomaly = node.mean_anomaly_deg
    axis_au = _km_to_au(node.separation_km)
    state = comet_orbital_state(
        "elliptical", axis_au * (1.0 - node.eccentricity), node.eccentricity, node.inclination_deg,
        node.arg_periapsis_deg, node.ascending_node_deg, node.mass_solar,
        mean_anomaly_rad=math.radians(anomaly % 360.0), orbital_period_years=node.period_years)
    return (state["position_x_au"], state["position_y_au"], state["position_z_au"])


def _add(a, b):
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def _scaled(vector, factor):
    return (vector[0] * factor, vector[1] * factor, vector[2] * factor)


def star_offsets_au(masses_solar, slots, nodes, unix, dark=()):
    """
    The offsets from the barycentre, in AU, of the stars of a system with a tree of `nodes`, in the order of
    `masses_solar` and `slots` (`(node_no, slot)` for each star: the child of that node it is), at `unix`. `dark` lists the
    members that are not stars, `(mass_solar, (node_no, slot))` each (GEN.1013): they weigh in the tree and are left out of
    the offsets returned.
    """
    members = [TreeMember(number, parent, slot, mass)
               for number, (mass, (parent, slot)) in enumerate(zip(masses_solar, slots))]
    members += [TreeMember(len(members) + index, parent, slot, mass) for index, (mass, (parent, slot)) in enumerate(dark)]
    offsets = member_offsets_au(members, nodes, unix)
    return [offsets[number] for number in range(len(masses_solar))]


def member_offsets_au(members, nodes, unix):
    """
    Every member's position at `unix` relative to the system's barycentre, in AU, as `{member_no: (x, y, z)}`.

    Args:
        members (iterable of TreeMember): All the members of the system.
        nodes (iterable of OrbitNode): All its nodes (empty for a single object: the member sits at the origin).
        unix (float): The time, in seconds since the Unix epoch.
    """
    members, nodes = list(members), list(nodes)
    if not nodes:
        return {member.member_no: (0.0, 0.0, 0.0) for member in members}
    children = {}                       # node_no -> {slot: ("member", obj) | ("node", obj)}
    for member in members:
        children.setdefault(member.parent_node_no, {})[member.slot] = ("member", member)
    for node in nodes:
        if node.parent_node_no is not None:
            children.setdefault(node.parent_node_no, {})[node.slot] = ("node", node)
    roots = [node for node in nodes if node.parent_node_no is None]
    if len(roots) != 1:
        raise ValueError(f"a tree has exactly one root node, found {len(roots)}")
    centres = {roots[0].node_no: (0.0, 0.0, 0.0)}   # each node's centre of mass relative to the barycentre
    offsets = {}
    pending = [roots[0]]
    while pending:
        node = pending.pop()
        centre = centres[node.node_no]
        slots = children.get(node.node_no, {})
        if node.kind == "group":
            for kind, child in slots.values():
                if kind == "node":
                    centres[child.node_no] = centre
                    pending.append(child)
                else:
                    fly = tuple(c * (unix - node.epoch_unix) * _KM_PER_S_TO_AU_PER_S for c in (child.velocity_kms or (0.0, 0.0, 0.0)))
                    offsets[child.member_no] = _add(centre, _add(child.offset_au or (0.0, 0.0, 0.0), fly))
            continue
        low, high = slots.get(0), slots.get(1)
        if low is None or high is None:
            raise ValueError(f"node {node.node_no} needs a child in slot 0 and in slot 1")
        mass_low, mass_high = low[1].mass_solar, high[1].mass_solar
        relative = node_relative_au(node, unix)
        total = mass_low + mass_high
        places = ((low, _scaled(relative, -mass_high / total)), (high, _scaled(relative, mass_low / total)))
        for (kind, child), shift in places:
            place = _add(centre, shift)
            if kind == "node":
                centres[child.node_no] = place
                pending.append(child)
            else:
                offsets[child.member_no] = place
    unplaced = [member.member_no for member in members if member.member_no not in offsets]
    if unplaced:
        raise ValueError(f"members {unplaced} hang from no node of the tree")
    return offsets
