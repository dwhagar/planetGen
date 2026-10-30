# planetgen/physics/planet_zones.py

"""
Where a planet can stay in a hierarchy of stars, and the light it gets there (GEN.314).

One pure function over a member tree (`OrbitNode` and `TreeMember` of `physics/system_tree.py`, or the generator's
`multistar.Node` and `multistar.Member`: any objects with the same fields), the one place the zones are worked out.
GEN.264 stage C, GEN.130 and GEN.275 call it and compute no zones of their own. The rules are
docs/design/systems/multistar-and-compact-systems.md 5.6 and 6.5:

- **S-type**, round one star X: the window runs from X's floor (the engulfment or Roche floor the caller gives, 0
  when it gives none) out to the Holman and Wiegert S-type limit of X against its sibling subtree, taken as a point
  mass at the separation and eccentricity of the node they share, and never past X's own Hill radius in that pair
  (the tighter of the two, as a wide pair's `Star.a_crit_au` has always been). A lone star's window is its floor out
  to the disk cap.
- **P-type**, round a pair node P (the two subtrees below it as points): from the Holman and Wiegert P-type limit
  (the lighter subtree's fraction, P's eccentricity) or the highest floor of a member below it, whichever is further
  out, up to the S-type limit of P as a body against its own sibling in its parent's orbit. The root has no sibling,
  so its window runs to the disk cap. The window is empty when the parent's axis is too near P's: a hierarchy needs
  an axis ratio of about 8 for a circumbinary zone.
- Everything is capped at the disk size (`tuning.PLANET_ZONE_DISK_CAP_AU`, 100 AU).

A `group` node (a young system that is not hierarchical) has no orbit to set a window by: its members get none and a
node under it is taken as a root.

`flux_at` is the other half: the time-averaged light on a planet of a zone's host, in units of what one solar
luminosity gives at 1 AU. The host's own light falls with the planet's distance squared; every member outside the host
adds `L <1/d^2>` with `d` the separation of the node that member and the host meet at and
`<1/d^2> = 1 / (a^2 sqrt(1 - e^2))` over an eccentric orbit (6.5). Both treat a subtree as a point at its centre of
mass, which holds where a planet is well inside its host's window and the hierarchy is as hierarchical as the
stability test makes it.

Distances are AU, masses and luminosities solar.
"""

import math
from dataclasses import dataclass

from planetgen import tuning
from planetgen.physics import unit_registry as units
from planetgen.physics.orbits import (
    hill_radius_in_pair_au, holman_wiegert_circumbinary_a_crit_au, holman_wiegert_critical_semimajor_axis)

_km_to_au = units.converter("km", "au")

S_TYPE = "S"
"""str: A zone round one star."""

P_TYPE = "P"
"""str: A zone round a pair node."""


@dataclass(frozen=True)
class Zone:
    """
    One window. `host_no` is the member's number for an S-type zone and the node's for a P-type one. `inner_au` and
    `outer_au` are the closest and farthest stable orbit; the window is empty when `outer_au <= inner_au`.
    """
    kind: str
    host_no: int
    inner_au: float
    outer_au: float

    @property
    def width_au(self):
        return max(0.0, self.outer_au - self.inner_au)

    @property
    def exists(self):
        return self.outer_au > self.inner_au


class _Tree:
    """The members and nodes indexed by number, with each subtree's mass and members."""

    def __init__(self, nodes, members):
        self.nodes = {node.node_no: node for node in nodes}
        self.members = {member.member_no: member for member in members}
        self.children = {}
        for member in self.members.values():
            if member.parent_node_no is not None:
                self.children.setdefault(member.parent_node_no, {})[member.slot] = ("member", member.member_no)
        for node in self.nodes.values():
            if node.parent_node_no is not None:
                self.children.setdefault(node.parent_node_no, {})[node.slot] = ("node", node.node_no)
        self._below = {}

    def below(self, child):
        """The member numbers in the subtree of `child`, a `("member" | "node", number)` pair."""
        if child not in self._below:
            kind, number = child
            if kind == "member":
                self._below[child] = frozenset((number,))
            else:
                found = set()
                for grandchild in self.children.get(number, {}).values():
                    found |= self.below(grandchild)
                self._below[child] = frozenset(found)
        return self._below[child]

    def mass(self, child):
        return sum(self.members[number].mass_solar for number in self.below(child))

    def parent_of(self, child):
        """The node `child` hangs from, or None for the root (or a member with no node)."""
        kind, number = child
        parent_no = (self.members if kind == "member" else self.nodes)[number].parent_node_no
        return None if parent_no is None else self.nodes[parent_no]

    def sibling(self, child, parent):
        slot = (self.members if child[0] == "member" else self.nodes)[child[1]].slot
        return self.children[parent.node_no].get(1 - slot)

    def chain(self, child):
        """The node numbers from `child`'s parent up to the root."""
        found = []
        parent = self.parent_of(child)
        while parent is not None:
            found.append(parent.node_no)
            parent = self.parent_of(("node", parent.node_no))
        return found


def _axis_au(node):
    return _km_to_au(node.separation_km)


def _s_limit_au(tree, child, parent):
    """The S-type limit of subtree `child` against its sibling in `parent`'s orbit: the Holman and Wiegert fit or the
    Hill radius in the pair, whichever is tighter."""
    sibling = tree.sibling(child, parent)
    if sibling is None:
        raise ValueError(f"node {parent.node_no} needs a child in slot 0 and in slot 1")
    own, other = tree.mass(child), tree.mass(sibling)
    total = own + other
    axis = _axis_au(parent)
    return min(holman_wiegert_critical_semimajor_axis(axis, other / total, parent.eccentricity),
               hill_radius_in_pair_au(axis, own / total, parent.eccentricity))


def planet_zones(nodes, members, floors_au=None, disk_cap_au=None):
    """
    Every S-type window (one per member) and every P-type window (one per pair node) of a tree.

    Args:
        nodes (iterable): The tree's nodes (empty for a lone object).
        members (iterable): Its members.
        floors_au (dict, optional): `{member_no: AU}` inside which a member has removed every planet (its
            engulfment or Roche floor); a member left out has none.
        disk_cap_au (float, optional): The disk size; `tuning.PLANET_ZONE_DISK_CAP_AU` when left out.

    Returns:
        list of Zone: Members' S-type zones in member order, then the P-type zones in node order. A member of a
            `group` node has no zone and a `group` node has none.

    Raises:
        ValueError: For a node without a child in each slot.
    """
    tree = _Tree(list(nodes), list(members))
    floors = floors_au or {}
    cap = tuning.PLANET_ZONE_DISK_CAP_AU if disk_cap_au is None else disk_cap_au
    zones = []
    for number in sorted(tree.members):
        parent = tree.parent_of(("member", number))
        floor = floors.get(number, 0.0)
        if parent is None:
            zones.append(Zone(S_TYPE, number, floor, cap))
        elif parent.kind != "group":
            zones.append(Zone(S_TYPE, number, floor, min(cap, _s_limit_au(tree, ("member", number), parent))))
    for number in sorted(tree.nodes):
        node = tree.nodes[number]
        if node.kind == "group":
            continue
        low, high = tree.children.get(number, {}).get(0), tree.children.get(number, {}).get(1)
        if low is None or high is None:
            raise ValueError(f"node {number} needs a child in slot 0 and in slot 1")
        mass_low, mass_high = tree.mass(low), tree.mass(high)
        inner = holman_wiegert_circumbinary_a_crit_au(
            _axis_au(node), min(mass_low, mass_high) / (mass_low + mass_high), node.eccentricity)
        inner = max([inner] + [floors.get(member_no, 0.0) for member_no in tree.below(low) | tree.below(high)])
        parent = tree.parent_of(("node", number))
        outer = cap if parent is None or parent.kind == "group" else min(
            cap, _s_limit_au(tree, ("node", number), parent))
        zones.append(Zone(P_TYPE, number, inner, outer))
    return zones


def ancestor_limits_au(zone, nodes, members):
    """
    Every S-type limit above `zone`'s host, nearest first: the host (a star, or a close pair taken as one body) against the stars beside it at
    its own node, then against those beside that node, up to the root. `planet_zones` takes only the nearest one; a planet must stay inside
    all of them (GEN.1014).
    """
    tree = _Tree(list(nodes), list(members))
    child = ("member", zone.host_no) if zone.kind == S_TYPE else ("node", zone.host_no)
    limits = []
    parent = tree.parent_of(child)
    while parent is not None and parent.kind != "group":
        limits.append(_s_limit_au(tree, child, parent))
        child = ("node", parent.node_no)
        parent = tree.parent_of(child)
    return limits


def host_mass_solar(zone, nodes, members):
    """The mass a body in `zone` goes round, in solar masses: the host star's, or the sum of the pair node's two sides'."""
    tree = _Tree(list(nodes), list(members))
    host = ("member", zone.host_no) if zone.kind == S_TYPE else ("node", zone.host_no)
    return tree.mass(host)


def pair_sides_solar(nodes, members, node_no):
    """The masses of the two sides of pair node `node_no`, in solar masses: `(slot 0, slot 1)`."""
    tree = _Tree(list(nodes), list(members))
    sides = tree.children[node_no]
    return tree.mass(sides[0]), tree.mass(sides[1])


def sibling_of(nodes, members, member_no):
    """
    The node a member hangs from and the numbers of the members in the subtree beside it, which its planets see as the
    companion (GEN.332): `(node, frozenset of member numbers)`.

    Raises:
        ValueError: For a member that hangs from no node.
    """
    tree = _Tree(list(nodes), list(members))
    parent = tree.parent_of(("member", member_no))
    if parent is None:
        raise ValueError(f"member {member_no} hangs from no node")
    return parent, tree.below(tree.sibling(("member", member_no), parent))


def flux_at(zone, distance_au, luminosities_sol, nodes, members):
    """
    The light on a planet that goes round `zone`'s host at `distance_au`, in units of one solar luminosity at 1 AU.

    Args:
        zone (Zone): The window the planet is in.
        distance_au (float): The planet's distance from the host (the star for an S-type zone, the pair's centre of
            mass for a P-type one).
        luminosities_sol (dict): `{member_no: luminosity in solar units}`; every member of the tree needs one.
        nodes (iterable): The tree's nodes.
        members (iterable): Its members.

    Returns:
        float: The host's own light plus, for every other member, its luminosity over the time-averaged square of
            its distance (`1 / (a^2 sqrt(1 - e^2))` at the node the member and the host meet at).

    Raises:
        ValueError: For a distance that is not positive, or a member that meets the host at a `group` node.
    """
    if not distance_au > 0.0:
        raise ValueError(f"flux_at: the planet's distance must be positive, got {distance_au!r}")
    tree = _Tree(list(nodes), list(members))
    host = ("member", zone.host_no) if zone.kind == S_TYPE else ("node", zone.host_no)
    own = tree.below(host)
    host_chain = set(tree.chain(host)) | ({zone.host_no} if zone.kind == P_TYPE else set())
    total = sum(luminosities_sol[number] for number in own) / (distance_au * distance_au)
    for number in tree.members:
        if number in own:
            continue
        meeting = next(node_no for node_no in tree.chain(("member", number)) if node_no in host_chain)
        node = tree.nodes[meeting]
        if node.kind == "group":
            raise ValueError(f"member {number} meets the host at group node {meeting}, which has no orbit")
        axis = _axis_au(node)
        total += luminosities_sol[number] / (axis * axis * math.sqrt(1.0 - node.eccentricity ** 2))
    return total
