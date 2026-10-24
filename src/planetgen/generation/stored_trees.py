# planetgen/generation/stored_trees.py

"""
A system built from a tree of stars (GEN.315, GEN.631).

A scatter row with companions (`scatter_companions`) stores a member tree: the placed object is member 0, the others
are drawn from their own seeds and birth masses, and the orbits between them are nodes. A system the fill makes from
nothing draws the same tree (`multistar.draw_member_tree`) around its primary. Either way the tree is the system's stars
and orbits: `plan` turns a stored tree, and `drawn` a tree just drawn, into what `StarSystem(member_tree=...)` takes, and
the build draws no companion of its own.

The pair classes (`binary.BinaryStarProxy`, `wide_binary.WideBinaryPair`) hold exactly two ordinary stars, so a system
built from a tree has one pair, the two stars of a node whose children are both stars (the node of the placed
object when it is one, else the heaviest), and every other star hangs in the tree as one of the system's `member_stars` (GEN.630). A tree
builds when each member is a star `Star` makes from a birth mass and an age. Any other stored tree (a neutron star, a
black hole or a brown dwarf among the members) builds the placed object alone and is marked `complete=False`, its stored
rows kept; a drawn tree has such members cut out before it is built (`multistar.without_members`).

The stored phases hold at the time the tree was drawn; `plan` moves the nodes on to the epoch the built system holds its
orbits at, so the nodes the system stores are the same orbits at that instant. A drawn tree holds its phases at
`GENERATED_EPOCH_UNIX`, which the store replaces with the database's orbit epoch when it saves the system.
"""

import collections
import dataclasses

from planetgen.physics import constants
from planetgen.physics.kepler import wrapped_turn_deg
from planetgen.physics.system_tree import OrbitNode
from planetgen.physics import unit_registry as units

_yr_to_s = units.converter("yr", "s")

PAIR_STAR_KINDS = ("main-sequence", "pre-main-sequence", "protostar", "giant", "white-dwarf")
"""tuple: The member kinds `Star` builds from a birth mass and an age, which the pair classes can hold."""

UNBUILT_KINDS = ("brown-dwarf", "neutron-star", "black-hole")
"""tuple: The member kinds of a drawn tree that `Star` cannot be built as: there is no brown-dwarf host class, and a
compact remnant is a `phenomena.compact_remnant` object, not a star of a pair (`multistar.without_members` cuts them out)."""

PAIR_NODE_KINDS = ("close", "wide")
"""tuple: The node kinds the pair classes stand for (a P-type pair is `close`, an S-type pair `wide`)."""

MASS_TOLERANCE = 1.0e-6
"""float: How far a built companion's mass may differ from its stored mass, as a share of it. A pair the evolution
changed in a way the stored mass records (mass transfer, a merger) does not build from a birth mass and an age."""


GENERATED_EPOCH_UNIX = 0.0
"""float: The epoch a tree drawn at generation holds its phases at, until the store saves the system and puts the database's
orbit epoch there (no real orbit epoch is the start of 1970)."""


@dataclasses.dataclass(frozen=True)
class PlannedStar:
    """One member of a plan: where it hangs in the plan's nodes (`parent_node_no`, `slot`) and what builds it (`seed`,
    `initial_mass_solar`; `mass_solar` is the mass the tree gave it)."""
    member_no: int
    parent_node_no: int
    slot: int
    kind: str
    mass_solar: float
    initial_mass_solar: float
    seed: int


@dataclasses.dataclass(frozen=True)
class StoredPlan:
    """
    What a build takes from a tree. `members` are all the stars of the tree in member order, the placed object (member 0,
    already built) first, and `nodes` its orbits, numbered from 0 in the tree's order with the phases held at the build's
    epoch; both are empty for a single star. `complete` says every stored member is in the system built, so the row's
    stored members and nodes can go; it is `False` for a tree the build cannot hold, which builds the placed object alone.
    `exact` says a built star must have the mass the tree gave it (a stored tree: a star that does not is a change the
    birth mass and age do not give); without it (a tree just drawn) the nodes take the masses the build made.
    """
    members: tuple = ()
    nodes: tuple = ()
    complete: bool = True
    exact: bool = True

    def pair_node_no(self):
        """The node of the system's pair (module docstring), `None` for a single star: the node of the placed star when its
        children are both stars, else the one of the heaviest two stars."""
        if not self.nodes:
            return None
        counts = collections.Counter(member.parent_node_no for member in self.members)
        placed_at = self.members[0].parent_node_no
        if counts[placed_at] == 2:
            return placed_at
        mass = collections.Counter()
        for member in self.members:
            mass[member.parent_node_no] += member.mass_solar
        return max((node_no for node_no, count in counts.items() if count == 2), key=lambda node_no: mass[node_no])


SINGLE = StoredPlan()
"""StoredPlan: A tree with no members: one star, and no companion drawn."""


def node_at(node, epoch_unix):
    """`node` with its mean anomaly carried on from its own epoch to `epoch_unix`; the same orbit at that instant."""
    if epoch_unix is None or epoch_unix == node.epoch_unix:
        return node
    turn = wrapped_turn_deg(epoch_unix - node.epoch_unix, _yr_to_s(node.period_years)) \
        if node.period_years > 0 else 0.0
    return dataclasses.replace(node, mean_anomaly_deg=(node.mean_anomaly_deg + turn) % 360.0,
                               epoch_unix=float(epoch_unix))


def _planned(members, nodes, node_for, exact):
    """The plan of `members` and `nodes` (anything with the fields of `PlannedStar`, and `OrbitNode`s), the node numbers
    taken from 0 in node order and each node made by `node_for`."""
    ordered = sorted(nodes, key=lambda node: node.node_no)
    number = {node.node_no: index for index, node in enumerate(ordered)}
    planned_nodes = tuple(dataclasses.replace(node_for(node), node_no=number[node.node_no],
                                              parent_node_no=number.get(node.parent_node_no)) for node in ordered)
    planned_members = tuple(PlannedStar(m.member_no, number[m.parent_node_no], m.slot, m.kind, m.mass_solar,
                                        m.initial_mass_solar, m.seed) for m in sorted(members, key=lambda m: m.member_no))
    return StoredPlan(planned_members, planned_nodes, True, exact)


def plan(tree, epoch_unix=None):
    """
    The `StoredPlan` of a row's stored tree (`db.scatter_systems.StoredTree`, or `None` for a row with no members), its
    nodes carried to `epoch_unix` (the orbit epoch the built system holds its phases at; `None`: the nodes' own).
    """
    if tree is None or not tree.members:
        return SINGLE
    if len(tree.members) < 2 or any(member.kind not in PAIR_STAR_KINDS for member in tree.members[1:]) \
            or any(node.kind not in PAIR_NODE_KINDS for node in tree.nodes):
        return StoredPlan(complete=False)
    return _planned(tree.members, tree.nodes, lambda node: node_at(node, epoch_unix), True)


def drawn(tree):
    """The `StoredPlan` of a `multistar.MemberTree` just drawn and pruned (`multistar.without_members`) of the members a
    build cannot hold: its phases held at `GENERATED_EPOCH_UNIX`, its masses the build's to correct."""
    if tree.member_count < 2:
        return SINGLE
    return _planned(tree.members, tree.nodes, lambda node: OrbitNode(
        node_no=node.node_no, parent_node_no=node.parent_node_no, slot=node.slot, mass_solar=node.mass_solar,
        separation_km=node.separation_km, eccentricity=node.eccentricity, inclination_deg=node.inclination_deg,
        ascending_node_deg=node.ascending_node_deg, arg_periapsis_deg=node.arg_periapsis_deg,
        mean_anomaly_deg=node.mean_anomaly_deg, epoch_unix=GENERATED_EPOCH_UNIX, period_years=node.period_years,
        kind=node.kind, interaction=node.interaction, flags=node.flags), False)


def matches(star, member):
    """Whether the built `star` is the stored `member`: its mass is the stored mass."""
    return abs(star.mass / constants.SOLAR_MASS_TO_KG - member.mass_solar) <= MASS_TOLERANCE * member.mass_solar
