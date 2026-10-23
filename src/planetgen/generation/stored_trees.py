# planetgen/generation/stored_trees.py

"""
A system built from the tree its scatter row stored (GEN.315).

A scatter row with companions (`scatter_companions`) stores a member tree: the placed object is member 0, the others
are drawn from their own seeds and birth masses, and the orbits between them are nodes. When a sector's fill (or the
uncharted-system route) builds the row into a system, that tree is the system's stars and orbits, and the build draws no
companion of its own: `plan` turns the stored tree into what `StarSystem(member_tree=...)` takes.

The pair classes (`binary.BinaryStarProxy`, `wide_binary.WideBinaryPair`) hold exactly two ordinary stars, so a tree
builds today when it is one star, or two stars joined by one `close` or `wide` node. Any other tree (three or more
members, a neutron star, a black hole or a brown dwarf among them) builds the placed object alone and is marked
`complete=False`, its stored rows kept for the N-star system that retires the pair classes (GEN.316). A system never
draws a companion of its own for a scattered row, whichever way its tree builds.

The stored phases hold at the time the tree was drawn; `plan` moves the node on to the epoch the built system holds its
orbits at, so the node the system stores is the same orbit at that instant.
"""

import dataclasses

from planetgen.physics import constants
from planetgen.physics.kepler import wrapped_turn_deg

PAIR_STAR_KINDS = ("main-sequence", "pre-main-sequence", "protostar", "giant", "white-dwarf")
"""tuple: The member kinds `Star` builds from a birth mass and an age, which the pair classes can hold."""

PAIR_NODE_KINDS = ("close", "wide")
"""tuple: The node kinds the pair classes stand for (a P-type pair is `close`, an S-type pair `wide`)."""

MASS_TOLERANCE = 1.0e-6
"""float: How far a built companion's mass may differ from its stored mass, as a share of it. A pair the evolution
changed in a way the stored mass records (mass transfer, a merger) does not build from a birth mass and an age."""


@dataclasses.dataclass(frozen=True)
class StoredPlan:
    """
    What a build takes from a stored tree. `companion` is the stored second member and `node` the orbit between the two
    (phases held at the build's epoch), or both `None` for a single star. `complete` says every stored member is in
    the system built, so the row's stored members and nodes can go; it is `False` for a tree the pair classes cannot
    hold, which builds the placed object alone.
    """
    companion: object = None
    node: object = None
    complete: bool = True


SINGLE = StoredPlan()
"""StoredPlan: A scattered row with no members: one star, and no companion drawn."""


def node_at(node, epoch_unix):
    """`node` with its mean anomaly carried on from its own epoch to `epoch_unix`; the same orbit at that instant."""
    if epoch_unix is None or epoch_unix == node.epoch_unix:
        return node
    turn = wrapped_turn_deg(epoch_unix - node.epoch_unix, node.period_years * constants.SECONDS_PER_YEAR) \
        if node.period_years > 0 else 0.0
    return dataclasses.replace(node, mean_anomaly_deg=(node.mean_anomaly_deg + turn) % 360.0,
                               epoch_unix=float(epoch_unix))


def plan(tree, epoch_unix=None):
    """
    The `StoredPlan` of a row's stored tree (`db.scatter_systems.StoredTree`, or `None` for a row with no members),
    its node carried to `epoch_unix` (the orbit epoch the built system holds its phases at; `None`: the node's own).
    """
    if tree is None or not tree.members:
        return SINGLE
    if len(tree.members) == 2 and len(tree.nodes) == 1 and tree.members[1].kind in PAIR_STAR_KINDS \
            and tree.nodes[0].kind in PAIR_NODE_KINDS:
        return StoredPlan(tree.members[1], node_at(tree.nodes[0], epoch_unix), True)
    return StoredPlan(complete=False)


def matches(star, member):
    """Whether the built `star` is the stored `member`: its mass is the stored mass."""
    return abs(star.mass / constants.SOLAR_MASS_TO_KG - member.mass_solar) <= MASS_TOLERANCE * member.mass_solar
