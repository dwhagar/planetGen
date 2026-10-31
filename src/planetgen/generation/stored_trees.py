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
black hole or a brown dwarf among the members) is still the system's tree (GEN.1013): those members are `dark` (`PlannedStar.built`
false), they hang in the nodes with their mass and radius so the orbits and the gravity are the drawn tree's, and the stars
around them are built. The pair is two built stars, so a tree with fewer than two of them builds the one star in its
tree.

The stored phases hold at the time the tree was drawn; `plan` moves the nodes on to the epoch the built system holds its
orbits at, so the nodes the system stores are the same orbits at that instant. A drawn tree holds its phases at
`GENERATED_EPOCH_UNIX`, which the store replaces with the database's orbit epoch when it saves the system.
"""

import collections
import dataclasses

from planetgen.physics.kepler import wrapped_turn_deg
from planetgen.physics.system_tree import OrbitNode
from planetgen.physics import unit_registry as units

_kg_to_msun = units.converter("kg", "Msun")
_km_to_au = units.converter("km", "au")

_yr_to_s = units.converter("yr", "s")

PAIR_STAR_KINDS = ("main-sequence", "pre-main-sequence", "protostar", "giant", "white-dwarf")
"""tuple: The member kinds `Star` builds from a birth mass and an age, which the pair classes can hold."""

UNBUILT_KINDS = ("brown-dwarf", "neutron-star", "black-hole", "stripped-helium")
"""tuple: The member kinds of a drawn tree that `Star` cannot be built as: there is no brown-dwarf host class, a compact
remnant is a `phenomena.compact_remnant` object, not a star of a pair, and a stripped helium star has no record a star is made
from. They stay in the tree as dark members (GEN.1013): mass and radius, no light and no planets."""

PAIR_NODE_KINDS = ("close", "wide")
"""tuple: The node kinds the pair classes stand for (a P-type pair is `close`, an S-type pair `wide`)."""

RECORD_TOLERANCE = 1.0e-9
"""float: How far a built companion's mass and radius may differ from the member's, as a share of it (GEN.1009). The build makes the
star from the member record the tree drew its orbits from, so the two are the same numbers to rounding; a member that is not
(a pair the evolution changed by mass transfer or a merger does not build from a birth mass and an age) is left out."""


GENERATED_EPOCH_UNIX = 0.0
"""float: The epoch a tree drawn at generation holds its phases at, until the store saves the system and puts the database's
orbit epoch there (no real orbit epoch is the start of 1970)."""


@dataclasses.dataclass(frozen=True)
class PlannedStar:
    """One member of a plan: where it hangs in the plan's nodes (`parent_node_no`, `slot`) and what builds it (`seed`,
    `initial_mass_solar`; `mass_solar` and `radius_au` are the mass and radius the tree gave it, which the star built from the
    record has to have too; `age_yr` is the age the record is made at, None for the primary's)."""
    member_no: int
    parent_node_no: int
    slot: int
    kind: str
    mass_solar: float
    initial_mass_solar: float
    seed: int
    radius_au: float | None = None
    age_yr: float | None = None

    @property
    def built(self):
        """Whether `Star` builds this member; a dark member (`UNBUILT_KINDS`) only hangs in the tree."""
        return self.kind not in UNBUILT_KINDS

    def dark(self):
        """This member as the `DarkMember` a built system keeps of it."""
        return DarkMember(self.member_no, self.parent_node_no, self.slot, self.kind, self.mass_solar, self.radius_au)


@dataclasses.dataclass(frozen=True)
class DarkMember:
    """A member of a built system's tree that is not a star of it (GEN.1013): a brown dwarf, a neutron star, a black hole
    (`UNBUILT_KINDS`). It hangs in `parent_node_no`, `slot` with the `mass_solar` and `radius_au` the tree drew; the system
    keeps no more of it, no light and no planets."""
    member_no: int
    parent_node_no: int
    slot: int
    kind: str
    mass_solar: float
    radius_au: float | None = None


@dataclasses.dataclass(frozen=True)
class StoredPlan:
    """
    What a build takes from a tree. `members` are all the stars of the tree in member order, the placed object (member 0,
    already built) first, and `nodes` its orbits, numbered from 0 in the tree's order with the phases held at the build's
    epoch; both are empty for a single star. `complete` says every stored member is in the system built, so the row's
    stored members and nodes can go; it is `False` for a tree the build cannot hold, which builds the placed object alone.
    A built star must have the mass and radius the tree gave its member (`matches`): one that does not is a change the birth mass
    and age do not give.
    """
    members: tuple = ()
    nodes: tuple = ()
    complete: bool = True

    @property
    def lit_members(self):
        """The members `Star` builds, in member order."""
        return tuple(member for member in self.members if member.built)

    @property
    def dark_members(self):
        """The members that only hang in the tree (`UNBUILT_KINDS`) as `DarkMember`s, in member order."""
        return tuple(member.dark() for member in self.members if not member.built)

    def pair_node_no(self):
        """
        The node of the system's pair (module docstring), `None` for fewer than two built stars: the node of the placed star
        when its children are both built stars, else the node of two built stars that holds the most mass. A tree with no
        such node (every pair of built stars has a dark member or a sub-pair between them) has the pair of its two heaviest
        built stars, which is the lowest node above both.
        """
        lit = self.lit_members
        if len(lit) < 2:
            return None
        counts = collections.Counter(member.parent_node_no for member in lit)
        placed_at = self.members[0].parent_node_no
        if counts[placed_at] == 2:
            return placed_at
        mass = collections.Counter()
        for member in self.members:
            mass[member.parent_node_no] += member.mass_solar
        pairs = [node_no for node_no, count in counts.items() if count == 2]
        if pairs:
            return max(pairs, key=lambda node_no: mass[node_no])
        first, second = sorted(lit, key=lambda member: -member.mass_solar)[:2]
        parent = {node.node_no: node.parent_node_no for node in self.nodes}
        above = set()
        node_no = first.parent_node_no
        while node_no is not None:
            above.add(node_no)
            node_no = parent[node_no]
        node_no = second.parent_node_no
        while node_no not in above:
            node_no = parent[node_no]
        return node_no

    def pair_members(self):
        """`(low, high)`, the member numbers of the pair, by slot of the pair node (`None` for no pair): the two built children
        of the pair node, else (a pair with a dark member or a sub-pair between them) the two heaviest built stars by the
        slot of the pair node's child each hangs under."""
        pair_no = self.pair_node_no()
        if pair_no is None:
            return None
        lit = self.lit_members
        direct = [member for member in lit if member.parent_node_no == pair_no]
        if len(direct) == 2:
            return tuple(member.member_no for member in sorted(direct, key=lambda member: member.slot))
        first, second = sorted(lit, key=lambda member: -member.mass_solar)[:2]
        return tuple(member.member_no for member in sorted((first, second), key=lambda m: self.slot_under(m, pair_no)))

    def slot_under(self, member, node_no):
        """The slot of node `node_no` that `member` hangs under, directly or through the nodes below it."""
        parent = {node.node_no: (node.parent_node_no, node.slot) for node in self.nodes}
        at, slot = member.parent_node_no, member.slot
        while at != node_no:
            at, slot = parent[at]
        return slot


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


def _planned(members, nodes, node_for, age_for):
    """The plan of `members` and `nodes` (anything with the fields of `PlannedStar`, and `OrbitNode`s), the node numbers
    taken from 0 in node order and each node made by `node_for`; `age_for(member)` is the age its record is made at."""
    ordered = sorted(nodes, key=lambda node: node.node_no)
    number = {node.node_no: index for index, node in enumerate(ordered)}
    planned_nodes = tuple(dataclasses.replace(node_for(node), node_no=number[node.node_no],
                                              parent_node_no=number.get(node.parent_node_no)) for node in ordered)
    planned_members = tuple(PlannedStar(m.member_no, number[m.parent_node_no], m.slot, m.kind, m.mass_solar,
                                        m.initial_mass_solar, m.seed, m.radius_au, age_for(m))
                            for m in sorted(members, key=lambda m: m.member_no))
    return StoredPlan(planned_members, planned_nodes, True)


def plan(tree, epoch_unix=None):
    """
    The `StoredPlan` of a row's stored tree (`db.scatter_systems.StoredTree`, or `None` for a row with no members), its
    nodes carried to `epoch_unix` (the orbit epoch the built system holds its phases at; `None`: the nodes' own).
    """
    if tree is None or not tree.members:
        return SINGLE
    if len(tree.members) < 2 or any(member.kind not in PAIR_STAR_KINDS + UNBUILT_KINDS for member in tree.members[1:]) \
            or any(node.kind not in PAIR_NODE_KINDS for node in tree.nodes):
        return StoredPlan(complete=False)
    return _planned(tree.members, tree.nodes, lambda node: node_at(node, epoch_unix), lambda member: None)


def drawn(tree):
    """The `StoredPlan` of a `multistar.MemberTree` just drawn, dark members and all (GEN.1013): its phases held at
    `GENERATED_EPOCH_UNIX`, each member's record made at the age the tree drew it at (a merger product's rejuvenated one)."""
    if tree.member_count < 2:
        return SINGLE
    return _planned(tree.members, tree.nodes, lambda node: OrbitNode(
        node_no=node.node_no, parent_node_no=node.parent_node_no, slot=node.slot, mass_solar=node.mass_solar,
        separation_km=node.separation_km, eccentricity=node.eccentricity, inclination_deg=node.inclination_deg,
        ascending_node_deg=node.ascending_node_deg, arg_periapsis_deg=node.arg_periapsis_deg,
        mean_anomaly_deg=node.mean_anomaly_deg, epoch_unix=GENERATED_EPOCH_UNIX, period_years=node.period_years,
        kind=node.kind, interaction=node.interaction, flags=node.flags),
        lambda member: member.effective_age_yr or tree.age_yr)


def matches(star, member):
    """Whether the built `star` is the member the tree drew: its mass and radius are the member's (GEN.1009), which they are
    when it was made from the member's record. A member with no stored radius (a row made before v184) is not."""
    if member.radius_au is None:
        return False
    radius_au = _km_to_au(star.radius)
    return (abs(_kg_to_msun(star.mass) - member.mass_solar) <= RECORD_TOLERANCE * member.mass_solar
            and abs(radius_au - member.radius_au) <= RECORD_TOLERANCE * member.radius_au)
