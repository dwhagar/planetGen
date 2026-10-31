# planetgen/generation/multistar.py

"""
The member-tree generator (GEN.261): the stars, brown dwarfs, remnants and
protostars that share a system, in a stable hierarchy of up to twelve
members, drawn for one placed object (a scattered star or compact object, or
any object the fill creates). One function, `draw_member_tree`, with two call
sites: the scatter passes and the fill. It draws the members and their orbits
only: planets, moons and belts stay at sector fill.

The laws and numbers are in docs/design/systems/multi-star-systems.md (sections 4, 6
and 13) and docs/design/systems/multistar-and-compact-systems.md (sections 3, 4, 5
and 7); the stability criteria, the Kozai-Lidov screen and the Peters time are
`planetgen.physics.stability`. In short:

1. **How many stars.** The `star-count` table (`phenomena/hosts.py`, GEN.552) by the placed star's
   spectral class (O to M, from its initial mass), a tail for 8 to 12 stars
   that continues the table's last column, and the `prevalence` knob on the
   five-or-more columns. A forced `n_stars` or a lower `max_stars` wins.
2. **Which members.** The other stars have initial masses `m1 * U(0.1, 1)`
   (the ratio no lower than keeps the star over the hydrogen-burning limit,
   GEN.899: the table counts stars) and take their present state from the
   system's age (main sequence, giant, white dwarf, pre-main-sequence,
   protostar). A brown dwarf is one more member, in `tuning.BROWN_DWARF_COMPANION_CHANCE`
   of the trees a directive does not count. A compact companion (white dwarf, neutron star, black hole)
   is added at the odds of `MEMBER_STATE_BY_PLACED`, one uniform draw per
   placed object, so the scatter's bulk pass can screen with it
   (`companion_state`) and call the tree only for the hits.
3. **The shape.** Random splitting of the members into a binary tree
   (2+1; 3+1 at 60% or 2+2 at 40%; k = 2 weighted 1.5 for five or more). A
   compact companion starts as a pair with the placed object.
4. **The orbits, top down.** Each pair gets a period from its law (the
   truncated log-normal for stars, a log-uniform one for a compact object
   with a star and for a compact pair), a thermal eccentricity (circular
   under 7 days) and an isotropic orientation. A pair inside a larger one has
   its axis bounded by `stability.max_inner_axis` and is then checked with
   `stability.is_triple_stable` and the Kozai screen; up to 60 tries a
   pair, 30 restarts a tree, then the lightest star is dropped.
5. **The compact filters.** A compact pair whose Peters time is shorter than
   the age is rejected; a bound companion beyond 1,000 AU of a neutron star
   or a black hole that did not collapse directly is rejected; at most three
   compact members to a tree.

Every draw comes from `rng` (a `draw.Stream` or the `draw` module) in a fixed
order, so the same stream gives the same tree. Orbit angles are in degrees,
separations in AU (and km), periods in years and times in years.
"""

import dataclasses
import math
from dataclasses import dataclass
from statistics import NormalDist

from planetgen import tuning
from planetgen.generation import member_record
from planetgen.generation import close_binary_trend
from planetgen.generation.phenomena import hosts, system_types
from planetgen.physics import (
    accretion, brown_dwarf, constants, helium_core_white_dwarf, imf, mass_boundaries, ms_lifetime, neutron_star, pair_history, remnant_mass, remnant_progenitors, roche,
    stability, stellar_abundances, stellar_composition)
from planetgen.physics import stellar_evolution as evolution
from planetgen.util import draw
from planetgen.physics import unit_registry as units

_m_to_au = units.converter("m", "au")
_km_to_m = units.converter("km", "m")
_km_to_au = units.converter("km", "au")
_m_to_km = units.converter("m", "km")
_au_to_m = units.converter("au", "m")
_yr_to_d = units.converter("yr", "d")
_d_to_yr = units.converter("d", "yr")

TAIL_RATIOS = (0.12, 0.10)
"""tuple: The chance a system of 7 stars has an 8th (0.12) and the chance each
further star follows the last (0.10), up to `max_stars` (multi-star-systems.md
13.1). The tail comes out of the 7-star column, so the table's totals stand."""

MAX_MEMBERS = tuning.MAX_STARS_DEFAULT
"""int: The most members in a tree unless a caller lowers or raises it (`SystemConfig.MAX_STARS`, GEN.267)."""

PERIOD_LAWS_LOG10_DAYS = {"solar": (5.03, 2.28), "massive": (3.5, 1.5)}
"""dict: The mean and sigma of the log-normal in log10 of the period in days: the
solar-type law (Raghavan 2010) and, for a primary of `MASSIVE_PRIMARY_MIN_MASS_SOL`
or more, the law that peaks at log P = 3.5 (Moe and Di Stefano 2017); the sigma of
the second is a design default."""

MASSIVE_PRIMARY_MIN_MASS_SOL = 2.1

CLOSE_COMPANION_AXIS_AU = 10.0
"""float: The axis inside which a star's [Fe/H] moves its chance of a companion (GEN.103; Moe, Kratter and Badenes 2019 measured
the trend for periods under 10^4 days, which is 10 AU for two solar masses)."""
ROOT_MAX_AXIS_AU = 20_000.0
"""float: The widest root pair (0.1 pc), a directive raises it."""

ROOT_ECCENTRICITY_CAP = 0.95
INNER_ECCENTRICITY_CAP = 0.8
CIRCULAR_BELOW_DAYS = 7.0
COMPACT_PAIR_AXIS_AU = (0.003, 100.0)
COMPACT_PAIR_CIRCULAR_BELOW_AU = 0.1
INTERACTED_PERIOD_DAYS = (0.3, 1.0e4)
"""tuple: The log-uniform period range of a neutron star or black hole with a star."""

EXTENT_AXIS_FACTOR = 2.5
"""float: A pair's axis is at least this many times its members' summed radii."""

PERICENTRE_EXTENT_FACTOR = 2.0
"""float: The pericentre must clear this many times the members' extent."""

MAX_FILL_FACTOR = 1.2
"""float: How far past its Roche lobe a star may reach at the pair's closest approach (1 is semi-detached,
both past 1 is contact). Pairs of two single members use this in place of `EXTENT_AXIS_FACTOR`."""

STABILITY_MARGIN = 1.15
"""float: How much wider than the bare criterion an inner orbit's outer partner must be, at every level of a tree (GEN.1010, GM15): the
Eggleton and Kiseleva criterion is reliable to about 20 percent, and the nested-triple test over-predicts stable trees
(Vynatheya, Mardling and Hamers 2023), so each triple of the tree has to clear the criterion by this factor. A stand-in."""


@dataclass(frozen=True)
class Side:
    """One side of a pair for the keep-out rule: its `mass` (solar), its `extent` (AU: a single member's radius, a node's greatest
    distance from its own barycentre to the surface of a star below it, `node_extent_au`) and its `radius`, the radius of a single
    member (None for a node)."""
    mass: float
    extent: float
    radius: float | None = None


def single_side(mass, radius_au):
    return Side(mass, radius_au, radius_au)


def node_extent_au(axis_au, eccentricity, low, high):
    """The extent of the node of the pair `low`, `high` (`Side`s): the farthest a star below it reaches from the node's barycentre,
    both sides at the apocentre. The one definition the build, the tree and the verifier (`physics/clearance.py`) share."""
    apocentre = axis_au * (1.0 + eccentricity)
    total = low.mass + high.mass
    return max(apocentre * high.mass / total + low.extent, apocentre * low.mass / total + high.extent)


def pair_floor_au(low, high):
    """The closest the pair may be (GEN.1010): the Roche-lobe limit for two single members, else the extent rule, at every depth of
    the tree."""
    if low.radius is None or high.radius is None:
        return EXTENT_AXIS_FACTOR * (low.extent + high.extent)
    return max(roche.separation_for_filling(low.radius, MAX_FILL_FACTOR, low.mass, high.mass),
               roche.separation_for_filling(high.radius, MAX_FILL_FACTOR, high.mass, low.mass))


def pair_pericentre_floor_au(low, high):
    """The closest the pair's pericentre may be: the Roche-lobe limit for two single members, else `PERICENTRE_EXTENT_FACTOR` times the
    extents (what `pair_clears` asks)."""
    if low.radius is None or high.radius is None:
        return PERICENTRE_EXTENT_FACTOR * (low.extent + high.extent)
    return pair_floor_au(low, high)


def pair_clears(axis_au, eccentricity, low, high):
    """Whether the pair's closest approach keeps its sides inside their limits: neither single member past `MAX_FILL_FACTOR` of its
    Roche lobe, else the pericentre at least `PERICENTRE_EXTENT_FACTOR` times the sides' extents. The one rule for the tree, the
    build (`built_tree_clears`) and, re-derived, the verifier."""
    peri = axis_au * (1.0 - eccentricity)
    if low.radius is None or high.radius is None:
        return peri >= PERICENTRE_EXTENT_FACTOR * (low.extent + high.extent)
    return (roche.filling_factor(low.radius, peri, low.mass, high.mass) <= MAX_FILL_FACTOR
            and roche.filling_factor(high.radius, peri, high.mass, low.mass) <= MAX_FILL_FACTOR)

ACCRETOR_KINDS = ("white-dwarf", "neutron-star", "black-hole")
"""tuple: The member kinds that receive mass from a lobe-filling partner (an `accreting` pair)."""

SUDDEN_MASS_LOSS_KINDS = ("neutron-star", "black-hole")
"""tuple: The member kinds whose mass left in a supernova, not slowly (`pair_history.Member.sudden`)."""

INTERACTIONS = ("detached", "semi-detached", "contact", "accreting", "post-common-envelope", "merger-product")
"""tuple: The `Node.interaction` values; the classifier gives the first four and `post-common-envelope`
from geometry, `merger-product` comes from the merger draw of the close-binary types."""

REJUVENATION_FRACTION_RANGE = (0.1, 0.6)
"""tuple: A merger product's effective age as a fraction of the system's age (a design default)."""

FRAGILE_AXIS_AU = 10_000.0
"""float: A root pair wider than this is ultra-wide and flagged `fragile` (passing stars unbind it)."""

PAIR_RETRIES = 60
"""int: How many times a pair's draw may be repeated, for the screens that depend on the draw itself and for an orientation that leaves
no room (GEN.1011). The period and eccentricity are drawn once from their cut laws; no draw is repeated for overlap or stability."""
TREE_RESTARTS = 30
MAX_COMPACT_MEMBERS = 3
WIDE_COMPACT_AXIS_AU = 1000.0
"""float: A bound companion beyond this of a neutron star or a black hole that did
not collapse directly does not survive the collapse."""

DIRECT_COLLAPSE_HOLE_SOL = 20.0
"""float: A black hole counts as a direct collapse when its progenitor failed to explode (a black hole's progenitor
lies in a window of `remnant_progenitors`, which is every one the generator draws), or, when the progenitor is not known,
when the hole is `DIRECT_COLLAPSE_HOLE_SOL` or more (a stand-in; the note leaves it open)."""

COMPACT_TERTIARY_SHARE = 0.10
"""float: A compact binary has a third member in the share of the `compact-triple` row (30%, GEN.272; before it, black holes
only), which is a compact object in 10% of those (a design default)."""

BROWN_DWARF_DESERT_AU = 5.0
BROWN_DWARF_DESERT_KEEP = 0.01
"""A brown dwarf inside `BROWN_DWARF_DESERT_AU` of a star is kept in 1% of draws."""

REMNANT_WINDOW_FRACTION = 1.0 / 3.0
"""float: A remnant companion of a living star exists only in the last third of the star's
life, so the table value (a mean over all ages) is divided by this inside the window."""

COMPACT_KINDS = ("neutron-star", "black-hole")
MEMBER_KINDS = ("main-sequence", "giant", "white-dwarf", "neutron-star", "black-hole", "brown-dwarf",
                "pre-main-sequence", "protostar", "stripped-helium", "thorne-zytkow", "wolf-rayet")
"""tuple: Every member kind. The last three exist for the type generators that follow (the evolution of a
Wolf-Rayet star, a stripped helium star and a Thorne-Zytkow object is theirs); `star_labels` names them."""
COMPANION_STATES = ("star", "black-widow", "redback")
"""tuple: The states of `MEMBER_STATE_BY_PLACED` that are not member kinds: the companion is a star of the tree's laws, or
one of a pulsar's two spiders (GEN.272)."""
PLACED_KINDS = ("star", "black-hole", "intermediate-mass-black-hole", "neutron-star", "white-dwarf")

MEMBER_DIRECT_COLLAPSE = 1
MEMBER_MASS_REVERSED = 2
MEMBER_MASS_GAINER = 4
MEMBER_MERGER_PRODUCT = 8
MEMBER_BE = 16
MEMBER_FRAGILE = 32
MEMBER_THEORETICAL = 64
MEMBER_RUNAWAY_SUPERNOVA = 128
MEMBER_RUNAWAY_DYNAMICAL = 256
MEMBER_R_CRB = 512
MEMBER_FK_COMAE = 1024
MEMBER_MAGNETIC = 2048
MEMBER_BARIUM = 4096
MEMBER_SYMBIOTIC = 8192
MEMBER_MICROQUASAR = 16384
MEMBER_ABLATING = 32768
MEMBER_FLAGS = {"direct_collapse": MEMBER_DIRECT_COLLAPSE, "mass_reversed": MEMBER_MASS_REVERSED,
                "mass_gainer": MEMBER_MASS_GAINER, "merger_product": MEMBER_MERGER_PRODUCT, "be": MEMBER_BE,
                "fragile": MEMBER_FRAGILE, "theoretical": MEMBER_THEORETICAL,
                "runaway_supernova": MEMBER_RUNAWAY_SUPERNOVA, "runaway_dynamical": MEMBER_RUNAWAY_DYNAMICAL,
                "r_crb": MEMBER_R_CRB, "fk_comae": MEMBER_FK_COMAE, "magnetic": MEMBER_MAGNETIC,
                "barium": MEMBER_BARIUM, "symbiotic": MEMBER_SYMBIOTIC, "microquasar": MEMBER_MICROQUASAR,
                "ablating": MEMBER_ABLATING}
"""dict: The bits of `Member.flags`. The classifier sets `mass_gainer` (the star that receives the mass of a
semi-detached or contact pair), `mass_reversed` (the donor of a semi-detached pair that was the heavier at birth: an Algol),
`fragile` and `theoretical`; the close-binary types (GEN.270) set `merger_product` (a massive star that merged, a field blue
straggler, an R Coronae Borealis or FK Comae star), `be`, `mass_gainer` (a massive star that has gained mass), `r_crb`
and `fk_comae`; and the runaway origin (none, or one of the two channels, `draw_runaway_origin`, GEN.329) is the two
`runaway_*` bits: a built system's runaway star carries the one its draw named in `stars.flags`. All sixteen bits of the
column are taken, so a member of a young group has no flag of its own (GEN.327): its node's kind `group` and its
`decay_age_yr` say it has not decayed yet."""

RUNAWAY_ORIGIN_FLAGS = {"supernova": MEMBER_RUNAWAY_SUPERNOVA, "dynamical": MEMBER_RUNAWAY_DYNAMICAL}
"""dict: A runaway's origin -> the flag bit that records it."""
NODE_INTERACTED = 1
NODE_KL_ACTIVE = 2
NODE_WIND_FED = 4
NODE_DISC_FED = 8
"""int: The node bits of an `accreting` pair of a neutron star or a black hole (GEN.272): the gas is a hot star's wind
(`NODE_WIND_FED`) or a Be star's decretion disc (`NODE_DISC_FED`); neither is a donor that fills its lobe. They are what
`MemberTree.x_ray_luminosities` reads."""

MEMBER_STATE_BY_PLACED = {
    "star-massive": {"neutron-star": (0.002, 0.00033, 0.0133), "black-hole": (0.001, 0.00017, 0.0067)},
    "star-intermediate": {"white-dwarf": (0.01, 0.002, 0.05), "neutron-star": (0.0003, 0.00005, 0.0015),
                          "black-hole": (0.0002, 0.00005, 0.0015)},
    "black-hole": {"star": (0.010, 0.001, 0.022), "black-hole": (0.002, 0.00017, 0.024),
                   "neutron-star": (0.0005, 0.000075, 0.0007)},
    "intermediate-mass-black-hole": {"star": (0.01, 0.001, 0.05)},
    "neutron-star": {"star": (0.02, 0.01, 0.03), "neutron-star": (1.8e-4, 6e-5, 4.3e-4),
                     "black-hole": (2.5e-5, 7e-6, 7e-5)},
    "recycled-pulsar": {"white-dwarf": (0.661, 0.55, 0.72), "black-widow": (0.0624, 0.03, 0.09),
                        "redback": (0.0312, 0.015, 0.05), "star": (0.0234, 0.01, 0.04),
                        "neutron-star": (2.0e-3, 6e-4, 4.3e-3), "black-hole": (2.5e-5, 7e-6, 7e-5)},
    "white-dwarf": {"star": (0.205, 0.15, 0.22), "white-dwarf": (0.095, 0.05, 0.12)},
    "planetary-nebula-star": {"star": (0.15, 0.10, 0.20), "white-dwarf": (0.05, 0.03, 0.10)},
}
"""dict: Per class of placed object, the chance of a companion of each state as
`(default, low, high)`, the mean over all ages and per object
(multi-star-systems.md 4.2 to 4.4). A star class (`star-massive` for an initial mass of 14
solar masses or more, `star-intermediate` for 2 to 14) adds a compact companion to its
tree of stars, inside the last third of its life; the other classes decide whether the
placed object has a companion at all, and which one. Other placed stars get no compact
companion here. The states a compact or white dwarf placed object can have are `star`, `brown-dwarf`, `white-dwarf`,
`neutron-star`, `black-hole`, and (GEN.272) `black-widow` and `redback`, a recycled pulsar's spider companions, which ablate
it. The low-mass X-ray binaries are part of a neutron star's and a black hole's `star` odds (`draw_low_mass_xray_pair`)."""


@dataclass(frozen=True)
class Placed:
    """The object a tree is drawn for. `mass_solar` is its present mass; `initial_mass_solar` its birth mass
    (defaults to the present mass). `recycled` marks a millisecond pulsar; `planetary_nebula` a
    planetary-nebula central star (a white dwarf). `fe_h` is the placed star's birth [Fe/H], when known. `record` is the
    `member_record.MemberRecord` of a placed star that is already built (a system's primary): the tree reads its radius
    and mass, so it clears the star the build has; None for a scattered object, whose tree is drawn from its birth mass."""
    kind: str
    mass_solar: float
    initial_mass_solar: float | None = None
    recycled: bool = False
    planetary_nebula: bool = False
    fe_h: float | None = None
    record: "member_record.MemberRecord | None" = None

    @property
    def birth_mass_solar(self):
        return self.mass_solar if self.initial_mass_solar is None else self.initial_mass_solar

    @property
    def tree_fe_h(self):
        """The birth [Fe/H] the tree's stars are evolved at (GEN.651): the placed star's own, which its companions
        share (a system is born from one cloud), else the reference composition, 0 (a remnant's row stores none)."""
        return 0.0 if self.fe_h is None else self.fe_h


@dataclass(frozen=True)
class Member:
    """One member of a tree: member 0 is the placed object. `radius_au` is its extent; `flags` are the bits of
    `MEMBER_FLAGS`; `effective_age_yr` is a merger product's rejuvenated age (None: the system's age). `state` is the
    Yerkes class of a star (`V`, `III`, `VII`, ...), `PMS` or `protostar` for a young one, the spectral type of a brown
    dwarf (`M8`, `L5`, `T2`, `Y0`; `physics/brown_dwarf.py`, GEN.252), and None for a remnant; it is stored in
    `scatter_members.state`, a VARCHAR(16). A member of a `group` node carries its `offset_au` from the group's centre, its
    `velocity_kms` and the group's `decay_age_yr` (GEN.327); every other member has none of the three. `radius_max_au` and
    `mass_max_solar` are the largest radius and mass the member has had (GEN.1009, `member_record`); a tree draws them for
    every member."""
    member_no: int
    parent_node_no: int | None
    slot: int
    kind: str
    state: str | None
    flags: int
    mass_solar: float
    initial_mass_solar: float
    radius_au: float
    seed: int
    effective_age_yr: float | None = None
    offset_au: tuple | None = None
    velocity_kms: tuple | None = None
    decay_age_yr: float | None = None
    radius_max_au: float | None = None
    mass_max_solar: float | None = None


@dataclass(frozen=True)
class Node:
    """One orbit: the relative orbit of the child in `slot` 1 about the child in `slot` 0. Node 1 is the
    root. `kind` is `close`, `wide` or `compact`; `interaction` is `detached` until the Roche classifier
    (GEN.269) sets it; `flags` bit 0 is an interacted pair (the tight law), bit 1 `kl_active`."""
    node_no: int
    parent_node_no: int | None
    slot: int
    mass_solar: float
    separation_au: float
    separation_km: float
    eccentricity: float
    inclination_deg: float
    ascending_node_deg: float
    arg_periapsis_deg: float
    mean_anomaly_deg: float
    period_years: float
    kind: str
    interaction: str
    flags: int


@dataclass(frozen=True)
class MemberTree:
    """The result: `members` (placed object first), `nodes` (empty for a single object) and `dropped`, the
    number of members the generator could not fit and left out. `ejected` holds the member a decayed young group threw out
    (GEN.327): the lightest, which left the system and is in none of `members`."""
    members: tuple
    nodes: tuple
    age_yr: float
    dropped: int = 0
    ejected: tuple = ()

    @property
    def member_count(self):
        return len(self.members)

    @property
    def system_mass_solar(self):
        return sum(member.mass_solar for member in self.members)

    @property
    def star_count(self):
        return sum(member.kind not in COMPACT_KINDS and member.kind != "brown-dwarf" for member in self.members)

    @property
    def compact_count(self):
        return sum(member.kind in COMPACT_KINDS for member in self.members)

    def eclipsing_nodes(self):
        """The node numbers of the pairs that eclipse (GEN.270; derived, never stored): a node whose two children are
        single members and whose inclination is within `asin((R1 + R2) / a)` of edge-on (`pair_eclipses`)."""
        found = []
        for node in self.nodes:
            pair = [member for member in self.members if member.parent_node_no == node.node_no]
            if len(pair) == 2 and pair_eclipses(node.separation_au, node.inclination_deg, pair[0].radius_au,
                                                pair[1].radius_au):
                found.append(node.node_no)
        return tuple(found)

    def spectroscopic_nodes(self):
        """The node numbers of the pairs closer than `tuning.SPECTROSCOPIC_MAX_AXIS_AU` (derived, never stored): the
        orbits whose radial-velocity swing a spectrograph resolves."""
        return tuple(node.node_no for node in self.nodes if node.separation_au < tuning.SPECTROSCOPIC_MAX_AXIS_AU)

    def _pair_members(self, node):
        """The two members of a node whose children are single members, else None."""
        pair = [member for member in self.members if member.parent_node_no == node.node_no]
        return pair if len(pair) == 2 else None

    def x_ray_luminosities(self):
        """
        `{node_no: watts}` for every `accreting` pair of a neutron star or a black hole and a single member (GEN.272; derived,
        never stored): `eta Mdot c^2` at most the accretor's Eddington limit (`physics.accretion.x_ray_luminosity_w`). The
        gas is a wind (`NODE_WIND_FED`), a Be star's disc (`NODE_DISC_FED`) or a donor that fills its lobe. A white dwarf's
        accretion is not X-ray and is not listed.
        """
        found = {}
        for node in self.nodes:
            pair = self._pair_members(node) if node.interaction == "accreting" else None
            accretor = None if pair is None else next((m for m in pair if m.kind in COMPACT_KINDS), None)
            if accretor is None:
                continue
            donor = pair[0] if accretor is pair[1] else pair[1]
            mode = "wind" if node.flags & NODE_WIND_FED else "disc" if node.flags & NODE_DISC_FED else "overflow"
            found[node.node_no] = accretion.x_ray_luminosity_w(
                mode, accretor.kind, accretor.mass_solar, accretor.radius_au, donor.kind, donor.mass_solar,
                donor.radius_au, node.separation_au, node.eccentricity)
        return found

    def be_xray_nodes(self):
        """The node numbers of the Be/X-ray binaries (GEN.272; derived): an `accreting` pair fed by a disc, whose donor is
        flagged `be`."""
        return tuple(node.node_no for node in self.nodes
                     if node.interaction == "accreting" and node.flags & NODE_DISC_FED)

    def spider_nodes(self):
        """`{node_no: 'black-widow' or 'redback'}` for the pairs of an ablating pulsar and its companion (GEN.272; derived):
        a companion under 0.1 solar masses is a black widow, a heavier one a redback (the mass gap is empty)."""
        found = {}
        for node in self.nodes:
            pair = self._pair_members(node)
            if pair is None or not any(m.flags & MEMBER_ABLATING for m in pair):
                continue
            companion = next(m for m in pair if not m.flags & MEMBER_ABLATING)
            found[node.node_no] = "black-widow" if companion.mass_solar < 0.1 else "redback"
        return found

    def dormant_nodes(self):
        """The node numbers of the dormant neutron stars and black holes with a star (GEN.272, the Gaia kind; derived): a
        `detached` pair of a neutron star or a black hole and a main-sequence star, a subgiant or a giant."""
        found = []
        for node in self.nodes:
            pair = self._pair_members(node) if node.interaction == "detached" else None
            if pair is None:
                continue
            kinds = {m.kind for m in pair}
            if kinds & set(COMPACT_KINDS) and kinds & {"main-sequence", "giant"} and not kinds - set(COMPACT_KINDS) - {
                    "main-sequence", "giant"}:
                found.append(node.node_no)
        return tuple(found)


def pair_eclipses(separation_au, inclination_deg, radius_a_au, radius_b_au):
    """Whether a pair of these radii (AU) on an orbit of this semi-major axis (AU) and inclination (degrees, 90 is
    edge-on) eclipses: the inclination lies within `asin((R1 + R2) / a)` of edge-on (a pair whose radii reach its axis
    always does)."""
    ratio = (radius_a_au + radius_b_au) / separation_au
    return ratio >= 1.0 or abs(90.0 - inclination_deg) <= math.degrees(math.asin(ratio))


@dataclass
class _Leaf:
    kind: str
    state: str | None
    mass: float
    initial: float
    radius_au: float
    flags: int = 0
    effective_age_yr: float | None = None
    orbit_au: tuple | None = None
    """tuple or None: `(low, high)` AU. A pair of two single members, one of which has this, is placed with an axis
    log-uniform over it (GEN.271: the close double white dwarf, the post-common-envelope pair, the cataclysmic
    variable, whose axis is the one that fills the donor's lobe), in place of the period law."""
    accretes: str | None = None
    """str or None: `wind`, `disc` or `overflow` on the neutron star or black hole of a pair that accretes (GEN.272): the
    pair is `accreting` whatever its geometry says, and a pair of a neutron star or black hole with no such leaf never is."""
    radius_max_au: float | None = None
    mass_max_solar: float | None = None
    """float or None: The largest radius and mass the member has had (`member_record.MemberRecord`); None for a leaf with no
    larger past than its present, which reads as the present radius and the larger of its birth and present mass."""
    seed: int | None = None
    """int or None: The seed the member's record was made under, which its `Member` keeps; None draws one from the tree's source."""
    circular: bool = False
    """bool: The pair this member is in has a circular orbit and is `post-common-envelope`: the companion of a helium-core white
    dwarf the star draw made (GEN.996, `helium_core_draw`), on the axis its giant filled its lobe at."""

    @property
    def largest_radius_au(self):
        return self.radius_au if self.radius_max_au is None else max(self.radius_max_au, self.radius_au)

    @property
    def largest_mass_solar(self):
        return max(self.initial, self.mass) if self.mass_max_solar is None else max(self.mass_max_solar, self.mass)

    def member_seed(self, rng):
        return rng.getrandbits(63) if self.seed is None else self.seed

    @property
    def compact(self):
        return self.kind in COMPACT_KINDS

    @property
    def direct_collapse(self):
        return bool(self.flags & MEMBER_DIRECT_COLLAPSE)


@dataclass(frozen=True)
class RunawayDraw:
    """A runaway star's draw (GEN.329): its `origin` (`supernova`: the companion exploded and the survivor was flung out;
    `dynamical`: a close encounter in a young group) and whether a compact companion stays bound to it (a supernova
    that left the pair bound)."""
    origin: str
    bound_companion: bool


def runaway_rate(initial_mass_sol):
    """The share of stars of this initial mass that are runaways where the catalogue has an observed one: 10% to 30% of
    O stars (default 20%) and 5% to 10% of B stars (default 7.5%, Gies 1987), else None (the field rate of
    `tuning.PHENOMENON_DENSITY_PC3['runaway-star']` applies at the fill)."""
    entry = system_types.SYSTEM_TYPE_RATES["runaway-star"]
    letter = spectral_class(initial_mass_sol)
    return entry.default if letter == "O" else entry.by_class["B"] if letter == "B" else None


def runaway_origin_at(unit):
    """The channel the uniform `unit` (in [0, 1)) names for a runaway: `supernova` for the first two thirds of the
    range, `dynamical` for the rest (two thirds and one third of the massive runaways, Dray et al.)."""
    return "supernova" if unit < system_types.SYSTEM_TYPE_RATES["runaway-star"].by_class["supernova channel"] \
        else "dynamical"


def bound_companion_share():
    """The share of supernova-channel runaways that keep a compact companion: the collapse leaves a binary bound in the
    17% that the 83% of non-interacting massive binaries unbound by the first collapse do not take
    (`unbound-by-collapse`); the unbound ones are not stored, the survivor is an ordinary row."""
    return 1.0 - system_types.SYSTEM_TYPE_RATES["unbound-by-collapse"].default


def draw_white_dwarf_pollution(teff_k, key):
    """
    Whether a white dwarf is polluted by accreted rocky debris (GEN.253, the `white-dwarf-pollution` row): a cool one is at
    the registry's default share (25 to 50 percent, the middle 0.375), a hot one never. The draw is the star's own
    (`stellar_composition.key_draws(key).u_pollution`, `key` from `stellar_composition.abundance_key`), so nothing is stored
    and the star's makeup (`star_makeup(...).detail["polluted"]`) says the same.
    """
    return stellar_composition.white_dwarf_polluted(
        teff_k, stellar_composition.key_draws(key).u_pollution,
        system_types.SYSTEM_TYPE_RATES["white-dwarf-pollution"].default)


def draw_runaway_origin(initial_mass_sol, runaway_unit, origin_unit, bound_unit):
    """
    The runaway draw of a star of `initial_mass_sol` (GEN.329), from three uniforms in [0, 1) so the scatter can derive
    them from a row's seed: `None` for a star that is not a runaway (`runaway_unit` at or above `runaway_rate`, or a
    class with no observed rate), else its `RunawayDraw`. Only the supernova channel can leave a compact companion
    bound (`bound_companion_share`); a hypervelocity star is never drawn here (the Hills mechanism flings out one star).
    """
    rate = runaway_rate(initial_mass_sol)
    if rate is None or runaway_unit >= rate:
        return None
    origin = runaway_origin_at(origin_unit)
    return RunawayDraw(origin, origin == "supernova" and bound_unit < bound_companion_share())


def spectral_class(initial_mass_sol):
    """The spectral class (O to M) of a star of this initial mass."""
    for letter, lower in mass_boundaries.SPECTRAL_CLASS_LOWER_MASS_SOL:
        if initial_mass_sol >= lower:
            return letter
    return "M"


def companion_class(initial_mass_sol, letter=None):
    """The key into the `star-count` table (`hosts.STAR_COUNT_PROBABILITIES`) of a star: its
    spectral letter (from the mass when `letter` is None), with the B stars of 5 to 16 solar masses split into B-mid and
    B-early (Moe and Di Stefano's bins, GEN.449)."""
    letter = letter or spectral_class(initial_mass_sol)
    if letter == "B":
        for name, lower in mass_boundaries.B_SUBCLASS_LOWER_MASS_SOL:
            if initial_mass_sol >= lower:
                return name
    return letter


def placed_class(placed):
    """The `MEMBER_STATE_BY_PLACED` key of a placed object, or None when it has no companion table
    (a star of under 2 solar masses)."""
    if placed.kind == "star":
        birth = placed.birth_mass_solar
        return "star-massive" if birth >= 14.0 else "star-intermediate" if birth >= 2.0 else None
    if placed.kind == "neutron-star":
        return "recycled-pulsar" if placed.recycled else "neutron-star"
    if placed.kind == "white-dwarf":
        return "planetary-nebula-star" if placed.planetary_nebula else "white-dwarf"
    if placed.kind in PLACED_KINDS:
        return placed.kind
    raise ValueError(f"unknown placed kind {placed.kind!r}; choose from {PLACED_KINDS}")


def companion_state_at(placed, pick):
    """
    The companion state the uniform `pick` (in [0, 1)) names against the cumulative odds of
    `MEMBER_STATE_BY_PLACED` for a placed object that is not a star: `star`, `white-dwarf`, `brown-dwarf`,
    `neutron-star` or `black-hole`, or None. `companion_state` draws the uniform; the scatter's bulk screen
    (`scatter_companions`) derives it from the row's seed.
    """
    table = MEMBER_STATE_BY_PLACED.get(placed_class(placed), {})
    for state, (default, _low, _high) in table.items():
        if pick < default:
            return state
        pick -= default
    return None


def companion_state(placed, rng=draw):
    """
    One uniform draw against the cumulative odds of `MEMBER_STATE_BY_PLACED` for a placed object that is
    not a star: the companion state it gets (`star`, `white-dwarf`, `brown-dwarf`, `neutron-star` or
    `black-hole`), or None. This is the bulk screen of the scatter (multi-star-systems.md 6.3):
    only a row it names needs `draw_member_tree`.
    """
    return companion_state_at(placed, rng.random())


def same_class_companions(kind, recycled_share=0.0):
    """
    The expected number of companions of the placed object's own kind per placed object of `kind`
    (`black-hole`, `neutron-star`, ...): the odds of a black hole's black-hole companion, a neutron star's
    neutron-star one (young and recycled pulsars weighted by `recycled_share`). The scatter divides a class's
    object rate by `1 + this`, since a row is a system keyed by its placed object and its same-class companions
    are objects of the class too (multi-star-systems.md 6.4).
    """
    def odds(placed_class_name, state):
        return MEMBER_STATE_BY_PLACED.get(placed_class_name, {}).get(state, (0.0,))[0]
    if kind == "neutron-star":
        return (1.0 - recycled_share) * odds("neutron-star", kind) + recycled_share * odds("recycled-pulsar", kind)
    return odds(kind, kind)


def star_count_weights(initial_mass_sol, prevalence=1.0):
    """The weights of 1 to 7 stars for a star of this initial mass (the `star-count` table, `phenomena/hosts.py`), the
    five-or-more columns scaled by `prevalence`."""
    return hosts.star_count_weights(companion_class(initial_mass_sol), prevalence)


def _tilted(u, share, weight):
    """`u`, a uniform draw, moved so the part of the range under `share` of it has the odds of `share` times `weight`: the
    inside part stretches to the new share, the outside part to what is left."""
    if not 0.0 < share < 1.0:
        return u
    odds = share / (1.0 - share) * weight
    moved = odds / (1.0 + odds)
    if u < moved:
        return u * share / moved
    return share + (u - moved) * (1.0 - share) / (1.0 - moved)


def _inside_share(initial_mass_sol):
    """The share of a star's companions that are inside `CLOSE_COMPANION_AXIS_AU`: the part of the period law of its mass
    under the period of that axis for the star and a companion of the mean mass ratio."""
    low, high = companion_ratio_range(initial_mass_sol)
    total = initial_mass_sol * (1.0 + 0.5 * (low + high))
    period_d = _yr_to_d(math.sqrt(CLOSE_COMPANION_AXIS_AU ** 3 / total))
    mu, sigma = PERIOD_LAWS_LOG10_DAYS["massive" if initial_mass_sol >= MASSIVE_PRIMARY_MIN_MASS_SOL else "solar"]
    return NormalDist(mu, sigma).cdf(math.log10(period_d))


def close_trend(initial_mass_sol, fe_h, prevalence=1.0):
    """
    The `star-count` weights of a star of this birth [Fe/H], and the weight its root pair's inside range carries (GEN.103,
    GEN.705). The table is the Sun's neighbourhood (`fe_h` 0); a metal-poorer star has more companions inside
    `CLOSE_COMPANION_AXIS_AU`, a metal-richer one fewer, by `close_binary_trend.scaled_close_chance` (the factor on the *odds*
    of a companion inside the axis), and the chance of one outside it stays. The chance of any companion moves with it, and
    the share of the companions that are inside moves with it too: the weight returned multiplies the odds of the root pair
    being inside (`_Builder.period_days`).

    Returns:
        tuple: `(weights, inside_weight)`; `(None, 1.0)` for a star with no [Fe/H] or at the factor's 1, where the table stands.
    """
    if close_binary_trend.close_companion_factor(fe_h) == 1.0:
        return None, 1.0
    weights = star_count_weights(initial_mass_sol, prevalence)
    inside = _inside_share(initial_mass_sol)
    total = sum(weights)
    chance = 1.0 - weights[0] / total
    if not (0.0 < chance < 1.0 and 0.0 < inside < 1.0):
        return None, 1.0
    inside_chance = chance * inside
    moved = close_binary_trend.scaled_close_chance(inside_chance, fe_h)
    new_chance = min(moved + chance - inside_chance, 1.0 - 1.0e-9)
    new_inside = moved / new_chance
    scaled = (total * (1.0 - new_chance),) + tuple(w * new_chance / chance for w in weights[1:])
    return scaled, (new_inside / (1.0 - new_inside)) / (inside / (1.0 - inside))


def expected_star_count(initial_mass_sol, max_stars=MAX_MEMBERS, prevalence=1.0):
    """The mean number of stars in the system of a star of this initial mass (GEN.545): the `star-count` table's 1 to 7
    stars, and the tail from 7 up to `max_stars` that `_draw_star_count` walks, in expectation."""
    weights = star_count_weights(initial_mass_sol, prevalence)
    total = sum(weights)
    mean = sum((count + 1) * weight for count, weight in enumerate(weights)) / total
    ratio, reach = TAIL_RATIOS[0], 1.0
    for _ in range(7, max_stars):
        reach *= min(1.0, ratio * prevalence)
        mean += weights[-1] / total * reach
        ratio = TAIL_RATIOS[1]
    return mean


def _draw_star_count(initial_mass_sol, max_stars, prevalence, rng, weights=None):
    """The number of stars: one draw of the `star-count` table, or of `weights` (`close_trend`'s) when the caller moved them."""
    if weights is None:
        count = hosts.draw_star_count(companion_class(initial_mass_sol), prevalence, rng)
    else:
        count = rng.choices(range(1, len(weights) + 1), weights=weights)[0]
    ratio = TAIL_RATIOS[0]
    while count >= 7 and count < max_stars and rng.random() < min(1.0, ratio * prevalence):
        count += 1
        ratio = TAIL_RATIOS[1]
    return min(count, max_stars)


def _radius_au(state):
    if state["radius_km"] is not None:
        return _m_to_au(_km_to_m(state["radius_km"]))
    radius_sol = math.sqrt(state["luminosity_sol"]) * (tuning.SUN_EFFECTIVE_TEMPERATURE_K / state["temperature_k"]) ** 2
    return _m_to_au(radius_sol * constants.SOLAR_RADIUS_M)


def _leaf_of(record):
    """The roster leaf of a member record: what the tree draws the orbits from."""
    return _Leaf(record.kind, record.state, record.mass_solar, record.initial_mass_solar, record.radius_au,
                 radius_max_au=record.radius_max_au, mass_max_solar=record.mass_max_solar, seed=record.seed)


def _star_leaf(initial_mass, age_yr, rng, lighten=True, fe_h=0.0, substellar=False):
    """A star-like member of this birth mass at the system's age: its record (`member_record.make`, under a seed drawn here, which
    the member keeps, so the build makes the same star). A `substellar` one was drawn from the brown dwarf mass function (the
    objects of 0.074 to 0.08 solar masses among them are the stable late-M and L dwarfs of the census, `physics/brown_dwarf.py`)."""
    return _leaf_of(member_record.make(rng.getrandbits(63), initial_mass, age_yr, fe_h, lighten, substellar))


def _schwarzschild_au(mass_solar):
    return (2.0 * constants.G * mass_solar * constants.SOLAR_MASS_TO_KG / constants.SPEED_OF_LIGHT_M_S ** 2
            / constants.AU_M)


def _remnant_leaf(kind, progenitor, mass, direct=False):
    if kind == "white-dwarf":
        radius = _m_to_au(_km_to_m(evolution.white_dwarf_radius_km(mass)))
        return _Leaf(kind, "VII", mass, progenitor, radius)
    if kind == "neutron-star":
        return _Leaf(kind, None, mass, progenitor,
                     _m_to_au(_km_to_m(neutron_star.radius_km(mass))))
    return _Leaf(kind, None, mass, progenitor, _schwarzschild_au(mass), MEMBER_DIRECT_COLLAPSE if direct else 0)


def _remnant(kind, age_yr, lowest, rng, limit_mass=None, fe_h=0.0):
    """A remnant of this kind that has had time to form by `age_yr`, from a progenitor no lighter than
    `lowest`, or None when no progenitor can have collapsed yet. A neutron star or a black hole has a progenitor from the
    windows of its kind (`remnant_progenitors`), and a black hole's is a star that failed to explode."""
    if kind == "white-dwarf":
        low, high = tuning.WHITE_DWARF_PROGENITOR_SOL
        low = max(low, lowest, evolution.youngest_progenitor_sol(age_yr, fe_h))
        if low >= high:
            return None
        progenitor = imf.sample(low, high, rng)
        record = member_record.make(rng.getrandbits(63), progenitor, age_yr, fe_h, lighten=False)
        return _leaf_of(record) if record.kind == "white-dwarf" else None
    progenitor = remnant_progenitors.sample_mass_sol(kind, max(lowest, evolution.youngest_progenitor_sol(age_yr, fe_h)),
                                                     rng)
    if progenitor is None:
        return None
    mass = remnant_mass.draw_mass_sol(kind, rng=rng)
    if limit_mass is not None:
        mass = min(mass, limit_mass)
    return _remnant_leaf(kind, progenitor, mass, direct=kind == "black-hole")


def _alive_star_mass(low, high, age_yr, rng, fe_h=0.0):
    """An IMF mass in `[low, high]` that is still on the main sequence at this age (the heaviest such mass
    caps the range) for a star of birth [Fe/H] `fe_h`."""
    cap = ms_lifetime.mass_for_lifetime_gy(age_yr / 1.0e9, fe_h)
    high = max(min(high, cap), low * 1.0001)
    return imf.sample(low, high, rng)


def _is_direct_collapse(leaf):
    if leaf.kind != "black-hole":
        return False
    return leaf.direct_collapse


class _Failed(Exception):
    """A pair could not be placed in `PAIR_RETRIES` tries (its draws all found no room or failed a screen, `draw.counting` says which);
    the tree restarts."""


UNFUSING_KINDS = ("brown-dwarf", "pre-main-sequence", "protostar")
"""tuple: The kinds that have not begun (or never will begin) hydrogen fusion in the core: a big radius of theirs is no
Algol and no contact binary (`_gate_forbids`)."""


def _gate_forbids(interaction, heavy_initial_sol, age_yr, kinds=()):
    """Whether the catalogue's `where` column rules a pair's `interaction` out (GEN.270; audit test 11): a contact binary
    of solar-type stars needs the angular momentum loss of the intermediate and old populations (the young population has
    only the massive over-contact ones), an Algol is a young or intermediate-age binary
    (`contact-binary` and `algol` in `SYSTEM_TYPE_RATES`), and neither is a pair that holds a star not yet fusing
    (`UNFUSING_KINDS` among the `kinds` of its two members: a contracting or substellar object fills its lobe without being
    an interacting binary). The generator knows the age, not the population, so the gates are the ages of the populations'
    windows (`tuning.YOUNG_POPULATION_MAX_AGE_GY`, `OLD_POPULATION_MIN_AGE_GY`)."""
    if interaction not in ("contact", "semi-detached"):
        return False
    if any(kind in UNFUSING_KINDS for kind in kinds):
        return True
    age_gy = age_yr / 1.0e9
    if interaction == "contact":
        return age_gy < tuning.YOUNG_POPULATION_MAX_AGE_GY and heavy_initial_sol < MASSIVE_PRIMARY_MIN_MASS_SOL
    return age_gy >= tuning.OLD_POPULATION_MIN_AGE_GY


def _history_member(leaf):
    """A pair member as `pair_history` reads it: present and largest mass (solar masses) and radius (AU); a neutron star or a black hole lost its mass in a supernova."""
    return pair_history.Member(leaf.mass, leaf.largest_mass_solar, leaf.radius_au, leaf.largest_radius_au,
                               sudden=leaf.kind in SUDDEN_MASS_LOSS_KINDS)


def _detached_or_post_common_envelope(a, b, axis_au, e):
    """`detached`, or `post-common-envelope` when the history rule (GEN.1012; `pair_history`) finds that a member's largest radius
    reached its Roche lobe, or its partner, at the orbit the pair had when its members were their heaviest."""
    if pair_history.must_be_post_common_envelope(axis_au, e, _history_member(a), _history_member(b)):
        return "post-common-envelope", ()
    return "detached", ()


def classify_pair(a, b, axis_au, e):
    """
    The `interaction` of a pair of two single members from their geometry (GEN.269): `(interaction, gainers)`, where
    `gainers` are the members that receive mass (a tuple of `a` and `b`). `a` and `b` have a `kind`, a `radius_au`, a
    `mass`, an `initial` (a white dwarf's progenitor mass, solar masses) and the largest radius and mass they have had
    (`largest_radius_au`, `largest_mass_solar`). A lobe-filling donor feeds an accretor (`accreting`); two stars past
    their lobes at closest approach are `contact`, one is `semi-detached`; a pair that is not filling a lobe now but whose
    history forced a common envelope (`pair_history.must_be_post_common_envelope`, any member and any depth, GEN.1012)
    is `post-common-envelope`; otherwise `detached`.
    """
    peri = axis_au * (1.0 - e)
    fill_a = roche.filling_factor(a.radius_au, peri, a.mass, b.mass)
    fill_b = roche.filling_factor(b.radius_au, peri, b.mass, a.mass)
    accretor_a, accretor_b = a.kind in ACCRETOR_KINDS, b.kind in ACCRETOR_KINDS
    if accretor_a != accretor_b:
        donor_fill = fill_b if accretor_a else fill_a
        if donor_fill >= 1.0:
            return "accreting", ()
        return _detached_or_post_common_envelope(a, b, axis_au, e)
    if accretor_a:
        return _detached_or_post_common_envelope(a, b, axis_au, e)
    if fill_a >= 1.0 and fill_b >= 1.0:
        return "contact", (a, b)
    if fill_a >= 1.0:
        return "semi-detached", (b,)
    if fill_b >= 1.0:
        return "semi-detached", (a,)
    return _detached_or_post_common_envelope(a, b, axis_au, e)


def period_days(law, mass, a_lo, a_hi, massive, inside_weight=1.0, u=0.5):
    """The period, days, of a pair of total mass `mass` (solar masses) from its law cut to the axes `a_lo` to `a_hi` (AU), drawn from
    the unit `u`; None when the cut leaves no period. `massive` and `inside_weight` are the primary's (the law of a massive primary;
    `close_trend`'s weight on an axis inside `CLOSE_COMPANION_AXIS_AU`)."""
    lo = _yr_to_d(math.sqrt(a_lo ** 3 / mass))
    hi = _yr_to_d(math.sqrt(a_hi ** 3 / mass)) if math.isfinite(a_hi) else math.inf
    if law == "compact-pair":
        lo = max(lo, _yr_to_d(math.sqrt(COMPACT_PAIR_AXIS_AU[0] ** 3 / mass)))
        hi = min(hi, _yr_to_d(math.sqrt(COMPACT_PAIR_AXIS_AU[1] ** 3 / mass)))
    elif law == "interacted":
        lo, hi = max(lo, INTERACTED_PERIOD_DAYS[0]), min(hi, INTERACTED_PERIOD_DAYS[1])
    elif law == "stripped":
        lo, hi = max(lo, tuning.STRIPPED_PAIR_PERIOD_DAYS[0]), min(hi, tuning.STRIPPED_PAIR_PERIOD_DAYS[1])
    if not lo < hi:
        return None
    log_lo, log_hi = math.log10(lo), math.log10(hi) if math.isfinite(hi) else 12.0
    if law != "star":
        return 10.0 ** (log_lo + u * (log_hi - log_lo))
    mu, sigma = PERIOD_LAWS_LOG10_DAYS["massive" if massive else "solar"]
    dist = NormalDist(mu, sigma)
    p_lo, p_hi = dist.cdf(log_lo), dist.cdf(log_hi)
    if p_hi - p_lo < 1.0e-12:
        return 10.0 ** (log_lo + u * (log_hi - log_lo))
    # The tight-tail weight (design note 13.7) multiplies the density under a day of the pairs of solar-type
    # primaries; at 1 this is the plain draw, and the massive law is never weighted.
    weight = 1.0 if massive else system_types.TUNING_KNOBS["tight_tail_weight"]
    split = min(max(dist.cdf(0.0), p_lo), p_hi)
    tail, rest = (split - p_lo) * weight, p_hi - split
    if inside_weight != 1.0:
        # The odds of an axis inside `CLOSE_COMPANION_AXIS_AU` follow the star's [Fe/H] (`close_trend`): the draw's
        # inside part grows or shrinks, the outside part keeps its shape.
        edge = min(max(math.log10(_yr_to_d(math.sqrt(CLOSE_COMPANION_AXIS_AU ** 3 / mass))), log_lo), log_hi)
        reach = min(max(dist.cdf(edge), p_lo), p_hi)
        inside = (reach - p_lo) * weight if reach <= split else tail + (reach - split)
        u = _tilted(u, inside / (tail + rest), inside_weight)
    target = u * (tail + rest)
    point = p_lo + target / weight if target < tail else split + (target - tail)
    return 10.0 ** dist.inv_cdf(min(max(point, 1.0e-15), 1.0 - 1.0e-15))


class _Builder:
    """Draws the orbits of one tree shape, top down."""

    def __init__(self, leaves, age_yr, rng, fixed=None, max_axis_au=math.inf, inside_weight=1.0, fe_h=0.0):
        self.max_axis_au = max_axis_au
        self.fe_h = fe_h
        """float: The birth [Fe/H] of the tree's members (the placed star's), which sets how large a giant of a given core is
        (`helium_core_white_dwarf`)."""
        self.inside_weight = inside_weight
        """float: The weight on the odds that the root pair is inside `CLOSE_COMPANION_AXIS_AU` (`close_trend`)."""
        self.leaves = [dataclasses.replace(leaf) for leaf in leaves]
        """list: This draw's own copies of the roster's leaves: a white dwarf that proves to be a helium core (`helium_core`) is
        replaced in them, and a draw that restarts leaves the roster as it was."""
        self.age_yr = age_yr
        self.rng = rng
        self.nodes = []
        self.by_sub = {}
        self.extra_flags = {}
        self.fixed = fixed or {}
        """dict: `id(pair) -> (period_days, axis_au, eccentricity)` of a pair whose size and shape are kept (the inner
        pair of a contact binary that gains a tertiary): its orientation is drawn and the stability screens run, its
        period is not."""

    def leaves_of(self, sub):
        return (sub,) if isinstance(sub, int) else self.leaves_of(sub[0]) + self.leaves_of(sub[1])

    def mass(self, sub):
        return sum(self.leaves[i].mass for i in self.leaves_of(sub))

    def extent(self, sub):
        return sum(self.leaves[i].radius_au for i in self.leaves_of(sub))

    def law(self, left, right):
        if isinstance(left, int) and isinstance(right, int):
            a, b = self.leaves[left], self.leaves[right]
            if a.compact and b.compact:
                return "compact-pair"
            if a.compact or b.compact:
                return "interacted"
            if "stripped-helium" in (a.kind, b.kind):
                return "stripped"
        return "star"

    def pair_leaves(self, left, right):
        """The two leaves of a pair of single members, else None."""
        if isinstance(left, int) and isinstance(right, int):
            return self.leaves[left], self.leaves[right]
        return None

    def side(self, sub, placed=None):
        """The `Side` of `sub` for the keep-out rule: a single member's radius, else the node's extent. Before the orbits below it are
        drawn (`placed` None) a node's extent is the sum of the radii of its stars, the least it can be; with `placed`
        (`id(sub) -> record` of the drawn nodes) it is the barycentric extent of the orbits as drawn (`node_extent_au`)."""
        if isinstance(sub, int):
            leaf = self.leaves[sub]
            return single_side(leaf.mass, leaf.radius_au)
        if placed is None:
            return self.min_side(sub)
        record = placed[id(sub)]
        low, high = self.side(record["left"], placed), self.side(record["right"], placed)
        return Side(low.mass + high.mass, node_extent_au(record["separation_au"], record["eccentricity"], low, high))

    def interaction(self, left, right, axis, e):
        """The pair's `interaction` and the flags it puts on the two single members: `(interaction, {leaf: flag})`."""
        pair = self.pair_leaves(left, right)
        if pair is None:
            return "detached", {}
        interaction, gainers = classify_pair(pair[0], pair[1], axis, e)
        flags = {(left if gainer is pair[0] else right): MEMBER_MASS_GAINER for gainer in gainers}
        if interaction == "semi-detached":
            # An Algol: the donor was the heavier at birth (it evolved first), which is how it came to be the lighter.
            gainer = gainers[0]
            donor, donor_index = (pair[1], right) if gainer is pair[0] else (pair[0], left)
            if donor.initial > gainer.initial:
                flags[donor_index] = MEMBER_MASS_REVERSED
        return interaction, flags

    def orbit_span(self, left, right):
        """The `orbit_au` `(low, high)` a pair of two single members carries on one of them (GEN.271), else None."""
        pair = self.pair_leaves(left, right)
        if pair is None:
            return None
        return next((leaf.orbit_au for leaf in pair if leaf.orbit_au is not None), None)

    def double_white_dwarf(self, left, right):
        """Whether the pair is two white dwarfs, which the Peters rule screens as it does a compact pair."""
        pair = self.pair_leaves(left, right)
        return pair is not None and pair[0].kind == "white-dwarf" and pair[1].kind == "white-dwarf"

    def wide_cap(self, sub):
        """The widest axis the pair can have: 1,000 AU if it holds a neutron star or a black hole that did not
        collapse directly, else unbounded."""
        if any(self.leaves[i].compact and not _is_direct_collapse(self.leaves[i]) for i in self.leaves_of(sub)):
            return WIDE_COMPACT_AXIS_AU
        return math.inf

    def period_days(self, law, mass, a_lo, a_hi, massive, inside_weight=1.0):
        return period_days(law, mass, a_lo, a_hi, massive, inside_weight, self.rng.random())

    def eccentricity(self, law, period_d, a_au, root, limit=math.inf):
        """The pair's eccentricity from the law of its kind cut at `limit`, the largest the room it has allows (GEN.1011): the
        thermal law `2 e` on `[0, cap]` cut at `limit` has the inverse `limit * sqrt(u)`, the exact conditional draw."""
        u = self.rng.random()
        if (law == "compact-pair" and a_au < COMPACT_PAIR_CIRCULAR_BELOW_AU) or period_d < CIRCULAR_BELOW_DAYS:
            return 0.0
        cap = min(1.0 - (period_d / 2.0) ** (-2.0 / 3.0),
                  ROOT_ECCENTRICITY_CAP if root else INNER_ECCENTRICITY_CAP, limit)
        return cap * math.sqrt(u)

    def min_side(self, sub):
        """The `Side` of `sub` at its least: its stars at the closest their own floors let them be, on circular orbits."""
        if isinstance(sub, int):
            leaf = self.leaves[sub]
            return single_side(leaf.mass, leaf.radius_au)
        low, high = self.min_side(sub[0]), self.min_side(sub[1])
        return Side(low.mass + high.mass, node_extent_au(pair_floor_au(low, high), 0.0, low, high))

    def _apocentre_limit_au(self, budget, low, high):
        """The largest apocentre the pair may have so its extent stays inside `budget` (AU), `inf` without one."""
        if not math.isfinite(budget):
            return math.inf
        total = low.mass + high.mass
        return min((budget - low.extent) * total / high.mass, (budget - high.extent) * total / low.mass)

    def _eccentricity_limit(self, axis, peri_floor, apocentre_limit, stable_at):
        """The largest eccentricity at `axis` that keeps the pair's pericentre clear, its extent inside its budget and the triple above it
        stable (`stable_at(e)`, monotone in e: bisected), negative when even a circular orbit does not fit."""
        limit = 1.0 - peri_floor / axis
        if math.isfinite(apocentre_limit):
            limit = min(limit, apocentre_limit / axis - 1.0)
        if limit < 0.0:
            return limit
        if stable_at is None or stable_at(limit):
            return limit
        if not stable_at(0.0):
            return -1.0
        low, high = 0.0, limit
        for _ in range(30):
            middle = 0.5 * (low + high)
            low, high = (middle, high) if stable_at(middle) else (low, middle)
        return low

    def place(self, sub, parent, slot, parent_no, budget=math.inf):
        """
        Places the orbit of the pair `sub`, then its child pairs. `parent` is `(axis, eccentricity, inclination, node,
        sibling_mass)` or None for the root. `budget` is the most extent (AU) the pair may have, from the keep-out rule of the
        pair above it (`inf` for the root).

        The period and eccentricity are drawn once, from their laws cut to the range the pair's room allows (GEN.1011): the axis
        between the keep-out floor and the least of the stability limit, the cap and the extent budget; the eccentricity
        up to what the pericentre floor, the budget and the stability of the triple above allow. A draw can then only fail
        the screens that depend on the draw itself, or find no room (`pair.no_allowed_period`; `draw.redrawn` counts which).
        """
        left, right = sub
        if self.mass(left) < self.mass(right):
            left, right = right, left
        m_left, m_right = self.mass(left), self.mass(right)
        total = m_left + m_right
        drawn = left, right, m_left, m_right, total
        core = None
        law = self.law(left, right)
        extent = self.extent(left) + self.extent(right)
        sides = self.side(left), self.side(right)
        a_lo = pair_floor_au(*sides)
        peri_floor = pair_pericentre_floor_au(*sides)
        apocentre_limit = self._apocentre_limit_au(budget, *sides)
        cap = min(self.wide_cap(sub), self.max_axis_au, apocentre_limit)
        heavy = max(self.leaves[i].initial for i in self.leaves_of(sub))
        massive = heavy >= MASSIVE_PRIMARY_MIN_MASS_SOL
        fixed = self.fixed.get(id(sub))
        span = self.orbit_span(left, right)
        double_dwarf = self.double_white_dwarf(left, right)
        for _ in range(PAIR_RETRIES):
            left, right, m_left, m_right, total = drawn
            core = None
            rng = self.rng
            inclination = math.acos(1.0 - 2.0 * rng.random())
            node = 2.0 * math.pi * rng.random()
            arg = 2.0 * math.pi * rng.random()
            anomaly = 2.0 * math.pi * rng.random()
            mutual = None
            a_hi = ROOT_MAX_AXIS_AU
            stable_at = None
            if parent is not None:
                mutual = stability.mutual_inclination(inclination, node, parent[2], parent[3])
                a_hi = stability.max_inner_axis(m_left, m_right, parent[4], parent[0], 0.0, parent[1], mutual,
                                                safety=STABILITY_MARGIN)

                def stable_at(e, mutual=mutual):
                    return stability.is_triple_stable(m_left, m_right, parent[4], axis, parent[0], e, parent[1], mutual,
                                                      safety=STABILITY_MARGIN)
            a_hi = min(a_hi, cap)
            if fixed is None and span is not None:
                low, high = max(a_lo, span[0]), min(a_hi, span[1])
                if low > high:                # a span of one axis (a cataclysmic variable's, an AM CVn's) has low == high
                    draw.redrawn("pair.axis_out_of_range")
                    continue
                axis = 10.0 ** (math.log10(low) + rng.random() * (math.log10(high) - math.log10(low)))
                period_d = 365.25 * math.sqrt(axis ** 3 / total)
            elif fixed is None:
                period_d = self.period_days(law, total, a_lo, a_hi, massive,
                                            self.inside_weight if parent is None else 1.0)
                if period_d is None:
                    draw.redrawn("pair.no_allowed_period")
                    continue
                axis = stability.axis_au(_d_to_yr(period_d), total)
            else:
                period_d, axis, e = fixed
                if axis > a_hi:
                    draw.redrawn("pair.axis_out_of_range")
                    continue
            if fixed is None:
                limit = self._eccentricity_limit(axis, peri_floor, apocentre_limit, stable_at)
                if limit < 0.0:
                    draw.redrawn("pair.no_allowed_period")
                    continue
                e = self.eccentricity(law, period_d, axis, parent is None, limit)
            if not pair_clears(axis, e, *sides):
                draw.redrawn("pair.overlap")         # cannot happen from a cut draw; counted so TEST.476 can hold it at zero
                continue
            circular = (lambda both: both is not None and any(leaf.circular for leaf in both))(self.pair_leaves(left, right))
            if circular:
                e = 0.0                      # the companion of a helium-core white dwarf the star draw made (GEN.996)
            interaction, member_flags = self.interaction(left, right, axis, e)
            if circular:
                interaction = "post-common-envelope"
            core = self.helium_core(left, right, axis, e, interaction)
            if core is not None:
                # The giant that left this white dwarf filled its lobe on the red giant branch (GEN.996): the white dwarf is the
                # lighter helium core, so the pair is lighter, may swap its heavier member, and has a longer period.
                index, replaced = core
                m_left, m_right = (m_left - self.leaves[index].mass + replaced.mass if left == index else m_left,
                                   m_right - self.leaves[index].mass + replaced.mass if right == index else m_right)
                total = m_left + m_right
                period_d = 365.25 * math.sqrt(axis ** 3 / total)
                e = 0.0                      # the overflow circularised the orbit: the lobe was filled at this axis
                if m_left < m_right:
                    left, right, m_left, m_right = right, left, m_right, m_left
            pair = self.pair_leaves(left, right)
            accretes = None if pair is None else next((leaf.accretes for leaf in pair if leaf.accretes), None)
            if accretes is not None:
                interaction = "accreting"  # the wind of a supergiant or the disc of a Be star feeds it without a full lobe
            elif interaction == "accreting" and any(leaf.compact for leaf in pair):
                draw.redrawn("pair.interaction_gate")
                continue  # a neutron star or a black hole accretes only in the X-ray binary draws, which fix the rate
            if interaction == "accreting" and any(leaf.kind == "white-dwarf" for leaf in pair) and (
                    span is None or span[0] < span[1]):
                draw.redrawn("pair.interaction_gate")
                continue  # a white dwarf accretes only in the cataclysmic variable and AM CVn draws, which fix the rate
            if _gate_forbids(interaction, heavy, self.age_yr, () if pair is None else (pair[0].kind, pair[1].kind)):
                draw.redrawn("pair.interaction_gate")
                continue
            flags = NODE_INTERACTED if law != "star" or interaction in ("accreting", "post-common-envelope") else 0
            flags |= {"wind": NODE_WIND_FED, "disc": NODE_DISC_FED}.get(accretes, 0)
            if parent is not None:
                if not stable_at(e):
                    draw.redrawn("pair.unstable_triple")   # likewise cannot happen from a cut draw
                    continue
                verdict = stability.kozai_verdict(m_left, m_right, parent[4], axis, parent[0], e, parent[1], mutual,
                                                  self.age_yr, extent)
                if verdict == "unsafe":
                    draw.redrawn("pair.kozai")
                    continue
                if verdict == "kl_active":
                    flags |= NODE_KL_ACTIVE
            if ((law == "compact-pair" or double_dwarf)
                    and not stability.survives_gravitational_waves(m_left, m_right, axis, e, self.age_yr)):
                draw.redrawn("pair.gravitational_waves")
                continue
            if self._in_desert(left, right, axis):
                draw.redrawn("pair.brown_dwarf_desert")
                continue
            break
        else:
            raise _Failed
        for leaf_index, flag in member_flags.items():
            self.extra_flags[leaf_index] = self.extra_flags.get(leaf_index, 0) | flag
        if core is not None:
            self._replace_helium_core(*core, parent_no)
        node_no = len(self.nodes) + 1
        kind = "compact" if law in ("compact-pair", "interacted") else "close" if axis < 50.0 else "wide"
        record = {"node_no": node_no, "parent_node_no": parent_no, "slot": slot, "mass_solar": total,
                  "separation_au": axis, "eccentricity": e, "inclination_deg": math.degrees(inclination),
                  "ascending_node_deg": math.degrees(node), "arg_periapsis_deg": math.degrees(arg),
                  "mean_anomaly_deg": math.degrees(anomaly), "period_years": _d_to_yr(period_d), "kind": kind,
                  "flags": flags, "left": left, "right": right, "sub": sub,
                  "interaction": interaction}
        self.nodes.append(record)
        self.by_sub[id(sub)] = record
        peri, apocentre = axis * (1.0 - e), axis * (1.0 + e)
        for child_slot, (child, sibling) in enumerate(((left, right), (right, left))):
            if isinstance(child, int):
                continue
            # The room the child has is what this pair's keep-out rule and its own budget leave once its sibling is
            # accounted for (the sibling's extent as drawn, or its least when it is not drawn yet).
            sibling_extent = self.side(sibling, self.by_sub).extent if isinstance(sibling, int) or id(sibling) in self.by_sub \
                else self.min_side(sibling).extent
            room = peri / PERICENTRE_EXTENT_FACTOR - sibling_extent
            if math.isfinite(budget):
                room = min(room, budget - apocentre * self.mass(sibling) / total)
            self.place(child, (axis, e, inclination, node, self.mass(sibling)), child_slot, node_no, room)
        return node_no

    def helium_core(self, left, right, axis, e, interaction):
        """
        `(index, replacement)` when the pair is a post-common-envelope pair whose white dwarf is a helium core (GEN.996), else None:
        the white dwarf's giant filled its Roche lobe at the pair's axis before it reached the helium flash, so the white dwarf is
        the core it had then (`helium_core_white_dwarf.overflow_core_mass_sol`), and `replacement` is the leaf made from the member's
        own record with that mass. Only a white dwarf the roster drew from a record (it has a seed), not the placed object, which is
        built before the tree. The companion's mass is the one it had when the giant overflowed: its birth mass when it evolves
        later, its present one when it is a remnant already. Of a pair of two such white dwarfs the one that evolved first.
        """
        pair = self.pair_leaves(left, right)
        if pair is None or interaction != "post-common-envelope":
            return None
        candidates = [(index, other) for index, other in ((left, right), (right, left))
                      if index != 0 and self.leaves[index].kind == "white-dwarf" and self.leaves[index].seed is not None
                      and self.leaves[index].orbit_au is None
                      and self.leaves[index].mass >= tuning.HELIUM_FLASH_CORE_MASS_SOL]
        if not candidates:
            return None
        index, other = max(candidates, key=lambda c: self.leaves[c[0]].initial)
        leaf, partner = self.leaves[index], self.leaves[other]
        partner_mass = partner.initial if partner.initial <= leaf.initial else partner.mass
        mass = helium_core_white_dwarf.overflow_core_mass_sol(leaf.initial, partner_mass, axis, 0.0, self.fe_h)
        if mass is None:
            return None
        record = member_record.make(leaf.seed, leaf.initial, self.age_yr, self.fe_h, lighten=False, white_dwarf_mass_sol=mass)
        if record.kind != "white-dwarf":
            return None
        replacement = _leaf_of(record)
        replacement.flags, replacement.effective_age_yr = leaf.flags, leaf.effective_age_yr
        return index, replacement

    def _replace_helium_core(self, index, replacement, parent_no):
        """Puts the helium core white dwarf `replacement` in place of leaf `index` and takes its mass off the nodes above it, whose
        periods follow (an orbit keeps its axis when the mass in it changes)."""
        lost = self.leaves[index].mass - replacement.mass
        self.leaves[index] = replacement
        while parent_no is not None:
            node = self.nodes[parent_no - 1]
            node["mass_solar"] -= lost
            node["period_years"] = math.sqrt(node["separation_au"] ** 3 / node["mass_solar"])
            parent_no = node["parent_node_no"]

    def _in_desert(self, left, right, axis):
        if axis >= BROWN_DWARF_DESERT_AU or not (isinstance(left, int) and isinstance(right, int)):
            return False
        kinds = {self.leaves[left].kind, self.leaves[right].kind}
        if "brown-dwarf" not in kinds or kinds == {"brown-dwarf"}:
            return False
        return self.rng.random() >= BROWN_DWARF_DESERT_KEEP

    def clears_extents(self):
        """The final check: every pair's closest approach clears what it holds (`pair_clears` with the extents of the orbits as
        drawn)."""
        by_sub = {id(record["sub"]): record for record in self.nodes}
        return all(pair_clears(record["separation_au"], record["eccentricity"], self.side(record["left"], by_sub),
                               self.side(record["right"], by_sub)) for record in self.nodes)

    def finish(self):
        """Flags the members of an ultra-wide system `fragile` and the theoretical kinds `theoretical`."""
        fragile = bool(self.nodes) and self.nodes[0]["separation_au"] > FRAGILE_AXIS_AU
        for index, leaf in enumerate(self.leaves):
            if fragile:
                self.extra_flags[index] = self.extra_flags.get(index, 0) | MEMBER_FRAGILE
            if leaf.kind == "thorne-zytkow":
                self.extra_flags[index] = self.extra_flags.get(index, 0) | MEMBER_THEORETICAL


def built_tree_clears(nodes, stars):
    """
    Whether the built stars of a tree clear its orbits: the same `pair_clears` rule `_Builder.clears_extents` runs on the drawn leaves,
    run on the radii and masses the build gave the stars (`nodes` the tree's `OrbitNode`s, `stars` `(parent_node_no, radius_au,
    mass_solar)` for each).
    """
    below = {}
    for parent, radius_au, mass_solar in stars:
        below.setdefault(parent, []).append(single_side(mass_solar, radius_au))
    by_parent = {}
    for node in nodes:
        by_parent.setdefault(node.parent_node_no, []).append(node)
    sides = {}

    def side_of(node):
        if node.node_no not in sides:
            low, high = children_of(node)
            sides[node.node_no] = Side(low.mass + high.mass, node_extent_au(
                _km_to_au(node.separation_km), node.eccentricity, low, high))
        return sides[node.node_no]

    def children_of(node):
        found = below.get(node.node_no, []) + [side_of(child) for child in by_parent.get(node.node_no, ())]
        return found[0], found[1]

    return all(pair_clears(_km_to_au(node.separation_km), node.eccentricity, *children_of(node)) for node in nodes)


def _shape(units, rng):
    """A random binary tree over `units` (indices or nested pairs of them)."""
    if len(units) == 1:
        return units[0]
    items = list(units)
    rng.shuffle(items)
    count = len(items)
    if count <= 3:
        split = 1
    elif count == 4:
        split = 1 if rng.random() < 0.6 else 2
    else:
        sizes = list(range(1, count // 2 + 1))
        split = rng.choices(sizes, [1.5 if size == 2 else 1.0 for size in sizes])[0]
    return (_shape(items[:split], rng), _shape(items[split:], rng))


STAR_TYPE_DRAWS = 10
"""int: The uniforms a placed star's tree draws first, in a fixed order, for the close-binary types (GEN.270), the giants'
white dwarf companions (GEN.271) and the high-mass X-ray binaries (GEN.272): the merger product, the Be star, its stripped
companion, the mass gainer, the exotic star, the merger's rejuvenation, the barium or symbiotic star, whether the compact
companion of a massive star accretes, whether it is fed by a Be star's disc, and whether it launches jets."""


def _young(age_yr):
    return age_yr / 1.0e9 < tuning.YOUNG_POPULATION_MAX_AGE_GY


def draw_merger_product(birth_mass_sol, age_yr, unit):
    """
    Whether a placed star is a merger product, and which (GEN.270, `massive-merger-product` and `blue-straggler-field` in
    `SYSTEM_TYPE_RATES`; the shares are the `merger_product_fraction` knob): `massive` for a young star of
    `tuning.MERGER_MASSIVE_MIN_MASS_SOL` or more (8%; its tree has one star fewer, the one it ate), `field` for any other
    star of the intermediate and old populations (a blue straggler, 0.5%), else None. `unit` is one uniform.
    """
    shares = system_types.TUNING_KNOBS["merger_product_fraction"]
    if birth_mass_sol >= tuning.MERGER_MASSIVE_MIN_MASS_SOL:
        return "massive" if _young(age_yr) and unit < shares["massive"] else None
    return "field" if not _young(age_yr) and unit < shares["field"] else None


def draw_mass_gainer(birth_mass_sol, age_yr, unit):
    """Whether a young massive star has gained mass from a companion (20%, `massive-merger-product`'s `mass_gainer`)."""
    share = system_types.SYSTEM_TYPE_RATES["massive-merger-product"].by_class["mass_gainer"]
    return birth_mass_sol >= tuning.MERGER_MASSIVE_MIN_MASS_SOL and _young(age_yr) and unit < share


def draw_be_star(birth_mass_sol, age_yr, unit, stripped_unit):
    """
    `(be, stripped_companion)` for a placed star (GEN.270, `be-star`): a young B star is a Be star at the row's 20%, and
    30% of Be stars have a stripped helium companion that spun them up. Two uniforms.
    """
    entry = system_types.SYSTEM_TYPE_RATES["be-star"]
    if spectral_class(birth_mass_sol) != "B" or not _young(age_yr) or unit >= entry.default:
        return False, False
    return True, stripped_unit < entry.by_class["stripped_companion"]


def draw_exotic_star(leaf_kind, age_yr, unit):
    """
    `r-crb`, `fk-comae`, `thorne-zytkow` or None for a placed star (GEN.270, `rcb-fk-comae` and `thorne-zytkow`): an old
    giant is one of them at the type's Galaxy count over the old giants of the Galaxy (`tuning.OLD_GIANT_SYSTEM_SHARE` of
    its systems, a stand-in), one uniform. A Thorne-Zytkow object is theoretical.
    """
    if leaf_kind != "giant" or age_yr / 1.0e9 < tuning.OLD_POPULATION_MIN_AGE_GY:
        return None
    rates = system_types.SYSTEM_TYPE_RATES
    share = tuning.OLD_GIANT_SYSTEM_SHARE
    chances = (("r-crb", system_types.fraction_of_systems("rcb-fk-comae") / share),
               ("fk-comae", rates["rcb-fk-comae"].by_class["FK Comae"] / system_types.GALAXY_SYSTEMS / share),
               ("thorne-zytkow", system_types.fraction_of_systems("thorne-zytkow") / share))
    for kind, chance in chances:
        if unit < chance:
            return kind
        unit -= chance
    return None


def _stripped_helium_leaf(be_birth_mass_sol, rng):
    """A Be star's stripped helium companion: 0.3 to 0.5 solar masses, born heavier than the Be star it spun up."""
    mass = rng.uniform(*tuning.BE_STRIPPED_COMPANION_MASS_RANGE_SOL)
    progenitor = be_birth_mass_sol * rng.uniform(*tuning.BE_STRIPPED_PROGENITOR_RATIO_RANGE)
    radius_sol = tuning.STRIPPED_HELIUM_RADIUS_SOL_AT_0_47 * (mass / 0.47) ** (-1.0 / 3.0)
    return _Leaf("stripped-helium", None, mass, progenitor, _m_to_au(radius_sol * constants.SOLAR_RADIUS_M))


def _swap_units(shape, a, b):
    """`shape` with the members `a` and `b` changing places."""
    if isinstance(shape, int):
        return b if shape == a else a if shape == b else shape
    return (_swap_units(shape[0], a, b), _swap_units(shape[1], a, b))


def _leaf_pairs(shape):
    """The nodes of `shape` whose two children are both single members, as `(left, right)`."""
    if isinstance(shape, int):
        return []
    if isinstance(shape[0], int) and isinstance(shape[1], int):
        return [shape]
    return _leaf_pairs(shape[0]) + _leaf_pairs(shape[1])


def _beside_placed(shape):
    """What hangs beside member 0 (the placed object) in `shape`: a member or a subtree, `None` if 0 is the whole tree."""
    if isinstance(shape, int):
        return None
    left, right = shape
    if isinstance(left, int) and left == 0:
        return right
    if isinstance(right, int) and right == 0:
        return left
    found = _beside_placed(left)
    return found if found is not None else _beside_placed(right)


def _pair_placed(shape):
    """`shape` with the placed object (member 0) beside a single member, when the tree has a pair to give it: the placed
    object trades places with a member of a pair of the subtree beside it (the orbits are drawn after, for the new shape)."""
    beside = _beside_placed(shape)
    if beside is None or isinstance(beside, int):
        return shape
    return _swap_units(shape, 0, _leaf_pairs(beside)[0][0])


def draw_white_dwarf_binary(class_name, state, age_yr, unit):
    """
    The binary type of a placed white dwarf's companion (GEN.271, `close-double-white-dwarf`, `white-dwarf-post-ce`,
    `cataclysmic-variable` and `am-cvn` in `SYSTEM_TYPE_RATES`), given the `state` (`star` or `white-dwarf`) the screen drew
    against `MEMBER_STATE_BY_PLACED[class_name]`: `double-white-dwarf` for a white dwarf companion; for a star, `am-cvn` (a
    system of the old population), `cataclysmic-variable` or `post-common-envelope` (neither of the young population) at
    the rows' rates per white dwarf over the class's odds of a star companion, else None (a pair on the period law). One
    uniform.
    """
    if state == "white-dwarf":
        return "double-white-dwarf"
    if state != "star":
        return None
    rates = system_types.SYSTEM_TYPE_RATES
    table = MEMBER_STATE_BY_PLACED[class_name]
    star_odds = table["star"][0]
    companions = sum(default for default, _low, _high in table.values())
    old = age_yr / 1.0e9 >= tuning.OLD_POPULATION_MIN_AGE_GY
    chances = (("am-cvn", rates["am-cvn"].default / star_odds if old else 0.0),
               ("cataclysmic-variable", rates["cataclysmic-variable"].default / star_odds if not _young(age_yr) else 0.0),
               ("post-common-envelope",
                system_types.TUNING_KNOBS["post_ce_fraction"] * companions / star_odds if not _young(age_yr) else 0.0))
    for kind, chance in chances:
        if unit < chance:
            return kind
        unit -= chance
    return None


def draw_stripped_helium(unit):
    """Whether the star of a post-common-envelope white dwarf pair is a stripped helium star, a hot subdwarf (GEN.271,
    `hot-subdwarf-stripped`; `tuning.HOT_SUBDWARF_POST_CE_SHARE`). One uniform."""
    return unit < tuning.HOT_SUBDWARF_POST_CE_SHARE


def draw_white_dwarf_merger(age_yr, unit):
    """Whether a placed white dwarf of the old population is the product of a double white dwarf merger (GEN.271,
    `double-white-dwarf-merger`; the `white-dwarf` share of the `merger_product_fraction` knob, 10%). One uniform."""
    shares = system_types.TUNING_KNOBS["merger_product_fraction"]
    return age_yr / 1.0e9 >= tuning.OLD_POPULATION_MIN_AGE_GY and unit < shares["white-dwarf"]


def draw_giant_companion(leaf_kind, birth_mass_sol, unit):
    """
    `symbiotic`, `barium` or None for a placed star (GEN.271, `symbiotic-barium`): a G or K giant (a giant born lighter
    than `tuning.BARIUM_GIANT_MAX_INITIAL_MASS_SOL`) is a symbiotic star at the row's 5e-6 and a barium or CH star at its
    1%; both have a white dwarf companion. One uniform.
    """
    if leaf_kind != "giant" or birth_mass_sol >= tuning.BARIUM_GIANT_MAX_INITIAL_MASS_SOL:
        return None
    entry = system_types.SYSTEM_TYPE_RATES["symbiotic-barium"]
    if unit < entry.default:
        return "symbiotic"
    return "barium" if unit - entry.default < entry.by_class["barium and CH (G and K giants)"] else None


def draw_accreting_compact(case, unit, mode_unit=0.0, age_yr=0.0):
    """
    How a neutron star or a black hole accretes (GEN.272), for a registry `case`; `None` for none. One or two uniforms.

    * `high-mass-xray-binary` (a living massive star's compact companion): it accretes at the `accreting_fraction` knob (40%,
      `unit`); 60% of those (`be_xray_share`, `mode_unit`) are fed by a Be star's disc (`disc`, the Be/X-ray binary, of the
      young population only: a halo or old star is never one), the rest by a supergiant's wind (`wind`).
    * `low-mass-xray-binary` (the low-mass star of a neutron star or a black hole): `overflow` for the knob's 1% (`unit`), of
      the intermediate and old populations only.
    * `microquasar` (an accreting neutron star or black hole): `jets` for `tuning.MICROQUASAR_SHARE_OF_XRAY_BINARIES`.
    """
    knobs = system_types.TUNING_KNOBS["accreting_fraction"]
    if case == "high-mass-xray-binary":
        if unit >= knobs[case]:
            return None
        be_share = system_types.SYSTEM_TYPE_RATES[case].by_class["be_xray_share"]
        return "disc" if mode_unit < be_share and _young(age_yr) else "wind"
    if case == "low-mass-xray-binary":
        return "overflow" if unit < knobs[case] and not _young(age_yr) else None
    if case == "microquasar":
        return "jets" if unit < tuning.MICROQUASAR_SHARE_OF_XRAY_BINARIES else None
    raise ValueError(f"unknown accreting-compact case {case!r}")


def draw_low_mass_xray_pair(placed, age_yr, unit):
    """
    Whether the `star` companion of a neutron star or a black hole is the low-mass donor of a low-mass X-ray binary pair
    (GEN.272, `low-mass-xray-binary`): the row's 5e-5 of neutron stars and 1e-3 of black holes are part of the `star` odds, so
    the share of the star companions is the row's rate over those odds (0.25% of a neutron star's, 10% of a black hole's).
    Not for a recycled pulsar (its companions are the spiders, white dwarfs and the rest) and not in the young population.
    One uniform.
    """
    if placed.kind not in COMPACT_KINDS or placed.recycled or _young(age_yr):
        return False
    row = system_types.SYSTEM_TYPE_RATES["low-mass-xray-binary"]
    rate = row.by_class["black hole"] if placed.kind == "black-hole" else row.default
    return unit < rate / MEMBER_STATE_BY_PLACED[placed_class(placed)]["star"][0]


def _xray_span(donor, compact, period_days):
    """The axis span `(low, high)` of an X-ray binary of this period range, days, log-uniform: the range's axes about the
    pair's mass, the low end no closer than keeps the donor inside `MAX_FILL_FACTOR` of its lobe and the high end at least
    `tuning.XRAY_SPAN_ROOM` times the low one."""
    total = donor.mass + compact.mass
    low, high = (stability.axis_au(days / 365.25, total) for days in period_days)
    low = max(low, roche.separation_for_filling(donor.radius_au, MAX_FILL_FACTOR, donor.mass, compact.mass))
    return low, max(high, low * tuning.XRAY_SPAN_ROOM)


def _feed_compact(donor, compact, mode, jet_unit):
    """Makes the pair of `donor` and `compact` an X-ray binary fed in `mode` (GEN.272): the pair carries the mode, its axis
    is log-uniform over the period range of the mode, a Be donor of a `disc` pair is flagged `be`, and the compact object
    launches jets (`microquasar`) at `tuning.MICROQUASAR_SHARE_OF_XRAY_BINARIES`."""
    days = tuning.BE_XRAY_PERIOD_DAYS if mode == "disc" else tuning.SUPERGIANT_XRAY_PERIOD_DAYS
    compact.orbit_au = _xray_span(donor, compact, days)
    compact.accretes = mode
    if mode == "disc":
        donor.flags |= MEMBER_BE
    if draw_accreting_compact("microquasar", jet_unit) is not None:
        compact.flags |= MEMBER_MICROQUASAR


def _low_mass_xray_companion(placed, base, age_yr, rng):
    """The low-mass star of a low-mass X-ray binary pair (`draw_low_mass_xray_pair`, GEN.272), or None when none fits. Its
    pair accretes at the knob's 1% (`draw_accreting_compact`): the donor, a main-sequence star, fills its lobe at a period
    of `tuning.LOW_MASS_XRAY_PERIOD_DAYS`, and the compact object is flagged `microquasar` at its share. The rest are
    detached pairs on the interacted period law."""
    mode = draw_accreting_compact("low-mass-xray-binary", rng.random(), age_yr=age_yr)
    for _ in range(PAIR_RETRIES):
        donor = _star_leaf(_alive_star_mass(*tuning.LOW_MASS_XRAY_DONOR_MASS_RANGE_SOL, age_yr, rng, placed.tree_fe_h), age_yr,
                           rng, fe_h=placed.tree_fe_h)
        if mode is None:
            return donor
        axis = _filling_axis_au(donor, placed.mass_solar, tuning.LOW_MASS_XRAY_PERIOD_DAYS) if (
            donor.kind == "main-sequence") else None
        if axis is not None:
            donor.orbit_au = (axis, axis)
            donor.accretes = mode
            if draw_accreting_compact("microquasar", rng.random()) is not None:
                base.flags |= MEMBER_MICROQUASAR
            return donor
        draw.redrawn("accreting.donor")
    return None


def _spider_companion(placed, base, state, age_yr, rng):
    """The companion of a recycled pulsar that is a black widow (`black-widow`: 0.01 to 0.05 solar masses) or a redback
    (`redback`: 0.2 to 0.4), on an orbit of `tuning.SPIDER_PERIOD_DAYS`, with the pulsar flagged `ablating`."""
    low, high = tuning.BLACK_WIDOW_MASS_RANGE_SOL if state == "black-widow" else tuning.REDBACK_MASS_RANGE_SOL
    companion = _star_leaf(rng.uniform(low, high), age_yr, rng, fe_h=placed.tree_fe_h)
    total = placed.mass_solar + companion.mass
    companion.orbit_au = tuple(stability.axis_au(days / 365.25, total) for days in tuning.SPIDER_PERIOD_DAYS)
    base.flags |= MEMBER_ABLATING
    return companion


def _filling_axis_au(donor, accretor_mass, period_range_days):
    """The axis at which `donor` fills its lobe by `tuning.ACCRETING_FILL_FACTOR` about an accretor of `accretor_mass`, or
    None when the period falls outside `period_range_days`."""
    axis = roche.separation_for_filling(donor.radius_au, tuning.ACCRETING_FILL_FACTOR, donor.mass, accretor_mass)
    period_d = 365.25 * math.sqrt(axis ** 3 / (donor.mass + accretor_mass))
    return axis if period_range_days[0] <= period_d <= period_range_days[1] else None


def _accreting_companion(placed_mass, kind, age_yr, rng, fe_h=0.0):
    """The donor of a cataclysmic variable (`kind` `cataclysmic-variable`: a main-sequence star) or an AM CVn binary
    (`am-cvn`: a degenerate helium star) and its lobe-filling axis, `(leaf, axis_au)`, or None after the draws that fit
    no period."""
    if kind == "cataclysmic-variable":
        window = tuple(hours / 24.0 for hours in tuning.CATACLYSMIC_PERIOD_HOURS)
    else:
        window = tuple(minutes / 1440.0 for minutes in tuning.AM_CVN_PERIOD_MINUTES)
    for _ in range(PAIR_RETRIES):
        if kind == "cataclysmic-variable":
            donor = _star_leaf(rng.uniform(*tuning.CATACLYSMIC_DONOR_MASS_RANGE_SOL), age_yr, rng,
                               fe_h=fe_h)
            if donor.kind != "main-sequence":
                draw.redrawn("accreting.donor")
                continue
        else:
            low, high = tuning.AM_CVN_DONOR_MASS_RANGE_SOL
            mass = math.exp(math.log(low) + rng.random() * (math.log(high) - math.log(low)))
            radius_sol = tuning.AM_CVN_DONOR_RADIUS_SOL_AT_0_05 * (mass / 0.05) ** (-1.0 / 3.0)
            donor = _Leaf("stripped-helium", None, mass, rng.uniform(1.0, 2.5),
                          radius_sol * constants.SOLAR_RADIUS_M / constants.AU_M)
        axis = _filling_axis_au(donor, placed_mass, window)
        if axis is not None:
            donor.orbit_au = (axis, axis)
            return donor
        draw.redrawn("accreting.donor")
    return None


def _white_dwarf_companion(placed, base, binary, age_yr, rng):
    """The companion leaf of a placed white dwarf for the `binary` type `draw_white_dwarf_binary` named (GEN.271), with
    its orbit span, or None when none fits; a cataclysmic variable's white dwarf is magnetic at
    `tuning.CATACLYSMIC_MAGNETIC_SHARE` (a flag on `base`)."""
    if binary == "double-white-dwarf":
        companion = _remnant("white-dwarf", age_yr, 0.0, rng, fe_h=placed.tree_fe_h)
        span = tuning.DOUBLE_WHITE_DWARF_AXIS_AU
    elif binary == "post-common-envelope":
        companion = (_stripped_helium_leaf(1.0, rng) if draw_stripped_helium(rng.random())
                     else _star_leaf(_alive_star_mass(0.1, mass_boundaries.sol("core_collapse_limit"), age_yr, rng, placed.tree_fe_h), age_yr, rng,
                                fe_h=placed.tree_fe_h))
        span = tuning.POST_COMMON_ENVELOPE_AXIS_AU
    else:
        companion = _accreting_companion(placed.mass_solar, binary, age_yr, rng, placed.tree_fe_h)
        if binary == "cataclysmic-variable" and companion is not None and rng.random() < tuning.CATACLYSMIC_MAGNETIC_SHARE:
            base.flags |= MEMBER_MAGNETIC
        span = None
    if companion is not None and span is not None:
        companion.orbit_au = span
    return companion


def _giant_white_dwarf(leaf, binary, age_yr, rng, fe_h=0.0):
    """The white dwarf of a barium or symbiotic giant (GEN.271), born from a star heavier than the giant, on the orbit of
    the type, or None when no white dwarf can have formed by the age. The giant carries the flag."""
    dwarf = _remnant("white-dwarf", age_yr, leaf.initial, rng, fe_h=fe_h)
    if dwarf is None:
        return None
    days = tuning.SYMBIOTIC_PERIOD_DAYS if binary == "symbiotic" else tuning.BARIUM_PERIOD_DAYS
    total = leaf.mass + dwarf.mass
    dwarf.orbit_au = tuple(stability.axis_au(period / 365.25, total) for period in days)
    leaf.flags |= MEMBER_SYMBIOTIC if binary == "symbiotic" else MEMBER_BARIUM
    return dwarf


_DRAW_COMPANION = object()
"""object: The default of `draw_member_tree`'s `companion`: the tree draws the compact object's companion
state itself."""


def companion_floor_solar(fe_h=None):
    """The lightest birth mass of another star of a tree: the hydrogen-burning limit (`mass_boundaries.sol("hydrogen_burning_limit")`),
    or, for a place of birth [Fe/H] `fe_h` whose own limit is higher (GEN.936: 0.0804 at -0.2, 0.09 at -1), that one. A lighter object
    is a brown dwarf (`_star_leaf`) and would use up one of the table's places (GEN.899)."""
    limit = mass_boundaries.sol("hydrogen_burning_limit")
    return limit if fe_h is None else max(limit, stellar_abundances.hydrogen_burning_limit_for_fe_h(fe_h))


def companion_ratio_range(birth, fe_h=None):
    """The range `(low, high)` the birth mass of another star of a tree is drawn over, as a share of the primary's `birth`
    solar masses, uniform: `tuning.BINARY_MASS_RATIO_RANGE`, the low end no lower than puts the star under the floor
    (`companion_floor_solar`): the `star-count` table counts stars, so a lighter object would use up one of its
    places (GEN.899). The brown dwarfs of a tree are a member of their own (`tuning.BROWN_DWARF_COMPANION_CHANCE`); lighter
    objects than they are planets, filled in as planets are, not members of a tree (GEN.596)."""
    low_ratio, high_ratio = tuning.BINARY_MASS_RATIO_RANGE
    return min(max(low_ratio, companion_floor_solar(fe_h) / birth), high_ratio), high_ratio


def _companion_mass(birth, rng, fe_h=None):
    """The initial mass of another star of a tree whose primary was born with `birth` solar masses: `birth * U(low, high)`
    over `companion_ratio_range`."""
    return birth * rng.uniform(*companion_ratio_range(birth, fe_h))


def _is_helium_core(leaf):
    """Whether a roster leaf is a helium-core white dwarf (GEN.996): a white dwarf lighter than the helium flash core, which the
    initial-final relation gives none of (`tuning.WD_MASS_RANGE_SOL`), so only a companion's overflow made it."""
    return leaf.kind == "white-dwarf" and leaf.mass < tuning.HELIUM_FLASH_CORE_MASS_SOL


def _bind_helium_core_companion(dwarf, companion, fe_h):
    """Puts `companion` on the circular orbit at which the giant of the helium-core white dwarf `dwarf` filled its lobe with that
    core (`helium_core_white_dwarf.overflow_axis_au`; the companion's mass is its birth mass when it evolves after the giant, its
    present one when it is a remnant already), and gives the white dwarf the giant's radius at the overflow as its largest."""
    partner = companion.initial if companion.initial <= dwarf.initial else companion.mass
    axis = helium_core_white_dwarf.overflow_axis_au(dwarf.mass, dwarf.initial, partner, fe_h)
    companion.orbit_au, companion.circular = (axis, axis), True
    dwarf.radius_max_au = max(dwarf.radius_au, _m_to_au(helium_core_white_dwarf.giant_radius_rsun(dwarf.mass, fe_h)
                                                          * constants.SOLAR_RADIUS_M))
    dwarf.mass_max_solar = max(dwarf.initial, dwarf.mass)


def _roster(placed, age_yr, rng, max_stars, n_stars, prevalence, companion=_DRAW_COMPANION, bound_runaway=False,
            star_weights=None):
    """The members and the shape units of a tree: `(leaves, units)`, leaf 0 being the placed object."""
    birth = placed.birth_mass_solar
    if placed.kind == "star":
        types = [rng.random() for _ in range(STAR_TYPE_DRAWS)]
        merger = draw_merger_product(birth, age_yr, types[0])
        be, stripped = draw_be_star(birth, age_yr, types[1], types[2])
        effective_age = age_yr * (REJUVENATION_FRACTION_RANGE[0] + types[5] * (
            REJUVENATION_FRACTION_RANGE[1] - REJUVENATION_FRACTION_RANGE[0])) if merger else None
        leaf = (_star_leaf(birth, effective_age or age_yr, rng, lighten=False, fe_h=placed.tree_fe_h) if placed.record is None
                else _leaf_of(placed.record))
        leaf.mass = placed.mass_solar if leaf.kind != "white-dwarf" else leaf.mass
        if placed.record is not None:
            leaf.seed = None                 # the built primary needs no seed of its own: the member draws one
        leaf.effective_age_yr = effective_age
        if merger:
            leaf.flags |= MEMBER_MERGER_PRODUCT
        if be:
            leaf.flags |= MEMBER_BE
        if draw_mass_gainer(birth, age_yr, types[3]):
            leaf.flags |= MEMBER_MASS_GAINER
        exotic = draw_exotic_star(leaf.kind, age_yr, types[4])
        if exotic == "thorne-zytkow":
            leaf.kind, leaf.state = "thorne-zytkow", None
            leaf.flags |= MEMBER_THEORETICAL
            return [leaf], [0]
        if exotic is not None:
            leaf.flags |= MEMBER_MERGER_PRODUCT | (MEMBER_R_CRB if exotic == "r-crb" else MEMBER_FK_COMAE)
        leaves, units = [leaf], [0]
        if n_stars is not None and n_stars <= 1:     # a directive for a single object (no binaries, GEN.719)
            return leaves, units
        compact = (_bound_runaway_companion(placed, age_yr, rng) if bound_runaway
                   else _star_compact_companion(placed, age_yr, rng))
        giant_type = draw_giant_companion(leaf.kind, birth, types[6])
        if compact is None and giant_type is not None:
            compact = _giant_white_dwarf(leaf, giant_type, age_yr, rng, placed.tree_fe_h)
        elif compact is not None and not bound_runaway and placed_class(placed) == "star-massive":
            mode = draw_accreting_compact("high-mass-xray-binary", types[7], types[8], age_yr)
            if mode is not None:
                _feed_compact(leaf, compact, mode, types[9])
        if n_stars is None:
            count = _draw_star_count(birth, max_stars - (compact is not None), prevalence, rng, star_weights)
            if merger == "massive":
                count = max(1, count - 1)
            if stripped:
                count = max(count, 2)
        else:
            count = max(1, min(n_stars, max_stars - (compact is not None)))
        if n_stars is None and compact is None and _is_helium_core(leaf):
            count = max(count, 2)            # the companion that stripped it (GEN.996)
        for _ in range(count - 1):
            leaves.append(_star_leaf(_companion_mass(birth, rng, placed.fe_h), age_yr, rng, fe_h=placed.tree_fe_h))
        # A brown dwarf is no star of the table's count: a drawn tree may hold one more member (a directive's count is
        # of stars, so it adds none).
        brown = (n_stars is None and birth > mass_boundaries.sol("hydrogen_burning_limit") and len(leaves) < max_stars - (compact is not None)
                 and rng.random() < tuning.BROWN_DWARF_COMPANION_CHANCE)
        if brown:
            leaves.append(_star_leaf(brown_dwarf.sample_mass_solar(rng), age_yr, rng, fe_h=placed.tree_fe_h, substellar=True))
        if stripped and count >= 2:
            leaves[1] = _stripped_helium_leaf(birth, rng)
            leaf.flags |= MEMBER_MASS_GAINER
        if compact is not None:
            leaves.append(compact)
            units[0] = (0, len(leaves) - 1)
        if compact is None and count >= 2 and _is_helium_core(leaf):
            # A helium-core white dwarf (GEN.996, `helium_core_draw`) was made by a companion: it is in a pair with the first star
            # of the roster, on the axis at which its giant filled its lobe with this core, on a circle.
            _bind_helium_core_companion(leaf, leaves[1], placed.tree_fe_h)
            units[0] = (0, 1)
            units.extend(range(2, count + brown))
            return leaves, units
        units.extend(range(1, count + brown))
        return leaves, units
    base = {"black-hole": lambda: _remnant_leaf("black-hole", birth, placed.mass_solar, _placed_direct_collapse(placed)),
            "intermediate-mass-black-hole": lambda: _remnant_leaf("black-hole", birth, placed.mass_solar, True),
            "neutron-star": lambda: _remnant_leaf("neutron-star", birth, placed.mass_solar),
            "white-dwarf": lambda: _remnant_leaf("white-dwarf", birth, placed.mass_solar)}[placed.kind]()
    leaves = [base]
    if placed.kind == "white-dwarf" and draw_white_dwarf_merger(age_yr, rng.random()):
        base.flags |= MEMBER_MERGER_PRODUCT
        base.effective_age_yr = rejuvenated_age_yr(age_yr, rng)
    if n_stars is not None and n_stars <= 1:     # a directive for a single object (no binaries, GEN.719)
        return leaves, [0]
    if companion is _DRAW_COMPANION and _is_helium_core(base):
        # A helium-core white dwarf (GEN.996, `helium_core_draw`) was made by its companion: a star of the tree's mass law on the
        # axis its giant filled its lobe at.
        stripper = _star_leaf(_companion_mass(birth, rng, placed.fe_h), age_yr, rng, fe_h=placed.tree_fe_h)
        _bind_helium_core_companion(base, stripper, placed.tree_fe_h)
        return [base, stripper], [(0, 1)]
    state = companion_state(placed, rng) if companion is _DRAW_COMPANION else companion
    if state is None:
        return leaves, [0]
    binary = (draw_white_dwarf_binary(placed_class(placed), state, age_yr, rng.random())
              if placed.kind == "white-dwarf" else None)
    if binary is not None:
        companion = _white_dwarf_companion(placed, base, binary, age_yr, rng)
    elif state == "star" and placed.kind in COMPACT_KINDS and draw_low_mass_xray_pair(placed, age_yr, rng.random()):
        companion = _low_mass_xray_companion(placed, base, age_yr, rng)
    elif state in ("black-widow", "redback"):
        companion = _spider_companion(placed, base, state, age_yr, rng)
    else:
        companion = _compact_pair_companion(placed, state, age_yr, rng)
    if companion is None:
        return leaves, [0]
    leaves.append(companion)
    units = [(0, 1)]
    if (placed.kind in COMPACT_KINDS and max_stars >= 3
            and rng.random() < system_types.SYSTEM_TYPE_RATES["compact-triple"].default):
        third = _compact_tertiary(placed, age_yr, rng)
        if third is not None and sum(leaf.compact for leaf in leaves + [third]) <= MAX_COMPACT_MEMBERS:
            leaves.append(third)
            units.append(2)
    return leaves, units


def _placed_direct_collapse(placed):
    """Whether a placed black hole collapsed directly: its progenitor, when known, else its own mass."""
    if placed.initial_mass_solar is not None:
        return remnant_progenitors.fails(placed.initial_mass_solar)
    return placed.mass_solar >= DIRECT_COLLAPSE_HOLE_SOL


def _star_compact_companion(placed, age_yr, rng):
    """The compact companion of a living star, or None (one draw)."""
    table = MEMBER_STATE_BY_PLACED.get(placed_class(placed), {})
    lifetime_yr = evolution.remnant_formation_time_gy(placed.birth_mass_solar, placed.tree_fe_h) * 1.0e9
    in_window = age_yr >= lifetime_yr * (1.0 - REMNANT_WINDOW_FRACTION)
    pick = rng.random()
    for kind, (default, _low, _high) in table.items():
        chance = default / REMNANT_WINDOW_FRACTION
        if pick < chance:
            return _remnant(kind, age_yr, placed.birth_mass_solar, rng, fe_h=placed.tree_fe_h) if in_window else None
        pick -= chance
    return None


def _bound_runaway_companion(placed, age_yr, rng):
    """The compact companion a supernova left bound to a runaway star (GEN.329): a neutron star or a black hole in the
    proportion of the placed class's table, from a progenitor heavier than the runaway that has had time to collapse
    by `age_yr`; None when none can have."""
    table = MEMBER_STATE_BY_PLACED.get(placed_class(placed), {})
    odds = [(kind, default) for kind, (default, _low, _high) in table.items() if kind in COMPACT_KINDS]
    pick = rng.random() * sum(default for _kind, default in odds)
    for kind, default in odds:
        if pick < default:
            return _remnant(kind, age_yr, placed.birth_mass_solar, rng, fe_h=placed.tree_fe_h)
        pick -= default
    return None


def _compact_pair_companion(placed, state, age_yr, rng):
    """The companion leaf of a compact (or white dwarf) placed object, drawn for `state`."""
    if state in ("star", "brown-dwarf"):
        if state == "brown-dwarf":
            return _star_leaf(rng.uniform(0.02, 0.07), age_yr, rng, fe_h=placed.tree_fe_h)
        low, high = (0.1, 1.0) if placed.recycled else (0.5, 16.0) if placed.kind == "neutron-star" else (0.1, 8.0)
        return _star_leaf(_alive_star_mass(low, high, age_yr, rng, placed.tree_fe_h), age_yr, rng,
                          fe_h=placed.tree_fe_h)
    limit = placed.mass_solar if state == "black-hole" else None
    return _remnant(state, age_yr, 0.0, rng, limit_mass=limit, fe_h=placed.tree_fe_h)


def _compact_tertiary(placed, age_yr, rng):
    """The third member of a compact binary's tree (GEN.272, `compact-triple`'s 30%): a compact object of the placed
    object's own kind at 80% and the other kind at 20%, in `COMPACT_TERTIARY_SHARE` of them (a black hole's lighter than
    the placed one), else a star."""
    if rng.random() < COMPACT_TERTIARY_SHARE:
        own = rng.random() < 0.8
        kind = placed.kind if own else ("neutron-star" if placed.kind == "black-hole" else "black-hole")
        limit = placed.mass_solar if kind == "black-hole" and placed.kind == "black-hole" else None
        return _remnant(kind, age_yr, 0.0, rng, limit_mass=limit, fe_h=placed.tree_fe_h)
    return _star_leaf(_alive_star_mass(0.1, mass_boundaries.sol("core_collapse_limit"), age_yr, rng, placed.tree_fe_h), age_yr, rng, fe_h=placed.tree_fe_h)


def _assemble(leaves, builder, age_yr, dropped, rng):
    """The result for the placed shape in `builder`: members renumbered by roster order."""
    used = sorted({i for record in builder.nodes for i in builder.leaves_of((record["left"], record["right"]))})
    number = {leaf_index: n for n, leaf_index in enumerate(used)}
    where = {}
    nodes = []
    for record in builder.nodes:
        for child_slot, child in enumerate((record["left"], record["right"])):
            if isinstance(child, int):
                where[child] = (record["node_no"], child_slot)
        nodes.append(Node(node_no=record["node_no"], parent_node_no=record["parent_node_no"], slot=record["slot"],
                          mass_solar=record["mass_solar"], separation_au=record["separation_au"],
                          separation_km=_m_to_km(_au_to_m(record["separation_au"])),
                          eccentricity=record["eccentricity"], inclination_deg=record["inclination_deg"],
                          ascending_node_deg=record["ascending_node_deg"],
                          arg_periapsis_deg=record["arg_periapsis_deg"], mean_anomaly_deg=record["mean_anomaly_deg"],
                          period_years=record["period_years"], kind=record["kind"], interaction=record["interaction"],
                          flags=record["flags"]))
    members = []
    for leaf_index in used:
        leaf = leaves[leaf_index]
        node_no, slot = where[leaf_index]
        members.append(Member(member_no=number[leaf_index], parent_node_no=node_no, slot=slot, kind=leaf.kind,
                              state=leaf.state, flags=leaf.flags | builder.extra_flags.get(leaf_index, 0), mass_solar=leaf.mass,
                              initial_mass_solar=leaf.initial, radius_au=leaf.radius_au, seed=leaf.member_seed(rng),
                              effective_age_yr=leaf.effective_age_yr, radius_max_au=leaf.largest_radius_au,
                              mass_max_solar=leaf.largest_mass_solar))
    return MemberTree(tuple(members), tuple(nodes), age_yr, dropped)


def _contact_pair_alone(leaves, builder, n_stars, rng):
    """
    The tree of a contact binary left without its other stars, or None (GEN.270, `contact-binary`'s `tertiary` share,
    60% of contact binaries). The stability screen leaves the inner pair of a hierarchy only short periods, so 84% of the
    contact binaries the draw makes have a tertiary, more than the catalogue's share. A tree with a contact pair that
    holds the placed star keeps its other stars with the odds `tuning.CONTACT_TERTIARY_KEEP_SHARE`; otherwise the pair is
    placed alone, on the orbit it was drawn on with its orientation redrawn (member 0 stays the placed object). Only for
    a tree of stars and a count the draw chose (a forced count is a directive). None when the tree is kept or has no such
    pair, else `(leaves, builder)` of the pair.
    """
    if n_stars is not None or len(leaves) < 3 or any(leaf.compact for leaf in leaves):
        return None
    contact = next((record for record in builder.nodes
                    if record["interaction"] == "contact" and 0 in (record["left"], record["right"])), None)
    if contact is None:
        return None
    if rng.random() < tuning.CONTACT_TERTIARY_KEEP_SHARE:
        return None
    pair = [leaves[0], leaves[contact["right"] if contact["left"] == 0 else contact["left"]]]
    shape = (0, 1)
    fixed = {id(shape): (_yr_to_d(contact["period_years"]), contact["separation_au"], contact["eccentricity"])}
    for _ in range(TREE_RESTARTS):
        alone = _Builder(pair, builder.age_yr, rng, fixed)
        try:
            alone.place(shape, None, 0, None)
        except _Failed:
            draw.redrawn("tree.restart_pair_failed")
            continue
        if alone.clears_extents():
            alone.finish()
            return pair, alone
    return None


def rejuvenated_age_yr(age_yr, rng=draw):
    """A merger product's effective age: a fraction (`REJUVENATION_FRACTION_RANGE`) of the system's age. The
    system's own age stays on the tree."""
    return age_yr * rng.uniform(*REJUVENATION_FRACTION_RANGE)


def _in_sphere(radius, rng):
    """A point drawn uniformly inside a sphere of this radius, as a 3-tuple."""
    while True:
        point = [rng.uniform(-1.0, 1.0) for _ in range(3)]
        if sum(c * c for c in point) <= 1.0:
            return tuple(c * radius for c in point)


def _group_member(number, leaf, rng, parent_node_no, slot, flags, **extra):
    """A `Member` of a young group's tree from a roster leaf."""
    return Member(member_no=number, parent_node_no=parent_node_no, slot=slot, kind=leaf.kind, state=leaf.state,
                  flags=leaf.flags | flags, mass_solar=leaf.mass, initial_mass_solar=leaf.initial,
                  radius_au=leaf.radius_au, seed=leaf.member_seed(rng), effective_age_yr=leaf.effective_age_yr,
                  radius_max_au=leaf.largest_radius_au, mass_max_solar=leaf.largest_mass_solar, **extra)


def _group_tree(leaves, age_yr, decay_age_yr, rng):
    """
    The tree of a young group that has not decayed (GEN.327): one `group` node and its members, each at an offset from the
    group's centre (AU) with a velocity (km/s) from the dispersion `tuning.YOUNG_GROUP_VELOCITY_DISPERSION_KMS`, the centre
    of mass at rest at the origin so the stored point is the barycentre. The node's period is the crossing time (the
    members fly straight; the node has no orbit), its separation the group's radius.
    """
    low, high = tuning.YOUNG_GROUP_RADIUS_AU_RANGE
    radius_au = math.exp(rng.uniform(math.log(low), math.log(high)))
    sigma = tuning.YOUNG_GROUP_VELOCITY_DISPERSION_KMS
    masses = [leaf.mass for leaf in leaves]
    total = sum(masses)
    offsets = [_in_sphere(radius_au, rng) for _ in leaves]
    velocities = [tuple(rng.gauss(0.0, sigma) for _ in range(3)) for _ in leaves]

    def at_rest(vectors):
        centre = [sum(m * v[axis] for m, v in zip(masses, vectors)) / total for axis in range(3)]
        return [tuple(v[axis] - centre[axis] for axis in range(3)) for v in vectors]

    offsets, velocities = at_rest(offsets), at_rest(velocities)
    crossing_yr = radius_au * constants.AU_TO_KM / sigma / constants.SECONDS_PER_YEAR
    node = Node(node_no=1, parent_node_no=None, slot=0, mass_solar=total, separation_au=radius_au,
                separation_km=_m_to_km(_au_to_m(radius_au)), eccentricity=0.0, inclination_deg=0.0,
                ascending_node_deg=0.0, arg_periapsis_deg=0.0, mean_anomaly_deg=0.0, period_years=crossing_yr,
                kind="group", interaction="detached", flags=0)
    members = tuple(_group_member(n, leaf, rng, 1, n, 0, offset_au=offsets[n], velocity_kms=velocities[n],
                                  decay_age_yr=decay_age_yr) for n, leaf in enumerate(leaves))
    return MemberTree(members, (node,), age_yr)


def _decayed_group(leaves, age_yr, rng, max_axis_au, pair_placed):
    """
    The result of a young group older than its decay age (GEN.327): the lightest member (never the placed star) is
    ejected, and the others harden into a hierarchy no wider than `tuning.YOUNG_GROUP_REMNANT_MAX_AXIS_AU`. The ejected
    star is on the tree's `ejected`, flagged `runaway_dynamical`; it is not a row of its own, because the runaway rows the
    scatter draws (`draw_runaway_origin`) already carry the catalogue's share of dynamical runaways.
    """
    lightest = min(range(1, len(leaves)), key=lambda i: leaves[i].mass)
    kept = [leaf for i, leaf in enumerate(leaves) if i != lightest]
    tree = _grow_tree(kept, list(range(len(kept))), age_yr, rng, len(kept),
                      min(max_axis_au, tuning.YOUNG_GROUP_REMNANT_MAX_AXIS_AU), pair_placed)
    thrown = _group_member(0, leaves[lightest], rng, None, 0, MEMBER_RUNAWAY_DYNAMICAL)
    return dataclasses.replace(tree, ejected=(thrown,))


def draw_young_group(placed, age_yr, rng=draw, max_stars=MAX_MEMBERS, prevalence=1.0, bound_runaway=False,
                     max_axis_au=math.inf, pair_placed=False):
    """
    The tree of a massive star born in a young non-hierarchical group (a trapezium, GEN.327), or None when it was not.

    An O or early-B star (`companion_class`) younger than `tuning.YOUNG_GROUP_AGE_MAX_YR` is in a group at the `young-group`
    share of `SYSTEM_TYPE_RATES` (5%). The group holds 3 to 5 stars (`tuning.YOUNG_GROUP_SIZE_WEIGHTS`) and a decay age
    drawn from `tuning.YOUNG_GROUP_DECAY_AGE_RANGE_YR`. A system younger than that has the group itself, a `group` node whose
    members carry an offset and a velocity; a group is never stored past its decay age: an older system has the group's
    result, one star fewer in a tighter hierarchy (`decayed-trapezium-remnant`).

    Args:
        placed (Placed): The star the tree is drawn for.
        age_yr (float): The system's age in years.
        rng: The random source (a `draw.Stream` or the `draw` module). A star that is not eligible draws nothing.
        max_stars, prevalence, bound_runaway, max_axis_au, pair_placed: As for `draw_member_tree`.

    Returns:
        MemberTree or None: None leaves the caller to draw the star's tree as for any other star.
    """
    if (placed.kind != "star" or age_yr >= tuning.YOUNG_GROUP_AGE_MAX_YR
            or companion_class(placed.birth_mass_solar) not in ("O", "B-early")):
        return None
    if rng.random() >= system_types.SYSTEM_TYPE_RATES["young-group"].default:
        return None
    sizes = list(tuning.YOUNG_GROUP_SIZE_WEIGHTS)
    count = min(rng.choices(sizes, weights=list(tuning.YOUNG_GROUP_SIZE_WEIGHTS.values()))[0], max_stars)
    leaves, _units = _roster(placed, age_yr, rng, max_stars, count, prevalence, bound_runaway=bound_runaway)
    if len(leaves) < 3 or any(leaf.compact for leaf in leaves):
        return None
    low, high = tuning.YOUNG_GROUP_DECAY_AGE_RANGE_YR
    decay_age_yr = math.exp(rng.uniform(math.log(low), math.log(high)))
    if age_yr < decay_age_yr:
        return _group_tree(leaves, age_yr, decay_age_yr, rng)
    return _decayed_group(leaves, age_yr, rng, max_axis_au, pair_placed)


def _single(leaf, age_yr, dropped, rng):
    member = Member(0, None, 0, leaf.kind, leaf.state, leaf.flags, leaf.mass, leaf.initial, leaf.radius_au,
                    leaf.member_seed(rng), leaf.effective_age_yr, radius_max_au=leaf.largest_radius_au,
                    mass_max_solar=leaf.largest_mass_solar)
    return MemberTree((member,), (), age_yr, dropped)


def draw_member_tree(placed, age_yr, rng=draw, max_stars=MAX_MEMBERS, n_stars=None, prevalence=1.0,
                     companion=_DRAW_COMPANION, bound_runaway=False, max_axis_au=math.inf, pair_placed=False,
                     groups=False):
    """
    The members of the system around `placed` and their orbits, as a `MemberTree`.

    Args:
        placed (Placed): The scattered or filled object.
        age_yr (float): The system's age in years; states, remnants and the Peters and
            Kozai-Lidov checks read it.
        rng: The random source (a `draw.Stream` or the `draw` module).
        max_stars (int): The most members in the tree (`GEN.267`'s cap; 12).
        n_stars (int or None): Forces the number of stars of a star's tree (a directive). 1 is a single object
            for any placed object, a living star or a white dwarf or other compact object (no companion of any
            kind, GEN.719). Other counts are ignored for a compact or white dwarf placed object, whose companions
            come from `MEMBER_STATE_BY_PLACED`.
        prevalence (float): Scales the chance of five or more stars (1 is the table).
        companion (str or None): For a compact or white dwarf placed object, the companion state the caller's
            screen already drew (`companion_state_at`): the tree takes it instead of drawing its own, `None`
            for no companion. Left out, the tree draws it (`companion_state`). Ignored for a star.
        bound_runaway (bool): The placed star is a supernova-channel runaway that kept a compact companion
            (`draw_runaway_origin`, GEN.329): the tree gives it one, a neutron star or a black hole that has had time
            to form by the age, in place of the draw at the odds of `MEMBER_STATE_BY_PLACED`. Ignored for any other
            placed object.
        max_axis_au (float): The widest axis any orbit of the tree may have: what the place keeps of a wide pair against
            passing stars and the tide (`place_terms.wide_axis_limit_au`, GEN.103). The scatter, which has no place,
            leaves it unbounded.
        pair_placed (bool): Draw the shape so the placed object is in a pair of two single members, when the tree has one
            (a system build, whose pair is the primary and secondary); the scatter leaves the shape as drawn.
        groups (bool): A massive star younger than `tuning.YOUNG_GROUP_AGE_MAX_YR` is, at the `young-group` share, born in a
            young group (`draw_young_group`, GEN.327): a `group` node of stars flying straight while the system is younger
            than the group's decay age, the hierarchy the group left once it is older. A system build cannot hold a group,
            so only the scatter asks for it.

    Returns:
        MemberTree: The placed object is member 0. A single object has one member and no nodes.
    """
    if groups and placed.kind == "star":
        group = draw_young_group(placed, age_yr, rng, max_stars, prevalence, bound_runaway, max_axis_au, pair_placed)
        if group is not None:
            return group
    weights, inside_weight = close_trend(placed.birth_mass_solar, placed.fe_h, prevalence) if placed.kind == "star" else (None, 1.0)
    leaves, units = _roster(placed, age_yr, rng, max_stars, n_stars, prevalence, companion, bound_runaway, weights)
    return _grow_tree(leaves, units, age_yr, rng, n_stars, max_axis_au, pair_placed, inside_weight, placed.tree_fe_h)


def _grow_tree(leaves, units, age_yr, rng, n_stars, max_axis_au, pair_placed, inside_weight=1.0, fe_h=0.0):
    """The `MemberTree` of the `leaves` in the shape `units` names (the placed leaf first): hierarchies are drawn until one
    clears the stability and extent screens, and when none does the lightest loose member is dropped and the draw tries again."""
    dropped = 0

    def shape_of(units, rng):
        shape = _shape(units, rng)
        return _pair_placed(shape) if pair_placed else shape

    while True:
        if len(units) == 1 and isinstance(units[0], int):
            return _single(leaves[units[0]], age_yr, dropped, rng)
        for _ in range(TREE_RESTARTS):
            builder = _Builder(leaves, age_yr, rng, max_axis_au=max_axis_au, inside_weight=inside_weight, fe_h=fe_h)
            try:
                builder.place(shape_of(units, rng), None, 0, None)
            except _Failed:
                draw.redrawn("tree.restart_pair_failed")
                continue
            if not builder.clears_extents():
                draw.redrawn("tree.restart_overlap")
                continue
            builder.finish()
            alone = _contact_pair_alone(builder.leaves, builder, n_stars, rng)
            if alone is not None:
                builder = alone[1]
            return _assemble(builder.leaves, builder, age_yr, dropped, rng)
        loose = [u for u in units[1:] if isinstance(u, int)]
        if loose:
            drop = min(loose, key=lambda i: leaves[i].mass)
            units.remove(drop)
        elif len(units) > 1:
            units.pop()
        else:
            units[0] = units[0][0]
        draw.redrawn("tree.member_dropped")
        dropped += 1
