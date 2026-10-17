# planetgen/generation/multistar.py

"""
The member-tree generator (GEN.261): the stars, brown dwarfs, remnants and
protostars that share a system, in a stable hierarchy of up to twelve
members, drawn for one placed object (a scattered star or compact object, or
any object the fill creates). One function, `draw_member_tree`, with two call
sites: the scatter passes and the fill. It draws the members and their orbits
only: planets, moons and belts stay at sector fill.

The laws and numbers are in docs/design/multi-star-systems.md (sections 4, 6
and 13) and docs/design/multistar-and-compact-systems.md (sections 3, 4, 5
and 7); the stability criteria, the Kozai-Lidov screen and the Peters time are
`planetgen.physics.stability`. In short:

1. **How many stars.** `STAR_COUNT_PROBABILITIES` by the placed star's
   spectral class (O to M, from its initial mass), a tail for 8 to 12 stars
   that continues the table's last column, and the `prevalence` knob on the
   five-or-more columns. A forced `n_stars` or a lower `max_stars` wins.
2. **Which members.** The other stars have initial masses `m1 * U(0.1, 1)`
   and take their present state from the system's age (main sequence, giant,
   white dwarf, pre-main-sequence, protostar; a brown dwarf below 0.08
   solar masses). A compact companion (white dwarf, neutron star, black hole)
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
from planetgen.physics import constants, pre_main_sequence, roche, stability
from planetgen.physics import stellar_evolution as evolution
from planetgen.util import draw

STAR_COUNT_PROBABILITIES = {
    "M": (0.73000, 0.23500, 0.02800, 0.00560, 0.00123, 0.00015, 0.00001),
    "K": (0.60000, 0.31000, 0.06500, 0.02002, 0.00440, 0.00053, 0.00005),
    "G": (0.54500, 0.33500, 0.08000, 0.03202, 0.00705, 0.00085, 0.00008),
    "F": (0.50000, 0.34000, 0.11000, 0.04003, 0.00881, 0.00106, 0.00011),
    "A": (0.40000, 0.40000, 0.15000, 0.03741, 0.01048, 0.00189, 0.00023),
    "B": (0.25000, 0.40000, 0.22000, 0.08962, 0.03137, 0.00784, 0.00118),
    "B-mid": (0.24000, 0.36000, 0.25143, 0.10241, 0.03585, 0.00896, 0.00135),
    "B-early": (0.16000, 0.32000, 0.32686, 0.13314, 0.04660, 0.01165, 0.00175),
    "O": (0.10000, 0.25000, 0.30000, 0.21354, 0.09610, 0.03363, 0.00673),
}
"""dict: The chance of 1 to 7 stars to a system, by the placed star's spectral
class (multistar-and-compact-systems.md 3.3). The B-mid and B-early rows (GEN.449) take the single and binary fractions
of `tuning.MULTIPLICITY_BY_PRIMARY_MASS_BIN` and split the rest across 3 to 7 stars in the proportions of the B row."""

SPECTRAL_CLASS_LOWER_MASS_SOL = (("O", 16.0), ("B", 2.1), ("A", 1.4), ("F", 1.04), ("G", 0.8), ("K", 0.45), ("M", 0.0))
"""tuple: The initial mass at which each spectral class begins, heaviest first."""

TAIL_RATIOS = (0.12, 0.10)
"""tuple: The chance a system of 7 stars has an 8th (0.12) and the chance each
further star follows the last (0.10), up to `max_stars` (multi-star-systems.md
13.1). The tail comes out of the 7-star column, so the table's totals stand."""

PERIOD_LAWS_LOG10_DAYS = {"solar": (5.03, 2.28), "massive": (3.5, 1.5)}
"""dict: The mean and sigma of the log-normal in log10 of the period in days: the
solar-type law (Raghavan 2010) and, for a primary of `MASSIVE_PRIMARY_MIN_MASS_SOL`
or more, the law that peaks at log P = 3.5 (Moe and Di Stefano 2017); the sigma of
the second is a design default."""

MASSIVE_PRIMARY_MIN_MASS_SOL = 2.1
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

ACCRETOR_KINDS = ("white-dwarf", "neutron-star", "black-hole")
"""tuple: The member kinds that receive mass from a lobe-filling partner (an `accreting` pair)."""

INTERACTIONS = ("detached", "semi-detached", "contact", "accreting", "post-common-envelope", "merger-product")
"""tuple: The `Node.interaction` values; the classifier gives the first four and `post-common-envelope`
from geometry, `merger-product` comes from the merger draw of the close-binary types."""

REJUVENATION_FRACTION_RANGE = (0.1, 0.6)
"""tuple: A merger product's effective age as a fraction of the system's age (a design default)."""

FRAGILE_AXIS_AU = 10_000.0
"""float: A root pair wider than this is ultra-wide and flagged `fragile` (passing stars unbind it)."""

PAIR_RETRIES = 60
TREE_RESTARTS = 30
MAX_COMPACT_MEMBERS = 3
WIDE_COMPACT_AXIS_AU = 1000.0
"""float: A bound companion beyond this of a neutron star or a black hole that did
not collapse directly does not survive the collapse."""

DIRECT_COLLAPSE_PROGENITOR_SOL = 40.0
DIRECT_COLLAPSE_HOLE_SOL = 20.0
"""float: A black hole counts as a direct collapse when its progenitor was
`DIRECT_COLLAPSE_PROGENITOR_SOL` or more, or, when the progenitor is not known, when the hole is
`DIRECT_COLLAPSE_HOLE_SOL` or more (a stand-in; the note leaves it open)."""

BLACK_HOLE_TERTIARY_CHANCE = 0.30
COMPACT_TERTIARY_SHARE = 0.10
"""float: A black hole with a companion has a third member in 30% of systems, which is
a neutron star or a black hole in 10% of those (a design default)."""

BROWN_DWARF_MAX_MASS_SOL = 0.08
BROWN_DWARF_DESERT_AU = 5.0
BROWN_DWARF_DESERT_KEEP = 0.01
"""A brown dwarf inside `BROWN_DWARF_DESERT_AU` of a star is kept in 1% of draws."""

BROWN_DWARF_RADIUS_SOL = 0.1
NEUTRON_STAR_RADIUS_KM = 12.0
REMNANT_WINDOW_FRACTION = 1.0 / 3.0
"""float: A remnant companion of a living star exists only in the last third of the star's
life, so the table value (a mean over all ages) is divided by this inside the window."""

NEUTRON_STAR_PROGENITOR_SOL = (8.0, 30.0)
BLACK_HOLE_PROGENITOR_SOL = (20.0, 100.0)
WHITE_DWARF_PROGENITOR_SOL = (0.8, 8.0)

COMPACT_KINDS = ("neutron-star", "black-hole")
MEMBER_KINDS = ("main-sequence", "giant", "white-dwarf", "neutron-star", "black-hole", "brown-dwarf",
                "pre-main-sequence", "protostar", "stripped-helium", "thorne-zytkow", "wolf-rayet")
"""tuple: Every member kind. The last three exist for the type generators that follow (the evolution of a
Wolf-Rayet star, a stripped helium star and a Thorne-Zytkow object is theirs); `star_labels` names them."""
PLACED_KINDS = ("star", "black-hole", "intermediate-mass-black-hole", "neutron-star", "white-dwarf")

MEMBER_DIRECT_COLLAPSE = 1
MEMBER_MASS_REVERSED = 2
MEMBER_MASS_GAINER = 4
MEMBER_MERGER_PRODUCT = 8
MEMBER_BE = 16
MEMBER_FRAGILE = 32
MEMBER_THEORETICAL = 64
MEMBER_FLAGS = {"direct_collapse": MEMBER_DIRECT_COLLAPSE, "mass_reversed": MEMBER_MASS_REVERSED,
                "mass_gainer": MEMBER_MASS_GAINER, "merger_product": MEMBER_MERGER_PRODUCT, "be": MEMBER_BE,
                "fragile": MEMBER_FRAGILE, "theoretical": MEMBER_THEORETICAL}
"""dict: The bits of `Member.flags`. The classifier sets `mass_gainer`, `fragile` and `theoretical`;
`mass_reversed`, `merger_product` and `be` are set by the close-binary types, and `runaway_origin` by the
runaway origin draw (its bits are reserved there)."""
NODE_INTERACTED = 1
NODE_KL_ACTIVE = 2

MEMBER_STATE_BY_PLACED = {
    "star-massive": {"neutron-star": (0.002, 0.00033, 0.0133), "black-hole": (0.001, 0.00017, 0.0067)},
    "star-intermediate": {"white-dwarf": (0.01, 0.002, 0.05), "neutron-star": (0.0003, 0.00005, 0.0015),
                          "black-hole": (0.0002, 0.00005, 0.0015)},
    "black-hole": {"star": (0.010, 0.001, 0.022), "black-hole": (0.002, 0.00017, 0.024),
                   "neutron-star": (0.0005, 0.000075, 0.0007)},
    "intermediate-mass-black-hole": {"star": (0.01, 0.001, 0.05)},
    "neutron-star": {"star": (0.02, 0.01, 0.03), "neutron-star": (1.8e-4, 6e-5, 4.3e-4),
                     "black-hole": (2.5e-5, 7e-6, 7e-5)},
    "recycled-pulsar": {"white-dwarf": (0.62, 0.50, 0.70), "star": (0.12, 0.05, 0.20),
                        "brown-dwarf": (0.04, 0.01, 0.08), "neutron-star": (2.0e-3, 6e-4, 4.3e-3),
                        "black-hole": (2.5e-5, 7e-6, 7e-5)},
    "white-dwarf": {"star": (0.25, 0.20, 0.30), "white-dwarf": (0.05, 0.03, 0.10)},
    "planetary-nebula-star": {"star": (0.15, 0.10, 0.20), "white-dwarf": (0.05, 0.03, 0.10)},
}
"""dict: Per class of placed object, the chance of a companion of each state as
`(default, low, high)`, the mean over all ages and per object
(multi-star-systems.md 4.2 to 4.4). A star class (`star-massive` for an initial mass of 14
solar masses or more, `star-intermediate` for 2 to 14) adds a compact companion to its
tree of stars, inside the last third of its life; the other classes decide whether the
placed object has a companion at all, and which one. Other placed stars get no compact
companion here."""


@dataclass(frozen=True)
class Placed:
    """The object a tree is drawn for. `mass_solar` is its present mass; `initial_mass_solar` its birth mass
    (defaults to the present mass). `recycled` marks a millisecond pulsar; `planetary_nebula` a
    planetary-nebula central star (a white dwarf)."""
    kind: str
    mass_solar: float
    initial_mass_solar: float | None = None
    recycled: bool = False
    planetary_nebula: bool = False

    @property
    def birth_mass_solar(self):
        return self.mass_solar if self.initial_mass_solar is None else self.initial_mass_solar


@dataclass(frozen=True)
class Member:
    """One member of a tree: member 0 is the placed object. `radius_au` is its extent; `flags` are the bits of
    `MEMBER_FLAGS`; `effective_age_yr` is a merger product's rejuvenated age (None: the system's age). `state` is the
    Yerkes class of a star (`V`, `III`, `VII`, ...), `PMS` or `protostar` for a young one, and None for a remnant or
    a brown dwarf; it is stored in `scatter_members.state`, a VARCHAR(16)."""
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
    number of members the generator could not fit and left out."""
    members: tuple
    nodes: tuple
    age_yr: float
    dropped: int = 0

    @property
    def member_count(self):
        return len(self.members)

    @property
    def system_mass_solar(self):
        return sum(member.mass_solar for member in self.members)

    @property
    def star_count(self):
        return sum(member.kind not in COMPACT_KINDS for member in self.members)

    @property
    def compact_count(self):
        return sum(member.kind in COMPACT_KINDS for member in self.members)


@dataclass
class _Leaf:
    kind: str
    state: str | None
    mass: float
    initial: float
    radius_au: float
    flags: int = 0

    @property
    def compact(self):
        return self.kind in COMPACT_KINDS

    @property
    def direct_collapse(self):
        return bool(self.flags & MEMBER_DIRECT_COLLAPSE)


def spectral_class(initial_mass_sol):
    """The spectral class (O to M) of a star of this initial mass."""
    for letter, lower in SPECTRAL_CLASS_LOWER_MASS_SOL:
        if initial_mass_sol >= lower:
            return letter
    return "M"


def companion_class(initial_mass_sol, letter=None):
    """The key into `STAR_COUNT_PROBABILITIES` and `tuning.BINARY_SYSTEM_PROBABILITY_BY_SPECTRAL_CLASS` of a star: its
    spectral letter (from the mass when `letter` is None), with the B stars of 5 to 16 solar masses split into B-mid and
    B-early (Moe and Di Stefano's bins, GEN.449)."""
    letter = letter or spectral_class(initial_mass_sol)
    if letter == "B":
        for name, lower in tuning.B_SUBCLASS_LOWER_MASS_SOL:
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
    """The weights of 1 to 7 stars for a star of this initial mass, the five-or-more columns scaled by
    `prevalence`."""
    row = STAR_COUNT_PROBABILITIES[companion_class(initial_mass_sol)]
    return tuple(weight * (prevalence if count >= 5 else 1.0) for count, weight in enumerate(row, start=1))


def _draw_star_count(initial_mass_sol, max_stars, prevalence, rng):
    count = rng.choices(range(1, 8), star_count_weights(initial_mass_sol, prevalence))[0]
    ratio = TAIL_RATIOS[0]
    while count >= 7 and count < max_stars and rng.random() < min(1.0, ratio * prevalence):
        count += 1
        ratio = TAIL_RATIOS[1]
    return min(count, max_stars)


def _radius_au(state):
    if state["radius_km"] is not None:
        return state["radius_km"] * constants.KM_TO_M_FACTOR / constants.AU_M
    radius_sol = math.sqrt(state["luminosity_sol"]) * (tuning.SUN_EFFECTIVE_TEMPERATURE_K / state["temperature_k"]) ** 2
    return radius_sol * constants.SOLAR_RADIUS_M / constants.AU_M


def _star_leaf(initial_mass, age_yr, rng, lighten=True):
    """A star-like member of this birth mass at the system's age."""
    if initial_mass < BROWN_DWARF_MAX_MASS_SOL:
        return _Leaf("brown-dwarf", None, initial_mass, initial_mass,
                     BROWN_DWARF_RADIUS_SOL * constants.SOLAR_RADIUS_M / constants.AU_M)
    young = pre_main_sequence.state(initial_mass, age_yr)
    if young is not None:
        # GEN.242: still contracting; the radius follows from the track's luminosity and temperature.
        kind = "protostar" if young["yerkes_class"] == pre_main_sequence.PROTOSTAR else "pre-main-sequence"
        return _Leaf(kind, young["yerkes_class"], initial_mass, initial_mass, _radius_au(young))
    mass = initial_mass
    while True:
        state = evolution.evolve_star(mass, age_yr / 1.0e9, rng)
        if state is not None or not lighten:
            break
        mass *= 0.8
    if state is None:
        raise ValueError(f"a star of {initial_mass:g} solar masses has collapsed by {age_yr:g} years")
    yerkes = state["yerkes_class"]
    kind = "main-sequence" if yerkes == "V" else "white-dwarf" if yerkes == "VII" else "giant"
    return _Leaf(kind, yerkes, state["mass_sol"], mass, _radius_au(state))


def _schwarzschild_au(mass_solar):
    return (2.0 * constants.G * mass_solar * constants.SOLAR_MASS_TO_KG / constants.SPEED_OF_LIGHT_M_S ** 2
            / constants.AU_M)


def _remnant_leaf(kind, progenitor, mass, direct=False):
    if kind == "white-dwarf":
        radius = evolution.white_dwarf_radius_km(mass) * constants.KM_TO_M_FACTOR / constants.AU_M
        return _Leaf(kind, "VII", mass, progenitor, radius)
    if kind == "neutron-star":
        return _Leaf(kind, None, mass, progenitor, NEUTRON_STAR_RADIUS_KM * constants.KM_TO_M_FACTOR / constants.AU_M)
    return _Leaf(kind, None, mass, progenitor, _schwarzschild_au(mass), MEMBER_DIRECT_COLLAPSE if direct else 0)


def _youngest_progenitor_sol(age_yr):
    """The lightest star whose giant phases have ended by this age (a star lighter than
    `HELIUM_FLASH_MASS_SOL` also burns helium on the horizontal branch and climbs the asymptotic giant branch)."""
    def mass(fraction):
        return (tuning.SOLAR_MS_LIFESPAN_GY * fraction * 1.0e9 / age_yr) ** (1.0 / -tuning.MS_LIFESPAN_MASS_EXPONENT)
    light = mass(evolution.remnant_age_fraction(0.0))
    if light < tuning.HELIUM_FLASH_MASS_SOL:
        return light
    return max(mass(tuning.GIANT_PHASE_END_MS_FRACTION), tuning.HELIUM_FLASH_MASS_SOL)


def _remnant(kind, age_yr, lowest, rng, limit_mass=None):
    """A remnant of this kind that has had time to form by `age_yr`, from a progenitor no lighter than
    `lowest`, or None when no progenitor can have collapsed yet."""
    low, high = {"neutron-star": NEUTRON_STAR_PROGENITOR_SOL, "black-hole": BLACK_HOLE_PROGENITOR_SOL,
                 "white-dwarf": WHITE_DWARF_PROGENITOR_SOL}[kind]
    low = max(low, lowest, _youngest_progenitor_sol(age_yr))
    if low >= high:
        return None
    progenitor = evolution.sample_imf_mass_sol(low, high, rng)
    if kind == "white-dwarf":
        mass = min(max(tuning.WD_IFMR_SLOPE * progenitor + tuning.WD_IFMR_INTERCEPT_SOL, tuning.WD_MASS_RANGE_SOL[0]),
                   tuning.WD_MASS_RANGE_SOL[1])
        return _remnant_leaf(kind, progenitor, mass)
    mass_range = tuning.NEUTRON_STAR_MASS_RANGE_SOLAR if kind == "neutron-star" else tuning.BLACK_HOLE_MASS_RANGE_SOLAR
    mass = rng.uniform(*mass_range)
    if limit_mass is not None:
        mass = min(mass, limit_mass)
    return _remnant_leaf(kind, progenitor, mass, direct=kind == "black-hole" and progenitor >= DIRECT_COLLAPSE_PROGENITOR_SOL)


def _alive_star_mass(low, high, age_yr, rng):
    """An IMF mass in `[low, high]` that is still on the main sequence at this age (the heaviest such mass
    caps the range)."""
    cap = (tuning.SOLAR_MS_LIFESPAN_GY * 1.0e9 / age_yr) ** (1.0 / -tuning.MS_LIFESPAN_MASS_EXPONENT)
    high = max(min(high, cap), low * 1.0001)
    return evolution.sample_imf_mass_sol(low, high, rng)


def _is_direct_collapse(leaf):
    if leaf.kind != "black-hole":
        return False
    return leaf.direct_collapse


class _Failed(Exception):
    """A pair could not be placed in `PAIR_RETRIES` tries; the tree restarts."""


class _Builder:
    """Draws the orbits of one tree shape, top down."""

    def __init__(self, leaves, age_yr, rng):
        self.leaves = leaves
        self.age_yr = age_yr
        self.rng = rng
        self.nodes = []
        self.extra_flags = {}

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
        return "star"

    def pair_leaves(self, left, right):
        """The two leaves of a pair of single members, else None."""
        if isinstance(left, int) and isinstance(right, int):
            return self.leaves[left], self.leaves[right]
        return None

    def axis_floor(self, left, right, extent):
        """The closest the pair may be: the Roche-lobe limit for two single members, else the extent rule."""
        pair = self.pair_leaves(left, right)
        if pair is None:
            return EXTENT_AXIS_FACTOR * extent
        a, b = pair
        return max(roche.separation_for_filling(a.radius_au, MAX_FILL_FACTOR, a.mass, b.mass),
                   roche.separation_for_filling(b.radius_au, MAX_FILL_FACTOR, b.mass, a.mass))

    def clears(self, left, right, axis, e, extent):
        """Whether the pair's closest approach keeps its members inside their limits."""
        pair = self.pair_leaves(left, right)
        if pair is None:
            return axis * (1.0 - e) >= PERICENTRE_EXTENT_FACTOR * extent
        a, b = pair
        peri = axis * (1.0 - e)
        return (roche.filling_factor(a.radius_au, peri, a.mass, b.mass) <= MAX_FILL_FACTOR
                and roche.filling_factor(b.radius_au, peri, b.mass, a.mass) <= MAX_FILL_FACTOR)

    def interaction(self, left, right, axis, e):
        """The pair's `interaction` and the flags it puts on the two single members: `(interaction, {leaf: flag})`."""
        pair = self.pair_leaves(left, right)
        if pair is None:
            return "detached", {}
        a, b = pair
        peri = axis * (1.0 - e)
        fill_a = roche.filling_factor(a.radius_au, peri, a.mass, b.mass)
        fill_b = roche.filling_factor(b.radius_au, peri, b.mass, a.mass)
        accretor_a, accretor_b = a.kind in ACCRETOR_KINDS, b.kind in ACCRETOR_KINDS
        if accretor_a != accretor_b:
            donor_fill = fill_b if accretor_a else fill_a
            if donor_fill >= 1.0:
                return "accreting", {}
            if (a.kind == "white-dwarf" or b.kind == "white-dwarf") and axis < tuning.WD_PROGENITOR_ENGULFMENT_AU:
                return "post-common-envelope", {}
            return "detached", {}
        if accretor_a:
            return "detached", {}
        if fill_a >= 1.0 and fill_b >= 1.0:
            return "contact", {left: MEMBER_MASS_GAINER, right: MEMBER_MASS_GAINER}
        if fill_a >= 1.0:
            return "semi-detached", {right: MEMBER_MASS_GAINER}
        if fill_b >= 1.0:
            return "semi-detached", {left: MEMBER_MASS_GAINER}
        return "detached", {}

    def wide_cap(self, sub):
        """The widest axis the pair can have: 1,000 AU if it holds a neutron star or a black hole that did not
        collapse directly, else unbounded."""
        if any(self.leaves[i].compact and not _is_direct_collapse(self.leaves[i]) for i in self.leaves_of(sub)):
            return WIDE_COMPACT_AXIS_AU
        return math.inf

    def period_days(self, law, mass, a_lo, a_hi, massive):
        lo = 365.25 * math.sqrt(a_lo ** 3 / mass)
        hi = 365.25 * math.sqrt(a_hi ** 3 / mass) if math.isfinite(a_hi) else math.inf
        u = self.rng.random()
        if law == "compact-pair":
            lo = max(lo, 365.25 * math.sqrt(COMPACT_PAIR_AXIS_AU[0] ** 3 / mass))
            hi = min(hi, 365.25 * math.sqrt(COMPACT_PAIR_AXIS_AU[1] ** 3 / mass))
        elif law == "interacted":
            lo, hi = max(lo, INTERACTED_PERIOD_DAYS[0]), min(hi, INTERACTED_PERIOD_DAYS[1])
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
        return 10.0 ** dist.inv_cdf(min(max(p_lo + u * (p_hi - p_lo), 1.0e-15), 1.0 - 1.0e-15))

    def eccentricity(self, law, period_d, a_au, root):
        u = self.rng.random()
        if (law == "compact-pair" and a_au < COMPACT_PAIR_CIRCULAR_BELOW_AU) or period_d < CIRCULAR_BELOW_DAYS:
            return 0.0
        cap = min(1.0 - (period_d / 2.0) ** (-2.0 / 3.0),
                  ROOT_ECCENTRICITY_CAP if root else INNER_ECCENTRICITY_CAP)
        return cap * math.sqrt(u)

    def place(self, sub, parent, slot, parent_no):
        """Places the orbit of the pair `sub`, then its child pairs. `parent` is
        `(axis, eccentricity, inclination, node, sibling_mass)` or None for the root."""
        left, right = sub
        if self.mass(left) < self.mass(right):
            left, right = right, left
        m_left, m_right = self.mass(left), self.mass(right)
        total = m_left + m_right
        law = self.law(left, right)
        extent = self.extent(left) + self.extent(right)
        a_lo = self.axis_floor(left, right, extent)
        cap = self.wide_cap(sub)
        heavy = max(self.leaves[i].initial for i in self.leaves_of(sub))
        massive = heavy >= MASSIVE_PRIMARY_MIN_MASS_SOL
        for _ in range(PAIR_RETRIES):
            rng = self.rng
            inclination = math.acos(1.0 - 2.0 * rng.random())
            node = 2.0 * math.pi * rng.random()
            arg = 2.0 * math.pi * rng.random()
            anomaly = 2.0 * math.pi * rng.random()
            mutual = None
            a_hi = ROOT_MAX_AXIS_AU
            if parent is not None:
                mutual = stability.mutual_inclination(inclination, node, parent[2], parent[3])
                a_hi = stability.max_inner_axis(m_left, m_right, parent[4], parent[0], 0.0, parent[1], mutual)
            period_d = self.period_days(law, total, a_lo, min(a_hi, cap), massive)
            if period_d is None:
                continue
            axis = stability.axis_au(period_d / 365.25, total)
            e = self.eccentricity(law, period_d, axis, parent is None)
            if not self.clears(left, right, axis, e, extent):
                continue
            flags = NODE_INTERACTED if law != "star" else 0
            if parent is not None:
                if not stability.is_triple_stable(m_left, m_right, parent[4], axis, parent[0], e, parent[1], mutual):
                    continue
                verdict = stability.kozai_verdict(m_left, m_right, parent[4], axis, parent[0], e, parent[1], mutual,
                                                  self.age_yr, extent)
                if verdict == "unsafe":
                    continue
                if verdict == "kl_active":
                    flags |= NODE_KL_ACTIVE
            if law == "compact-pair" and not stability.survives_gravitational_waves(m_left, m_right, axis, e, self.age_yr):
                continue
            if self._in_desert(left, right, axis):
                continue
            break
        else:
            raise _Failed
        interaction, member_flags = self.interaction(left, right, axis, e)
        for leaf_index, flag in member_flags.items():
            self.extra_flags[leaf_index] = self.extra_flags.get(leaf_index, 0) | flag
        node_no = len(self.nodes) + 1
        kind = "compact" if law != "star" else "close" if axis < 50.0 else "wide"
        record = {"node_no": node_no, "parent_node_no": parent_no, "slot": slot, "mass_solar": total,
                  "separation_au": axis, "eccentricity": e, "inclination_deg": math.degrees(inclination),
                  "ascending_node_deg": math.degrees(node), "arg_periapsis_deg": math.degrees(arg),
                  "mean_anomaly_deg": math.degrees(anomaly), "period_years": period_d / 365.25, "kind": kind,
                  "flags": flags, "left": left, "right": right, "sub": sub,
                  "interaction": interaction}
        self.nodes.append(record)
        for child_slot, (child, sibling) in enumerate(((left, right), (right, left))):
            if not isinstance(child, int):
                self.place(child, (axis, e, inclination, node, self.mass(sibling)), child_slot, node_no)
        return node_no

    def _in_desert(self, left, right, axis):
        if axis >= BROWN_DWARF_DESERT_AU or not (isinstance(left, int) and isinstance(right, int)):
            return False
        kinds = {self.leaves[left].kind, self.leaves[right].kind}
        if "brown-dwarf" not in kinds or kinds == {"brown-dwarf"}:
            return False
        return self.rng.random() >= BROWN_DWARF_DESERT_KEEP

    def clears_extents(self):
        """The final check: every pair's closest approach clears what it holds (the Roche limit for two single
        members, the apocentres of the inner pairs otherwise)."""
        by_sub = {id(record["sub"]): record for record in self.nodes}
        for record in self.nodes:
            extent = 0.0
            for child in (record["left"], record["right"]):
                if isinstance(child, int):
                    extent += self.leaves[child].radius_au
                else:
                    inner = by_sub[id(child)]
                    extent += inner["separation_au"] * (1.0 + inner["eccentricity"])
            if not self.clears(record["left"], record["right"], record["separation_au"], record["eccentricity"],
                               extent):
                return False
        return True

    def finish(self):
        """Flags the members of an ultra-wide system `fragile` and the theoretical kinds `theoretical`."""
        fragile = bool(self.nodes) and self.nodes[0]["separation_au"] > FRAGILE_AXIS_AU
        for index, leaf in enumerate(self.leaves):
            if fragile:
                self.extra_flags[index] = self.extra_flags.get(index, 0) | MEMBER_FRAGILE
            if leaf.kind == "thorne-zytkow":
                self.extra_flags[index] = self.extra_flags.get(index, 0) | MEMBER_THEORETICAL


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


_DRAW_COMPANION = object()
"""object: The default of `draw_member_tree`'s `companion`: the tree draws the compact object's companion
state itself."""


def _roster(placed, age_yr, rng, max_stars, n_stars, prevalence, companion=_DRAW_COMPANION):
    """The members and the shape units of a tree: `(leaves, units)`, leaf 0 being the placed object."""
    birth = placed.birth_mass_solar
    if placed.kind == "star":
        leaf = _star_leaf(birth, age_yr, rng, lighten=False)
        leaf.mass = placed.mass_solar if leaf.kind != "white-dwarf" else leaf.mass
        leaves, units = [leaf], [0]
        compact = _star_compact_companion(placed, age_yr, rng)
        if n_stars is None:
            count = _draw_star_count(birth, max_stars - (compact is not None), prevalence, rng)
        else:
            count = max(1, min(n_stars, max_stars - (compact is not None)))
        for _ in range(count - 1):
            leaves.append(_star_leaf(birth * rng.uniform(0.1, 1.0), age_yr, rng))
        if compact is not None:
            leaves.append(compact)
            units[0] = (0, len(leaves) - 1)
        units.extend(range(1, count))
        return leaves, units
    base = {"black-hole": _remnant_leaf("black-hole", birth, placed.mass_solar, _placed_direct_collapse(placed)),
            "intermediate-mass-black-hole": _remnant_leaf("black-hole", birth, placed.mass_solar, True),
            "neutron-star": _remnant_leaf("neutron-star", birth, placed.mass_solar),
            "white-dwarf": _remnant_leaf("white-dwarf", birth, placed.mass_solar)}[placed.kind]
    leaves = [base]
    state = companion_state(placed, rng) if companion is _DRAW_COMPANION else companion
    if state is None:
        return leaves, [0]
    companion = _compact_pair_companion(placed, state, age_yr, rng)
    if companion is None:
        return leaves, [0]
    leaves.append(companion)
    units = [(0, 1)]
    if placed.kind == "black-hole" and rng.random() < BLACK_HOLE_TERTIARY_CHANCE and max_stars >= 3:
        third = _black_hole_tertiary(placed, age_yr, rng)
        if third is not None and sum(leaf.compact for leaf in leaves + [third]) <= MAX_COMPACT_MEMBERS:
            leaves.append(third)
            units.append(2)
    return leaves, units


def _placed_direct_collapse(placed):
    """Whether a placed black hole collapsed directly: its progenitor, when known, else its own mass."""
    if placed.initial_mass_solar is not None:
        return placed.initial_mass_solar >= DIRECT_COLLAPSE_PROGENITOR_SOL
    return placed.mass_solar >= DIRECT_COLLAPSE_HOLE_SOL


def _star_compact_companion(placed, age_yr, rng):
    """The compact companion of a living star, or None (one draw)."""
    table = MEMBER_STATE_BY_PLACED.get(placed_class(placed), {})
    lifetime_yr = (evolution.main_sequence_lifetime_gy(placed.birth_mass_solar)
                   * evolution.remnant_age_fraction(placed.birth_mass_solar) * 1.0e9)
    in_window = age_yr >= lifetime_yr * (1.0 - REMNANT_WINDOW_FRACTION)
    pick = rng.random()
    for kind, (default, _low, _high) in table.items():
        chance = default / REMNANT_WINDOW_FRACTION
        if pick < chance:
            return _remnant(kind, age_yr, placed.birth_mass_solar, rng) if in_window else None
        pick -= chance
    return None


def _compact_pair_companion(placed, state, age_yr, rng):
    """The companion leaf of a compact (or white dwarf) placed object, drawn for `state`."""
    if state in ("star", "brown-dwarf"):
        if state == "brown-dwarf":
            return _star_leaf(rng.uniform(0.02, 0.07), age_yr, rng)
        low, high = (0.5, 16.0) if placed.kind == "neutron-star" else (0.1, 8.0)
        return _star_leaf(_alive_star_mass(low, high, age_yr, rng), age_yr, rng)
    limit = placed.mass_solar if state == "black-hole" else None
    return _remnant(state, age_yr, 0.0, rng, limit_mass=limit)


def _black_hole_tertiary(placed, age_yr, rng):
    if rng.random() < COMPACT_TERTIARY_SHARE:
        kind = "black-hole" if rng.random() < 0.8 else "neutron-star"
        return _remnant(kind, age_yr, 0.0, rng, limit_mass=placed.mass_solar)
    return _star_leaf(_alive_star_mass(0.1, 8.0, age_yr, rng), age_yr, rng)


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
                          separation_km=record["separation_au"] * constants.AU_M / constants.KM_TO_M_FACTOR,
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
                              initial_mass_solar=leaf.initial, radius_au=leaf.radius_au, seed=rng.getrandbits(63)))
    return MemberTree(tuple(members), tuple(nodes), age_yr, dropped)


def rejuvenated_age_yr(age_yr, rng=draw):
    """A merger product's effective age: a fraction (`REJUVENATION_FRACTION_RANGE`) of the system's age. The
    system's own age stays on the tree."""
    return age_yr * rng.uniform(*REJUVENATION_FRACTION_RANGE)


def mark_merger_product(member, age_yr, rng=draw):
    """The member as a merger product: flagged `merger_product`, with a rejuvenated `effective_age_yr`."""
    return dataclasses.replace(member, flags=member.flags | MEMBER_MERGER_PRODUCT,
                               effective_age_yr=rejuvenated_age_yr(age_yr, rng))


def _single(leaf, age_yr, dropped, rng):
    member = Member(0, None, 0, leaf.kind, leaf.state, leaf.flags, leaf.mass, leaf.initial, leaf.radius_au,
                    rng.getrandbits(63))
    return MemberTree((member,), (), age_yr, dropped)


def draw_member_tree(placed, age_yr, rng=draw, max_stars=12, n_stars=None, prevalence=1.0,
                     companion=_DRAW_COMPANION):
    """
    The members of the system around `placed` and their orbits, as a `MemberTree`.

    Args:
        placed (Placed): The scattered or filled object.
        age_yr (float): The system's age in years; states, remnants and the Peters and
            Kozai-Lidov checks read it.
        rng: The random source (a `draw.Stream` or the `draw` module).
        max_stars (int): The most members in the tree (`GEN.267`'s cap; 12).
        n_stars (int or None): Forces the number of stars of a star's tree (a directive).
            Ignored for a compact or white dwarf placed object, whose companions come from
            `MEMBER_STATE_BY_PLACED`.
        prevalence (float): Scales the chance of five or more stars (1 is the table).
        companion (str or None): For a compact or white dwarf placed object, the companion state the caller's
            screen already drew (`companion_state_at`): the tree takes it instead of drawing its own, `None`
            for no companion. Left out, the tree draws it (`companion_state`). Ignored for a star.

    Returns:
        MemberTree: The placed object is member 0. A single object has one member and no nodes.
    """
    leaves, units = _roster(placed, age_yr, rng, max_stars, n_stars, prevalence, companion)
    dropped = 0
    while True:
        if len(units) == 1 and isinstance(units[0], int):
            return _single(leaves[units[0]], age_yr, dropped, rng)
        for _ in range(TREE_RESTARTS):
            builder = _Builder(leaves, age_yr, rng)
            try:
                builder.place(_shape(units, rng), None, 0, None)
            except _Failed:
                continue
            if not builder.clears_extents():
                continue
            builder.finish()
            return _assemble(leaves, builder, age_yr, dropped, rng)
        loose = [u for u in units[1:] if isinstance(u, int)]
        if loose:
            drop = min(loose, key=lambda i: leaves[i].mass)
            units.remove(drop)
        elif len(units) > 1:
            units.pop()
        else:
            units[0] = units[0][0]
        dropped += 1
