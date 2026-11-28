# tests/audit/catalogue.py

"""
The catalogue type audit (TEST.164, docs/design/systems/multi-star-systems.md section 15, tests 10 to 14) on the census harness
(`tests.audit.census`, `tests.audit.multiplicity`).

A seeded sample of `SYSTEMS` systems of each region class of the representation audit (`tests.audit.representation.REGION_POINTS`: the
bulge, the thin disk in an arm and between arms, the thick disk, the halo, a satellite galaxy, a globular cluster and intergalactic space) is drawn
the way the scatters draw it: a star of the region's population mix (`galaxy.profile.profile_at(...).mix`, the densities of the populations
at the region's point), its mass, age and metallicity from `star_population.sample_mass_band_stars` for its population, and its member tree from
`scatter_companions.star_tree`. A second sample of white dwarfs (`SYSTEMS` of each region class, aged from the window of their population) and a
third of young massive stars (`GROUP_ROWS`, aged from 10 kyr to 10 Myr) give the reference classes a field of living stars holds too few of. A
fourth sample (`SYSTEMS` living stars of a disk of constant star formation, aged 0.05 to 10 Gy, the registry's own reference of the Algol row:
"disk stars", design note 13.2, `TUNING_KNOBS["tight_tail_weight"]`, GEN.270, GEN.1855) is the Algol's: the share of the region classes' own population mixes
depends on how old their stars are, which the registry's share does not.
Every row and tree is classified and tested:

- **Test 10, type counts**: each type of `TYPES` (the rows of `generation.phenomena.system_types.SYSTEM_TYPE_RATES` that are a property of a
  drawn tree) is counted in its reference class (the FGK dwarfs, the disk stars, the white dwarfs, the O stars, ...) and the share is within four
  Poisson sigma of the registry's default, or of its low to high range where it has one (the range widened by the sigma), over the region classes
  together and in each region class that holds expected counts of at least `MIN_EXPECTED`.
- **Test 11, region gates**: no type is found in a population its registry row does not list (`SystemType.populations`, with the satellite twins
  of the disk's); the pair gates of `multistar._gate_forbids` are checked on the stored pairs too (no contact binary of solar-type stars in the
  young population, no Algol in the old).
- **Test 12, the N tail**: the systems with N = 8 to 12 stars against the columns of the `eight-to-twelve-stars` row; the systems of 5, 6 and 7
  stars of the living stars against what the `star-count` table gives them (the registry's `quintuple-to-septuple` columns come from surveys of
  living primaries, so the living stars are the comparison, and the report sets the registry's columns beside them).
- **Test 13, interaction consistency**: every pair of two single members is classified again from its stored masses, radii and orbit
  (`multistar.classify_pair`): a `semi-detached` pair has a donor at or past its Roche lobe at closest approach, a `contact` pair both, an
  `accreting` pair a lobe-filling donor and a receiver that is a white dwarf, a neutron star or a black hole (or a star flagged `mass_gainer`).
- **Test 14, group decay**: every stored `group` node is younger than its members' `decay_age_yr`; a group older than it is not stored.

The types that are not a property of a drawn tree (the planets of a pair, a Galaxy count, a disk) are listed in the report with the stage that
draws them (`NOT_A_TREE_PROPERTY`). The sample is fixed by the region and the chunk number, so every run draws the same systems.
"""

import collections
import concurrent.futures
import dataclasses
import functools
import math
import multiprocessing
import os
import random

from planetgen import tuning
from planetgen.galaxy import profile
from planetgen.generation import multistar, scatter_companions, star_population
from planetgen.generation.phenomena import system_types
from planetgen.physics import remnant_mass, roche
from planetgen.physics import stellar_evolution as evolution
from planetgen.util import draw
from tests.audit import multiplicity
from tests.audit import representation

SYSTEMS = 1_000_000
"""int: Systems of each region class, of its stars and of its white dwarfs (TEST.164: at least 1e6 per region class)."""

GROUP_ROWS = 1_000_000
"""int: Young massive stars of the group sample."""

CHUNK = 20_000
"""int: The rows of one worker's job."""

MIN_EXPECTED = 30
"""float: The fewest expected counts a region class needs to have its own check of a type; below it only the regions together are checked."""

SIGMA = multiplicity.SIGMA_COUNT
"""float: Four Poisson sigma, the harness's rule."""

EXAMPLES = 5
"""int: How many examples of a broken rule a report shows."""

STAR, WHITE_DWARF, YOUNG_MASSIVE = "star", "white-dwarf", "young-massive"
"""str: Three of the four samples (the fourth is `CONSTANT_DISK`)."""

GROUP_AGE_YR = (1.0e4, 1.0e7)
"""tuple: The ages, years, the young massive stars are drawn over (log-uniform): from the first stars to a group's longest decay age."""

GROUP_MASS_RANGE_SOL = (9.0, 60.0)
"""tuple: The birth masses of the young massive stars: O and early B stars, the mass range of `young-group`."""

WHITE_DWARF_MIN_WINDOW_GY = 0.5
"""float: A white dwarf is drawn in the populations whose age window reaches this old (a younger population has made none)."""

REGIONS = tuple(representation.REGION_POINTS)
"""tuple: The region classes, by name."""

DISK_AGE_GY = (0.05, 10.0)
"""tuple: The ages, Gy, of the disk of constant star formation of the Algol row: the registry's share is "of a disk of constant star formation (0.05 to 10 Gy)" (design note
13.2, GEN.270's measure), and `CONSTANT_DISK` is drawn so."""

CONSTANT_DISK = "constant-disk"
"""str: The fourth sample: living stars of a disk of constant star formation (`DISK_AGE_GY`), which the Algol row's reference class "disk stars" is."""


@functools.lru_cache(maxsize=None)
def region_mix(region):
    """`([population], [weight])` of a region class: the densities of the populations at its point (the sector profile's `mix`)."""
    point = representation.REGION_POINTS[region][0]
    mix = profile.profile_at(point, representation.WITH_COMPONENTS, representation.EDGE, representation.E1).mix
    names = sorted(name for name, weight in mix.items() if weight > 0.0 and name in tuning.STELLAR_POPULATION_AGE_RANGES_GY)
    return names, [mix[name] for name in names]


def rows_of(sample, region, chunk, count):
    """The seeds of one chunk: a function of the sample, the region and the chunk number."""
    rng = random.Random(f"catalogue:{sample}:{region}:{chunk}")
    return [rng.getrandbits(63) for _ in range(count)]


def jobs(systems=None, group_rows=None):
    """`[(sample, region, chunk number, rows)]` for the whole sample."""
    systems = SYSTEMS if systems is None else systems
    group_rows = GROUP_ROWS if group_rows is None else group_rows
    found = []
    for sample, region, total in ([(STAR, region, systems) for region in REGIONS] + [(WHITE_DWARF, region, systems) for region in REGIONS]
                                  + [(YOUNG_MASSIVE, "young", group_rows), (CONSTANT_DISK, CONSTANT_DISK, systems)]):
        number = 0
        while total > 0:
            found.append((sample, region, number, min(CHUNK, total)))
            total -= CHUNK
            number += 1
    return found


@dataclasses.dataclass(frozen=True)
class Row:
    """What the audit knows of one drawn system: its sample, region class, population, the placed star's birth mass, age and Yerkes class."""
    sample: str
    region: str
    population: str
    birth_mass_solar: float
    age_gy: float
    yerkes: str
    fe_h: float
    seed: int


def draw_row(sample, region, seed):
    """`(Row, bright_stars row)` of the system with this seed."""
    rng = random.Random(seed ^ 0x5DEECE66D)
    if sample == YOUNG_MASSIVE:
        birth = math.exp(rng.uniform(math.log(GROUP_MASS_RANGE_SOL[0]), math.log(GROUP_MASS_RANGE_SOL[1])))
        age_gy = min(math.exp(rng.uniform(math.log(GROUP_AGE_YR[0]), math.log(GROUP_AGE_YR[1]))) / 1.0e9,
                     0.95 * evolution.main_sequence_lifetime_gy(birth))     # a star is drawn alive
        info = Row(sample, region, "young", birth, age_gy, "V", 0.0, seed)
        return info, multiplicity._star_row(seed, birth, age_gy, birth, "V")
    if sample == CONSTANT_DISK:
        age_gy = rng.uniform(*DISK_AGE_GY)
        birth = multistar._alive_star_mass(0.1, 60.0, age_gy * 1.0e9, rng)         # a star alive at its age, over 0.1 to 60 solar masses (GEN.270's draw)
        info = Row(sample, region, "disk", birth, age_gy, "V", 0.0, seed)
        return info, multiplicity._star_row(seed, birth, age_gy, birth, "V")
    names, weights = region_mix(region)
    if sample == WHITE_DWARF:
        names, weights = zip(*[(n, w) for n, w in zip(names, weights) if tuning.STELLAR_POPULATION_AGE_RANGES_GY[n][1] >= WHITE_DWARF_MIN_WINDOW_GY])
    population = rng.choices(names, weights)[0]
    if sample == WHITE_DWARF:
        low, high = tuning.STELLAR_POPULATION_AGE_RANGES_GY[population]
        for _ in range(256):
            age_gy = rng.uniform(low, max(high, low + 1.0e-6))
            birth = (0.85 ** -1.35 - rng.random() * (0.85 ** -1.35 - 8.0 ** -1.35)) ** (-1.0 / 1.35)    # a Salpeter slope over 0.85 to 8
            if evolution.main_sequence_lifetime_gy(birth) < age_gy:
                break
        else:
            birth, age_gy = 1.0, max(high, 10.0)
        info = Row(sample, region, population, birth, age_gy, "VII", 0.0, seed)
        return info, multiplicity._star_row(seed, birth, age_gy, remnant_mass.white_dwarf_mass_sol(birth), "VII")      # the product's initial-final relation: the merger share and the pair types read the mass
    star = star_population.sample_mass_band_stars(1, multistar.companion_floor_solar(), None, population, rng=draw.Stream(seed))[0]
    row = list(multiplicity._star_row(seed, star["initial_mass_sol"], star["age_gy"], star["mass_kg"] / multiplicity.constants.SOLAR_MASS_TO_KG,
                                      star["yerkes_class"]))
    row[multiplicity._BRIGHT["metallicity_fe_h"]] = star["metallicity_fe_h"]
    info = Row(sample, region, population, star["initial_mass_sol"], star["age_gy"], star["yerkes_class"], star["metallicity_fe_h"], seed)
    return info, tuple(row)


def tree_of(sample, region, seed):
    """`(Row, tree)` of the system with this seed: the tree the draw gives (`scatter_companions.star_tree` with `keep_single`), a single object's one member
    included, so that the flags of a single object (a Be star, a runaway, a merger product) are counted too."""
    info, row = draw_row(sample, region, seed)
    return info, scatter_companions.star_tree(row, keep_single=True)


# --- The reference classes and the types ----------------------------------------------------------------------------


OWN_SAMPLE = {"young-group": "young", "algol": CONSTANT_DISK}
"""dict: The types counted in a sample of their own, not region by region: the young group in the young massive stars, the Algol in the constant-star-formation disk."""


def _is_star(info):
    return info.sample == STAR


REFERENCES = {
    "FGK dwarf systems": lambda info, tree: _is_star(info) and info.yerkes == "V" and 0.45 <= info.birth_mass_solar < 1.4,
    "disk stars": lambda info, tree: info.sample == CONSTANT_DISK,
    "Sun-like stars": lambda info, tree: _is_star(info) and info.yerkes == "V" and 0.8 <= info.birth_mass_solar < 1.2,
    "young B stars": lambda info, tree: (_is_star(info) and 2.1 <= info.birth_mass_solar < 16.0
                                         and info.age_gy < tuning.YOUNG_POPULATION_MAX_AGE_GY),
    "O stars": lambda info, tree: _is_star(info) and info.birth_mass_solar >= 16.0,
    "white dwarfs": lambda info, tree: info.sample == WHITE_DWARF,
    "old white dwarfs": lambda info, tree: info.sample == WHITE_DWARF and info.age_gy >= tuning.OLD_POPULATION_MIN_AGE_GY,
    "white-dwarf companions": lambda info, tree: info.sample == WHITE_DWARF and bool(tree.nodes),
    "old and intermediate stars": lambda info, tree: (_is_star(info) and info.age_gy >= tuning.YOUNG_POPULATION_MAX_AGE_GY
                                                      and info.yerkes != "VII"),
    "massive main-sequence stars": lambda info, tree: _is_star(info) and info.yerkes == "V" and info.birth_mass_solar >= 8.0,
    "giants": lambda info, tree: _is_star(info) and info.yerkes in ("III", "II", "I") and tree.members[0].kind == "giant",
    "living massive stars with a compact companion": lambda info, tree: (
        _is_star(info) and info.birth_mass_solar >= 14.0 and tree.members[0].kind in ("main-sequence", "giant") and tree.compact_count >= 1),
    "massive systems younger than 0.3 Myr": lambda info, tree: info.sample == YOUNG_MASSIVE and info.age_gy * 1.0e9 < tuning.YOUNG_GROUP_AGE_MAX_YR,
    "systems": lambda info, tree: info.sample == STAR,
}
"""dict: Reference class -> the predicate over a `Row` and its tree that puts a system in it (the registry's `reference` text, which the type
rows name; a class that reads the tree is one the registry defines by the system's contents)."""

REFERENCE_OF = {"FGK dwarf systems": "FGK dwarf systems", "disk stars": "disk stars", "Sun-like stars": "Sun-like stars", "B stars": "young B stars",
                "O stars": "O stars", "white dwarfs": "white dwarfs", "systems": "systems",
                "massive systems younger than 0.3 Myr": "massive systems younger than 0.3 Myr",
                "old and intermediate stars": "old and intermediate stars", "massive main-sequence stars": "massive main-sequence stars",
                "white-dwarf companions": "white-dwarf companions", "giants": "giants",
                "compact companions of living stars of 14 solar masses or more": "living massive stars with a compact companion"}
"""dict: The registry's `reference` text -> the key of `REFERENCES` the audit counts it in (the Be star's "B stars" are the young ones that host it)."""


def _pairs(tree):
    """`[(node, (member, member))]` of the nodes of the tree whose two children are single members."""
    found = []
    for node in tree.nodes:
        pair = [m for m in tree.members if m.parent_node_no == node.node_no]
        if len(pair) == 2:
            found.append((node, (min(pair, key=lambda m: m.slot), max(pair, key=lambda m: m.slot))))
    return found


def _any_flag(tree, flag):
    return any(m.flags & flag for m in tree.members)


def _accreting_with(tree, kinds_a, kinds_b):
    return any(node.interaction == "accreting" and ((a.kind in kinds_a and b.kind in kinds_b) or (b.kind in kinds_a and a.kind in kinds_b))
               for node, (a, b) in _pairs(tree))


def _brown_dwarf_axis(tree):
    """The axes, AU, of the nodes a brown dwarf member sits on directly."""
    by_no = {node.node_no: node for node in tree.nodes}
    return [by_no[m.parent_node_no].separation_au for m in tree.members if m.kind == "brown-dwarf" and m.parent_node_no in by_no]


TYPES = {
    "contact-binary": ("FGK dwarf systems", lambda tree: any(n.interaction == "contact" for n in tree.nodes)),
    "algol": ("disk stars", lambda tree: any(node.interaction == "semi-detached" and any(m.flags & multistar.MEMBER_MASS_REVERSED for m in pair)
                                             for node, pair in _pairs(tree))),
    "be-star": ("young B stars", lambda tree: _any_flag(tree, multistar.MEMBER_BE)),
    "runaway-star": ("O stars", lambda tree: _any_flag(tree, multistar.MEMBER_RUNAWAY_SUPERNOVA | multistar.MEMBER_RUNAWAY_DYNAMICAL)),
    "brown-dwarf-close": ("Sun-like stars", lambda tree: any(axis < 10.0 for axis in _brown_dwarf_axis(tree))),
    "brown-dwarf-wide": ("Sun-like stars", lambda tree: any(28.0 <= axis <= 1590.0 for axis in _brown_dwarf_axis(tree))),
    "white-dwarf-companion": ("white dwarfs", lambda tree: bool(tree.nodes)),
    "close-double-white-dwarf": ("white dwarfs", lambda tree: any(m.kind == "white-dwarf" for m in tree.members[1:])),
    "double-white-dwarf-merger": ("old white dwarfs", lambda tree: bool(tree.members[0].flags & multistar.MEMBER_MERGER_PRODUCT)),
    "cataclysmic-variable": ("white dwarfs", lambda tree: _accreting_with(tree, ("white-dwarf",), ("main-sequence",))),
    "am-cvn": ("white dwarfs", lambda tree: _accreting_with(tree, ("white-dwarf",), ("stripped-helium",))),
    "blue-straggler-field": ("old and intermediate stars", lambda tree: bool(tree.members[0].flags & multistar.MEMBER_MERGER_PRODUCT)),
    "massive-merger-product": ("massive main-sequence stars", lambda tree: bool(tree.members[0].flags & multistar.MEMBER_MERGER_PRODUCT)),
    "white-dwarf-post-ce": ("white-dwarf companions", lambda tree: any(
        node.interaction == "post-common-envelope" and any(m.kind == "white-dwarf" for m in pair) for node, pair in _pairs(tree))),
    "symbiotic-barium": ("giants", lambda tree: bool(tree.members[0].flags & multistar.MEMBER_SYMBIOTIC)),
    "high-mass-xray-binary": ("living massive stars with a compact companion", lambda tree: any(
        node.interaction == "accreting" and any(m.kind in multistar.COMPACT_KINDS for m in pair) for node, pair in _pairs(tree))),
    "young-group": ("massive systems younger than 0.3 Myr", lambda tree: any(n.kind == "group" for n in tree.nodes) or bool(tree.ejected)),
}
"""dict: Registry type -> `(reference class key of REFERENCES, predicate over a tree)`: the types that are a property of a drawn tree.
The predicate is asked of the systems that have a tree; a system with none (a single star) has none of these types."""

NOT_A_TREE_PROPERTY = {
    "close-binary": "a period range of a pair of the architecture (the closest pair's axis), counted by the multiplicity audit's tables",
    "intermediate-binary": "a period range of the architecture", "wide-binary": "a period range of the architecture",
    "ultra-wide-pair": "a period range of the architecture",
    "triple-2plus1": "the shape of the hierarchy, drawn at 60 and 40 percent (`multistar` shape draw, TEST.162)",
    "quadruple-3plus1": "the shape of the hierarchy", "quadruple-2plus2": "the shape of the hierarchy",
    "multiple-with-disk": "the disk descriptor drawn at fill (`draw_disk_descriptor`)",
    "luminous-blue-variable": "a count of the Galaxy, not a rate", "rcb-fk-comae": "a count of the Galaxy, not a rate",
    "thorne-zytkow": "a count of the Galaxy, not a rate", "microquasar": "a count of the Galaxy, not a rate",
    "hypervelocity-star": "a count of the Galaxy, drawn by the scatter of phenomena",
    "supermassive-black-hole-binary": "an extragalactic nucleus flag", "magnetar-binary": "a theoretical rate of 0",
    "rogue-planet-pair": "a rogue planet row", "s-type-planet": "the planet stage at fill", "circumbinary-planet": "the planet stage at fill",
    "circumbinary-disk": "the disk descriptor drawn at fill", "white-dwarf-pollution": "the planet stage at fill",
    "planetary-mass-companion": "a theoretical rate of the substellar draw at fill",
}
"""dict: Registry types the audit does not count in a tree, with why (the report lists them so that none is silently dropped)."""


COUNTED_ELSEWHERE = frozenset({
    "black-hole-black-hole", "dormant-black-hole-star", "double-neutron-star", "neutron-star-black-hole", "intermediate-mass-black-hole-star",
    "pulsar-massive-companion", "millisecond-pulsar-companion", "dormant-neutron-star-star",
    "quintuple-to-septuple", "eight-to-twelve-stars"})
"""frozenset: The registry types another check of the audit covers: the compact companions are the companion-state tables of the multiplicity audit
(`tests.audit.multiplicity`, TEST.162), the N columns are `tail_checks` here."""

NOT_AUDITED_HERE = {
    "brown-dwarf-pair": "its reference class, brown dwarfs, is not a sample of this audit (the sample's primaries are stars and white dwarfs)",
    "cluster-compact-binary": "its reference class, the black holes of clusters, is a phenomenon scatter class (a follow-up of TEST.162's sample)",
    "compact-triple": "theoretical, with no sample of compact binaries (a follow-up of TEST.162's sample)",
    "decayed-trapezium-remnant": "a decayed group is a tree with an ejected star; the audit counts those (`ejected` of the young massive sample) and that no group is older than its decay age",
    "hot-subdwarf-stripped": "no default of its own (a share of the post-common-envelope pairs, `tuning.HOT_SUBDWARF_POST_CE_SHARE`)",
    "low-mass-xray-binary": "its reference class, neutron stars with a low-mass star, is a phenomenon scatter class (a follow-up of TEST.162's sample)",
    "unbound-by-collapse": "the origin of a runaway, drawn with it (`draw_runaway_origin`); the runaway rate is audited, its two channels are not",
    "white-dwarf-brown-dwarf": "no default of its own (part of the brown dwarf companion rates)",
    "wolf-rayet-binary": "a Wolf-Rayet star is a late phase of a massive star and the sample's living stars are not drawn in it",
}
"""dict: Registry types this audit does not count, with why: each needs a sample the audit does not draw or has no rate of its own. The report lists them."""


REPORT_ONLY = {
    "brown-dwarf-close": "its generator (`draw_substellar_companion`, GEN.275) is not built: the brown dwarfs of a tree come from the mass ratio law and the brown dwarf draw",
    "brown-dwarf-wide": "its generator (`draw_substellar_companion`, GEN.275) is not built: the brown dwarfs of a tree come from the mass ratio law and the brown dwarf draw",
}
"""dict: Types the audit counts and reports but does not check: the registry row names a generator that has not landed, so the rate is not yet the row's."""

AT_LEAST = {
    "white-dwarf-post-ce": "the draw gives the share (`post_ce_fraction`) and the history rule (`pair_history`, GEN.1012) makes every pair whose giant reached its lobe "
                           "one too, so the registry's 10% is a floor",
}
"""dict: Types whose registry rate is the share the draw gives, to which the geometry or the history adds pairs: the audit checks that the count is at least the share."""


def registry_rate(name):
    """`(default, low, high)` of a registry type (low and high are the default when it has no range)."""
    entry = system_types.SYSTEM_TYPE_RATES[name]
    low = entry.default if entry.low is None else entry.low
    high = entry.default if entry.high is None else entry.high
    return entry.default, low, high


def populations_of(name):
    """The populations the registry row allows (its `populations`, which carry the satellite twins)."""
    return set(system_types.SYSTEM_TYPE_RATES[name].populations)


# --- The checks of one tree ---------------------------------------------------------------------------------------------


def interaction_problems(tree):
    """`[str]`: the pairs of two single members whose stored `interaction` the geometry does not give (test 13)."""
    found = []
    for node, (a, b) in _pairs(tree):
        if node.interaction not in ("contact", "semi-detached", "accreting"):
            continue
        peri = node.separation_au * (1.0 - node.eccentricity)
        fill_a = roche.filling_factor(a.radius_au, peri, a.mass_solar, b.mass_solar)
        fill_b = roche.filling_factor(b.radius_au, peri, b.mass_solar, a.mass_solar)
        if node.interaction == "contact" and not (fill_a >= 1.0 and fill_b >= 1.0):
            found.append(f"contact pair at node {node.node_no} has a member inside its lobe ({fill_a:.3g}, {fill_b:.3g})")
        if node.interaction == "semi-detached" and not (fill_a >= 1.0 or fill_b >= 1.0):
            found.append(f"semi-detached pair at node {node.node_no} has no donor at its lobe ({fill_a:.3g}, {fill_b:.3g})")
        if node.interaction == "accreting":
            receivers = [(m, fill) for m, other, fill in ((a, b, fill_b), (b, a, fill_a)) if m.kind in multistar.ACCRETOR_KINDS or m.flags & multistar.MEMBER_MASS_GAINER]
            if not receivers:
                found.append(f"accreting pair at node {node.node_no} has no receiver that is degenerate or a mass gainer")
            elif not any(fill >= 1.0 for _m, fill in receivers) and not node.flags & (multistar.NODE_WIND_FED | multistar.NODE_DISK_FED):
                found.append(f"accreting pair at node {node.node_no} has no donor at its lobe ({fill_a:.3g}, {fill_b:.3g})")
    return found


def gate_problems(tree, info):
    """
    `[str]`: the types of the tree that its age's gate forbids (test 11). The generator knows a system's age, not its population, so a type
    is allowed at the ages of the windows of the populations its registry row lists (`tuning.STELLAR_POPULATION_AGE_RANGES_GY`; a young
    group's massive stars are the young window's): a type is forbidden when its reference class holds the system, the tree has the type
    and the age is in none of the windows. The pair gates of `multistar._gate_forbids` are checked on the stored pairs as well.
    """
    found = []
    for name, (ref_name, predicate) in TYPES.items():
        if not REFERENCES[ref_name](info, tree):
            continue
        windows = [tuning.STELLAR_POPULATION_AGE_RANGES_GY[p] for p in populations_of(name) if p in tuning.STELLAR_POPULATION_AGE_RANGES_GY]
        if any(low <= info.age_gy <= high for low, high in windows):
            continue
        if predicate(tree):
            found.append(f"{name} at an age of no population it lists")
    for node, (a, b) in _pairs(tree):
        if node.interaction == "contact" and info.age_gy < tuning.YOUNG_POPULATION_MAX_AGE_GY and max(a.initial_mass_solar, b.initial_mass_solar) < multistar.MASSIVE_PRIMARY_MIN_MASS_SOL:
            found.append("a contact binary of solar-type stars in the young population")
        if node.interaction == "semi-detached" and info.age_gy >= tuning.OLD_POPULATION_MIN_AGE_GY:
            found.append("an Algol in the old population")
    return found


def group_problems(tree):
    """`[str]`: the group nodes of the tree older than their decay age (test 14)."""
    found = []
    for node in tree.nodes:
        if node.kind != "group":
            continue
        decays = [m.decay_age_yr for m in tree.members if m.decay_age_yr is not None]
        if not decays or any(tree.age_yr >= decay for decay in decays):
            found.append(f"group node {node.node_no} at {tree.age_yr:.3g} yr, decay ages {sorted(set(round(d) for d in decays))[:3]}")
    return found


# --- One chunk ----------------------------------------------------------------------------------------------------------


@dataclasses.dataclass
class Tally:
    """What the rows of one or more chunks add up to, kept apart by `(sample, region)`."""
    rows: collections.Counter = dataclasses.field(default_factory=collections.Counter)
    references: collections.Counter = dataclasses.field(default_factory=collections.Counter)
    expected_stars: collections.Counter = dataclasses.field(default_factory=collections.Counter)
    alive_stars: collections.Counter = dataclasses.field(default_factory=collections.Counter)
    types: collections.Counter = dataclasses.field(default_factory=collections.Counter)
    stars: collections.Counter = dataclasses.field(default_factory=collections.Counter)
    pairs_checked: collections.Counter = dataclasses.field(default_factory=collections.Counter)
    interaction_bad: collections.Counter = dataclasses.field(default_factory=collections.Counter)
    gate_bad: collections.Counter = dataclasses.field(default_factory=collections.Counter)
    gate_kinds: collections.Counter = dataclasses.field(default_factory=collections.Counter)
    groups: int = 0
    groups_bad: int = 0
    decayed_rows: int = 0
    ejected: int = 0
    examples: list = dataclasses.field(default_factory=list)

    def add(self, other):
        for field in dataclasses.fields(self):
            mine, theirs = getattr(self, field.name), getattr(other, field.name)
            if isinstance(mine, collections.Counter):
                mine.update(theirs)
            elif isinstance(mine, list):
                mine[:] = (mine + theirs)[:EXAMPLES * 4]
            else:
                setattr(self, field.name, mine + theirs)
        return self


def run_chunk(job):
    """The `Tally` of one job `(sample, region, chunk number, rows)`."""
    sample, region, chunk, count = job
    tally = Tally()
    for seed in rows_of(sample, region, chunk, count):
        info, tree = tree_of(sample, region, seed)
        tally.rows[(sample, region)] += 1
        for ref_name, predicate in REFERENCES.items():
            if predicate(info, tree):
                tally.references[(region, ref_name)] += 1
        if sample == STAR:
            tally.stars[(region, tree.star_count)] += 1
        if sample == YOUNG_MASSIVE and info.age_gy * 1.0e9 >= tuning.YOUNG_GROUP_AGE_MAX_YR:
            tally.decayed_rows += 1
        if sample == YOUNG_MASSIVE and tree.ejected:
            tally.ejected += 1
        for name, (ref_name, predicate) in TYPES.items():
            if REFERENCES[ref_name](info, tree) and predicate(tree):
                tally.types[(region, name)] += 1
        if sample == STAR and info.yerkes != "VII":
            tally.alive_stars[(region, "rows")] += 1
            tally.alive_stars[(region, tree.star_count)] += 1
            for count, share in star_count_shares(info.birth_mass_solar, info.fe_h).items():
                tally.expected_stars[(region, count)] += share
        pairs = [node for node, _pair in _pairs(tree) if node.interaction in ("contact", "semi-detached", "accreting")]
        tally.pairs_checked[region] += len(pairs)
        for problem in interaction_problems(tree):
            tally.interaction_bad[region] += 1
            if len(tally.examples) < EXAMPLES * 4:
                tally.examples.append(f"{sample}/{region} seed {seed}: {problem}")
        for problem in gate_problems(tree, info):
            tally.gate_bad[region] += 1
            tally.gate_kinds[problem] += 1
            if len(tally.examples) < EXAMPLES * 4:
                tally.examples.append(f"{sample}/{region} seed {seed}: {problem}")
        groups = sum(1 for n in tree.nodes if n.kind == "group")
        tally.groups += groups
        for problem in group_problems(tree):
            tally.groups_bad += 1
            if len(tally.examples) < EXAMPLES * 4:
                tally.examples.append(f"{sample}/{region} seed {seed}: {problem}")
    return tally


def run(workers=None, sample_jobs=None):
    """The `Tally` of the whole sample, the chunks drawn by forked workers (up to 4 by default, 1 to run them one after another)."""
    sample_jobs = jobs() if sample_jobs is None else sample_jobs
    workers = min(4, os.cpu_count() or 1) if workers is None else workers
    if workers <= 1:
        done = [run_chunk(job) for job in sample_jobs]
    else:
        with concurrent.futures.ProcessPoolExecutor(max_workers=workers, mp_context=multiprocessing.get_context("fork")) as pool:
            done = list(pool.map(run_chunk, sample_jobs, chunksize=2))
    total = Tally()
    for tally in done:
        total.add(tally)
    return total


# --- The expectations ---------------------------------------------------------------------------------------------------


def rate_check(name, observed, reference_count, default, low, high, label=""):
    """A `multiplicity.Check` of a type's count in its reference class: within `SIGMA` Poisson sigma of the default, or, where the row has a
    range, of the nearest end of the range (the range widened by the sigma)."""
    text = f"{label}{name}: {observed:,} of {reference_count:,}"
    if name in AT_LEAST:
        expected = reference_count * default
        return multiplicity.Check(text + " (a floor)", observed, expected, math.sqrt(max(expected, 1.0)), SIGMA, "at least")
    if low == high:
        return multiplicity.share_check(text, observed, reference_count, default)
    nearest = min(max(observed, reference_count * low), reference_count * high)
    return multiplicity.Check(text + f" (range {low:g} to {high:g})", observed, nearest, math.sqrt(max(nearest, 1.0)), SIGMA)


def type_checks(tally):
    """Test 10: every type of `TYPES` in its reference class, over the regions together and in each region with enough expected counts."""
    found = []
    for name, (ref_name, _predicate) in TYPES.items():
        if name in REPORT_ONLY:
            continue
        default, low, high = registry_rate(name)
        reference_total = sum(tally.references[(region, ref_name)] for region in REGIONS)
        if name in OWN_SAMPLE:
            where = OWN_SAMPLE[name]
            found.append(rate_check(name, tally.types[(where, name)], tally.references[(where, ref_name)], default, low, high))
            continue
        observed = sum(tally.types[(region, name)] for region in REGIONS)
        found.append(rate_check(name, observed, reference_total, default, low, high, "all regions: "))
        for region in REGIONS:
            count = tally.references[(region, ref_name)]
            if count * (low if low == high else max(default, low)) >= MIN_EXPECTED:
                found.append(rate_check(name, tally.types[(region, name)], count, default, low, high, f"{region}: "))
    return found


def gate_checks(tally):
    """Test 11: no type in a population its row does not allow, no gated pair."""
    return [multiplicity.Check(f"{region}: stored types or pairs a gate forbids", tally.gate_bad[region], 0, 1.0, 0.5) for region in REGIONS + ("young", CONSTANT_DISK)]


def star_count_shares(initial_mass_solar, fe_h):
    """`{N: chance}` of a system of N stars (1 to `multistar.MAX_MEMBERS`) of a star of this birth mass and [Fe/H]: the `star-count` table as the tree draws
    it (`multistar.close_trend` moves it for a star of another [Fe/H] than the Sun's), and the tail from 7 stars up (`multistar.TAIL_RATIOS`)."""
    weights, _inside = multistar.close_trend(initial_mass_solar, fe_h, 1.0)
    weights = weights or multistar.star_count_weights(initial_mass_solar, 1.0)
    total = sum(weights)
    shares = {count + 1: weight / total for count, weight in enumerate(weights)}
    top, reach, ratio = shares[7], 1.0, multistar.TAIL_RATIOS[0]
    shares[7] = 0.0
    for count in range(7, multistar.MAX_MEMBERS):
        step = min(1.0, ratio)
        shares[count] += top * reach * (1.0 - step)
        reach *= step
        ratio = multistar.TAIL_RATIOS[1]
        shares[count + 1] = 0.0
    shares[multistar.MAX_MEMBERS] += top * reach
    return shares


def tail_checks(tally):
    """
    Test 12, the N tail. The systems with 8 to 12 stars against the tail table's columns of the registry, over the regions together (the design
    note's test). The systems with 5, 6 and 7 stars against what the `star-count` table gives the living stars of the sample (`star_count_shares`,
    the tree's own draw, [Fe/H] trend included). The registry's columns for them come from surveys of living primaries, so living stars are the
    proper comparison; a row whose star has died (a white dwarf, neutron star or hole) takes its companions from the tables of its remnant, not the
    star count's, and is left out of it. `tail_report` sets the registry's columns against the living stars as well.
    """
    found = []
    systems = sum(tally.rows[(STAR, region)] for region in REGIONS)
    counts = collections.Counter()
    for (region, n), count in tally.stars.items():
        counts[n] += count
    tail = system_types.SYSTEM_TYPE_RATES["eight-to-twelve-stars"].by_class
    observed = sum(c for n, c in counts.items() if 8 <= n <= 12)
    share = sum(tail[f"N{n}"] for n in range(8, 13))
    found.append(multiplicity.share_check(f"systems with N = 8 to 12 stars (of {systems:,}, over all regions)", observed, systems, share))
    alive = sum(tally.alive_stars[(region, "rows")] for region in REGIONS)
    for column in (5, 6, 7):
        observed = sum(tally.alive_stars[(region, column)] for region in REGIONS)
        expected = sum(tally.expected_stars[(region, column)] for region in REGIONS)
        found.append(multiplicity.Check(f"systems with N = {column} stars (of {alive:,} living stars, over all regions, the star-count table's)",
                                        observed, expected, math.sqrt(max(expected, 1.0)), SIGMA))
    return found


def tail_report(tally):
    """The lines of the report that set the registry's N = 5, 6 and 7 columns against the living stars of the sample (not checked: `tail_checks` holds them to the table)."""
    alive = sum(tally.alive_stars[(region, "rows")] for region in REGIONS)
    five = system_types.SYSTEM_TYPE_RATES["quintuple-to-septuple"].by_class
    lines = []
    for column in (5, 6, 7):
        count = sum(tally.alive_stars[(region, column)] for region in REGIONS)
        share = count / alive
        lines.append(f"  N = {column}: {count:,} of {alive:,} living stars ({share:.3g}) against the registry's {five[f'N{column}']:.3g} "
                     f"(a ratio of {share / five[f'N{column}']:.2f})")
    return lines


def consistency_checks(tally):
    """Test 13: no pair whose stored interaction the geometry does not give."""
    return [multiplicity.Check(f"{region}: pairs failing the Roche check (of {tally.pairs_checked[region]:,} interacting)", tally.interaction_bad[region], 0, 1.0, 0.5)
            for region in REGIONS + ("young", CONSTANT_DISK)]


def decay_checks(tally):
    """Test 14: no stored group older than its decay age, and the sample holds groups and decayed rows (the test saw both sides)."""
    return [multiplicity.Check("stored groups older than their decay age", tally.groups_bad, 0, 1.0, 0.5),
            multiplicity.Check(f"young massive rows with a stored group (of {tally.groups:,} groups)", min(tally.groups, 1), 1, 1.0, 0.5),
            multiplicity.Check(f"young massive rows older than 0.3 Myr (of {tally.decayed_rows:,})", min(tally.decayed_rows, 1), 1, 1.0, 0.5)]


def all_checks(tally):
    """Every check of the audit over a sample."""
    return type_checks(tally) + gate_checks(tally) + tail_checks(tally) + consistency_checks(tally) + decay_checks(tally)


def format_report(tally, checks=None):
    """The audit's report: the sample, every check, the types that are not a tree's, and the examples of what broke."""
    checks = all_checks(tally) if checks is None else checks
    lines = ["Catalogue type audit (TEST.164)", "", "Sample (systems):"]
    for (sample, region), rows in sorted(tally.rows.items()):
        lines.append(f"  {sample:<14} {region:<22} {rows:>10,}")
    lines += ["", "Checks:"] + ["  " + check.line() for check in checks]
    bad = [check for check in checks if not check.ok]
    lines += ["", f"{len(checks)} checks, {len(bad)} failed.", "", "Types of the registry that are not a property of a drawn tree:"]
    lines += [f"  {name}: {why}" for name, why in NOT_A_TREE_PROPERTY.items()]
    lines += ["", "The registry's N = 5, 6 and 7 columns against the living stars of the sample:"]
    lines += tail_report(tally)
    lines += ["", "Counted and reported, not checked:"]
    for name, why in REPORT_ONLY.items():
        ref_name = TYPES[name][0]
        count = sum(tally.references[(region, ref_name)] for region in REGIONS)
        lines.append(f"  {name}: {sum(tally.types[(region, name)] for region in REGIONS):,} of {count:,} {ref_name} against the row's "
                     f"{system_types.SYSTEM_TYPE_RATES[name].default:g}; {why}")
    lines += ["", "Types of the registry the audit does not count here (they need a sample it does not draw, or have no rate of their own):"]
    lines += [f"  {name}: {why}" for name, why in NOT_AUDITED_HERE.items()]
    lines += ["", "Types the audit counts elsewhere: " + ", ".join(sorted(COUNTED_ELSEWHERE))]
    unaccounted = set(system_types.SYSTEM_TYPE_RATES) - set(TYPES) - set(NOT_A_TREE_PROPERTY) - set(NOT_AUDITED_HERE) - COUNTED_ELSEWHERE
    lines += ["Types of the registry that are in none of these lists: " + (", ".join(sorted(unaccounted)) or "none")]
    lines += ["", "Types a gate forbids, by kind:"] + [f"  {count:>8,}  {kind}" for kind, count in tally.gate_kinds.most_common(20)]
    lines += ["  example: " + example for example in tally.examples[:EXAMPLES * 2]]
    return "\n".join(lines) + "\n"


def report(workers=None):
    """`(text, failed checks)` for the whole sample."""
    tally = run(workers)
    checks = all_checks(tally)
    return format_report(tally, checks), [check for check in checks if not check.ok]


def main(argv=None):
    """`python -m tests.audit.catalogue [report.txt]`: runs the audit and prints (and writes) the report; exit code 1 when a check failed."""
    import sys
    argv = list(sys.argv[1:] if argv is None else argv)
    text, failed = report()
    print(text, end="")
    if argv:
        with open(argv[0], "w", encoding="utf-8") as handle:
            handle.write(text)
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
