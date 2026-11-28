# tests/audit/multiplicity.py

"""
The multiplicity audit (TEST.162, docs/design/systems/multi-star-systems.md section 9) on the census harness
(`tests.audit.census`): a seeded scatter of at least `ROWS` rows of each placed class gets its member trees drawn
(`generation.scatter_companions`, the function the scatter passes call) and every tree and every count is tested.

The harness counts sectors; a scatter row is not a sector. So this audit uses the harness's rule (`census.within_sigma`, four
Poisson sigma) and its report, with its own seeded sample: rows whose seeds are fixed by the class and the chunk number
(`rows_of`), so every run draws the same scatter, and chunks drawn at the same time by forked workers add up to the one-by-one
run.

What a chunk checks, for each row it holds:

- the companion state of the row (none, star, black hole, neutron star, white dwarf, a pulsar's spider) against its table
  (`multistar.MEMBER_STATE_BY_PLACED`), a star row's number of stars against the `star-count` table (`multistar.expected_star_count`);
- the bookkeeping of a row: the same-class companions per row against `multistar.same_class_companions`, so that rows and companions
  of a class keep its object density (`phenomenon_scatter.row_density_divisor`); nobody is both a row and a member (a member's seed is
  never the row's, no two members of a tree share a seed);
- every stored tree: each level passes the stability screen of `physics/hierarchy_stability.py`, no pair's Kozai cycle brings its
  pericenter inside four times its extents, and no compact pair has a Peters time shorter than the system's age;
- the same seed draws the same members (a sample of the trees is drawn twice);
- nothing a planet, a moon or a belt: every member's kind is one of `multistar.MEMBER_KINDS`, which is also the stored kinds.

The counts of black hole pairs, the nightly update and the fill that builds the stored members are database checks:
`tests.audit.multiplicity_db`.
"""

import collections
import concurrent.futures
import dataclasses
import math
import multiprocessing
import os
import random

from planetgen.db import orbit_trees, store
from planetgen.generation import multistar, phenomenon_scatter, scatter_companions
from planetgen.physics import constants, hierarchy_stability, stability
from planetgen.physics import stellar_evolution as evolution
from planetgen.physics.system_tree import TreeMember
from planetgen.util import draw
from tests.audit import census as audit

ROWS = 1_000_000
"""int: How many rows of each placed class of the phenomenon scatter, and in all of the star classes, the audit draws (TEST.162: at least 1e6)."""

CHUNK = 20_000
"""int: The rows of one worker's job."""

DETERMINISM_PER_CHUNK = 20
"""int: How many of a chunk's trees are drawn a second time and compared."""

EXAMPLES = 5
"""int: How many examples of a broken rule a report shows."""

SIGMA_COUNT = audit.DEFAULT_SIGMA
SIGMA_MEAN = audit.DEFAULT_SIGMA
"""float: The sigma of a share and of a mean number of stars: four, the harness's rule (the design note's two sigma for the means would fail
about one run in three over the nine star classes at a fixed seed)."""

_BRIGHT = {name: i for i, name in enumerate(store.BRIGHT_STAR_COLUMNS)}
_PHENOMENON = {name: i for i, name in enumerate(store.PHENOMENON_SCATTER_COLUMNS)}


@dataclasses.dataclass(frozen=True)
class PhenomenonClass:
    """One class of the phenomenon scatter whose rows have a companion screen: its name, the `(kind, subtype)` of the row, the placed object's
    mass (solar) and its `phenomenon_scatter.PLACED_KINDS` name."""
    name: str
    kind: str
    subtype: str | None
    mass_solar: float | None


PHENOMENON_CLASSES = (
    PhenomenonClass("black-hole", "black-hole", "stellar", 15.0),
    PhenomenonClass("intermediate-mass-black-hole", "black-hole", "intermediate", 1.0e3),
    PhenomenonClass("neutron-star", "neutron-star", None, 1.4),
    PhenomenonClass("planetary-nebula-star", "planetary-nebula", None, None),
)

STAR_CLASSES = (("O", 16.0, 60.0), ("B-early", 9.0, 16.0), ("B-mid", 5.0, 9.0), ("B", 2.1, 5.0), ("A", 1.4, 2.1), ("F", 1.04, 1.4),
                ("G", 0.8, 1.04), ("K", 0.45, 0.8), ("M", 0.1, 0.45))
"""tuple: The `star-count` classes with the birth masses (solar) a star row of the class is drawn over, log-uniform."""

WHITE_DWARF = "white-dwarf"
"""str: The name of the white dwarf rows of the star scatter (class VII)."""

STAR_AGE_CAP_GY = 10.0
"""float: The oldest a star row is drawn (the disk's age)."""

OLD_DENSITY_AGE_FRACTION = 0.5
"""float: A star row is this far through its main-sequence life, as the rows of the tests of the scatter companions are."""


def classes():
    """Every class the audit draws, by name: the phenomenon classes, the star classes by their spectral class and the white dwarf rows."""
    return [c.name for c in PHENOMENON_CLASSES] + [name for name, _lo, _hi in STAR_CLASSES] + [WHITE_DWARF]


def rows_in(name):
    """How many rows of the class the audit draws: `ROWS` of each phenomenon class and of the white dwarf rows, and `ROWS` over all the star
    classes (the same number of each, rounded up)."""
    return -(-ROWS // len(STAR_CLASSES)) if name in {s[0] for s in STAR_CLASSES} else ROWS


def jobs():
    """`[(class name, chunk number, rows in the chunk)]` for the whole sample."""
    found = []
    for name in classes():
        total, number = rows_in(name), 0
        while total > 0:
            found.append((name, number, min(CHUNK, total)))
            total -= CHUNK
            number += 1
    return found


def rows_of(name, chunk, count):
    """The seeds of one chunk: a function of the class's name and the chunk number, so the scatter is the same on every run."""
    rng = random.Random(f"multiplicity:{name}:{chunk}")
    return [rng.getrandbits(63) for _ in range(count)]


def _phenomenon_row(cls, seed):
    row = [None] * len(store.PHENOMENON_SCATTER_COLUMNS)
    for key, value in (("shell_index", 3), ("band_index", 0), ("wedge_index", 0), ("kind", cls.kind), ("subtype", cls.subtype),
                       ("position_x_mpc", 1000), ("position_y_mpc", 2000), ("position_z_mpc", 3000), ("seed", seed),
                       ("mass_solar", cls.mass_solar)):
        row[_PHENOMENON[key]] = value
    return tuple(row)


def _star_row(seed, initial_mass, age_gy, mass, yerkes):
    row = [None] * len(store.BRIGHT_STAR_COLUMNS)
    for key, value in (("shell_index", 3), ("band_index", 0), ("wedge_index", 0), ("position_x_mpc", 1000), ("position_y_mpc", 2000),
                       ("position_z_mpc", 3000), ("population", "thin-disk"), ("star_type", "audit"), ("yerkes_class", yerkes),
                       ("mass_kg", mass * constants.SOLAR_MASS_TO_KG), ("radius_km", 1.0e7), ("temperature_k", 1.0e4),
                       ("luminosity_w", 1.0e27), ("age_gy", age_gy), ("lifespan_gy", age_gy * 2.0), ("initial_mass_sol", initial_mass),
                       ("phase_end_age_gy", age_gy * 2.0), ("seed", seed)):
        row[_BRIGHT[key]] = value
    return tuple(row)


def star_row(name, seed):
    """The `bright_stars` row of one star of class `name` (a spectral class or `WHITE_DWARF`): its birth mass log-uniform over the class's
    range and its age half way through its main-sequence life (at most `STAR_AGE_CAP_GY`), both drawn from the seed; a white dwarf of
    0.6 solar masses from a 1 to 6 solar mass star, older than the star's life."""
    rng = random.Random(seed ^ 0x5DEECE66D)
    if name == WHITE_DWARF:
        birth = math.exp(rng.uniform(math.log(1.0), math.log(6.0)))
        life = evolution.main_sequence_lifetime_gy(birth)
        return _star_row(seed, birth, min(STAR_AGE_CAP_GY, 1.5 * life), 0.6, "VII")
    low, high = next((lo, hi) for n, lo, hi in STAR_CLASSES if n == name)
    birth = math.exp(rng.uniform(math.log(low), math.log(high)))
    life = evolution.main_sequence_lifetime_gy(birth)
    return _star_row(seed, birth, min(STAR_AGE_CAP_GY, OLD_DENSITY_AGE_FRACTION * life), birth, "V")


def tree_of(name, seed):
    """`(row, tree)` of the row of class `name` with this seed, `tree` None for a single object (what the scatter passes store)."""
    for cls in PHENOMENON_CLASSES:
        if cls.name == name:
            row = _phenomenon_row(cls, seed)
            return row, scatter_companions.phenomenon_tree(row)
    row = star_row(name, seed)
    return row, scatter_companions.star_tree(row)


# --- The checks of one tree ---------------------------------------------------------------------------------------


def stored_form(tree):
    """`(OrbitNode list, TreeMember list)` of a drawn tree, as the stored rows hold them."""
    nodes = orbit_trees.nodes_from_member_tree(tree, 0.0)
    members = [TreeMember(m.member_no, m.parent_node_no, m.slot, m.mass_solar) for m in tree.members if m.parent_node_no is not None]
    return nodes, members


def _child_mass(tree, node_no, slot):
    """The mass of the child of a node in a slot: a member's, or the subtree's of a child node."""
    for member in tree.members:
        if member.parent_node_no == node_no and member.slot == slot:
            return member.mass_solar
    for node in tree.nodes:
        if node.parent_node_no == node_no and node.slot == slot:
            return node.mass_solar
    return None


def _leaves_under(tree, node_no):
    """The members below a node."""
    found = []
    for member in tree.members:
        if member.parent_node_no == node_no:
            found.append(member)
    for node in tree.nodes:
        if node.parent_node_no == node_no:
            found += _leaves_under(tree, node.node_no)
    return found


def stability_problems(tree):
    """`[str]`: the levels of the tree that fail the stability criteria (`hierarchy_stability.verdicts`), the pairs whose Kozai cycle
    brings the pericenter inside four times the extents of what they hold (`stability.kozai_verdict`, with the radii of the members
    under the pair, which is at most the extent the generator used), and the compact pairs whose Peters time is under the system's
    age."""
    found = []
    if not tree.nodes:
        return found
    nodes, members = stored_form(tree)
    for verdict in hierarchy_stability.verdicts([n for n in nodes], members):
        if not verdict.stable:
            found.append(f"level at node {verdict.node_no} ({verdict.shape}) fails {verdict.criterion}: {verdict.detail}")
    by_no = {node.node_no: node for node in tree.nodes}
    for node in tree.nodes:
        if node.kind == "group":
            continue
        parent = by_no.get(node.parent_node_no)
        m_left, m_right = _child_mass(tree, node.node_no, 0), _child_mass(tree, node.node_no, 1)
        if parent is not None and parent.kind != "group" and m_left is not None and m_right is not None:
            sibling = _child_mass(tree, parent.node_no, 1 - node.slot)
            mutual = stability.mutual_inclination(math.radians(node.inclination_deg), math.radians(node.ascending_node_deg),
                                                  math.radians(parent.inclination_deg), math.radians(parent.ascending_node_deg))
            extent = sum(member.radius_au for member in _leaves_under(tree, node.node_no))
            if sibling is not None and stability.kozai_verdict(m_left, m_right, sibling, node.separation_au, parent.separation_au,
                                                              node.eccentricity, parent.eccentricity, mutual, tree.age_yr, extent) == "unsafe":
                found.append(f"node {node.node_no} about node {parent.node_no}: a Kozai cycle brings the pericenter inside 4 extents")
        leaves = [m for m in tree.members if m.parent_node_no == node.node_no]
        if len(leaves) == 2 and all(m.kind in multistar.COMPACT_KINDS + ("white-dwarf",) for m in leaves):
            if not stability.survives_gravitational_waves(leaves[0].mass_solar, leaves[1].mass_solar, node.separation_au,
                                                          node.eccentricity, tree.age_yr):
                found.append(f"compact pair at node {node.node_no} merges by gravitational waves before the system's age")
    return found


# --- One chunk ----------------------------------------------------------------------------------------------------


@dataclasses.dataclass
class Tally:
    """What the rows of one or more chunks add up to."""
    rows: int = 0
    trees: int = 0
    states: collections.Counter = dataclasses.field(default_factory=collections.Counter)
    star_total: int = 0
    star_squares: int = 0
    stars: int = 0
    draw_total: int = 0
    draw_squares: int = 0
    members: int = 0
    same_class: int = 0
    same_kind_members: int = 0
    both_row_and_member: int = 0
    duplicate_member_seeds: int = 0
    not_a_member_kind: collections.Counter = dataclasses.field(default_factory=collections.Counter)
    determinism_checked: int = 0
    determinism_broken: int = 0
    unstable: int = 0
    unstable_examples: list = dataclasses.field(default_factory=list)
    black_hole_pairs: int = 0
    young_groups_old: int = 0

    def add(self, other):
        for field in dataclasses.fields(self):
            mine, theirs = getattr(self, field.name), getattr(other, field.name)
            if isinstance(mine, list):
                mine[:] = (mine + theirs)[:EXAMPLES]
            elif isinstance(mine, collections.Counter):
                mine.update(theirs)
            else:
                setattr(self, field.name, mine + theirs)
        return self


def _placed_kind_of(name):
    for cls in PHENOMENON_CLASSES:
        if cls.name == name:
            return phenomenon_scatter.PLACED_KINDS.get((cls.kind, cls.subtype))
    return "white-dwarf" if name == WHITE_DWARF else "star"


def _state_of(name, seed, tree):
    """The companion state of a row: `(recycled, state)` for a neutron star row, else `state`, where state is None for a single object."""
    if name == "neutron-star":
        recycled = scatter_companions.screen_unit(seed, "recycled") < phenomenon_scatter.NEUTRON_STAR_RECYCLED_SHARE
        hit = scatter_companions.phenomenon_companion("neutron-star", None, seed)
        return recycled, None if hit is None else hit[2]
    if tree is None or tree.member_count < 2:
        return None
    return tree.members[1].kind if tree.members[1].kind in multistar.COMPACT_KINDS + ("white-dwarf",) else "star"


def run_chunk(job):
    """The `Tally` of one job `(class name, chunk number, rows)`."""
    name, chunk, count = job
    tally = Tally()
    kind = _placed_kind_of(name)
    for seed in rows_of(name, chunk, count):
        tally.rows += 1
        row, tree = tree_of(name, seed)
        state = _state_of(name, seed, tree)
        if name in {c.name for c in PHENOMENON_CLASSES}:
            hit = scatter_companions.phenomenon_companion(next(c.kind for c in PHENOMENON_CLASSES if c.name == name),
                                                          next(c.subtype for c in PHENOMENON_CLASSES if c.name == name), seed)
            tally.states[state if name == "neutron-star" else (None if hit is None else hit[2])] += 1
        else:
            tally.states[state] += 1
        if name not in {c.name for c in PHENOMENON_CLASSES} and name != WHITE_DWARF:
            count = multistar._draw_star_count(row[_BRIGHT["initial_mass_sol"]], multistar.MAX_MEMBERS, 1.0, draw.Stream(f"{seed}:starcount"))
            tally.draw_total += count
            tally.draw_squares += count * count
        if tree is None:
            if name not in {c.name for c in PHENOMENON_CLASSES}:
                tally.star_total += 1
                tally.star_squares += 1
                tally.stars += 1
            continue
        tally.trees += 1
        tally.members += tree.member_count
        if name not in {c.name for c in PHENOMENON_CLASSES}:
            drawn = tree.star_count + tree.dropped
            tally.star_total += drawn
            tally.star_squares += drawn * drawn
            tally.stars += 1
        placed = tree.members[0]
        tally.same_kind_members += sum(1 for m in tree.members[1:] if m.kind == _member_kind_of(kind))
        if kind in ("black-hole", "neutron-star") and _same_class_state(kind, state):
            tally.same_class += 1
        if any(m.seed == seed for m in tree.members[1:]):
            tally.both_row_and_member += 1
        seeds = [m.seed for m in tree.members]
        tally.duplicate_member_seeds += len(seeds) - len(set(seeds))
        for member in tree.members:
            if member.kind not in multistar.MEMBER_KINDS:
                tally.not_a_member_kind[member.kind] += 1
        if kind == "black-hole" and placed.kind == "black-hole":
            tally.black_hole_pairs += sum(1 for m in tree.members[1:] if m.kind == "black-hole")
        problems = stability_problems(tree)
        if problems:
            tally.unstable += 1
            if len(tally.unstable_examples) < EXAMPLES:
                tally.unstable_examples.append(f"{name} seed {seed}: {problems[0]}")
        if tally.trees % max(1, 1000 // DETERMINISM_PER_CHUNK) == 0 and tally.determinism_checked < DETERMINISM_PER_CHUNK:
            tally.determinism_checked += 1
            if tree_of(name, seed)[1] != tree:
                tally.determinism_broken += 1
    return tally


def _same_class_state(kind, state):
    """Whether the companion state of a row is its own class (a black hole's black hole, a neutron star's neutron star, whichever table)."""
    return (state[1] if isinstance(state, tuple) else state) == kind


def _member_kind_of(placed_kind):
    """The member kind a same-class companion of a placed kind has."""
    return {"intermediate-mass-black-hole": "black-hole"}.get(placed_kind, placed_kind)


# --- The whole sample ---------------------------------------------------------------------------------------------


def run(workers=None, sample_jobs=None):
    """`{class name: Tally}` for the whole sample, the chunks drawn by forked workers (`workers` up to 4 by default, 1 to run them one after
    another); the chunks are fixed by the seeds so the answer is the one-by-one run's."""
    sample_jobs = jobs() if sample_jobs is None else sample_jobs
    workers = min(4, os.cpu_count() or 1) if workers is None else workers
    if workers <= 1:
        done = [run_chunk(job) for job in sample_jobs]
    else:
        with concurrent.futures.ProcessPoolExecutor(max_workers=workers, mp_context=multiprocessing.get_context("fork")) as pool:
            done = list(pool.map(run_chunk, sample_jobs, chunksize=2))
    by_class = {name: Tally() for name in classes()}
    for job, tally in zip(sample_jobs, done):
        by_class[job[0]].add(tally)
    return by_class


# --- The expectations ---------------------------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class Check:
    """One comparison: `observed` against `expected` within `sigma_count` sigma of `sigma`."""
    name: str
    observed: float
    expected: float
    sigma: float
    k: float
    side: str = "both"

    @property
    def ok(self):
        """`both`: within k sigma; `at most`: not above expected by more than k sigma; `at least`: not below it by more than k sigma."""
        gap = self.observed - self.expected
        if self.side == "at most":
            return gap <= self.k * self.sigma + 1.0e-9
        if self.side == "at least":
            return -gap <= self.k * self.sigma + 1.0e-9
        return abs(gap) <= self.k * self.sigma + 1.0e-9

    def line(self):
        return (f"{'ok  ' if self.ok else 'FAIL'} {self.name}: observed {self.observed:,.6g}, expected {self.expected:,.6g}"
                f" ({abs(self.observed - self.expected) / self.sigma if self.sigma else 0.0:.2f} sigma of {self.sigma:,.4g}, limit {self.k:g})")


def share_check(name, observed, rows, share, k=SIGMA_COUNT):
    return Check(name, observed, rows * share, math.sqrt(max(rows * share * (1.0 - share), 1.0e-12)), k)


def state_table_of(placed_class_name):
    """`{state: odds}` of one `multistar.MEMBER_STATE_BY_PLACED` entry (`"white-dwarf"`, `"star-massive"`, ...)."""
    return {state: values[0] for state, values in multistar.MEMBER_STATE_BY_PLACED[placed_class_name].items()}


def state_table(cls_name, recycled=None):
    """`{state: odds}` of a phenomenon class (and, for a neutron star, of its recycled or young rows)."""
    if cls_name == "neutron-star":
        key = "recycled-pulsar" if recycled else "neutron-star"
    else:
        key = {"black-hole": "black-hole", "intermediate-mass-black-hole": "intermediate-mass-black-hole",
               "planetary-nebula-star": "planetary-nebula-star"}[cls_name]
    return {state: values[0] for state, values in multistar.MEMBER_STATE_BY_PLACED[key].items()}


def companion_state_checks(by_class):
    """The share of rows with each companion state of each phenomenon class against its table, within four sigma."""
    found = []
    for cls in PHENOMENON_CLASSES:
        tally = by_class[cls.name]
        if cls.name == "neutron-star":
            recycled_odds = phenomenon_scatter.NEUTRON_STAR_RECYCLED_SHARE
            young, recycled = (rows * 1.0 for rows in (tally.rows * (1.0 - recycled_odds), tally.rows * recycled_odds))
            young_table, recycled_table = state_table(cls.name, False), state_table(cls.name, True)
            expected = collections.Counter()
            for state, odds in young_table.items():
                expected[state] += odds * (1.0 - recycled_odds)
            for state, odds in recycled_table.items():
                expected[state] += odds * recycled_odds
            counts = collections.Counter()
            for key, count in tally.states.items():
                counts[key[1] if isinstance(key, tuple) else key] += count
            table = expected
        else:
            table = state_table(cls.name)
            counts = tally.states
        for state, odds in sorted(table.items(), key=lambda item: str(item[0])):
            found.append(share_check(f"{cls.name}: rows with a {state} companion", counts[state], tally.rows,
                                      odds if cls.name != "neutron-star" else table[state]))
        none = 1.0 - sum(table.values())
        found.append(share_check(f"{cls.name}: rows with no companion", counts[None], tally.rows, none))
        stray = set(counts) - set(table) - {None}
        found.append(Check(f"{cls.name}: companion states the table has no entry for", len(stray), 0, 1.0, 0.5))
    return found


def star_count_checks(by_class):
    """
    The number of stars of each star class against `multistar.expected_star_count` (the `star-count` table and its tail), the mean over the
    birth masses the class was drawn at:

    - the count each row's seed draws (`multistar._draw_star_count`, the draw the tree makes first) is within `SIGMA_MEAN` of the mean;
    - the stars a stored tree holds are never more than that draw's mean (within the sigma) and are at least `TREE_KEEPS` of it: a tree
      keeps fewer stars than were drawn where the generator takes some away on purpose (a contact binary keeps its other stars with
      the odds of `tuning.CONTACT_TERTIARY_KEEP_SHARE`, a merger leaves one star for two).
    """
    found = []
    for name, low, high in STAR_CLASSES:
        tally = by_class[name]
        masses = [math.exp(math.log(low) + (i + 0.5) / 64.0 * (math.log(high) - math.log(low))) for i in range(64)]
        expected = sum(multistar.expected_star_count(m) for m in masses) / len(masses)
        drawn = tally.draw_total / tally.rows
        sigma = math.sqrt(max(tally.draw_squares / tally.rows - drawn * drawn, 0.0) / tally.rows)
        found.append(Check(f"{name} stars: the draw's mean stars per system", drawn, expected, sigma, SIGMA_MEAN))
        mean = tally.star_total / tally.rows
        sigma = math.sqrt(max(tally.star_squares / tally.rows - mean * mean, 0.0) / tally.rows)
        found.append(Check(f"{name} stars: stars a tree holds, at most the draw's mean", mean, expected, sigma, SIGMA_MEAN, "at most"))
        found.append(Check(f"{name} stars: stars a tree holds, at least {TREE_KEEPS:.0%} of the draw's mean", mean, TREE_KEEPS * expected, sigma,
                           SIGMA_MEAN, "at least"))
    return found


TREE_KEEPS = 0.9
"""float: The least share of the drawn stars a stored tree keeps on average (the contact binaries and mergers take at most a tenth, measured at 3 to 4 percent)."""


def bookkeeping_checks(by_class):
    """Same-class companions per row against the table (so a class's rows and their companions keep its object density), the black hole
    pairs against the stock's range, and that nobody is both a row and a member."""
    found = []
    for cls in PHENOMENON_CLASSES:
        tally = by_class[cls.name]
        placed = phenomenon_scatter.PLACED_KINDS.get((cls.kind, cls.subtype))
        if placed is None or cls.kind == "planetary-nebula":
            continue
        per_row = multistar.same_class_companions(placed, phenomenon_scatter.NEUTRON_STAR_RECYCLED_SHARE if placed == "neutron-star" else 0.0)
        found.append(share_check(f"{cls.name}: same-class companions (rows over 1 + {per_row:.4g} keep the class density)",
                                  tally.same_class, tally.rows, per_row))
        found.append(Check(f"{cls.name}: members of the class's own kind in the trees (the companion states, plus stars that have collapsed"
                           " and tertiaries, so never fewer)", tally.same_kind_members, tally.same_class, math.sqrt(max(tally.same_class, 1)),
                           SIGMA_COUNT, "at least"))
        divisor = phenomenon_scatter.row_density_divisor(cls.kind, cls.subtype)
        found.append(Check(f"{cls.name}: the rows' divisor is 1 plus the table's same-class odds", divisor, 1.0 + per_row, 1.0e-9, 1.0))
    for name in classes():
        tally = by_class[name]
        found.append(Check(f"{name}: rows that are also a member of their own tree", tally.both_row_and_member, 0, 1.0, 0.5))
        found.append(Check(f"{name}: members of one tree that share a seed", tally.duplicate_member_seeds, 0, 1.0, 0.5))
    return found


def tree_checks(by_class):
    """The checks every stored tree passes: stability, the Kozai and Peters screens, determinism, and no member that is not a member kind."""
    found = []
    for name in classes():
        tally = by_class[name]
        found.append(Check(f"{name}: trees failing the stability, Kozai or Peters screen (of {tally.trees:,})", tally.unstable, 0, 1.0, 0.5))
        found.append(Check(f"{name}: redrawn trees that differ (of {tally.determinism_checked:,})", tally.determinism_broken, 0, 1.0, 0.5))
        found.append(Check(f"{name}: members that are a planet, a moon, a belt or any kind not stored", sum(tally.not_a_member_kind.values()), 0,
                           1.0, 0.5))
    return found


def black_hole_pair_checks(by_class):
    """The black hole pairs: the share of stellar black holes with a black hole companion lies between the gravitational-wave stock
    (1.7e-4) and the Lamberts count (2.4e-2) of holes (design note section 9, test 9)."""
    tally = by_class["black-hole"]
    share = tally.black_hole_pairs / max(1, tally.rows)
    low, high = BLACK_HOLE_PAIR_SHARE
    sigma = math.sqrt(max(share * (1.0 - share), 1.0e-12) / max(1, tally.rows))
    return [Check("black holes with a black hole companion lie between the stock and the Lamberts count", share,
                  min(max(share, low), high), sigma, SIGMA_COUNT)]


BLACK_HOLE_PAIR_SHARE = (1.7e-4, 2.4e-2)
"""tuple: The gravitational-wave stock's share of black holes with a black hole companion and the Lamberts count's (multi-star-systems.md section 9)."""


def all_checks(by_class):
    """Every check of the audit over a sample."""
    return (companion_state_checks(by_class) + star_count_checks(by_class) + bookkeeping_checks(by_class) + tree_checks(by_class)
            + black_hole_pair_checks(by_class))


def format_report(by_class, checks=None):
    """The audit's report: the sample, every check and the examples of what broke."""
    checks = all_checks(by_class) if checks is None else checks
    lines = ["Multiplicity audit (TEST.162)", "", "Sample (rows, trees drawn, members):"]
    for name in classes():
        tally = by_class[name]
        lines.append(f"  {name:<32} {tally.rows:>10,} rows {tally.trees:>10,} trees {tally.members:>10,} members")
    lines += ["", "Checks:"] + ["  " + check.line() for check in checks]
    bad = [check for check in checks if not check.ok]
    lines += ["", f"{len(checks)} checks, {len(bad)} failed."]
    for name in classes():
        for example in by_class[name].unstable_examples:
            lines.append("  example: " + example)
    return "\n".join(lines) + "\n"


def report(workers=None):
    """`(text, failed checks)` for the whole sample."""
    by_class = run(workers)
    checks = all_checks(by_class)
    return format_report(by_class, checks), [check for check in checks if not check.ok]


def main(argv=None):
    """`python -m tests.audit.multiplicity [report.txt]`: prints the audit's report and writes it to the file when one is named; exits 1 when a
    check failed."""
    import argparse
    import sys
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument("path", nargs="?")
    options = parser.parse_args(argv)
    text, failed = report()
    sys.stdout.write(text)
    if options.path:
        with open(options.path, "w", encoding="utf-8") as handle:
            handle.write(text)
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
