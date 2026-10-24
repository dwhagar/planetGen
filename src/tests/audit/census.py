# tests/audit/census.py

"""
The census audit harness (TEST.186): fill seeded sectors in each region class, count what the fills made by stage and
kind, compare the counts with an expectation within k Poisson sigma, and print one report.

The audits that read it add checks to it and no counters of their own: the representation audit (TEST.160), the
fill-stage and system-pass tests (TEST.161), the multiplicity audit (TEST.162), the catalogue type audit (TEST.164),
the whole-sphere tests (TEST.184, TEST.185) and the expected-numbers checks (TEST.196, TEST.206). The harness
knows nothing about which kinds should appear where; a caller hands it the expected number per sector for each region
and kind, and it does the counting, the Poisson test, the "expected 20 or more appears at least once" rule and the report.

The pieces:

- `RegionCase`: one region class to fill (a position in a galaxy shape, how many systems, which seeds).
- `fill_region_sector(case, seed, args)`: the default fill of one sector of a region, DB-free, through
  `run_sector.generate_sector` with a `FillContext` of the position's profile. A caller that needs another fill
  (scatter rows, a stored skeleton) passes its own `fill` to `run_census`.
- `run_census(cases, ...)`: fills every seed of every case and returns a `Census`.
- `within_sigma`, `must_appear`: the two rules.
- `compare(census, expectations)`: the rules over the whole census, as a list of `Finding`.
- `format_report(census, expectations, findings)`: the one report.

A region's seeds are fixed by its name and the sector number (`seed_for`), so a census is the same on every run and
the cases of a census do not move each other's draws. Counts come from `sector.stage_counts` (the fill's stage I,
`generation.fill_stages`), keyed by `(stage, kind)`.

Counts are tested as Poisson numbers: the sectors of a region add up to a Poisson count whose mean is the sector
mean times the number of sectors. A clumped kind (a cluster's members) has more variance than that, so the
caller gives such a kind the expectation of its parents and counts the parents.
"""

import collections
import copy
import math
from dataclasses import dataclass, field

from planetgen.galaxy.profile import profile_at
from planetgen.generation import bright_stars as brightStars
from planetgen.generation import run_galaxy, run_sector
from planetgen.physics.units import pc_to_ly
from planetgen.util import draw

DEFAULT_SIGMA = 4.0
"""float: How many Poisson sigma a count may stray from its expectation."""

APPEAR_AT = 20.0
"""float: A kind expected this many times or more over a region's sectors appears at least once."""

OUT_OF_BOUNDS = "out of bounds"
NEVER_APPEARED = "never appeared"
FORBIDDEN = "forbidden"
"""str: The kinds of `Finding.problem`."""


@dataclass(frozen=True)
class RegionCase:
    """
    One region class to fill.

    Attributes:
        name (str): The region class (`"thin-disk-arm"`, `"halo"`); also the root of its seeds.
        center_pc (tuple): The sector's centre in the galaxy frame, parsecs.
        shape (GalaxyShape): The galaxy shape its profile is read from.
        edge_pc (float): The sector edge, parsecs.
        expected_at_density_1 (float): The systems a sector of this edge holds at relative density 1
            (`SectorProfile.expected_systems` scales it).
        sectors (int): How many seeded sectors to fill.
        num_systems (int or None): A fixed system count per sector; `None` draws it from the density.
        density (float or None): The relative density of the draw; `None` takes the profile's.
        fill_kwargs (dict): Extra arguments of `bright_stars.FillContext` (`phenomenon_rows`, `below_cut` and the
            rest), for a case that fills with a scatter.
    """

    name: str
    center_pc: tuple
    shape: object
    edge_pc: float
    expected_at_density_1: float
    sectors: int = 8
    num_systems: int | None = None
    density: float | None = None
    fill_kwargs: dict = field(default_factory=dict)


def seed_for(region, index):
    """The seed of sector `index` of `region`: fixed by the two, so every run fills the same sectors."""
    return f"census:{region}:{index}"


@dataclass
class RegionCensus:
    """
    What the fills of one region made.

    Attributes:
        name (str): The region class.
        sectors (int): The sectors filled.
        per_sector (list): One `Counter` of `(stage, kind) -> n` per sector, in seed order.
        totals (Counter): `(stage, kind) -> n` over all the sectors.
    """

    name: str
    sectors: int = 0
    per_sector: list = field(default_factory=list)
    totals: collections.Counter = field(default_factory=collections.Counter)

    def add(self, counts):
        """Adds one sector's `(stage, kind) -> n` counts."""
        counts = collections.Counter(counts)
        self.per_sector.append(counts)
        self.totals.update(counts)
        self.sectors += 1

    def kinds(self):
        """The kinds that were made, sorted."""
        return sorted({kind for _stage, kind in self.totals})

    def by_kind(self):
        """`{kind: n}` over all the stages."""
        out = collections.Counter()
        for (_stage, kind), n in self.totals.items():
            out[kind] += n
        return out

    def observed(self, key):
        """The count of `key`: a `(stage, kind)` pair, or a kind alone to add its stages."""
        if isinstance(key, tuple):
            return self.totals.get(key, 0)
        return self.by_kind().get(key, 0)

    def sectors_with(self, key):
        """In how many sectors `key` (a `(stage, kind)` pair or a kind) was made at least once."""
        if isinstance(key, tuple):
            return sum(1 for counts in self.per_sector if counts.get(key, 0) > 0)
        return sum(1 for counts in self.per_sector if any(k == key and n > 0 for (_s, k), n in counts.items()))


class Census:
    """The regions of a census by name, in the order they were filled."""

    def __init__(self, regions=()):
        self.regions = {region.name: region for region in regions}

    def __getitem__(self, name):
        return self.regions[name]

    def __iter__(self):
        return iter(self.regions.values())

    def __len__(self):
        return len(self.regions)


def _distance_ly(center_pc):
    return pc_to_ly(math.sqrt(sum(c * c for c in center_pc)))


def default_args():
    """The shared generation arguments of a fill (`run_galaxy._default_generation_args`)."""
    return run_galaxy._default_generation_args()


def fill_region_sector(case, seed, args=None):
    """
    Fills one sector of `case` on `seed`: the sector's profile at the case's position, then
    `run_sector.generate_sector` with it, the way a galaxy sector is filled (no database).

    Args:
        case (RegionCase): The region.
        seed: What the fill's draws are bound to (`seed_for`).
        args (argparse.Namespace, optional): The shared generation arguments; built when not given.

    Returns:
        SpaceSector: The filled sector; its `stage_counts` are the census.
    """
    args = copy.copy(args if args is not None else default_args())
    profile = profile_at(case.center_pc, case.shape, case.edge_pc, case.expected_at_density_1)
    if case.num_systems is not None:
        args.num_systems, args.density = case.num_systems, None
    else:
        args.num_systems, args.density = None, profile.density if case.density is None else case.density
    context = brightStars.FillContext(profile, **case.fill_kwargs)
    with draw.bound(seed):
        _name, sector = run_sector.generate_sector(args, galactic_center_dist_ly=_distance_ly(case.center_pc),
                                                   fill=context)
    return sector


def run_census(cases, fill=fill_region_sector, args=None):
    """
    Fills every seeded sector of every case and counts them.

    Args:
        cases (iterable): `RegionCase`s, or anything with `name` and `sectors`, for a custom `fill`.
        fill (callable): `fill(case, seed, args)` -> a sector with `stage_counts`.
        args (argparse.Namespace, optional): The shared generation arguments, built once when not given.

    Returns:
        Census: One `RegionCensus` per case.
    """
    args = args if args is not None else default_args()
    regions = []
    for case in cases:
        region = RegionCensus(case.name)
        for index in range(case.sectors):
            sector = fill(case, seed_for(case.name, index), args)
            region.add(sector.stage_counts)
        regions.append(region)
    return Census(regions)


def sigma_distance(observed, expected):
    """How many Poisson sigma `observed` lies from `expected` (positive above); infinite for a count where none was
    expected, zero when both are zero."""
    if expected <= 0.0:
        return 0.0 if observed == 0 else math.inf
    return (observed - expected) / math.sqrt(expected)


def within_sigma(observed, expected, k=DEFAULT_SIGMA):
    """Whether `observed` is within `k` Poisson sigma of `expected`; where `expected` is zero, only zero is."""
    if expected <= 0.0:
        return observed == 0
    return abs(observed - expected) <= k * math.sqrt(expected)


def must_appear(expected, observed, at=APPEAR_AT):
    """Whether the rule 'a kind expected `at` times or more appears at least once' is broken."""
    return expected >= at and observed == 0


@dataclass(frozen=True)
class Finding:
    """
    A count the census rules refuse.

    Attributes:
        region (str): The region class.
        key: The kind, or the `(stage, kind)` pair the expectation named.
        observed (int): The count over the region's sectors.
        expected (float): The expectation over the region's sectors.
        problem (str): `OUT_OF_BOUNDS`, `NEVER_APPEARED` or `FORBIDDEN` (a count where the gate expects none).
    """

    region: str
    key: object
    observed: int
    expected: float
    problem: str

    def label(self):
        """The key as text: `"scattered/black-hole"` or `"black-hole"`."""
        return "/".join(self.key) if isinstance(self.key, tuple) else str(self.key)


def expected_total(per_sector, sectors):
    """The expectation over `sectors` sectors of a kind that averages `per_sector` per sector."""
    return per_sector * sectors


def compare(census, expectations, k=DEFAULT_SIGMA, appear_at=APPEAR_AT):
    """
    The census rules over every region and every expected key.

    Args:
        census (Census): What the fills made.
        expectations (dict): `{region: {key: per_sector_mean}}`; a key is a kind (its stages added) or a
            `(stage, kind)` pair. A mean of zero says the region's gate forbids the kind. A region with no entry
            is not compared.
        k (float): The Poisson sigma allowed.
        appear_at (float): The expected total from which a kind must appear.

    Returns:
        list: The `Finding`s, in region and key order; empty when every count passes.
    """
    findings = []
    for region in census:
        for key, mean in sorted(expectations.get(region.name, {}).items(), key=lambda item: str(item[0])):
            expected = expected_total(mean, region.sectors)
            observed = region.observed(key)
            if expected <= 0.0:
                if observed:
                    findings.append(Finding(region.name, key, observed, expected, FORBIDDEN))
            elif must_appear(expected, observed, appear_at):
                findings.append(Finding(region.name, key, observed, expected, NEVER_APPEARED))
            elif not within_sigma(observed, expected, k):
                findings.append(Finding(region.name, key, observed, expected, OUT_OF_BOUNDS))
    return findings


def unlisted(census, expectations):
    """`{region: [kinds made that no expectation names]}`: the report lists them; they are not failures."""
    out = {}
    for region in census:
        if region.name not in expectations:
            continue
        named = {key[1] if isinstance(key, tuple) else key for key in expectations.get(region.name, {})}
        extra = [kind for kind in region.kinds() if kind not in named]
        if extra:
            out[region.name] = extra
    return out


def format_report(census, expectations=None, findings=None, k=DEFAULT_SIGMA):
    """
    The one report: for each region its sectors and, per expected key, the observed and expected counts with the
    distance in sigma and the verdict; then the kinds nothing expected; then a summary line.

    Args:
        census (Census): What the fills made.
        expectations (dict, optional): As `compare`. Without it the report lists the counts only.
        findings (list, optional): The result of `compare`; computed when not given.
        k (float): The sigma used when `findings` is computed.

    Returns:
        str: The report text.
    """
    expectations = expectations or {}
    findings = compare(census, expectations, k) if findings is None else findings
    bad = {(f.region, f.key): f.problem for f in findings}
    lines = ["Census audit", ""]
    for region in census:
        lines.append(f"{region.name}: {region.sectors} sectors")
        wanted = expectations.get(region.name, {})
        if wanted:
            width = max(len(Finding(region.name, key, 0, 0.0, "").label()) for key in wanted)
            for key in sorted(wanted, key=str):
                expected = expected_total(wanted[key], region.sectors)
                observed = region.observed(key)
                distance = sigma_distance(observed, expected)
                shown = "inf" if math.isinf(distance) else f"{distance:+.1f}"
                verdict = bad.get((region.name, key), "ok")
                label = Finding(region.name, key, 0, 0.0, "").label()
                lines.append(f"  {label:<{width}}  observed {observed:>8}  expected {expected:>10.2f}  "
                             f"{shown:>6} sigma  {verdict}")
        else:
            for kind, n in sorted(region.by_kind().items()):
                lines.append(f"  {kind}  observed {n:>8}")
        extra = unlisted(census, expectations).get(region.name)
        if extra and wanted:
            lines.append("  made, with no expectation: " + ", ".join(extra))
        lines.append("")
    lines.append(f"{len(findings)} of the expected counts broke a rule" if findings
                 else "Every expected count is within bounds")
    return "\n".join(lines) + "\n"


def write_report(path, census, expectations=None, findings=None, k=DEFAULT_SIGMA):
    """Writes `format_report` to `path` and returns the text."""
    text = format_report(census, expectations, findings, k)
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(text)
    return text
