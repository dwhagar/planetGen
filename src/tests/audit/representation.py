# tests/audit/representation.py

"""
The regions and the expectations the whole-sphere audits share (TEST.160, TEST.184, TEST.185), on the census harness
(`tests.audit.census`).

A real sector holds a handful of systems, and a per-star kind (a comet, a supernova remnant) is rolled at one in a
thousand to one in a billion stars, so a few dozen seeded fills would never show most kinds at all. The audit therefore
fills with the rates dialed up to `BOOSTED_RATE_PER_STAR` for every kind (`boosted_rates`): what is audited is that each
kind's count follows its rate, the stars of the sector and the profile's `kind_factors` of the place, in every region
class, and that a kind the place's gate forbids never appears. The rates themselves, the factors and the registry's
expected numbers have tests of their own (`test_gen_kind_factors.py`, `test_phenomena_registry.py`).

A kind that needs a free spot of its own (a neutron star, a black hole, a brown dwarf, a planetary nebula's central
star) is dropped when the star already there has filled the cell with its Hill sphere, which in a halo or a cluster is
most of it (the tidal sphere is wide where the galaxy pulls little). Those `PLACEMENT_LIMITED` kinds are audited from
above only: never more than the rate allows. In a whole galaxy their rows are placed by the scatters at their own points
and the fill builds them where they are (`phenomenon_scatter`), which `test_phenomenon_scatter.py` covers.
"""

import collections
import math

from planetgen import tuning
from planetgen.galaxy import component_field, profile, remnant_density
from planetgen.generation import run_sector
from tests import test_gen_component_pass as clusters
from tests import test_satellites as satellites
from tests.audit import census as audit

BOOSTED_RATE_PER_STAR = 0.3
"""float: What every per-star kind is rolled at in an audit fill, before the place's factor."""

SECTORS = 100
"""int: How many seeded sectors each region is filled with."""

SYSTEMS_PER_SECTOR = 2
"""int: Systems in each audit sector. A cube holds about seven before the Hill spheres leave no room, and a kind that
adds a system of its own (a planetary nebula, a brown dwarf) needs the room."""

PER_STAR_KINDS = ("rogue-planet", "comet", "brown-dwarf", "molecular-cloud", "supernova-remnant", "planetary-nebula",
                  "neutron-star")
"""tuple: The kinds a fill rolls per star that the audit follows."""

PLACEMENT_LIMITED = frozenset({"brown-dwarf", "planetary-nebula", "neutron-star", "black-hole"})
"""frozenset: The kinds that need a free spot of their own, audited from above only (module docstring)."""

HOSTED = "nebula"
"""str: The kind of an emission or reflection nebula grown around a hot star."""

EDGE = satellites.EDGE
E1 = satellites.E1

CLUSTER = clusters._cluster(name="Audit Cluster", mass=3.0e6, half=5.0, w0=7, center=(12000.0, 3000.0, 9000.0))
"""Component: An old globular cluster in the halo, away from the disk's gas."""

SHAPE = satellites.SHAPE
"""GalaxyShape: The Milky Way's, the one the audits read profiles from."""

_PLACED = satellites._placed()
_LMC = next(c for c in _PLACED if c.name == "Large Magellanic Cloud")
WITH_COMPONENTS = SHAPE._replace(components=component_field.ComponentField(_PLACED + [CLUSTER], EDGE, E1))


def _arm_points(radius_pc=6000.0, height_pc=100.0, steps=72):
    """The crest and the trough of the thin disk's arms on a circle: the angles of the highest and lowest young share."""
    best, worst = None, None
    for step in range(steps):
        angle = 2.0 * math.pi * step / steps
        point = (radius_pc * math.cos(angle), radius_pc * math.sin(angle), height_pc)
        mix = profile.profile_at(point, SHAPE, EDGE, E1).mix
        share = mix.get("young", 0.0) / sum(mix.values())
        if best is None or share > best[0]:
            best = (share, point)
        if worst is None or share < worst[0]:
            worst = (share, point)
    return best[1], worst[1]


_ARM, _INTERARM = _arm_points()

REGION_POINTS = {
    "bulge": ((300.0, 0.0, 50.0), "bulge"),
    "thin-disk-arm": (_ARM, "old"),
    "thin-disk-interarm": (_INTERARM, "old"),
    "thick-disk": ((8128.0, 0.0, 900.0), "old"),
    "halo": ((9000.0, 9000.0, 20000.0), "halo"),
    "satellite": (_LMC.center_pc, "satellite"),
    "globular-cluster": (CLUSTER.center_pc, "cluster"),
    "intergalactic": ((150000.0, 100000.0, 150000.0), "halo"),
}
"""dict: Region class -> (a point of it in the galaxy frame, the population its sector mix is dominated by)."""


def region_cases(names=None, sectors=SECTORS):
    """The `RegionCase`s of the audit, each a sector of `SYSTEMS_PER_SECTOR` systems at the region's point."""
    names = names or tuple(REGION_POINTS)
    return [audit.RegionCase(name, REGION_POINTS[name][0], WITH_COMPONENTS, EDGE, E1, sectors=sectors,
                             num_systems=SYSTEMS_PER_SECTOR) for name in names]


def factors_at(name):
    """The `kind_factors` of the sector profile of region `name`."""
    return profile.profile_at(REGION_POINTS[name][0], WITH_COMPONENTS, EDGE, E1).kind_factors


def dominant_population(name):
    """The population with the largest share of region `name`'s sector mix."""
    mix = profile.profile_at(REGION_POINTS[name][0], WITH_COMPONENTS, EDGE, E1).mix
    return max(mix, key=mix.get)


def boosted_rates(monkeypatch, rate=BOOSTED_RATE_PER_STAR):
    """Sets every per-star kind's base rate to `rate` for the length of the test (the place's factors still apply)."""
    for kind in tuning.PHENOMENON_RATE_SCALE:
        density = tuning.PHENOMENON_DENSITY_PC3.get(kind)
        if density:
            monkeypatch.setitem(tuning.PHENOMENON_RATE_SCALE, kind,
                                rate * tuning.REFERENCE_STELLAR_DENSITY_PC3 / density)
    monkeypatch.setattr(remnant_density, "base_rate_per_star", lambda kind=None, *args, **kwargs: rate)


class Fill:
    """
    The fill callable of a census that also records, per region, the stars each sector held when its phenomena were
    rolled (`run_sector.sector_star_count`, which every per-star rate multiplies).
    """

    def __init__(self, monkeypatch):
        self.stars = collections.Counter()
        self.sectors = collections.Counter()
        self._current = None
        real = run_sector.sector_star_count

        def counting(sector):
            n = real(sector)
            if self._current is not None:
                self.stars[self._current] += n
            return n

        monkeypatch.setattr(run_sector, "sector_star_count", counting)

    def __call__(self, case, seed, args):
        self._current = case.name
        self.sectors[case.name] += 1
        try:
            return audit.fill_region_sector(case, seed, args)
        finally:
            self._current = None


def expectations(fill, names, rate=BOOSTED_RATE_PER_STAR):
    """
    `{region: {kind: per-sector mean}}`: the boosted rate times the stars the sectors held times the region's factor,
    for the kinds that are not `PLACEMENT_LIMITED`, plus the hosted nebula at zero where the region forbids it.
    """
    out = {}
    for name in names:
        factors = factors_at(name)
        per_sector = fill.stars[name] / fill.sectors[name]
        out[name] = {kind: rate * per_sector * factors.get(kind, 1.0)
                     for kind in PER_STAR_KINDS if kind not in PLACEMENT_LIMITED}
    return out


def limited_findings(census, fill, names, rate=BOOSTED_RATE_PER_STAR, k=audit.DEFAULT_SIGMA):
    """The `PLACEMENT_LIMITED` kinds that exceed what the rate allows (more than `k` sigma above it)."""
    findings = []
    for name in names:
        factors = factors_at(name)
        region = census[name]
        for kind in sorted(PLACEMENT_LIMITED):
            expected = rate * fill.stars[name] * factors.get(kind, 1.0)
            observed = region.observed(kind)
            if observed > expected + k * math.sqrt(max(expected, 1.0)):
                findings.append(audit.Finding(name, kind, observed, expected, audit.OUT_OF_BOUNDS))
    return findings


def run(monkeypatch, names=None, sectors=SECTORS):
    """Fills the audit regions with boosted rates; returns `(census, fill, names)`."""
    names = tuple(names or REGION_POINTS)
    boosted_rates(monkeypatch)
    fill = Fill(monkeypatch)
    census = audit.run_census(region_cases(names, sectors), fill=fill)
    return census, fill, names


def report(sectors=SECTORS):
    """
    Runs the audit and returns `(text, findings)`: the census report of `audit.format_report` for every region, and the
    kinds that broke a rule, the `PLACEMENT_LIMITED` ones included (listed under the report).
    """
    import pytest
    patch = pytest.MonkeyPatch()
    try:
        census, fill, names = run(patch, sectors=sectors)
    finally:
        patch.undo()
    wanted = expectations(fill, names)
    findings = audit.compare(census, wanted)
    limited = limited_findings(census, fill, names)
    text = audit.format_report(census, wanted, findings)
    text += "\nPlacement-limited kinds (never more than the rate allows): " + (
        "all within the rate\n" if not limited else "; ".join(
            f"{f.region}/{f.label()} observed {f.observed}, at most {f.expected:.1f}" for f in limited) + "\n")
    return text, findings + limited


def main(argv=None):
    """`python -m tests.audit.representation [report.txt] [--sectors N]`: prints the representation audit's report and
    writes it to the file when one is named; exits 1 when a count broke a rule."""
    import argparse
    import sys
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument("path", nargs="?")
    parser.add_argument("--sectors", type=int, default=SECTORS)
    options = parser.parse_args(argv)
    text, findings = report(options.sectors)
    sys.stdout.write(text)
    if options.path:
        with open(options.path, "w", encoding="utf-8") as handle:
            handle.write(text)
    return 1 if findings else 0


if __name__ == "__main__":
    raise SystemExit(main())
