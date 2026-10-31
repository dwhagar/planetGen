# tests/audit/prevalence.py

"""
The prevalence census harness (TEST.430): seeded draws of planets and moons per star type and zone, counted into class and
overall-colour shares, and compared with a stored file of expected ranges (`tests/data/prevalence_expectations.json`).

Two kinds of draw feed it:

- **Bare planet draws** (`draw_cell`): one `Planet` at a chosen orbit of a star of a given type, class and moons drawn as the
  system draw would; this is the cheap, many-sample part (a cell is a star type by a zone).
- **Full systems** (`draw_systems`): `StarSystem` draws of one star type, which add what only a whole system decides (belts,
  comets, the planet count, the habitable-world retry).

Both are seeded from the cell's name (`seed_for`), so a run repeats and a cell does not move another cell's draws.

The expectations file holds, per measure, a `min` and a `max`, a `kind` (`observed` rate with its source, `norm` researched, or
`stand-in` with the tier and tweak level of docs/design/stand-in-values.md) and, where the generator is known to be off, an
`anomaly` item number. `compare` returns a list of `Finding`; `format_report` writes the latest-run report.

Reading the figures: shares are fractions of the cell's bodies (0 to 1). A class absent from a cell's expectations must not
appear there (`barred`); a class with a `min` above zero must appear (`missing`).
"""

import collections
import json
import multiprocessing
import os
import pathlib
import signal
import zlib
from dataclasses import dataclass

STAR_TYPES = ("M2V", "K5V", "G2V", "F5V", "A0V", "K0III")
"""tuple: The census stars: four dwarfs, an A dwarf and a giant (Boss, 7:21 PM EDT 2026-10-28)."""

ZONES = ("h", "e", "c")
"""tuple: The zones, hot, ecosphere and cold (`Planet.zone`)."""

ZONE_NAMES = {"h": "hot", "e": "ecosphere", "c": "cold"}

TIERS = ("Blue", "Green", "Yellow", "Red")

CEILINGS = ("unlimited", "microbial", "none")
"""tuple: What the life class of a scored body allows (`ceiling_of`): complex life allowed (no ceiling), microbial life only, or no life."""

GROUPS = (("classes", "class"), ("tiers", "tier"), ("ceilings", "ceiling"))
"""tuple: The measures a cell's ranges hold: `(key in the file, name in a finding)`."""

EXPECTATIONS_PATH = pathlib.Path(__file__).resolve().parent.parent / "data" / "prevalence_expectations.json"

REPORT_NAME = "class-census-latest.md"

DRAW_TIME_LIMIT_S = 45
"""int: Seconds one planet draw may take before the census abandons it and counts it (`abandoned`); a draw that never returns
would stall the whole census, and is itself a finding (a generation hang, see `Finding` kind `hang`)."""


class _Abandoned(BaseException):
    """Raised by the alarm when one draw ran past `DRAW_TIME_LIMIT_S`."""


def seed_for(*parts):
    """The seed of a cell from its name: the same on every run and different for each cell."""
    return zlib.crc32(":".join(str(part) for part in parts).encode()) or 1


def _orbit_au(zone, habitable_zone, rng):
    """A census orbit in a zone: log-uniform from a fifth of the inner edge to the inner edge for hot, uniform inside the
    ecosphere, log-uniform from the outer edge to 20 times it for cold."""
    low, high = habitable_zone
    if zone == "h":
        return low * 0.2 * (0.9 / 0.2) ** rng.random()
    if zone == "e":
        return low + (high - low) * rng.random()
    return high * 1.1 * (20.0 / 1.1) ** rng.random()


def ceiling_of(body):
    """`CEILINGS` entry of a scored body: no life class is `none`, a class whose complexity index stays under the multicellular threshold is
    `microbial`, any other `unlimited`."""
    from planetgen import tuning

    if body.life_class is None:
        return "none"
    if body.life_complexity is not None and body.life_complexity < tuning.LIFE_COMPLEXITY["thresholds"]["multicellularity"]:
        return "microbial"
    return "unlimited"


def _body_record(body, star_type, zone, kind, life):
    from planetgen.physics import habitability

    score = getattr(body, "phi4", None)
    return {
        "star": star_type, "zone": zone, "kind": kind, "cls": body.planet_class,
        "tier": habitability.tier(score) if score is not None else None,
        "ceiling": ceiling_of(body) if score is not None else None,
        "habitable": habitability.is_habitable(score),
    }


def finish_bodies(system, planets):
    """Runs on bare planets the passes `StarSystem.__init__` runs on a finished system (snow line, rings, insolation, interior heat,
    life data), so each body has the score, ceiling and life a drawn system's bodies have. Mirrors the "snow line pass" to the
    "life data pass" of `generation/system.py`; a pass added there that changes a score belongs here too."""
    from planetgen.generation import life
    from planetgen.physics import insolation, interior_heat, magnetism, rings

    system.planets = list(planets)
    system.secondary_planets = []
    system.member_planets = []
    for planet in planets:
        planet.born_beyond_snow_line = planet.distance >= planet.star.snow_line_au
        for moon in planet.moons:
            moon.born_beyond_snow_line = planet.born_beyond_snow_line
    for planet in planets:
        rings.apply(planet)
        for moon in planet.moons:
            rings.apply(moon, planet)
    insolation.refresh_system(system)
    for planet in planets:
        interior_heat.update_heat(planet)
        life.apply_life_data(planet)
        for moon in planet.moons:
            interior_heat.update_heat(moon)
            magnetism.generate_induced_field(moon, planet)
            life.apply_life_data(moon)


def draw_cell(star_type, zone, count):
    """`count` bare planets (and their moons) at orbits of `zone` around a `star_type` star, finished as a system's bodies are
    (`finish_bodies`), as a list of body records."""
    import random

    from planetgen.util import draw
    from planetgen.generation import life
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.planet import Planet
    from planetgen.generation.system import StarSystem

    draw.set_run_seed(seed_for("cell", star_type, zone))
    rng = random.Random(seed_for("orbit", star_type, zone))
    config = SystemConfig()
    config.STAR_TYPE = star_type
    config.PLANETS = False
    system = StarSystem(system_config=config)
    star = system.primary_star
    zone_edges = star.habitable_zone
    planets, abandoned, raised = [], 0, []

    def _alarm(*_):
        raise _Abandoned()

    previous = signal.signal(signal.SIGALRM, _alarm)
    try:
        for _ in range(count):
            orbit = _orbit_au(zone, zone_edges, rng)
            signal.alarm(DRAW_TIME_LIMIT_S)
            try:
                planets.append(Planet(config, star, zone_edges, orbit))
            except _Abandoned:
                abandoned += 1
            except Exception as error:  # a crash is a finding, not the end of the census
                raised.append(f"{type(error).__name__}: {error}")
            finally:
                signal.alarm(0)
        signal.alarm(DRAW_TIME_LIMIT_S * max(len(planets), 1))
        try:
            finish_bodies(system, planets)
        finally:
            signal.alarm(0)
    finally:
        signal.signal(signal.SIGALRM, previous)
    records = [{"star": star_type, "zone": zone, "kind": "abandoned", "cls": None, "tier": None, "ceiling": None,
                "habitable": False}] * abandoned
    records += [{"star": star_type, "zone": zone, "kind": "error", "cls": message[:120], "tier": None, "ceiling": None,
                 "habitable": False} for message in raised]
    for planet in planets:
        records.append(_body_record(planet, star_type, planet.zone, "planet", life))
        for moon in planet.moons or []:
            records.append(_body_record(moon, star_type, planet.zone, "moon", life))
    return records


def _cell_job(job):
    return draw_cell(*job)


def draw_cells(count, star_types=STAR_TYPES, zones=ZONES, workers=None):
    """Every cell's records, drawn on `workers` processes (default: the cores, at most the cells)."""
    jobs = [(star, zone, count) for star in star_types for zone in zones]
    workers = min(workers or os.cpu_count() or 1, len(jobs))
    # always in worker processes: the time guard is an alarm, which a test run's own timeout (also an alarm) must not share
    with multiprocessing.get_context("fork").Pool(workers) as pool:
        results = pool.map(_cell_job, jobs, chunksize=1)
    return [record for cell in results for record in cell]


def tally(records):
    """Counts by `(star, zone, kind)`: `{"n": bodies, "cls": Counter, "tier": Counter, "habitable_by_cls": Counter,
    "ceiling": Counter}`. Bodies with no score (gas and ice giants) are in `cls` and not in `tier`."""
    cells = {}
    for rec in records:
        cell = cells.setdefault((rec["star"], rec["zone"], rec["kind"]), {
            "n": 0, "cls": collections.Counter(), "tier": collections.Counter(),
            "habitable_by_cls": collections.Counter(), "ceiling": collections.Counter()})
        cell["n"] += 1
        cell["cls"][rec["cls"]] += 1
        if rec["tier"] is not None:
            cell["tier"][rec["tier"]] += 1
            cell["ceiling"][rec["ceiling"] or "unlimited"] += 1
        if rec["habitable"]:
            cell["habitable_by_cls"][rec["cls"]] += 1
    return cells


def pooled(cells, star=None, zone=None, kind=None):
    """The cells matching `star`, `zone` and `kind` (None for any) added together, as `tally` shapes one."""
    total = {"n": 0, "cls": collections.Counter(), "tier": collections.Counter(),
             "habitable_by_cls": collections.Counter(), "ceiling": collections.Counter()}
    for (s, z, k), cell in cells.items():
        if (star is None or s == star) and (zone is None or z == zone) and (kind is None or k == kind):
            total["n"] += cell["n"]
            for key in ("cls", "tier", "habitable_by_cls", "ceiling"):
                total[key].update(cell[key])
    return total


# --- The expectations ------------------------------------------------------------------------------------------------


def load_expectations(path=EXPECTATIONS_PATH):
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)


@dataclass
class Finding:
    """One flag: `kind` is `out-of-range`, `missing` or `barred`; `known` names the item that explains it, if any."""
    kind: str
    cell: str
    measure: str
    value: float
    low: float
    high: float
    source: str = ""
    known: str = ""

    def line(self):
        known = f" (known: {self.known})" if self.known else ""
        return (f"{self.kind.upper()} {self.cell} {self.measure}: {self.value:.4f} not in "
                f"[{self.low:.4f}, {self.high:.4f}] [{self.source}]{known}")


def _shares(counter, whole):
    return {key: count / whole for key, count in counter.items()} if whole else {}


def cell_key(star, zone, kind):
    return f"{star}/{ZONE_NAMES[zone]}/{kind}"


def compare(cells, expectations):
    """Every flag over the census: for each expected cell, each class share and each tier share against its range."""
    findings = []
    for key, spec in expectations.get("cells", {}).items():
        star, zone_name, kind = key.split("/")
        zone = next(z for z, name in ZONE_NAMES.items() if name == zone_name)
        cell = cells.get((star, zone, kind))
        if not cell or not cell["n"]:
            findings.append(Finding("missing", key, "cell", 0.0, 1.0, 1.0, "no bodies drawn"))
            continue
        scored = sum(cell["tier"].values())
        shares = {"classes": _shares(cell["cls"], cell["n"]), "tiers": _shares(cell["tier"], scored),
                  "ceilings": _shares(cell["ceiling"], scored)}
        barred = set(spec.get("barred", ()))
        for cls, share in shares["classes"].items():
            if cls in barred:
                findings.append(Finding("barred", key, f"class {cls}", share, 0.0, 0.0, spec.get("source", "")))
        for group, prefix in GROUPS:
            for name, bounds in spec.get(group, {}).items():
                observed = shares[group].get(name, 0.0)
                low, high = bounds["min"], bounds["max"]
                if observed < low or observed > high:
                    kind_flag = "missing" if observed == 0.0 and low > 0 else "out-of-range"
                    findings.append(Finding(kind_flag, key, f"{prefix} {name}", observed, low, high,
                                            f"{bounds.get('kind', '?')}: {bounds.get('source', '')}", bounds.get("anomaly", "")))
    for (star, zone, kind), cell in sorted(cells.items()):
        if kind == "error" and cell["n"]:
            findings.append(Finding("error", cell_key(star, zone, kind), "draws that raised", float(cell["n"]), 0.0, 0.0,
                                    "; ".join(sorted(cell["cls"]))[:300], expectations.get("errors", {}).get("anomaly", "")))
        if kind == "abandoned" and cell["n"]:
            findings.append(Finding("hang", cell_key(star, zone, kind), "draws over the time limit", float(cell["n"]), 0.0, 0.0,
                                    f"a planet draw ran past {DRAW_TIME_LIMIT_S} s",
                                    expectations.get("hangs", {}).get("anomaly", "")))
    return findings


def stale_anomalies(expectations, findings):
    """The ranges tagged with a known anomaly that the census no longer breaks (`cell group name (item)`): the fault is fixed, so
    the tag can go."""
    broken = {(f.cell, f.measure) for f in findings}
    stale = []
    for key, spec in expectations.get("cells", {}).items():
        for group, prefix in GROUPS:
            for name, bounds in spec.get(group, {}).items():
                if bounds.get("anomaly") and (key, f"{prefix} {name}") not in broken:
                    stale.append(f"{key} {prefix} {name} ({bounds['anomaly']})")
    return stale


def format_report(cells, findings, header=""):
    """The latest-run report: one table per star type and zone, the flags, the known anomalies apart."""
    lines = ["# Class census (TEST.430), latest run", ""]
    if header:
        lines += [header, ""]
    unexplained = [f for f in findings if not f.known]
    explained = [f for f in findings if f.known]
    lines += [f"Flags: {len(unexplained)} new, {len(explained)} known anomalies.", ""]
    for finding in unexplained:
        lines.append(f"- {finding.line()}")
    if explained:
        lines += ["", "Known anomalies:", ""]
        lines += [f"- {finding.line()}" for finding in explained]
    lines += ["", "## Shares", ""]
    for star in STAR_TYPES:
        for zone in ZONES:
            for kind in ("planet", "moon"):
                cell = cells.get((star, zone, kind))
                if not cell or not cell["n"]:
                    continue
                classes = ", ".join(f"{c} {n / cell['n']:.1%}" for c, n in sorted(cell["cls"].items(), key=lambda kv: -kv[1]))
                scored = sum(cell["tier"].values())
                tiers = ", ".join(f"{t} {cell['tier'][t] / scored:.1%}" for t in TIERS) if scored else "no scores"
                life = ", ".join(f"{c} {cell['ceiling'][c] / scored:.1%}" for c in CEILINGS) if scored else "no scores"
                lines.append(f"- {cell_key(star, zone, kind)} (n={cell['n']}): {classes}; colour: {tiers}; life: {life}")
                if cell["habitable_by_cls"]:
                    habitable = ", ".join(f"{c} {n} of {cell['cls'][c]}" for c, n in sorted(cell["habitable_by_cls"].items()))
                    lines.append(f"  - Blue or Green by class: {habitable}")
    return "\n".join(lines) + "\n"


def write_report(text, directory=None, name=REPORT_NAME):
    """Writes the report as `name` (the class census's by default) in `directory`, `$PLANETGEN_REPORT_DIR` or `.test-reports/`."""
    directory = pathlib.Path(directory or os.environ.get("PLANETGEN_REPORT_DIR") or ".test-reports")
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / name
    path.write_text(text, encoding="utf-8")
    return path
