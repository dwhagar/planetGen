# stellarObjects/plausibility.py

"""
Physical-plausibility anomaly finder.

See TODO.md ("Future ideas" -> "Physical-plausibility test suite (anomaly
finder)") for the original ask. This module is the reusable engine; the
`src/tests/physical_plausibility_cli.py` script is the CLI wrapper for running a
large batch and reading a human report, and `tests/test_physical_plausibility.py`
is the opt-in (`@pytest.mark.slow`) pytest wrapper around the hard-invariant
half of it. Splitting it this way answers the TODO's open question directly:

    "whether this lives in tests/ as a slow/opt-in suite or as a separate
    standalone script"

...with "both, for different halves of the job": a large statistical batch
run producing a report meant for a human to read doesn't fit naturally into
a pass/fail assertion, but a set of always-true hard invariants absolutely
should be a regression-gating test. See the two-tier design below.

Two-tier "sane range" design
-----------------------------
1. Hard physical invariants -- always true, and either analytically derived
   from the generator's own declared per-class data
   (`program_constants.PLANET_CLASSES`, `physical_constants.PLANET_DENSITY`,
   etc), or simple non-negotiable physical facts (a temperature or pressure
   can't be negative or non-finite). Surface gravity is a multilinear
   function of (radius, density) for terrestrial classes, or, for gas
   giants, a function of (radius, core density, atmosphere density,
   core/atmosphere ratio) that is monotonic (though not multilinear -- see
   `theoretical_gravity_bounds_g`'s docstring) in each of those variables
   individually -- either property guarantees a function's extrema over a
   box occur at the box's corners, so `theoretical_gravity_bounds_g` below
   evaluates every corner of the class's declared ranges and takes the
   min/max, rather than hand-picking a number. Any generated value outside
   these bounds is unambiguously a bug (in the generator, or in this
   derivation -- either way worth reporting). These gate the pytest test.

2. Statistical outliers -- for surface_temperature and atmospheric_pressure,
   no such closed-form bound is tractable (both depend on the host star's
   luminosity, orbital distance, a randomized albedo, and a greenhouse
   factor derived from randomized atmospheric composition). Sane ranges for
   these are instead derived from a large generated sample's own
   distribution, per planet class, using Tukey's fences (values beyond
   `k * IQR` outside the inter-quartile range, k=3.0 -- the conservative
   "far outlier" convention rather than the more trigger-happy k=1.5 --
   since host stars are deliberately sampled across the whole main-sequence
   spectral range and some real, non-buggy spread is expected). These are
   reported for human review, not asserted to be zero: doing that would
   reintroduce exactly the kind of band-aid the now-disabled Class M/P
   clamps used to be (see planetPhysics.py).

Host star sampling ("per star class where relevant")
------------------------------------------------------
Rather than a single fixed G2V host star (as tests/test_planets.py uses,
reasonably, for deterministic regression checks), each sample draws its host
star's spectral type from a broad main-sequence (Yerkes V) grid spanning
O..M. This is what makes the tool able to catch a bug that only manifests
for certain host star classes -- e.g. a greenhouse-effect blowup that only
shows up for a hot, luminous primary. Evolved (non-V) Yerkes classes are
deliberately excluded here: their far larger luminosity swings would
dominate a class's temperature distribution with real (not buggy) variance,
and tests/test_star_matrix.py already exhaustively sanity-checks star
physical properties across every spectral/Yerkes combination.
"""

import itertools
import math
import random
import statistics
from collections import defaultdict

from . import physical_constants as pc
from . import program_constants as prog_c
from .config import SystemConfig
from .planetData import Planet
from .starData import Star

ZONE_CHARS = "hec"

# Main-sequence-only, spanning the full spectral range, several subclasses
# per letter so temperature/luminosity aren't fixed to one boundary value.
SPECTRAL_CLASSES = ("O", "B", "A", "F", "G", "K", "M")
SUBCLASSES = (0, 3, 5, 7, 9)
HOST_STAR_TYPES = tuple(f"{spec}{sub}V" for spec in SPECTRAL_CLASSES for sub in SUBCLASSES)

# (class, zone) pairs PLANET_CLASSES actually declares valid -- same
# derivation tests/test_planets.py uses, kept independent here on purpose
# (importing it from the test module would make this importable-anywhere
# module depend on the tests package).
VALID_CLASS_ZONE_PAIRS = [
    (cls, zone)
    for cls, data in prog_c.PLANET_CLASSES.items()
    for zone in ZONE_CHARS
    if data[zone]
]

# The physical quantities this tool checks. surface_temperature and
# atmospheric_pressure are the two the atmospheric-pressure unit bug
# ([5.3.0], see CHANGELOG.md) actually broke; gravity is included because
# it's the other quantity the now-disabled Class M/P clamps used to paper
# over; scale_height is included as a cheap extra since it feeds directly
# into the pressure calculation; density is included because gas giants'
# core/atmosphere density blend (generate_planet_properties in
# planetPhysics.py) is a separate source of implausible values from the
# gravity/pressure issues above (see the manual habitability/atmosphere
# sanity review -- docs/analysis/habitability-atmosphere-sanity-review.md
# on branch worktree-agent-a75772018725a8694 as of this writing, not yet
# merged to main), and gravity outliers alone don't distinguish "gravity is
# high because the planet is big" from "gravity is off because density is
# off."
STATISTICAL_METRICS = ("surface_temperature", "atmospheric_pressure", "gravity", "scale_height", "density")


def distance_for_zone(star, zone):
    """
    Picks an orbital distance (AU) that lands in `zone` for `star`, mirroring
    tests/test_planets.py's `distance_for_zone` so generation here matches
    already-proven-valid usage of `Planet`.
    """
    inner, outer = star.habitable_zone
    if zone == 'h':
        return inner * 0.5
    if zone == 'c':
        return outer * 2.0
    return (inner + outer) / 2.0


def _corners(*ranges):
    return itertools.product(*ranges)


def theoretical_gravity_bounds_g(cls):
    """
    Analytically-derived (min, max) surface gravity in Earth g's for planet
    class `cls`, computed directly from the generator's own declared radius
    and density ranges (not hand-authored).

    Derivation: `calculate_surface_gravity` computes
    `g = G * mass / radius_m^2`, and `mass = volume(radius) * density`, so
    `g = (4/3) * pi * G * density_kg_m3 * radius_m` -- linear in radius and
    (for terrestrial classes) linear in density. For gas giants,
    `generate_planet_properties` blends the core and atmosphere densities via
    a mass-weighted harmonic mean: `1/density = ratio/rock_density +
    (1 - ratio) / atm_density` (see planetPhysics.py) -- physically correct
    for combining two densities via a mass fraction, but NOT multilinear
    (the ratio and atm_density terms invert `density`, not sum linearly).
    It is, however, monotonic in each of (rock_density, ratio, atm_density)
    individually holding the others fixed: `1/density` is a sum of terms
    each moving monotonically (one increasing, one decreasing) as `ratio`
    varies, and each of `rock_density`/`atm_density` only ever appears with
    a fixed-sign coefficient on its own reciprocal. A function monotonic in
    each variable over a box also has its extrema at the box's corners (same
    conclusion multilinearity would give, via a different property), so
    evaluating all corners and taking min/max still gives the exact
    theoretical range -- no approximation involved.

    Returns:
        tuple: (min_gravity_g, max_gravity_g).
    """
    data = prog_c.PLANET_CLASSES[cls]
    min_r_km, max_r_km = data["radius_range"]
    radii_m = (min_r_km * pc.KM_TO_M_FACTOR, max_r_km * pc.KM_TO_M_FACTOR)
    ptype = data["type"]

    def gravity_g(radius_m, density_kg_m3):
        g_ms2 = (4 / 3) * math.pi * pc.G * density_kg_m3 * radius_m
        return g_ms2 / pc.EARTH_GRAVITY

    # Class-specific density range if declared (e.g. a brown-dwarf-like
    # sub-stellar class -- see program_constants.PLANET_CLASSES), else the
    # default range shared by every other class of this body type. Must
    # mirror planetPhysics.py's own per-class density_range lookup exactly,
    # or this bound would flag legitimate override values as violations.
    values = []
    if ptype == "t":
        min_d, max_d = data.get("density_range", pc.PLANET_DENSITY["t"])  # g/cm^3
        densities_kgm3 = (min_d * 1000, max_d * 1000)
        for radius_m, density_kgm3 in _corners(radii_m, densities_kgm3):
            values.append(gravity_g(radius_m, density_kgm3))
    elif "density_range" in data:
        # A class with its own density_range skips the core/envelope blend
        # entirely in planetPhysics.py -- see
        # generate_planet_properties' matching `if planet.body_type == 'g'
        # and "density_range" not in class_data` guard -- so its declared
        # range alone bounds gravity here too.
        min_d, max_d = data["density_range"]  # g/cm^3
        densities_kgm3 = (min_d * 1000, max_d * 1000)
        for radius_m, density_kgm3 in _corners(radii_m, densities_kgm3):
            values.append(gravity_g(radius_m, density_kgm3))
    else:
        min_rock, max_rock = pc.PLANET_DENSITY["g"]  # g/cm^3
        min_ratio, max_ratio = prog_c.GAS_GIANT_CORE_ATMOSPHERE_RATIO
        # Must match planetPhysics.generate_planet_properties' envelope-side
        # blend input exactly: physical_constants.GAS_ENVELOPE_BULK_DENSITY
        # (already g/cm^3), not ATMOSPHERE_DENSITY["g"] (a different,
        # ~1000x-lighter physical layer -- see that constant's docstring).
        atm_gcm3_range = pc.GAS_ENVELOPE_BULK_DENSITY
        rock_range = (min_rock, max_rock)
        ratio_range = (min_ratio, max_ratio)
        for radius_m, rock_gcm3, ratio, atm_gcm3 in _corners(radii_m, rock_range, ratio_range, atm_gcm3_range):
            # Mass-weighted harmonic mean -- must match the blend formula in
            # planetPhysics.generate_planet_properties exactly.
            density_gcm3 = 1 / (ratio / rock_gcm3 + (1 - ratio) / atm_gcm3)
            density_kgm3 = density_gcm3 * 1000
            values.append(gravity_g(radius_m, density_kgm3))

    return min(values), max(values)


def generate_sample(cls, zone, n, include_moons=True):
    """
    Generates `n` independent planets of class `cls` in zone `zone`, each
    orbiting a freshly-generated host star whose spectral type is drawn
    uniformly from `HOST_STAR_TYPES`.

    Note: `stellarObjects.utils.reseed_rng()` reseeds the global `random`
    module from `secrets` at the start of most generation calls, so
    generation here (like everywhere else in this codebase) cannot be made
    deterministic via a seed -- successive runs of this tool will always
    see fresh random draws.

    Args:
        cls (str): Planet class code.
        zone (str): Zone character ('h', 'e', 'c'); must be valid for `cls`.
        n (int): Number of independent planets to generate.
        include_moons (bool): If True (default), also generates each
            planet's moons (via the normal random chance) and includes
            their metrics in the returned records, tagged `is_moon=True`.

    Returns:
        list[dict]: One record per generated body (planet, plus any moons).
    """
    records = []
    for _ in range(n):
        star_type = random.choice(HOST_STAR_TYPES)
        cfg = SystemConfig()
        cfg.STAR_TYPE = star_type
        star = Star(cfg)
        distance = distance_for_zone(star, zone)
        planet = Planet(
            cfg, star, star.habitable_zone, distance,
            planet_class=cls, moon_count=None if include_moons else 0,
        )
        records.append(_extract_record(planet, star_type, is_moon=False))
        for moon in planet.moons:
            records.append(_extract_record(moon, star_type, is_moon=True))
    return records


def _extract_record(body, star_type, is_moon):
    return {
        "planet_class": body.planet_class,
        "zone": body.zone,
        "star_type": star_type,
        "star_spectral_class": star_type[0],
        "is_moon": is_moon,
        "has_atmosphere": body.atmosphere != "None",
        "gravity": body.gravity,
        "surface_temperature": body.surface_temperature,
        "atmospheric_pressure": body.atmospheric_pressure,
        "scale_height": body.scale_height,
        "density": body.density,
        "mass": body.mass,
        "radius": body.radius,
    }


def check_hard_invariants(record):
    """
    Checks the always-true physical invariants for one record.

    Returns:
        list[str]: Human-readable descriptions of any violated invariant
            (empty if none).
    """
    issues = []

    gravity = record["gravity"]
    if not math.isfinite(gravity) or gravity <= 0:
        issues.append(f"gravity={gravity!r} is not a finite, positive value")
    else:
        lo, hi = theoretical_gravity_bounds_g(record["planet_class"])
        # 1% margin for floating-point slack in the corner-evaluation above,
        # not a physical fudge factor.
        margin = 0.01
        if not (lo * (1 - margin) <= gravity <= hi * (1 + margin)):
            issues.append(
                f"gravity={gravity:.4f}g outside theoretical range "
                f"[{lo:.4f}, {hi:.4f}]g analytically derived for class "
                f"{record['planet_class']!r}"
            )

    temperature = record["surface_temperature"]
    if not math.isfinite(temperature) or temperature <= 0:
        issues.append(f"surface_temperature={temperature!r} is not a finite, positive value (Kelvin)")

    pressure = record["atmospheric_pressure"]
    if not math.isfinite(pressure) or pressure < 0:
        issues.append(f"atmospheric_pressure={pressure!r} is not a finite, non-negative value (Pa)")

    if record["has_atmosphere"]:
        scale_height = record["scale_height"]
        if not math.isfinite(scale_height) or scale_height <= 0:
            issues.append(f"scale_height={scale_height!r} is not finite/positive despite having an atmosphere")
    else:
        if pressure != 0.0:
            issues.append(f"atmospheric_pressure={pressure!r} but class {record['planet_class']!r} has no atmosphere")

    return issues


def iqr_bounds(values, k=3.0):
    """
    Tukey's-fences bounds (Q1 - k*IQR, Q3 + k*IQR) for `values`.

    k=3.0 (the conventional "far outlier" threshold, vs. the more
    trigger-happy k=1.5 "mild outlier" one) is used by default since host
    stars are deliberately sampled across the whole main-sequence spectral
    range, and some real spread from that is expected, not buggy.

    Returns:
        tuple: (lower_bound, upper_bound). `(-inf, inf)` if there aren't
            enough values (< 4) to compute quartiles meaningfully.
    """
    if len(values) < 4:
        return (-math.inf, math.inf)
    q1, _, q3 = statistics.quantiles(values, n=4, method="inclusive")
    iqr = q3 - q1
    return (q1 - k * iqr, q3 + k * iqr)


def analyze(records, metrics=STATISTICAL_METRICS, k=3.0):
    """
    Groups `records` by planet class and computes, per metric: summary
    statistics and a list of statistical (IQR-based) outliers, plus the
    hard-invariant violations for every record regardless of grouping.

    Returns:
        dict: `{planet_class: {"n", "hard_violations", "metrics": {metric:
            {"n", "min", "max", "mean", "median", "stdev", "iqr_bounds",
            "outlier_count", "outlier_fraction", "outlier_examples"}}}}`.
    """
    by_class = defaultdict(list)
    for record in records:
        by_class[record["planet_class"]].append(record)

    report = {}
    for cls, recs in sorted(by_class.items()):
        hard_violations = []
        for record in recs:
            issues = check_hard_invariants(record)
            if issues:
                hard_violations.append((record, issues))

        metric_report = {}
        for metric in metrics:
            # scale_height is None for classes with no atmosphere (never set
            # by calculate_atmospheric_conditions) -- exclude those records
            # from this metric rather than treating None as an anomaly.
            values = [r[metric] for r in recs if r[metric] is not None and math.isfinite(r[metric])]
            if not values:
                continue
            lo, hi = iqr_bounds(values, k=k)
            outliers = [
                r for r in recs
                if r[metric] is not None and math.isfinite(r[metric]) and not (lo <= r[metric] <= hi)
            ]
            metric_report[metric] = {
                "n": len(values),
                "min": min(values),
                "max": max(values),
                "mean": statistics.mean(values),
                "median": statistics.median(values),
                "stdev": statistics.stdev(values) if len(values) > 1 else 0.0,
                "iqr_bounds": (lo, hi),
                "outlier_count": len(outliers),
                "outlier_fraction": len(outliers) / len(recs),
                "outlier_examples": outliers[:5],
            }

        report[cls] = {
            "n": len(recs),
            "hard_violations": hard_violations,
            "metrics": metric_report,
        }
    return report


def run_full_sample(n_per_class_zone=100, include_moons=True):
    """
    Generates `n_per_class_zone` planets for every valid (class, zone) pair
    `PLANET_CLASSES` declares, across the full `HOST_STAR_TYPES` grid.

    Returns:
        list[dict]: All generated records (planets and, if `include_moons`,
            their moons).
    """
    records = []
    for cls, zone in VALID_CLASS_ZONE_PAIRS:
        records.extend(generate_sample(cls, zone, n_per_class_zone, include_moons=include_moons))
    return records


def format_report(report, show_outlier_examples=True):
    """
    Renders `analyze`'s output as a human-readable multi-line report string.
    """
    lines = []
    total_hard_violations = sum(len(v["hard_violations"]) for v in report.values())
    total_n = sum(v["n"] for v in report.values())
    lines.append(f"Physical-plausibility report: {total_n} bodies across {len(report)} planet classes.")
    lines.append(f"Hard-invariant violations: {total_hard_violations}")
    lines.append("")

    for cls, data in report.items():
        lines.append(f"=== Class {cls} (n={data['n']}) ===")
        if data["hard_violations"]:
            lines.append(f"  HARD INVARIANT VIOLATIONS: {len(data['hard_violations'])}")
            for record, issues in data["hard_violations"][:5]:
                lines.append(
                    f"    - star={record['star_type']} zone={record['zone']} is_moon={record['is_moon']}: "
                    + "; ".join(issues)
                )
        for metric, stats in data["metrics"].items():
            flag = " <-- outliers present" if stats["outlier_count"] else ""
            lines.append(
                f"  {metric}: min={stats['min']:.4g} max={stats['max']:.4g} "
                f"mean={stats['mean']:.4g} median={stats['median']:.4g} "
                f"stdev={stats['stdev']:.4g} "
                f"outliers={stats['outlier_count']}/{data['n']} "
                f"({stats['outlier_fraction']:.1%}){flag}"
            )
            if show_outlier_examples and stats["outlier_examples"]:
                for record in stats["outlier_examples"]:
                    lines.append(
                        f"      e.g. {metric}={record[metric]:.4g} "
                        f"(star={record['star_type']}, zone={record['zone']}, is_moon={record['is_moon']})"
                    )
        lines.append("")

    return "\n".join(lines)
