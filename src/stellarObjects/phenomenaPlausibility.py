# stellarObjects/phenomenaPlausibility.py

"""
Exotic-phenomena plausibility anomaly finder.

The seven-phenomenon counterpart to `plausibility.py`'s planet-focused
engine (see that module's docstring for the full two-tier design
rationale, which this module follows exactly rather than inventing a new
pattern): `src/tests/test_phenomena_plausibility.py` gates the hard-
invariant half in CI, and `src/tests/phenomena_plausibility_cli.py` is the
CLI wrapper for a large, human-reviewed statistical batch.

Two-tier "sane range" design (mirrors `plausibility.py` exactly):
1. Hard physical invariants -- always true, analytically derived directly
   from each phenomenon's own declared `program_constants` ranges/formulas
   (e.g. a black hole's event horizon radius must match the Schwarzschild
   formula for its own generated mass; a neutron star's luminosity must
   match the Stefan-Boltzmann law for its own radius/temperature; a
   supernova remnant's radius must match the Sedov-Taylor relation for its
   own age). Any generated value violating one of these is unambiguously a
   bug. These gate the pytest test.
2. Statistical checks -- IQR-based (`plausibility.iqr_bounds`, reused
   directly rather than reimplemented) outlier detection on continuous
   metrics, plus category-frequency comparisons against each phenomenon's
   own configured chance (e.g. black hole accretion-disk proportion vs.
   `program_constants.BLACK_HOLE_ACCRETION_DISK_CHANCE`) -- reported for
   human review via `format_report`, never asserted to be exactly the
   configured value: this project has no chi-square/KS-test machinery
   anywhere, and a fixed-size sample's own real sampling noise means a
   observed proportion should land NEAR its configured chance, not
   exactly on it. Tolerance bands are simple, documented margins (the same
   empirically-calibrated-threshold convention
   `test_planet_physics_fixes.py`'s correlation checks already use), not a
   formal significance test.
"""

import math
import statistics
from collections import Counter, defaultdict

from . import physical_constants as pc
from . import program_constants as prog_c
from .asteroidFieldData import AsteroidField
from .compactRemnant import BlackHole, NeutronStar
from .config import SystemConfig
from .nebulaData import Nebula
from .plausibility import iqr_bounds
from .roguePlanetData import InterstellarComet, RoguePlanet
from .supernovaRemnantData import SupernovaRemnant

PHENOMENON_TYPES = prog_c.PHENOMENON_TYPE_CHOICES
"""tuple: The seven phenomenon types this module can generate/check --
identical to `phenomenonGen.py`'s own `--type` choices."""

STATISTICAL_METRICS_BY_TYPE = {
    "black-hole": ("mass_solar", "spin", "event_horizon_radius_km"),
    "neutron-star": ("mass_solar", "radius_km", "surface_temperature_k", "spin_period_ms", "magnetic_field_gauss"),
    "nebula": ("radius_ly",),
    "supernova-remnant": ("age_years", "radius_ly"),
    "rogue-planet": ("mass_kg", "radius_km"),
    "comet": ("nucleus_diameter_km", "velocity_kms"),
    "asteroid-field": ("radius_ly",),
}
"""dict: `phenomenon_type -> (metric_name, ...)` -- the continuous metrics
`analyze` runs IQR outlier detection on for that type, mirroring
`plausibility.STATISTICAL_METRICS`'s role but per-type rather than one
shared list, since these seven types share almost no metric names."""

# Rogue planets have no single configured "chance" for planet_type the way
# e.g. an accretion disk does -- it falls out of where
# ROGUE_PLANET_GAS_GIANT_MASS_THRESHOLD_JUPITER sits within the uniform
# ROGUE_PLANET_MASS_RANGE_JUPITER draw. Computed here (not hand-picked) so
# it stays correct if either constant ever changes.
_ROGUE_MASS_LO, _ROGUE_MASS_HI = prog_c.ROGUE_PLANET_MASS_RANGE_JUPITER
_ROGUE_GAS_GIANT_FRACTION = (
    (_ROGUE_MASS_HI - prog_c.ROGUE_PLANET_GAS_GIANT_MASS_THRESHOLD_JUPITER)
    / (_ROGUE_MASS_HI - _ROGUE_MASS_LO)
)

CATEGORY_EXPECTATIONS = {
    "black-hole": {
        "has_accretion_disk": {
            True: prog_c.BLACK_HOLE_ACCRETION_DISK_CHANCE,
            False: 1 - prog_c.BLACK_HOLE_ACCRETION_DISK_CHANCE,
        },
        "is_intermediate_mass": {
            True: prog_c.BLACK_HOLE_INTERMEDIATE_MASS_CHANCE,
            False: 1 - prog_c.BLACK_HOLE_INTERMEDIATE_MASS_CHANCE,
        },
    },
    "neutron-star": {
        "pulsar_type": {
            "millisecond": prog_c.NEUTRON_STAR_PULSAR_CHANCE * prog_c.PULSAR_MILLISECOND_CHANCE,
            "young": prog_c.NEUTRON_STAR_PULSAR_CHANCE * (1 - prog_c.PULSAR_MILLISECOND_CHANCE),
            "non-pulsing": 1 - prog_c.NEUTRON_STAR_PULSAR_CHANCE,
        },
    },
    "nebula": {
        "nebula_type": {t: 1 / len(prog_c.NEBULA_TYPES) for t in prog_c.NEBULA_TYPES},
    },
    "supernova-remnant": {
        "morphology": {m: 1 / len(prog_c.SUPERNOVA_REMNANT_MORPHOLOGIES) for m in prog_c.SUPERNOVA_REMNANT_MORPHOLOGIES},
        "progenitor_type": {
            "Type Ia": prog_c.SUPERNOVA_PROGENITOR_TYPE_IA_CHANCE,
            "core-collapse": 1 - prog_c.SUPERNOVA_PROGENITOR_TYPE_IA_CHANCE,
        },
    },
    "rogue-planet": {
        "planet_type": {
            "g": _ROGUE_GAS_GIANT_FRACTION,
            "t": 1 - _ROGUE_GAS_GIANT_FRACTION,
        },
    },
    "comet": {
        "is_active": {
            True: prog_c.INTERSTELLAR_COMET_ACTIVE_CHANCE,
            False: 1 - prog_c.INTERSTELLAR_COMET_ACTIVE_CHANCE,
        },
    },
    "asteroid-field": {
        "density": {"dense": 1 / 3, "sparse": 1 / 3, "typical": 1 / 3},
    },
}
"""dict: `phenomenon_type -> {field: {value: expected_proportion, ...}}`,
each inner dict's proportions summing to 1 -- the discrete-choice
counterpart to `STATISTICAL_METRICS_BY_TYPE`'s continuous metrics.
Reported by `format_report`, never asserted exactly (see module
docstring)."""


def _build_phenomenon(phenomenon_type, system_config):
    """
    Instantiates one phenomenon of `phenomenon_type`, always standalone
    (never anchored to a `StarSystem` -- that integration path is already
    covered by `test_phenomena.py`'s own `StarSystem(compact_remnant=...)`
    tests, not this module's concern).

    Args:
        phenomenon_type (str): One of `PHENOMENON_TYPES`.
        system_config (SystemConfig): The config to generate with.

    Returns:
        The generated phenomenon object.

    Raises:
        ValueError: If `phenomenon_type` isn't recognized.
    """
    if phenomenon_type == "black-hole":
        return BlackHole(system_config)
    if phenomenon_type == "neutron-star":
        return NeutronStar(system_config)
    if phenomenon_type == "nebula":
        return Nebula(system_config)
    if phenomenon_type == "supernova-remnant":
        return SupernovaRemnant(system_config)
    if phenomenon_type == "rogue-planet":
        return RoguePlanet(system_config)
    if phenomenon_type == "comet":
        return InterstellarComet(system_config)
    if phenomenon_type == "asteroid-field":
        return AsteroidField(system_config)
    raise ValueError(f"Unknown phenomenon type: {phenomenon_type!r}")


def _extract_record(phenomenon_type, obj):
    """
    Flattens one generated phenomenon object into a plain dict record --
    the shape every other function in this module (`check_hard_invariants`,
    `analyze`) operates on, mirroring `plausibility._extract_record`'s
    identical role.
    """
    record = {
        "phenomenon_type": phenomenon_type,
        "galactic_orbital_speed_kms": obj.galactic_orbital_speed_kms,
        "galactic_orbital_period_gy": obj.galactic_orbital_period_gy,
        "galactic_orbital_phase_deg": obj.galactic_orbital_phase_deg,
        "galactic_min_update_interval_years": obj.galactic_min_update_interval_years,
    }

    if phenomenon_type == "black-hole":
        record.update({
            "mass_solar": obj.mass_solar,
            "spin": obj.spin,
            "event_horizon_radius_km": obj.event_horizon_radius_km,
            "has_accretion_disk": obj.has_accretion_disk,
            "luminosity_w": obj.luminosity,
            "is_intermediate_mass": obj.mass_solar >= prog_c.BLACK_HOLE_INTERMEDIATE_MASS_RANGE_SOLAR[0],
        })
    elif phenomenon_type == "neutron-star":
        record.update({
            "mass_solar": obj.mass_solar,
            "radius_km": obj.radius,
            "spin_period_ms": obj.spin_period_ms,
            "magnetic_field_gauss": obj.magnetic_field_gauss,
            "pulsar_type": obj.pulsar_type,
            "surface_temperature_k": obj.surface_temperature_k,
            "luminosity_w": obj.luminosity,
        })
    elif phenomenon_type == "nebula":
        record.update({"nebula_type": obj.nebula_type, "radius_ly": obj.radius_ly})
    elif phenomenon_type == "supernova-remnant":
        record.update({
            "morphology": obj.morphology,
            "age_years": obj.age_years,
            "radius_ly": obj.radius_ly,
            "progenitor_type": obj.progenitor_type,
            "has_compact_remnant": obj.compact_remnant is not None,
        })
    elif phenomenon_type == "rogue-planet":
        record.update({
            "mass_kg": obj.mass_kg,
            "radius_km": obj.radius_km,
            "planet_type": obj.planet_type,
        })
    elif phenomenon_type == "comet":
        record.update({
            "nucleus_diameter_km": obj.nucleus_diameter_km,
            "velocity_kms": obj.velocity_kms,
            "is_active": obj.is_active,
        })
    elif phenomenon_type == "asteroid-field":
        record.update({"density": obj.density, "radius_ly": obj.radius_ly})

    return record


def generate_sample(phenomenon_type, n):
    """
    Generates `n` independent standalone phenomena of `phenomenon_type`.

    Args:
        phenomenon_type (str): One of `PHENOMENON_TYPES`.
        n (int): Number of independent bodies to generate.

    Returns:
        list[dict]: One record per generated body.
    """
    records = []
    for _ in range(n):
        cfg = SystemConfig()
        obj = _build_phenomenon(phenomenon_type, cfg)
        records.append(_extract_record(phenomenon_type, obj))
    return records


def run_full_sample(n_per_type=200):
    """
    Generates `n_per_type` bodies for every phenomenon type.

    Args:
        n_per_type (int): Number of independent bodies to generate per type.

    Returns:
        list[dict]: All generated records, across every type.
    """
    records = []
    for phenomenon_type in PHENOMENON_TYPES:
        records.extend(generate_sample(phenomenon_type, n_per_type))
    return records


def _isclose_or_flag(issues, label, actual, expected, rel_tol=1e-6):
    """Shared helper: appends an issue to `issues` if `actual` isn't
    `math.isclose` to `expected` -- used by every "must match a derived
    formula exactly" check below, so each one reads as a single line."""
    if not math.isclose(actual, expected, rel_tol=rel_tol):
        issues.append(f"{label}={actual!r} does not match the derived value {expected!r}")


def check_hard_invariants(record):
    """
    Checks the always-true physical invariants for one record, both the
    universal galactic-orbit fields every phenomenon type shares (see
    `nebulaData.Nebula`'s `galactic_orbital_*` docstring) and the
    type-specific ones.

    Returns:
        list[str]: Human-readable descriptions of any violated invariant
            (empty if none).
    """
    issues = []
    phenomenon_type = record["phenomenon_type"]

    phase = record["galactic_orbital_phase_deg"]
    if not (math.isfinite(phase) and 0 <= phase < 360):
        issues.append(f"galactic_orbital_phase_deg={phase!r} is not finite and in [0, 360)")
    period = record["galactic_orbital_period_gy"]
    if not (math.isfinite(period) and period > 0):
        issues.append(f"galactic_orbital_period_gy={period!r} is not a finite, positive value")
    speed = record["galactic_orbital_speed_kms"]
    if not (math.isfinite(speed) and speed > 0):
        issues.append(f"galactic_orbital_speed_kms={speed!r} is not a finite, positive value")
    guard = record["galactic_min_update_interval_years"]
    if not (math.isfinite(guard) and guard > 0):
        issues.append(f"galactic_min_update_interval_years={guard!r} is not a finite, positive value")

    if phenomenon_type == "black-hole":
        mass = record["mass_solar"]
        stellar_lo, stellar_hi = prog_c.BLACK_HOLE_MASS_RANGE_SOLAR
        inter_lo, inter_hi = prog_c.BLACK_HOLE_INTERMEDIATE_MASS_RANGE_SOLAR
        if not (stellar_lo <= mass <= stellar_hi or inter_lo <= mass <= inter_hi):
            issues.append(f"mass_solar={mass!r} outside both stellar-mass and intermediate-mass ranges")
        spin_lo, spin_hi = prog_c.BLACK_HOLE_SPIN_RANGE
        if not (spin_lo <= record["spin"] <= spin_hi):
            issues.append(f"spin={record['spin']!r} outside {prog_c.BLACK_HOLE_SPIN_RANGE}")
        expected_radius_km = (
            2 * pc.G * (mass * pc.SOLAR_MASS_TO_KG) / pc.SPEED_OF_LIGHT_M_S ** 2
        ) / 1000
        _isclose_or_flag(issues, "event_horizon_radius_km", record["event_horizon_radius_km"], expected_radius_km)
        if not record["has_accretion_disk"] and record["luminosity_w"] != 0.0:
            issues.append("luminosity_w is nonzero despite has_accretion_disk being False")
        if record["has_accretion_disk"] and record["luminosity_w"] <= 0.0:
            issues.append("luminosity_w is not positive despite has_accretion_disk being True")

    elif phenomenon_type == "neutron-star":
        mass_lo, mass_hi = prog_c.NEUTRON_STAR_MASS_RANGE_SOLAR
        if not (mass_lo <= record["mass_solar"] <= mass_hi):
            issues.append(f"mass_solar={record['mass_solar']!r} outside {prog_c.NEUTRON_STAR_MASS_RANGE_SOLAR}")
        radius_lo, radius_hi = prog_c.NEUTRON_STAR_RADIUS_RANGE_KM
        if not (radius_lo <= record["radius_km"] <= radius_hi):
            issues.append(f"radius_km={record['radius_km']!r} outside {prog_c.NEUTRON_STAR_RADIUS_RANGE_KM}")
        temp_lo, temp_hi = prog_c.NEUTRON_STAR_SURFACE_TEMPERATURE_RANGE_K
        if not (temp_lo <= record["surface_temperature_k"] <= temp_hi):
            issues.append(
                f"surface_temperature_k={record['surface_temperature_k']!r} outside "
                f"{prog_c.NEUTRON_STAR_SURFACE_TEMPERATURE_RANGE_K}"
            )
        if record["pulsar_type"] not in ("young", "millisecond", "non-pulsing"):
            issues.append(f"pulsar_type={record['pulsar_type']!r} is not a recognized value")
        radius_m = record["radius_km"] * pc.KM_TO_M_FACTOR
        expected_luminosity_w = (
            pc.STEFAN_BOLTZMANN_CONSTANT * 4 * math.pi * radius_m ** 2
            * record["surface_temperature_k"] ** 4
        )
        _isclose_or_flag(issues, "luminosity_w", record["luminosity_w"], expected_luminosity_w)

    elif phenomenon_type == "nebula":
        type_data = prog_c.NEBULA_TYPES.get(record["nebula_type"])
        if type_data is None:
            issues.append(f"nebula_type={record['nebula_type']!r} is not a recognized value")
        else:
            lo, hi = type_data["radius_range_ly"]
            if not (lo <= record["radius_ly"] <= hi):
                issues.append(f"radius_ly={record['radius_ly']!r} outside {record['nebula_type']}'s own {(lo, hi)}")

    elif phenomenon_type == "supernova-remnant":
        if record["morphology"] not in prog_c.SUPERNOVA_REMNANT_MORPHOLOGIES:
            issues.append(f"morphology={record['morphology']!r} is not a recognized value")
        age_lo, age_hi = prog_c.SUPERNOVA_REMNANT_AGE_RANGE_YEARS
        if not (age_lo <= record["age_years"] <= age_hi):
            issues.append(f"age_years={record['age_years']!r} outside {prog_c.SUPERNOVA_REMNANT_AGE_RANGE_YEARS}")
        expected_radius_ly = (
            prog_c.SEDOV_TAYLOR_RADIUS_COEFFICIENT_LY
            * (record["age_years"] ** prog_c.SEDOV_TAYLOR_TIME_EXPONENT)
        )
        _isclose_or_flag(issues, "radius_ly", record["radius_ly"], expected_radius_ly)
        if record["progenitor_type"] not in ("Type Ia", "core-collapse"):
            issues.append(f"progenitor_type={record['progenitor_type']!r} is not a recognized value")
        if record["progenitor_type"] == "Type Ia" and record["has_compact_remnant"]:
            issues.append("a Type Ia progenitor should never leave a compact remnant behind")

    elif phenomenon_type == "rogue-planet":
        mass_jupiter = record["mass_kg"] / pc.JUPITER_MASS_TO_KG
        lo, hi = prog_c.ROGUE_PLANET_MASS_RANGE_JUPITER
        if not (lo <= mass_jupiter <= hi):
            issues.append(f"mass (in Jupiter masses)={mass_jupiter!r} outside {prog_c.ROGUE_PLANET_MASS_RANGE_JUPITER}")
        if record["planet_type"] not in ('t', 'g'):
            issues.append(f"planet_type={record['planet_type']!r} is not 't' or 'g'")
        elif record["planet_type"] == 'g' and mass_jupiter < prog_c.ROGUE_PLANET_GAS_GIANT_MASS_THRESHOLD_JUPITER:
            issues.append("planet_type is 'g' but mass is below the gas-giant threshold")
        elif record["planet_type"] == 't' and mass_jupiter >= prog_c.ROGUE_PLANET_GAS_GIANT_MASS_THRESHOLD_JUPITER:
            issues.append("planet_type is 't' but mass is at/above the gas-giant threshold")
        if not (math.isfinite(record["radius_km"]) and record["radius_km"] > 0):
            issues.append(f"radius_km={record['radius_km']!r} is not a finite, positive value")

    elif phenomenon_type == "comet":
        lo, hi = prog_c.INTERSTELLAR_COMET_NUCLEUS_DIAMETER_RANGE_KM
        if not (lo <= record["nucleus_diameter_km"] <= hi):
            issues.append(
                f"nucleus_diameter_km={record['nucleus_diameter_km']!r} outside "
                f"{prog_c.INTERSTELLAR_COMET_NUCLEUS_DIAMETER_RANGE_KM}"
            )
        lo, hi = prog_c.INTERSTELLAR_OBJECT_SPEED_KMS_RANGE
        if not (lo <= record["velocity_kms"] <= hi):
            issues.append(f"velocity_kms={record['velocity_kms']!r} outside {prog_c.INTERSTELLAR_OBJECT_SPEED_KMS_RANGE}")

    elif phenomenon_type == "asteroid-field":
        lo, hi = prog_c.ASTEROID_FIELD_RADIUS_RANGE_LY
        if not (lo <= record["radius_ly"] <= hi):
            issues.append(f"radius_ly={record['radius_ly']!r} outside {prog_c.ASTEROID_FIELD_RADIUS_RANGE_LY}")
        if record["density"] not in ("dense", "sparse", "typical"):
            issues.append(f"density={record['density']!r} is not a recognized value")

    return issues


def analyze(records, k=3.0):
    """
    Groups `records` by `phenomenon_type` and computes, per type: summary
    statistics and IQR-based outliers for that type's own continuous
    metrics (`STATISTICAL_METRICS_BY_TYPE`), observed-vs-expected category
    frequencies for its own discrete choices (`CATEGORY_EXPECTATIONS`), and
    hard-invariant violations for every record regardless of grouping.
    Mirrors `plausibility.analyze`'s shape/role.

    Returns:
        dict: `{phenomenon_type: {"n", "hard_violations", "metrics":
            {metric: {...same shape as plausibility.analyze's per-metric
            dict...}}, "categories": {field: {"counts", "proportions",
            "expected"}}}}`.
    """
    by_type = defaultdict(list)
    for record in records:
        by_type[record["phenomenon_type"]].append(record)

    report = {}
    for phenomenon_type, recs in sorted(by_type.items()):
        hard_violations = []
        for record in recs:
            issues = check_hard_invariants(record)
            if issues:
                hard_violations.append((record, issues))

        metric_report = {}
        for metric in STATISTICAL_METRICS_BY_TYPE.get(phenomenon_type, ()):
            values = [r[metric] for r in recs if r.get(metric) is not None and math.isfinite(r[metric])]
            if not values:
                continue
            lo, hi = iqr_bounds(values, k=k)
            outliers = [
                r for r in recs
                if r.get(metric) is not None and math.isfinite(r[metric]) and not (lo <= r[metric] <= hi)
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

        category_report = {}
        for field, expected_proportions in CATEGORY_EXPECTATIONS.get(phenomenon_type, {}).items():
            counts = Counter(r[field] for r in recs if field in r)
            n_categorized = sum(counts.values())
            category_report[field] = {
                "counts": dict(counts),
                "proportions": {value: count / n_categorized for value, count in counts.items()} if n_categorized else {},
                "expected": dict(expected_proportions),
            }

        report[phenomenon_type] = {
            "n": len(recs),
            "hard_violations": hard_violations,
            "metrics": metric_report,
            "categories": category_report,
        }
    return report


def format_report(report, show_outlier_examples=True, category_flag_margin=0.1):
    """
    Renders `analyze`'s output as a human-readable multi-line report
    string, mirroring `plausibility.format_report`'s shape, plus a
    category-frequency section per type.

    Args:
        report (dict): `analyze`'s return value.
        show_outlier_examples (bool): Whether to print individual outlier
            example lines under each metric.
        category_flag_margin (float): A category's observed proportion is
            flagged (not failed) when it differs from its expected
            proportion by more than this -- a simple, documented margin
            (see module docstring), not a formal significance test.

    Returns:
        str: The formatted report.
    """
    lines = []
    total_hard_violations = sum(len(v["hard_violations"]) for v in report.values())
    total_n = sum(v["n"] for v in report.values())
    lines.append(f"Exotic-phenomena plausibility report: {total_n} bodies across {len(report)} phenomenon types.")
    lines.append(f"Hard-invariant violations: {total_hard_violations}")
    lines.append("")

    for phenomenon_type, data in report.items():
        lines.append(f"=== {phenomenon_type} (n={data['n']}) ===")
        if data["hard_violations"]:
            lines.append(f"  HARD INVARIANT VIOLATIONS: {len(data['hard_violations'])}")
            for record, issues in data["hard_violations"][:5]:
                lines.append(f"    - {record}: " + "; ".join(issues))

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
                    lines.append(f"      e.g. {metric}={record[metric]:.4g}")

        for field, cat_data in data["categories"].items():
            lines.append(f"  {field} (observed vs. expected):")
            for value, expected_proportion in sorted(cat_data["expected"].items(), key=lambda item: str(item[0])):
                observed_proportion = cat_data["proportions"].get(value, 0.0)
                flag = (
                    " <-- deviates from expected"
                    if abs(observed_proportion - expected_proportion) > category_flag_margin else ""
                )
                lines.append(
                    f"    {value!r}: observed {observed_proportion:.1%} vs. expected "
                    f"{expected_proportion:.1%}{flag}"
                )
        lines.append("")

    return "\n".join(lines)
