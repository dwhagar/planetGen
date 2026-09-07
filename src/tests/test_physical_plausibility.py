"""
Physical-plausibility anomaly-finder regression tests.

See `stellarObjects/plausibility.py` for the full two-tier design rationale,
and TODO.md's "Physical-plausibility test suite (anomaly finder)" future
idea for the original ask. This module only gates the *hard-invariant*
half of that design (unambiguous bugs -- see `check_hard_invariants` and
`theoretical_gravity_bounds_g`): a value outside an analytically-derived
bound, or a non-finite/wrong-sign temperature or pressure, is always a bug
and should fail CI.

The *statistical* half (Tukey's-fences outlier detection across a large,
host-star-spectral-type-spanning batch) is deliberately NOT asserted here --
some real, non-buggy spread across that batch is expected, and asserting an
outlier count of zero would just reintroduce the kind of band-aid the
now-disabled Class M/P clamps used to be. That half is meant for a human to
read via the `physical_plausibility_cli.py` CLI script, not to gate a test run.
This is the module's answer to the TODO's open question ("whether this lives
in tests/ as a slow/opt-in suite or as a separate standalone script"): both,
split by tier. The generation batch below is deliberately small (a handful
of bodies per (class, zone) pair, not the CLI's default 150) precisely
because only the fast, deterministic hard-invariant checks run here --
keeping this a normal (not opt-in/slow) test, consistent with how
test_full_matrix.py already runs its own large combinatorial sweep in the
regular suite rather than behind a marker.

Run with: pytest tests/test_physical_plausibility.py
"""
import math

import pytest

from stellarObjects import plausibility
from stellarObjects import program_constants as prog_c

# Small on purpose -- see module docstring. Large enough to exercise each
# (class, zone) pair's random draws a few times over; still fast (well
# under a second for the whole grid, including moons).
N_PER_PAIR = 15


@pytest.fixture(scope="module")
def full_sample():
    """One shared batch of records (planets + moons) for the whole module,
    covering every valid (class, zone) pair across the full host-star grid."""
    return plausibility.run_full_sample(n_per_class_zone=N_PER_PAIR, include_moons=True)


def test_valid_class_zone_pairs_is_nonempty():
    # Sanity check on the derivation itself: if PLANET_CLASSES ever loses
    # every zone flag for every class, the rest of this module would
    # silently pass on zero bodies without this.
    assert plausibility.VALID_CLASS_ZONE_PAIRS
    for cls, zone in plausibility.VALID_CLASS_ZONE_PAIRS:
        assert prog_c.PLANET_CLASSES[cls][zone]


@pytest.mark.parametrize("cls", sorted(prog_c.PLANET_CLASSES))
def test_theoretical_gravity_bounds_are_finite_and_ordered(cls):
    # Corner-evaluation should never itself blow up or invert, regardless
    # of whether the class has a valid zone (a class with no valid zone
    # still has radius/density ranges declared).
    lo, hi = plausibility.theoretical_gravity_bounds_g(cls)
    assert math.isfinite(lo) and math.isfinite(hi)
    assert 0 < lo <= hi


@pytest.mark.parametrize("cls,zone", plausibility.VALID_CLASS_ZONE_PAIRS,
                          ids=[f"{c}-{z}" for c, z in plausibility.VALID_CLASS_ZONE_PAIRS])
def test_generated_sample_has_no_hard_invariant_violations(cls, zone):
    records = plausibility.generate_sample(cls, zone, N_PER_PAIR, include_moons=True)
    assert records  # generation actually produced bodies (plus possibly moons)
    violations = []
    for record in records:
        issues = plausibility.check_hard_invariants(record)
        if issues:
            violations.append((record["planet_class"], record["is_moon"], record["star_type"], issues))
    assert not violations, f"Hard-invariant violations found: {violations}"


def test_analyze_report_structure(full_sample):
    report = plausibility.analyze(full_sample)
    assert set(report) == {r["planet_class"] for r in full_sample}
    for cls, data in report.items():
        assert data["n"] > 0
        assert data["hard_violations"] == []
        for metric, stats in data["metrics"].items():
            assert metric in plausibility.STATISTICAL_METRICS
            assert stats["n"] > 0
            assert stats["min"] <= stats["median"] <= stats["max"]
            assert stats["iqr_bounds"][0] <= stats["iqr_bounds"][1]
            assert 0 <= stats["outlier_fraction"] <= 1


def test_format_report_runs_without_error(full_sample):
    report = plausibility.analyze(full_sample)
    text = plausibility.format_report(report)
    assert "Physical-plausibility report:" in text
    assert "Hard-invariant violations: 0" in text


# --- Unit tests on the checking logic itself, independent of generation ---

def _base_record(**overrides):
    record = {
        "planet_class": "M",
        "zone": "e",
        "star_type": "G5V",
        "star_spectral_class": "G",
        "is_moon": False,
        "has_atmosphere": True,
        "gravity": 1.0,
        "surface_temperature": 288.0,
        "atmospheric_pressure": 101325.0,
        "scale_height": 8.5,
        "density": 5.0,
        "mass": 5.9e24,
        "radius": 6371.0,
    }
    record.update(overrides)
    return record


def test_check_hard_invariants_accepts_plausible_record():
    assert plausibility.check_hard_invariants(_base_record()) == []


@pytest.mark.parametrize("bad_gravity", [0.0, -1.0, math.inf, math.nan])
def test_check_hard_invariants_flags_bad_gravity(bad_gravity):
    issues = plausibility.check_hard_invariants(_base_record(gravity=bad_gravity))
    assert any("gravity" in issue for issue in issues)


def test_check_hard_invariants_flags_gravity_outside_theoretical_bounds():
    lo, hi = plausibility.theoretical_gravity_bounds_g("M")
    issues = plausibility.check_hard_invariants(_base_record(gravity=hi * 10))
    assert any("outside theoretical range" in issue for issue in issues)


@pytest.mark.parametrize("bad_temperature", [0.0, -10.0, math.inf, math.nan])
def test_check_hard_invariants_flags_bad_temperature(bad_temperature):
    issues = plausibility.check_hard_invariants(_base_record(surface_temperature=bad_temperature))
    assert any("surface_temperature" in issue for issue in issues)


@pytest.mark.parametrize("bad_pressure", [-1.0, math.inf, math.nan])
def test_check_hard_invariants_flags_bad_pressure(bad_pressure):
    issues = plausibility.check_hard_invariants(_base_record(atmospheric_pressure=bad_pressure))
    assert any("atmospheric_pressure" in issue for issue in issues)


def test_check_hard_invariants_flags_pressure_without_atmosphere():
    issues = plausibility.check_hard_invariants(
        _base_record(has_atmosphere=False, atmospheric_pressure=1000.0)
    )
    assert any("no atmosphere" in issue for issue in issues)


def test_check_hard_invariants_flags_missing_scale_height_with_atmosphere():
    issues = plausibility.check_hard_invariants(
        _base_record(has_atmosphere=True, scale_height=0.0)
    )
    assert any("scale_height" in issue for issue in issues)


def test_iqr_bounds_too_few_values_is_unbounded():
    assert plausibility.iqr_bounds([1.0, 2.0, 3.0]) == (-math.inf, math.inf)


def test_iqr_bounds_flags_a_far_outlier():
    values = [10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 1000.0]
    lo, hi = plausibility.iqr_bounds(values, k=3.0)
    assert lo <= 10.0 and 15.0 <= hi < 1000.0


def test_run_full_sample_covers_every_valid_pair():
    records = plausibility.run_full_sample(n_per_class_zone=2, include_moons=False)
    seen_classes = {r["planet_class"] for r in records}
    expected_classes = {cls for cls, _ in plausibility.VALID_CLASS_ZONE_PAIRS}
    assert seen_classes == expected_classes
