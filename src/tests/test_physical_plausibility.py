"""
Physical-plausibility anomaly-finder regression tests.

See `planetgen/generation/plausibility.py` for the full two-tier design rationale,
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

from planetgen.generation import plausibility

# Small on purpose -- see module docstring. Large enough to exercise each
# zone's random draws a few times over; still fast (well under a second for
# the whole grid, including moons).
N_PER_ZONE = 30


@pytest.fixture(scope="module")
def full_sample():
    """One shared batch of records (planets + moons) for the whole module,
    covering every zone across the full host-star grid."""
    return plausibility.run_full_sample(n_per_zone=N_PER_ZONE, include_moons=True)


def test_zones_are_nonempty():
    assert plausibility.ZONE_CHARS


@pytest.mark.parametrize("body_type", ["t", "g"])
def test_theoretical_gravity_bounds_are_finite_and_ordered(body_type):
    lo, hi = plausibility.theoretical_gravity_bounds_g(body_type)
    assert math.isfinite(lo) and math.isfinite(hi)
    assert 0 < lo <= hi


@pytest.mark.parametrize("zone", plausibility.ZONE_CHARS)
def test_generated_sample_has_no_hard_invariant_violations(zone):
    records = plausibility.generate_sample(zone, N_PER_ZONE, include_moons=True)
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
        "body_type": "t",
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
    lo, hi = plausibility.theoretical_gravity_bounds_g("t")
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


def test_run_full_sample_covers_every_zone():
    records = plausibility.run_full_sample(n_per_zone=2, include_moons=False)
    assert {r["zone"] for r in records} == set(plausibility.ZONE_CHARS)
