"""
Exotic-phenomena plausibility anomaly-finder regression tests.

See `stellarObjects/phenomenaPlausibility.py` for the full two-tier design
rationale -- the seven-phenomenon counterpart to `plausibility.py`/
`test_physical_plausibility.py`. This module only gates the *hard-
invariant* half (unambiguous bugs -- see `check_hard_invariants`): a value
outside an analytically-derived bound or formula (e.g. an event horizon
radius that doesn't match the Schwarzschild formula for its own mass) is
always a bug and should fail CI.

The *statistical* half (IQR outlier detection, category-frequency
comparisons against each phenomenon's own configured chance) is
deliberately NOT asserted here -- some real spread/sampling noise across a
batch is expected, the same reasoning `test_physical_plausibility.py`'s own
module docstring gives. That half is meant for a human to read via the
`phenomena_plausibility_cli.py` CLI script, not to gate a test run.

The generation batch below is deliberately small (a handful of bodies per
phenomenon type, not the CLI's default 200) precisely because only the
fast, deterministic hard-invariant checks run here -- keeping this a
normal (not opt-in/slow) test.

Run with: pytest src/tests/test_phenomena_plausibility.py
"""
import math

import pytest

from stellarObjects import phenomenaPlausibility as pp
from stellarObjects import physical_constants as pc
from stellarObjects import program_constants as prog_c

# Small on purpose -- see module docstring. Large enough to exercise each
# type's random draws (including its rarer branches, e.g. an accretion
# disk or an intermediate-mass black hole) a few times over; still fast.
N_PER_TYPE = 30


@pytest.fixture(scope="module")
def full_sample():
    """One shared batch of records for the whole module, covering every
    phenomenon type."""
    return pp.run_full_sample(n_per_type=N_PER_TYPE)


def test_phenomenon_types_is_nonempty():
    assert pp.PHENOMENON_TYPES
    assert len(pp.PHENOMENON_TYPES) == 7


@pytest.mark.parametrize("phenomenon_type", pp.PHENOMENON_TYPES)
def test_generated_sample_has_no_hard_invariant_violations(phenomenon_type):
    records = pp.generate_sample(phenomenon_type, N_PER_TYPE)
    assert records  # generation actually produced bodies
    violations = []
    for record in records:
        issues = pp.check_hard_invariants(record)
        if issues:
            violations.append((phenomenon_type, issues))
    assert not violations, f"Hard-invariant violations found: {violations}"


def test_analyze_report_structure(full_sample):
    report = pp.analyze(full_sample)
    assert set(report) == {r["phenomenon_type"] for r in full_sample}
    for phenomenon_type, data in report.items():
        assert data["n"] > 0
        assert data["hard_violations"] == []
        for metric, stats in data["metrics"].items():
            assert metric in pp.STATISTICAL_METRICS_BY_TYPE[phenomenon_type]
            assert stats["n"] > 0
            assert stats["min"] <= stats["median"] <= stats["max"]
            assert stats["iqr_bounds"][0] <= stats["iqr_bounds"][1]
            assert 0 <= stats["outlier_fraction"] <= 1
        for field, cat_data in data["categories"].items():
            assert field in pp.CATEGORY_EXPECTATIONS[phenomenon_type]
            assert math.isclose(sum(cat_data["expected"].values()), 1.0, rel_tol=1e-9)
            assert math.isclose(sum(cat_data["proportions"].values()), 1.0, rel_tol=1e-9)


def test_format_report_runs_without_error(full_sample):
    report = pp.analyze(full_sample)
    text = pp.format_report(report)
    assert "Exotic-phenomena plausibility report:" in text
    assert "Hard-invariant violations: 0" in text


def test_run_full_sample_covers_every_type():
    records = pp.run_full_sample(n_per_type=2)
    seen_types = {r["phenomenon_type"] for r in records}
    assert seen_types == set(pp.PHENOMENON_TYPES)


# --- Unit tests on the checking logic itself, independent of generation ---

_BLACK_HOLE_MASS_SOLAR = 10.0
_BLACK_HOLE_EVENT_HORIZON_KM = (
    2 * pc.G * (_BLACK_HOLE_MASS_SOLAR * pc.SOLAR_MASS_TO_KG) / pc.SPEED_OF_LIGHT_M_S ** 2
) / 1000

_NEUTRON_STAR_RADIUS_KM = 11.0
_NEUTRON_STAR_SURFACE_TEMPERATURE_K = 5e5
_NEUTRON_STAR_LUMINOSITY_W = (
    pc.STEFAN_BOLTZMANN_CONSTANT * 4 * math.pi * (_NEUTRON_STAR_RADIUS_KM * pc.KM_TO_M_FACTOR) ** 2
    * _NEUTRON_STAR_SURFACE_TEMPERATURE_K ** 4
)

_SUPERNOVA_AGE_YEARS = 1000.0


def _base_record(phenomenon_type, **overrides):
    common = {
        "phenomenon_type": phenomenon_type,
        "galactic_orbital_speed_kms": 205.7,
        "galactic_orbital_period_gy": 0.236,
        "galactic_orbital_phase_deg": 180.0,
        "galactic_min_update_interval_years": 1e-9,
    }
    by_type = {
        "black-hole": {
            "mass_solar": _BLACK_HOLE_MASS_SOLAR, "spin": 0.5,
            "event_horizon_radius_km": _BLACK_HOLE_EVENT_HORIZON_KM,
            "has_accretion_disk": False, "luminosity_w": 0.0, "is_intermediate_mass": False,
        },
        "neutron-star": {
            "mass_solar": 1.4, "radius_km": _NEUTRON_STAR_RADIUS_KM, "spin_period_ms": 100.0,
            "magnetic_field_gauss": 1e11, "pulsar_type": "young",
            "surface_temperature_k": _NEUTRON_STAR_SURFACE_TEMPERATURE_K,
            "luminosity_w": _NEUTRON_STAR_LUMINOSITY_W,
        },
        "nebula": {"nebula_type": "emission", "radius_ly": 50.0},
        "supernova-remnant": {
            "morphology": "shell", "age_years": _SUPERNOVA_AGE_YEARS,
            "radius_ly": (
                prog_c.SEDOV_TAYLOR_RADIUS_COEFFICIENT_LY
                * (_SUPERNOVA_AGE_YEARS ** prog_c.SEDOV_TAYLOR_TIME_EXPONENT)
            ),
            "progenitor_type": "core-collapse", "has_compact_remnant": True,
        },
        "rogue-planet": {"mass_kg": 1e27, "radius_km": 70000.0, "planet_type": "g"},
        "comet": {"nucleus_diameter_km": 1.0, "velocity_kms": 40.0, "is_active": True},
        "asteroid-field": {"density": "typical", "radius_ly": 0.5},
    }
    record = {**common, **by_type[phenomenon_type]}
    record.update(overrides)
    return record


@pytest.mark.parametrize("phenomenon_type", pp.PHENOMENON_TYPES)
def test_check_hard_invariants_accepts_plausible_record(phenomenon_type):
    assert pp.check_hard_invariants(_base_record(phenomenon_type)) == []


@pytest.mark.parametrize("bad_phase", [-1.0, 360.0, math.inf, math.nan])
def test_check_hard_invariants_flags_bad_galactic_phase(bad_phase):
    issues = pp.check_hard_invariants(_base_record("nebula", galactic_orbital_phase_deg=bad_phase))
    assert any("galactic_orbital_phase_deg" in issue for issue in issues)


@pytest.mark.parametrize("bad_period", [0.0, -1.0, math.inf, math.nan])
def test_check_hard_invariants_flags_bad_galactic_period(bad_period):
    issues = pp.check_hard_invariants(_base_record("nebula", galactic_orbital_period_gy=bad_period))
    assert any("galactic_orbital_period_gy" in issue for issue in issues)


def test_check_hard_invariants_flags_black_hole_event_horizon_mismatch():
    issues = pp.check_hard_invariants(_base_record("black-hole", event_horizon_radius_km=999.0))
    assert any("event_horizon_radius_km" in issue for issue in issues)


def test_check_hard_invariants_flags_black_hole_mass_out_of_range():
    # 50 solar masses falls in the gap between the stellar-mass range
    # (5-20) and the intermediate-mass range (100-1000) -- its own
    # Schwarzschild radius is computed so only the mass-range check (not
    # also a radius mismatch) fires.
    mass_solar = 50.0
    radius_km = (2 * pc.G * (mass_solar * pc.SOLAR_MASS_TO_KG) / pc.SPEED_OF_LIGHT_M_S ** 2) / 1000
    issues = pp.check_hard_invariants(
        _base_record("black-hole", mass_solar=mass_solar, event_horizon_radius_km=radius_km)
    )
    assert any("mass_solar" in issue for issue in issues)


def test_check_hard_invariants_flags_black_hole_disk_luminosity_mismatch():
    issues = pp.check_hard_invariants(_base_record("black-hole", has_accretion_disk=True, luminosity_w=0.0))
    assert any("luminosity_w" in issue for issue in issues)


def test_check_hard_invariants_flags_neutron_star_luminosity_mismatch():
    issues = pp.check_hard_invariants(_base_record("neutron-star", luminosity_w=1.0))
    assert any("luminosity_w" in issue for issue in issues)


def test_check_hard_invariants_flags_nebula_radius_outside_its_type_range():
    issues = pp.check_hard_invariants(_base_record("nebula", nebula_type="planetary", radius_ly=50.0))
    assert any("radius_ly" in issue for issue in issues)


def test_check_hard_invariants_flags_supernova_remnant_sedov_taylor_mismatch():
    issues = pp.check_hard_invariants(_base_record("supernova-remnant", radius_ly=999.0))
    assert any("radius_ly" in issue for issue in issues)


def test_check_hard_invariants_flags_type_ia_with_compact_remnant():
    issues = pp.check_hard_invariants(
        _base_record("supernova-remnant", progenitor_type="Type Ia", has_compact_remnant=True)
    )
    assert any("Type Ia" in issue for issue in issues)


def test_check_hard_invariants_flags_rogue_planet_type_mismatch():
    issues = pp.check_hard_invariants(_base_record("rogue-planet", mass_kg=1e20, planet_type="g"))
    assert any("planet_type" in issue for issue in issues)


def test_check_hard_invariants_flags_comet_velocity_out_of_range():
    issues = pp.check_hard_invariants(_base_record("comet", velocity_kms=500.0))
    assert any("velocity_kms" in issue for issue in issues)


def test_check_hard_invariants_flags_asteroid_field_bad_density():
    issues = pp.check_hard_invariants(_base_record("asteroid-field", density="chunky"))
    assert any("density" in issue for issue in issues)
