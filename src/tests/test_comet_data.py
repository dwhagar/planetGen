"""
Star-bound comet regression tests
==================================

Covers `stellarObjects/cometData.Comet` -- a comet gravitationally bound
to a star (contrast `roguePlanetData.InterstellarComet`, covered by
`test_phenomena.py`, which is always standalone/unbound). See
docs/design/comet-orbital-realism.md for the design this implements.

Run with: pytest src/tests/test_comet_data.py
"""

import math

import pytest

from stellarObjects import physical_constants, program_constants
from stellarObjects.cometData import Comet
from stellarObjects.config import SystemConfig

TRIALS = 20


def make_config(**overrides):
    cfg = SystemConfig()
    for attr, value in overrides.items():
        setattr(cfg, attr, value)
    return cfg


# ---------------------------------------------------------------------------
# Generation sanity.
# ---------------------------------------------------------------------------

def test_elliptical_comet_within_configured_ranges():
    for _ in range(TRIALS):
        comet = Comet(make_config(), primary_mass_solar=1.0, orbit_type="elliptical")
        assert comet.orbit_type == "elliptical"
        assert comet.period_class in program_constants.COMET_PERIOD_CLASSES
        class_data = program_constants.COMET_PERIOD_CLASSES[comet.period_class]
        assert class_data["eccentricity_range"][0] <= comet.eccentricity <= class_data["eccentricity_range"][1]
        assert 0 <= comet.inclination_deg <= class_data["inclination_max_deg"]
        assert program_constants.BOUND_COMET_NUCLEUS_DIAMETER_RANGE_KM[0] <= comet.nucleus_diameter_km <= program_constants.BOUND_COMET_NUCLEUS_DIAMETER_RANGE_KM[1]
        assert program_constants.COMET_PERIHELION_DISTANCE_RANGE_AU[0] <= comet.perihelion_distance_au <= program_constants.COMET_PERIHELION_DISTANCE_RANGE_AU[1]
        assert comet.orbital_period_years is not None and comet.orbital_period_years > 0
        assert comet.parabolic_mean_anomaly is None
        assert 0 <= comet.mean_anomaly_deg < 360
        assert comet.min_update_interval_years is not None and comet.min_update_interval_years > 0
        assert 0 < len(comet.composition) <= 3
        assert all(component in program_constants.COMET_COMPOSITION for component in comet.composition)


def test_parabolic_comet_within_configured_ranges():
    for _ in range(TRIALS):
        comet = Comet(make_config(), primary_mass_solar=1.0, orbit_type="parabolic")
        assert comet.orbit_type == "parabolic"
        assert comet.period_class is None
        assert program_constants.PARABOLIC_COMET_ECCENTRICITY_RANGE[0] <= comet.eccentricity <= program_constants.PARABOLIC_COMET_ECCENTRICITY_RANGE[1]
        assert 0 <= comet.inclination_deg <= program_constants.PARABOLIC_COMET_INCLINATION_MAX_DEG
        assert comet.orbital_period_years is None
        assert comet.mean_anomaly_deg is None
        assert comet.parabolic_mean_anomaly is not None
        assert comet.min_update_interval_years is None


def test_comet_never_closer_than_its_own_perihelion():
    # distance_au is this comet's CURRENT distance, sampled at some random
    # point along its orbit -- it must never be less than the orbit's own
    # closest-approach distance, regardless of orbit_type.
    for orbit_type in ("elliptical", "parabolic"):
        for _ in range(TRIALS):
            comet = Comet(make_config(), primary_mass_solar=1.0, orbit_type=orbit_type)
            assert comet.distance_au >= comet.perihelion_distance_au - 1e-9


def test_comet_orbital_speed_and_position_are_finite_and_positive():
    for orbit_type in ("elliptical", "parabolic"):
        comet = Comet(make_config(), primary_mass_solar=1.0, orbit_type=orbit_type)
        assert math.isfinite(comet.orbital_speed_kms) and comet.orbital_speed_kms > 0
        assert math.isfinite(comet.distance_au) and comet.distance_au > 0
        for coord in (comet.position_x_au, comet.position_y_au, comet.position_z_au):
            assert math.isfinite(coord)


def test_orbit_type_defaults_to_a_random_roll_without_override():
    # Without an explicit orbit_type, generation should be capable of
    # producing both kinds across enough trials (a smoke test against the
    # override path silently becoming the only path).
    seen_types = {Comet(make_config(), primary_mass_solar=1.0).orbit_type for _ in range(60)}
    assert seen_types == {"elliptical", "parabolic"}


# ---------------------------------------------------------------------------
# Serialization round-trips (to_dict/from_dict).
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("orbit_type", ["elliptical", "parabolic"])
def test_to_dict_from_dict_round_trips_render_identically(orbit_type):
    cfg = make_config()
    comet = Comet(cfg, primary_mass_solar=1.0, orbit_type=orbit_type)
    reloaded = Comet.from_dict(comet.to_dict(), cfg)
    assert str(reloaded) == str(comet)
    assert reloaded.to_dict() == comet.to_dict()


# ---------------------------------------------------------------------------
# Paragraph rendering.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("orbit_type", ["elliptical", "parabolic"])
@pytest.mark.parametrize("markdown", [True, False])
def test_paragraph_list_renders_section_header(orbit_type, markdown):
    cfg = make_config(MARKDOWN=markdown)
    comet = Comet(cfg, primary_mass_solar=1.0, orbit_type=orbit_type)
    paragraphs = comet.to_paragraph_list()
    assert paragraphs
    assert all(isinstance(p, str) and p for p in paragraphs)
    header = paragraphs[0]
    if markdown:
        assert header.startswith("##")
        assert not header.endswith("==")
    else:
        assert header.startswith("==")
        assert header.endswith("==")
    assert comet.name in header
    assert str(comet) == "\n\n".join(paragraphs)


def test_elliptical_paragraph_mentions_period_and_orbit_type():
    comet = Comet(make_config(), primary_mass_solar=1.0, orbit_type="elliptical")
    description = comet.to_paragraph_list()[1]
    assert "elliptical" in description
    assert "period" in description


def test_parabolic_paragraph_mentions_it_will_not_return():
    comet = Comet(make_config(), primary_mass_solar=1.0, orbit_type="parabolic")
    description = comet.to_paragraph_list()[1]
    assert "will not return" in description
