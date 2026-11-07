# tests/test_field_host.py

"""The asteroid field a young supernova remnant hosts (GEN.1104, `generation/field_host.py`): only a core-collapse remnant, young enough that the field has not passed a light year, hosts one; its
columns are the rule's and its radius follows from them; the same remnant hosts the same field; a scattered remnant of a sector brings the field with it."""

import pytest

from planetgen.galaxy.sector import SpaceSector
from planetgen.generation import field_host, run_sector
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.asteroid_field import RING_RADIUS_AU, radius_ly_at
from planetgen.generation.phenomena.supernova_remnant import SupernovaRemnant
from planetgen.physics import constants, field_dispersal, field_mass, field_rates
from planetgen.util import draw


def _remnant(seed):
    with draw.bound(f"remnant:{seed}"):
        return SupernovaRemnant(SystemConfig())


def _hosts(count):
    """The first `count` core-collapse remnants that host a field, as `(remnant, hosted)`."""
    found = []
    for seed in range(1, 4000):
        remnant = _remnant(seed)
        hosted = field_host.hosted_field(remnant, SystemConfig(), f"k:{seed}")
        if hosted is not None:
            found.append((remnant, hosted))
            if len(found) == count:
                return found
    pytest.fail(f"fewer than {count} hosts in 4000 remnants")


def test_a_type_ia_remnant_hosts_nothing():
    for seed in range(1, 200):
        remnant = _remnant(seed)
        if remnant.progenitor_type == "Type Ia":
            assert field_host.hosted_field(remnant, SystemConfig(), f"k:{seed}") is None
            return
    pytest.fail("no Type Ia remnant in 200")


def test_a_hosted_field_has_the_columns_of_the_rule_and_the_radius_they_give():
    for remnant, hosted in _hosts(6):
        field = hosted.field
        assert remnant.progenitor_type == "core-collapse"
        assert field.shape == "annulus" and field.origin == "unbound belt" and "unbound" in field.history
        assert field.birth_age_yr == remnant.age_years and field.rms_speed_kms == pytest.approx(hosted.speed_ms / 1000.0)
        assert field.radius_ly == pytest.approx(radius_ly_at(field.birth_age_yr, field.rms_speed_kms, "annulus"), rel=1e-12)
        assert field.radius_ly < 1.0                  # not past a light year: it would be retired
        assert not field_dispersal.is_retired(field.birth_age_yr, hosted.speed_ms, hosted.belt_radius_m)
        assert field.mass_kg == hosted.mass.mass_kg >= field_mass.MIN_FIELD_MASS_EARTH * constants.EARTH_MASS_TO_KG
        assert 0.0 < field.largest_body_kg <= field.mass_kg


def test_the_same_remnant_hosts_the_same_field():
    for seed in range(1, 4000):
        remnant = _remnant(seed)
        first = field_host.hosted_field(remnant, SystemConfig(), f"k:{seed}")
        if first is not None:
            again = field_host.hosted_field(remnant, SystemConfig(), f"k:{seed}")
            assert first.field.to_dict() == again.field.to_dict() and first.mass == again.mass
            return
    pytest.fail("no hosting remnant in 4000")


def test_a_remnant_older_than_the_field_takes_to_pass_a_light_year_hosts_nothing():
    remnant, hosted = _hosts(1)[0]
    remnant.age_years = 1.0e9
    for key in range(200):
        assert field_host.hosted_field(remnant, SystemConfig(), f"old:{key}") is None


def test_the_share_of_core_collapse_remnants_that_host_a_field_is_under_the_share_that_kept_a_belt():
    core_collapse = hosted = 0
    for seed in range(1, 600):
        remnant = _remnant(seed)
        if remnant.progenitor_type != "core-collapse":
            continue
        core_collapse += 1
        hosted += field_host.hosted_field(remnant, SystemConfig(), f"k:{seed}") is not None
    assert 0 < hosted < core_collapse * max(share for _age, share in field_rates.PROGENITOR_BELT_SHARE_BY_AGE)


def test_a_ring_spreads_faster_in_radius_than_a_ball_at_young_ages_and_alike_when_old():
    young = radius_ly_at(100.0, 10.0, "annulus"), radius_ly_at(100.0, 10.0, "ball")
    assert young[0] > young[1] and young[0] >= RING_RADIUS_AU * constants.AU_M / constants.LIGHTYEAR_M
    old = radius_ly_at(1.0e5, 10.0, "annulus"), radius_ly_at(1.0e5, 10.0, "ball")
    assert old[0] == pytest.approx(old[1], rel=1e-3)


def test_a_scattered_remnant_of_a_sector_brings_its_field_with_it():
    sector = SpaceSector("Host Sector", edge_ly=10.0)
    for seed in range(1, 4000):
        row = {"seed": f"scatter:{seed}"}
        with draw.bound(row["seed"]):
            remnant = SupernovaRemnant(SystemConfig())
        if field_host.hosted_field(remnant, remnant.system_config, str(row["seed"])) is None:
            continue
        entry = sector.add_phenomenon(remnant, "supernova-remnant", position=(1.0, 2.0, 3.0))
        made = run_sector._hosted_field_entry(sector, entry, row)
        assert made is not None and made.phenomenon_type == "asteroid-field" and made.position == entry.position
        assert made.phenomenon.origin == "unbound belt" and made in sector.phenomena
        return
    pytest.fail("no hosting remnant in 4000")
