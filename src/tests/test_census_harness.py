"""
TEST.186: the census audit harness (`tests.audit.census`): the Poisson rule, the 'expected 20 or more appears' rule,
the comparison over a census, the one report, and a real fill of seeded sectors counted by stage and kind.
"""

import math
import random

import pytest

from planetgen import tuning
from planetgen.galaxy.density import build_galaxy_shape
from planetgen.generation import fill_stages
from planetgen.physics.units import ly_to_pc
from tests.toy_shape import toy_shape
from tests.audit import census as audit

EDGE_PC = ly_to_pc(tuning.DEFAULT_SECTOR_EDGE_LY)
SHAPE = toy_shape(disk_scale_length_pc=40.0, disk_scale_height_pc=12.0, bulge_scale_radius_pc=10.0,
                           bulge_amplitude=2.0, arm_count=2, pitch_angle_rad=math.radians(15), arm_amplitude=0.4)


def _region(name, totals, sectors=10):
    region = audit.RegionCensus(name)
    for i in range(sectors):
        region.add({key: n // sectors + (1 if i < n % sectors else 0) for key, n in totals.items()})
    return region


# --- The rules ----------------------------------------------------------------------------------------------------


def test_within_sigma_is_k_poisson_sigma_either_side_of_the_expectation():
    assert audit.within_sigma(100, 100)
    assert audit.within_sigma(140, 100)         # 4 sigma
    assert not audit.within_sigma(141, 100)
    assert audit.within_sigma(60, 100)
    assert not audit.within_sigma(59, 100)
    assert audit.within_sigma(141, 100, k=4.2)
    assert not audit.within_sigma(130, 100, k=2)


def test_where_none_is_expected_only_none_passes():
    assert audit.within_sigma(0, 0.0)
    assert not audit.within_sigma(1, 0.0)
    assert audit.sigma_distance(0, 0.0) == 0.0
    assert math.isinf(audit.sigma_distance(1, 0.0))
    assert audit.sigma_distance(110, 100) == pytest.approx(1.0)


def test_a_kind_expected_twenty_times_must_appear():
    assert audit.must_appear(20.0, 0)
    assert not audit.must_appear(19.9, 0)
    assert not audit.must_appear(20.0, 1)
    assert audit.must_appear(5.0, 0, at=5)


def test_poisson_counts_pass_four_sigma_almost_always():
    rng = random.Random(186)

    def poisson(mean):
        limit, product, n = math.exp(-mean), 1.0, 0
        while True:
            product *= rng.random()
            if product <= limit:
                return n
            n += 1

    passed = sum(audit.within_sigma(poisson(60.0), 60.0) for _ in range(400))
    assert passed >= 398


# --- The comparison -----------------------------------------------------------------------------------------------


def _census():
    return audit.Census([
        _region("disk", {("scattered", "black-hole"): 30, ("systems", "star-system"): 200,
                         ("flags", "fast-star"): 1}),
        _region("halo", {("scattered", "black-hole"): 0, ("systems", "star-system"): 3,
                         ("nebulae", "molecular-cloud"): 4}),
    ])


def test_compare_passes_counts_within_bounds_and_names_the_ones_that_break_a_rule():
    expectations = {
        "disk": {"black-hole": 3.0, ("systems", "star-system"): 20.0, "wolf-rayet": 2.0},
        "halo": {"black-hole": 5.0, ("systems", "star-system"): 0.0, "molecular-cloud": 0.0},
    }
    findings = audit.compare(_census(), expectations)
    got = {(f.region, f.label(), f.problem): (f.observed, f.expected) for f in findings}
    assert got == {
        ("disk", "wolf-rayet", audit.NEVER_APPEARED): (0, 20.0),
        ("halo", "black-hole", audit.NEVER_APPEARED): (0, 50.0),
        ("halo", "systems/star-system", audit.FORBIDDEN): (3, 0.0),
        ("halo", "molecular-cloud", audit.FORBIDDEN): (4, 0.0),
    }


def test_compare_flags_a_count_out_of_bounds_and_skips_regions_with_no_expectation():
    findings = audit.compare(_census(), {"disk": {"black-hole": 1.0}})
    assert [(f.region, f.problem, f.observed, f.expected) for f in findings] == [
        ("disk", audit.OUT_OF_BOUNDS, 30, 10.0)]
    assert audit.compare(_census(), {}) == []
    assert audit.compare(_census(), {"disk": {"black-hole": 1.0}}, k=10) == []


def test_a_kind_key_adds_its_stages_and_a_stage_kind_key_reads_one():
    region = _region("disk", {("scattered", "black-hole"): 6, ("below-cut", "black-hole"): 4})
    assert region.observed("black-hole") == 10
    assert region.observed(("below-cut", "black-hole")) == 4
    assert region.observed("neutron-star") == 0
    assert region.kinds() == ["black-hole"]
    assert region.sectors_with("black-hole") == 6
    assert region.sectors_with(("below-cut", "black-hole")) == 4


def test_kinds_made_with_no_expectation_are_listed_and_not_failed():
    assert audit.unlisted(_census(), {"disk": {"black-hole": 3.0, ("systems", "star-system"): 20.0}}) == {
        "disk": ["fast-star"]}


# --- The report ---------------------------------------------------------------------------------------------------


def test_the_report_is_one_text_with_every_expected_count_and_its_verdict():
    expectations = {"disk": {"black-hole": 3.0, "wolf-rayet": 2.0}, "halo": {"black-hole": 5.0}}
    findings = audit.compare(_census(), expectations)
    text = audit.format_report(_census(), expectations, findings)
    assert text.startswith("Census audit")
    assert "disk: 10 sectors" in text and "halo: 10 sectors" in text
    assert "wolf-rayet" in text and "never appeared" in text
    assert "made, with no expectation: fast-star, star-system" in text
    assert "2 of the expected counts broke a rule" in text
    assert text == audit.format_report(_census(), expectations, findings)


def test_the_report_without_expectations_lists_the_counts_and_a_clean_run_says_so():
    assert "black-hole  observed       30" in audit.format_report(_census())
    clean = audit.format_report(_census(), {"disk": {"black-hole": 3.0}})
    assert clean.rstrip().endswith("Every expected count is within bounds")


def test_the_report_is_written_to_a_file(tmp_path):
    path = tmp_path / "census.txt"
    text = audit.write_report(path, _census(), {"disk": {"black-hole": 3.0}})
    assert path.read_text(encoding="utf-8") == text


# --- A real fill --------------------------------------------------------------------------------------------------


CENTERS = {"inner": (100.0, 0.0, 0.0), "outer": (900.0, 0.0, 300.0)}


def _cases(*names, sectors=3):
    return [audit.RegionCase(name, CENTERS[name], SHAPE, EDGE_PC, 2000.0, sectors=sectors, num_systems=2)
            for name in names]


def test_a_census_counts_each_seeded_sector_by_stage_and_kind():
    census = audit.run_census(_cases("inner", "outer"))
    assert list(census.regions) == ["inner", "outer"]
    for region in census:
        assert region.sectors == len(region.per_sector) == 3
        assert region.totals == sum(region.per_sector, audit.collections.Counter())
        assert region.observed((fill_stages.SYSTEMS, "star-system")) == 6
        assert set(stage for stage, _kind in region.totals) <= set(fill_stages.STAGES)


def test_a_census_is_the_same_on_every_run_and_regions_do_not_move_each_other():
    both = audit.run_census(_cases("inner", "outer"))
    again = audit.run_census(_cases("inner", "outer"))
    alone = audit.run_census(_cases("outer"))
    assert [r.per_sector for r in both] == [r.per_sector for r in again]
    assert alone["outer"].per_sector == both["outer"].per_sector
    assert audit.seed_for("outer", 1) != audit.seed_for("outer", 2) != audit.seed_for("inner", 2)


def test_a_case_may_bring_its_own_fill():
    seen = []

    class Sector:
        stage_counts = {("systems", "star-system"): 1}

    def fill(case, seed, args):
        seen.append(seed)
        return Sector()

    census = audit.run_census(_cases("inner", sectors=2), fill=fill, args=object())
    assert seen == [audit.seed_for("inner", 0), audit.seed_for("inner", 1)]
    assert census["inner"].observed("star-system") == 2
