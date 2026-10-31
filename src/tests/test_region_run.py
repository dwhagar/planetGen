# tests/test_region_run.py

"""GEN.1128: the plan of a region run (`planetgen/galaxy/region_run.py`): a sphere from a start point as whole sectors, units of a block along one
(shell, band) arc, nearest first, resume by skipping what is charted, the barrier between stages 0 to 2 and the top-up, and the cost."""

import math

import pytest

from planetgen.galaxy import region as region_module
from planetgen.galaxy import region_run as rr
from planetgen.galaxy import spherical

EDGE = 4.0
CENTRES = [(0.0, 0.0, 0.0), (30.0, -12.0, 5.0), (8.0, 0.0, 0.0)]


def _addresses(plan):
    out = []
    for unit in plan.units:
        out.extend(tuple(rr.unit_items(unit, spherical.wedge_count)[i:i + 3]) for i in range(0, 3 * unit.count, 3))
    return out


def _ball(center, radius, grow=0):
    return region_module.Region(center, rr.reach_pc(radius, EDGE, grow), EDGE)


@pytest.mark.parametrize("center", CENTRES)
def test_the_units_hold_every_sector_of_the_region_once(center):
    plan = rr.plan(center, 20.0, EDGE)
    addresses = _addresses(plan)
    assert len(addresses) == len(set(addresses)) == plan.sectors
    assert set(addresses) == set(_ball(center, 20.0).addresses())


@pytest.mark.parametrize("block", [1, 5, 27, 200])
def test_a_unit_is_a_run_of_consecutive_wedges_of_one_band_no_larger_than_the_block(block):
    plan = rr.plan(CENTRES[1], 24.0, EDGE, block_sectors=block)
    assert plan.units and all(1 <= unit.count <= block for unit in plan.units)
    for unit in plan.units:
        wedges = spherical.wedge_count(unit.shell, unit.band)
        triples = rr.unit_items(unit, spherical.wedge_count)
        assert triples[0::3] == [unit.shell] * unit.count and triples[1::3] == [unit.band] * unit.count
        assert triples[2::3] == [(unit.first + k) % wedges for k in range(unit.count)]


def test_the_units_run_nearest_first_and_the_plan_is_the_same_every_time():
    plan = rr.plan(CENTRES[1], 24.0, EDGE)
    nearest = [unit.nearest_pc for unit in plan.units]
    assert nearest == sorted(nearest)
    assert rr.plan(CENTRES[1], 24.0, EDGE) == plan


def test_a_unit_never_changes_what_a_sector_is_so_two_regions_cut_the_same_sector_the_same_way_in_every_block_size():
    """The sector is the unit of randomness: the sectors a plan holds depend on the sphere only, never on the cut."""
    small = set(_addresses(rr.plan(CENTRES[1], 20.0, EDGE, block_sectors=3)))
    large = set(_addresses(rr.plan(CENTRES[1], 20.0, EDGE, block_sectors=500)))
    assert small == large


def test_a_stopped_run_resumes_with_exactly_the_sectors_not_yet_charted():
    plan = rr.plan(CENTRES[0], 20.0, EDGE)
    charted = set()
    for unit in plan.units[: len(plan.units) // 2]:
        charted.update(tuple(rr.unit_items(unit, spherical.wedge_count)[i:i + 3]) for i in range(0, 3 * unit.count, 3))
    resumed = rr.plan(CENTRES[0], 20.0, EDGE, skip=charted)
    addresses = set(_addresses(resumed))
    assert addresses.isdisjoint(charted)
    assert addresses | charted == set(_addresses(plan))
    assert resumed.skipped == len(charted) and resumed.sectors == plan.sectors - len(charted)


def test_a_second_overlapping_region_costs_only_its_new_sectors():
    first = set(_addresses(rr.plan(CENTRES[0], 16.0, EDGE)))
    second = rr.plan((6.0, 0.0, 0.0), 16.0, EDGE, skip=first)
    assert set(_addresses(second)) == set(_ball((6.0, 0.0, 0.0), 16.0).addresses()) - first
    assert 0 < second.sectors < len(first)


def test_a_skipped_sector_splits_the_run_it_sat_in_and_no_unit_spans_it():
    plan = rr.plan(CENTRES[0], 20.0, EDGE, block_sectors=500)
    victim = next(unit for unit in plan.units if unit.count >= 3)
    wedges = spherical.wedge_count(victim.shell, victim.band)
    middle = (victim.shell, victim.band, (victim.first + 1) % wedges)
    cut = rr.plan(CENTRES[0], 20.0, EDGE, block_sectors=500, skip={middle})
    for unit in cut.units:
        assert middle not in {tuple(rr.unit_items(unit, spherical.wedge_count)[i:i + 3]) for i in range(0, 3 * unit.count, 3)}


def test_growing_shells_adds_sectors_in_the_new_shells_only():
    base = set(_addresses(rr.plan(CENTRES[0], 20.0, EDGE)))
    grown = set(_addresses(rr.plan(CENTRES[0], 20.0, EDGE, grow_shells=2)))
    assert base < grown
    assert math.isclose(rr.plan(CENTRES[0], 20.0, EDGE, grow_shells=2).reach_pc - rr.plan(CENTRES[0], 20.0, EDGE).reach_pc, 2 * EDGE)
    farthest = max(math.dist(spherical.sector_position_pc(*a, EDGE), CENTRES[0]) for a in base)
    assert all(math.dist(spherical.sector_position_pc(*a, EDGE), CENTRES[0]) >= farthest - 3 * EDGE for a in grown - base)


def test_a_shell_costs_about_the_surface_it_adds():
    """Section 9: growing by shells costs about r squared, roughly 3,100 sectors for the first 4 pc shell around a 50 pc sphere."""
    inner = rr.plan(CENTRES[0], 50.0, EDGE)
    outer = rr.plan(CENTRES[0], 50.0, EDGE, grow_shells=1)
    added = outer.sectors - inner.sectors
    assert 2000 < added < 5000, added


def test_the_sector_count_is_the_sphere_volume_over_a_sector_s_within_the_allowance():
    cost = rr.cost(rr.plan(CENTRES[0], 50.0, EDGE))
    assert 8000 < cost["volume_sectors"] < 8400
    assert cost["volume_sectors"] <= cost["sectors"] < 1.35 * cost["volume_sectors"]
    assert cost["units"] >= cost["sectors"] / rr.BLOCK_SECTORS and cost["largest_unit"] <= rr.BLOCK_SECTORS
    assert "seconds" not in cost
    assert rr.cost(rr.plan(CENTRES[0], 10.0, EDGE), seconds_per_sector=2.0)["seconds"] == 2.0 * rr.plan(CENTRES[0], 10.0, EDGE).sectors


def test_the_barrier_holds_the_top_up_until_every_unit_of_the_set_is_generated():
    plan = rr.plan(CENTRES[0], 12.0, EDGE)
    keys = [rr.key(unit) for unit in plan.units]
    stage, pending = rr.next_units(plan, topup=True)
    assert stage == rr.STAGE_GENERATE and pending == list(plan.units)
    stage, pending = rr.next_units(plan, done_generate=keys[:-1], topup=True)
    assert stage == rr.STAGE_GENERATE and [rr.key(u) for u in pending] == keys[-1:]
    stage, pending = rr.next_units(plan, done_generate=keys, topup=True)
    assert stage == rr.STAGE_TOPUP and len(pending) == len(keys)
    stage, pending = rr.next_units(plan, done_generate=keys, done_topup=keys[:2], topup=True)
    assert stage == rr.STAGE_TOPUP and len(pending) == len(keys) - 2
    assert rr.next_units(plan, done_generate=keys, done_topup=keys, topup=True) == (rr.STAGE_DONE, [])


def test_a_run_with_no_top_up_stage_is_done_when_every_unit_is_generated():
    plan = rr.plan(CENTRES[0], 12.0, EDGE)
    assert rr.next_units(plan, done_generate=[rr.key(u) for u in plan.units]) == (rr.STAGE_DONE, [])


def test_a_top_up_unit_done_before_the_set_was_generated_does_not_open_the_stage():
    plan = rr.plan(CENTRES[0], 12.0, EDGE)
    keys = [rr.key(u) for u in plan.units]
    stage, _pending = rr.next_units(plan, done_generate=keys[1:], done_topup=keys, topup=True)
    assert stage == rr.STAGE_GENERATE


def test_bad_arguments_are_refused():
    with pytest.raises(ValueError):
        rr.plan(CENTRES[0], -1.0, EDGE)
    with pytest.raises(ValueError):
        rr.plan(CENTRES[0], 10.0, EDGE, block_sectors=0)
    with pytest.raises(ValueError):
        rr.plan(CENTRES[0], 10.0, EDGE, grow_shells=-1)
    with pytest.raises(ValueError):
        rr.plan(CENTRES[0], 10.0, EDGE, grow_shells=1.5)


def test_a_zero_radius_holds_the_start_sector_and_its_face_neighbours_by_the_allowance():
    plan = rr.plan(CENTRES[0], 0.0, EDGE)
    assert plan.sectors >= 1 and plan.units[0].nearest_pc <= EDGE
