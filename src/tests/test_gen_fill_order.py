"""
GEN.101: the order a neighbourhood fill generates its sectors in (`planetgen.generation.fill_order`): the sector at
the centre first, then outward in blocks of four sectors a side along a Hilbert curve, deterministic whatever order the
candidates come in, and the steps that are not between neighbours counted and reported.
"""

import itertools
import random

import pytest

from planetgen import tuning
from planetgen.galaxy.spherical import enumerate_sectors_within_radius
from planetgen.generation import fill_order, run_galaxy
from planetgen.physics.units import ly_to_pc
from planetgen.queue import wrapper
from planetgen.util import log

EDGE_PC = ly_to_pc(tuning.DEFAULT_SECTOR_EDGE_LY)
CENTRES = {"the galactic centre": (0.0, 0.0, 0.0), "off the axis": (1000.0, 200.0, 50.0),
           "far out on the axis": (0.0, 0.0, 9000.0)}


def _cells(center, radius_pc):
    return list(enumerate_sectors_within_radius(center, radius_pc, EDGE_PC))


@pytest.mark.parametrize("bits", [1, 2, 3])
def test_the_hilbert_curve_visits_every_point_once_and_every_step_is_face_adjacent(bits):
    side = 1 << bits
    index = {point: fill_order.hilbert_index(point, bits) for point in itertools.product(range(side), repeat=3)}
    assert sorted(index.values()) == list(range(side ** 3))
    at = {value: point for point, value in index.items()}
    for i in range(side ** 3 - 1):
        assert sum(abs(a - b) for a, b in zip(at[i], at[i + 1])) == 1


def test_a_point_outside_the_lattice_is_refused():
    with pytest.raises(ValueError):
        fill_order.hilbert_index((4, 0, 0), 2)


@pytest.mark.parametrize("name", CENTRES)
def test_the_order_holds_every_candidate_once_and_starts_at_the_sector_nearest_the_centre(name):
    center = CENTRES[name]
    cells = _cells(center, 40.0)
    ordered = fill_order.fill_order(cells, center, EDGE_PC)
    assert sorted(ordered) == sorted(cells)
    nearest = min(cells, key=lambda cell: (round(cell[6], 9), cell[:3]))
    assert ordered[0] == nearest


@pytest.mark.parametrize("name", CENTRES)
def test_the_order_does_not_depend_on_the_order_the_candidates_came_in(name):
    center = CENTRES[name]
    cells = _cells(center, 40.0)
    shuffled = cells[:]
    random.Random(7).shuffle(shuffled)
    assert fill_order.fill_order(shuffled, center, EDGE_PC) == fill_order.fill_order(cells, center, EDGE_PC)
    assert fill_order.fill_order(reversed(cells), center, EDGE_PC) == fill_order.fill_order(cells, center, EDGE_PC)


@pytest.mark.parametrize("name", CENTRES)
def test_at_least_95_percent_of_the_steps_of_a_50_pc_ball_are_between_neighbours(name):
    center = CENTRES[name]
    ordered = fill_order.fill_order(_cells(center, 50.0), center, EDGE_PC)
    steps, jumps, mean, longest = fill_order.step_statistics(ordered, EDGE_PC)
    assert steps == len(ordered) - 1 > 5000
    assert (steps - jumps) / steps >= 0.95, (jumps, steps)
    assert mean > fill_order.JUMP_EDGES and longest >= mean       # a jump is longer than a neighbour step, by definition


def test_the_ball_is_filled_from_the_centre_outward_in_blocks():
    center = CENTRES["off the axis"]
    cells = _cells(center, 50.0)
    ordered = fill_order.fill_order(cells, center, EDGE_PC)
    # the nearest tenth of the sectors is, for the most part, among the first tenth to be filled
    tenth = len(cells) // 10
    first = {cell[:3] for cell in ordered[:tenth]}
    nearest = sorted(cells, key=lambda cell: cell[6])[:tenth]
    assert sum(cell[:3] in first for cell in nearest) >= 0.65 * tenth


def test_a_single_candidate_and_no_candidate_are_orders_of_their_own():
    assert fill_order.fill_order([], (0.0, 0.0, 0.0), EDGE_PC) == []
    one = [(0, 0, 0, 0.0, 0.0, 0.0, 0.0)]
    assert fill_order.fill_order(one, (0.0, 0.0, 0.0), EDGE_PC) == one
    assert fill_order.step_statistics(one, EDGE_PC) == (0, 0, 0.0, 0.0)
    assert fill_order.step_statistics([], EDGE_PC) == (0, 0, 0.0, 0.0)


def test_a_step_is_a_jump_only_beyond_two_and_a_half_edges():
    def cell(x):
        return (0, 0, 0, x * EDGE_PC, 0.0, 0.0, 0.0)
    ordered = [cell(0), cell(1), cell(3.4), cell(9.4)]       # steps of 1, 2.4 and 6 edges
    assert fill_order.step_statistics(ordered, EDGE_PC) == (3, 1, pytest.approx(6.0), pytest.approx(6.0))


def test_the_lattice_puts_the_centre_in_the_middle_of_its_block():
    assert fill_order.lattice_offset((0.0, 0.0, 0.0), (0.0, 0.0, 0.0), EDGE_PC) == (0, 0, 0)
    assert fill_order.lattice_offset((EDGE_PC * 1.4, -EDGE_PC * 0.6, 0.0), (0.0, 0.0, 0.0), EDGE_PC) == (1, -1, 0)
    block, local = fill_order._block_of((0, 0, 0))
    assert block == (0, 0, 0) and local == (2, 2, 2)
    assert fill_order._block_of((-2, -2, -2)) == ((0, 0, 0), (0, 0, 0))
    assert fill_order._block_of((2, 2, 2))[0] == (1, 1, 1)


# --- In the run ---------------------------------------------------------------------------------------------------


def test_a_local_neighborhood_generates_its_sectors_in_the_fill_order_and_says_how_it_walks(mysql_config, monkeypatch):
    from tests.test_galaxy_gen import _all_sectors, _mysql_argv, _plan_wide_galaxy, _run_cli
    _plan_wide_galaxy(mysql_config)
    _run_cli(["--shell", "1", "--band", "0", "--limit", "1", "--num-systems", "1"] + _mysql_argv(mysql_config))
    (seed,) = _all_sectors(mysql_config)
    center = (seed["center_x_pc"], seed["center_y_pc"], seed["center_z_pc"])
    submitted = []
    monkeypatch.setattr(run_galaxy, "_submit_batch",
                        lambda args, batch, *rest, **kwargs: submitted.append((args, list(batch), kwargs)))
    messages = []
    monkeypatch.setattr(log, "normal", lambda message, *args, **kwargs: messages.append(message))
    radius = 14.0
    _run_cli(["--center-sector", str(seed["id"]), "--radius-pc", str(radius), "--num-systems", "1"]
             + _mysql_argv(mysql_config))
    ((args, batch, options),) = submitted
    expected = [cell for cell in fill_order.fill_order(_cells(center, radius), center, EDGE_PC)
                if cell[:3] != (seed["shell_index"], seed["band_index"], seed["wedge_index"])]
    assert batch == [] and options["planned"] == "neighborhood" and options["count"] == len(expected)
    items, suffixes = run_galaxy.planned_items({
        "database": wrapper.database_args(mysql_config), "args": run_galaxy._sector_args_json(args),
        "edge_pc": EDGE_PC, "items": [], "planned": "neighborhood"})
    assert [tuple(items[place:place + 3]) for place in range(0, len(items), 3)] == [cell[:3] for cell in expected]
    assert len(suffixes) == len(expected) and suffixes[0].startswith(", ") and f"pc from sector_id={seed['id']}" in suffixes[0]
    (line,) = [m for m in messages if m.startswith("Fill order: ")]
    assert f"{len(expected) + 1:,} sector(s)" in line and "between neighbouring sectors" in line
    assert "jump(s)" in line


def test_a_caller_that_only_wants_the_set_gets_the_enumeration_unsorted(mysql_config):
    from tests.test_galaxy_gen import _plan_wide_galaxy
    _plan_wide_galaxy(mysql_config)
    bounds = run_galaxy._BatchDensity(mysql_config).bounds
    center = (0.0, 0.0, 0.0)
    plain, _occupied, _outside = run_galaxy._neighborhood_candidates(center, 12.0, EDGE_PC, mysql_config, bounds,
                                                                      ordered=False)
    ordered, _occupied, _outside = run_galaxy._neighborhood_candidates(center, 12.0, EDGE_PC, mysql_config, bounds)
    assert sorted(plain) == sorted(ordered) and plain != ordered
    assert ordered == fill_order.fill_order(plain, center, EDGE_PC)
