"""OPS.202 (part 2): `scatter-backfill`, `refresh-run-blocks`, `generate-start-sector`, `generate-neighborhood`, `generate-sectors` and `generate-system` run as pieces under the engine's loop."""

import os
import types

import pytest

from planetgen.generation import run_galaxy, run_plan, run_system
from planetgen.queue import commands

CONVERTED = ("scatter-backfill", "refresh-run-blocks", "generate-start-sector", "generate-neighborhood", "generate-sectors", "generate-system")


@pytest.mark.parametrize("kind", CONVERTED)
def test_the_galaxy_commands_have_pieces_and_no_loop_of_their_own(kind):
    spec = commands.REGISTRY[kind]
    assert spec.pieces is not None and spec.do is not None and spec.gather is not None and spec.unit is None
    for function in (commands.pieces_function(kind), commands.do_function(kind), commands.gather_function(kind)):
        assert callable(function)


def test_a_backfill_unit_has_a_piece_for_each_chunk_of_its_sectors_in_sorted_order():
    size = run_plan.BACKFILL_CHUNK_SECTORS
    items = [[2, 0, wedge, 1.5] for wedge in range(2 * size + 3)][::-1]       # out of order on purpose
    pieces = run_galaxy.backfill_pieces({"items": items}, 0)
    assert [len(piece["items"]) for piece in pieces] == [size, size, 3]
    flat = [item for piece in pieces for item in piece["items"]]
    assert flat == sorted(items, key=lambda item: tuple(item[:3]))
    assert run_galaxy.backfill_pieces({"items": items}, 1) == pieces[1:]
    assert run_galaxy.backfill_pieces({"items": []}, 0) == []


def test_a_sky_unit_is_one_piece_and_its_gather_is_that_pieces_result():
    unit = {"items": [], "sky": {"center": [0, 0, 0], "limit": 7.5}}
    assert run_galaxy.backfill_pieces(unit, 0) == [{"sky": True}] and run_galaxy.backfill_pieces(unit, 1) == []
    result = {"sectors": 0, "stars": 3, "shells": {}, "sky": {"stars": 3}}
    assert run_galaxy.backfill_gather(unit, [result], None) is result


def test_the_backfill_result_adds_the_pieces_in_the_order_the_kinds_first_appeared():
    unit = {"items": [[1, 0, 0, 1.0]] * 5}
    first = {"sectors": 2, "stars": 5, "shells": {"3": [["O star", "O stars", 2], ["B star", "B stars", 3]]}}
    second = {"sectors": 1, "stars": 4, "shells": {"3": [["B star", "B stars", 1], ["A star", "A stars", 1]], "4": [["O star", "O stars", 2]]}}
    result = run_galaxy.backfill_gather(unit, [first, second], None)
    assert result == {"sectors": 3, "stars": 9, "visited": 5,
                      "shells": {"3": [["O star", "O stars", 2], ["B star", "B stars", 4], ["A star", "A stars", 1]], "4": [["O star", "O stars", 2]]}}
    assert run_galaxy.backfill_gather({"items": []}, [], None) == {"sectors": 0, "stars": 0, "visited": 0, "shells": {}}


def test_a_refresh_run_blocks_unit_has_the_nebulae_piece_then_the_sectors_piece():
    assert run_galaxy.refresh_run_blocks_pieces({}, 0) == [{"run": "nebulae"}, {"run": "sectors"}]
    assert run_galaxy.refresh_run_blocks_pieces({}, 1) == [{"run": "sectors"}]
    gathered = run_galaxy.refresh_run_blocks_gather({}, [{"nebulae": 4, "warnings": []}, {"sectors": 7, "warnings": ["a warning"]}], None)
    assert gathered == {"nebulae": 4, "sectors": 7, "warnings": ["a warning"]}


def test_a_sector_unit_has_a_piece_for_each_sector_and_none_for_the_ones_saved(monkeypatch):
    monkeypatch.setattr(run_galaxy, "_sector_args_of", lambda unit: None)       # (the unit's arguments are read only to find the saved sectors)
    unit = {"args": {}, "items": [1, 0, 0, 1, 0, 1, 1, 0, 2], "suffixes": ["a", "b", "c"]}
    pieces = run_galaxy.generate_sectors_pieces(unit, 0)
    assert pieces == [{"address": [1, 0, 0], "suffix": "a"}, {"address": [1, 0, 1], "suffix": "b"}, {"address": [1, 0, 2], "suffix": "c"}]
    assert run_galaxy.generate_sectors_pieces(unit, 2) == pieces[2:]
    monkeypatch.setattr(run_galaxy, "_saved_addresses", lambda args: {(1, 0, 1)})
    assert run_galaxy.generate_sectors_pieces({**unit, "again": True}, 0) == [pieces[0], pieces[2]]       # the suffix goes by place, so a skipped sector leaves the others theirs
    assert run_galaxy.generate_sectors_pieces({"args": {}, "items": [1, 0, 0]}, 0) == [{"address": [1, 0, 0], "suffix": ""}]


def test_the_sector_unit_result_is_what_each_sector_reported_and_this_process():
    results = [{"address": [1, 0, 0]}, {"address": [1, 0, 1]}]
    assert run_galaxy.generate_sectors_gather({}, results, None) == {"sectors": results, "pid": os.getpid()}


def test_a_system_unit_generates_then_finishes_and_its_result_is_the_finishing_pieces():
    assert run_system.system_pieces({}, 0) == [{"run": "generate"}, {"run": "finish"}] and run_system.system_pieces({}, 1) == [{"run": "finish"}]
    assert run_system.system_gather({}, [{}, {"name": "S", "star_system_id": 3}], None) == {"name": "S", "star_system_id": 3}


def test_a_system_unit_draws_its_tries_from_one_stream_across_its_two_pieces(monkeypatch):
    """The system is made on the stream seeded from the run's, and rendered on the same one, as the one unit function did."""
    from planetgen.queue import engine
    from planetgen.util import draw
    seen = []
    system = types.SimpleNamespace(unmet_requirements=[], name="Test")
    monkeypatch.setattr(run_system, "StarSystem", lambda system_config: seen.append(draw.current()) or system)
    monkeypatch.setattr(run_system, "build_system_config", lambda args: types.SimpleNamespace(STAR_TYPE=None))
    monkeypatch.setattr(run_system, "class_request_report", lambda system, config: None)
    monkeypatch.setattr(run_system, "render_star_system", lambda system, fmt: "page")
    monkeypatch.setattr(run_galaxy, "_sector_args_of", lambda unit: types.SimpleNamespace(output="x", format="text", strict=False))
    monkeypatch.setattr(run_system.run_common, "_strict", lambda args: False)
    unit = {"run_seed": 42}
    ctx = engine.Context("generate-system", unit, lambda done, total: None)
    assert run_system.system_do(unit, {"run": "generate"}, ctx) == {}
    result = run_system.system_do(unit, {"run": "finish"}, ctx)
    assert result == {"warning": None, "name": "Test", "text": "page", "star_system_id": None, "asked": None}
    assert len(seen) == 1 and seen[0] is ctx.state["stream"]
