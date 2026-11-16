"""OPS.198: the scatters and fields (`scatter-mass`, `scatter-luminosity`, `scatter-band`, `scatter-phenomena`, `draw-nebula-field`, `fill-map-luminosity`) run as pieces under the engine's loop."""

import pytest

from planetgen.generation import bright_stars as brightStars
from planetgen.generation import run_plan
from planetgen.queue import commands

CONVERTED = ("scatter-mass", "scatter-luminosity", "scatter-band", "scatter-phenomena", "draw-nebula-field", "fill-map-luminosity")


@pytest.mark.parametrize("kind", CONVERTED)
def test_the_scatters_and_fields_have_pieces_and_no_loop_of_their_own(kind):
    spec = commands.REGISTRY[kind]
    assert spec.pieces is not None and spec.do is not None and spec.gather is not None and spec.unit is None
    for function in (commands.pieces_function(kind), commands.do_function(kind), commands.gather_function(kind)):
        assert callable(function)


def test_a_shell_scatter_unit_has_a_piece_for_each_shell_then_the_far_shells_then_the_components():
    args = {"shells": [[3, 2], [4, 1]], "far": [20, 30], "component_parts": [["cluster", "A", 0, 1]]}
    pieces = run_plan.scatter_shells_pieces(args, 0)
    assert pieces == [{"run": "shell", "shell": 3, "high": 2}, {"run": "shell", "shell": 4, "high": 1},
                      {"run": "far", "first": 20, "last": 30}, {"run": "components"}]
    assert run_plan.scatter_shells_pieces(args, 2) == pieces[2:]
    assert run_plan.scatter_shells_pieces({"shells": [[3, 2]], "far": None, "component_parts": None}, 0) == [{"run": "shell", "shell": 3, "high": 2}]
    assert run_plan.scatter_shells_pieces({"shells": [], "far": [20, 30]}, 0) == [{"run": "far", "first": 20, "last": 30}]


def test_a_phenomenon_scatter_unit_has_the_same_pieces_in_the_same_order():
    args = {"shells": [[3, 2], [4, -1]], "far": [20, 30], "component_parts": [["cluster", "A", 0, 1]]}
    assert run_plan.scatter_phenomena_pieces(args, 0) == [{"run": "shell", "shell": 3, "high": 2}, {"run": "shell", "shell": 4, "high": -1},
                                                         {"run": "far", "first": 20, "last": 30}, {"run": "components"}]
    assert run_plan.scatter_phenomena_pieces(args, 3) == [{"run": "components"}]


def test_a_nebula_field_unit_has_a_piece_for_each_slab_after_the_cursor():
    assert run_plan.draw_nebula_field_pieces({"i": -3, "count": 4}, 0) == [-3, -2, -1, 0]
    assert run_plan.draw_nebula_field_pieces({"i": -3, "count": 4}, 2) == [-1, 0]
    assert run_plan.draw_nebula_field_pieces({"i": 5, "count": 1}, 1) == []


@pytest.mark.parametrize("low, high", [(1, 10), (1, run_plan.MAP_LUMINOSITY_PIECE_IDS), (7, 3 * run_plan.MAP_LUMINOSITY_PIECE_IDS + 5), (100, 100)])
def test_the_map_luminosity_pieces_cover_every_id_of_the_unit_once(low, high):
    pieces = run_plan.fill_map_luminosity_pieces({"low": low, "high": high}, 0)
    assert pieces[0]["low"] == low and pieces[-1]["high"] == high
    assert all(a["high"] + 1 == b["low"] for a, b in zip(pieces, pieces[1:]))
    assert all(1 <= piece["high"] - piece["low"] + 1 <= run_plan.MAP_LUMINOSITY_PIECE_IDS for piece in pieces)
    assert all(piece["weight"] == piece["high"] - piece["low"] + 1 for piece in pieces)       # OPS.243: a piece weighs its ids
    assert run_plan.fill_map_luminosity_pieces({"low": low, "high": high}, 1) == pieces[1:]


def test_the_shell_scatter_result_is_the_pieces_added_up_in_order():
    one, two = brightStars.POPULATIONS[:2]
    first = {"counts": {one: 5, two: 1}, "types": [["O star", "O stars", 2], ["B star", "B stars", 4]], "shells": [[3, 6]], "systems": 2, "seconds": 1.5}
    second = {"counts": {one: 2, two: 0}, "types": [["B star", "B stars", 1], ["A star", "A stars", 1]], "shells": [[4, 2]], "systems": 1, "seconds": 0.5}
    args = {"speed": {"kind": "mass", "density": 1.0}}
    result = run_plan.scatter_shells_gather(args, [first, second], None)
    assert result["counts"][one] == 7 and result["counts"][two] == 1 and set(result["counts"]) == set(brightStars.POPULATIONS)
    assert result["types"] == [["O star", "O stars", 2], ["B star", "B stars", 5], ["A star", "A stars", 1]]
    assert result["shells"] == [[3, 6], [4, 2]] and result["systems"] == 3 and result["seconds"] == pytest.approx(2.0)
    assert result["samples"] == [{"kind": "mass", "density": 1.0, "seconds": pytest.approx(2.0), "systems": 8, "stars": 8}]
    assert result["row"] == {"objects": 8, "touched": [3, 4]}
    assert run_plan.scatter_shells_gather({}, [{**first, "counts": {}, "shells": [[3, 0]]}], None)["samples"] == []        # no star, no time per star (PERF.53)


def test_the_phenomenon_result_is_the_pieces_added_up_with_their_write_time():
    results = [{"counts": {"neutron star": 3}, "write_seconds": 0.25, "systems": 1, "seconds": 1.0},
               {"counts": {"neutron star": 1, "black hole": 2}, "write_seconds": 0.5, "systems": 0, "seconds": 2.0}]
    result = run_plan.scatter_phenomena_gather({"weight": 7.0, "key": 4}, results, None)
    assert result["counts"] == {"neutron star": 4, "black hole": 2} and result["systems"] == 1
    assert result["seconds"] == pytest.approx(3.0) and result["write_seconds"] == pytest.approx(0.75)
    assert [sample["kind"] for sample in result["samples"]] == ["phenomena", "phenomena-write"]
    assert result["row"] == {"objects": 6, "touched": [4]}
    assert run_plan.scatter_phenomena_gather({"weight": 7.0, "key": 4}, [{"counts": {}, "write_seconds": 0.0, "systems": 0, "seconds": 0.1}], None)["row"] == {"objects": 0, "touched": []}


def test_the_nebula_and_fill_results_keep_the_shape_the_followers_read():
    drawn = [{"written": 3, "counts": {"M": 3}, "shells": (1, 2), "slab": -1}, {"written": 0, "counts": {}, "shells": None, "slab": 0}]
    assert run_plan.draw_nebula_field_gather({"i": -1, "count": 2}, drawn, None) == {"written": 3, "slab": -1, "count": 2, "drawn": drawn}
    assert run_plan.fill_map_luminosity_gather({"low": 5, "high": 14}, [{"filled": 4}, {"filled": 2}], None) == {"filled": 6, "ids": 10}
