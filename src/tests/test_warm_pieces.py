# tests/test_warm_pieces.py
"""OPS.205: the warm commands and `end-map-rebuild-maintenance` give the engine `pieces` and `do` (`queue/engine.py`) in place of a `unit` function that loops for itself."""

import pytest

from planetgen.galaxy import field_job
from planetgen.queue import commands, engine, unitbeat
from planetgen.web import warmup

CONVERTED = ("warm-cluster-field", "warm-object-index", "warm-galaxy-map", "warm-galaxy-nodes", "warm-galaxy-deep-nodes", "warm-lists", "warm-sectors", "warm-systems",
             "warm-phenomena", "warm-search", "end-map-rebuild-maintenance")


@pytest.mark.parametrize("kind", CONVERTED)
def test_every_warm_command_runs_through_the_engines_loop(kind):
    spec = commands.spec(kind)
    assert spec.unit is None and spec.pieces is not None and spec.do is not None and spec.gather is not None
    assert callable(commands.pieces_function(kind)) and callable(commands.do_function(kind)) and callable(commands.gather_function(kind))


def test_an_object_slice_is_cut_in_pieces_of_fifty_and_the_first_carries_what_the_cap_left_out():
    args = {"database": "x", "cap": None, "offset": 1000, "count": 130, "left_out": 7}
    pieces = warmup.pieces_warm_systems(args, 0)
    assert [(p["offset"], p["count"], p["left_out"]) for p in pieces] == [(1000, 50, 7), (1050, 50, 0), (1100, 30, 0)]
    assert warmup.pieces_warm_systems(args, 2) == pieces[2:]
    assert sum(p["count"] for p in pieces) == args["count"]


def test_an_empty_slice_and_a_unit_with_no_slice_are_one_piece_each():
    assert warmup.pieces_warm_sectors({"database": "x", "offset": 0, "count": 0, "left_out": 4}, 0) == [{"offset": 0, "count": 0, "left_out": 4}]
    assert warmup.pieces_warm_phenomena({"database": "x", "cap": None}, 0) == [{}]
    assert warmup.pieces_warm_phenomena({"database": "x", "cap": None}, 1) == []


def test_node_keys_are_cut_in_chunks_of_the_most_one_request_holds():
    keys = [[0, i, 0] for i in range(warmup.NODE_CHUNK * 2 + 2)]
    chunks = warmup.pieces_warm_galaxy_nodes({"database": "x", "keys": keys}, 0)
    assert [len(chunk) for chunk in chunks] == [warmup.NODE_CHUNK, warmup.NODE_CHUNK, 2]
    assert warmup.pieces_warm_galaxy_deep_nodes({"database": "x", "keys": keys}, 1) == chunks[1:]
    assert warmup.pieces_warm_galaxy_nodes({"database": "x", "keys": []}, 0) == []


def test_the_object_index_is_one_piece_per_kind():
    from planetgen.db import objectindex
    assert warmup.pieces_warm_object_index({"database": "x"}, 0) == list(objectindex.KINDS)
    assert warmup.pieces_warm_object_index({"database": "x"}, 1) == list(objectindex.KINDS)[1:]


def test_a_page_family_counts_what_it_visits_as_the_units_second_level(monkeypatch):
    steps = []
    monkeypatch.setattr(unitbeat, "step", lambda name, total=None: steps.append((name, total)))
    advanced = []
    monkeypatch.setattr(unitbeat, "advance", lambda amount=1: advanced.append(amount))

    def warm_family(db, family, cap=None, app=None, config=None, window=None, report=None):
        report(0, 3)
        report(2, 3)
        report(3, 3)
        return {"family": family, "warmed": 3, "noun": "pages", "cached": 1, "built": 2, "left_out": 0, "seconds": 0.5, "line": f"{family}: ok"}

    monkeypatch.setattr(warmup, "warm_family", warm_family)
    beat = []
    result = commands.run_alone("warm-lists", {"database": "x", "cap": None}, lambda done, estimate: beat.append((done, estimate)))
    assert result["family"] == "lists" and result["warmed"] == 3
    assert beat == [(0, 1), (1, 1)]                                       # the engine's count: one piece of one
    assert steps == [("Warming the lists", 3)] and sum(advanced) == 3     # the second level: the pages of the lists, from the family's own report


def test_an_object_unit_warms_each_part_of_its_slice_and_adds_the_parts_up(monkeypatch):
    seen = []

    def warm_family(db, family, cap=None, app=None, config=None, window=None, report=None):
        seen.append(dict(window))
        return {"family": family, "warmed": window["count"], "noun": "systems", "cached": 1, "built": window["count"] - 1, "left_out": window["left_out"],
                "chunk": window["count"], "seconds": 1.0, "line": "x"}

    monkeypatch.setattr(warmup, "warm_family", warm_family)
    args = {"database": "x", "cap": None, "offset": 0, "count": 120, "left_out": 9}
    result = commands.run_alone("warm-systems", args)
    assert [(w["offset"], w["count"], w["left_out"]) for w in seen] == [(0, 50, 9), (50, 50, 0), (100, 20, 0)]
    assert (result["warmed"], result["cached"], result["built"], result["left_out"], result["seconds"]) == (120, 3, 117, 9, 3.0)
    assert result["line"].startswith("systems: 120 systems warmed") and "9 left out by the cap" in result["line"]


class _Ctx:
    def __init__(self):
        self.stepped = []

    def step(self, name, total=None):
        self.stepped.append((name, total))

    def advance(self, amount=1):
        pass


def _nodes(monkeypatch, full_after=None, cached=0):
    """The node pieces with the site and the cache faked: `full_after` chunks are built, then the cache reads full."""
    built = []

    def is_full():
        return full_after is not None and len(built) >= full_after

    class _App:
        def test_request_context(self, path):
            import contextlib
            return contextlib.nullcontext()

    monkeypatch.setattr(warmup.tilecache, "cache_is_full", is_full)
    monkeypatch.setattr(warmup.tilecache, "max_cache_bytes", lambda: 100 * 1048576)
    monkeypatch.setattr(warmup, "app_for", lambda config=None: _App())
    monkeypatch.setattr(warmup, "fetch_nodes", lambda db, keys, ahead=False: built.append(list(keys)) or {"cached": min(cached, len(keys))})
    return built


def test_a_node_unit_builds_every_chunk_and_says_what_was_stored(monkeypatch):
    built = _nodes(monkeypatch, cached=10)
    keys = [[0, i, 0] for i in range(warmup.NODE_CHUNK + 5)]
    result = commands.run_alone("warm-galaxy-nodes", {"database": "x", "keys": keys, "stored": 3})
    assert [len(chunk) for chunk in built] == [warmup.NODE_CHUNK, 5]
    assert (result["warmed"], result["left_out"], result["cached"], result["built"]) == (len(keys) + 3, 0, 10 + 5 + 3, len(keys) - 15)
    assert result["family"] == "galaxy-nodes" and result["line"].startswith("galaxy-nodes: ")


def test_a_node_unit_leaves_out_what_the_cache_has_no_room_for(monkeypatch):
    built = _nodes(monkeypatch, full_after=1)
    keys = [[0, i, 0] for i in range(warmup.NODE_CHUNK * 3)]
    result = commands.run_alone("warm-galaxy-deep-nodes", {"database": "x", "keys": keys})
    assert len(built) == 1 and result["left_out"] == warmup.NODE_CHUNK * 2 and result["warmed"] == warmup.NODE_CHUNK
    assert "left out as the cache filled" in result["line"] and result["family"] == "galaxy-deep-nodes"


def test_a_node_unit_that_finds_the_cache_full_at_the_start_leaves_all_out(monkeypatch):
    built = _nodes(monkeypatch, full_after=0)
    result = commands.run_alone("warm-galaxy-nodes", {"database": "x", "keys": [[0, 1, 0], [0, 2, 0]], "stored": 4})
    assert built == [] and (result["left_out"], result["warmed"], result["built"]) == (2, 4, 0)
    assert "the cache is full" in result["line"]


def test_a_node_unit_with_no_key_says_what_the_stage_found_stored(monkeypatch):
    _nodes(monkeypatch)
    assert commands.run_alone("warm-galaxy-nodes", {"database": "x", "keys": [], "stored": 12})["line"].startswith("galaxy-nodes: 12 nodes warmed (0 built, 12 already stored)")
    assert "no node to warm, no block is stored yet" in commands.run_alone("warm-galaxy-nodes", {"database": "x", "keys": []})["line"]
    assert commands.run_alone("warm-galaxy-deep-nodes", {"database": "x", "keys": []})["line"] == "galaxy-deep-nodes: no node to warm"


def test_the_map_rebuild_reopens_the_site_in_one_piece(monkeypatch):
    from planetgen.admin import maintenance
    calls = []
    monkeypatch.setattr(maintenance, "turn_off_for", lambda cause, by: calls.append((cause, by)) or True)
    assert commands.run_alone("end-map-rebuild-maintenance", {}) == {"reopened": True, "line": "Maintenance mode is off: the Galaxy Map is rebuilt."}
    assert calls == [(maintenance.SECTOR_SWITCH, "map rebuild")]
    monkeypatch.setattr(maintenance, "turn_off_for", lambda cause, by: False)
    assert commands.run_alone("end-map-rebuild-maintenance", {})["reopened"] is False


def test_a_cluster_field_unit_with_nothing_to_make_has_no_piece():
    assert field_job.pieces_field({"database": "x", "state": "empty"}, 0) == []
    assert field_job.pieces_field({"database": "x", "state": "ready"}, 0) == []
    assert field_job.gather_field({"database": "x", "state": "empty"}, [], None)["line"].startswith("cluster-field: ")


def test_the_engine_counts_the_pieces_of_a_node_unit(monkeypatch):
    _nodes(monkeypatch)
    beat = []
    keys = [[0, i, 0] for i in range(warmup.NODE_CHUNK + 1)]
    engine.run("warm-galaxy-nodes", {"database": "x", "keys": keys}, lambda done, total: beat.append((done, total)))
    assert beat == [(0, 2), (1, 2), (2, 2)]
