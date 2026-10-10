# tests/test_layer_groups.py

"""PERF.57: after 5 empty layers in a row the scatter passes draw the next 10 as one group, doubling while it is empty."""

import statistics

import pytest

from planetgen.generation import bright_stars as brightStars
from planetgen.generation import layer_groups, phenomenon_scatter, run_plan
from planetgen.util import log
from planetgen import tuning
from tests.test_bright_star_scatter import EDGE_PC, E_VALUE, SHAPE, THRESHOLD, _plan_args, _seed_galaxy


class _InlineQueue:
    def wait_until(self, condition):
        assert condition()


def _walk(results, **kwargs):
    """Drives a `Walk` over layers that take `results[position]` objects when drawn singly (groups take the sum
    of theirs); returns the order of the units it drew."""
    walk = layer_groups.Walk([0.0] * len(results), **kwargs)
    drawn = []

    def single(position):
        drawn.append(("single", position))
        walk.record(position, results[position])

    def group(start, size):
        drawn.append(("group", start, size))
        walk.record_group(start, size, sum(results[start:start + size]))

    walk.drive(_InlineQueue(), single, group)
    return walk, drawn


def test_five_empty_layers_in_a_row_start_a_group_of_ten():
    _walk_obj, drawn = _walk([3, 0, 0, 0, 0, 0] + [0] * 6 + [0] * 4)
    assert drawn[:6] == [("single", position) for position in range(6)]
    assert drawn[6] == ("group", 6, 10)


def test_a_group_that_takes_nothing_is_followed_by_one_twice_the_size():
    _walk_obj, drawn = _walk([0] * 100)
    sizes = [unit[2] for unit in drawn if unit[0] == "group"]
    assert sizes == [10, 20, 40, 25]   # the last one is what is left of the 100 layers
    assert [unit for unit in drawn if unit[0] == "single"] == [("single", position) for position in range(5)]


def test_a_placement_takes_the_walk_back_to_single_layers_and_the_first_group_size():
    results = [0] * 5 + [0] * 10 + [0] * 20      # five empty, an empty group of 10, then a group of 20 that places
    results[20] = 4                              # inside the group of 20
    results += [0] * 10
    walk, drawn = _walk(results)
    assert drawn[5] == ("group", 5, 10) and drawn[6] == ("group", 15, 20)
    assert drawn[7] == ("single", 35)            # single layers again
    assert walk.group_layers == 10 or walk.group_sizes[:2] == [10, 20]


def test_nothing_is_grouped_when_the_limit_is_zero():
    _walk_obj, drawn = _walk([0] * 30, empty_before=0)
    assert all(unit[0] == "single" for unit in drawn) and len(drawn) == 30


def test_a_layer_is_not_started_while_results_to_come_could_put_the_walk_in_a_group():
    walk = layer_groups.Walk([0.0] * 20, certain=30.0)
    for position in range(2):
        walk.record(position, 0)
    # Layers 2 to 4 are out: if they were all empty the walk would be in a group by layer 5.
    assert not walk.must_wait(4) and walk.must_wait(5)
    walk.record(4, 1)                            # one of them took something, so layer 5 is safe
    assert not walk.must_wait(5)


def test_a_layer_expected_to_hold_many_is_counted_as_certain_to_take_some():
    walk = layer_groups.Walk([0.0, 0.0, 0.0, 500.0, 0.0, 0.0, 0.0, 0.0], certain=30.0)
    assert not walk.must_wait(7)


# -- the group draw -------------------------------------------------------------------------------

LAYERS = [(0, 8), (-1, 6), (1, 6)]


def test_a_group_places_stars_in_its_layers_with_about_the_layers_expected_count():
    fractions = brightStars.band_fractions(THRESHOLD)
    expected = sum(brightStars.layer_expected_stars(SHAPE, layer, ring, EDGE_PC, E_VALUE, fractions)
                   for layer, ring in LAYERS)
    counts = [len(list(brightStars.scatter_group(SHAPE, LAYERS, EDGE_PC, E_VALUE, THRESHOLD, seed)))
              for seed in range(12)]
    per_layer = [sum(len(list(brightStars.scatter_layer(SHAPE, layer, ring, EDGE_PC, E_VALUE, THRESHOLD, seed)))
                     for layer, ring in LAYERS) for seed in range(12)]
    assert expected > 20
    assert statistics.mean(counts) == pytest.approx(statistics.mean(per_layer), rel=0.15)
    rows = list(brightStars.scatter_group(SHAPE, LAYERS, EDGE_PC, E_VALUE, THRESHOLD, 3))
    assert {row[1] for row in rows} <= {layer for layer, _ring in LAYERS}


def test_a_group_places_more_stars_in_the_denser_layers():
    rows = [row for seed in range(10) for row in brightStars.scatter_group(SHAPE, LAYERS, EDGE_PC, E_VALUE,
                                                                           THRESHOLD, seed)]
    by_layer = {layer: sum(1 for row in rows if row[1] == layer) for layer, _ring in LAYERS}
    assert by_layer[0] > by_layer[1] and by_layer[0] > by_layer[-1]


def test_the_same_group_gives_the_same_stars():
    first = list(brightStars.scatter_group(SHAPE, LAYERS, EDGE_PC, E_VALUE, THRESHOLD, 5))
    again = list(brightStars.scatter_group(SHAPE, LAYERS, EDGE_PC, E_VALUE, THRESHOLD, 5))
    assert first == again and first


def test_a_group_leaves_filled_sectors_out():
    rows = list(brightStars.scatter_group(SHAPE, LAYERS, EDGE_PC, E_VALUE, THRESHOLD, 5))
    skip = {(row[0], row[1], row[2]) for row in rows}
    assert list(brightStars.scatter_group(SHAPE, LAYERS, EDGE_PC, E_VALUE, THRESHOLD, 5, skip_addresses=skip)) == []


def test_a_phenomena_group_places_about_the_layers_expected_count():
    expected = sum(phenomenon_scatter.layer_expected(SHAPE, layer, ring, EDGE_PC, E_VALUE, 5.0)
                   for layer, ring in LAYERS)
    counts = [len(list(phenomenon_scatter.scatter_group(SHAPE, LAYERS, EDGE_PC, E_VALUE, seed, min_mass_solar=5.0)))
              for seed in range(25)]
    per_layer = [sum(len(list(phenomenon_scatter.scatter_layer(SHAPE, layer, ring, EDGE_PC, E_VALUE, seed,
                                                                 min_mass_solar=5.0)))
                     for layer, ring in LAYERS) for seed in range(25)]
    assert expected > 1
    assert statistics.mean(counts) == pytest.approx(statistics.mean(per_layer), rel=0.3)


# -- the passes ---------------------------------------------------------------------------------

def _empty_layers_but_the_middle(monkeypatch):
    real = brightStars.scatter_layer
    monkeypatch.setattr(brightStars, "scatter_layer",
                        lambda shape, layer_index, *args, **kwargs:
                        real(shape, layer_index, *args, **kwargs) if layer_index == 0 else iter(()))


def test_the_star_pass_draws_empty_layers_as_groups_and_says_so(mysql_config, monkeypatch):
    _seed_galaxy(mysql_config)
    # The toy galaxy has three layers, so shrink the trigger to make a group of what is left.
    monkeypatch.setattr(tuning, "SCATTER_EMPTY_LAYERS_BEFORE_GROUP", 1)
    monkeypatch.setattr(tuning, "SCATTER_GROUP_LAYERS", 10)
    _empty_layers_but_the_middle(monkeypatch)
    messages = []
    monkeypatch.setattr(log, "normal", lambda message, *args, **kwargs: messages.append(message))
    run_plan.scatter_bright_stars(_plan_args(mysql_config, "--bright-stars-only", "--workers", "1"))
    # Layer 0 takes stars; layer -1 is empty, so layer 1 is drawn as a group of what is left (one layer).
    assert any("layers drawn in 1 groups (sizes 1)" in message for message in messages)


def test_the_phenomena_pass_draws_empty_layers_as_groups_and_keeps_the_special_phenomena(mysql_config, monkeypatch):
    from tests import test_phenomenon_scatter as phen
    phen._seed_galaxy(mysql_config)
    monkeypatch.setattr(tuning, "SCATTER_EMPTY_LAYERS_BEFORE_GROUP", 1)
    monkeypatch.setattr(phenomenon_scatter, "scatter_layer", lambda *args, **kwargs: iter(()))
    messages = []
    monkeypatch.setattr(log, "normal", lambda message, *args, **kwargs: messages.append(message))
    summary = run_plan.scatter_phenomena(phen._plan_args(mysql_config, "--workers", "1"))
    assert any(message.startswith("Phenomena: ") and "layers drawn in" in message and "groups" in message
               for message in messages)
    assert summary["total"] > 0   # the group places some, and the nucleus and hypervelocity stars are not layers


def test_a_normal_galaxy_is_drawn_layer_by_layer(mysql_config, monkeypatch):
    """A galaxy of a few layers never has five empty ones in a row: no groups."""
    from tests import test_phenomenon_scatter as phen
    phen._seed_galaxy(mysql_config)
    messages = []
    monkeypatch.setattr(log, "normal", lambda message, *args, **kwargs: messages.append(message))
    run_plan.scatter_phenomena(phen._plan_args(mysql_config, "--workers", "1"))
    assert not any("drawn in" in message and "groups" in message for message in messages)


WIDE_EXTENTS = [(layer, 6) for layer in range(-40, 41)]
"""A toy galaxy of 81 layers, nearly all of them empty, so the walk groups the outer ones."""


def _stars_with(mysql_config, workers):
    from planetgen.db import store
    store.save_galaxy_shape(SHAPE, edge_pc=EDGE_PC, outer_ring_index=8, expected_system_count_at_density_1=E_VALUE,
                            config=mysql_config)
    store.replace_galaxy_layers(WIDE_EXTENTS, config=mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        conn.execute("DELETE FROM bright_stars")
        conn.commit()
    finally:
        conn.close()
    messages = []
    real = log.normal
    log.normal = lambda message, *args, **kwargs: messages.append(message)
    try:
        run_plan.scatter_bright_stars(_plan_args(mysql_config, "--bright-stars-only", "--workers", str(workers),
                                                 "--bright-star-min-luminosity", "400000"))
    finally:
        log.normal = real
    conn = store.get_connection(mysql_config)
    try:
        rows = conn.execute("SELECT ring_index, layer_index, ring_slot_index, position_x_mpc, position_y_mpc,"
                            " position_z_mpc FROM bright_stars ORDER BY layer_index, ring_index, ring_slot_index,"
                            " position_x_mpc").fetchall()
    finally:
        conn.close()
    return [tuple(row.values()) if hasattr(row, "values") else tuple(row) for row in rows], messages


def test_any_number_of_workers_draws_the_same_stars_with_the_same_groups(mysql_config):
    one, one_log = _stars_with(mysql_config, 1)
    two, two_log = _stars_with(mysql_config, 2)
    assert any("groups" in message for message in one_log), one_log
    assert one == two and one
    assert [m for m in one_log if "groups" in m] == [m for m in two_log if "groups" in m]
