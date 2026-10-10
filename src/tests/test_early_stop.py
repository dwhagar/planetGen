# tests/test_early_stop.py

"""The layer-walking scatter passes stop once enough layers in a row (in walk order) produced nothing."""

import pytest

from planetgen.db import store
from planetgen.generation import bright_stars as brightStars
from planetgen.generation import early_stop, run_plan
from planetgen.util import log
from planetgen import tuning
from tests.test_bright_star_scatter import _plan_args, _seed_galaxy


def test_the_default_is_a_hundred_layers():
    assert tuning.SCATTER_DRY_LAYERS == 100
    assert early_stop.DryStreak(500).limit == 100


def test_a_run_of_empty_layers_in_a_row_stops_the_walk():
    streak = early_stop.DryStreak(10, limit=3)
    for position, produced in enumerate([5, 0, 0]):
        streak.record(position, produced)
    assert not streak.stopped
    streak.record(3, 0)
    assert streak.stopped and streak.stopped_at == 4


def test_a_productive_layer_starts_the_count_again():
    streak = early_stop.DryStreak(10, limit=3)
    for position, produced in enumerate([0, 0, 1, 0, 0]):
        streak.record(position, produced)
    assert not streak.stopped


def test_layers_finishing_out_of_order_count_in_walk_order():
    streak = early_stop.DryStreak(10, limit=3)
    streak.record(3, 0)
    streak.record(2, 0)
    streak.record(1, 0)
    assert not streak.stopped   # layer 0 hasn't reported: it might be productive
    streak.record(0, 4)
    assert streak.stopped_at == 4


def test_a_limit_of_zero_never_stops():
    streak = early_stop.DryStreak(10, limit=0)
    for position in range(10):
        streak.record(position, 0)
    assert not streak.stopped


def _empty_outer_layers(monkeypatch):
    real = brightStars.scatter_layer
    monkeypatch.setattr(brightStars, "scatter_layer",
                        lambda shape, layer_index, *args, **kwargs:
                        iter(()) if layer_index != 0 else real(shape, layer_index, *args, **kwargs))


def _stars(config):
    conn = store.get_connection(config)
    try:
        return conn.execute("SELECT COUNT(*) AS n FROM bright_stars").fetchone()["n"]
    finally:
        conn.close()


def test_stopping_early_gives_the_rows_a_full_walk_gives(mysql_config, monkeypatch):
    _empty_outer_layers(monkeypatch)
    _seed_galaxy(mysql_config)
    monkeypatch.setattr(tuning, "SCATTER_DRY_LAYERS", 0)
    full = run_plan.scatter_bright_stars(_plan_args(mysql_config, "--bright-stars-only", "--workers", "1"))
    full_rows = _stars(mysql_config)

    monkeypatch.setattr(tuning, "SCATTER_DRY_LAYERS", 1)
    messages = []
    monkeypatch.setattr(log, "normal", lambda message, *args, **kwargs: messages.append(message))
    early = run_plan.scatter_bright_stars(_plan_args(mysql_config, "--bright-stars-only", "--workers", "1"))
    assert early["total"] == full["total"] and _stars(mysql_config) == full_rows > 0
    assert any("stopped after 1 layers in a row with nothing; 1 of 3 layers not walked" in m for m in messages)


def test_the_phenomena_pass_gives_the_same_rows_on_a_normal_galaxy(mysql_config, monkeypatch):
    """A galaxy of a few layers never has 100 empty ones in a row, so the default walks them all."""
    from tests import test_phenomenon_scatter as phen
    phen._seed_galaxy(mysql_config)
    monkeypatch.setattr(tuning, "SCATTER_DRY_LAYERS", 0)
    full = run_plan.scatter_phenomena(phen._plan_args(mysql_config, "--workers", "1"))
    monkeypatch.setattr(tuning, "SCATTER_DRY_LAYERS", 100)
    early = run_plan.scatter_phenomena(phen._plan_args(mysql_config, "--workers", "1"))
    assert early["total"] == full["total"] > 0


def test_the_phenomena_pass_stops_after_empty_layers_but_keeps_the_special_phenomena(mysql_config, monkeypatch):
    from tests import test_phenomenon_scatter as phen
    from planetgen.generation import phenomenon_scatter
    phen._seed_galaxy(mysql_config)
    monkeypatch.setattr(phenomenon_scatter, "scatter_layer", lambda *args, **kwargs: iter(()))
    monkeypatch.setattr(tuning, "SCATTER_DRY_LAYERS", 1)
    messages = []
    monkeypatch.setattr(log, "normal", lambda message, *args, **kwargs: messages.append(message))
    summary = run_plan.scatter_phenomena(phen._plan_args(mysql_config, "--workers", "1"))
    assert any("Phenomena: stopped after 1 layers in a row with nothing; 2 of 3 layers not walked" in m
               for m in messages)
    assert summary["total"] > 0   # the nucleus and hypervelocity stars are not layers
