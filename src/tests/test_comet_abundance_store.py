# tests/test_comet_abundance_store.py

"""GEN.1775: the comet abundance tier round-trips through a stored galaxy (the plan's pins) and reaches the config of a sector's systems."""

import types

import pytest

from planetgen.db import store
from planetgen.galaxy import dials, shape_draw
from planetgen.generation import bright_stars
from planetgen.generation.config import SystemConfig

pytestmark = pytest.mark.fresh_database

SEED = bytes.fromhex("00112233445566778899AABBCCDDEEFF")


def _save(mysql_config, **given):
    args = types.SimpleNamespace(**{**{row.dest: None for row in dials.LEVERS}, **given})
    shape, drawn = shape_draw.plan_shape(args, SEED)
    store.save_galaxy_shape(shape, edge_pc=4.0, outer_shell_index=10, expected_system_count_at_density_1=1.0, config=mysql_config, galaxy_seed=SEED,
                            pinned_stellar_mass_solar=shape_draw.pinned_stellar_mass(args, drawn), pins=shape_draw.pins_of(args))


def _load(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        return store.get_galaxy_shape(conn)
    finally:
        conn.close()


def test_a_plan_with_no_choice_reads_back_science_and_a_chosen_tier_reads_back_as_chosen(mysql_config):
    _save(mysql_config)
    assert _load(mysql_config).shape.comet_tier == "science"
    _save(mysql_config, comet_abundance="abundant")
    info = _load(mysql_config)
    assert info.shape.comet_tier == "abundant" and dials.choice_pinned(info.pins, "comet_abundance") == "abundant"
    assert info.shape.abundance_tier == "abundant", "the Class M tier is its own dial: its box is ticked when no pin says otherwise (GEN.1832)"
    _save(mysql_config, class_m="science")
    info = _load(mysql_config)
    assert info.shape.comet_tier == "science" and info.shape.abundance_tier == "science"


def test_a_sectors_config_takes_the_tier_only_when_it_is_not_science(mysql_config, monkeypatch):
    monkeypatch.setattr(bright_stars, "pick_population", lambda densities, rng: "thin_disk")
    monkeypatch.setattr(bright_stars, "age_laws_of", lambda shape: {})
    for tier, expected in (("science", None), ("sparse", "sparse")):
        _save(mysql_config, comet_abundance=tier)
        context = types.SimpleNamespace(shape=_load(mysql_config).shape, densities={}, place_terms=None, profile=None)
        config = bright_stars.FillContext.apply_population(context, SystemConfig())
        assert config.COMET_TIER == expected
