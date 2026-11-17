# tests/test_comet_abundance.py

"""GEN.1775: the comet abundance tier of a galaxy (`generation/comet_abundance.py`): Science is the physics (comets only where a star's own outer disk was scattered), each tier scales the registered stand-in chance of the
stars without one, and the choice goes through the plan, the stored galaxy and a sector's config."""

import types

import pytest

from planetgen import tuning
from planetgen.galaxy import dials, shape_draw
from planetgen.generation import comet_abundance
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import comet_supply
from planetgen.util import draw
from tests.test_systems import make_config

SEED = bytes.fromhex("00112233445566778899AABBCCDDEEFF")
WEIGHTS = {"jupiter_family": 0.5, "halley_type": 0.3, "long_period": 0.2}


def _config(tier=None):
    config = SystemConfig()
    config.COMET_TIER = tier
    return config


def test_science_is_the_default_and_no_factor_acts_on_it():
    assert comet_abundance.DEFAULT == "science" == comet_abundance.TIERS[0]
    assert comet_abundance.tier_of(SystemConfig()) == "science" and comet_abundance.multiplier(SystemConfig()) == 0.0
    assert comet_abundance.tier_of(_config("science")) == "science"


def test_every_other_tier_has_a_factor_that_rises_and_plentiful_is_the_old_chance():
    factors = [tuning.COMET_TIER_MULTIPLIERS[tier] for tier in comet_abundance.TIERS[1:]]
    assert factors == sorted(factors) and len(set(factors)) == 4 and set(tuning.COMET_TIER_MULTIPLIERS) == set(comet_abundance.TIERS[1:])
    assert tuning.COMET_TIER_MULTIPLIERS["plentiful"] == 1.0, "the tier that keeps the chance the generator held before the tier"


def test_a_tier_that_does_not_exist_is_refused_and_the_tier_is_not_saved_with_the_config():
    with pytest.raises(ValueError, match="not one of science"):
        comet_abundance.tier_of(_config("common"))
    config = _config("abundant")
    assert "comet_tier" not in config.to_dict() and SystemConfig.from_dict(config.to_dict()).COMET_TIER is None


def test_a_star_the_physics_gives_comets_has_them_on_every_tier():
    for tier in comet_abundance.TIERS:
        assert comet_abundance.presence_chance(_config(tier), True) == 1.0


def test_science_is_the_physics_and_each_tier_scales_the_stand_in_for_the_stars_without_a_scattered_disk():
    assert comet_abundance.presence_chance(_config("science"), False) == 0.0
    chances = [comet_abundance.presence_chance(_config(tier), False) for tier in comet_abundance.TIERS]
    assert chances == sorted(chances) and chances[0] == 0.0
    for tier in comet_abundance.TIERS[1:]:
        assert comet_abundance.presence_chance(_config(tier), False) == pytest.approx(min(1.0, tuning.SYSTEM_COMET_CHANCE * tuning.COMET_TIER_MULTIPLIERS[tier]))


def _share_with_comets(tier, reservoir, monkeypatch, trials=300):
    monkeypatch.setattr(comet_supply, "family_weights", lambda belts: WEIGHTS if reservoir else None)
    system = StarSystem(system_config=make_config("G2V", PLANETS=False))
    system.system_config.COMET_TIER = tier
    hits = 0
    for seed in range(trials):
        with draw.bound(f"comet-abundance:{seed}"):
            hits += bool(system._generate_comets(system.star, None, []))
    return hits / trials


def test_in_the_generator_science_gives_no_comet_without_a_disk_and_every_one_with_it(monkeypatch):
    assert _share_with_comets("science", False, monkeypatch) == 0.0
    assert _share_with_comets("science", True, monkeypatch) == 1.0
    assert _share_with_comets("abundant", True, monkeypatch) == 1.0


def test_in_the_generator_each_tier_gives_the_stars_without_a_disk_its_share(monkeypatch):
    shares = {tier: _share_with_comets(tier, False, monkeypatch) for tier in comet_abundance.TIERS}
    assert shares["science"] == 0.0 and list(shares.values()) == sorted(shares.values())
    for tier in comet_abundance.TIERS[1:]:
        assert shares[tier] == pytest.approx(comet_abundance.presence_chance(_config(tier), False), abs=0.08)


def test_a_forced_request_still_gives_comets_on_science(monkeypatch):
    monkeypatch.setattr(comet_supply, "family_weights", lambda belts: None)
    system = StarSystem(system_config=make_config("G2V", PLANETS=False, COMETS=True))
    with draw.bound("comet-abundance-forced"):
        assert system._generate_comets(system.star, None, [])


def test_the_lever_is_a_census_dial_with_the_same_choices_as_the_class_m_abundance():
    row = dials.lever("comet_abundance")
    assert row.kind == "enum" and row.dial == "census" and row.choices == dials.lever("class_m").choices == comet_abundance.TIERS
    assert row.flag == "--comet-abundance" and row.parse("abundant") == "abundant"


def test_the_plan_pins_the_choice_and_a_sector_config_takes_it_only_off_science():
    args = types.SimpleNamespace(**{row.dest: None for row in dials.LEVERS})
    assert "comet_abundance" not in shape_draw.pins_of(args)
    args.comet_abundance = "plentiful"
    pins = shape_draw.pins_of(args)
    assert dials.choice_pinned(pins, "comet_abundance") == "plentiful"
    shape, _drawn = shape_draw.plan_shape(args, SEED)
    assert shape.comet_tier == "science", "the shape carries the tier only once the stored pins are read back"
