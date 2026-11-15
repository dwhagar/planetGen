# tests/test_abundance_tier.py

"""GEN.1704: the Class M abundance tier of a galaxy (docs/design/habitability/class-m-abundance-tiers.md): what a tier's dials are, and that Science moves nothing."""

from types import SimpleNamespace

import pytest

from planetgen import tuning
from planetgen.generation import abundance_tier
from planetgen.generation.config import SystemConfig


def _planet(tier, zone="e", body_type="t", is_moon=False):
    config = SystemConfig()
    config.ABUNDANCE_TIER = tier
    return SimpleNamespace(system_config=config, zone=zone, body_type=body_type, is_moon=is_moon)


def test_science_is_the_default_and_every_dial_is_off():
    config = SystemConfig()
    assert config.ABUNDANCE_TIER is None and abundance_tier.tier_of(config) == "science"
    assert "ABUNDANCE_TIER" not in vars(config)      # a class attribute: a config on science holds no field of it (the golden digests list every field)
    assert abundance_tier.dials_of(config) == abundance_tier.OFF == (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    assert not abundance_tier.active(config) and abundance_tier.request_share(config) == 0.0
    for tier in (None, "science"):
        planet = _planet(tier)
        assert abundance_tier.water_share(planet) == 0.0 and abundance_tier.bare_share(planet) == 0.0 and abundance_tier.dry_part(planet) == 0.0
        assert abundance_tier.quiet_share(planet) == 0.0 and not abundance_tier.quiet_interior(planet) and not abundance_tier.desert_sea(planet)
        assert abundance_tier.cycle_share(planet) == 0.0


def test_every_other_tier_has_a_row_of_six_shares_that_never_fall_from_one_tier_to_the_next():
    assert set(tuning.ABUNDANCE_TIER_DIALS) == set(abundance_tier.TIERS) - {"science"}
    assert set(abundance_tier.TARGET_SYSTEMS_PER_CLASS_M) == set(abundance_tier.TIERS)
    rows = [abundance_tier.dials_of(_planet(tier)) for tier in abundance_tier.TIERS]
    for lower, higher in zip(rows, rows[1:]):
        assert all(0.0 <= a <= b <= 1.0 for a, b in zip(lower, higher))
    assert all(len(row) == 6 and 0.0 < row.dry_part < 0.5 and row.cycle_share > 0.0 for row in rows[1:])        # desert worlds beside the temperate one, which stays the larger part (GEN.1721)
    targets = [abundance_tier.TARGET_SYSTEMS_PER_CLASS_M[tier] for tier in abundance_tier.TIERS[1:]]
    assert targets == sorted(targets, reverse=True)


def test_the_targets_run_from_one_in_500_to_one_in_10_in_even_log_steps():
    targets = [abundance_tier.TARGET_SYSTEMS_PER_CLASS_M[tier] for tier in abundance_tier.TIERS[1:]]
    assert targets[0] == 500 and targets[-1] == 10          # Boss, 9:28 and 9:30 PM EST 2026-11-14
    step = (targets[0] / targets[-1]) ** (1.0 / (len(targets) - 1))
    for lower, higher in zip(targets, targets[1:]):
        assert abs(lower / higher / step - 1.0) < 0.02       # whole numbers, so even to within 2 percent


def test_the_tier_is_not_saved_with_a_system_config():
    config = SystemConfig()
    config.ABUNDANCE_TIER = "abundant"
    assert "abundance_tier" not in config.to_dict()
    assert SystemConfig.from_dict(config.to_dict()).ABUNDANCE_TIER is None


def test_a_tier_that_does_not_exist_is_refused():
    with pytest.raises(ValueError, match="not one of science"):
        abundance_tier.tier_of(_planet("common"))


@pytest.mark.parametrize("changes", [{"zone": "h"}, {"zone": "c"}, {"zone": None}, {"body_type": "j"}, {"body_type": "a"}, {"is_moon": True}])
def test_the_water_and_envelope_dials_act_on_rocky_planets_in_the_ecosphere_only(changes):
    inside = _planet("abundant")
    assert abundance_tier.candidate(inside)
    assert abundance_tier.water_share(inside) == tuning.ABUNDANCE_TIER_DIALS["abundant"][0] > 0.0
    assert abundance_tier.bare_share(inside) == tuning.ABUNDANCE_TIER_DIALS["abundant"][1] > 0.0
    assert abundance_tier.dry_part(inside) == tuning.ABUNDANCE_TIER_DIALS["abundant"][3] > 0.0
    assert abundance_tier.quiet_share(inside) == tuning.ABUNDANCE_TIER_DIALS["abundant"][4] > 0.0
    outside = _planet("abundant", **changes)
    assert not abundance_tier.candidate(outside)
    assert abundance_tier.water_share(outside) == 0.0 and abundance_tier.bare_share(outside) == 0.0
    assert abundance_tier.quiet_share(outside) == 0.0


def test_a_desert_world_and_a_quiet_interior_are_read_from_what_the_tier_set_on_the_body():
    plain = _planet("abundant")
    assert not abundance_tier.desert_sea(plain) and not abundance_tier.quiet_interior(plain)
    plain._water_dialed = "low"
    assert not abundance_tier.desert_sea(plain)          # Earth's class is not a desert world
    plain._water_dialed, plain._quiet_interior = "dry", True
    assert abundance_tier.desert_sea(plain) and abundance_tier.quiet_interior(plain)


def test_the_desert_water_range_is_inside_the_dry_class_and_the_quiet_sink_is_inside_earths_published_range():
    from planetgen.physics import biogenic_oxygen, delivered_water
    low, high = tuning.ABUNDANCE_DESERT_WATER_RANGE
    dry_low, dry_high = delivered_water.bracket(delivered_water.DESERT_CLASS)
    assert dry_low <= low < high <= dry_high
    assert 1.0e12 / biogenic_oxygen.EARTH_REDUCTANT_MOL_YR <= tuning.ABUNDANCE_QUIET_SINK_FACTOR < 1.0        # registry row BO2: 1e12 to 5e12 mol a year

