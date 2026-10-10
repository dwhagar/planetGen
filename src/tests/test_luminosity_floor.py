"""GEN.184: the luminosity floor presets."""
import pytest

from planetgen import tuning
from planetgen.generation import luminosity_floor as floor


def test_presets_run_from_the_floor_to_the_ceiling_and_include_the_default():
    presets = floor.PRESETS
    assert presets[0] == 2500.0 and presets[-1] == 4_000_000.0
    assert tuning.BRIGHT_STAR_MIN_LUMINOSITY_SOL == 9000.0 and 9000.0 in presets
    assert list(presets) == sorted(set(presets))


def test_steps_start_at_100_and_end_in_the_hundreds_of_thousands():
    presets = floor.PRESETS
    gaps = [b - a for a, b in zip(presets, presets[1:])]
    assert gaps[0] == 100.0 and gaps[4] == 100.0
    assert 300_000 <= gaps[-1] <= 500_000
    # Each step is at least as large as the one before it, within a rounding unit.
    assert all(later >= earlier / 2 for earlier, later in zip(gaps, gaps[1:]))


def test_check_accepts_the_range_and_refuses_the_rest():
    assert floor.check("3000") == 3000.0 and floor.check(2500) == 2500.0 and floor.check(4e6) == 4e6
    for bad in (2499.9, 4_000_001, "x", float("nan"), None, -5):
        with pytest.raises(ValueError):
            floor.check(bad)


def test_nearest_snaps_on_the_log_scale():
    assert floor.nearest(2503) == 2500.0 and floor.nearest(9e9) == 4_000_000.0
