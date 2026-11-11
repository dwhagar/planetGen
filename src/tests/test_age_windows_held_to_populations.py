"""
GEN.1525 (row 4 of the GEN.945 audit): the age windows of the default star-formation history and of free compact remnants held against the thin disk's population windows
(`tuning.STELLAR_POPULATION_AGE_RANGES_GY`, `THIN_DISK_OLD_WINDOW_GY`).
"""

import pytest

from planetgen import tuning


def _thin_disk_windows():
    ranges = tuning.STELLAR_POPULATION_AGE_RANGES_GY
    return [ranges["young"], ranges["intermediate"], tuning.THIN_DISK_OLD_WINDOW_GY]


def test_the_default_star_formation_window_is_the_three_thin_disk_populations_laid_end_to_end():
    windows = _thin_disk_windows()
    assert all(a[1] == b[0] for a, b in zip(windows, windows[1:]))      # no gap and no overlap between young, intermediate and the thin disk's old window
    assert (windows[0][0], windows[-1][1]) == tuning.STAR_FORMATION_AGE_RANGE_GY
    assert tuning.STAR_FORMATION_AGE_RANGE_GY[1] < tuning.UNIVERSE_AGE_GY


def test_free_compact_remnants_are_as_old_as_the_thin_disk_and_older_than_a_magnetar_lives():
    low, high = tuning.COMPACT_REMNANT_AGE_RANGE_GY
    assert high == tuning.STAR_FORMATION_AGE_RANGE_GY[1]
    assert low * 1.0e9 == pytest.approx(tuning.MAGNETAR_ACTIVE_LIFETIME_YR)      # 10 kyr: the youngest a free neutron star is, so none is a magnetar (GEN.278)
