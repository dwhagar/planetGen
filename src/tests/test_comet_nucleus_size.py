"""
GEN.1765: a star-bound comet's nucleus diameter follows the cumulative size law of the Jupiter family ([R4139], [R5296]) between the ends of
`BOUND_COMET_NUCLEUS_DIAMETER_RANGE_KM`, not a flat draw: small nuclei are the common ones.
"""

import math

import pytest

from planetgen import tuning
from planetgen.generation.comet import Comet
from planetgen.generation.config import SystemConfig
from planetgen.util import draw

COUNT = 1500


def _diameters():
    config = SystemConfig()
    out = []
    for index in range(COUNT):
        with draw.bound(draw.Stream(f"comet-size:{index}")):
            out.append(Comet(config, 1.0, orbit_type="parabolic").nucleus_diameter_km)
    return out


def _share_above(diameter_km):
    """The share of the draw above a diameter by the law cut to the ends: (D^-q - hi^-q) / (lo^-q - hi^-q)."""
    low, high = tuning.BOUND_COMET_NUCLEUS_DIAMETER_RANGE_KM
    q = tuning.COMET_NUCLEUS_CUMULATIVE_INDEX
    return (diameter_km ** -q - high ** -q) / (low ** -q - high ** -q)


def test_the_cumulative_index_is_the_one_the_papers_measure():
    assert tuning.COMET_NUCLEUS_CUMULATIVE_INDEX == pytest.approx(1.92, abs=0.2)          # [R4139] 1.92 +- 0.20; [R5296] about 1.9


def test_the_nuclei_stay_between_the_ends_and_most_are_small():
    sizes = _diameters()
    low, high = tuning.BOUND_COMET_NUCLEUS_DIAMETER_RANGE_KM
    assert all(low <= size <= high for size in sizes)
    assert sorted(sizes)[COUNT // 2] < math.sqrt(low * high)          # the median sits under the geometric middle; a flat draw puts it at the arithmetic one (11 km)
    assert sorted(sizes)[COUNT // 2] < 5.0


@pytest.mark.parametrize("diameter_km", [4.0, 6.0, 10.0, 15.0])
def test_the_share_of_nuclei_above_a_size_follows_the_law(diameter_km):
    sizes = _diameters()
    share = sum(1 for size in sizes if size > diameter_km) / COUNT
    assert share == pytest.approx(_share_above(diameter_km), abs=0.035)
