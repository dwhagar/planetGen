# tests/test_rogue_density.py

"""GEN.1115: the smoothed star density a rogue planet is drawn from (`galaxy/rogue_density.py`; rogue-planet-drift.md sections 3.3 and 5)."""

import math
import random

import numpy as np
import pytest

from planetgen.galaxy import rogue_density as rd
from planetgen.physics import rogue_drift

SCALE_LENGTH, SCALE_HEIGHT = 2600.0, 300.0


def disc(x, y, z):
    """An exponential disc with a sech^2 height: no azimuth dependence."""
    r = np.hypot(x, y)
    return np.exp(-r / SCALE_LENGTH) / np.cosh(z / (2.0 * SCALE_HEIGHT)) ** 2


def spiral(x, y, z, arms=2, pitch=math.radians(15.0), amplitude=0.4):
    """The disc with logarithmic spiral arms: a cosine in the arm angle."""
    r = np.maximum(np.hypot(x, y), 1.0)
    phase = arms * (np.arctan2(y, x) - np.log(r / 8000.0) / math.tan(pitch))
    return disc(x, y, z) * (1.0 + amplitude * np.cos(phase))


def test_the_ring_average_removes_arms_and_bars_exactly():
    for arms in (1, 2, 3, 4, 7):
        mean = rd.ring_average(lambda x, y, z: spiral(x, y, z, arms=arms), 8000.0, 0.0)
        assert mean == pytest.approx(float(disc(np.array([8000.0]), np.array([0.0]), np.array([0.0]))[0]), rel=1e-9)


def test_the_ring_average_of_a_disc_is_its_value():
    assert rd.ring_average(disc, 5000.0, 120.0) == pytest.approx(math.exp(-5000.0 / SCALE_LENGTH) / math.cosh(120.0 / 600.0) ** 2, rel=1e-12)


def test_smoothing_leaves_a_smooth_density_nearly_unchanged_and_a_linear_one_exactly():
    uniform = rd.old_density(lambda x, y, z: np.full_like(x, 2.5), 8000.0, 10.0)
    assert uniform == pytest.approx(2.5)
    linear = rd.old_density(lambda x, y, z: 1.0 + 1e-3 * np.hypot(x, y) + 2e-3 * z, 8000.0, 15.0)
    assert linear == pytest.approx(1.0 + 1e-3 * 8000.0 + 2e-3 * 15.0, rel=1e-12)
    smooth = rd.old_density(disc, 8000.0, 0.0)
    assert smooth == pytest.approx(float(disc(np.array([8000.0]), np.array([0.0]), np.array([0.0]))[0]), rel=0.01)


def test_the_smoothing_agrees_with_a_brute_force_average_over_the_gaussians():
    rng = random.Random(3)
    total, count = 0.0, 40000
    for _ in range(count):
        r = abs(5000.0 + rd.RADIAL_SMOOTH_PC * rng.gauss(0.0, 1.0))
        z = 100.0 + rd.VERTICAL_SMOOTH_PC * rng.gauss(0.0, 1.0)
        total += float(disc(np.array([r]), np.array([0.0]), np.array([z]))[0])
    assert rd.old_density(disc, 5000.0, 100.0) == pytest.approx(total / count, rel=0.01)


def test_an_old_rogue_does_not_see_the_arm_the_sector_is_in():
    """Two sectors on the same ring, one on an arm crest and one between arms: the old density is the same, the sector's own is not."""
    r = 8000.0
    for turn in np.linspace(0.0, math.pi / 2, 9):
        position = (r * math.cos(turn), r * math.sin(turn), 0.0)
        assert rd.old_density(spiral, math.hypot(position[0], position[1]), 0.0) == pytest.approx(rd.old_density(spiral, r, 0.0), rel=1e-12)
    values = [float(spiral(np.array([r * math.cos(t)]), np.array([r * math.sin(t)]), np.array([0.0]))[0]) for t in np.linspace(0.0, math.pi, 25)]
    assert max(values) > 1.5 * min(values)


def test_the_young_share_follows_the_donor_and_the_rest_is_smooth():
    position = (8000.0, 0.0, 0.0)
    local = float(spiral(np.array([8000.0]), np.array([0.0]), np.array([0.0]))[0])
    smooth = rd.old_density(spiral, 8000.0, 0.0)
    share = rd.young_share(10000.0)
    assert share == pytest.approx(0.01)
    assert rd.rogue_density(spiral, position, 10000.0) == pytest.approx(share * local + (1.0 - share) * smooth, rel=1e-12)
    assert rd.young_share(50.0) == 1.0
    assert rd.rogue_density(spiral, position, 50.0) == pytest.approx(local, rel=1e-12)


def test_the_sectors_and_the_map_read_the_same_function():
    """The map reads `old_density(R, z)`; the sector's mean with no young share is that number at the sector's radius and height."""
    position = (3000.0, 4000.0, 25.0)
    assert rd.rogue_density(spiral, position, 1.0e12) == pytest.approx(rd.old_density(spiral, 5000.0, 25.0), rel=1e-6)


def test_a_sector_on_the_axis_is_reflected_not_negative():
    assert rd.old_density(disc, 0.0, 0.0) > 0.0
    assert rd.old_density(disc, 0.0, 0.0) <= 1.0


def test_no_smoothing_is_the_ring_average():
    assert rd.old_density(spiral, 8000.0, 40.0, 0.0, 0.0) == pytest.approx(rd.ring_average(spiral, 8000.0, 40.0), rel=1e-12)


def test_bad_input_is_refused():
    with pytest.raises(ValueError):
        rd.old_density(disc, -1.0, 0.0)
    with pytest.raises(ValueError):
        rd.old_density(disc, 100.0, 0.0, -1.0, 0.0)
    with pytest.raises(ValueError):
        rd.young_share(0.0)


def test_the_smoothing_widths_are_the_drift_modules_rms_offsets():
    """The stand-in kick distribution (median 3 km/s, ln-sigma 0.8, isotropic) at 1 to 10 Gyr, through `rogue_drift.drift`: the rms radial and vertical offsets are the module's widths."""
    rng = random.Random(1)
    xs, zs = [], []
    for _ in range(20000):
        speed = 3.0 * math.exp(0.8 * rng.gauss(0.0, 1.0))
        cosine, azimuth = rng.uniform(-1.0, 1.0), rng.uniform(0.0, 2.0 * math.pi)
        sine = math.sqrt(1.0 - cosine * cosine)
        offset = rogue_drift.drift((speed * sine * math.cos(azimuth), speed * sine * math.sin(azimuth), speed * cosine), rng.uniform(1000.0, 10000.0))
        xs.append(offset.x_pc)
        zs.append(offset.z_pc)
    assert math.sqrt(sum(x * x for x in xs) / len(xs)) == pytest.approx(rd.RADIAL_SMOOTH_PC, rel=0.15)
    assert math.sqrt(sum(z * z for z in zs) / len(zs)) == pytest.approx(rd.VERTICAL_SMOOTH_PC, rel=0.15)


def test_the_real_galaxy_density_works_as_the_density_function():
    from planetgen.galaxy import density
    from planetgen.physics import mathcheck
    try:
        shape = mathcheck._reference_galaxy_shape()
    except Exception as error:      # the isochrone data are fetched at setup (update.sh); a machine without them cannot build the shape
        pytest.skip(f"no reference galaxy shape: {error}")

    def at(x, y, z):
        return density.total_density_array(x, y, z, shape)

    old = rd.old_density(at, 8000.0, 20.0)
    local = float(at(np.array([8000.0]), np.array([0.0]), np.array([20.0]))[0])
    assert old > 0.0
    assert 0.3 * local < old < 3.0 * local
    far = rd.old_density(at, 20000.0, 20.0)
    assert far < old
