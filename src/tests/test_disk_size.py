"""The gas disk's outer edge (GEN.1520, physics/disk_size.py): 100 au for the Sun, the star's mass to the 0.5, a factor of 2 of scatter drawn once per star."""

import math

import pytest

from planetgen.generation.star import Star
from planetgen.generation.config import SystemConfig
from planetgen.physics import disk_size
from planetgen.util import draw


def test_the_median_edge_is_one_hundred_au_for_the_sun_and_follows_the_mass_to_the_half():
    assert disk_size.median_edge_au(1.0) == pytest.approx(100.0)
    assert disk_size.median_edge_au(0.25) == pytest.approx(50.0)
    assert disk_size.median_edge_au(0.2) == pytest.approx(44.7, rel=0.01)
    assert disk_size.median_edge_au(1.4) == pytest.approx(118.3, rel=0.01)
    with pytest.raises(ValueError):
        disk_size.median_edge_au(0.0)


def test_the_scatter_is_a_factor_of_two_a_sigma_and_cut_at_two_sigma():
    assert disk_size.edge_au_from(1.0, 1.0) == pytest.approx(200.0, rel=0.01)
    assert disk_size.edge_au_from(1.0, -1.0) == pytest.approx(50.0, rel=0.01)
    draws = [disk_size.draw_edge_au(1.0, draw.Stream(f"t{index}")) for index in range(400)]
    assert min(draws) >= 100.0 * 10.0 ** (-0.3 * disk_size.EDGE_SCATTER_MAX_SIGMA) - 1e-9
    assert max(draws) <= 100.0 * 10.0 ** (0.3 * disk_size.EDGE_SCATTER_MAX_SIGMA) + 1e-9
    logs = sorted(math.log10(value / 100.0) for value in draws)
    assert abs(logs[len(logs) // 2]) < 0.06                    # the median of the draws is the median edge


def test_a_star_gets_the_same_edge_every_time_and_it_moves_no_other_draw():
    config = SystemConfig()
    with draw.bound(5):
        star = Star(system_config=config)
    with draw.bound(6):
        first = draw.random()
    with draw.bound(6):
        edge = disk_size.edge_au(star)
        assert draw.random() == first
    with draw.bound(99):
        assert disk_size.edge_au(star) == edge
    assert disk_size.median_edge_au(disk_size.formation_mass_solar(star)) / 4.0 <= edge <= disk_size.median_edge_au(disk_size.formation_mass_solar(star)) * 4.0


def test_the_growth_limit_is_the_suns_tuned_value_scaled_with_the_disk_and_the_core_taper_starts_there():
    from planetgen import tuning
    from planetgen.physics import planet_mass
    assert disk_size.growth_limit_au(1.0) == pytest.approx(tuning.PLANET_OUTER_GROWTH_AU)
    assert disk_size.growth_limit_au(0.25) == pytest.approx(0.5 * tuning.PLANET_OUTER_GROWTH_AU)
    assert disk_size.growth_limit_au(0.38) < disk_size.growth_limit_au(0.7) < tuning.PLANET_OUTER_GROWTH_AU < disk_size.growth_limit_au(1.3)
    # a 0.25 Msun star's cores fall off from half the Sun's distance: at the same distance beyond it they are lighter relative to the untapered mass than the Sun's
    star_kg, far = 0.25 * 1.98847e30, 3.0 * disk_size.growth_limit_au(0.25)
    tapered = planet_mass.oligarch_core_kg(far, star_kg, 0.8, 1.0, 1.0, 10.0)
    inside = planet_mass.oligarch_core_kg(disk_size.growth_limit_au(0.25), star_kg, 0.8, 1.0, 1.0, 10.0)
    assert tapered < inside
