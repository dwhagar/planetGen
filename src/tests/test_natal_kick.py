# tests/test_natal_kick.py

"""GEN.554: the share of a galaxy's neutron stars its natal kicks send out of it follows from the kick law and the galaxy's
escape speed (`planetgen.galaxy.natal_kick`, docs/design/systems/natal-kick.md), and the old 20% to 40% draw is a test target."""

import math

import numpy as np
import pytest

from planetgen import tuning
from planetgen.galaxy import mass_chain, natal_kick, structure as galaxyStructure, structure_draw
from tests.wide_masses import wide_mass

MILKY_WAY_STARS = 5.0e10
MILKY_WAY_HALO = 1.15e12           # the published halo mass the repository's tests quote (test_sector_gen.py, TEST.245)


def _chain(stellar_mass, halo_mass):
    return mass_chain.MassChain(stellar_mass, halo_mass, mass_chain.virial_radius_kpc(halo_mass),
                                mass_chain.median_concentration(halo_mass))


@pytest.mark.parametrize("speed,threshold", [(0.8, 2.4), (0.5, 3.0), (1.5, 2.0), (0.0, 2.5), (0.3, 0.5)])
def test_the_survival_matches_the_kick_drawn_at_random(speed, threshold):
    rng = np.random.default_rng(554)
    kick = rng.standard_normal((400_000, 3))
    kick[:, 0] += speed
    drawn = float((np.linalg.norm(kick, axis=1) > threshold).mean())
    assert natal_kick.survival(speed, threshold) == pytest.approx(drawn, abs=0.004)


def test_the_survival_falls_with_the_threshold_and_rises_with_the_progenitors_speed():
    thresholds = [natal_kick.survival(0.8, b) for b in (0.5, 1.0, 2.0, 3.0, 4.0)]
    assert thresholds == sorted(thresholds, reverse=True) and thresholds[0] <= 1.0
    speeds = [natal_kick.survival(a, 2.5) for a in (0.0, 0.4, 0.8, 1.6)]
    assert speeds == sorted(speeds)
    assert natal_kick.survival(0.5, 0.0) == 1.0


def test_the_birth_radii_divide_the_disc_into_equal_shares():
    radii = natal_kick.birth_radii_over_scale_length()
    assert len(radii) == tuning.NEUTRON_STAR_KICK_BIRTH_RADII
    assert list(radii) == sorted(radii)
    assert sum(radii) / len(radii) == pytest.approx(2.0, rel=0.02)       # the mean radius of R exp(-R / L) is 2 L
    share = lambda x: 1.0 - (1.0 + x) * math.exp(-x)
    assert share(radii[len(radii) // 2]) == pytest.approx(0.5 + 0.5 / len(radii), abs=1e-6)


def test_the_escape_speed_falls_outward_and_a_deeper_halo_holds_more():
    chain = _chain(MILKY_WAY_STARS, MILKY_WAY_HALO)
    speeds = [natal_kick.escape_speed_kms(chain, 2.6, r) for r in (1.0, 2.6, 8.2, 20.0, 100.0)]
    assert speeds == sorted(speeds, reverse=True)
    heavier = _chain(MILKY_WAY_STARS, 2.0 * MILKY_WAY_HALO)
    assert natal_kick.escape_speed_kms(heavier, 2.6, 8.2) > speeds[2]


def test_a_galaxy_of_the_milky_ways_masses_has_the_escape_speed_the_milky_way_is_measured_to_have():
    """About 550 km/s at the solar circle ([R3523], 533 +54/-41), a derived check and no input."""
    speed = natal_kick.escape_speed_kms(_chain(MILKY_WAY_STARS, MILKY_WAY_HALO), 2.6, 8.2)
    assert 490.0 < speed < 600.0


def test_the_milky_ways_neutron_stars_are_ejected_inside_the_old_target_range():
    low, high = tuning.NEUTRON_STAR_UNBOUND_FRACTION_TARGET_RANGE
    share = natal_kick.neutron_star_unbound_share(_chain(MILKY_WAY_STARS, MILKY_WAY_HALO), 220.0, 2.6)
    assert low <= share <= high


def test_the_reference_structure_holds_the_median_halos_share_and_shows_the_gap():
    """The median halo of 5e10 solar masses of stars (1.6e12) ejects fewer than the old 20% to 40% target: the gap is shown."""
    low, high = tuning.NEUTRON_STAR_UNBOUND_FRACTION_TARGET_RANGE
    share = galaxyStructure.REFERENCE_STRUCTURE.neutron_star_unbound_share
    assert 0.5 * low <= share <= high
    halo = mass_chain.median_halo_mass(MILKY_WAY_STARS)
    assert share == pytest.approx(natal_kick.neutron_star_unbound_share(_chain(MILKY_WAY_STARS, halo), 220.0, 2.6))


def test_a_faster_kick_ejects_more_and_a_faster_disc_or_deeper_well_ejects_fewer():
    chain = _chain(MILKY_WAY_STARS, MILKY_WAY_HALO)
    base = natal_kick.neutron_star_unbound_share(chain, 220.0, 2.6)
    assert natal_kick.neutron_star_unbound_share(chain, 220.0, 2.6, sigma_kms=1.3 * tuning.NEUTRON_STAR_KICK_SIGMA_KMS) > base
    assert natal_kick.neutron_star_unbound_share(chain, 220.0, 2.6, sigma_kms=0.7 * tuning.NEUTRON_STAR_KICK_SIGMA_KMS) < base
    assert natal_kick.neutron_star_unbound_share(_chain(MILKY_WAY_STARS, 3.0 * MILKY_WAY_HALO), 220.0, 2.6) < base


def test_a_seed_draws_the_share_its_galaxy_has_and_the_same_every_time():
    for n in range(20):
        seed = n.to_bytes(16, "big")
        share = structure_draw.draw_neutron_star_unbound_share(seed)
        assert share == structure_draw.draw_structure(seed).neutron_star_unbound_share
        assert share == structure_draw.draw_neutron_star_unbound_share(seed)
        assert 0.0 <= share <= 1.0


def test_small_galaxies_lose_more_of_their_neutron_stars_than_large_ones():
    seeds = [n.to_bytes(16, "big") for n in range(120)]
    rows = sorted((wide_mass(seed), structure_draw.draw_neutron_star_unbound_share(seed, wide_mass(seed))) for seed in seeds)
    small, large = [share for _m, share in rows[:30]], [share for _m, share in rows[-30:]]
    assert sum(small) / 30 > sum(large) / 30 + 0.2
    assert min(share for _m, share in rows) < 0.1 and max(share for _m, share in rows) > 0.5


def test_the_ejected_speeds_are_ascending_quantiles_that_rise_with_the_kick_and_fall_with_the_well():
    chain = _chain(MILKY_WAY_STARS, MILKY_WAY_HALO)
    speeds = natal_kick.ejecta_speeds_kms(chain, 220.0, 2.6)
    assert len(speeds) == tuning.EJECTA_SPEED_QUANTILES and list(speeds) == sorted(speeds) and speeds[0] > 0.0
    assert sum(speeds) / len(speeds) < 700.0
    faster = natal_kick.ejecta_speeds_kms(chain, 220.0, 2.6, sigma_kms=350.0)
    assert sum(faster) / len(faster) > sum(speeds) / len(speeds)
    deeper = natal_kick.ejecta_speeds_kms(_chain(MILKY_WAY_STARS, 4.0 * MILKY_WAY_HALO), 220.0, 2.6)
    assert sum(deeper) / len(deeper) < sum(speeds) / len(speeds)


def test_the_ejected_speed_law_integrates_to_the_unbound_share():
    """The weights over the speeds above the escape speed add up to `survival` at that radius."""
    chain = _chain(MILKY_WAY_STARS, MILKY_WAY_HALO)
    sigma = tuning.NEUTRON_STAR_KICK_SIGMA_KMS
    a = 220.0 / sigma
    escape = natal_kick.escape_speed_kms(chain, 2.6, 8.2) / sigma
    steps = int(natal_kick._EJECTA_SPEED_SPAN / natal_kick._EJECTA_SPEED_STEP)
    total = sum(natal_kick._speed_density(escape + (i + 0.5) * natal_kick._EJECTA_SPEED_STEP, a)
                for i in range(steps)) * natal_kick._EJECTA_SPEED_STEP
    assert total == pytest.approx(natal_kick.survival(a, escape), rel=2e-3)


def test_a_galaxy_with_a_well_too_deep_to_leave_ejects_nothing_and_has_no_speeds():
    deep = _chain(MILKY_WAY_STARS, 1.0e15)
    assert natal_kick.ejecta_speeds_kms(deep, 220.0, 2.6) == ()
    assert galaxyStructure.REFERENCE_STRUCTURE.neutron_star_ejecta_speeds_kms == natal_kick.ejecta_speeds_kms(
        galaxyStructure._MILKY_WAY_CHAIN, 220.0, 2.6)
