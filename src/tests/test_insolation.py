# tests/test_insolation.py

"""GEN.148: the summed insolation of a planet in a multi-star system and its `flux_variation`."""

import math

import pytest

from planetgen.physics import insolation


def test_the_mean_inverse_square_of_a_circular_orbit_is_one_over_a_squared():
    assert insolation.mean_inverse_square(50.0, 0.0) == pytest.approx(1 / 2500.0)
    assert insolation.mean_inverse_square(50.0, 0.6) == pytest.approx(1 / (2500.0 * 0.8))


def test_a_companion_at_50_au_adds_about_4e_minus_4_of_its_light_to_a_planet_at_1_au():
    """The design note's figure: a companion of the same luminosity at 50 AU (circular) is 4e-4 of the host's light."""
    swing = insolation.companion_swing(1.0, 1.0, 1.0, 50.0, 50.0)
    assert swing == 0.0  # circular: the companion's light never changes
    eccentric = insolation.companion_swing(1.0, 1.0, 1.0, 25.0, 75.0)
    assert eccentric == pytest.approx((1 / 625.0 - 1 / 5625.0) / (1.0 + 1 / (2500.0 * math.sqrt(1 - 0.25))), rel=1e-9)
    assert 1e-3 < eccentric < 2e-3


def test_a_dim_host_beside_a_bright_eccentric_companion_swings_by_no_more_than_the_whole_light():
    """GEN.597: the peak to peak over the mean passes 1 there; the fraction is capped at 1 as the circumbinary one is."""
    assert insolation.companion_swing(1.0e-4, 5.8, 1.0, 420.0, 1580.0) == 1.0
    assert insolation.companion_swing(1.0e-4, 0.1, 1.0, 1.2, 6.4) == 1.0
    assert 0.0 < insolation.companion_swing(1.0, 1.0, 1.0, 25.0, 75.0) < 1.0


def test_a_closer_and_brighter_companion_swings_the_light_more():
    far = insolation.companion_swing(1.0, 1.0, 1.0, 100.0, 300.0)
    near = insolation.companion_swing(1.0, 1.0, 1.0, 30.0, 90.0)
    brighter = insolation.companion_swing(1.0, 1.0, 50.0, 30.0, 90.0)
    assert 0.0 < far < near < brighter


def test_a_circumbinary_planets_light_swings_by_the_luminosity_centroid_offset():
    """Equal stars of equal mass have the centroid at the barycentre: no swing. A dim light companion makes one."""
    assert insolation.circumbinary_swing(1.0, 1.0, 1.0, 1.0, 0.2, 1.0) == 0.0
    # 1 and 0.1 Msun, 1 and 0.01 Lsun, 0.2 AU apart, planet at 1 AU: r1 = 0.0182, r2 = 0.182,
    # delta = |1 * 0.0182 - 0.01 * 0.182| / 1.01 = 0.0162; 4 delta / d = 0.065.
    swing = insolation.circumbinary_swing(1.0, 0.01, 1.0, 0.1, 0.2, 1.0)
    assert swing == pytest.approx(0.065, abs=0.001)
    assert insolation.circumbinary_swing(1.0, 0.01, 1.0, 0.1, 0.2, 4.0) == pytest.approx(swing / 4)


def test_a_lone_star_swings_only_by_the_orbit_and_a_generated_system_sets_it_on_every_body():
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from planetgen.util import draw
    from planetgen.physics.orbit_elements import flux_swing
    from tests.fuzz_support import deterministic_entropy

    seen = {"close": 0, "wide": 0, None: 0}
    swung = 0
    for seed in range(240):   # GEN.597: seeds 36, 106, 225 and 245 once gave a swing above 1
        with deterministic_entropy(seed):
            draw.set_run_seed(seed)
            system = StarSystem(SystemConfig(), galactic_center_dist_ly=26000.0)
        seen[system.binary_type] += 1
        for planet in system.planets + system.secondary_planets:
            if not hasattr(planet, "moons"):
                continue
            assert planet.flux_variation is not None and 0.0 <= planet.flux_variation <= 1.0    # a giant's too (GEN.784)
            if planet.body_type != "g":
                swung += planet.flux_variation > 0.0
            for moon in planet.moons:
                assert moon.flux_variation is not None and 0.0 <= moon.flux_variation <= 1.0
                if system.binary_type is None:
                    assert moon.flux_variation == pytest.approx(flux_swing(planet.eccentricity))
            if system.binary_type is None:
                assert planet.flux_variation == pytest.approx(flux_swing(planet.eccentricity))
    assert seen[None] > 0 and seen["close"] + seen["wide"] > 0
    assert swung > 0
