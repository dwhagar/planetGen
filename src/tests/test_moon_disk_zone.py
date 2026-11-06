"""GEN.1244: the regular moons of a giant form inside the centrifugal radius of its disk, outside the corotation radius of its spin."""

import math

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.physics import constants, moon_disk_zone
from tests import test_fuzz_system_generation as fuzz

JUPITER_MASS_KG = 1.898e27
JUPITER_RADIUS_KM = 69911.0
JUPITER_SPIN = 2.0 * math.pi / (9.925 * 3600.0)


def test_the_corotation_radius_of_jupiters_spin_is_two_and_a_quarter_radii():
    radius = moon_disk_zone.corotation_km(JUPITER_MASS_KG, JUPITER_SPIN)
    assert radius / JUPITER_RADIUS_KM == pytest.approx(2.25, rel=0.03)


def test_a_planet_that_does_not_spin_has_its_corotation_at_infinity():
    assert moon_disk_zone.corotation_km(JUPITER_MASS_KG, 0.0) == math.inf
    assert moon_disk_zone.corotation_km(JUPITER_MASS_KG, None) == math.inf


def test_the_zone_runs_from_corotation_to_thirty_radii():
    inner, outer = moon_disk_zone.zone_km(JUPITER_MASS_KG, JUPITER_RADIUS_KM, JUPITER_SPIN, 1.5 * JUPITER_RADIUS_KM, 400.0 * JUPITER_RADIUS_KM)
    assert outer == moon_disk_zone.CENTRIFUGAL_RADII * JUPITER_RADIUS_KM
    assert inner == pytest.approx(moon_disk_zone.corotation_km(JUPITER_MASS_KG, JUPITER_SPIN))


def test_a_planet_that_spins_too_slowly_has_no_cavity():
    low = 1.5 * JUPITER_RADIUS_KM
    inner, outer = moon_disk_zone.zone_km(JUPITER_MASS_KG, JUPITER_RADIUS_KM, JUPITER_SPIN / 1000.0, low, 400.0 * JUPITER_RADIUS_KM)
    assert inner == low and inner < outer


def test_the_zone_never_passes_the_stable_limit_of_the_orbit():
    _inner, outer = moon_disk_zone.zone_km(JUPITER_MASS_KG, JUPITER_RADIUS_KM, JUPITER_SPIN, 1.5 * JUPITER_RADIUS_KM, 10.0 * JUPITER_RADIUS_KM)
    assert outer == 10.0 * JUPITER_RADIUS_KM


def test_the_disk_moons_of_generated_giants_sit_inside_the_zone(monkeypatch):
    """Twelve systems of a sun-like star: every moon born in a giant's disk was born inside the zone it was drawn for, whose outer edge is 30 radii of the giant
    as it was when its moons were made (it was 46 to 130 before GEN.1244; the tide may carry a moon out later). TEST.657: the giant then contracts to its radius at
    the star's final age (`_settle_to_final_age`), a few percent smaller, so the test reads the zone each giant's moons were drawn in, not its final radius."""
    zones = {}
    real = moon_disk_zone.zone_km

    def recording(mass_kg, radius_km, spin_rate_rad_s, low_km, high_km):
        zone = real(mass_kg, radius_km, spin_rate_rad_s, low_km, high_km)
        zones[mass_kg] = zone
        return zone

    monkeypatch.setattr(moon_disk_zone, "zone_km", recording)
    found = 0
    for seed in range(12):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        cfg.BINARY_SYSTEM = False
        cfg.MOONS = True
        system = fuzz._generate(cfg, seed)
        for planet in system.planets:
            if planet.body_type != "g":
                continue
            for moon in planet.moons:
                if moon.moon_origin == "disk":
                    found += 1
                    outer_km = zones[planet.mass][1]
                    assert getattr(moon, "_birth_distance_au", moon.distance) * constants.AU_TO_KM <= outer_km * 1.001
                    assert outer_km <= moon_disk_zone.CENTRIFUGAL_RADII * planet.radius * 1.06     # and the zone is 30 radii of the giant, whose radius has since shrunk by under 6 percent
    assert found > 0
