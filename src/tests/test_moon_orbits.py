# tests/test_moon_orbits.py
"""MAP.308: a moon's orbit about its planet from the stored state, and the moon system scene drawn from it."""

import math
import re

import pytest

from planetgen.physics import constants
from planetgen.physics.state_vectors import state_from_elements
from planetgen.web.maps import moonorbits, systemmap

PLANET_MASS_KG = 5.97e24
PLANET_RADIUS_KM = 6371.0
MOON_MASS_KG = 7.3e22
MU = constants.G * 1.0e-9 * (PLANET_MASS_KG + MOON_MASS_KG)


def _planet(moons):
    return {"id": 7, "name": "Aster", "planet_class": "M", "body_type": "r", "radius_km": PLANET_RADIUS_KM,
            "mass_kg": PLANET_MASS_KG, "moons": moons, "atmosphere": "None", "distance_km": 1.5e8,
            "zone": "e", "period_years": 1.0, "gravity_g": 1.0, "life_class": None}


def _moon(index, periapsis_km, eccentricity, inclination_deg=10.0, radius_km=1737.0, true_anomaly=1.0):
    position, velocity = state_from_elements(MU, periapsis_km, eccentricity, math.radians(inclination_deg),
                                             0.7, 0.4, true_anomaly)
    return {"id": 100 + index, "name": f"Moon {index}", "planet_class": "C", "body_type": "r", "radius_km": radius_km,
            "mass_kg": MOON_MASS_KG, "zone": "e", "period_years": 0.07, "gravity_g": 0.16,
            "life_class": None, "distance_km": periapsis_km / (1 - eccentricity),
            "position_x_km": position[0], "position_y_km": position[1], "position_z_km": position[2],
            "velocity_x_kms": velocity[0], "velocity_y_kms": velocity[1], "velocity_z_kms": velocity[2],
            "atmosphere": "None"}


def test_a_moon_orbit_is_the_ellipse_its_state_gives():
    moon = _moon(1, 300_000.0, 0.3)
    orbit = moonorbits.moon_orbit(_planet([moon]), moon)
    assert not orbit["circular"]
    assert orbit["eccentricity"] == pytest.approx(0.3, abs=2e-3)
    assert orbit["periapsis_km"] == pytest.approx(300_000.0, rel=2e-3)
    assert orbit["apoapsis_km"] == pytest.approx(300_000.0 * 1.3 / 0.7, rel=2e-3)
    assert len(orbit["points"]) == moonorbits.SAMPLES
    # The moon is on its own orbit: the sampled path passes close to its projected position.
    nearest = min(math.hypot(x - moon["position_x_km"], y - moon["position_y_km"]) for x, y in orbit["points"])
    assert nearest < 0.1 * orbit["apoapsis_km"]
    # No point is farther from the planet than the apoapsis.
    assert max(math.hypot(x, y) for x, y in orbit["points"]) <= orbit["apoapsis_km"] * 1.0001


def test_a_moon_without_a_state_is_the_circle_of_its_distance():
    moon = {"id": 5, "name": "Lone", "distance_km": 250_000.0}
    orbit = moonorbits.moon_orbit(_planet([moon]), moon)
    assert orbit["circular"] and orbit["eccentricity"] == 0.0
    assert orbit["apoapsis_km"] == orbit["periapsis_km"] == 250_000.0
    assert all(math.hypot(x, y) == pytest.approx(250_000.0) for x, y in orbit["points"])
    assert moonorbits.moon_orbit(_planet([]), {"id": 6, "name": "None"})["points"] == []


def test_a_hyperbolic_state_falls_back_to_a_circle():
    moon = _moon(2, 300_000.0, 0.2)
    moon["velocity_x_kms"] *= 10.0
    moon["velocity_y_kms"] *= 10.0
    moon["velocity_z_kms"] *= 10.0
    orbit = moonorbits.moon_orbit(_planet([moon]), moon)
    assert orbit["circular"] and orbit["points"]


def _scenes(moons):
    markup = systemmap._render_moon_scene(_planet(moons))
    return {match.group(1): match.group(0) for match in re.finditer(
        r'<svg class="sysmap-svg(?! sysmap-orbits-layer)[^"]*" data-scene="([^"]+)".*?</svg>', markup, re.S)}


def test_a_planet_with_moons_gets_two_scenes_on_one_linear_scale():
    moons = [_moon(1, 380_000.0, 0.05), _moon(2, 900_000.0, 0.2, true_anomaly=2.5)]
    scenes = _scenes(moons)
    assert set(scenes) == {"planet-7", "planet-7-true"}
    for scene_id, svg in scenes.items():
        assert "data-pxperkm=" in svg and "data-lokm" not in svg, scene_id
    assert 'data-scale-twin="planet-7-true"' in scenes["planet-7"]
    assert 'data-scale-twin="planet-7"' in scenes["planet-7-true"]
    assert 'data-scale-true="true"' in scenes["planet-7-true"] and 'data-scale-true="enlarged"' in scenes["planet-7"]


def _radii(svg, klass):
    return [float(m) for m in re.findall(rf'<g class="sysmap-body {klass}"[^>]*>.*?<circle class="sysmap-body-fill"[^>]*? r="([\d.]+)"', svg, re.S)]


def test_the_true_scale_scene_draws_the_planet_to_the_scale_of_its_orbits():
    moons = [_moon(1, 380_000.0, 0.05), _moon(2, 900_000.0, 0.2, true_anomaly=2.5)]
    svg = _scenes(moons)["planet-7-true"]
    per_km = float(re.search(r'data-pxperkm="([\d.e+-]+)"', svg).group(1))
    orbits = [moonorbits.moon_orbit(_planet(moons), m) for m in moons]
    # The outermost apoapsis is the scene's reach from the middle.
    assert max(o["apoapsis_km"] for o in orbits) * per_km == pytest.approx(systemmap._MOON_SCENE_REACH_PX, rel=1e-4)
    planet_px = _radii(svg, "sysmap-planet")[0]
    assert planet_px == pytest.approx(max(PLANET_RADIUS_KM * per_km, systemmap._TRUE_PLANET_MIN_PX), abs=0.06)
    # The planet's radius to an orbit's radius is the real ratio: within a pixel.
    orbit_px = orbits[1]["periapsis_km"] * per_km
    assert abs(planet_px / orbit_px - PLANET_RADIUS_KM / orbits[1]["periapsis_km"]) * orbit_px < 1.0 + systemmap._TRUE_PLANET_MIN_PX


def test_the_default_scene_enlarges_every_body_by_one_factor_and_says_so():
    moons = [_moon(1, 380_000.0, 0.05, radius_km=1737.0), _moon(2, 900_000.0, 0.2, radius_km=2400.0, true_anomaly=2.5)]
    scenes = _scenes(moons)
    svg, true_svg = scenes["planet-7"], scenes["planet-7-true"]
    per_km = float(re.search(r'data-pxperkm="([\d.e+-]+)"', svg).group(1))
    assert per_km == pytest.approx(float(re.search(r'data-pxperkm="([\d.e+-]+)"', true_svg).group(1)), rel=1e-5)
    enlarge = float(re.search(r'data-enlarge="([\d.]+)"', svg).group(1))
    assert enlarge > 3, "an Earth at this scale is a few pixels: the default enlarges it"
    planet_px = _radii(svg, "sysmap-planet")[0]
    assert planet_px == pytest.approx(PLANET_RADIUS_KM * per_km * enlarge, abs=0.06)
    assert 14.0 <= planet_px <= 40.0
    # It stays well clear of the nearest moon's closest approach.
    assert planet_px <= systemmap._ENLARGED_CLEARANCE * 380_000.0 * per_km * 1.0001 + 0.06
    moon_px = _radii(svg, "sysmap-moon")
    assert moon_px[0] == pytest.approx(max(1737.0 * per_km * enlarge, systemmap._ENLARGED_MOON_MIN_PX), abs=0.06)
    assert moon_px[1] > moon_px[0], "the larger moon is drawn larger"
    assert "times their true size" in svg and "true size and distance" in true_svg
    # A planet that is already large against its moons' orbits is not enlarged past the clearance.
    close_scenes = _scenes([_moon(1, 40_000.0, 0.01)])
    assert set(close_scenes) == {"planet-7"}, "nothing is enlarged, so there is no second scene to switch to"
    assert "data-scale-twin" not in close_scenes["planet-7"]
    close = close_scenes["planet-7"]
    close_per_km = float(re.search(r'data-pxperkm="([\d.e+-]+)"', close).group(1))
    assert _radii(close, "sysmap-planet")[0] <= max(systemmap._ENLARGED_CLEARANCE * 40_000.0 * close_per_km, PLANET_RADIUS_KM * close_per_km) + 0.06
    # Each moon's orbit is an ellipse path, not a circle.
    assert svg.count('<path class="sysmap-orbit"') == 0     # the orbits live in the orbits layer


def test_every_moon_orbit_is_a_path_in_the_orbits_layer():
    moons = [_moon(1, 380_000.0, 0.3), _moon(2, 900_000.0, 0.1, true_anomaly=2.5)]
    markup = systemmap._render_moon_scene(_planet(moons))
    layers = re.findall(r'<svg class="sysmap-svg sysmap-orbits-layer[^"]*" data-scene="([^"]+)".*?</svg>', markup, re.S)
    assert layers == ["planet-7", "planet-7-true"]
    assert markup.count('<path class="sysmap-orbit"') == 4 and 'class="sysmap-orbit" cx=' not in markup
