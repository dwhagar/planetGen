"""
GEN.1454 part 2: `planets.settle_obliquity` gives a free-spinning planet in a chaotic zone of the spin-orbit resonances a tilt drawn from the zone, and leaves every other planet alone.
"""

import math
from types import SimpleNamespace

import pytest

from planetgen.physics import obliquity_variation as ov, planets

SUN = 1.98847e30
EARTH = 5.972e24


def _planet(distance, tilt=23.0, rotation_hours=24.0, moons=(), inclination=1.0, node=10.0, body_type="t"):
    return SimpleNamespace(body_type=body_type, is_moon=False, period=distance ** 1.5, mass=EARTH if body_type == "t" else 300 * EARTH, distance=distance, radius=6371.0,
                           orbital_inclination_deg=inclination, orbital_ascending_node_deg=node, axial_tilt_deg=tilt, rotation_period_hours=rotation_hours,
                           spin_axis_x=math.sin(math.radians(tilt)), spin_axis_y=0.0, spin_axis_z=math.cos(math.radians(tilt)), eccentricity=0.02, moons=list(moons),
                           star=SimpleNamespace(), name=f"P{distance}")


@pytest.fixture
def quiet(monkeypatch):
    """The things a changed tilt re-makes are not under test here."""
    calls = []
    monkeypatch.setattr(planets.daylight, "generate", lambda body: calls.append(("daylight", body.name)))
    monkeypatch.setattr(planets.sky_suns, "suns_of", lambda *args, **kwargs: ())
    monkeypatch.setattr(planets.sky_suns, "encode", lambda suns: "")
    monkeypatch.setattr(planets, "orient_moons", lambda body: calls.append(("moons", body.name)))
    return calls


def _chaotic(monkeypatch, minimum=10.0, maximum=60.0):
    monkeypatch.setattr(ov, "analyse", lambda alpha, tilt, modes: ov.Variation(minimum, maximum, True, True, 0.0, (), alpha))


def test_a_planet_in_a_chaotic_zone_takes_a_tilt_from_the_zone_and_the_sky_and_moons_follow(monkeypatch, quiet):
    _chaotic(monkeypatch)
    moon = SimpleNamespace(mass=7.3e22, distance=0.00257, eccentricity=0.05)
    inner, outer = _planet(1.0, moons=[moon]), _planet(5.0, body_type="g")
    assert planets.settle_obliquity([inner, outer], SUN) >= 1
    assert 10.0 <= inner.axial_tilt_deg <= 60.0 and inner.axial_tilt_deg != 23.0
    assert ("daylight", inner.name) in quiet and ("moons", inner.name) in quiet
    axis = (inner.spin_axis_x, inner.spin_axis_y, inner.spin_axis_z)
    assert math.sqrt(sum(c * c for c in axis)) == pytest.approx(1.0)
    normal = planets.spin.orbit_normal(1.0, 10.0)
    assert math.degrees(math.acos(sum(a * b for a, b in zip(axis, normal)))) == pytest.approx(inner.axial_tilt_deg, abs=1e-6)


def test_the_tilt_is_the_same_for_the_same_planet_and_comes_from_the_planets_own_stream(monkeypatch, quiet):
    _chaotic(monkeypatch)
    first, second = _planet(1.0), _planet(1.0)
    planets.settle_obliquity([first, _planet(5.0, body_type="g")], SUN)
    planets.settle_obliquity([second, _planet(5.0, body_type="g")], SUN)
    assert first.axial_tilt_deg == second.axial_tilt_deg


def test_a_regular_planet_a_lone_planet_a_locked_planet_and_one_with_no_axis_keep_their_tilt(monkeypatch, quiet):
    keep = _planet(1.0)
    assert planets.settle_obliquity([keep], SUN) == 0                                   # a lone planet has no modes
    monkeypatch.setattr(ov, "analyse", lambda alpha, tilt, modes: ov.Variation(20.0, 26.0, False, False, 1.0, (), alpha))
    regular = _planet(1.0)
    assert planets.settle_obliquity([regular, _planet(5.0, body_type="g")], SUN) == 0 and regular.axial_tilt_deg == 23.0
    _chaotic(monkeypatch)
    locked = _planet(1.0, rotation_hours=1.0 * 1.0 ** 1.5 * 365.25 * 24.0)               # one turn an orbit
    assert planets.spin.spin_state(locked.rotation_period_hours, locked.period)[0] == "locked"
    planets.settle_obliquity([locked, _planet(5.0, body_type="g")], SUN)
    assert locked.axial_tilt_deg == 23.0
    bare = _planet(1.0)
    bare.spin_axis_x = None
    planets.settle_obliquity([bare, _planet(5.0, body_type="g")], SUN)
    assert bare.axial_tilt_deg == 23.0
    assert quiet == [] or all(entry[1] != locked.name for entry in quiet)


def test_the_inputs_of_the_analysis_are_the_planets_own_moons_orbits_and_spin(monkeypatch, quiet):
    seen = []

    def spy(alpha, tilt, modes):
        seen.append((alpha, tilt, modes))
        return ov.Variation(tilt, tilt, False, False, 0.0, (), alpha)

    monkeypatch.setattr(ov, "analyse", spy)
    moon = SimpleNamespace(mass=7.3e22, distance=0.00257, eccentricity=0.05)
    with_moon, bare = _planet(1.0, moons=[moon]), _planet(1.0)
    giant = _planet(5.0, body_type="g", inclination=3.0)
    planets.settle_obliquity([with_moon, giant], SUN)
    steady, tilt, modes = seen[0]
    seen.clear()
    planets.settle_obliquity([bare, giant], SUN)
    assert steady > 2.5 * seen[0][0] and tilt == 23.0                                      # the Moon-sized moon raises the precession constant
    assert len(modes) >= 1 and all(isinstance(m, ov.Mode) for m in modes)
