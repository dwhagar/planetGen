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


def test_a_regular_planet_a_lone_planet_and_one_with_no_axis_keep_their_tilt(monkeypatch, quiet):
    keep = _planet(1.0)
    assert planets.settle_obliquity([keep], SUN) == 0                                   # a lone planet has no modes
    monkeypatch.setattr(ov, "analyse", lambda alpha, tilt, modes: ov.Variation(20.0, 26.0, False, False, 1.0, (), alpha))
    regular = _planet(1.0)
    assert planets.settle_obliquity([regular, _planet(5.0, body_type="g")], SUN) == 0 and regular.axial_tilt_deg == 23.0
    _chaotic(monkeypatch)
    bare = _planet(1.0)
    bare.spin_axis_x = None
    planets.settle_obliquity([bare, _planet(5.0, body_type="g")], SUN)
    assert bare.axial_tilt_deg == 23.0


def _locked(distance=0.4, **kwargs):
    planet = _planet(distance, rotation_hours=distance ** 1.5 * 365.25 * 24.0, inclination=7.0, **kwargs)
    assert planets.spin.spin_state(planet.rotation_period_hours, planet.period)[0] == "locked"
    return planet


def test_a_locked_rocky_planet_takes_the_cassini_tilt_of_its_figure_and_a_giant_or_a_free_one_does_not(monkeypatch, quiet):
    """GEN.1536: the tilt is the first Cassini state against the nodal modes of the orbit, not the early share of the inclination; the sky and moons follow."""
    monkeypatch.setattr(planets, "internal_heat_flux_w_m2", lambda body: 0.02)
    locked = _locked(moons=[SimpleNamespace(mass=1.0e20, distance=0.001, eccentricity=0.0)])
    giant = _locked(1.5, body_type="g")
    outer = _planet(5.0, body_type="g")
    assert planets.settle_obliquity([locked, giant, outer], SUN) >= 1
    assert 0.0 <= locked.axial_tilt_deg < 7.0 and locked.axial_tilt_deg != 23.0
    assert giant.axial_tilt_deg == 23.0 and outer.axial_tilt_deg == 23.0
    assert ("daylight", locked.name) in quiet and ("moons", locked.name) in quiet
    axis = (locked.spin_axis_x, locked.spin_axis_y, locked.spin_axis_z)
    normal = planets.spin.orbit_normal(7.0, 10.0)
    assert math.degrees(math.acos(sum(a * b for a, b in zip(axis, normal)))) == pytest.approx(locked.axial_tilt_deg, abs=1e-6)


def test_the_cassini_tilt_is_the_same_for_the_same_planet_and_a_hot_one_keeps_the_hydrostatic_figure(monkeypatch, quiet):
    monkeypatch.setattr(planets, "internal_heat_flux_w_m2", lambda body: 0.02)
    first, second = _locked(), _locked()
    planets.settle_obliquity([first, _planet(5.0, body_type="g")], SUN)
    planets.settle_obliquity([second, _planet(5.0, body_type="g")], SUN)
    assert first.axial_tilt_deg == second.axial_tilt_deg
    monkeypatch.setattr(planets, "internal_heat_flux_w_m2", lambda body: 5.0)
    hot = _locked()
    planets.settle_obliquity([hot, _planet(5.0, body_type="g")], SUN)
    assert hot.axial_tilt_deg > first.axial_tilt_deg      # the relaxed figure is the fossil one divided by its factor (22.7 to 94): the axis is held less to the orbit


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


def test_the_range_of_the_tilt_is_kept_for_a_planet_that_is_worked_and_not_for_a_locked_one(monkeypatch, quiet):
    """GEN.1535: `obliquity_min_deg` and `obliquity_max_deg` are the analysed range, set whether or not the zone is chaotic; a locked planet takes the Cassini tilt (GEN.1536) and has no range."""
    monkeypatch.setattr(ov, "analyse", lambda alpha, tilt, modes: ov.Variation(18.0, 31.0, False, False, 0.0, (), alpha))
    free, locked = _planet(1.0), _planet(5.0, rotation_hours=5.0 ** 1.5 * 8766.0)
    free.obliquity_min_deg = free.obliquity_max_deg = locked.obliquity_min_deg = locked.obliquity_max_deg = None
    assert planets.settle_obliquity([free, locked], SUN) == 1       # the locked one
    assert (free.obliquity_min_deg, free.obliquity_max_deg) == (18.0, 31.0)
    assert locked.obliquity_min_deg is None and locked.obliquity_max_deg is None


def test_a_lone_planet_with_no_companion_is_given_the_range_of_its_own_tilt():
    lone = _planet(1.0, tilt=23.0)
    lone.obliquity_min_deg = lone.obliquity_max_deg = None
    assert planets.settle_obliquity([lone], SUN) == 0
    assert lone.obliquity_min_deg == lone.obliquity_max_deg == 23.0


def _companion(normal, axis_au=50.0, mass=SUN * 0.5, eccentricity=0.3):
    from planetgen.physics import radiation_field
    return radiation_field.Companion(SimpleNamespace(mass=mass), 1, axis_au, eccentricity, (axis_au, 0.0, 0.0), normal)


def test_a_companion_in_a_plane_of_its_own_drives_a_mode_in_the_orbit_of_a_planet():
    """GEN.1535: the mode's amplitude is the mutual inclination and its rate is negative (a prograde companion regresses the node)."""
    from planetgen.physics import radiation_field
    planet = _planet(1.0, inclination=0.0, node=0.0)
    tilt = math.radians(30.0)
    radiation_field.attach(planet.star, [_companion((0.0, -math.sin(tilt), math.cos(tilt)))])
    modes = planets.companion_modes(planet, SUN)
    assert len(modes) == 1 and modes[0].frequency_rad_yr < 0.0 and math.degrees(modes[0].amplitude_rad) == pytest.approx(30.0, abs=0.01)


def test_a_coplanar_or_unplaced_or_nearby_companion_drives_no_mode():
    from planetgen.physics import radiation_field
    planet = _planet(1.0, inclination=0.0, node=0.0)
    radiation_field.attach(planet.star, [_companion((0.0, 0.0, 1.0)), _companion(None), _companion((0.0, -0.5, 0.866), axis_au=3.0)])
    assert planets.companion_modes(planet, SUN) == ()
    assert planets.companion_modes(_planet(1.0), SUN) == ()      # no companions attached


def test_a_planet_in_a_chaotic_zone_is_named_by_the_wandering_tilt_oddity_and_a_regular_one_is_not(monkeypatch, quiet):
    """GEN.1535: oddity W6 (W61 for a giant) names the planet `settle_obliquity` found chaotic; the rows draw nothing (rate 1), so the table names exactly those."""
    _chaotic(monkeypatch)
    inner, outer = _planet(1.0), _planet(5.0, body_type="g")
    inner.oddities = outer.oddities = ()
    planets.settle_obliquity([inner, outer], SUN)
    assert "W6" in inner.oddities and "W61" in outer.oddities
    monkeypatch.setattr(ov, "analyse", lambda alpha, tilt, modes: ov.Variation(20.0, 26.0, False, False, 3.0, (), alpha))
    regular = _planet(1.0)
    regular.oddities = ()
    planets.settle_obliquity([regular, _planet(5.0)], SUN)
    assert regular.oddities == ()
