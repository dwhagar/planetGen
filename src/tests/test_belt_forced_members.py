# tests/test_belt_forced_members.py

"""
GEN.1492 slice 1: a belt's members take the secular forcing of a companion star (`physics/belt_forcing.py`, wired into `AsteroidBelt._draw_members`). A belt of a lone star must give the
members of before the change to the last digit; a companion raises the mean eccentricity of a belt inside the theory; a giant planet beside the belt dilutes it.
"""

import math

import pytest

from planetgen.generation.belt import AsteroidBelt
from planetgen.generation.system import StarSystem
from planetgen.physics import belt_roster, radiation_field as rf, spin, vec3
from planetgen.util import draw
from tests.test_radiation_field import star as plain_star
from tests.test_systems import make_config

JUPITER = 1.898e27


@pytest.fixture()
def system():
    return StarSystem(system_config=make_config("G2V", PLANETS=False))


@pytest.fixture()
def host(system):
    star = system.star
    rf.attach(star, ())
    yield star
    rf.attach(star, ())


def _members(system, host):
    with draw.bound(5):
        belt = AsteroidBelt(system.system_config, host, 2.5, 2.2, 3.4)
    return belt, belt.members


def _companion(host, axis_au, eccentricity, mass_ratio=0.5):
    companion = plain_star(host.luminosity)
    companion.mass = mass_ratio * host.mass
    rf.attach(host, [rf.Companion(companion, 1, axis_au, eccentricity)])


def _unforced(system, host, monkeypatch):
    """The members the belt draws under the companion's light with the forcing switched off."""
    with monkeypatch.context() as patch:
        patch.setattr(AsteroidBelt, "_force_by_companions", lambda self, star, rows: rows)
        _belt, rows = _members(system, host)
        return [dict(r) for r in rows]


def _mean_e(rows):
    return sum(r["eccentricity"] for r in rows) / len(rows)


def test_the_members_of_a_lone_star_are_what_they_were(system, host):
    belt, rows = _members(system, host)
    assert rows and belt.members is rows
    belt.give_secular_planets([5.2], [JUPITER])
    assert belt.members == rows


def test_a_companion_raises_the_mean_eccentricity_and_keeps_every_other_field(system, host, monkeypatch):
    _companion(host, 20.0, 0.5)
    lone = _unforced(system, host, monkeypatch)
    _belt, forced = _members(system, host)
    assert len(forced) == len(lone) and _mean_e(forced) > _mean_e(lone) + 0.005
    for before, after in zip(lone, forced):
        assert {k: v for k, v in before.items() if k not in ("eccentricity", "periapsis_deg")} == {k: v for k, v in after.items() if k not in ("eccentricity", "periapsis_deg")}
        assert 0.0 <= after["eccentricity"] <= belt_roster.E_MAX and 0.0 <= after["periapsis_deg"] < 360.0


def test_a_giant_beside_the_belt_dilutes_the_forcing(system, host):
    _companion(host, 20.0, 0.5)
    belt, bare = _members(system, host)
    bare = [dict(r) for r in bare]
    belt.give_secular_planets([5.2], [JUPITER])
    assert _mean_e(belt.members) < _mean_e(bare)


def test_a_companion_beyond_the_theory_leaves_the_belt_alone(system, host, monkeypatch):
    _companion(host, 5.0, 0.5)               # alpha 0.44 and more
    lone = _unforced(system, host, monkeypatch)
    _belt, rows = _members(system, host)
    assert rows == lone


def test_the_forced_periapsis_follows_the_companions_direction(system, host):
    _companion(host, 200.0, 0.8, mass_ratio=1.0)
    _belt, rows = _members(system, host)
    direction = math.degrees(__import__("planetgen.physics.belt_forcing", fromlist=["x"]).periapsis_rad(host.mass, 1, 200.0))
    spread = [abs((r["periapsis_deg"] - direction + 180.0) % 360.0 - 180.0) for r in rows]
    assert sum(spread) / len(spread) < 90.0


TILT = spin.orbit_normal(20.0, 90.0)             # a companion's plane 20 degrees from the planets'


def _tilted_companion(host, normal, axis_au=20.0, eccentricity=0.3):
    companion = plain_star(host.luminosity)
    companion.mass = 0.5 * host.mass
    rf.attach(host, [rf.Companion(companion, 1, axis_au, eccentricity, normal=normal)])


def _mean_i(rows):
    return sum(r["inclination_deg"] for r in rows) / len(rows)


def test_a_companion_with_no_plane_keeps_every_inclination_and_node(system, host, monkeypatch):
    _companion(host, 20.0, 0.5)
    lone = _unforced(system, host, monkeypatch)
    _belt, rows = _members(system, host)
    assert [(r["inclination_deg"], r["node_deg"]) for r in rows] == [(r["inclination_deg"], r["node_deg"]) for r in lone]


def test_a_coplanar_companion_keeps_every_inclination_and_node(system, host, monkeypatch):
    _tilted_companion(host, (0.0, 0.0, 1.0))
    lone = _unforced(system, host, monkeypatch)
    _belt, rows = _members(system, host)
    assert [(r["inclination_deg"], r["node_deg"]) for r in rows] == [(r["inclination_deg"], r["node_deg"]) for r in lone]


def test_a_tilted_companion_tilts_a_bare_belt_toward_its_plane(system, host, monkeypatch):
    _tilted_companion(host, TILT)
    lone = _unforced(system, host, monkeypatch)
    _belt, rows = _members(system, host)
    assert _mean_i(rows) > _mean_i(lone) + 1.0
    for before, after in zip(lone, rows):
        assert 0.0 <= after["inclination_deg"] <= 180.0 and 0.0 <= after["node_deg"] < 360.0
        assert {k: v for k, v in before.items() if k not in ("inclination_deg", "node_deg", "eccentricity", "periapsis_deg")} == {k: v for k, v in after.items() if k not in ("inclination_deg", "node_deg", "eccentricity", "periapsis_deg")}


def test_a_giant_beside_the_belt_holds_it_flatter_than_a_bare_belt(system, host):
    _tilted_companion(host, TILT)
    belt, bare = _members(system, host)
    bare = [dict(r) for r in bare]
    belt.give_secular_planets([5.2], [JUPITER])
    assert _mean_i(belt.members) < _mean_i(bare)


def test_the_tilt_keeps_the_angle_to_the_forced_axis(system, host, monkeypatch):
    _tilted_companion(host, TILT)
    lone = _unforced(system, host, monkeypatch)
    _belt, rows = _members(system, host)
    assert any(a["inclination_deg"] != b["inclination_deg"] for a, b in zip(lone, rows))
    for after in rows:
        assert math.isclose(vec3.norm(spin.orbit_normal(after["inclination_deg"], after["node_deg"])), 1.0, abs_tol=1e-9)
