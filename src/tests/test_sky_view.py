# tests/test_sky_view.py

"""GEN.1443: the bright sources in a body's sky (`planetgen.physics.sky_view`) and the sky text of the description (`describe.build._sky_block`)."""

import math
from types import SimpleNamespace

import pytest

from planetgen.describe import build
from planetgen.physics import constants, daylight, radiation_field as rf, sky_suns, sky_view

SUN_W = constants.SOLAR_LUMINOSITY


def star(luminosity_w=SUN_W, temperature=5772.0, radius_km=6.957e5, state="main-sequence", kind="G2V"):
    return SimpleNamespace(luminosity=luminosity_w, temperature=temperature, radius=radius_km, state=state, type=kind, mass=1.989e30)


def with_companion(host, other, axis_au):
    rf.attach(host, [rf.Companion(other, 1, axis_au, 0.0, (axis_au, 0.0, 0.0))])
    return host


def test_a_lone_star_gives_the_host_alone_and_no_text():
    host = star()
    (body,) = sky_view.bodies(host, 1.0)
    assert body.role == "host" and body.share == pytest.approx(1.0)
    assert body.angular_deg == pytest.approx(0.53, abs=0.02)         # the Sun from Earth
    assert build._sky_sentences([body]) == []


def test_a_bright_companion_is_a_second_sun_that_makes_a_day():
    host = with_companion(star(), star(luminosity_w=0.5 * SUN_W, temperature=4200.0, kind="K5V"), 6.0)
    host_body, second = sky_view.bodies(host, 1.0)
    assert second.role == "sun" and second.light_w_m2 >= daylight.DAY_SUN_W_M2
    assert host_body.share + second.share == pytest.approx(1.0)
    (sentence,) = build._sky_sentences([host_body, second])
    assert "second sun" in sentence and "K-type star" in sentence and "a day of its own" in sentence


def test_a_dim_companion_lights_the_night_and_a_dimmer_one_does_not_show():
    sun = with_companion(star(), star(luminosity_w=1.0e-2 * SUN_W), 10.0)             # 0.14 W/m2: over the twilight of a day
    assert [b.role for b in sky_view.bodies(sun, 1.0)] == ["host", "sun"]
    night = with_companion(star(), star(luminosity_w=1.3e-2 * SUN_W), 60.0)           # 5 mW/m2: two full Moons, under the day
    host_body, light = sky_view.bodies(night, 1.0)
    assert light.role == "night light" and sky_view.NIGHT_LIGHT_W_M2 <= light.light_w_m2 < daylight.DAY_SUN_W_M2
    (sentence,) = build._sky_sentences([host_body, light])
    assert "lights its nights" in sentence and "full Moon" in sentence
    gone = with_companion(star(), star(luminosity_w=1.0e-4 * SUN_W), 60.0)            # 0.04 mW/m2: not shown
    assert [b.role for b in sky_view.bodies(gone, 1.0)] == ["host"]


def test_a_white_dwarf_companion_is_named_by_its_kind_and_a_white_dwarf_host_is_named():
    host = with_companion(star(), star(luminosity_w=0.01 * SUN_W, temperature=12000.0, radius_km=7.0e3, state="white-dwarf", kind="DA"), 2.0)
    second = sky_view.bodies(host, 0.8)[1]
    assert second.kind == "white dwarf" and second.angular_deg < 0.5
    lone = star(luminosity_w=0.003 * SUN_W, temperature=9000.0, radius_km=7.0e3, state="white-dwarf", kind="DA")
    (body,) = sky_view.bodies(lone, 0.05)
    (sentence,) = build._sky_sentences([body])
    assert sentence.startswith("Its sun is a white dwarf")


def test_the_light_is_the_one_the_temperature_already_reads():
    host = with_companion(star(), star(luminosity_w=0.3 * SUN_W), 5.0)
    bodies = sky_view.bodies(host, 1.3)
    assert sum(b.light_w_m2 for b in bodies) == pytest.approx(rf.starlight_w_m2(host, 1.3), rel=1e-9)
    assert bodies[1].light_w_m2 == pytest.approx(sky_suns.light_of(host, rf.companions_of(host)[0], 1.3))


def test_the_sky_text_of_a_moon_follows_its_planets_distance():
    host = with_companion(star(), star(luminosity_w=0.5 * SUN_W), 6.0)
    planet = SimpleNamespace(star=host, distance=1.0, eccentricity=0.0)
    moon = SimpleNamespace(star=host, distance=0.002, eccentricity=0.0)
    assert build._sky_block(moon, planet)["text"] == build._sky_block(planet, None)["text"]
    assert build._sky_block(SimpleNamespace(star=star(), distance=1.0), None)["text"] == []


def _sky_texts(document):
    found = []
    stack = list(document["sections"]) if isinstance(document, dict) and "sections" in document else list(document)
    while stack:
        node = stack.pop()
        for block in node.get("blocks", ()):
            if block["type"] == "prose" and block["role"] == "sky":
                found.append(" ".join(str(part) for sentence in block["text"] for part in sentence))
        stack.extend(node.get("children", ()))
    return found


def test_a_generated_white_dwarf_system_names_its_sun_on_every_planet():
    from planetgen.describe import describe
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from planetgen.util import draw
    draw.set_run_seed(1443)
    for _ in range(6):
        config = SystemConfig()
        config.BINARY_SYSTEM = False
        config.STAR_TYPE = "M2VII"
        system = StarSystem(system_config=config)
        if system.planets:
            texts = _sky_texts(describe(system))
            assert texts and all("Its sun is a white dwarf" in text for text in texts)
            return
    pytest.skip("no white dwarf system with planets in six draws")
