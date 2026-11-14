# tests/test_gen1620_class_request.py

"""GEN.1620: `planetgen system --force-class LETTER` and `--force-life-class ID` ask for a class (reusing FORCE_CLASS and a new FORCE_LIFE_CLASS),
the result says what came of it, and with neither asked a seeded system is byte-identical."""

import argparse
import json

import pytest

from planetgen import describe
from planetgen.generation import run_system
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.util import draw


def _system(config, seed=11):
    with draw.bound(seed):
        return StarSystem(system_config=config)


def _args(**given):
    base = dict(star_type="G2V", num_orbits=None, name=None, age=None, system_file=None, stars=None, max_stars=None,
                flavor_chance_system=None, flavor_chance_planet=None, max_planet_flavor=False, debug_areas=None,
                force_class=None, force_life_class=None)
    for name, _attr, _description in run_system.TRISTATE_OPTIONS:
        base[name] = None
    base.update(given)
    return argparse.Namespace(**base)


def test_with_no_request_a_seeded_system_is_byte_identical():
    plain = SystemConfig()
    plain.STAR_TYPE = "G2V"
    built = run_system.build_system_config(_args())
    assert built.FORCE_CLASS is None and built.FORCE_LIFE_CLASS is None
    assert json.dumps(describe.describe(_system(plain))) == json.dumps(describe.describe(_system(built)))


def test_the_options_set_the_config():
    built = run_system.build_system_config(_args(force_class="o", force_life_class="6"))
    assert (built.FORCE_CLASS, built.FORCE_LIFE_CLASS) == ("O", "6")


def test_a_life_class_is_found_by_number_or_by_number_and_variant():
    config = SystemConfig()
    config.STAR_TYPE = "G2V"
    system = _system(config)
    held = [(body.life_class, body.life_class_variant) for planet in system.all_planets for body in (planet, *getattr(planet, "moons", ()))
            if getattr(body, "life_class", None) is not None]
    if not held:
        pytest.skip("this seed holds no life")
    number, variant = held[0]
    assert system.bodies_with_life_class(str(number))
    assert system.bodies_with_life_class(f"{number}{variant or ''}")
    assert not system.bodies_with_life_class("0z")


def test_the_report_says_what_came_of_the_request():
    config = SystemConfig()
    config.STAR_TYPE = "G2V"
    system = _system(config)
    assert run_system.class_request_report(system, config) is None
    classes = {getattr(body, "planet_class", None) for planet in system.all_planets for body in (planet, *getattr(planet, "moons", ()))}
    classes.discard(None)
    config.FORCE_CLASS = sorted(classes)[0]
    assert "holds it" in run_system.class_request_report(system, config)
    config.FORCE_CLASS = next(one for one in "ABCDEFGHIJKLMNOPRSTUVWXYZ" if one not in classes)
    assert "no body of this system came out as class" in run_system.class_request_report(system, config)
    config.FORCE_CLASS, config.FORCE_LIFE_CLASS = None, "9"
    assert "Life class 9 asked for" in run_system.class_request_report(system, config)


def test_a_life_class_that_does_not_come_out_is_an_unmet_requirement():
    config = SystemConfig()
    config.STAR_TYPE = "M5V"
    config.FORCE_LIFE_CLASS = "9"
    system = _system(config)
    assert system.bodies_with_life_class("9") or "a world with life class 9" in system.unmet_requirements


def test_the_cli_refuses_a_class_that_does_not_exist():
    assert run_system.class_request_problem("O", "4c") is None and run_system.class_request_problem(None, "6") is None
    assert "not a planet class" in run_system.class_request_problem("Q", None)
    assert "not a life class" in run_system.class_request_problem(None, "12")
