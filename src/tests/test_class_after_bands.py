# tests/test_class_after_bands.py

"""
GEN.1587: a planet or moon is classed after the latitude bands exist and after its water is made again at the final age. The classes that read the hottest and coldest region, the
night side and the wet bulb (M, O, H, P and R) see the bands that `generate_spin` ends by making (`habitability_world.generate`); classed before them they score 0, so a request for
one of them could never be met and a moon's class was judged on the water it was drawn with, not the water it ends with.
"""

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.planet import Planet
from planetgen.generation.system import StarSystem
from planetgen.physics import planets as planet_physics
from planetgen.util import draw

from tests.fuzz_support import first_seed_with


def _star():
    cfg = SystemConfig()
    cfg.STAR_TYPE, cfg.PLANETS, cfg.BINARY_SYSTEM = "G2V", False, False
    with draw.bound(1):
        return cfg, StarSystem(system_config=cfg).star


def test_a_planet_built_for_a_request_is_classed_after_its_spin_made_its_bands(monkeypatch):
    calls = []
    for name in ("generate_spin", "classify", "prepare_class"):
        original = getattr(planet_physics, name)
        monkeypatch.setattr(planet_physics, name, lambda *args, _o=original, _n=name, **kwargs: (calls.append(_n), _o(*args, **kwargs))[1])
    cfg, star = _star()
    with draw.bound(7):
        planet = Planet(cfg, star, star.habitable_zone, 1.0, moon_count=0, mass=5.972e24)
    assert calls.count("classify") == 1 and calls.index("prepare_class") < calls.index("generate_spin") < calls.index("classify")
    assert planet.band1_mean_k is not None            # the bands were there when the class was judged
    assert planet.planet_class == planet_physics.planet_classes.classify(planet)


def test_a_moon_is_classed_after_its_water_is_made_again_at_the_final_age(monkeypatch):
    log = []
    refresh, apply_class = planet_physics.refresh_at_final_age, planet_physics.apply_class
    monkeypatch.setattr(planet_physics, "refresh_at_final_age", lambda body: (log.append(("refresh", id(body))), refresh(body))[1])
    monkeypatch.setattr(planet_physics, "apply_class", lambda body: (log.append(("class", id(body))), apply_class(body))[1])

    def build(seed):
        cfg = SystemConfig()
        cfg.STAR_TYPE, cfg.PLANETS, cfg.MAX_PLANETS, cfg.MOONS, cfg.BINARY_SYSTEM = "G2V", True, True, True, False
        with draw.bound(seed):
            return StarSystem(system_config=cfg)

    _seed, system = first_seed_with(build, lambda s: any(m.body_type != 'a' for p in s.planets for m in getattr(p, "moons", ())), seeds=range(1, 40), what="a system with a moon")
    moons = [moon for planet in system.planets for moon in getattr(planet, "moons", ())]
    assert moons
    for moon in moons:
        kinds = [kind for kind, ident in log if ident == id(moon)]
        if "refresh" in kinds:
            assert kinds.index("class") > kinds.index("refresh"), "the moon was classed before its water was made again"
