# tests/test_class_after_bands.py

"""
GEN.1587: the class of a body is judged when its latitude bands exist. The classes that read the hottest and coldest region, the night side and the wet bulb (M, O, H, P and R,
`class_domains.derived`) saw no bands at build time and so never matched: a request for one of them could not be met, and a moon was classed before its water was remade at the final age.
"""

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.planet import Planet
from planetgen.generation.star import Star
from planetgen.physics import constants, planet_classes, planets as planet_physics
from planetgen.util import draw


def _star():
    config = SystemConfig()
    config.STAR_TYPE = "G2V"
    config.WEIRDNESS_DIAL = 0.0
    with draw.bound("class-after-bands-star"):
        return config, Star(config)


def _earth_like(seed):
    config, star = _star()
    inner, outer = star.habitable_zone
    with draw.bound(seed):
        return Planet(config, star, star.habitable_zone, 0.5 * (inner + outer), mass=constants.EARTH_MASS_TO_KG, moon_count=0)


def test_a_rocky_planet_is_classed_with_its_latitude_bands_made(monkeypatch):
    seen = []
    real = planet_classes.classify_detail

    def recording(body):
        seen.append((body.body_type, getattr(body, "band1_mean_k", None)))
        return real(body)

    monkeypatch.setattr(planet_classes, "classify_detail", recording)
    for seed in range(4):
        seen.clear()
        _earth_like(seed)
        assert seen and all(band is not None for _type, band in seen), seen        # every classification of the build sees the bands


def test_a_class_that_reads_the_bands_is_met_by_a_body_as_it_is_built(monkeypatch):
    """The Earth twin without bands is F at match 0, with bands M at 1 (GEN.1587): a class that reads the bands must be what the build gives, or a request for it can never be met."""
    real = planet_classes.classify_detail

    def reads_the_bands(body):
        detail = real(body)
        if detail is None or getattr(body, "band1_mean_k", None) is None:
            return detail
        return type(detail)(**{**detail.__dict__, "letter": "M"}) if hasattr(detail, "__dict__") else detail._replace(letter="M")

    monkeypatch.setattr(planet_classes, "classify_detail", reads_the_bands)
    planet = _earth_like(0)
    assert planet.planet_class == "M"


def test_the_class_a_body_is_built_with_is_the_class_of_its_finished_conditions():
    for seed in range(4):
        planet = _earth_like(seed)
        assert planet.planet_class == planet_classes.classify(planet)


def test_a_body_that_is_prepared_and_not_classed_has_no_class_and_one_that_is_classed_has_one():
    planet = _earth_like(0)
    letter = planet.planet_class
    planet_physics.prepare_class(planet)
    assert planet.planet_class is None and planet.class_match is None
    planet_physics.classify(planet)
    assert planet.planet_class == letter


def test_a_moon_is_classed_after_its_water_is_remade_at_the_final_age(monkeypatch):
    from planetgen.generation.system import StarSystem

    events = []
    real_refresh, real_class = planet_physics.refresh_at_final_age, planet_physics.apply_class
    monkeypatch.setattr(planet_physics, "refresh_at_final_age", lambda body: (events.append(("refresh", id(body))), real_refresh(body))[1])
    monkeypatch.setattr(planet_physics, "apply_class", lambda body: (events.append(("class", id(body))), real_class(body))[1])
    config = SystemConfig()
    config.PLANETS = True
    moons = []
    for seed in range(40):
        events.clear()
        with draw.bound(seed):
            system = StarSystem(system_config=config)
        moons = [moon for planet in system.planets for moon in getattr(planet, "moons", ())]
        if moons:
            break
    assert moons, "no moon in the first 40 seeds"
    for moon in moons:
        sequence = [kind for kind, ident in events if ident == id(moon)]
        assert "class" in sequence and "refresh" in sequence
        last_refresh = max(i for i, kind in enumerate(sequence) if kind == "refresh")
        last_class = max(i for i, kind in enumerate(sequence) if kind == "class")
        assert last_class > last_refresh, sequence      # the class is the last word on a moon whose water was remade
