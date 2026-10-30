"""
GEN.206: a planet's description and evolutionary timeline never claim more life than its conditions allow. No supported life class
means no class, timeline or claim of life; a class whose complexity index stays under the multicellular threshold means a timeline that
stops before any complex flora or fauna.
"""
import re

import pytest

from planetgen import tuning
from planetgen.describe import build as describe_build, life_facts
from planetgen.describe.quantity import inline_text
from planetgen.generation import life as planetLife
from planetgen.generation.config import SystemConfig
from planetgen.generation.star import Star
from planetgen.physics import life_complexity
from planetgen.generation.system import StarSystem

from tests.fuzz_support import deterministic_entropy
from tests import class_support, render_support

LIFE_WORDS = re.compile(r"life|bacteri|vegetation|flora|fauna|evolution|alga|microb|organism|biosphere|civiliz", re.I)
COMPLEX_WORDS = re.compile(r"vegetation|macro-algae|moss|ferns|forests|sapient|civilization|invertebrate", re.I)

# Planets drawn in the ecosphere ("e" zone), the only place a timeline is drawn; whatever class each comes out as (GEN.646).
ECOSPHERE_SEEDS = range(6)


@pytest.fixture(scope="module")
def host_star():
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    return Star(cfg)


def _planet(host_star, seed=0):
    inner, outer = host_star.habitable_zone
    return class_support.draw_planet(host_star, (inner + outer) / 2.0, seed)


def _description_and_life_text(planet):
    """Everything a reader sees about the world's nature and life, except the index's own explanation
    (which says what the scores are)."""
    return render_support.prose_text(describe_build._planet_section(planet))


def test_a_star_too_young_for_a_crust_and_oceans_has_no_life_anywhere(host_star, monkeypatch):
    """GEN.997: the class letter no longer bars the young star's worlds, so the life pass does: under `LIFE_MIN_STAR_AGE_GY` no body holds a class,
    though the same world around an old star does."""
    for seed in range(60):
        planet = _planet(host_star, seed)
        planetLife.apply_life_data(planet)
        if planet.life_class is not None:
            break
    else:
        pytest.skip("no world of the first 60 seeds is alive around this star")
    monkeypatch.setattr(host_star, "age", tuning.LIFE_MIN_STAR_AGE_GY * 0.5)
    planetLife.apply_life_data(planet)
    assert planet.life_class is None and not planet.has_timeline
    for moon in planet.moons:
        planetLife.apply_life_data(moon)
        assert moon.life_class is None


@pytest.mark.parametrize("seed", ECOSPHERE_SEEDS)
def test_a_world_no_class_supports_makes_no_claim_of_life(host_star, seed, monkeypatch):
    monkeypatch.setattr(planetLife.life_classes, "score_classes", lambda *args, **kwargs: {})
    planet = _planet(host_star, seed)
    planet.moons = []  # a moon is scored and described on its own
    planetLife.apply_life_data(planet)
    assert planet.life_class is None
    assert not planet.has_timeline and planet.evolutionary_pace is None
    assert planet.reflection_spectrum_visible is None and planet.reflection_spectrum_non_visible is None
    assert not LIFE_WORDS.search(planet.description), planet.description
    text = _description_and_life_text(planet)
    # The planet's table row gives the class and the like; the sentences must not talk of life.
    sentences = " ".join(s for s in re.split(r"(?<=[.!?])\s+", text) if not s.startswith("{{"))
    assert not LIFE_WORDS.search(sentences), sentences


@pytest.mark.parametrize("seed", ECOSPHERE_SEEDS)
def test_a_world_whose_class_stays_microbial_stops_its_timeline_before_complex_life(host_star, seed, monkeypatch):
    # the complexity index is set by hand: a class that never passes the complex-cell threshold
    monkeypatch.setattr(planetLife.life_complexity, "assess", lambda *args, **kwargs: life_complexity.Complexity(1e-6, 1.0, 1.0, 0, "photosynthesis"))
    for attempt in range(4):  # the pace and age draws differ; the cap holds on every one
        planet = _planet(host_star, 10 * seed + attempt)
        planetLife.apply_life_data(planet)
        assert planet.life_stage in (None, "abiogenesis", "photosynthesis")
        assert not COMPLEX_WORDS.search(planet.description), planet.description
        life = life_facts.resolve(planet)
        timeline_text = " ".join(inline_text(sentence) for sentence in describe_build._timeline_sentences(life)) if life and life.timeline else ""
        assert not COMPLEX_WORDS.search(timeline_text), timeline_text


def test_no_class_description_claims_life_or_names_a_planet():
    """GEN.221: the class text is the same whatever the index allows, so it claims no life, no flora, no percentage and
    no "planet" (a moon reads it as it stands); the life sentences say what the index allows."""
    claims = re.compile(r"life|bacteri|vegetation|flora|Earth-like|planet|%")
    for cls, found in tuning.PLANET_CLASSES.items():
        assert not claims.search(found["description"]), (cls, found["description"])


def test_a_class_description_stands_whatever_the_conditions_allow(host_star, monkeypatch):
    for supported in (True, False):
        if not supported:
            monkeypatch.setattr(planetLife.life_classes, "score_classes", lambda *args, **kwargs: {})
        for is_moon in (False, True):
            planet = _planet(host_star)
            planet.is_moon = is_moon
            planetLife.apply_life_data(planet)
            stands = tuning.PLANET_CLASSES[planet.planet_class]["description"]
            if is_moon:
                stands = stands.replace("planets", "moons").replace("planet", "moon")
            assert planet.description == stands


def test_a_body_with_no_scores_is_left_to_its_class(host_star):
    planet = _planet(host_star)
    planet.phi_bio = planet.phi_cpx = None
    planetLife.apply_life_data(planet)
    assert planet.description == tuning.PLANET_CLASSES[planet.planet_class]["description"]


def test_no_generated_description_contradicts_its_own_habitability_index():
    """A sample of generated systems: wherever the index allows less life than the text would claim, it doesn't."""
    lifeless = microbial_only = 0
    for seed in range(40):
        with deterministic_entropy(1000 + seed):
            system = StarSystem(system_config=SystemConfig())
        for planet in system.planets + system.secondary_planets:
            for body in [planet] + list(getattr(planet, "moons", ())):
                if body.body_type == "a" or getattr(body, "phi_bio", None) is None:
                    continue
                if body.life_class is None:
                    lifeless += 1
                    assert not body.has_timeline, body.name
                    assert not LIFE_WORDS.search(body.description), (body.name, body.description)
                elif body.life_complexity < tuning.LIFE_COMPLEXITY["thresholds"]["multicellularity"]:
                    microbial_only += 1
                    assert body.life_stage in (None, "abiogenesis", "photosynthesis", "complex_cells")
                    assert not COMPLEX_WORDS.search(body.description), (body.name, body.description)
    assert lifeless > 0 and microbial_only >= 0
