"""
GEN.206: a planet's description and evolutionary timeline never claim more life than its habitability index
allows. No microbial life (`phi_bio` shows 0.00) means no chemistry, timeline or claim of life; microbial life
only (`phi_cpx` shows 0.00) means a timeline that stops before any complex flora or fauna.
"""
import re

import pytest

from planetgen import tuning
from planetgen.generation import life as planetLife
from planetgen.generation.config import SystemConfig
from planetgen.generation.evolution import life_stage_from_paragraphs
from planetgen.generation.planet import Planet
from planetgen.generation.star import Star
from planetgen.generation.system import StarSystem

from tests.fuzz_support import deterministic_entropy

LIFE_WORDS = re.compile(r"life|bacteri|vegetation|flora|fauna|evolution|alga|microb|organism|biosphere|civiliz", re.I)
COMPLEX_WORDS = re.compile(r"vegetation|macro-algae|moss|ferns|forests|sapient|civilization|invertebrate", re.I)

# Habitable classes that can sit in the ecosphere ("e" zone), the only place a timeline is drawn.
ECOSPHERE_HABITABLE = [cls for cls in tuning.HABITABLE_PLANET_CLASSES if tuning.PLANET_CLASSES[cls]["e"]]


@pytest.fixture(scope="module")
def host_star():
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    return Star(cfg)


def _planet(host_star, cls):
    inner, outer = host_star.habitable_zone
    return Planet(SystemConfig(), host_star, host_star.habitable_zone, (inner + outer) / 2.0, planet_class=cls)


def _description_and_life_text(planet):
    """Everything a reader sees about the world's nature and life, except the index's own explanation
    (which says what the scores are)."""
    habitability = set(planet._habitability_paragraphs())
    paragraphs = [p for p in planet.to_paragraph_list() if p not in habitability and not p.lstrip().startswith(("=", "#", "{"))]
    return " ".join(paragraphs)


@pytest.mark.parametrize("cls", ECOSPHERE_HABITABLE)
def test_a_world_the_index_gives_no_microbial_life_makes_no_claim_of_life(host_star, cls):
    planet = _planet(host_star, cls)
    planet.moons = []  # a moon is scored and described on its own
    planet.phi_bio, planet.phi_cpx = 0.0, 0.0
    planetLife.apply_life_data(planet)
    assert planet.life_chemical is None
    assert planet.evolutionary_data == [] and planet.evolutionary_speed is None
    assert planet.reflection_spectrum_visible is None and planet.reflection_spectrum_non_visible is None
    assert not LIFE_WORDS.search(planet.description), planet.description
    text = _description_and_life_text(planet)
    # The planet's table row gives the class and the like; the sentences must not talk of life.
    sentences = " ".join(s for s in re.split(r"(?<=[.!?])\s+", text) if not s.startswith("{{"))
    assert not LIFE_WORDS.search(sentences), sentences


@pytest.mark.parametrize("cls", ECOSPHERE_HABITABLE)
def test_a_world_the_index_allows_only_microbial_life_stops_its_timeline_before_complex_life(host_star, cls):
    for _ in range(12):  # the pace and age draws differ; the cap holds on every one
        planet = _planet(host_star, cls)
        planet.phi_bio, planet.phi_cpx = 0.6, 0.0
        planetLife.apply_life_data(planet)
        assert life_stage_from_paragraphs(planet.evolutionary_data) in (None, "abiogenesis", "photosynthesis")
        assert not COMPLEX_WORDS.search(planet.description), planet.description
        assert not COMPLEX_WORDS.search(" ".join(planet.evolutionary_data)), planet.evolutionary_data


def test_a_class_description_that_claims_life_is_restored_when_the_index_allows_it(host_star):
    planet = _planet(host_star, "L")
    planet.phi_bio, planet.phi_cpx = 0.0, 0.0
    planetLife.apply_life_data(planet)
    assert planet.description == tuning.NO_LIFE_CLASS_DESCRIPTIONS["L"]
    planet.phi_bio, planet.phi_cpx = 0.5, 0.0
    planetLife.apply_life_data(planet)
    assert planet.description == tuning.MICROBIAL_ONLY_CLASS_DESCRIPTIONS["L"]
    planet.phi_bio, planet.phi_cpx = 0.5, 0.5
    planetLife.apply_life_data(planet)
    assert planet.description == tuning.PLANET_CLASSES["L"]["description"]


def test_a_moon_keeps_its_own_wording_when_the_index_gives_no_life(host_star):
    planet = _planet(host_star, "G")
    planet.is_moon = True
    planet.phi_bio, planet.phi_cpx = 0.0, 0.0
    planetLife.apply_life_data(planet)
    assert planet.description == tuning.NO_LIFE_CLASS_DESCRIPTIONS["G"]  # no "planet" in it to reword
    assert planet.life_chemical is None


def test_a_body_with_no_scores_is_left_to_its_class(host_star):
    planet = _planet(host_star, "M")
    planet.phi_bio = planet.phi_cpx = None
    assert planetLife.life_ceiling(planet) is None
    planetLife.apply_life_data(planet)
    assert planet.description == tuning.PLANET_CLASSES["M"]["description"]


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
                ceiling = planetLife.life_ceiling(body)
                if ceiling == "none":
                    lifeless += 1
                    assert body.life_chemical is None and body.evolutionary_data == [], body.name
                    assert not LIFE_WORDS.search(body.description), (body.name, body.description)
                elif ceiling == "microbial":
                    microbial_only += 1
                    assert life_stage_from_paragraphs(body.evolutionary_data) in (None, "abiogenesis", "photosynthesis")
                    assert not COMPLEX_WORDS.search(body.description), (body.name, body.description)
    assert lifeless > 0 and microbial_only >= 0
