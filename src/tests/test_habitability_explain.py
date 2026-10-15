"""
UX.91: a planet's and moon's description explains its PHI-4 colours: the score and the equipment, then each of
the four factors with its colour, this world's input values and which input sets the colour. The explanation reads
the same world and band edges as the stored score, so its colours always agree with the chip.
"""

import random

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import habitability as hab
from planetgen.physics import habitability_explain as explain


def _bold(text):
    return f"**{text}**"


def _earth():
    return hab.World(pressure_kpa=101.0, temperature_c=15.0, gases={"O2": 0.21, "CO2": 0.0004}, gravity_ms2=9.8,
                     life_stage="biosphere")


def test_an_earth_like_world_is_one_paragraph_of_six_blue_categories_with_their_values():
    paragraphs = explain.paragraphs(_earth(), _bold)
    assert len(paragraphs) == 1
    text = paragraphs[0]
    assert text.startswith("**Habitability** is 1.00 (Blue).") and text.endswith("A visitor needs no equipment.")
    for domain in ("Pressure", "Temperature", "Air", "Water", "Radiation", "Suitability"):
        assert f"**{domain}: Blue**." in text
    assert "At 101 kPa no gear is needed." in text and "At 15 C no gear is needed." in text


def test_a_factor_names_the_input_that_sets_its_colour_and_where_it_sits():
    world = hab.World(pressure_kpa=101.0, temperature_c=15.0, gases={"O2": 0.21, "CO2": 0.03}, gravity_ms2=9.8)
    air = next(f for f in explain.factors(world) if f["domain"] == "Air")
    assert air["tier"] == "Yellow" and air["limiting"]["label"] == "carbon dioxide"
    assert "outside the Green band (up to 2 kPa), inside the Yellow band (up to 5 kPa)" in air["limiting"]["phrase"]
    text = explain.paragraphs(world, _bold)[0]
    assert "a mask with a scrubber" in text and "**Air: Yellow**" in text


def test_a_planet_of_a_compact_host_says_its_dose_is_rated_lethal():
    world = hab.World(pressure_kpa=101.0, temperature_c=15.0, gases={"O2": 0.21}, gravity_ms2=9.8,
                      surface_dose_sv_yr=1000.0)
    text = explain.paragraphs(world, _bold, compact_host="neutron star")[0]
    assert "1,000 Sv a year" in text and "A planet of a neutron star is rated as a lethal dose" in text
    assert "**Radiation: Red**" in text


def _worlds(count=400):
    rng = random.Random(91)
    for _ in range(count):
        o2 = rng.choice([0.0, 0.05, 0.21, 0.4, 0.9])
        co2 = rng.choice([0.0, 0.0004, 0.01, 0.05, 0.5])
        yield hab.World(
            pressure_kpa=10 ** rng.uniform(-1, 5), temperature_c=rng.uniform(-150, 500),
            gases={"O2": o2, "CO2": co2, "CO": rng.choice([0.0, 1e-6, 1e-3]), "H2S": rng.choice([0.0, 1e-6, 1e-3])},
            gravity_ms2=rng.uniform(1, 25), relative_humidity=rng.uniform(0.05, 0.99), solvent="water",
            ph=rng.choice([None, 2.0, 7.0, 10.0, 13.0]), water_activity=rng.choice([1.0, 0.8, 0.5]),
            chaotropicity_kj_kg=rng.choice([0.0, 50.0, 90.0]), surface_dose_sv_yr=10 ** rng.uniform(-3, 3),
            life_stage=rng.choice(hab.LIFE_STAGES), water_source=rng.choice(["liquid", "ice", "hydrated", "vapour", None]),
            life_on_surface=rng.random() < 0.8)


def test_the_explained_colours_always_equal_the_stored_ones():
    for world in _worlds():
        stored = hab.phi4(world)
        for factor in explain.factors(world):
            score, colour = stored[factor["domain"]]
            assert (factor["score"], factor["tier"]) == (score, colour)
            assert factor["limiting"]["score"] == pytest.approx(min(i["score"] for i in factor["items"]))
            # The domain's own score is the worst of its inputs (the temperature input may be capped at Yellow).
            assert factor["score"] == pytest.approx(factor["limiting"]["score"], abs=1e-3), factor["domain"]


def _system_text(markdown):
    for seed in range(1, 60):
        random.seed(seed)
        config = SystemConfig()
        config.MARKDOWN = markdown
        text = str(StarSystem(config))
        if "Habitability** is" in text or "Habitability''' is" in text:
            return text
    raise AssertionError("no rocky body in 60 systems")


def test_a_generated_planet_or_moon_description_carries_the_explanation():
    text = _system_text(True)
    assert text.count("**Habitability** is") >= 1
    for domain in ("Pressure", "Temperature", "Air", "Water", "Radiation", "Suitability"):
        assert f"**{domain}: " in text


def test_the_wikitext_page_uses_wikitext_emphasis():
    text = _system_text(False)
    assert "'''Habitability''' is" in text and "**Habitability" not in text
