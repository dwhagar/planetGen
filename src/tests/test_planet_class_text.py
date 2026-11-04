"""
The description of each planet class (UX.252, GEN.1212, `planetgen.physics.planet_class_text`).
"""

import pytest

from planetgen import tuning
from planetgen.physics import planet_class_text
from planetgen.util import plainwords


def test_every_planet_class_has_a_description():
    assert sorted(planet_class_text.DESCRIPTIONS) == sorted(tuning.PLANET_CLASSES)


def test_a_moon_reads_moon_where_the_text_says_planet():
    for letter, text in planet_class_text.DESCRIPTIONS.items():
        moon_text = planet_class_text.description(letter, moon=True)
        assert "planet" not in moon_text.replace("planetary", "")
        if "planet" not in text:
            assert moon_text == text


def test_an_unknown_class_is_refused():
    with pytest.raises(planet_class_text.PlanetClassTextError):
        planet_class_text.description("Z")


def test_descriptions_follow_the_plain_words_guide():
    for letter, text in planet_class_text.DESCRIPTIONS.items():
        sentence = f"It is {text}."
        assert plainwords.check_sentence(sentence) == [], (letter, sentence)


def test_the_asteroid_type_has_a_description_and_no_letter():
    for moon in (False, True):
        text = planet_class_text.description(None, moon=moon)
        assert text.startswith("an asteroid")
        assert plainwords.check_sentence(f"It is {text}.") == []
