"""
The description of each planet class (UX.252, GEN.1212, `planetgen.physics.planet_class_text`).
"""

import re

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
        planet_class_text.description("Q")


def test_descriptions_follow_the_plain_words_guide():
    for letter, text in planet_class_text.DESCRIPTIONS.items():
        sentence = f"It is {text}."
        assert plainwords.check_sentence(sentence) == [], (letter, sentence)


def test_the_asteroid_type_has_a_description_and_no_letter():
    for moon in (False, True):
        text = planet_class_text.description(None, moon=moon)
        assert text.startswith("an asteroid")
        assert plainwords.check_sentence(f"It is {text}.") == []


def test_a_class_has_one_row_with_its_phrase_and_a_sentence_without_a_number():
    """GEN.1432: the phrase after "It is" and the class page sentence are in one row of planet_classes.json; the sentence holds no number, because the page lists the conditions from class_fit under it."""
    from planetgen.physics import planet_classes
    for cls in planet_classes.CLASSES:
        assert cls.phrase and cls.phrase == planet_class_text.DESCRIPTIONS[cls.letter]
        assert not any(character.isdigit() for character in cls.sentence), (cls.letter, cls.sentence)
    assert not any(character.isdigit() for character in planet_classes.ASTEROID.sentence)
    assert planet_class_text.ASTEROID_PHRASE == planet_classes.ASTEROID.phrase


def test_no_class_text_states_a_feature_of_one_planetary_system():
    """GEN.1432: rings on a gas giant, a tilted field on an ice giant, Venus-like and Titan-like are Solar System facts; the words say what the conditions say."""
    from planetgen.physics import planet_classes
    texts = [t for cls in planet_classes.CLASSES for t in (cls.phrase, cls.sentence)] + [planet_classes.ASTEROID.phrase]
    for text in texts:
        assert not re.search(r"\b(rings?|tilted|magnetic|Venus|Titan|Mars|Jupiter|Saturn)\b", text), text


def test_a_phrase_or_sentence_that_breaks_the_rules_is_refused(tmp_path, monkeypatch):
    import json
    from planetgen.physics import planet_classes
    good = json.load(open(planet_classes.DATA_FILE, encoding="utf-8"))
    for change in (lambda d: d["rows"][0].update(sentence="A class of 3 Earth masses."), lambda d: d["rows"][0].update(phrase=""),
                   lambda d: d["rows"][0].update(phrase="Upper case")):
        document = json.loads(json.dumps(good))
        change(document)
        path = tmp_path / "bad.json"
        path.write_text(json.dumps(document), encoding="utf-8")
        monkeypatch.setattr(planet_classes, "DATA_FILE", str(path))
        with pytest.raises(planet_classes.ClassError):
            planet_classes._load()
