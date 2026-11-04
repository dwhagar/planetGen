"""
The words of the life classes (GEN.892, `planetgen.physics.life_class_text`).
"""

import json
import re

import pytest

from planetgen.physics import life_class_text, life_classes

NUMBERS = sorted({c.number for c in life_classes.CLASSES})


def _every_string(node):
    if isinstance(node, str):
        yield node
    elif isinstance(node, dict):
        for value in node.values():
            yield from _every_string(value)
    elif isinstance(node, list):
        for value in node:
            yield from _every_string(value)


def test_every_class_has_words_for_its_first_two_rungs_and_the_generic_rungs_cover_the_rest():
    for number in NUMBERS:
        own = life_class_text.DOCUMENT["classes"][str(number)]
        assert "1" in own and "2" in own
    for rung in life_class_text.RUNGS:
        for number in NUMBERS:
            name, description = life_class_text.stage(number, rung)
            assert name and description and "{n}" not in description


def test_a_generic_rung_names_the_class():
    assert "Class 7" in life_class_text.stage_description("7", 3)
    assert "Class 2" in life_class_text.stage_description("2", 5)


def test_a_class_with_words_of_its_own_uses_them_and_the_variant_letter_makes_no_difference():
    assert life_class_text.stage_name("1a", 2) == life_class_text.stage_name("1c", 2) == "Oxygenating mats"
    assert life_class_text.stage_name("4c", 5) == "Vent and seep animals"


def test_animal_grades_have_their_own_names_and_keep_the_class_sentence():
    assert [life_class_text.stage_name("4", 5, grade=g) for g in life_class_text.GRADES] == \
        ["Sessile animals", "Small mobile animals", "Large active animals"]
    assert "tubeworms" in life_class_text.stage_description("4", 5, grade="A1")
    assert life_class_text.stage_description("1", 5, grade="A3").startswith("Animals have appeared among Class 1 life.")


@pytest.mark.parametrize("call", [lambda: life_class_text.stage("10", 1), lambda: life_class_text.stage("x", 1),
                                  lambda: life_class_text.stage("1", 7), lambda: life_class_text.stage("1", 4, grade="A1"),
                                  lambda: life_class_text.stage("1", 5, grade="A9"), lambda: life_class_text.lines("0")])
def test_asking_for_something_that_does_not_exist_is_refused(call):
    with pytest.raises(life_class_text.LifeClassTextError):
        call()


def test_the_four_lines_are_the_rule_files_and_a_bare_number_reads_its_first_variant():
    said = life_class_text.lines("7")
    assert said["drinks"] == "liquid methane and ethane" and set(said) == {"drinks", "eats", "breathes", "gives_off"}
    assert life_class_text.lines("4") == life_class_text.lines("4a")
    assert life_class_text.lines("4b")["breathes"] != life_class_text.lines("4a")["breathes"]


def test_the_sentence_reads_in_plain_words():
    assert life_class_text.sentence("7") == \
        "It drinks liquid methane and ethane, eats hydrogen and acetylene (hydrogenation), breathes nothing and gives off methane."


def test_the_text_file_follows_the_house_rules():
    with open(life_class_text.DATA_FILE, encoding="utf-8") as handle:
        strings = list(_every_string(json.load(handle)))
    assert len(strings) > 60
    for text in strings:
        assert "—" not in text and "–" not in text, text
        assert not re.search(r"\b[A-Z]{2,5}\.\d+\b", text), f"a TODO code in shown text: {text}"
        assert text == text.strip() and "  " not in text, text


def test_a_text_file_with_a_missing_part_is_refused():
    document = json.loads(json.dumps(life_class_text.DOCUMENT))
    del document["classes"]["9"]
    with pytest.raises(life_class_text.LifeClassTextError):
        life_class_text._check(document)
    document = json.loads(json.dumps(life_class_text.DOCUMENT))
    document["classes"]["1"]["4"]["description"] = ""
    with pytest.raises(life_class_text.LifeClassTextError):
        life_class_text._check(document)
    document = json.loads(json.dumps(life_class_text.DOCUMENT))
    del document["grades"]["A2"]
    with pytest.raises(life_class_text.LifeClassTextError):
        life_class_text._check(document)


def test_a_stage_sentence_names_the_whole_class_not_its_family():
    """UX.251: "Class 3a life has grown ...", never "Class 3 life" for a class with letters; a bare number stays a number."""
    said = [life_class_text.stage_description("3a", rung) for rung in life_class_text.RUNGS]
    assert any("Class 3a life" in text for text in said) and not any("Class 3 life" in text for text in said)
    assert any("Class 2 life" in life_class_text.stage_description("2", rung) for rung in life_class_text.RUNGS)
    assert any("Class 3 life" in life_class_text.stage_description(3, rung) for rung in life_class_text.RUNGS)
    assert "Class 3a life" in life_class_text.stage_description("3a", 5, "A1")
