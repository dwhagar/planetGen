"""GEN.1220: the class code is a letter plus a digit; digit 0 is not written."""

import pytest

from planetgen.physics import class_code, planet_classes


def test_digit_zero_is_not_written():
    assert class_code.format_code("M", 0) == "M" and class_code.format_code("M", None) == "M"
    assert class_code.format_code("F", 2) == "F2"
    assert class_code.format_code(None, 0) is None


@pytest.mark.parametrize("text, expected", [("M", ("M", 0)), ("F2", ("F", 2)), ("f4", ("F", 4)), (" y9 ", ("Y", 9))])
def test_a_code_parses_to_its_letter_and_digit(text, expected):
    assert class_code.parse_code(text) == expected


@pytest.mark.parametrize("text", ["", "M10", "1M", "MM", None])
def test_a_bad_code_is_refused(text):
    with pytest.raises(class_code.ClassCodeError):
        class_code.parse_code(text)


def test_a_digit_outside_zero_to_nine_is_refused():
    with pytest.raises(class_code.ClassCodeError):
        class_code.format_code("M", 10)


def test_format_and_parse_round_trip_for_every_letter_and_digit():
    for letter in class_code.LETTERS:
        for digit in class_code.DIGITS:
            assert class_code.parse_code(class_code.format_code(letter, digit)) == (letter, digit)


def test_codes_sort_by_letter_then_digit():
    codes = [("M", 0), ("F", 2), ("F", 0), ("A", 0)]
    assert sorted(codes, key=lambda code: class_code.sort_key(*code)) == [("A", 0), ("F", 0), ("F", 2), ("M", 0)]


def test_the_new_letters_are_classes_now_and_the_confirmed_meanings_hold():
    """Boss, 11:59 AM EDT 2026-10-31: Y sulfuric acid, U methane-ethane, X ammonia-water, Z hot thin-air, R Hothouse (classes since GEN.1207)."""
    classes = {cls.letter for cls in planet_classes.CLASSES}
    for letter in ("R", "Y", "U", "X", "Z"):
        assert letter in classes and letter not in planet_classes.RESERVE
    assert "acid" in class_code.LETTERS["Y"] and "methane" in class_code.LETTERS["U"] and "ammonia" in class_code.LETTERS["X"]
    assert "Hothouse" in class_code.LETTERS["R"] and "thin" in class_code.LETTERS["Z"]
