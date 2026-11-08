# tests/test_class_code_filter.py

"""GEN.1208: the class filters of the database queries read a written class code (`class_code.parse_code`): a letter alone matches all its digits, a letter and a digit one position."""

from planetgen.db import query


def _clause(alias, codes):
    clauses, params = [], []
    query._class_clause(alias, codes, clauses, params)
    return clauses, params


def test_letters_alone_filter_as_they_always_did():
    assert _clause("p", {"O", "M"}) == (["p.planet_class IN (?,?)"], ["M", "O"])
    assert _clause("p", set()) == ([], [])


def test_a_letter_and_a_digit_filter_that_position_only():
    clauses, params = _clause("m", {"F2"})
    assert clauses == ["(m.planet_class = ? AND COALESCE(m.class_variant, 0) = ?)"]
    assert params == ["F", 2]
    # Digit 0 is the ideal of the family and a row with no stored digit has it.
    assert _clause("p", {"M0"})[1] == ["M", 0]


def test_letters_and_codes_together_match_either():
    clauses, params = _clause("p", {"M", "F2"})
    assert clauses == ["(p.planet_class IN (?) OR (p.planet_class = ? AND COALESCE(p.class_variant, 0) = ?))"]
    assert params == ["M", "F", 2]


def test_a_tag_that_is_not_a_class_code_matches_nothing_as_before():
    assert _clause("p", {"F22"}) == (["p.planet_class IN (?)"], ["F22"])
