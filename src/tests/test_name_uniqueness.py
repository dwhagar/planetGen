# tests/test_name_uniqueness.py

"""
Pure unit tests for `stellarObjects.nameUniqueness` -- no database, no
generation, just the collision-resolution state machines and their
inverse (`strip_decoration`). See that module's own docstring for the
full sector > system > planet/moon hierarchy these feed into; the actual
database integration (`_db.py`'s `_reserve_*`/`_confirm_*` functions) is
covered separately by `test_db_persistence.py`.

Run with: pytest tests/test_name_uniqueness.py
"""

from stellarObjects import nameUniqueness as nu
from stellarObjects.names import (
    COMPANION_SUFFIXES, DIMINUTIVE_PREFIXES, GREEK_LETTERS, ROMAN_NUMERAL_VALUES,
)

BASE = "Sol"


# ---------------------------------------------------------------------------
# resolve_greek_roman_collision
# ---------------------------------------------------------------------------

def test_no_collision_returns_the_bare_base_name():
    new_name, rename = nu.resolve_greek_roman_collision(BASE, 0)
    assert new_name == BASE
    assert rename is None


def test_first_collision_renames_the_existing_row_to_alpha_and_new_row_to_beta():
    new_name, rename = nu.resolve_greek_roman_collision(BASE, 1)
    assert new_name == f"{GREEK_LETTERS[1]} {BASE}"  # "Beta Sol"
    assert rename == (BASE, f"{GREEK_LETTERS[0]} {BASE}")  # ("Sol", "Alpha Sol")


def test_third_collision_uses_the_next_greek_letter_with_no_rename():
    new_name, rename = nu.resolve_greek_roman_collision(BASE, 2)
    assert new_name == f"{GREEK_LETTERS[2]} {BASE}"  # "Gamma Sol"
    assert rename is None


def test_greek_tier_runs_through_every_letter_with_no_further_renames():
    for existing_count in range(2, len(GREEK_LETTERS)):
        new_name, rename = nu.resolve_greek_roman_collision(BASE, existing_count)
        assert new_name == f"{GREEK_LETTERS[existing_count]} {BASE}"
        assert rename is None
    # The 24th (last) Greek letter, Omega, is reached at existing_count == 23.
    new_name, _rename = nu.resolve_greek_roman_collision(BASE, len(GREEK_LETTERS) - 1)
    assert new_name == f"Omega {BASE}"


def test_greek_tier_exhaustion_converts_alpha_to_roman_numeral_one_and_two():
    n_greek = len(GREEK_LETTERS)
    new_name, rename = nu.resolve_greek_roman_collision(BASE, n_greek)
    alpha = f"Alpha {BASE}"
    assert rename == (alpha, f"{alpha} I")
    assert new_name == f"{alpha} II"


def test_roman_tier_advances_through_every_remaining_value_with_no_rename():
    n_greek = len(GREEK_LETTERS)
    alpha = f"Alpha {BASE}"
    # existing_count == n_greek + 1 mints roman index 2 ("III"), since the
    # n_greek transition step above already minted indices 0 and 1 in one go.
    for offset, expected_index in enumerate(range(2, len(ROMAN_NUMERAL_VALUES)), start=1):
        new_name, rename = nu.resolve_greek_roman_collision(BASE, n_greek + offset)
        assert rename is None
        assert new_name == f"{alpha} {nu.ROMAN_NUMERALS_BY_VALUE[ROMAN_NUMERAL_VALUES[expected_index]]}"


def test_exhausting_every_slot_reports_none_and_none():
    # The Alpha row is renamed into (rather than adding a distinct new
    # row for) the roman tier's first slot, so the last row actually
    # created lands one short of the full 24-Greek + 12-roman label
    # count -- GREEK_ROMAN_CAPACITY - 1 rows total, existing_count 0..
    # (GREEK_ROMAN_CAPACITY - 2) inclusive.
    last_valid_count = nu.GREEK_ROMAN_CAPACITY - 2
    new_name, rename = nu.resolve_greek_roman_collision(BASE, last_valid_count)
    assert new_name == f"Alpha {BASE} C"  # the last (12th) roman slot

    new_name, rename = nu.resolve_greek_roman_collision(BASE, last_valid_count + 1)
    assert new_name is None
    assert rename is None


def test_capacity_matches_greek_plus_roman_list_lengths():
    assert nu.GREEK_ROMAN_CAPACITY == len(GREEK_LETTERS) + len(ROMAN_NUMERAL_VALUES)


# ---------------------------------------------------------------------------
# resolve_diminutive / resolve_companion -- identical shape, different lists
# ---------------------------------------------------------------------------

def test_resolve_diminutive_first_use_returns_the_first_prefix():
    prefix, next_index = nu.resolve_diminutive(None)
    assert prefix == DIMINUTIVE_PREFIXES[0]
    assert next_index == 0


def test_resolve_diminutive_advances_on_repeat():
    prefix, next_index = nu.resolve_diminutive(0)
    assert prefix == DIMINUTIVE_PREFIXES[1]
    assert next_index == 1


def test_resolve_diminutive_exhausts_after_the_last_prefix():
    last_index = len(DIMINUTIVE_PREFIXES) - 1
    prefix, next_index = nu.resolve_diminutive(last_index)
    assert prefix is None
    assert next_index is None


def test_resolve_companion_first_use_returns_the_first_suffix():
    suffix, next_index = nu.resolve_companion(None)
    assert suffix == COMPANION_SUFFIXES[0]
    assert next_index == 0


def test_resolve_companion_advances_on_repeat():
    suffix, next_index = nu.resolve_companion(0)
    assert suffix == COMPANION_SUFFIXES[1]
    assert next_index == 1


def test_resolve_companion_exhausts_after_the_last_suffix():
    last_index = len(COMPANION_SUFFIXES) - 1
    suffix, next_index = nu.resolve_companion(last_index)
    assert suffix is None
    assert next_index is None


# ---------------------------------------------------------------------------
# strip_decoration -- inverse of all three resolvers above
# ---------------------------------------------------------------------------

def test_strip_decoration_is_a_no_op_on_an_undecorated_name():
    assert nu.strip_decoration(BASE) == BASE


def test_strip_decoration_is_a_no_op_on_a_two_word_undecorated_name():
    # Sector names are always two words (generate_sector_name) -- neither
    # word should ever falsely match a decoration token.
    assert nu.strip_decoration("Vorhibe Drabria") == "Vorhibe Drabria"


def test_strip_decoration_strips_every_greek_roman_variant():
    count = 0
    while True:
        new_name, _rename = nu.resolve_greek_roman_collision(BASE, count)
        if new_name is None:
            break
        assert nu.strip_decoration(new_name) == BASE, new_name
        count += 1
    assert count == nu.GREEK_ROMAN_CAPACITY - 1


def test_strip_decoration_strips_a_diminutive_prefix():
    prefix, _index = nu.resolve_diminutive(None)
    assert nu.strip_decoration(f"{prefix} {BASE}") == BASE


def test_strip_decoration_strips_a_companion_suffix():
    suffix, _index = nu.resolve_companion(None)
    assert nu.strip_decoration(f"{BASE} {suffix}") == BASE
