"""
ASCII-printable regression tests for stellarObjects/names.py.

Every generated name ultimately gets pasted into a wiki page and stored as
TEXT in the database -- it should never depend on non-ASCII characters
rendering correctly somewhere downstream. This walks every public
list-of-strings constant in names.py and asserts each string is 7-bit
printable ASCII, catching any future addition automatically instead of
relying on manual review each time a list is extended.

Run with: pytest tests/test_names.py
"""
import pytest

from stellarObjects import names

# Every public module-level list-of-strings constant in names.py that feeds
# generated names. Enumerated explicitly (rather than introspecting the
# module) so an unrelated new list constant doesn't silently get skipped or
# unexpectedly included.
NAME_LIST_CONSTANTS = [
    "BAD_CONSONANTS",
    "STAR_NAMES",
    "PLANET_NAMES",
    "MOON_NAMES",
    "STAR_PREFIXES",
    "STAR_SUFFIXES",
    "PLANET_PREFIXES",
    "PLANET_SUFFIXES",
    "MOON_PREFIXES",
    "MOON_SUFFIXES",
    "SECTOR_NAMES",
    "SECTOR_PREFIXES",
    "SECTOR_SUFFIXES",
    "UNIVERSAL_PHONEMES",
]


def is_ascii_printable(s):
    return all(32 <= ord(c) < 127 for c in s)


@pytest.mark.parametrize("constant_name", NAME_LIST_CONSTANTS)
def test_name_list_is_ascii_printable(constant_name):
    string_list = getattr(names, constant_name)
    offenders = [s for s in string_list if not is_ascii_printable(s)]
    assert not offenders, (
        f"{constant_name} contains non-ASCII-printable strings: {offenders!r}"
    )


def test_name_list_constants_cover_every_list_of_strings_in_module():
    """Guards against a future new list constant being forgotten above."""
    actual_list_constants = {
        attr for attr in vars(names)
        if attr.isupper()
        and isinstance(getattr(names, attr), list)
        and all(isinstance(item, str) for item in getattr(names, attr))
    }
    assert actual_list_constants == set(NAME_LIST_CONSTANTS)
