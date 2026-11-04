# tests/test_rewrite_bibliographies.py

"""
`scripts/rewrite_bibliographies.py` (DOC.98 follow-up): an APA-form bibliography entry of an imported note is read, tied to a source of the
reference list by its first author, year and title words, and written as `[R0123]`. No network is used here.
"""

import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

import rewrite_bibliographies as rb  # noqa: E402

ENTRY = "* Example, A. B., & Other, C. (2011). Report on the tidal heating of icy moons: a worked case. _Journal of Examples_, 109(2), 101-135."


def test_an_apa_entry_is_read_into_its_first_author_year_and_rest():
    m = rb.APA.match(ENTRY)
    assert (m.group("lead"), m.group("first"), m.group("year")) == ("* ", "Example", "2011")
    assert m.group("rest").startswith("Report on the tidal heating")


def test_a_sentence_is_not_an_entry():
    assert rb.APA.match("The tidal heating was computed in 2011 (see above).") is None
    assert rb.APA.match("* [R0123]") is None


def test_title_words_ignore_case_accents_and_short_words():
    assert rb.title_words("Évolution of the Stars, a Study") == {"evolution", "stars", "study"}


def test_a_record_covers_an_entry_when_most_of_its_title_is_in_it():
    rest = rb.APA.match(ENTRY).group("rest")
    assert rb.covers(rest, "Report on the tidal heating of icy moons: a worked case")
    assert rb.covers(rest, "Report on the tidal heating of icy moons")
    assert not rb.covers(rest, "Thermal management of lunar landers and sublimators")


def test_a_short_title_never_ties_by_accident():
    assert not rb.covers(ENTRY, "Moons")
