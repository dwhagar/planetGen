# tests/test_sky_help.py

"""
VIEW.38: the help page about the sky behind a system (`/help/sky`). It renders, is in the help registry and the Help index, reads the numbers
that live in Python from the code that uses them, keeps the plain-words style guide, and is linked from the System View's sky switch.
"""

from planetgen.api import system_sky
from planetgen.galaxy import sky, sky_band
from planetgen.util import plainwords
from planetgen.web.lib import help_registry, sky_help
from planetgen.web.lib.skypack import MAX_LIMIT, MIN_LIMIT
from tests.test_class_pages import _html, app, client  # noqa: F401


def _texts(page):
    found = list(page["intro"])
    for section in page["sections"]:
        found += section["paragraphs"]
    return found


def test_the_sky_help_page_renders_and_is_in_the_help_index(client):
    html = _html(client, "/help/sky")
    assert 'id="sky-help"' in html and "The sky behind a system" in html
    for section in sky_help.page()["sections"]:
        assert f'id="{section["id"]}"' in html
    assert "/help/sky" in _html(client, "/help")
    assert "sky_help_page" in {topic.endpoint for topic in help_registry.topics()}


def test_the_numbers_of_the_page_are_the_codes(client):
    text = " ".join(_texts(sky_help.page()))
    assert f"{sky.DEFAULT_LIMIT:g}" in text
    assert f"from {MIN_LIMIT:g} to {MAX_LIMIT:g}" in text
    assert f"{sky_band.CELL_PC:g} parsec cell" in text
    assert f"{system_sky.MAX_SPRITES:,}" in text


def test_the_page_keeps_the_plain_words_style_guide():
    broken = [(text, rule) for text in _texts(sky_help.page()) for sentence in plainwords.sentences(text)
              for rule, _detail in plainwords.check_sentence(sentence)]
    assert broken == []
    assert all("\u2014" not in text for text in _texts(sky_help.page()))


def test_the_page_links_to_the_stand_in_registry(client):
    html = _html(client, "/help/sky")
    assert "/references/stand-ins" in html
