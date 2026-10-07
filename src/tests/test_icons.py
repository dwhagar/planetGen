# tests/test_icons.py

"""The site's icon sprite (`static/icons.svg`) and `icon(name)` (UX.28)."""

import os
import re
import xml.etree.ElementTree as ET

import pytest

from planetgen.web import ICON_NAMES, icon
from planetgen.web.app import create_app

SPRITE = os.path.join(os.path.dirname(__file__), "..", "html", "static", "icons.svg")


def test_the_sprite_holds_exactly_the_named_icons():
    root = ET.parse(SPRITE).getroot()
    symbols = root.findall("{http://www.w3.org/2000/svg}symbol")
    assert {symbol.get("id") for symbol in symbols} == ICON_NAMES
    assert all(symbol.get("viewBox") == "0 0 24 24" for symbol in symbols)


def test_icon_is_hidden_markup_pointing_into_the_sprite():
    with create_app().test_request_context("/"):
        markup = str(icon("edit"))
        with pytest.raises(ValueError):
            icon("no-such-icon")
    assert re.fullmatch(r'<svg class="icon" aria-hidden="true" focusable="false">'
                        r'<use href="/static/icons\.svg\?v=[^"]+#edit"></use></svg>', markup)
