# tests/test_static_css_wellformed.py

"""
The site's stylesheet is one file that many branches append to; a bad merge leaves a conflict marker or an unclosed block
that silently swallows every rule after it. This test catches both.
"""

import pathlib
import re

STYLE = pathlib.Path(__file__).resolve().parents[1] / "html" / "static" / "style.css"


def _code():
    return re.sub(r"/\*.*?\*/", "", STYLE.read_text(encoding="utf-8"), flags=re.S)


def test_the_stylesheet_has_no_merge_conflict_markers():
    assert not re.search(r"^(<{7}|={7}|>{7})", _code(), flags=re.M)


def test_the_stylesheet_closes_every_block_it_opens():
    code = _code()
    depth = 0
    for char in code:
        depth += (char == "{") - (char == "}")
        assert depth >= 0, "a closing brace with nothing open"
    assert depth == 0, f"{depth} block(s) never closed"
