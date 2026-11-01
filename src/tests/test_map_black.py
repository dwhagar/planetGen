# tests/test_map_black.py

"""Every map is drawn on black in either theme (Boss, 2026-10-31): the viewports' background is the `--map-bg` token (black), they take the dark tokens so their
labels and controls are light on black, and no scene paints a navy of its own over it."""

import re
from pathlib import Path

from planetgen.web import STATIC_DIR

CSS = (Path(STATIC_DIR) / "style.css").read_text(encoding="utf-8")
VIEWPORTS = (".starmap-viewport", ".galaxymap-viewport", ".phenomenonmap-viewport", ".phenomrender-viewport", ".navmap-viewport", ".sysview-viewport")


def _rule(selector, wanted):
    """The body of the rule for `selector` that has the property `wanted` (a selector can be the last of a list as well)."""
    for match in re.finditer(r"(?m)^" + re.escape(selector) + r" \{\n(.*?)\n\}", CSS, re.S):
        if wanted in match.group(1):
            return match.group(1)
    raise AssertionError(f"{selector} has no {wanted}")


def test_the_map_token_is_black_and_every_viewport_paints_it():
    assert re.search(r"--map-bg: #000000;", CSS)
    for selector in VIEWPORTS:
        assert _rule(selector, "background: var(--map-bg);"), selector


def _tokens(block):
    return dict(re.findall(r"--([a-z-]+):\s*([^;]+);", block))


def test_the_viewports_take_the_dark_tokens_and_the_text_color_from_them():
    dark = re.search(r':root\[data-theme="dark"\] \{\n(.*?)^\}', CSS, re.S | re.M).group(1)
    maps = re.search(r"(?m)^\.starmap-viewport,\n\.galaxymap-viewport,\n\.phenomenonmap-viewport,\n\.phenomrender-viewport,\n\.navmap-viewport \{\n(.*?)^\}", CSS, re.S).group(1)
    assert _tokens(maps) == _tokens(dark) and len(_tokens(maps)) > 15
    assert "color: var(--text);" in re.search(r"(?m)^\.starmap-viewport,\n\.galaxymap-viewport,\n\.phenomenonmap-viewport,\n\.phenomrender-viewport,\n\.navmap-viewport \{\n  color: var\(--text\);\n\}", CSS).group(0)


def test_no_map_surface_keeps_a_navy_of_its_own():
    assert "#05070c" not in CSS.lower() and "#05070d" not in CSS.lower()
    for name in ("nebulaview.js", "phenomenonrender.js"):
        text = (Path(STATIC_DIR) / name).read_text(encoding="utf-8")
        assert "scene.background = new THREE.Color(0x000000)" in text, name
    assert "#05070c" not in (Path(__file__).resolve().parents[1] / "planetgen" / "web" / "maps" / "phenomenonrender.py").read_text(encoding="utf-8").lower()


def test_the_lines_of_a_map_use_the_map_line_token_and_the_body_view_is_on_black():
    for selector in (".galaxymap-ring", ".galaxymap-axis", ".navmap-compass-ring", ".sysmap-orbit"):
        assert _rule(selector, "stroke: var(--map-line);"), selector
    assert _rule(".body-view-canvas", "background: var(--map-bg);")
