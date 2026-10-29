# tests/test_poi_icons.py

"""
MAP.320: the point-of-interest symbols. Seven vector icons, one SVG file each under `static/poi/`, listed once on the Python
side (`web/lib/poi_icons.py`) and once on the map's (`static/poiicons.js`), coloured by the tokens of `static/style.css`.

The files, the two lists and the stylesheet are held to each other here; `tests/test_poi_icons_browser.py` checks what a
browser draws from them.
"""

import colorsys
import html
import os
import re
import xml.etree.ElementTree as ET

import pytest
from flask import render_template_string

from planetgen.physics.habitability_world import TIER_NAMES
from planetgen.web import poi_icon
from planetgen.web.app import create_app
from planetgen.web.lib import poi_icons

STATIC = os.path.join(os.path.dirname(__file__), "..", "html", "static")
SVG = "{http://www.w3.org/2000/svg}"
IDS = ("habitable-world", "uninhabitable-world", "inhabited-world", "gas-giant", "asteroid-field", "sector", "region")
MIN_LIGHTNESS_GAP = 12.0
"""float: How far, in CIE L* (0 to 100), brown and grey sit from each habitability colour, so they differ in lightness
from all four of them (the Blue and Green a viewer who cannot tell hues apart may take for each other included)."""


def _file(icon_id):
    return os.path.join(STATIC, "poi", f"{icon_id}.svg")


def _css():
    with open(os.path.join(STATIC, "style.css"), encoding="utf-8") as handle:
        return handle.read()


def _js_icons():
    """`[(id, label, token or None)]` of the `POI_ICONS` list of `static/poiicons.js`."""
    with open(os.path.join(STATIC, "poiicons.js"), encoding="utf-8") as handle:
        text = handle.read()
    body = text[text.index("export const POI_ICONS = ["):text.index("];")]
    return [(icon_id, label, None if token == "null" else token.strip('"'))
            for icon_id, label, token in re.findall(r'id: "([a-z-]+)", label: "([^"]+)", token: (null|"[^"]+")', body)]


def _themes():
    """`{"light": tokens, "dark": tokens, "os-dark": tokens}` of the hex custom properties each theme block defines."""
    css = _css()
    light = re.search(r"^:root \{(.*?)^\}", css, re.S | re.M).group(1)
    dark = re.search(r':root\[data-theme="dark"\] \{(.*?)^\}', css, re.S | re.M).group(1)
    auto = re.search(r':root:not\(\[data-theme="light"\]\) \{(.*?)^  \}', css, re.S | re.M).group(1)

    def parse(block):
        return dict(re.findall(r"--([a-z-]+):\s*(#[0-9a-fA-F]{6})\b", block))
    return {"light": parse(light), "dark": parse(dark), "os-dark": parse(auto)}


def _linear(channel):
    return channel / 12.92 if channel <= 0.04045 else ((channel + 0.055) / 1.055) ** 2.4


def _luminance(colour):
    red, green, blue = (_linear(int(colour[index:index + 2], 16) / 255) for index in (1, 3, 5))
    return 0.2126 * red + 0.7152 * green + 0.0722 * blue


def _lightness(colour):
    """CIE L* of a `#rrggbb` colour."""
    y = _luminance(colour)
    return 116 * y ** (1 / 3) - 16 if y > 0.008856 else 903.3 * y


def _contrast(a, b):
    high, low = sorted((_luminance(a), _luminance(b)), reverse=True)
    return (high + 0.05) / (low + 0.05)


def _hue_and_saturation(colour):
    red, green, blue = (int(colour[index:index + 2], 16) / 255 for index in (1, 3, 5))
    hue, _lightness_, saturation = colorsys.rgb_to_hls(red, green, blue)
    return hue * 360, saturation


# --- The files ---------------------------------------------------------------------------------------------------

def test_there_is_one_file_per_symbol_and_no_other():
    assert tuple(icon.id for icon in poi_icons.ICONS) == IDS
    assert sorted(os.listdir(os.path.join(STATIC, "poi"))) == sorted(f"{icon_id}.svg" for icon_id in IDS)


@pytest.mark.parametrize("icon_id", IDS)
def test_each_file_is_a_plain_square_silhouette_that_names_itself(icon_id):
    root = ET.parse(_file(icon_id)).getroot()
    assert root.tag == f"{SVG}svg" and root.get("viewBox") == "0 0 24 24" and root.get("id") == poi_icons.SPRITE_ID
    assert root.find(f"{SVG}title").text == poi_icons.BY_ID[icon_id].label
    shapes = [child for child in root if child.tag != f"{SVG}title"]
    assert shapes and all(child.tag == f"{SVG}path" and child.get("d") for child in shapes)
    # The paint is the page's (`.poi-icon`, a token's colour), never the file's: a colour or a style written into the
    # file would ignore the theme, and the Content-Security-Policy refuses an inline style besides.
    text = open(_file(icon_id), encoding="utf-8").read()
    assert not re.search(r"\b(fill|stroke|style|class|opacity)=", text.replace('fill-rule="evenodd"', ""))
    assert "<style" not in text and "<script" not in text and "#" not in text.replace('id="poi"', "")


# --- The two lists -------------------------------------------------------------------------------------------------

def test_the_python_and_the_map_lists_name_the_same_symbols_in_the_same_order():
    assert [(icon.id, icon.label) for icon in poi_icons.ICONS] == [(icon_id, label) for icon_id, label, _token in _js_icons()]
    scored = [icon.scored for icon in poi_icons.ICONS]
    assert scored == [token is None for _icon_id, _label, token in _js_icons()]     # no token of its own: the score colours it


def test_every_token_the_map_names_is_a_property_of_the_stylesheet():
    css = _css()
    for _icon_id, _label, token in _js_icons():
        if token:
            assert re.search(rf"^\s*{re.escape(token)}:\s*#", css, re.M), token
    with open(os.path.join(STATIC, "poiicons.js"), encoding="utf-8") as handle:
        tiers = re.search(r"TIER_TOKENS = \[(.*?)\]", handle.read()).group(1)
    assert re.findall(r'"(--tier-[a-z]+)"', tiers) == [f"--tier-{name.lower()}" for name in TIER_NAMES]


def test_the_stylesheet_colours_every_symbol_and_every_tier():
    css = _css()
    for icon in poi_icons.ICONS:
        if not icon.scored and icon.id not in ("sector", "region"):
            assert f".poi-icon.poi-{icon.id} {{" in css, icon.id
    for tier in poi_icons.TIER_CLASSES:
        assert f".poi-icon.{tier} {{" in css, tier


# --- The tokens ------------------------------------------------------------------------------------------------------

@pytest.mark.parametrize("theme", ["light", "dark", "os-dark"])
def test_brown_and_grey_clear_three_to_one_on_the_page_and_stand_off_every_habitability_colour(theme):
    colours = _themes()[theme]
    for token in ("poi-brown", "poi-grey", "poi-gas"):
        for surface in ("bg", "bg-elevated"):
            assert _contrast(colours[token], colours[surface]) >= 3.0, (theme, token, surface)    # WCAG 1.4.11
    for token in ("poi-brown", "poi-grey"):
        for name in TIER_NAMES:
            gap = abs(_lightness(colours[token]) - _lightness(colours[f"tier-{name.lower()}"]))
            assert gap >= MIN_LIGHTNESS_GAP, (theme, token, name, gap)


@pytest.mark.parametrize("theme", ["light", "dark", "os-dark"])
def test_brown_is_brown_and_never_red_and_grey_is_grey(theme):
    colours = _themes()[theme]
    hue, saturation = _hue_and_saturation(colours["poi-brown"])
    assert 20 <= hue <= 50 and saturation >= 0.3, (theme, hue, saturation)        # between red (0) and yellow (60), muted
    assert _hue_and_saturation(colours["poi-grey"])[1] <= 0.2, theme


def test_both_dark_blocks_carry_the_same_tokens():
    themes = _themes()
    assert themes["dark"] == themes["os-dark"]


# --- The markup ------------------------------------------------------------------------------------------------------

def _static_url(path):
    return f"/static/{path}?v=test"


def test_a_symbol_is_an_svg_that_names_itself_and_points_into_its_file():
    markup = poi_icons.symbol_html("gas-giant", _static_url)
    assert markup == ('<svg class="poi-icon poi-gas-giant" role="img" aria-label="Gas giant" data-poi="gas-giant" '
                      'focusable="false"><use href="/static/poi/gas-giant.svg?v=test#poi"></use></svg>')
    root = ET.fromstring(markup.replace("<svg ", f'<svg xmlns="{SVG[1:-1]}" ', 1))
    assert root.get("role") == "img" and root.get("aria-label") == "Gas giant"


def test_a_world_the_score_colours_takes_its_tier_and_says_it():
    for tier, name in enumerate(TIER_NAMES):
        markup = poi_icons.symbol_html("inhabited-world", _static_url, tier)
        assert f'poi-tier-{tier}' in markup and f'aria-label="Inhabited world: {name}"' in markup
    assert "poi-tier-1" in poi_icons.symbol_html("habitable-world", _static_url)     # no tier yet: Green
    assert "poi-tier" not in poi_icons.symbol_html("sector", _static_url, 2)          # a symbol with a token of its own ignores it
    assert "poi-tier" not in poi_icons.symbol_html("uninhabitable-world", _static_url, 3)    # brown, never red


def test_an_unknown_symbol_or_tier_is_refused():
    with pytest.raises(ValueError):
        poi_icons.symbol_html("star", _static_url)
    with pytest.raises(ValueError):
        poi_icons.symbol_html("habitable-world", _static_url, 4)
    with pytest.raises(ValueError):
        poi_icons.symbol_html("habitable-world", _static_url, -1)


def test_the_template_global_and_the_key_draw_every_symbol_once_and_the_scored_ones_in_four_colours():
    with create_app().test_request_context("/"):
        assert 'href="/static/poi/region.svg?v=' in str(poi_icon("region"))
        key = render_template_string('{% from "partials/poi_key.html" import poi_key with context %}{{ poi_key("h") }}')
    labels = re.findall(r'aria-label="([^"]+)"', key)
    expected = []
    for icon in poi_icons.ICONS:
        expected += [f"{icon.label}: {name}" for name in TIER_NAMES] if icon.scored else [icon.label]
    assert labels == [label for label in labels if label in expected] and sorted(labels) == sorted(expected)
    for icon in poi_icons.ICONS:
        assert f'id="poi-key-{icon.id}"' in key and html.escape(icon.meaning) in key and icon.shape in key
