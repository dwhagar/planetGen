# tests/test_habitability_symbols.py

"""The PHI-6 habitability symbols (UX.101): the sprite shapes, the helper's markup, every rating colour in both themes."""

import os
import re
import xml.etree.ElementTree as ET

import pytest

from planetgen.physics.habitability_world import DOMAINS, TIER_NAMES
from flask import render_template_string

from planetgen.web import habitability_symbol, habitability_symbol_row
from planetgen.web.app import create_app
from planetgen.web.lib import habitability_symbols as symbols

STATIC = os.path.join(os.path.dirname(__file__), "..", "html", "static")
SVG = "{http://www.w3.org/2000/svg}"
SPRITE_URL = "/static/icons.svg?v=test"


def _parse(markup):
    """The markup as XML in the SVG namespace, as a browser reads inline `<svg>`."""
    return ET.fromstring(re.sub(r"^<(\w+)", rf'<\1 xmlns="{SVG[1:-1]}"', markup, count=1))


def _sprite():
    return {symbol.get("id"): symbol for symbol in ET.parse(os.path.join(STATIC, "icons.svg")).getroot()
            .findall(f"{SVG}symbol")}


def _css():
    with open(os.path.join(STATIC, "style.css"), encoding="utf-8") as handle:
        return handle.read()


def test_there_are_seven_symbols_six_domains_and_the_overall():
    assert symbols.SYMBOL_NAMES == (*DOMAINS, "overall")
    assert set(symbols.LABELS) == set(symbols.SHAPES) == set(symbols.SPRITE_IDS) == set(symbols.SYMBOL_NAMES)


def test_every_symbol_is_in_the_sprite_as_an_outline_drawn_in_the_text_colour():
    sprite = _sprite()
    for name, sprite_id in symbols.SPRITE_IDS.items():
        symbol = sprite[sprite_id]
        assert symbol.get("viewBox") == "0 0 24 24", name
        shapes = list(symbol)
        assert shapes, name
        for shape in shapes:
            assert shape.get("fill") in (None, "currentColor"), name
            assert shape.get("stroke") is None, f"{name} sets its own stroke, which would ignore the rating colour"
            assert "#" not in "".join(shape.attrib.values()), f"{name} has a fixed colour"


def test_the_shapes_are_the_ones_boss_drew():
    sprite = _sprite()

    def tags(name):
        return [shape.tag.replace(SVG, "") for shape in sprite[symbols.SPRITE_IDS[name]]]
    assert tags("radiation") == ["circle", "path"]
    assert sprite["hab-radiation"].find(f"{SVG}path").get("d").count("M") == 3  # three arcs
    assert tags("temperature") == ["path", "path"]  # a triangle and its vertical line
    assert tags("pressure") == ["rect", "path", "path"]  # a square and two triangles
    assert tags("air") == ["path", "path"]  # a hexagon and three lines
    assert sprite["hab-air"].findall(f"{SVG}path")[1].get("d").count("M") == 3
    assert tags("water") == ["path", "circle"]  # a downward triangle and a circle
    assert tags("suitability") == ["circle", "path", "circle"]  # a circle, a diamond and a dot
    assert tags("overall") == ["circle", "path", "circle"]  # a ring, a triangle and a solid circle
    assert sprite["hab-overall"].findall(f"{SVG}circle")[1].get("fill") == "currentColor"


@pytest.mark.parametrize("name", symbols.SYMBOL_NAMES)
@pytest.mark.parametrize("tier", range(4))
def test_a_symbol_carries_its_rating_in_its_class_and_in_words(name, tier):
    markup = symbols.symbol_html(name, tier, SPRITE_URL)
    words = f"{symbols.LABELS[name]}: {TIER_NAMES[tier]}"
    root = _parse(markup)
    assert root.get("class") == f"icon hab-symbol tier-{TIER_NAMES[tier].lower()}"
    assert root.get("role") == "img"
    assert root.get("aria-label") == words
    assert root.find(f"{SVG}title").text == words
    assert root.find(f"{SVG}use").get("href") == f"{SPRITE_URL}#hab-{name}"


def test_detail_follows_the_rating_and_is_escaped():
    markup = symbols.symbol_html("pressure", 2, SPRITE_URL, '101 kPa <b>"x"</b>')
    assert _parse(markup).get("aria-label") == 'Pressure: Yellow, 101 kPa <b>"x"</b>'
    assert "<b>" not in markup


@pytest.mark.parametrize("name,tier", [("noise", 0), ("pressure", 4), ("pressure", -1), ("pressure", None),
                                        ("pressure", True), ("pressure", "Blue")])
def test_a_name_or_rating_the_symbols_do_not_have_is_refused(name, tier):
    with pytest.raises(ValueError):
        symbols.symbol_html(name, tier, SPRITE_URL)


def _body(**tiers):
    return {"phi4": 0.5, **{f"tier_{domain}": tiers.get(domain, 0) for domain in DOMAINS}}


def test_a_row_is_the_overall_symbol_then_the_six_domains_in_order():
    body = _body(radiation=3, temperature=2)
    body["phi4"] = 0.2  # Red
    row = _parse(symbols.symbol_row(body, SPRITE_URL))
    labels = [svg.get("aria-label") for svg in row.findall(f"{SVG}svg")]
    assert labels == ["Overall habitability: Red", "Pressure: Blue", "Temperature: Yellow", "Air chemistry: Blue",
                      "Water chemistry: Blue", "Radiation: Red", "Suitability: Blue"]
    assert len(_parse(symbols.symbol_row(body, SPRITE_URL, with_overall=False))) == 6


def test_the_overall_symbol_takes_the_colour_of_the_gated_score_not_the_worst_domain():
    from planetgen.physics import habitability
    for score in (0.95, 0.8, 0.5, 0.1):
        body = _body(suitability=3)  # Suitability Red never makes the overall worse than Yellow
        body["phi4"] = score
        assert symbols.overall_tier(body) == TIER_NAMES.index(habitability.tier(score))


def test_a_body_that_was_not_scored_has_no_row():
    gas_giant = {"phi4": None, **{f"tier_{domain}": None for domain in DOMAINS}}
    assert symbols.overall_tier(gas_giant) is None
    assert symbols.symbol_row(gas_giant, SPRITE_URL) == ""
    assert symbols.symbol_row({}, SPRITE_URL) == ""


def test_the_template_globals_point_into_the_sprite_and_are_not_escaped_again():
    app = create_app()
    with app.test_request_context("/"):
        single = habitability_symbol("water", 1, "pH 7")
        row = habitability_symbol_row(_body())
        rendered = render_template_string("{{ habitability_symbol('air', 3) }}")
    assert re.search(r'href="/static/icons\.svg\?v=[^"]+#hab-water"', str(single))
    assert 'aria-label="Water chemistry: Green, pH 7"' in str(single)
    assert str(row).count("<svg") == 7
    assert 'aria-label="Air chemistry: Red"' in rendered and "&lt;" not in rendered


def _luminance(colour):
    channels = [int(colour[index:index + 2], 16) / 255 for index in (1, 3, 5)]
    linear = [value / 12.92 if value <= 0.03928 else ((value + 0.055) / 1.055) ** 2.4 for value in channels]
    return 0.2126 * linear[0] + 0.7152 * linear[1] + 0.0722 * linear[2]


def _contrast(a, b):
    high, low = sorted((_luminance(a), _luminance(b)), reverse=True)
    return (high + 0.05) / (low + 0.05)


def _theme_blocks():
    """`{"light": {...}, "dark": {...}}` of the custom properties each theme defines in `style.css`."""
    css = _css()
    light = re.search(r"^:root \{(.*?)^\}", css, re.S | re.M).group(1)
    dark = re.search(r':root\[data-theme="dark"\] \{(.*?)^\}', css, re.S | re.M).group(1)
    auto = re.search(r':root:not\(\[data-theme="light"\]\) \{(.*?)^  \}', css, re.S | re.M).group(1)

    def parse(block):
        return dict(re.findall(r"--([a-z-]+):\s*(#[0-9a-fA-F]{6})\b", block))
    return {"light": parse(light), "dark": parse(dark), "os-dark": parse(auto)}


@pytest.mark.parametrize("theme", ["light", "dark", "os-dark"])
@pytest.mark.parametrize("tier", TIER_NAMES)
def test_every_rating_colour_reads_against_every_surface_in_both_themes(theme, tier):
    colours = _theme_blocks()[theme]
    ink = colours[f"tier-{tier.lower()}"]
    for surface in ("bg", "bg-elevated"):
        # WCAG 1.4.11: a graphic needs 3:1 against what is behind it.
        assert _contrast(ink, colours[surface]) >= 3.0, (theme, tier, surface)


def test_the_tier_classes_the_symbols_use_are_drawn_by_the_stylesheet():
    css = _css()
    for tier in TIER_NAMES:
        assert re.search(rf"^\.tier-{tier.lower()} \{{ color: var\(--tier-{tier.lower()}\); \}}", css, re.M)
    assert ".hab-symbol {" in css
