# tests/test_bughunt_render_edges.py

"""
Tier 1 bug-hunt coverage: `mdconvert.py`/`systemmap.py` fed
adversarial or degenerate data -- the two rendering layers most exposed to
whatever the database actually holds (mdconvert renders the stored
wikitext/Markdown write-up; systemmap renders raw `planets`/`stars`/
`asteroid_belts` rows).

Tier 1: `markdown_to_html` never raises and always HTML-escapes dangerous
input (confirmed: `<script>...` comes back as `&lt;script&gt;...`, never
literal), across a battery of malformed/adversarial markdown/wikitext
(broken tables, unterminated formatting, huge input, control characters,
`None`).

Known bug (strict xfail, TEST.4; was a Tier 2 report):
`render_system_map_panel` given a `NaN`/`Inf` position (never
produced by normal generation -- every physics function that could yield
one either already guards against it or was hardened by
`test_bughunt_physics_edges.py`/`test_bughunt_serialization_fuzz.py`, so
this is a defense-in-depth check at the render boundary, not a live path)
does not crash, but does silently embed the literal string "nan"/"inf"
into an SVG numeric attribute -- invalid SVG a browser will just fail to
draw that one element for, not a page-level crash. Kept as a strict xfail
(the fix would touch every `:.1f`-style
coordinate format across 5 separate map-rendering files -- `starmap.py`/
`systemmap.py`/`galaxymap.py`/`navmap.py`/`phenomenonmap.py` -- a wider
change than this pass's scope). A `None` or negative radius renders clean.
"""

import math
import os
import re
import sys

import pytest

_SRC_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_SRC_DIR, "html", "lib"))

import mdconvert  # noqa: E402
import systemmap as sm  # noqa: E402

ADVERSARIAL_MARKDOWN = [
    "",
    "|||broken|table|||",
    "# " * 1000,
    "```unterminated code block",
    "**unterminated bold",
    "| a | b |\n|---|\n| only one cell |",
    "a" * 100_000,
    "\x00\x01\x02 control chars",
    "<script>alert(1)</script>",
    "***" * 500,
    "\n\n\n\n\n",
    "|" * 500,
]


@pytest.mark.parametrize("text", ADVERSARIAL_MARKDOWN)
def test_markdown_to_html_never_crashes_on_adversarial_input(text):
    html = mdconvert.markdown_to_html(text)
    assert isinstance(html, str)


def test_markdown_to_html_script_tag_is_escaped_not_executable():
    html = mdconvert.markdown_to_html("<script>alert(document.cookie)</script>")
    assert "<script>" not in html
    assert "&lt;script&gt;" in html


def test_markdown_to_html_event_handler_attribute_is_escaped():
    html = mdconvert.markdown_to_html('<img src=x onerror="alert(1)">')
    assert "<img" not in html
    assert "onerror=" not in html or "&quot;" in html


def _star(id_, **kw):
    d = {
        "id": id_, "mass_kg": 1.989e30, "radius_km": 696_000, "star_type": "G2V",
        "temperature_k": 5778.0, "luminosity_w": 3.828e26,
    }
    d.update(kw)
    return d


def _planet(id_, star_id, x_km, y_km, **kw):
    d = {
        "id": id_, "star_id": star_id, "name": "P", "position_x_km": x_km, "position_y_km": y_km,
        "position_z_km": 0.0, "distance_km": math.hypot(x_km, y_km) if math.isfinite(x_km) and math.isfinite(y_km) else x_km,
        "radius_km": 6371.0, "planet_class": "G", "body_type": "t", "zone": "e", "period_years": 1.0,
        "gravity_g": 1.0, "life_chemical": None, "moons": [], "atmosphere": None, "composition": None,
        "surface_temperature_k": None,
    }
    d.update(kw)
    return d


# A number-valued attribute holding NaN or infinity (`cx="nan"`,
# `data-xkm="inf"`): matched as a whole attribute value, so class names like
# `sysmap-info` don't count (the old Tier 2 check did match them, which is
# why it flagged even the None and negative radius cases).
_NON_FINITE_ATTR = re.compile(r'[\w-]+="-?(?:nan|inf)(?: km)?"', re.IGNORECASE)


def _non_finite_attrs(planet):
    html = sm.render_system_map_panel({"name": "Test", "binary_configuration": None}, [_star(1)], [planet], [])
    return _NON_FINITE_ATTR.findall(html)


@pytest.mark.parametrize("planet", [
    _planet(1, 1, 1e8, 0.0, radius_km=None),
    _planet(1, 1, 1e8, 0.0, radius_km=-100.0),
    _planet(1, 1, 1e8, 0.0, radius_km=0.0),
], ids=["None radius", "negative radius", "zero radius"])
def test_system_map_odd_radius_renders_finite_numbers(planet):
    """TEST.4: was half of a Tier 2 report; these already render clean."""
    assert _non_finite_attrs(planet) == []


@pytest.mark.xfail(strict=True, reason="systemmap writes NaN/inf positions straight into SVG attributes "
                                       "(cx=\"nan\", data-xkm=\"inf\"); generation never stores one, so "
                                       "this is defense in depth at the render boundary")
@pytest.mark.parametrize("x_km", [float("nan"), float("inf")], ids=["NaN", "inf"])
def test_system_map_non_finite_position_renders_finite_numbers(x_km):
    """TEST.4: the other half of the old Tier 2 report, now a strict xfail.
    The page still renders (no exception); only the bad body's own numbers
    are broken."""
    assert _non_finite_attrs(_planet(1, 1, x_km, 0.0)) == []
