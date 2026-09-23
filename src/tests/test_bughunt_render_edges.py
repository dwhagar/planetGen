# tests/test_bughunt_render_edges.py

"""
Tier 1/Tier 2 bug-hunt coverage: `mdconvert.py`/`systemmap.py` fed
adversarial or degenerate data -- the two rendering layers most exposed to
whatever the database actually holds (mdconvert renders the stored
wikitext/Markdown write-up; systemmap renders raw `planets`/`stars`/
`asteroid_belts` rows).

Tier 1: `markdown_to_html` never raises and always HTML-escapes dangerous
input (confirmed: `<script>...` comes back as `&lt;script&gt;...`, never
literal), across a battery of malformed/adversarial markdown/wikitext
(broken tables, unterminated formatting, huge input, control characters,
`None`).

Tier 2 (not fixed -- see this test's own docstring below for why):
`render_system_map_panel` given a `NaN`/`Inf` position or radius (never
produced by normal generation -- every physics function that could yield
one either already guards against it or was hardened by
`test_bughunt_physics_edges.py`/`test_bughunt_serialization_fuzz.py`, so
this is a defense-in-depth check at the render boundary, not a live path)
does not crash, but does silently embed the literal string "nan"/"inf"
into an SVG numeric attribute -- invalid SVG a browser will just fail to
draw that one element for, not a page-level crash. Reported here as a
known, tracked soft finding (the fix would touch every `:.1f`-style
coordinate format across 5 separate map-rendering files -- `starmap.py`/
`systemmap.py`/`galaxymap.py`/`navmap.py`/`phenomenonmap.py` -- a wider
change than this pass's scope) rather than silently left uncovered.
"""

import math
import os
import sys

import pytest

from tests.bughunt_support import Tier2Report, tier2

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


def test_markdown_to_html_none_input_does_not_crash():
    html = mdconvert.markdown_to_html(None)
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


@tier2
def test_system_map_degenerate_values_do_not_crash_soft_nan_report():
    """See this file's module docstring -- Tier 2, tracked-not-fixed."""
    report = Tier2Report()
    system = {"name": "Test", "binary_configuration": None}

    cases = {
        "None radius_km": _planet(1, 1, 1e8, 0.0, radius_km=None),
        "NaN position": _planet(1, 1, float("nan"), 0.0),
        "Inf position": _planet(1, 1, float("inf"), 0.0),
        "negative radius_km": _planet(1, 1, 1e8, 0.0, radius_km=-100.0),
    }
    for label, planet in cases.items():
        try:
            html = sm.render_system_map_panel(system, [_star(1)], [planet], [])
        except Exception as exc:  # noqa: BLE001
            pytest.fail(f"{label}: render_system_map_panel raised {type(exc).__name__}: {exc}")
        if "nan" in html.lower() or "inf" in html.lower():
            report.add(f"{label}: rendered SVG contains a literal 'nan'/'inf' numeric attribute")

    report.flush("systemmap-degenerate-values")
