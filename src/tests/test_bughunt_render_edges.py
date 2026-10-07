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

Fixed (MAP.57; was a strict xfail, TEST.4, and before that a Tier 2 report):
`render_system_map_panel` given a `NaN`/`Inf` number (never
produced by normal generation, so this is defense in depth at the render
boundary) used to embed the literal string "nan"/"inf" into an SVG
numeric attribute. It now copies every row with such numbers set to
`None` first: a body with a lost position is drawn at its orbit distance
with a note, or left out when it has no distance either.
A `None` or negative radius renders clean.
"""

import math
import re

import pytest

from planetgen.web.lib import mdconvert  # noqa: E402
from planetgen.web.maps import systemmap as sm  # noqa: E402

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


@pytest.mark.parametrize("x_km", [float("nan"), float("inf"), float("-inf")], ids=["NaN", "inf", "-inf"])
def test_system_map_non_finite_position_renders_finite_numbers(x_km):
    """TEST.4, MAP.57: a body whose stored position is NaN or infinite
    never puts "nan"/"inf" into the SVG."""
    assert _non_finite_attrs(_planet(1, 1, x_km, 0.0, distance_km=1e8)) == []


def test_system_map_lost_position_is_drawn_at_its_orbit_distance_with_a_note():
    """MAP.57: drawn due east at its orbit distance, with a note in its
    info panel."""
    html = sm.render_system_map_panel({"name": "Test", "binary_configuration": None}, [_star(1)],
                                      [_planet(1, 1, float("nan"), 0.0, distance_km=1e8)], [])
    marker = re.search(r'<g class="sysmap-body sysmap-planet"[^>]*>', html).group(0)
    assert 'data-note="Position not recorded' in marker
    assert 'data-xkm="100000000.0"' in marker and 'data-ykm="0.0"' in marker


def test_system_map_leaves_out_a_body_with_no_position_or_distance():
    """MAP.57: nothing left to place it by, so it isn't drawn."""
    planet = _planet(1, 1, float("nan"), 0.0, distance_km=float("nan"))
    planet["name"] = "Lost"
    html = sm.render_system_map_panel({"name": "Test", "binary_configuration": None}, [_star(1)], [planet], [])
    assert "Lost" not in html
    assert not _NON_FINITE_ATTR.findall(html)


def test_system_map_non_finite_numbers_anywhere_render_finite():
    """MAP.57: NaN or infinity in a star, moon, belt, facility or the
    binary offset never reaches the SVG either."""
    nan, inf = float("nan"), float("inf")
    moon = _planet(5, None, inf, nan, radius_km=nan, gravity_g=inf, period_years=nan)
    moon["name"] = "Moonlet"
    planet = _planet(2, 1, 1.5e8, 0.0, radius_km=inf, mass_kg=nan, moons=[moon])
    stars = [_star(1, radius_km=nan, temperature_k=nan), _star(3, mass_kg=inf, luminosity_w=nan)]
    belts = [{"id": 9, "star_id": None, "distance_km": nan, "lower_limit_km": inf, "upper_limit_km": nan,
              "density": "sparse", "composition_summary": "rock"}]
    facilities = [{"id": 4, "name": "Dock", "kind": "station", "host_type": "planet", "host_id": 2,
                   "placement": "orbital", "orbit_distance_km": nan, "orbit_period_years": inf,
                   "orbital_speed_kms": nan, "orbit_phase_deg": nan}]
    system = {"name": "Test", "binary_configuration": "close", "binary_mutual_position_x_km": nan,
              "binary_mutual_position_y_km": inf}
    html = sm.render_system_map_panel(system, stars, [planet], belts, facilities)
    assert not _NON_FINITE_ATTR.findall(html)
    assert not re.search(r'="[^"]*\b(?:nan|inf)\b', html, re.IGNORECASE)
