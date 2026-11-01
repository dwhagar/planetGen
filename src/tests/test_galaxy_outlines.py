# tests/test_galaxy_outlines.py

"""
UX.135: the galaxy family outlines. Six vector drawings, one SVG file each under `static/galaxy/`, listed once on the Python
side (`web/lib/galaxy_outlines.py`). The files and the list are held to each other and to the family keys of the database here;
`tests/test_galaxy_outlines_browser.py` checks what a browser draws from them in both themes.
"""

import os
import re
import xml.etree.ElementTree as ET

import pytest

from planetgen.db import models
from planetgen.web import galaxy_outline
from planetgen.web.app import create_app
from planetgen.web.lib import galaxy_outlines

STATIC = os.path.join(os.path.dirname(__file__), "..", "html", "static")
SVG = "{http://www.w3.org/2000/svg}"
MAX_BYTES = 4096
IDS = tuple(outline.id for outline in galaxy_outlines.OUTLINES)


def _file(outline_id):
    return os.path.join(STATIC, "galaxy", f"{outline_id}.svg")


def _text(outline_id):
    with open(_file(outline_id), encoding="utf-8") as handle:
        return handle.read()


def test_the_outline_list_is_the_six_the_design_names_and_the_files_are_exactly_those():
    assert IDS == ("spiral", "barred_spiral", "lenticular", "elliptical", "irregular", "dwarf_spheroid")
    assert sorted(name[:-4] for name in os.listdir(os.path.join(STATIC, "galaxy"))) == sorted(IDS)
    assert len(set(IDS)) == len(IDS)


def test_every_family_key_of_the_database_has_an_outline_and_only_those_keys_are_drawn():
    constraint = next(item for item in models.galaxy_shape.constraints if item.name == "chk_galaxy_shape_family")
    keys = tuple(re.findall(r"'(\w+)'", str(constraint.sqltext).replace("_utf8mb4", "")))
    assert sorted(keys) == sorted(galaxy_outlines.FAMILY_KEYS)
    assert {outline.family for outline in galaxy_outlines.OUTLINES} == set(keys)
    for family in keys:
        assert galaxy_outlines.outline_id(family) in IDS
    assert galaxy_outlines.outline_id("spiral", bar=True) == "barred_spiral"
    assert galaxy_outlines.outline_id("elliptical", bar=True) == "elliptical"
    with pytest.raises(ValueError):
        galaxy_outlines.outline_id("quasar")


@pytest.mark.parametrize("outline_id", IDS)
def test_each_file_is_a_plain_outline_with_a_name_a_viewbox_and_no_script_or_external_reference(outline_id):
    text = _text(outline_id)
    assert len(text.encode("utf-8")) <= MAX_BYTES
    root = ET.fromstring(text)
    assert root.tag == f"{SVG}svg" and root.get("id") == galaxy_outlines.SPRITE_ID
    assert re.fullmatch(r"0 0 \d+ \d+", root.get("viewBox") or ""), root.get("viewBox")
    assert root.findtext(f"{SVG}title") == galaxy_outlines.BY_ID[outline_id].label
    tags = {element.tag.replace(SVG, "") for element in root.iter()}
    assert tags <= {"svg", "title", "g", "path", "circle", "ellipse", "rect", "line", "polyline", "polygon"}, tags
    lowered = text.replace('xmlns="http://www.w3.org/2000/svg"', "").lower()      # the one namespace is the only address
    for forbidden in ("<script", "<style", "<image", "<foreignobject", "href=", "url(", "http:", "https:", "onload", "javascript:"):
        assert forbidden not in lowered, forbidden.replace("<", "")
    assert text.count("xmlns=") == 1


@pytest.mark.parametrize("outline_id", IDS)
def test_each_outline_is_one_stroke_weight_in_the_text_colour_with_no_fill(outline_id):
    root = ET.fromstring(_text(outline_id))
    group = root.find(f"{SVG}g")
    assert group.get("fill") == "none" and group.get("stroke") == "currentColor" and group.get("stroke-width") == "2"
    for element in root.iter():
        for attribute in ("fill", "stroke", "stroke-width"):
            if element is not group and element.get(attribute) is not None:
                pytest.fail(f"{element.tag} sets its own {attribute}: one stroke weight, one paint")
    assert not re.search(r"#[0-9a-fA-F]{3,6}|rgb\(", _text(outline_id))


def test_the_outline_markup_names_the_drawing_and_points_into_its_file():
    with create_app().test_request_context("/"):
        markup = str(galaxy_outline("barred_spiral"))
        with pytest.raises(ValueError):
            galaxy_outline("no-such-galaxy")
    assert re.fullmatch(r'<svg class="galaxy-outline" role="img" aria-label="Barred spiral galaxy" data-galaxy-outline="barred_spiral" '
                        r'focusable="false"><use href="/static/galaxy/barred_spiral\.svg\?v=[^"]+#galaxy"></use></svg>', markup)


def test_the_stylesheet_paints_the_outline_in_the_text_colour():
    with open(os.path.join(STATIC, "style.css"), encoding="utf-8") as handle:
        rule = re.search(r"^\.galaxy-outline \{\n(.*?)^\}", handle.read(), re.S | re.M).group(1)
    assert "color: var(--text);" in rule and "stroke: currentColor;" in rule and "fill: none;" in rule
