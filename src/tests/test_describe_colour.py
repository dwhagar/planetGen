"""GEN.420: the description's "Atmosphere Colour" row and the swatch inline item in every output."""

import pytest

from planetgen.describe import describe, model
from planetgen.describe.html import inline_html
from planetgen.describe.quantity import display_value, inline_text
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.util import draw

SWATCH = model.swatch("#99bfff", "pale blue")


def test_a_swatch_is_written_by_every_format():
    assert inline_text([SWATCH], "text") == "pale blue (#99bfff)"
    assert inline_text(["Sky: ", SWATCH], "markdown") == "Sky: pale blue (#99bfff)"
    assert inline_text([SWATCH], "wikitext") == '<span style="background:#99bfff">&emsp;</span> pale blue'
    html = inline_html([SWATCH])
    assert 'class="swatch"' in html and 'data-colour="#99bfff"' in html and "style=" not in html and html.endswith(" pale blue")


def test_a_swatch_is_checked_and_found_by_the_walker():
    document = model.document("system", "system:1", "S", [model.section("planet", "planet:1", "P", [
        model.facts("Planet", [model.fact("atmosphere_colour", "Atmosphere Colour", SWATCH)])])])
    model.check(document)
    assert SWATCH in list(model.walk_inline(document))
    for bad in (model.swatch("blue", "blue"), model.swatch("#99bfff", ""), model.swatch("#99BFFF", "pale blue")):
        broken = model.document("system", "system:1", "S", [model.section("planet", "planet:1", "P", [
            model.facts("Planet", [model.fact("atmosphere_colour", "Atmosphere Colour", bad)])])])
        with pytest.raises(ValueError):
            model.check(broken)


def test_a_solid_body_with_air_gets_the_row_and_an_airless_one_or_a_giant_does_not():
    with_air = airless = giants = 0
    for seed in range(12):
        draw.set_run_seed(seed)
        system = StarSystem(system_config=SystemConfig())
        for section in model.sections(describe(system)):
            if section["kind"] not in ("planet", "moon"):
                continue
            rows = {row["key"]: row for block in section["blocks"] if block["type"] == "facts" for row in block["rows"]}
            if "pressure" not in rows:                      # a giant has an internal pressure instead
                assert "atmosphere_colour" not in rows
                giants += 1
            elif display_value(rows["pressure"]["value"][0]["qty"]) > 0:
                swatch = rows["atmosphere_colour"]["value"][0]["swatch"]
                assert len(swatch["hex"]) == 7 and swatch["word"]
                with_air += 1
            else:
                assert "atmosphere_colour" not in rows
                airless += 1
    assert with_air > 0 and airless > 0 and giants > 0
