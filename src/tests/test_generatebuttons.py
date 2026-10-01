"""
Tests for `html/static/generatebuttons.js`, the admin Generate buttons
shared by the Sector Map and the Galaxy Map, run under node (skipped
where node isn't installed): its light-year to parsec conversion must
match the Python one, and its "up to about N sectors" estimate must bound
a real count of the sector slots a neighborhood reaches.
"""

import json
import math
import os
import shutil
import subprocess

import pytest

from stellarObjects import program_constants
from stellarObjects.generationLimits import MAX_GENERATE_RADIUS_LY
from stellarObjects.utils import ly_to_pc, pc_to_ly

NODE = shutil.which("node")
MODULE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "html", "static", "generatebuttons.js")

pytestmark = pytest.mark.skipif(NODE is None, reason="node is not installed")


def _run(script):
    source = (
        f"import * as G from {json.dumps('file://' + os.path.abspath(MODULE))};\n"
        f"{script}\n"
    )
    result = subprocess.run([NODE, "--input-type=module", "-e", source], capture_output=True, text=True, check=True)
    return json.loads(result.stdout)


@pytest.mark.parametrize("ly", [13, 100, 250.5, 652])
def test_light_years_convert_the_same_way_python_does(ly):
    got = _run(f"console.log(JSON.stringify(G.lyToPc({ly})));")
    assert got == pytest.approx(ly_to_pc(ly), rel=1e-12)


def test_the_radius_bounds_match_the_generation_limits():
    out = _run("console.log(JSON.stringify({min: G.NEIGHBORHOOD_MIN_LY, max: G.NEIGHBORHOOD_MAX_LY, "
               "fallback: G.NEIGHBORHOOD_DEFAULT_LY, confirm: G.NEIGHBORHOOD_CONFIRM_SECTORS}));")
    assert out["fallback"] == round(pc_to_ly(program_constants.DEFAULT_GENERATE_RADIUS_PC))
    assert out["confirm"] == 5000
    # Never past what the Generate page accepts, and never under one sector.
    assert out["max"] <= MAX_GENERATE_RADIUS_LY
    assert out["min"] >= 13.0 - 1e-9


@pytest.mark.parametrize("radius_ly, edge_ly", [(100, 13.046), (652, 13.046), (100, 3.0), (40, 13.046)])
def test_the_sector_estimate_bounds_a_real_count(radius_ly, edge_ly):
    """The estimate is the sphere's volume over a sector's, so it is close
    to (and never far under) the cubes of that size whose centers the
    sphere holds."""
    got = _run(f"console.log(JSON.stringify(G.neighborhoodSectors({radius_ly}, {edge_ly})));")
    steps = int(radius_ly // edge_ly) + 1
    counted = 0
    for i in range(-steps, steps + 1):
        for j in range(-steps, steps + 1):
            for k in range(-steps, steps + 1):
                if math.dist((i * edge_ly, j * edge_ly, k * edge_ly), (0, 0, 0)) <= radius_ly:
                    counted += 1
    assert counted * 0.8 <= got <= counted * 1.3


def test_the_smallest_radius_still_estimates_a_few_sectors():
    """At the minimum radius, barely more than the sector itself."""
    got = _run("console.log(JSON.stringify(G.neighborhoodSectors(13, 13.046)));")
    assert 1 <= got <= 10


def test_an_unknown_sector_size_has_no_estimate():
    out = _run("console.log(JSON.stringify([G.neighborhoodSectors(100, undefined), "
               "G.neighborhoodSectors(100, 0), G.neighborhoodSectors(0, 13)]));")
    assert out == [None, None, None]


# A bare-bones `document` for the DOM-building functions: elements keep
# their properties and children, which is all the forms need.
_FAKE_DOM = """
class El {
  constructor(tag) { this.tag = tag; this.children = []; this.listeners = {}; }
  appendChild(c) { this.children.push(c); return c; }
  insertBefore(c, ref) { const i = this.children.indexOf(ref); this.children.splice(i < 0 ? 0 : i, 0, c); return c; }
  get firstChild() { return this.children[0]; }
  get lastChild() { return this.children[this.children.length - 1]; }
  addEventListener(name, fn) { this.listeners[name] = fn; }
}
globalThis.document = { createElement: (tag) => new El(tag) };
function forms(box) {
  return box.children.filter((c) => c.tag === "form").map((f) => ({
    action: f.action,
    label: f.children.find((c) => c.tag === "button").textContent,
    fields: Object.fromEntries(f.children.filter((c) => c.tag === "input").map((i) => [i.name, i.value])),
  }));
}
const TARGET = {url: "/admin/generate", csrfField: "csrf_token", csrfToken: "tok"};
"""


def test_block_buttons_post_the_generate_pages_block_mode():
    out = _run(_FAKE_DOM + "console.log(JSON.stringify(forms(G.blockGenerateButtons(TARGET, "
               "{block: '3.40.7.0', layer: -1, wholeBlock: true}))));")
    common = {"csrf_token": "tok", "action": "galaxy", "mode": "block", "block": "3.40.7.0"}
    assert out == [
        {"action": "/admin/generate", "label": "Generate this layer", "fields": {**common, "block_layer": "-1"}},
        {"action": "/admin/generate", "label": "Generate this block", "fields": {**common, "whole_block": "1"}},
    ]


def test_block_buttons_leave_out_what_the_offer_leaves_out():
    out = _run(_FAKE_DOM + "console.log(JSON.stringify(["
               "forms(G.blockGenerateButtons(TARGET, {block: '3.1.2.0', layer: null, wholeBlock: true})).map((f) => f.label),"
               "forms(G.blockGenerateButtons(TARGET, {block: '3.1.2.0', layer: 1, wholeBlock: false})).map((f) => f.label)]));")
    assert out == [["Generate this block"], ["Generate this layer"]]
