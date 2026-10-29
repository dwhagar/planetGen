"""
Tests for `html/static/generatebuttons.js`, the admin Generate buttons
shared by the Sector Map and the Galaxy Map, run under node (skipped
where node isn't installed): its light-year to parsec conversion must
match the Python one, and its "up to about N sectors" estimate must bound
a real count of the sectors a neighborhood reaches.
"""

import json
import math
import os
import shutil
import subprocess

import pytest

from planetgen import tuning
from planetgen.generation.limits import MAX_GENERATE_RADIUS_LY
from planetgen.physics.units import ly_to_pc, pc_to_ly

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
    assert out["fallback"] == round(pc_to_ly(tuning.DEFAULT_GENERATE_RADIUS_PC))
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


def test_sector_buttons_post_the_generate_pages_address_and_shell_modes():
    out = _run(_FAKE_DOM + "console.log(JSON.stringify(forms(G.generateButtons(TARGET, 40, -1, 7, 13.046))));")
    base = {"csrf_token": "tok", "action": "galaxy"}
    address = {"address_shell": "40", "address_band": "-1", "address_wedge": "7"}
    assert [(f["label"], f["fields"]) for f in out] == [
        ("Generate this sector", {**base, "mode": "address", **address}),
        ("Generate neighborhood", {**base, "mode": "address", **address, "address_radius_pc": "11.957"}),
        ("Generate the entire shell (not recommended)",
         {**base, "mode": "shell", "shell_index": "40", "whole_shell": "1"}),
    ]


def test_the_block_button_posts_the_generate_pages_block_mode():
    """UX.127: a block of the ladder posts `mode=block` with its key, and asks first when it is large or its size unknown."""
    out = _run(_FAKE_DOM + """
const box = G.blockGenerateButtons(TARGET, "9.40.2.5", 600);
const big = G.blockGenerateButtons(TARGET, "243.0.0.0", 90000);
console.log(JSON.stringify({forms: forms(box), hint: box.children[1].textContent, bigConfirms: !!big.children[0].listeners.submit,
  smallConfirms: !!box.children[0].listeners.submit,
  unknownConfirms: !!G.blockGenerateButtons(TARGET, "9.40.2.5").children[0].listeners.submit}));""")
    assert out["forms"] == [{"action": "/admin/generate", "label": "Generate this block", "fields": {
        "csrf_token": "tok", "action": "galaxy", "mode": "block", "block_key": "9.40.2.5", "whole_block": "1"}}]
    assert out["hint"] == "Starts a background job on the Generate page."
    assert out["bigConfirms"] and out["unknownConfirms"] and not out["smallConfirms"]


def test_the_block_button_is_disabled_with_its_reason_while_a_job_runs():
    out = _run(_FAKE_DOM + """
const target = {...TARGET, busy: "Generate sectors", jobUrl: "/admin/generate#current-job"};
const box = G.blockGenerateButtons(target, "9.40.2.5", 600);
const button = box.children[0].children.find((c) => c.tag === "button");
const hint = box.children[1];
console.log(JSON.stringify({disabled: button.disabled, title: button.title, hint: hint.textContent,
  link: hint.children.map((a) => [a.href, a.textContent])}));""")
    assert out["disabled"] is True
    assert out["title"] == "Another job is running"
    assert out["hint"].startswith("Not now: Generate sectors is still running.")
    assert out["link"] == [["/admin/generate#current-job", "See its progress"]]
