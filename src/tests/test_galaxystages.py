"""
Tests for `html/static/galaxystages.js`, the Galaxy Map drill-down's
stage rules, run under node (skipped where node isn't installed): the
galaxy's outline and each block's sector total must match a sector-by-
sector count, the galaxy must list every block that can hold a sector, a
window must be the blocks around its centre, stage URLs must round-trip,
malformed ones must open the galaxy with a reason, a sector's designation
must lead to the window that shows it, the breadcrumb must name every
step, and the flight must start and end where it says (MAP.150, MAP.151,
MAP.205, MAP.206: the shell, band and wedge grid).
"""

import json
import math
import os
import shutil
import subprocess

import pytest

from planetgen.galaxy.density import build_galaxy_shape, shape_with_terms

NODE = shutil.which("node")
STATIC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "html", "static")
STAGES = os.path.join(STATIC, "galaxystages.js")
PRISMS = os.path.join(STATIC, "galaxyprisms.js")

pytestmark = pytest.mark.skipif(NODE is None, reason="node is not installed")

SHAPE = build_galaxy_shape(2800.0, 300.0, 200.0, 1.0, 2, math.radians(15.0), 0.4)
EDGE_PC = 4.0
GALAXY_RADIUS_PC = 15000.0


def _shape():
    """The test shape with the galaxy's sector threshold, as the page
    embeds it (lib/galaxymap3d._density_shape)."""
    from planetgen.galaxy.skeleton import expected_system_count_at_density_1
    from planetgen.physics.units import pc_to_ly

    return {**shape_with_terms(SHAPE), "sector_min_density": 1.0 / expected_system_count_at_density_1(pc_to_ly(EDGE_PC))}


def _run(script):
    """Runs `script` as an ES module with the stages module as `S`, the
    prisms module as `P`, the shape as `shape` and the outline as
    `outline`; returns its JSON output."""
    source = (
        f"import * as S from {json.dumps('file://' + os.path.abspath(STAGES))};\n"
        f"import * as P from {json.dumps('file://' + os.path.abspath(PRISMS))};\n"
        f"const shape = {json.dumps(_shape())};\n"
        f"const edge = {EDGE_PC};\n"
        f"const outline = S.galaxyOutline(edge, shape, {GALAXY_RADIUS_PC});\n"
        f"{script}\n"
    )
    result = subprocess.run([NODE, "--input-type=module"], input=source, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr[-3000:]
    return json.loads(result.stdout)


def test_outline_matches_every_band():
    """The outline is a safe superset: no band the density allows lies above
    a shell's extent."""
    out = _run("""
const bad = [];
for (let shell = 0; shell <= outline.outerShell + 2; shell += 37) {
  const high = P.outlineHigh(outline, shell);
  const range = P.bandRange(Math.min(shell, 4000));
  for (let band = 0; band <= Math.min(Math.max(high, 0) + 3, range[1]); band++) {
    const allowed = P.bandAllowed(shell, band, edge, shape);
    if (allowed && band > high) bad.push([shell, band, high, allowed]);
  }
}
console.log(JSON.stringify({bad, outer: outline.outerShell}));
""")
    assert out["bad"] == []
    assert out["outer"] > 100


def test_block_totals_match_a_sector_by_sector_count():
    """drillBlockTotal agrees with counting each sector in the outline that
    joins the block, for level-3 and level-9 blocks across the galaxy."""
    out = _run("""
function brute(block) {
  let n = 0;
  for (let shell = block.shell * block.m; shell < block.shell * block.m + block.m; shell++) {
    const range = P.bandRange(shell);
    for (let band = range[0]; band <= range[1]; band++) {
      for (let wedge = 0; wedge < P.wedgeCount(shell, band); wedge++) {
        if (!P.sectorInOutline(shell, band, wedge, outline)) continue;
        if (P.sameBlock(P.drillChainOf(shell, band, wedge).find(c => c.m === block.m), block)) n++;
      }
    }
  }
  return n;
}
const checks = [];
for (const [shell, band, wedge] of [[5, 0, 3], [30, 2, 10], [40, -10, 20], [57, 30, 90], [59, -2, 5]]) {
  const chain = P.drillChainOf(shell, band, wedge);
  for (const block of [chain.find(c => c.m === 3), chain.find(c => c.m === 9)]) {
    checks.push([block.m, S.drillBlockTotal(block, outline), S.blockTotal(block, outline), brute(block), S.drillBlockHolds(block, outline)]);
  }
}
console.log(JSON.stringify(checks));
""")
    for m, total, cached, brute, holds in out:
        assert total == cached == brute, (m, total, cached, brute)
        assert holds == (brute > 0)
    assert any(total > 0 for _m, total, *_ in out)


def test_the_galaxy_lists_every_block_that_can_hold_a_sector():
    """The galaxy is the level-243 blocks that can hold a sector, each with
    its bounds, one choice per block; a click on one opens the window of the
    level-81 blocks around its middle."""
    out = _run("""
const r = S.resolveStage(S.galaxyStage(), outline, edge);
const brute = [];
for (let shell = 0; shell <= Math.floor(outline.outerShell / 243); shell++) {
  const rows = P.drillBandCount(243, shell);
  for (let band = -((rows - 1) / 2); band <= (rows - 1) / 2; band++) {
    for (let wedge = 0; wedge < P.drillWedgeCount(243, shell, band); wedge++) {
      if (S.drillBlockTotal({m: 243, shell, band, wedge}, outline) > 0) brute.push([shell, band, wedge].join("."));
    }
  }
}
const listed = r.view.blocks.map(b => [b.shell, b.band, b.wedge].join("."));
console.log(JSON.stringify({
  problem: r.problem, kind: r.kind, win: r.view.win, options: r.options.length, blocks: r.view.blocks.length,
  same: JSON.stringify(listed.slice().sort()) === JSON.stringify(brute.slice().sort()),
  bounded: r.view.blocks.every(b => b.bounds && S.blockTotal(b, outline) > 0 && b.m === 243),
  picks: r.options.every(o => o.pick.kind === "cell" && o.blocks.length === 1),
  into: S.windowInto(S.galaxyStage(), r.view.blocks[0], null, edge),
  sameTwice: S.resolveStage(S.galaxyStage(), outline, edge).view.blocks === r.view.blocks,
}));
""")
    assert out["problem"] is None and out["kind"] == "cell" and out["win"] is None
    assert out["options"] == out["blocks"] > 50
    assert out["same"] and out["bounded"] and out["picks"] and out["sameTwice"]
    assert out["into"]["win"]["centre"]["m"] == 81 and out["into"]["win"]["g"] == 9


def test_camera_urls_round_trip_and_bad_ones_are_refused():
    """MAP.228: the address holds the camera, to ten digits, and a bad one is refused with a reason."""
    out = _run("""
const cameras = [
  {target: [1234.5678901234, -98765.4321, 0.0001234567], dist: 0.0000123456789, quat: [0.5, -0.5, 0.5, 0.5]},
  {target: [0, 0, 0], dist: 262144, quat: [0, 0, 0, 1]},
];
const queries = cameras.map(c => "?" + S.cameraQuery(c));
const parsed = queries.map(q => S.parseStageQuery(q));
const bad = ["?cam=", "?cam=1,2,3", "?cam=1,2,3,4,5,6,7", "?cam=1,2,3,0,0,0,0,1", "?cam=1,2,3,-4,0,0,0,1", "?cam=1,2,3,4,0,0,0,0",
  "?cam=a,2,3,4,0,0,0,1", "?cam=1,2,3,4,0,0,0,Infinity", "?cam=1,2,,4,0,0,0,1"].map(q => S.parseStageQuery(q));
console.log(JSON.stringify({queries, cameras: parsed.map(p => p.camera), problems: parsed.map(p => p.problem),
  bad: bad.map(p => [p.camera, p.problem]), unit: S.parseCamera("1,2,3,4,0,0,0,2").quat}));
""")
    assert out["queries"][1] == "?cam=0,0,0,262144,0,0,0,1"
    assert out["problems"] == [None, None]
    first = out["cameras"][0]
    assert first["target"] == [1234.567890, -98765.4321, 0.0001234567] and first["dist"] == 0.0000123456789
    assert first["quat"] == [0.5, -0.5, 0.5, 0.5]
    for camera, problem in out["bad"]:
        assert camera is None and problem and "no camera" in problem
    assert out["unit"] == [0, 0, 0, 1]


@pytest.mark.parametrize("query", ["?at=27.14.30.-1", "?p=a1.90", "?slab=2", "?at=3.1.1.0&p=L0,s9.4", "?ring=3&layer=1",
                                   "?win=27.40.-1.11", "?win=27.40.-1.11&g=5"])
def test_the_links_from_before_the_camera_address_open_the_galaxy(query):
    """MAP.150 retired the arc, slab and segment picks and MAP.228 the window tokens (?win=, &g=); their links
    have no backward-compatible reading and open the galaxy."""
    out = _run(f"console.log(JSON.stringify(S.parseStageQuery({json.dumps(query)})));")
    assert out["camera"] is None
    assert out["sector"] is None


@pytest.mark.parametrize("query", ["?sector=ZZZ", "?sector=", "?sector=123456789012"])
def test_a_malformed_sector_link_opens_the_galaxy_with_a_reason(query):
    out = _run(f"console.log(JSON.stringify(S.parseStageQuery({json.dumps(query)})));")
    assert out["camera"] is None and out["sector"] is None
    assert out["problem"]


def test_a_window_beyond_the_outline_holds_its_cells_and_one_past_the_reach_is_refused():
    """MAP.243: the halo and the satellite galaxies lie beyond the outline; a window there holds every cell the grid
    names. Past the grid's last shell nothing can be named, and the window is refused."""
    out = _run("""
const beyond = S.windowStage({m: 243, shell: 20, band: 0, wedge: 0});
const high = S.windowStage({m: 3, shell: 1400, band: 0, wedge: 0});
const past = S.windowStage({m: 243, shell: 400, band: 0, wedge: 0});
const r = [beyond, high].map(stage => S.resolveStage(stage, outline, edge));
console.log(JSON.stringify({
  problems: r.map(x => x.problem), blocks: r.map(x => x.view.blocks.length),
  past: S.resolveStage(past, outline, edge).problem,
  outside: P.sectorInOutline(4860, 0, 0, outline),
}));
""")
    assert out["problems"] == [None, None]
    assert all(count > 0 for count in out["blocks"])
    assert out["past"] and "Nothing of the galaxy" in out["past"]
    assert out["outside"] is False  # the first window lies beyond the outline


def test_a_sector_beyond_the_outline_has_a_stage_and_one_past_the_reach_has_none():
    out = _run("""
const code = S.sectorDesignation(4000, 7, 3);
const parsed = S.parseStageQuery("?sector=" + code);
const s = parsed.sector;
console.log(JSON.stringify({sector: s, stage: s ? S.sectorStage(s.shell, s.band, s.wedge) : "none",
  past: S.sectorStage(90000, 0, 0), inOutline: P.sectorInOutline(4000, 7, 3, outline)}));
""")
    assert out["sector"] == {"shell": 4000, "band": 7, "wedge": 3}
    assert out["inOutline"] is False
    assert out["stage"] == {"win": {"centre": {"m": 1, "shell": 4000, "band": 7, "wedge": 3}, "g": 9}}
    assert out["past"] is None


def test_sector_links_open_the_window_that_shows_the_sector():
    """MAP.26: a sector opens among its neighbours, the window of sectors
    centred on it."""
    out = _run("""
const shell = 700, band = -20, wedge = Math.floor(P.wedgeCount(shell, band) / 3);
const code = S.sectorDesignation(shell, band, wedge);
const parsed = S.parseStageQuery("?sector=" + code);
const stage = S.sectorStage(shell, band, wedge);
const r = S.resolveStage(stage, outline, edge);
console.log(JSON.stringify({shell, band, wedge, code, back: S.parseSectorDesignation(code), parsed, stage,
  sectors: r.view.blocks.every(b => b.m === 1),
  holds: r.view.blocks.some(b => b.shell === shell && b.wedge === wedge && b.band === band),
  crumbs: S.crumbs(stage, outline, edge).map(c => c.label)}));
""")
    sector = {"shell": out["shell"], "band": out["band"], "wedge": out["wedge"]}
    assert out["back"] == sector
    assert out["parsed"]["sector"] == sector
    assert out["parsed"]["camera"] is None
    assert out["stage"] == {"win": {"centre": {"m": 1, **sector}, "g": 9}}
    assert out["sectors"] and out["holds"]
    crumbs = out["crumbs"]
    assert crumbs[0] == "Galaxy"
    assert [c.split(" ")[0] for c in crumbs[1:]] == ["Block"] * 5 + ["Sector"]
    assert crumbs[-1] == f"Sector 700·-20·{out['wedge']:,} (1)"


def test_parent_stages_walk_back_to_the_galaxy():
    out = _run("""
let stage = S.sectorStage(600, 4, 20);
const crumbs = S.crumbs(stage, outline, edge).length;
let steps = 0;
let repeated = false;
for (let parent = S.parentStage(stage, outline, edge); parent; parent = S.parentStage(stage, outline, edge)) {
  repeated = repeated || S.sameStage(parent, stage);
  stage = parent;
  steps++;
  if (steps > 50) break;
}
console.log(JSON.stringify({steps, crumbs, repeated, last: stage}));
""")
    assert not out["repeated"]
    assert out["steps"] == out["crumbs"] - 1 == 6
    assert out["last"] == {"win": None}


def test_a_window_is_the_blocks_around_its_centre():
    out = _run("""
const stage = S.windowStage({m: 27, shell: 10, band: 0, wedge: 11});
const r = S.resolveStage(stage, outline, edge);
const middle = P.cellPositionPc(stage.win.centre, edge);
const reach = 9 * 27 * edge / 2;
console.log(JSON.stringify({
  problem: r.problem, kind: r.kind, blocks: r.view.blocks.length, options: r.options.length,
  bounded: r.view.blocks.every(b => b.bounds && S.blockTotal(b, outline) > 0 && b.m === 27),
  holdsCentre: r.view.blocks.some(b => S.sameBlock(b, stage.win.centre)),
  pick: r.options[0].pick.kind,
  near: r.view.blocks.every(b => Math.hypot(...P.cellPositionPc(b, edge).map((v, k) => v - middle[k])) <= reach + 1e-6),
}));
""")
    assert out["problem"] is None and out["kind"] == "cell" and out["pick"] == "cell"
    assert out["blocks"] == out["options"] and 100 < out["blocks"] <= 729
    assert out["bounded"] and out["holdsCentre"] and out["near"]


def test_a_click_opens_the_window_one_size_down_and_the_way_up_is_a_step_back():
    out = _run("""
const stage = S.windowStage({m: 243, shell: 3, band: 0, wedge: 4});
const r = S.resolveStage(stage, outline, edge);
const rows = [];
for (const block of r.view.blocks.slice(0, 40)) {
  const into = S.windowInto(stage, block, null, edge);
  const back = S.parentStage(into, outline, edge);
  const parent = P.drillParent(into.win.centre);
  rows.push([into.win.centre.m, P.formatDrillKey(parent) === P.formatDrillKey(block), back.win.centre.m, back.win.g, into.win.g]);
}
const sectorWindow = S.windowStage({m: 1, shell: 50, band: 3, wedge: 3});
const sectorBlock = S.resolveStage(sectorWindow, outline, edge).view.blocks[0];
const fromGalaxy = S.windowInto(S.galaxyStage(), S.resolveStage(S.galaxyStage(), outline, edge).view.blocks[5], null, edge);
const three = S.windowChildCentre({m: 3, shell: 20, band: 0, wedge: 6}, null, edge);
console.log(JSON.stringify({rows, sector: S.windowInto(sectorWindow, sectorBlock, null, edge), top: S.parentStage(S.galaxyStage(), outline, edge),
  topWindow: S.parentStage(stage, outline, edge), fromGalaxy, three,
  threeHeld: S.sameBlock(P.drillParent(three), {m: 3, shell: 20, band: 0, wedge: 6})}));
""")
    assert out["rows"] and all(m == 81 and inside and back == 243 and g == 9 and into_g == 9 for m, inside, back, g, into_g in out["rows"])
    assert out["sector"] is None
    assert out["top"] is None
    assert out["topWindow"] == {"win": None}
    assert out["fromGalaxy"]["win"]["centre"]["m"] == 81
    # A level-3 block opens a window of the sectors round the one at its middle.
    assert out["three"]["m"] == 1 and out["threeHeld"]


def test_a_click_opens_the_window_round_the_cell_under_the_click():
    """The window centres on the child cell under the click, kept inside the
    block clicked; a click off the block's sides goes to the middle."""
    out = _run("""
const block = S.resolveStage(S.galaxyStage(), outline, edge).view.blocks.find(b => b.shell === 3 && b.band === 0);
const b = block.bounds;
const point = (rf, cf, tf) => {
  const r = b.r0 + (b.r1 - b.r0) * rf, c = b.c0 + (b.c1 - b.c0) * cf, t = b.t0 + (b.t1 - b.t0) * tf;
  return {x: r * Math.sin(c) * Math.cos(t), y: r * Math.sin(c) * Math.sin(t), z: r * Math.cos(c)};
};
const rows = [];
for (const [rf, cf, tf] of [[0.05, 0.05, 0.05], [0.95, 0.95, 0.95], [0.5, 0.5, 0.5], [0.05, 0.95, 0.5], [0.95, 0.05, 0.1]]) {
  const p = point(rf, cf, tf);
  const centre = S.windowChildCentre(S.stripBlock(block), p, edge);
  const box = P.drillBlockBounds(centre, edge);
  rows.push({parent: P.formatDrillKey(P.drillParent(centre)) === P.formatDrillKey(block), holds: P.boundsContain(box, p.x, p.y, p.z)});
}
const far = S.windowChildCentre(block, {x: -1e6, y: 1e6, z: 5e5}, edge);
const middle = S.windowChildCentre(block, null, edge);
console.log(JSON.stringify({rows, farParent: P.formatDrillKey(P.drillParent(far)) === P.formatDrillKey(block),
  middleParent: P.formatDrillKey(P.drillParent(middle)) === P.formatDrillKey(block)}));
""")
    assert all(row["parent"] and row["holds"] for row in out["rows"]), out
    assert out["farParent"] and out["middleParent"]


def test_a_click_opens_the_window_round_the_middle_of_what_it_went_through():
    """A ray that meets a block's top face opens the window round the middle of its stretch inside the block,
    not round the face it met."""
    out = _run("""
const b = P.drillBlockBounds({m: 243, shell: 0, band: 0, wedge: 0}, edge);
// Straight down through the core block from above: enters at the top, leaves at the bottom.
const top = {x: 100, y: 50, z: Math.sqrt(b.r1 * b.r1 - 100 * 100 - 50 * 50) - 1e-6};
const mid = S.chordMiddle(b, top, {x: 100, y: 50, z: 1e6});
const miss = S.chordMiddle(b, {x: 1e5, y: 0, z: 0}, {x: 1e5, y: 0, z: 1e6});
console.log(JSON.stringify({mid, top, miss, none: S.chordMiddle(b, null, null) === null}));
""")
    assert out["mid"]["x"] == pytest.approx(100, abs=1e-6) and out["mid"]["y"] == pytest.approx(50, abs=1e-6)
    assert abs(out["mid"]["z"]) < 5.0 and out["top"]["z"] > 400
    assert out["miss"] == {"x": 1e5, "y": 0, "z": 0} and out["none"]


def test_the_course_stage_is_the_smallest_window_holding_every_stop():
    """courseStage (section 9.4) climbs only as far as it must: the same
    sector stays at the sector level, neighbours share a window that holds
    both, and stops across the galaxy fall back to the galaxy."""
    out = _run("""
const wedges = P.wedgeCount(700, -20);
const home = {shell: 700, band: -20, wedge: Math.floor(wedges / 3)};
const near = {shell: 701, band: -20, wedge: Math.floor(P.wedgeCount(701, -20) / 3)};
const mid = {shell: 714, band: -24, wedge: Math.floor(P.wedgeCount(714, -24) / 3)};
const core = {shell: 5, band: 0, wedge: 1};
const across = {shell: 700, band: -20, wedge: (Math.floor(wedges / 3) + Math.floor(wedges / 2)) % wedges};
const high = {shell: 5000, band: 0, wedge: 1};
const holds = (stage, s) => {
  const chain = P.drillChainOf(s.shell, s.band, s.wedge);
  const cell = chain.find(c => c.m === stage.win.centre.m);
  return S.resolveStage(stage, outline, edge).view.blocks.some(b => S.sameBlock(b, cell));
};
const nearStage = S.courseStage([home, near], outline, edge);
const midStage = S.courseStage([home, mid], outline, edge);
console.log(JSON.stringify({
  same: S.sameStage(S.courseStage([home, home], outline, edge), S.sectorStage(home.shell, home.band, home.wedge)),
  nearHolds: holds(nearStage, home) && holds(nearStage, near),
  nearLevel: nearStage.win.centre.m,
  midHolds: holds(midStage, home) && holds(midStage, mid),
  midLevel: midStage.win.centre.m,
  far: S.courseStage([home, across], outline, edge),
  apart: S.courseStage([core, home], outline, edge),
  outside: S.courseStage([core, high], outline, edge),
  none: S.courseStage([], outline, edge),
}));
""")
    assert out["same"]
    assert out["nearHolds"] and out["nearLevel"] == 1
    assert out["midHolds"] and out["midLevel"] > 1
    assert out["far"]["win"] is None or out["far"]["win"]["centre"]["m"] >= 3
    # A core sector and one 2,800 pc out need a window of the larger blocks.
    assert out["apart"]["win"] is None or out["apart"]["win"]["centre"]["m"] >= 81
    assert out["outside"] == {"win": None}
    assert out["none"] == {"win": None}


def test_a_window_stage_has_a_breadcrumb_down_from_the_galaxy():
    out = _run("""
const centre = P.drillChainOf(600, 4, 20).find(c => c.m === 3);
const stage = S.windowStage(centre);
const crumbs = S.crumbs(stage, outline, edge);
let steps = 0, current = stage, seen = [];
for (let parent = S.parentStage(current, outline, edge); parent; parent = S.parentStage(current, outline, edge)) {
  seen.push(parent.win ? parent.win.centre.m : "galaxy");
  current = parent;
  if (++steps > 10) break;
}
console.log(JSON.stringify({labels: crumbs.map(c => c.label), last: crumbs.map(c => c.last), seen, current}));
""")
    assert out["labels"][0] == "Galaxy" and len(out["labels"]) == 6
    assert [label.split(" ")[0] for label in out["labels"][1:]] == ["Block"] * 5
    assert out["labels"][-1].endswith("(3)") and out["last"] == [False] * 5 + [True]
    assert out["seen"] == [9, 27, 81, 243, "galaxy"]
    assert out["current"] == {"win": None}


def test_a_cells_edges_lie_on_its_faces():
    """MAP.52: the outline of a cell is arcs and lines on its own faces: every end point has its radius at a bound,
    or its colatitude, or its longitude, and lies inside the bounds; a cap has no southern or meridian edges to
    draw where it has none."""
    out = _run("""
const cells = [P.sectorBounds(30, 2, 7, edge), P.sectorBounds(30, 0, 0, edge), P.sectorBounds(30, P.bandRange(30)[1], 0, edge),
  P.drillBlockBounds({m: 27, shell: 3, band: 1, wedge: 2}, edge)];
console.log(JSON.stringify(cells.map(b => {
  const pts = S.cellEdgePoints(b);
  const bad = [];
  for (let k = 0; k + 2 < pts.length; k += 3) {
    const x = pts[k], y = pts[k + 1], z = pts[k + 2];
    const r = Math.hypot(x, y, z), c = Math.atan2(Math.hypot(x, y), z);
    let t = Math.atan2(y, x); if (t < 0) t += 2 * Math.PI;
    const near = (a, v) => Math.abs(a - v) < 1e-6 * Math.max(1, Math.abs(v));
    const onFace = near(r, b.r0) || near(r, b.r1) || near(c, b.c0) || near(c, b.c1) || near(t, b.t0) || near(t, b.t1);
    const inside = r >= b.r0 - 1e-6 && r <= b.r1 + 1e-6 && c >= b.c0 - 1e-6 && c <= b.c1 + 1e-6;
    if (!onFace || !inside) bad.push([x, y, z]);
  }
  return {points: pts.length / 3, even: (pts.length / 3) % 2 === 0, bad: bad.length};
})));
""")
    for cell in out:
        assert cell["bad"] == 0 and cell["even"] and cell["points"] > 0, cell


def test_extent_and_column_follow_the_cells_bounds():
    out = _run("""
const blocks = S.resolveStage(S.windowStage({m: 27, shell: 10, band: 0, wedge: 11}), outline, edge).view.blocks;
const e = S.extentOf(blocks);
const reach = blocks.every(b => S.boundsPoints([b]).every((v, k, all) => k % 3 !== 0
  || Math.hypot(v - e.center[0], all[k + 1] - e.center[1], all[k + 2] - e.center[2]) <= e.radius + 1e-6));
// A column: the points of a cell stand over the cell's own plane footprint, and nothing over a gap.
const cell = P.sectorBounds(300, 5, 11, edge);
const inside = [], outside = [];
for (let n = 0; n < 200; n++) {
  const f = 0.02 + 0.96 * ((n * 0.6180339887) % 1), g = 0.02 + 0.96 * ((n * 0.7548776662) % 1), h = 0.02 + 0.96 * ((n * 0.5698402909) % 1);
  const r = cell.r0 + (cell.r1 - cell.r0) * f, c = cell.c0 + (cell.c1 - cell.c0) * g, t = cell.t0 + (cell.t1 - cell.t0) * h;
  inside.push(S.columnContains(cell, r * Math.sin(c) * Math.cos(t), r * Math.sin(c) * Math.sin(t)));
}
const centre = P.cellCenter(cell);
const R = Math.hypot(centre[0], centre[1]);
outside.push(S.columnContains(cell, centre[0] * 2, centre[1] * 2), S.columnContains(cell, centre[0] / 2, centre[1] / 2),
  S.columnContains(cell, -centre[0], -centre[1]));
const cap = P.sectorBounds(300, P.bandRange(300)[1], 0, edge);
const core = P.sectorBounds(0, 0, 0, edge);
console.log(JSON.stringify({reach, radius: e.radius, inside, outside,
  cap: S.columnContains(cap, 0.5, 0.5), capFar: S.columnContains(cap, 1000, 0),
  core: S.columnContains(core, 1, 1), coreFar: S.columnContains(core, 3, 0)}));
""")
    assert out["reach"] and out["radius"] > 0
    assert all(out["inside"]) and not any(out["outside"])
    assert out["cap"] and not out["capFar"]
    assert out["core"] and not out["coreFar"]


def test_labels_and_ranges_read_in_shells_bands_and_wedges():
    out = _run("""
console.log(JSON.stringify({
  label: S.blockLabel({shell: 1705, band: -20, wedge: 3225}),
  win: S.windowLabel({centre: {m: 27, shell: 3, band: 0, wedge: 2}, g: 9}),
  win1: S.windowLabel({centre: {m: 1, shell: 30, band: -1, wedge: 7}, g: 9}),
  key: S.parseCellKey("27.3.-1.2"), keyBad: [S.parseCellKey("27.3"), S.parseCellKey("1.2.3.x"), S.parseCellKey(null)],
  bearing: S.bearingRange({t0: 3 * Math.PI / 2, t1: 2 * Math.PI}), whole: S.bearingRange({t0: 0, t1: 2 * Math.PI}),
  latitude: S.latitudeRange({c0: Math.PI / 3, c1: Math.PI / 2}),
}));
""")
    assert out["label"] == "1,705·-20·3,225"
    assert out["win"] == "Block 3·0·2 (27)" and out["win1"] == "Sector 30·-1·7 (1)"
    assert out["key"] == {"m": 27, "shell": 3, "band": -1, "wedge": 2}
    assert out["keyBad"] == [None, None, None]
    assert out["bearing"] == pytest.approx([270.0, 360.0])
    assert out["whole"] == pytest.approx([0.0, 360.0])
    assert out["latitude"] == pytest.approx([30.0, 0.0])


def test_the_address_bar_reads_an_address():
    out = _run("""
const wedge = Math.floor(P.wedgeCount(312, -3) / 3);
const texts = ["312/-3/" + wedge, "shell 312 band -3 wedge " + wedge, "  SHELL 312, BAND -3, WEDGE " + wedge + " ", "0/0/0"];
console.log(JSON.stringify({wedge, parsed: texts.map(t => S.parseAddress(t, edge))}));
""")
    wedge = out["wedge"]
    expected = [{"shell": 312, "band": -3, "wedge": wedge}] * 3 + [{"shell": 0, "band": 0, "wedge": 0}]
    assert [p["sector"] for p in out["parsed"]] == expected
    assert all("problem" not in p for p in out["parsed"])


def test_the_address_bar_reads_a_point_and_a_designation():
    out = _run("""
console.log(JSON.stringify({
  point: S.parseAddress("(100.5, -20, 3) pc", edge),
  loose: S.parseAddress("100.5 -20 3", edge),
  code: S.parseAddress(S.sectorDesignation(1705, -20, 3225), edge),
  cell: P.sectorAddressAt(100.5, -20, 3, edge),
}));
""")
    assert out["point"]["sector"] == out["cell"]
    assert out["point"]["point"] == [100.5, -20, 3]
    assert out["loose"]["sector"] == out["cell"]
    assert out["code"]["sector"] == {"shell": 1705, "band": -20, "wedge": 3225}


@pytest.mark.parametrize("text, kind", [
    ("Belcana Scaonon", "name"),
    ("Bead", "name"),          # short enough that it is not a designation
    ("5/0/9999", "problem"),   # shell 5 has no wedge 9999
    ("5/99/1", "problem"),     # and no band 99
    ("99999/0/1", "problem"),  # and there is no shell 99999
    ("   ", "problem"),
])
def test_the_address_bar_tells_names_from_bad_addresses(text, kind):
    out = _run(f"console.log(JSON.stringify(S.parseAddress({json.dumps(text)}, edge)));")
    assert kind in out and out[kind]
    assert "sector" not in out


def test_flight_starts_and_ends_on_its_views():
    out = _run("""
const path = S.flightPath([0, 0], 30000, [8000, -2000], 1200);
const a = path.at(0), b = path.at(path.S);
const still = S.flightPath([5, 5], 100, [5, 5], 10);
console.log(JSON.stringify({S: path.S, a, b, ms: S.flightMs(path.S), still: [still.at(0).w, still.at(still.S).w],
  mid: path.at(path.S / 2).w}));
""")
    assert out["a"]["center"] == pytest.approx([0, 0], abs=1e-6)
    assert out["a"]["w"] == pytest.approx(30000)
    assert out["b"]["center"] == pytest.approx([8000, -2000], abs=1e-6)
    assert out["b"]["w"] == pytest.approx(1200)
    assert 500 <= out["ms"] <= 1600
    assert out["still"] == pytest.approx([100, 10])
    assert out["S"] > 0
