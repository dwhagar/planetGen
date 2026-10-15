"""
Tests for `html/static/galaxystages.js`, the Galaxy Map drill-down's
stage rules, run under node (skipped where node isn't installed): the
galaxy's outline and each block's sector total must match a sector-by-
sector count, the galaxy must list every block that can hold a sector, a
window must be the blocks around its centre, stage URLs must round-trip,
malformed ones must open the galaxy with a reason, a sector's designation
must lead to the window that shows it, the breadcrumb must name every
step, and the flight must start and end where it says.
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
    result = subprocess.run([NODE, "--input-type=module", "-e", source], capture_output=True, text=True, check=True)
    return json.loads(result.stdout)


def test_outline_matches_every_layer():
    """Each ring's extent is the highest drawable |layer|, and every layer
    up to it is drawable (the bound only falls away from the plane)."""
    out = _run("""
const bad = [];
for (let ring = 0; ring <= outline.maxRing + 2 && ring < outline.extents.length; ring += 7) {
  const reach = outline.extents[ring];
  for (let layer = 0; layer <= Math.max(reach, 0) + 3; layer++) {
    const drawable = P.sectorDrawable(ring, layer, edge, shape);
    if (drawable !== (layer <= reach)) bad.push([ring, layer, reach, drawable]);
  }
}
console.log(JSON.stringify({bad, maxRing: outline.maxRing, maxLayer: outline.maxLayer}));
""")
    assert out["bad"] == []
    assert out["maxRing"] > 100 and out["maxLayer"] > 10


def test_block_totals_match_a_sector_by_sector_count():
    """drillBlockTotal agrees with counting each sector the skeleton
    draws, for level-3 and level-27 blocks across the galaxy."""
    out = _run("""
function brute(block) {
  const half = (block.m - 1) / 2;
  let n = 0;
  for (let ring = block.ring * block.m; ring < block.ring * block.m + block.m; ring++) {
    const slots = P.drillSlotRange(block, ring);
    for (let layer = block.slab * block.m - half; layer <= block.slab * block.m + half; layer++) {
      if (!P.sectorDrawable(ring, layer, edge, shape)) continue;
      n += Math.max(0, slots.last - slots.first + 1);
    }
  }
  return n;
}
const checks = [];
for (const [ring, layer, slot] of [[5, 0, 3], [300, 2, 1000], [700, -10, 4000], [1500, 30, 9000], [40, -60, 200]]) {
  const chain = P.drillChainOf(ring, layer, slot);
  for (const block of [chain[1], chain[2]]) {
    checks.push([block.m, S.drillBlockTotal(block, outline), brute(block)]);
  }
}
console.log(JSON.stringify(checks));
""")
    for m, total, brute in out:
        assert total == brute, (m, total, brute)
    assert any(total > 0 for _m, total, _b in out)



def test_the_galaxy_lists_every_block_that_can_hold_a_sector():
    """The galaxy is the level-243 blocks that can hold a sector, each with
    its bounds and total, one choice per block; a click on one opens the
    window of the level-27 blocks around its middle."""
    out = _run("""
const r = S.resolveStage(S.galaxyStage(), outline, edge);
const brute = [];
for (let ring = 0; ring <= Math.floor(outline.maxRing / 243); ring++) {
  for (let slab = -10; slab <= 10; slab++) {
    for (let wedge = 0; wedge < P.drillWedgeCount(243, ring); wedge++) {
      if (S.drillBlockTotal({m: 243, ring, wedge, slab}, outline) > 0) brute.push([ring, wedge, slab].join("."));
    }
  }
}
const listed = r.view.blocks.map(b => [b.ring, b.wedge, b.slab].join("."));
console.log(JSON.stringify({
  problem: r.problem, kind: r.kind, win: r.view.win, options: r.options.length, blocks: r.view.blocks.length,
  same: JSON.stringify(listed.slice().sort()) === JSON.stringify(brute.slice().sort()),
  bounded: r.view.blocks.every(b => b.bounds && b.total > 0 && b.m === 243),
  picks: r.options.every(o => o.pick.kind === "cell" && o.blocks.length === 1),
  into: S.windowInto(S.galaxyStage(), r.view.blocks[0]),
  sameTwice: S.resolveStage(S.galaxyStage(), outline, edge).view.blocks === r.view.blocks,
  query: S.stageQuery(S.galaxyStage()),
}));
""")
    assert out["problem"] is None and out["kind"] == "cell" and out["win"] is None
    assert out["options"] == out["blocks"] > 50
    assert out["same"] and out["bounded"] and out["picks"] and out["sameTwice"]
    assert out["into"]["win"]["centre"]["m"] == 27 and out["into"]["win"]["g"] == 9
    assert out["query"] == ""


def test_stage_urls_round_trip_and_bad_ones_are_refused():
    out = _run("""
const centre = {m: 27, ring: 40, wedge: 11, slab: -1};
const stages = [S.galaxyStage(), S.windowStage(centre), S.windowStage(centre, 5), S.windowStage({m: 1, ring: 5, wedge: 7, slab: 0})];
const queries = stages.map(S.stageQuery);
const parsed = queries.map(q => S.parseStageQuery(q));
const bad = ["?win=81.0.0.0", "?win=27.3.99999.0", "?win=x", "?win=27.3.5.0&g=4", "?win=27.3.5.0&g=x", "?win=27.3.5.0&g=33"]
  .map(q => S.parseStageQuery(q));
console.log(JSON.stringify({queries, same: parsed.map((p, n) => S.sameStage(p.stage, stages[n]) && p.problem === null),
  bad: bad.map(p => [p.problem, p.stage])}));
""")
    assert out["queries"] == ["", "?win=27.40.11.-1", "?win=27.40.11.-1&g=5", "?win=1.5.7.0"]
    assert all(out["same"])
    for problem, stage in out["bad"]:
        assert problem and stage == {"win": None}


@pytest.mark.parametrize("query", ["?at=27.14.30.-1", "?p=a1.90", "?slab=2", "?at=3.1.1.0&p=L0,s9.4"])
def test_the_links_from_before_the_window_stages_open_the_galaxy(query):
    """MAP.150 retired the arc, slab and segment picks; their links have no
    backward-compatible reading and open the galaxy."""
    out = _run(f"console.log(JSON.stringify(S.parseStageQuery({json.dumps(query)})));")
    assert out["stage"] == {"win": None}
    assert out["sector"] is None


@pytest.mark.parametrize("query", ["?sector=ZZZ"])
def test_a_malformed_sector_link_opens_the_galaxy_with_a_reason(query):
    out = _run(f"console.log(JSON.stringify(S.parseStageQuery({json.dumps(query)})));")
    assert out["stage"] == {"win": None}
    assert out["problem"]


def test_a_window_outside_the_galaxy_is_refused():
    out = _run("""
const far = S.windowStage({m: 243, ring: 15, wedge: 0, slab: 30});
const high = S.windowStage({m: 3, ring: 300, wedge: 0, slab: 400});
console.log(JSON.stringify([S.resolveStage(far, outline, edge).problem, S.resolveStage(high, outline, edge).problem]));
""")
    assert all(problem and "Nothing of the galaxy" in problem for problem in out)


def test_a_sector_far_outside_the_galaxy_parses_but_has_no_stage():
    out = _run("""
const parsed = S.parseStageQuery("?sector=FFFFFFFFFFFFFFFF");
const s = parsed.sector;
console.log(JSON.stringify({sector: s, stage: s ? S.sectorStage(s.ring, s.layer, s.slot, outline) : "none"}));
""")
    assert out["sector"]
    assert out["stage"] is None


def test_sector_links_open_the_window_that_shows_the_sector():
    """MAP.26: a sector opens among its neighbours, the window of sectors
    centred on it."""
    out = _run("""
const ring = 1705, layer = -20, slot = 3225;
const code = S.sectorDesignation(ring, layer, slot);
const parsed = S.parseStageQuery("?sector=" + code);
const stage = S.sectorStage(ring, layer, slot, outline);
const r = S.resolveStage(stage, outline, edge);
console.log(JSON.stringify({code, back: S.parseSectorDesignation(code), parsed, stage,
  sectors: r.view.blocks.every(b => b.m === 1),
  holds: r.view.blocks.some(b => b.ring === ring && b.wedge === slot && b.slab === layer),
  crumbs: S.crumbs(stage, outline, edge).map(c => c.label)}));
""")
    assert out["back"] == {"ring": 1705, "layer": -20, "slot": 3225}
    assert out["parsed"]["sector"] == out["back"]
    assert out["parsed"]["stage"] == {"win": None}
    assert out["stage"] == {"win": {"centre": {"m": 1, "ring": 1705, "wedge": 3225, "slab": -20}, "g": 9}}
    assert out["sectors"] and out["holds"]
    crumbs = out["crumbs"]
    assert crumbs[0] == "Galaxy"
    assert [c.split(" ")[0] for c in crumbs[1:]] == ["Block", "Block", "Block", "Sector"]
    assert crumbs[-1] == "Sector 1,705·-20·3,225 (1)"


def test_parent_stages_walk_back_to_the_galaxy():
    out = _run("""
let stage = S.sectorStage(900, 4, 2000, outline);
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
    assert out["steps"] == out["crumbs"] - 1 == 4
    assert out["last"] == {"win": None}


def test_a_window_is_the_blocks_around_its_centre():
    out = _run("""
const stage = S.windowStage({m: 27, ring: 40, wedge: 11, slab: 0});
const r = S.resolveStage(stage, outline, edge);
console.log(JSON.stringify({
  problem: r.problem, kind: r.kind, blocks: r.view.blocks.length, options: r.options.length,
  bounded: r.view.blocks.every(b => b.bounds && b.total > 0 && b.m === 27),
  holdsCentre: r.view.blocks.some(b => S.sameBlock(b, stage.win.centre)),
  pick: r.options[0].pick.kind, mid: (r.view.a0 + r.view.a1) / 2,
  midBearing: (11 + 0.5) / P.drillWedgeCount(27, 40) * 2 * Math.PI,
}));
""")
    assert out["problem"] is None and out["kind"] == "cell" and out["pick"] == "cell"
    assert out["blocks"] == out["options"] and 300 < out["blocks"] <= 729
    assert out["bounded"] and out["holdsCentre"]
    # Bearings are counted within half a turn of the centre's.
    assert out["mid"] == pytest.approx(out["midBearing"])


def test_a_click_opens_the_window_one_size_down_and_the_way_up_is_a_step_back():
    out = _run("""
const stage = S.windowStage({m: 243, ring: 3, wedge: 4, slab: 0});
const r = S.resolveStage(stage, outline, edge);
const rows = [];
for (const block of r.view.blocks.slice(0, 40)) {
  const into = S.windowInto(stage, block);
  const back = S.parentStage(into, outline, edge);
  const parent = P.drillParent(into.win.centre);
  rows.push([into.win.centre.m, P.formatDrillKey(parent) === P.formatDrillKey(block), back.win.centre.m, back.win.g, into.win.g]);
}
const sectorWindow = S.windowStage({m: 1, ring: 50, wedge: 3, slab: 0});
const sectorBlock = S.resolveStage(sectorWindow, outline, edge).view.blocks[0];
const fromGalaxy = S.windowInto(S.galaxyStage(), S.resolveStage(S.galaxyStage(), outline, edge).view.blocks[5]);
console.log(JSON.stringify({rows, sector: S.windowInto(sectorWindow, sectorBlock), top: S.parentStage(S.galaxyStage(), outline, edge),
  topWindow: S.parentStage(stage, outline, edge), fromGalaxy,
  three: S.windowChildCentre({m: 3, ring: 20, wedge: 6, slab: 2})}));
""")
    assert all(m == 27 and inside and back == 243 and g == 9 and into_g == 9 for m, inside, back, g, into_g in out["rows"])
    assert out["sector"] is None
    assert out["top"] is None
    assert out["topWindow"] == {"win": None}
    assert out["fromGalaxy"]["win"]["centre"]["m"] == 27
    three = out["three"]
    assert three["m"] == 1 and three["slab"] == 6


def test_a_click_opens_the_window_round_the_cell_under_the_click():
    """The window centres on the child cell under the click, kept inside the
    block clicked; a click off the block's sides goes to the nearest."""
    out = _run("""
const block = S.resolveStage(S.galaxyStage(), outline, edge).view.blocks.find(b => b.ring === 3 && b.slab === 0);
const b = block.bounds;
const point = (rf, tf) => { const r = b.r0 + (b.r1 - b.r0) * rf, t = b.t0 + (b.t1 - b.t0) * tf; return {x: r * Math.cos(t), y: r * Math.sin(t)}; };
const rows = [];
for (const [rf, tf] of [[0.05, 0.05], [0.95, 0.95], [0.5, 0.5], [0.05, 0.95], [0.95, 0.05]]) {
  const centre = S.windowChildCentre(S.stripBlock(block), point(rf, tf), edge);
  const box = P.drillBlockBounds(centre, edge);
  const p = point(rf, tf), r = Math.hypot(p.x, p.y);
  rows.push({parent: P.formatDrillKey(P.drillParent(centre)) === P.formatDrillKey(block), rIn: r >= box.r0 - 1e-6 && r <= box.r1 + 1e-6,
    tIn: (() => { const t = Math.atan2(p.y, p.x); const d = (((t - box.t0) % (2 * Math.PI)) + 2 * Math.PI) % (2 * Math.PI); return d <= box.t1 - box.t0 + 1e-9; })()});
}
const far = S.windowChildCentre(block, {x: -1e6, y: 1e6}, edge);
const middle = S.windowChildCentre(block);
console.log(JSON.stringify({rows, farParent: P.formatDrillKey(P.drillParent(far)) === P.formatDrillKey(block),
  middleParent: P.formatDrillKey(P.drillParent(middle)) === P.formatDrillKey(block)}));
""")
    assert all(row["parent"] and row["rIn"] and row["tIn"] for row in out["rows"]), out
    assert out["farParent"] and out["middleParent"]


def test_the_course_stage_is_the_smallest_window_holding_every_stop():
    """courseStage (section 9.4) climbs only as far as it must: the same
    sector stays at the sector level, neighbours share a window that holds
    both, and stops across the galaxy fall back to the galaxy."""
    out = _run("""
const home = {ring: 1705, layer: -20, slot: 3225};
const near = {ring: 1706, layer: -20, slot: 3226};
const mid = {ring: 1725, layer: -22, slot: 3240};
const core = {ring: 5, layer: 0, slot: 1};
const across = {ring: 1705, layer: -20, slot: (3225 + P.ringSectorCount(1705) / 2) % P.ringSectorCount(1705)};
const high = {ring: 5, layer: 2000, slot: 1};
const holds = (stage, s) => {
  const chain = P.drillChainOf(s.ring, s.layer, s.slot);
  const cell = chain.find(c => c.m === stage.win.centre.m);
  return S.resolveStage(stage, outline, edge).view.blocks.some(b => S.sameBlock(b, cell));
};
const nearStage = S.courseStage([home, near], outline, edge);
const midStage = S.courseStage([home, mid], outline, edge);
console.log(JSON.stringify({
  same: S.sameStage(S.courseStage([home, home], outline, edge), S.sectorStage(1705, -20, 3225, outline)),
  nearHolds: holds(nearStage, home) && holds(nearStage, near),
  nearLevel: nearStage.win.centre.m,
  midHolds: holds(midStage, home) && holds(midStage, mid),
  midLevel: midStage.win.centre.m,
  far: S.courseStage([home, across], outline, edge),
  apart: S.courseStage([core, high], outline, edge),
  none: S.courseStage([], outline, edge),
}));
""")
    assert out["same"]
    assert out["nearHolds"] and out["nearLevel"] == 1
    assert out["midHolds"] and out["midLevel"] > 1
    assert out["far"] == {"win": None}
    assert out["apart"] == {"win": None}
    assert out["none"] == {"win": None}


def test_a_window_stage_has_a_breadcrumb_down_from_the_galaxy():
    out = _run("""
const centre = P.drillChainOf(900, 4, 2000)[2];
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
    assert out["labels"][0] == "Galaxy" and len(out["labels"]) == 4
    assert [label.split(" ")[0] for label in out["labels"][1:]] == ["Block", "Block", "Block"]
    assert out["labels"][-1].endswith("(3)") and out["last"] == [False, False, False, True]
    assert out["seen"] == [27, 243, "galaxy"]
    assert out["current"] == {"win": None}


def test_an_outline_runs_along_its_blocks_sides_and_closes():
    """MAP.52: an outline lies on its blocks' own sides (a radial edge at a
    block's t0 or t1, a circle at a block's r0 or r1), leaves out the sides
    blocks share, and closes: every corner meets an even number of edges."""
    out = _run("""
const top = S.resolveStage(S.galaxyStage(), outline, edge);
const half = (from, to) => top.view.blocks.filter(b => b.ring >= from && b.ring <= to && b.bounds.t1 <= Math.PI);
const groups = [half(2, 4), half(5, 9), top.view.blocks.filter(b => b.slab === 0 && b.ring > 3 && b.bounds.t0 >= Math.PI / 2 && b.bounds.t1 <= Math.PI)]
  .map(blocks => ({ a0: 0, a1: Math.PI, blocks }));
const results = groups.map(o => {
  const mid = (o.a0 + o.a1) / 2;
  const e = S.outlineEdges(o.blocks, mid);
  const wrap = t => mid + ((((t - mid + Math.PI) % (2 * Math.PI)) + 2 * Math.PI) % (2 * Math.PI)) - Math.PI;
  const near = (a, b) => Math.abs(a - b) < 1e-6;
  const onSide = e.radial.every(r => o.blocks.some(b => (near(wrap(b.bounds.t0), r.t) || near(wrap(b.bounds.t1), r.t))
    && near(b.bounds.r0, r.r0) && near(b.bounds.r1, r.r1)));
  const onRing = e.circles.every(c => o.blocks.some(b => near(b.bounds.r0, c.r) || near(b.bounds.r1, c.r)));
  const ends = new Map();
  // Every bearing meets at the center.
  const q = x => String(Math.round(x * 1e6) + 0);
  const add = (r, t) => { const k = r < 1e-9 ? "center" : q(r / 1e3) + "/" + q(Math.cos(t)) + "/" + q(Math.sin(t)); ends.set(k, (ends.get(k) || 0) + 1); };
  e.radial.forEach(r => { add(r.r0, r.t); add(r.r1, r.t); });
  e.circles.forEach(c => { if (c.t1 - c.t0 < 2 * Math.PI - 1e-9) { add(c.r, c.t0); add(c.r, c.t1); } });
  const closed = Array.from(ends.values()).every(n => n % 2 === 0);
  return {onSide, onRing, closed, radial: e.radial.length, circles: e.circles.length};
});
console.log(JSON.stringify(results));
""")
    assert out
    for result in out:
        assert result["onSide"] and result["onRing"] and result["closed"], result


@pytest.mark.parametrize("text, expected", [
    ("312/-3/1042", {"ring": 312, "layer": -3, "slot": 1042}),
    ("ring 312 layer -3 slot 1042", {"ring": 312, "layer": -3, "slot": 1042}),
    ("  RING 312, LAYER -3, SLOT 1042 ", {"ring": 312, "layer": -3, "slot": 1042}),
    ("0/0/0", {"ring": 0, "layer": 0, "slot": 0}),
])
def test_the_address_bar_reads_an_address(text, expected):
    out = _run(f"console.log(JSON.stringify(S.parseAddress({json.dumps(text)}, edge)));")
    assert out["sector"] == expected
    assert "problem" not in out


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
    assert out["code"]["sector"] == {"ring": 1705, "layer": -20, "slot": 3225}


@pytest.mark.parametrize("text, kind", [
    ("Belcana Scaonon", "name"),
    ("Bead", "name"),          # short enough that it is not a designation
    ("5/0/999", "problem"),    # ring 5 has no slot 999
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

