"""
Tests for `html/static/galaxystages.js`, the Galaxy Map drill-down's
stage rules, run under node (skipped where node isn't installed): the
galaxy's outline and each block's sector total must match a sector-by-
sector count, stage URLs must round-trip, malformed ones must open the
galaxy with a reason, a sector's designation must lead to the stage that
shows it, the breadcrumb must name every step, and the flight must start
and end where it says.
"""

import json
import math
import os
import shutil
import subprocess

import pytest

from stellarObjects.galaxyDensity import build_galaxy_shape

NODE = shutil.which("node")
STATIC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "html", "static")
STAGES = os.path.join(STATIC, "galaxystages.js")
PRISMS = os.path.join(STATIC, "galaxyprisms.js")

pytestmark = pytest.mark.skipif(NODE is None, reason="node is not installed")

SHAPE = build_galaxy_shape(2800.0, 350.0, 200.0, 1.0, 2, math.radians(15.0), 0.4)
EDGE_PC = 4.0
GALAXY_RADIUS_PC = 15000.0


def _shape():
    """The test shape with the galaxy's sector threshold, as the page
    embeds it (lib/galaxymap3d._density_shape)."""
    from stellarObjects.galaxySkeleton import expected_system_count_at_density_1
    from stellarObjects.utils import pc_to_ly

    return {**SHAPE._asdict(), "sector_min_density": 1.0 / expected_system_count_at_density_1(pc_to_ly(EDGE_PC))}


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


def test_children_fit_their_container_and_hold_sectors():
    """stageChildren lists only blocks that can hold sectors, grouped by
    slab, each inside its container's bounds."""
    out = _run("""
const at = P.drillChainOf(400, 1, 1500)[0];
const outer = P.drillBlockBounds(at, edge);
const groups = S.stageChildren(at, outline, edge);
const bad = [];
let blocks = 0;
groups.forEach(g => g.blocks.forEach(b => {
  blocks++;
  const c = b.bounds;
  if (!(b.total > 0) || b.slab !== g.slab || b.m !== at.m / 9) bad.push(b);
  if (c.r0 < outer.r0 - 1e-6 || c.r1 > outer.r1 + 1e-6 || c.z0 < outer.z0 - 1e-6 || c.z1 > outer.z1 + 1e-6) bad.push(b);
}));
const sum = groups.reduce((n, g) => n + g.blocks.reduce((k, b) => k + b.total, 0), 0);
console.log(JSON.stringify({bad, blocks, sum, total: S.drillBlockTotal(at, outline),
  slabs: groups.map(g => g.slab)}));
""")
    assert out["bad"] == []
    assert out["blocks"] > 0
    assert out["sum"] == out["total"]
    assert out["slabs"] == sorted(out["slabs"])


@pytest.mark.parametrize("stage", [
    {"at": None, "picks": []},
    {"at": None, "picks": [{"kind": "quadrant", "n": 1}, {"kind": "layer", "lo": -4, "hi": -2}]},
    {"at": {"m": 243, "ring": 1, "wedge": 2, "slab": 0}, "picks": []},
    {"at": {"m": 27, "ring": 14, "wedge": 30, "slab": -1}, "picks": [{"kind": "layer", "lo": 0, "hi": 0},
                                                                     {"kind": "region", "n": 4}]},
    {"at": {"m": 3, "ring": 130, "wedge": 500, "slab": 2}, "picks": [{"kind": "layer", "lo": 7, "hi": 7}]},
])
def test_stage_urls_round_trip(stage):
    out = _run(f"""
const stage = {json.dumps(stage)};
const query = S.stageQuery(stage);
const parsed = S.parseStageQuery(query);
console.log(JSON.stringify({{query, parsed, same: S.sameStage(parsed.stage, stage)}}));
""")
    assert out["same"], out
    assert out["parsed"]["problem"] is None


def test_an_older_slab_link_reads_as_a_layer_pick():
    out = _run('console.log(JSON.stringify(S.parseStageQuery("?at=27.14.30.-1&slab=2")));')
    assert out["stage"]["picks"] == [{"kind": "layer", "lo": 2, "hi": 2}]
    assert out["problem"] is None


@pytest.mark.parametrize("query", ["?at=5.1.1.0", "?at=junk", "?slab=x", "?sector=ZZZ", "?p=x9"])
def test_malformed_urls_open_the_galaxy_with_a_reason(query):
    out = _run(f"console.log(JSON.stringify(S.parseStageQuery({json.dumps(query)})));")
    assert out["stage"] == {"at": None, "picks": []}
    assert out["problem"]


def test_a_block_or_pick_outside_the_galaxy_is_refused():
    out = _run("""
const high = {m: 3, ring: 300, wedge: 0, slab: 400};
console.log(JSON.stringify([
  S.resolveStage({at: high, picks: []}, outline, edge).problem,
  S.resolveStage({at: null, picks: [{kind: "quadrant", n: 9}]}, outline, edge).problem,
  S.resolveStage({at: null, picks: [{kind: "quadrant", n: 1}, {kind: "layer", lo: 400, hi: 401}]}, outline, edge).problem,
  S.resolveStage({at: null, picks: [{kind: "quadrant", n: 1}]}, outline, edge).problem]));
""")
    assert out[0] and "outside" in out[0]
    assert out[1] and "no quadrant" in out[1]
    assert out[2] and "no layer" in out[2]
    assert out[3] is None


def test_a_sector_far_outside_the_galaxy_parses_but_has_no_stage():
    out = _run("""
const parsed = S.parseStageQuery("?sector=FFFFFFFFFFFFFFFF");
const s = parsed.sector;
console.log(JSON.stringify({sector: s, stage: s ? S.sectorStage(s.ring, s.layer, s.slot, outline, edge) : "none"}));
""")
    assert out["sector"]
    assert out["stage"] is None


def test_the_ladder_is_quarter_then_layer_then_arc():
    """MAP.17/MAP.19 (Boss's "Layer + arc"): the galaxy offers four
    quarters, then a layer (at most three choices, lowest first), then an
    arc of the ring band in view (at most a 3 by 3 grid, sorted into thirds
    of the view), and the choices together hold every block in view. A
    choice's bearings reach as far as its blocks do: near the core a wedge
    is wider than a quarter, and the quarter zooms into that wedge."""
    out = _run("""
const top = S.settleStage({at: null, picks: []}, outline, edge);
const quarter = S.settleStage({at: null, picks: [top.options[1].pick]}, outline, edge);
const layer = S.settleStage({at: null, picks: [top.options[1].pick, quarter.options[1].pick]}, outline, edge);
const covers = r => r.options.reduce((n, o) => n + o.blocks.length, 0) === r.view.blocks.length;
const span = r => r.view.a1 - r.view.a0;
const inside = (o, t) => { const d = ((t - o.a0 + 1e-9) % (2 * Math.PI) + 2 * Math.PI) % (2 * Math.PI); return d <= o.a1 - o.a0 + 2e-9; };
const holds = o => o.blocks.every(b => inside(o, b.bounds.t0) && inside(o, b.bounds.t1));
const tight = o => [o.a0, o.a1].every(a => o.blocks.some(b => Math.abs(Math.cos(b.bounds.t0) - Math.cos(a)) + Math.abs(Math.sin(b.bounds.t0) - Math.sin(a)) < 1e-9 || Math.abs(Math.cos(b.bounds.t1) - Math.cos(a)) + Math.abs(Math.sin(b.bounds.t1) - Math.sin(a)) < 1e-9));
console.log(JSON.stringify({
  wedges: top.options.concat(layer.options).map(o => holds(o) && tight(o)),
  top: [top.kind, top.options.length, covers(top)],
  quarter: [quarter.kind, quarter.options.length, covers(quarter), span(quarter)],
  layers: quarter.options.map(o => [o.pick.lo, o.pick.hi]),
  layer: [layer.kind, layer.options.length, covers(layer), span(layer)],
  arcs: layer.options.map(o => (o.a1 - o.a0) / span(layer)),
}));
""")
    assert out["top"] == ["quadrant", 4, True]
    kind, count, covers, quarter_span = out["quarter"]
    assert kind == "layer" and 2 <= count <= 3 and covers
    assert quarter_span > math.pi / 2
    assert all(out["wedges"])
    los = [lo for lo, _hi in out["layers"]]
    assert los == sorted(los)
    kind, count, covers, _span = out["layer"]
    assert kind == "region" and 2 <= count <= 9 and covers
    assert all(share >= 1 / 3 - 1e-9 for share in out["arcs"])


def test_sector_links_open_the_layer_that_shows_the_sector():
    """MAP.26: a sector opens among its own layer's neighbours, at the
    sector level of the slice holding it."""
    out = _run("""
const ring = 1705, layer = -20, slot = 3225;
const code = S.sectorDesignation(ring, layer, slot);
const parsed = S.parseStageQuery("?sector=" + code);
const stage = S.sectorStage(ring, layer, slot, outline, edge);
const r = S.resolveStage(stage, outline, edge);
const kinds = stage.picks.map(p => p.kind);
console.log(JSON.stringify({code, back: S.parseSectorDesignation(code), parsed, at: stage.at, kinds,
  sectors: r.view.blocks.every(b => b.m === 1 && b.slab === layer),
  holds: r.view.blocks.some(b => b.ring === ring && b.wedge === slot),
  crumbs: S.crumbs(stage, outline, edge).map(c => c.label)}));
""")
    assert out["back"] == {"ring": 1705, "layer": -20, "slot": 3225}
    assert out["parsed"]["sector"] == out["back"]
    assert out["parsed"]["stage"] == {"at": None, "picks": []}
    assert out["at"]["m"] == 3
    assert out["kinds"][-1] == "layer"
    assert out["sectors"] and out["holds"]
    crumbs = out["crumbs"]
    assert crumbs[0] == "Galaxy"
    assert crumbs[1].startswith("Quarter ")
    assert crumbs[-1] == "Layer -20"


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


def test_parent_stages_walk_back_to_the_galaxy():
    out = _run("""
let stage = S.sectorStage(900, 4, 2000, outline, edge);
const crumbs = S.crumbs(stage, outline, edge).length;
let steps = 0;
let repeated = false;
for (let parent = S.parentStage(stage, outline, edge); parent; parent = S.parentStage(stage, outline, edge)) {
  repeated = repeated || S.sameStage(parent, stage);
  stage = parent;
  steps++;
  if (steps > 50) break;
}
const top = S.settleStage({at: null, picks: []}, outline, edge).stage;
console.log(JSON.stringify({steps, crumbs, repeated, last: stage, top}));
""")
    assert not out["repeated"]
    assert out["steps"] == out["crumbs"] - 1
    assert out["last"] == out["top"]


def test_the_course_stage_is_the_smallest_one_holding_every_stop():
    """courseStage (section 9.4) climbs only as far as it must: the same
    sector stays at its layer's sector level, neighbours share a stage
    that holds both, and stops across the galaxy fall back toward the
    galaxy."""
    out = _run("""
const home = {ring: 1705, layer: -20, slot: 3225};
const near = {ring: 1706, layer: -20, slot: 3226};
const core = {ring: 5, layer: 0, slot: 1};
const high = {ring: 5, layer: 900, slot: 1};
const holds = (stage, s) => S.resolveStage(stage, outline, edge).view.blocks.some(b => {
  const chain = P.drillChainOf(s.ring, s.layer, s.slot);
  return b.m === 1 ? b.ring === s.ring && b.wedge === s.slot && b.slab === s.layer : chain.some(c => S.sameBlock(c, b));
});
const nearStage = S.courseStage([home, near], outline, edge);
console.log(JSON.stringify({
  same: S.sameStage(S.courseStage([home, home], outline, edge), S.sectorStage(1705, -20, 3225, outline, edge)),
  nearHolds: holds(nearStage, home) && holds(nearStage, near),
  nearDeep: !!nearStage.at,
  far: S.courseStage([home, core], outline, edge),
  apart: S.courseStage([core, high], outline, edge),
  none: S.courseStage([], outline, edge),
  top: S.settleStage({at: null, picks: []}, outline, edge).stage,
}));
""")
    assert out["same"]
    assert out["nearHolds"] and out["nearDeep"]
    assert out["far"]["at"] is None
    assert out["apart"] == out["top"]
    assert out["none"] == out["top"]


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


def test_a_block_one_sector_tall_shows_its_sectors():
    """A level-27 block whose sectors all lie in one layer skips its
    level-3 blocks: its view is the sectors, picked by arcs."""
    out = _run("""
let found = null;
for (let ring = 0; ring * 27 <= outline.maxRing && !found; ring += 3) {
  for (let slab = 0; slab <= Math.ceil(outline.maxLayer / 27) + 1 && !found; slab++) {
    const at = {m: 27, ring, wedge: 0, slab};
    if (S.drillBlockTotal(at, outline) > 0 && S.thinSectors(at, outline, edge)) found = at;
  }
}
const r = found && S.settleStage({at: found, picks: []}, outline, edge);
const plane = S.thinSectors({m: 27, ring: 10, wedge: 0, slab: 0}, outline, edge);
console.log(JSON.stringify({found, plane, kind: r && r.kind,
  sectors: r && r.view.blocks.every(b => b.m === 1),
  layers: r && new Set(r.view.blocks.map(b => b.slab)).size,
  total: r && r.view.blocks.length, expected: found && S.drillBlockTotal(found, outline)}));
""")
    assert out["found"], "no thin block in the test galaxy"
    assert out["plane"] is None
    assert out["sectors"] and out["layers"] == 1
    assert out["kind"] == "region"
    assert out["total"] == out["expected"]
