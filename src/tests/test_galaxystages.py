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
    {"at": None, "slab": None},
    {"at": None, "slab": 0},
    {"at": {"m": 243, "ring": 1, "wedge": 2, "slab": 0}, "slab": None},
    {"at": {"m": 27, "ring": 14, "wedge": 30, "slab": -1}, "slab": 0},
    {"at": {"m": 3, "ring": 130, "wedge": 500, "slab": 2}, "slab": 7},
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


@pytest.mark.parametrize("query", ["?at=5.1.1.0", "?at=junk", "?slab=x", "?sector=ZZZ"])
def test_malformed_urls_open_the_galaxy_with_a_reason(query):
    out = _run(f"console.log(JSON.stringify(S.parseStageQuery({json.dumps(query)})));")
    assert out["stage"] == {"at": None, "slab": None}
    assert out["problem"]


def test_a_block_outside_the_galaxy_is_refused():
    out = _run("""
const far = {m: 243, ring: outline.maxRing, wedge: 0, slab: 0};
const high = {m: 3, ring: 300, wedge: 0, slab: 400};
console.log(JSON.stringify([S.validStage({at: high, slab: null}, outline, edge),
  S.validStage({at: null, slab: 99}, outline, edge), S.validStage({at: null, slab: 0}, outline, edge)]));
""")
    assert out[0] and "outside" in out[0]
    assert out[1] and "empty" in out[1]
    assert out[2] is None


def test_a_sector_far_outside_the_galaxy_parses_but_is_refused():
    out = _run("""
const parsed = S.parseStageQuery("?sector=FFFFFFFFFFFFFFFF");
console.log(JSON.stringify(S.validStage(parsed.stage, outline, edge)));
""")
    assert "outside" in out


def test_sector_links_open_the_layer_that_shows_the_sector():
    out = _run("""
const ring = 1705, layer = -20, slot = 3225;
const code = S.sectorDesignation(ring, layer, slot);
const parsed = S.parseStageQuery("?sector=" + code);
const children = S.stageChildren(parsed.stage.at, outline, edge)
  .filter(g => g.slab === parsed.stage.slab)
  .flatMap(g => g.blocks)
  .some(b => b.ring === ring && b.wedge === slot && b.slab === layer);
console.log(JSON.stringify({code, back: S.parseSectorDesignation(code), parsed, children,
  number: S.stageNumber(parsed.stage), crumbs: S.crumbs(parsed.stage).map(c => c.label)}));
""")
    assert out["back"] == {"ring": 1705, "layer": -20, "slot": 3225}
    assert out["parsed"]["sector"] == out["back"]
    assert out["number"] == 8
    assert out["children"]
    crumbs = out["crumbs"]
    assert len(crumbs) == 8
    assert crumbs[0] == "Galaxy"
    assert crumbs[-1] == "Layer -20"
    assert crumbs[1].startswith("Slab ") and crumbs[2].startswith("Block ")


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
const stage = S.sectorStage(900, 4, 2000);
const chain = S.stageChain(stage);
console.log(JSON.stringify({numbers: chain.map(S.stageNumber), top: S.parentStage({at: null, slab: null})}));
""")
    assert out["numbers"] == [1, 2, 3, 4, 5, 6, 7, 8]
    assert out["top"] is None


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
