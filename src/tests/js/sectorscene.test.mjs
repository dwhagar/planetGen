// tests/js/sectorscene.test.mjs -- static/sectorscene.js, the Sector Map's
// scene (TEST.58, moved from the old sectormap.js tests with MAP.68): the
// entries it draws, the "Show on map" keys, the screen-reader labels and
// the rogue-planet markers (rings, a bigger point, an easier pick). The
// scene is lib/starmap.py's real JSON (test_js_unit.py's fixtures); the
// picking and the page around it are the browser tests' job
// (test_web_browser_fixture_maps.py).

import { test } from "node:test";
import assert from "node:assert/strict";

import { installDom, STATIC_URL } from "./fakedom.mjs";

const F = JSON.parse(process.env.PLANETGEN_JS_FIXTURES || "null");
if (!F) throw new Error("run through test_js_unit.py (PLANETGEN_JS_FIXTURES is not set)");

installDom("http://localhost/sector/1");
const SS = await import(new URL("sectorscene.js", STATIC_URL).href);

function build() {
  const data = structuredClone(F.sectorMap.data);
  const sector = SS.buildSectorScene(data, { origin: [0, 0, 0], unit: 1, flipY: false, pixelRatio: 1 });
  return { data, sector };
}

test("every star, cloud and neighbor is an entry, in list order", () => {
  const { data, sector } = build();
  assert.equal(sector.entries.length, data.stars.length + data.clouds.length + data.neighbors.length);
  assert.deepEqual(sector.entries.slice(0, data.stars.length).map((e) => e.name), data.stars.map((s) => s.name));
});

test("an entry is found by its key, only if it is on the map", () => {
  const { data, sector } = build();
  const nebula = data.clouds.find((c) => c.key === "nebula:1");
  assert.equal(sector.entryByKey.get("nebula:1").name, nebula.name);
  assert.equal(sector.entryByKey.get("nebula:999"), undefined);
});

test("every entry has a screen-reader label", () => {
  const { sector } = build();
  for (const entry of sector.entries) assert.ok(SS.entryLabel(entry), "a label");
  const neighbor = sector.entries.find((e) => e.isNeighbor);
  assert.equal(SS.entryLabel(neighbor), neighbor.exists ? neighbor.name || "Unnamed sector" : "Not yet generated (" + neighbor.designation + ")");
});

test("a rogue planet is one dim point with no ring or rogue-only mark; the kind's Mark grows it (MAP.184)", () => {
  const { data, sector } = build();
  const all = [];
  sector.group.traverse((o) => all.push(o));
  assert.equal(all.filter((o) => o.isSprite && o.material.sizeAttenuation === false).length, 0, "no ring");
  assert.equal(sector.setRoguesMarked, undefined);
  assert.equal(sector.roguesMarked, undefined);
  const points = all.find((o) => o.isPoints);
  const rogue = data.clouds.find((c) => c.kind === "roguePlanet");
  assert.equal(rogue.markedLight, undefined);
  const index = data.stars.length + data.clouds.filter((c) => c.light).indexOf(rogue);
  const size = () => points.geometry.getAttribute("pointSize").array[index];
  assert.equal(size(), rogue.light.sizePx);
  sector.setKindMarked("roguePlanet", true);
  assert.ok(size() > rogue.light.sizePx, "marked as a kind, it grows");
  sector.setKindMarked("roguePlanet", false);
  assert.equal(size(), rogue.light.sizePx);
});

test("the kinds in a sector are listed in order, and a hidden kind goes from the map (MAP.79)", () => {
  const { data, sector } = build();
  const kinds = sector.kinds();
  assert.deepEqual(kinds.map((k) => k.kind), ["star", "nebula", "roguePlanet", "neighbor"].filter((k) => kinds.some((x) => x.kind === k)));
  assert.equal(kinds.find((k) => k.kind === "star").count, data.stars.length);
  assert.equal(SS.kindOf(data.stars[0]), "star");
  assert.equal(SS.kindOf({ kind: "blackHoleQuiescent" }), "blackHole");

  const all = [];
  sector.group.traverse((o) => all.push(o));
  const points = all.find((o) => o.isPoints);
  const sizes = () => Array.from(points.geometry.getAttribute("pointSize").array);
  const before = sizes();
  const pointLayer = sector.layers.find((l) => l.name === "sector-points");
  const bodies = sector.layers.find((l) => l.name === "sector-bodies");
  const volumes = sector.layers.find((l) => l.name === "sector-volumes");
  const pickable = () => pointLayer.points().length + bodies.meshes().length + volumes.meshes().length;
  const total = pickable();

  sector.setKindHidden("star", true);
  assert.equal(sector.kindHidden("star"), true);
  data.stars.forEach((_, i) => assert.equal(sizes()[i], 0, "a hidden star draws nothing"));
  assert.ok(pointLayer.points().every((e) => SS.kindOf(e) !== "star"), "and can't be picked");
  sector.setKindHidden("star", false);
  assert.deepEqual(sizes(), before, "shown again, as it was");

  for (const kind of ["nebula", "neighbor", "roguePlanet"]) {
    if (!kinds.some((k) => k.kind === kind)) continue;
    sector.setKindHidden(kind, true);
    assert.ok(pickable() < total, kind + " is no longer pickable");
    sector.entries.filter((e) => SS.kindOf(e) === kind).forEach((entry) => {
      assert.equal(pointLayer.points().includes(entry), false);
    });
    sector.setKindHidden(kind, false);
    assert.equal(pickable(), total, kind + " is back");
  }
});

test("a hidden rogue planet is gone from the map and comes back as it was (MAP.79)", () => {
  const { data, sector } = build();
  const all = [];
  sector.group.traverse((o) => all.push(o));
  const points = all.find((o) => o.isPoints);
  const rogue = data.clouds.find((c) => c.kind === "roguePlanet");
  const index = data.stars.length + data.clouds.filter((c) => c.light).indexOf(rogue);
  sector.setKindMarked("roguePlanet", true);
  const marked = points.geometry.getAttribute("pointSize").array[index];
  sector.setKindHidden("roguePlanet", true);
  assert.equal(points.geometry.getAttribute("pointSize").array[index], 0);
  sector.setKindHidden("roguePlanet", false);
  assert.equal(points.geometry.getAttribute("pointSize").array[index], marked, "still marked");
});

test("an object left in an uncharted sector says so and has no page, bookmark or course buttons (MAP.162)", () => {
  const star = { name: "B star (uncharted)", starType: "B", temp: "20000 K", uncharted: true };
  const cloud = { name: "Black hole (uncharted)", kind: "blackHoleQuiescent", typeLabel: "Black Hole (Stellar)",
    distanceText: "1 ly from sector center", key: "black_hole:s5", uncharted: true };
  const nav = { active: () => false, actionsFor: () => [{ label: "Start here" }] };
  [star, cloud].forEach((entry) => {
    const spec = SS.infoSpec(entry, {}, nav);
    assert.equal(spec.title, entry.name);
    assert.deepEqual(spec.links, []);
    assert.deepEqual(spec.nav, []);
    assert.ok(!spec.bookmark);
    assert.ok(spec.fields.some((field) => field[0] === "Status" && /Uncharted/.test(field[1])));
  });
});
