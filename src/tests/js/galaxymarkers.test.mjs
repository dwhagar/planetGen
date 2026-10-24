// tests/js/galaxymarkers.test.mjs -- static/galaxymarkers.js (MAP.219): the Galaxy Map's marker kinds and what they choose to draw.

import { test } from "node:test";
import assert from "node:assert/strict";

import { STATIC_URL } from "./fakedom.mjs";

const M = await import(new URL("galaxymarkers.js", STATIC_URL).href);

const hole = (id, mass, x = 0) => ({ type: "black_hole", id, name: "BH " + id, map_luminosity_sol: mass, x, y: 0, z: 0 });
const star = (id) => ({ type: "neutron_star", id, x: 1, y: 1, z: 1 });

function kinds() {
  const registry = M.createMarkerKinds();
  registry.register(M.BLACK_HOLE_KIND);
  registry.register({ id: "marker-test", label: "Test", glyph: "habitable", color: "#0f0", minLevel: 3,
    source: (tile) => (tile.test || []).map((t) => ({ key: t.id, ...t })), rank: (item) => item.rank || 0 });
  return registry;
}

test("a kind needs an id, a label and a known glyph, and is registered once", () => {
  const registry = M.createMarkerKinds();
  assert.throws(() => registry.register({ id: "x", label: "X", glyph: "nope" }));
  assert.throws(() => registry.register({ label: "X", glyph: "habitable" }));
  registry.register({ id: "x", label: "X", glyph: "habitable" });
  assert.throws(() => registry.register({ id: "x", label: "X", glyph: "habitable" }));
  assert.equal(registry.get("x").defaultOn, true);
  assert.equal(registry.get("y"), null);
});

test("black holes come from the tiles' points, and only black holes", () => {
  const items = M.BLACK_HOLE_KIND.source({ points: [hole(1, 5), star(2), hole("s3", 9)] });
  assert.deepEqual(items.map((i) => i.key), ["1", "s3"]);
  assert.equal(M.BLACK_HOLE_KIND.describe(items[0]), "BH 1, Black Hole");
  assert.deepEqual(M.BLACK_HOLE_KIND.source({}), []);
});

test("a kind that is off draws nothing and a tile coarser than its minimum level lists nothing for it", () => {
  const registry = kinds();
  const tiles = [{ key: "1/0/0/0", tile: { points: [hole(1, 5)], test: [{ id: "a", x: 0, y: 0, z: 0 }] } },
                 { key: "4/0/0/0", tile: { test: [{ id: "b", x: 0, y: 0, z: 0 }] } }];
  const all = M.selectMarkers({ kinds: registry.kinds(), isOn: () => true, tiles });
  assert.deepEqual(all.map((c) => c.kind.id + ":" + c.item.key).sort(), ["marker-black-hole:1", "marker-test:b"]);
  const holesOff = M.selectMarkers({ kinds: registry.kinds(), isOn: (id) => id !== "marker-black-hole", tiles });
  assert.deepEqual(holesOff.map((c) => c.item.key), ["b"]);
});

test("an item listed by several tiles is drawn once and a missing tile is skipped", () => {
  const registry = kinds();
  const tiles = [{ key: "2/0/0/0", tile: { points: [hole(1, 5)] } }, { key: "3/0/0/0", tile: { points: [hole(1, 5)] } },
                 { key: "3/0/0/1", tile: undefined }];
  assert.equal(M.selectMarkers({ kinds: registry.kinds(), isOn: () => true, tiles }).length, 1);
});

test("items out of view are left out, and the budget keeps the highest rank", () => {
  const registry = kinds();
  const points = [hole(1, 1, 0), hole(2, 9, 0), hole(3, 5, 0), hole(4, 7, 100)];
  const tiles = [{ key: "3/0/0/0", tile: { points } }];
  const chosen = M.selectMarkers({ kinds: registry.kinds(), isOn: () => true, tiles, budget: 2, inView: (item) => item.x < 50 });
  assert.deepEqual(chosen.map((c) => c.item.key), ["2", "3"]);
});

test("kinds share the budget evenly and a kind with less hands the rest on", () => {
  const registry = kinds();
  const many = Array.from({ length: 10 }, (_, i) => ({ id: "t" + i, rank: i, x: 0, y: 0, z: 0 }));
  const tiles = [{ key: "4/0/0/0", tile: { points: [hole(1, 1), hole(2, 2)], test: many } }];
  const chosen = M.selectMarkers({ kinds: registry.kinds(), isOn: () => true, tiles, budget: 6 });
  const count = (id) => chosen.filter((c) => c.kind.id === id).length;
  assert.equal(count("marker-black-hole"), 2);
  assert.equal(count("marker-test"), 4);
  const tight = M.selectMarkers({ kinds: registry.kinds(), isOn: () => true, tiles, budget: 4 });
  assert.equal(tight.filter((c) => c.kind.id === "marker-black-hole").length, 2);
  assert.equal(tight.filter((c) => c.kind.id === "marker-test").length, 2);
});

test("an item with its own size is hidden while it is under two pixels", () => {
  const registry = M.createMarkerKinds();
  registry.register({ id: "marker-cloud", label: "Clouds", glyph: "habitable", extentPc: (item) => item.radius,
    source: (tile) => tile.clouds.map((c) => ({ key: c.id, ...c })) });
  const tiles = [{ key: "0/0/0/0", tile: { clouds: [{ id: "big", radius: 10, x: 0, y: 0, z: 0 }, { id: "small", radius: 1, x: 0, y: 0, z: 0 }] } }];
  const chosen = M.selectMarkers({ kinds: registry.kinds(), isOn: () => true, tiles, pxOf: (item) => item.radius * 0.5 });
  assert.deepEqual(chosen.map((c) => c.item.key), ["big"]);
  assert.equal(M.MIN_EXTENT_PX, 2);
});

test("the legend line carries the kind's color", () => {
  assert.deepEqual(M.legendOf(M.BLACK_HOLE_KIND), { text: "Black holes", swatch: M.BLACK_HOLE_KIND.color, swatchClass: "map-marker-swatch" });
});

test("every glyph draws without error on a recording context", () => {
  const calls = [];
  const ctx = new Proxy({}, { get: (_t, name) => (...args) => { calls.push(name); }, set: () => true });
  Object.keys(M.GLYPHS).forEach((name) => M.GLYPHS[name](ctx, 64, "#fff", "#000"));
  assert.ok(calls.includes("arc") && calls.includes("stroke"));
});

test("the habitable world kind lists a tile's systems, most worlds first, and is off until asked for", () => {
  const kind = M.HABITABLE_KIND;
  assert.equal(kind.defaultOn, false);
  const tile = { worlds: [{ id: 5, name: "Eden", worlds: 2, x: 1, y: 2, z: 3 }, { id: 6, name: "Rest", worlds: 1, x: 0, y: 0, z: 0 }] };
  assert.deepEqual(kind.source(tile).map((i) => i.key), ["5", "6"]);
  assert.deepEqual(kind.source({}), []);
  assert.ok(kind.rank({ worlds: 3 }) > kind.rank({ worlds: 1 }));
  assert.equal(kind.describe({ name: "Eden", worlds: 2 }), "Eden, 2 habitable worlds");
  assert.equal(kind.describe({ name: "Rest", worlds: 1 }), "Rest, 1 habitable world");
  const registry = M.createMarkerKinds();
  registry.register(kind);
  const chosen = M.selectMarkers({ kinds: registry.kinds(), isOn: (id) => id === "marker-habitable", tiles: [{ key: "9/0/0/0", tile }] });
  assert.deepEqual(chosen.map((c) => c.item.name), ["Rest", "Eden"].reverse());
});
