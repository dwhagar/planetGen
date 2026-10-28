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
  registry.register({ id: "marker-test", label: "Test", glyph: "black-hole", color: "#0f0", maxBlock: 9,
    source: (node) => (node.test || []).map((t) => ({ key: t.id, ...t })), rank: (item) => item.rank || 0 });
  return registry;
}

test("a kind needs an id, a label and a known glyph, and is registered once", () => {
  const registry = M.createMarkerKinds();
  assert.throws(() => registry.register({ id: "x", label: "X", glyph: "nope" }));
  assert.throws(() => registry.register({ label: "X", glyph: "black-hole" }));
  registry.register({ id: "x", label: "X", glyph: "black-hole" });
  assert.throws(() => registry.register({ id: "x", label: "X", glyph: "black-hole" }));
  assert.equal(registry.get("x").defaultOn, true);
  assert.equal(registry.get("y"), null);
});

test("black holes come from the nodes' points, and only black holes", () => {
  const items = M.BLACK_HOLE_KIND.source({ points: [hole(1, 5), star(2), hole("s3", 9)] });
  assert.deepEqual(items.map((i) => i.key), ["1", "s3"]);
  assert.equal(M.BLACK_HOLE_KIND.describe(items[0]), "BH 1, Black Hole");
  assert.deepEqual(M.BLACK_HOLE_KIND.source({}), []);
});

test("a kind that is off draws nothing and a node of a bigger block than its maximum lists nothing for it", () => {
  const registry = kinds();
  const nodes = [{ key: "81.0.0.0", node: { points: [hole(1, 5)], test: [{ id: "a", x: 0, y: 0, z: 0 }] } },
                 { key: "3.0.0.0", node: { test: [{ id: "b", x: 0, y: 0, z: 0 }] } }];
  const all = M.selectMarkers({ kinds: registry.kinds(), isOn: () => true, nodes });
  assert.deepEqual(all.map((c) => c.kind.id + ":" + c.item.key).sort(), ["marker-black-hole:1", "marker-test:b"]);
  const holesOff = M.selectMarkers({ kinds: registry.kinds(), isOn: (id) => id !== "marker-black-hole", nodes });
  assert.deepEqual(holesOff.map((c) => c.item.key), ["b"]);
});

test("an item listed by several nodes is drawn once and a missing node is skipped", () => {
  const registry = kinds();
  const nodes = [{ key: "27.0.0.0", node: { points: [hole(1, 5)] } }, { key: "9.0.0.0", node: { points: [hole(1, 5)] } },
                 { key: "9.0.0.1", node: undefined }];
  assert.equal(M.selectMarkers({ kinds: registry.kinds(), isOn: () => true, nodes }).length, 1);
});

test("items out of view are left out, and the budget keeps the highest rank", () => {
  const registry = kinds();
  const points = [hole(1, 1, 0), hole(2, 9, 0), hole(3, 5, 0), hole(4, 7, 100)];
  const nodes = [{ key: "9.0.0.0", node: { points } }];
  const chosen = M.selectMarkers({ kinds: registry.kinds(), isOn: () => true, nodes, budget: 2, inView: (item) => item.x < 50 });
  assert.deepEqual(chosen.map((c) => c.item.key), ["2", "3"]);
});

test("kinds share the budget evenly and a kind with less hands the rest on", () => {
  const registry = kinds();
  const many = Array.from({ length: 10 }, (_, i) => ({ id: "t" + i, rank: i, x: 0, y: 0, z: 0 }));
  const nodes = [{ key: "3.0.0.0", node: { points: [hole(1, 1), hole(2, 2)], test: many } }];
  const chosen = M.selectMarkers({ kinds: registry.kinds(), isOn: () => true, nodes, budget: 6 });
  const count = (id) => chosen.filter((c) => c.kind.id === id).length;
  assert.equal(count("marker-black-hole"), 2);
  assert.equal(count("marker-test"), 4);
  const tight = M.selectMarkers({ kinds: registry.kinds(), isOn: () => true, nodes, budget: 4 });
  assert.equal(tight.filter((c) => c.kind.id === "marker-black-hole").length, 2);
  assert.equal(tight.filter((c) => c.kind.id === "marker-test").length, 2);
});

test("an item with its own size is hidden while it is under two pixels", () => {
  const registry = M.createMarkerKinds();
  registry.register({ id: "marker-cloud", label: "Clouds", glyph: "black-hole", extentPc: (item) => item.radius,
    source: (node) => node.clouds.map((c) => ({ key: c.id, ...c })) });
  const nodes = [{ key: "243.0.0.0", node: { clouds: [{ id: "big", radius: 10, x: 0, y: 0, z: 0 }, { id: "small", radius: 1, x: 0, y: 0, z: 0 }] } }];
  const chosen = M.selectMarkers({ kinds: registry.kinds(), isOn: () => true, nodes, pxOf: (item) => item.radius * 0.5 });
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

const SYSTEMS_OF_INTEREST = [
  { id: 5, name: "Eden", x: 1, y: 2, z: 3, reasons: { habitable: { n: 2, score: 0.82, tier: "Green" }, "gas-giant": { n: 1 } } },
  { id: 6, name: "Rest", x: 0, y: 0, z: 0, reasons: { habitable: { n: 1, score: 0.55, tier: "Yellow" } } },
  { id: 7, name: "Jove", x: 4, y: 0, z: 0, reasons: { "gas-giant": { n: 3 } } },
];

test("the habitable world kind lists a node's systems that have the reason, most worlds first, and is off until asked for", () => {
  const kind = M.HABITABLE_KIND;
  assert.equal(kind.defaultOn, false);
  const node = { interests: SYSTEMS_OF_INTEREST };
  const items = kind.source(node);
  assert.deepEqual(items.map((i) => i.key), ["5", "6"]);
  assert.deepEqual(items.map((i) => [i.n, i.tier]), [[2, 1], [1, 2]]);      // the tier is an index into the colours: Green 1, Yellow 2
  assert.deepEqual(kind.source({}), []);
  assert.ok(kind.rank({ n: 3 }) > kind.rank({ n: 1 }));
  assert.equal(kind.describe(items[0]), "Eden, 2 habitable worlds, best Green");
  assert.equal(kind.describe(items[1]), "Rest, 1 habitable world, best Yellow");
  const registry = M.createMarkerKinds();
  registry.register(kind);
  const chosen = M.selectMarkers({ kinds: registry.kinds(), isOn: (id) => id === "marker-habitable", nodes: [{ key: "9.0.0.0", node }] });
  assert.deepEqual(chosen.map((c) => c.item.name), ["Eden", "Rest"]);
});

test("the gas giant kind lists the systems with a gas giant", () => {
  const kind = M.GAS_GIANT_KIND;
  assert.equal(kind.id, "marker-gas-giant");
  assert.equal(kind.icon, "gas-giant");
  assert.deepEqual(kind.source({ interests: SYSTEMS_OF_INTEREST }).map((i) => [i.name, i.n]), [["Eden", 1], ["Jove", 3]]);
  assert.equal(kind.describe({ name: "Jove", n: 3 }), "Jove, 3 gas giants");
  assert.equal(kind.describe({ name: "Eden", n: 1 }), "Eden, 1 gas giant");
});

test("a system with two reasons is drawn once, by the first reason that is on, and a reason that is off does not claim it", () => {
  const registry = M.createMarkerKinds();
  M.INTEREST_KINDS.forEach((kind) => registry.register(kind));
  const nodes = [{ key: "9.0.0.0", node: { interests: SYSTEMS_OF_INTEREST } }];
  const both = M.selectMarkers({ kinds: registry.kinds(), isOn: () => true, nodes });
  assert.deepEqual(both.map((c) => [c.kind.id, c.item.name]).sort(), [
    ["marker-gas-giant", "Jove"], ["marker-habitable", "Eden"], ["marker-habitable", "Rest"]]);
  const giantsOnly = M.selectMarkers({ kinds: registry.kinds(), isOn: (id) => id === "marker-gas-giant", nodes });
  assert.deepEqual(giantsOnly.map((c) => c.item.name).sort(), ["Eden", "Jove"]);
});

test("the reasons are in the order of the best world's symbol and a barren world or a belt is a reason of its own (MAP.321)", () => {
  assert.deepEqual(M.INTEREST_KINDS.map((k) => k.interest), ["inhabited", "habitable", "uninhabitable", "gas-giant", "asteroid-field"]);
  assert.deepEqual(M.INTEREST_KINDS.map((k) => k.icon), ["inhabited-world", "habitable-world", "uninhabitable-world", "gas-giant", "asteroid-field"]);
  assert.ok(M.INTEREST_KINDS.every((k) => k.defaultOn === false));
  const systems = [
    { id: 1, name: "Rock", x: 0, y: 0, z: 0, reasons: { uninhabitable: { n: 2, score: 0.1, tier: "Red" } } },
    { id: 2, name: "Belt", x: 1, y: 0, z: 0, reasons: { "asteroid-field": { n: 1 } } },
    { id: 3, name: "Mixed", x: 2, y: 0, z: 0, reasons: { "asteroid-field": { n: 1 }, "gas-giant": { n: 1 }, uninhabitable: { n: 1, score: 0.2, tier: "Red" } } },
  ];
  const registry = M.createMarkerKinds();
  M.INTEREST_KINDS.forEach((kind) => registry.register(kind));
  const nodes = [{ key: "9.0.0.0", node: { interests: systems } }];
  const chosen = M.selectMarkers({ kinds: registry.kinds(), isOn: () => true, nodes });
  assert.deepEqual(chosen.map((c) => [c.item.name, c.kind.id]).sort(), [
    ["Belt", "marker-asteroid-field"], ["Mixed", "marker-uninhabitable"], ["Rock", "marker-uninhabitable"]]);
  const onlyBelts = M.selectMarkers({ kinds: registry.kinds(), isOn: (id) => id === "marker-asteroid-field", nodes });
  assert.deepEqual(onlyBelts.map((c) => c.item.name).sort(), ["Belt", "Mixed"]);
  const settled = { id: 4, name: "Home", x: 3, y: 0, z: 0, reasons: { habitable: { n: 2, score: 0.8, tier: "Green" }, inhabited: { n: 1, score: 0.8, tier: "Green" } } };
  assert.deepEqual(M.selectMarkers({ kinds: registry.kinds(), isOn: () => true, nodes: [{ key: "9.0.0.0", node: { interests: [settled] } }] })
    .map((c) => c.kind.id), ["marker-inhabited"], "an inhabited world's symbol comes before a habitable one's");
  assert.equal(M.interestLines(settled).map((l) => l.text).join("; "), "1 inhabited world, best Green; 2 habitable worlds, best Green");
  assert.equal(M.interestLines(systems[2]).map((l) => l.text).join("; "), "1 uninhabitable world; 1 gas giant; 1 asteroid belt");
});

test("the infobox lines of a system name every reason it has, in the registry's order", () => {
  assert.deepEqual(M.interestLines(SYSTEMS_OF_INTEREST[0]).map((l) => [l.id, l.text]), [
    ["habitable", "2 habitable worlds, best Green"], ["gas-giant", "1 gas giant"]]);
  assert.deepEqual(M.interestLines({ reasons: {} }), []);
  assert.deepEqual(M.interestLines(null), []);
});

