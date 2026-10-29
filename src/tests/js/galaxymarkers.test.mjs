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

test("the habitable world kind lists a node's systems that have the reason, most worlds first, and is on until turned off (MAP.303)", () => {
  const kind = M.HABITABLE_KIND;
  assert.equal(kind.defaultOn, true);
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
  const worlds = M.INTEREST_KINDS.filter((k) => k.family !== "class");
  assert.deepEqual(worlds.map((k) => k.interest), ["inhabited", "civilization", "sentient", "life", "habitable", "uninhabitable", "gas-giant", "asteroid-field"]);
  assert.deepEqual(worlds.map((k) => k.icon), ["inhabited-world", undefined, undefined, undefined, "habitable-world", "uninhabitable-world", "gas-giant", "asteroid-field"]);
  assert.deepEqual(worlds.map((k) => k.glyph), [undefined, "civilization", "sentient", "life", undefined, undefined, undefined, undefined]);
  assert.deepEqual(M.INTEREST_KINDS.filter((k) => k.defaultOn).map((k) => k.interest), ["inhabited", "habitable", "gas-giant"]);
  assert.deepEqual(M.INTEREST_KINDS.slice(worlds.length).map((k) => k.family), M.PLANET_CLASS_LETTERS.map(() => "class"),
    "the planet classes come after every other reason");
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

test("each planet class is an optional marker kind of its own: a glyph with its letter, off until chosen, listed only by small blocks (MAP.220)", () => {
  const classes = M.INTEREST_KINDS.filter((k) => k.family === "class");
  assert.deepEqual(classes.map((k) => k.interest), M.PLANET_CLASS_LETTERS.map((letter) => "class-" + letter));
  assert.ok(classes.every((k) => k.defaultOn === false && !k.icon && M.GLYPHS[k.glyph] && k.maxBlock === M.CLASS_MAX_BLOCK));
  assert.equal(classes[0].id, "marker-class-A");
  assert.equal(classes[0].label, "Class A worlds");
  assert.equal(classes[0].describe({ name: "Rho", n: 1 }), "Rho, 1 class A world");
  assert.equal(classes[0].describe({ n: 3 }), "3 class A worlds");
  const registry = M.createMarkerKinds();
  M.INTEREST_KINDS.forEach((kind) => registry.register(kind));
  const systems = [
    { id: 1, name: "Ocean", x: 0, y: 0, z: 0, reasons: { "class-O": { n: 2 }, "class-G": { n: 1 }, habitable: { n: 1, score: 0.8, tier: "Green" } } },
    { id: 2, name: "Rock", x: 1, y: 0, z: 0, reasons: { "class-G": { n: 3 } } },
  ];
  const small = [{ key: "1.0.0.0", node: { interests: systems } }];
  const big = [{ key: "9.0.0.0", node: { interests: systems } }];
  const only = (letter, nodes) => M.selectMarkers({ kinds: registry.kinds(), isOn: (id) => id === "marker-class-" + letter, nodes });
  assert.deepEqual(only("G", small).map((c) => [c.kind.id, c.item.name]).sort(), [["marker-class-G", "Ocean"], ["marker-class-G", "Rock"]]);
  assert.deepEqual(only("O", small).map((c) => c.item.name), ["Ocean"]);
  assert.deepEqual(only("G", big), [], "a block bigger than the classes' limit does not list them, so no mark is drawn from it");
  const off = M.selectMarkers({ kinds: registry.kinds(), isOn: (id) => id === "marker-habitable", nodes: small });
  assert.deepEqual(off.map((c) => [c.kind.id, c.item.name]), [["marker-habitable", "Ocean"]], "classes that are off claim nothing");
  assert.equal(M.interestLines(systems[1]).map((l) => l.text).join("; "), "3 class G worlds");
  assert.deepEqual(M.legendOf(classes[0]), { text: "Class A worlds", swatch: classes[0].color, swatchClass: "map-marker-swatch" });
});

test("worlds with life, with sentient life and with a civilization are three optional kinds with their own sign, beside the habitable worlds (MAP.311)", () => {
  const kinds = ["civilization", "sentient", "life"].map((reason) => M.INTEREST_KINDS.find((k) => k.interest === reason));
  assert.deepEqual(kinds.map((k) => [k.id, k.label, k.defaultOn, k.maxBlock]), [
    ["marker-civilization", "Civilizations", false, 81], ["marker-sentient", "Sentient life", false, 81], ["marker-life", "Life", false, 81]]);
  assert.ok(kinds.every((k) => !k.icon && M.GLYPHS[k.glyph] && k.color));
  assert.equal(new Set(kinds.map((k) => k.color)).size, 3);
  const registry = M.createMarkerKinds();
  M.INTEREST_KINDS.forEach((kind) => registry.register(kind));
  const systems = [
    { id: 1, name: "Green", x: 0, y: 0, z: 0, reasons: { life: { n: 2 }, habitable: { n: 1, score: 0.8, tier: "Green" } } },
    { id: 2, name: "Minds", x: 1, y: 0, z: 0, reasons: { life: { n: 1 }, sentient: { n: 1 } } },
    { id: 3, name: "Cities", x: 2, y: 0, z: 0, reasons: { life: { n: 1 }, sentient: { n: 1 }, civilization: { n: 1 } } },
    { id: 4, name: "Dead", x: 3, y: 0, z: 0, reasons: { uninhabitable: { n: 1, score: 0.1, tier: "Red" } } },
  ];
  const nodes = [{ key: "9.0.0.0", node: { interests: systems } }];
  const only = (kind) => M.selectMarkers({ kinds: registry.kinds(), isOn: (id) => id === "marker-" + kind, nodes }).map((c) => c.item.name).sort();
  assert.deepEqual(only("life"), ["Cities", "Green", "Minds"]);
  assert.deepEqual(only("sentient"), ["Cities", "Minds"]);
  assert.deepEqual(only("civilization"), ["Cities"]);
  const both = M.selectMarkers({ kinds: registry.kinds(), isOn: (id) => id === "marker-life" || id === "marker-sentient", nodes });
  assert.deepEqual(both.map((c) => [c.item.name, c.kind.id]).sort(), [["Cities", "marker-sentient"], ["Green", "marker-life"], ["Minds", "marker-sentient"]],
    "a system is drawn once, by the strongest sign that is on");
  assert.equal(kinds[0].describe({ name: "Cities", n: 2 }), "Cities, 2 worlds with a civilization");
  assert.equal(M.interestLines(systems[2]).map((l) => l.text).join("; "), "1 world with a civilization; 1 world with sentient life; 1 world with life");
});

test("the infobox lines of a system name every reason it has, in the registry's order", () => {
  assert.deepEqual(M.interestLines(SYSTEMS_OF_INTEREST[0]).map((l) => [l.id, l.text]), [
    ["habitable", "2 habitable worlds, best Green"], ["gas-giant", "1 gas giant"]]);
  assert.deepEqual(M.interestLines({ reasons: {} }), []);
  assert.deepEqual(M.interestLines(null), []);
});


const sys = (id, name, x, tier, n) => ({ key: String(id), id, name, x, y: 0, z: 0, n: n || 1, interest: "x", tier });
const mark = (kindId, item) => ({ kind: M.INTEREST_KINDS.find((k) => k.id === kindId), item: Object.assign({}, item, { interest: kindId.slice(7) }) });
const atX = (item) => ({ x: item.x * 10, y: 0 });

test("a crowd of marks becomes one mark with the best reason's symbol, a count and every system", () => {
  const chosen = [
    mark("marker-gas-giant", sys(1, "Jove", 0, undefined, 3)),
    mark("marker-habitable", sys(2, "Eden", 1, 1, 2)),
    mark("marker-inhabited", sys(3, "Rest", 2, 2, 1)),
    mark("marker-habitable", sys(4, "Far", 50, 0, 1)),
  ];
  const merged = M.mergeMarks(chosen, atX, 22);
  assert.equal(merged.length, 2);
  const crowd = merged.find((e) => e.item.merged);
  assert.equal(crowd.kind.id, "marker-inhabited");
  assert.equal(crowd.item.count, 3);
  assert.equal(crowd.item.tier, 2);
  assert.deepEqual(crowd.item.members.map((m) => m.name).sort(), ["Eden", "Jove", "Rest"]);
  assert.deepEqual([crowd.item.x, crowd.item.n], [1, 6]);
  assert.equal(crowd.kind.describe(crowd.item), "3 points of interest");
  assert.equal(merged.find((e) => !e.item.merged).item.name, "Far", "a lone mark stays as it was");
});

test("zooming in spreads the marks and the crowd splits again, in the same order each time", () => {
  const chosen = [mark("marker-habitable", sys(1, "A", 0, 1)), mark("marker-habitable", sys(2, "B", 1, 1))];
  assert.equal(M.mergeMarks(chosen, atX, 22).length, 1);
  const spread = M.mergeMarks(chosen, (item) => ({ x: item.x * 100, y: 0 }), 22);
  assert.deepEqual(spread.map((e) => e.item.name), ["A", "B"]);
  assert.deepEqual(M.mergeMarks(chosen.slice().reverse(), atX, 22)[0].item.key, M.mergeMarks(chosen, atX, 22)[0].item.key);
});

test("a body's symbol and a region's token never merge, and a black hole never merges with a star's mark", () => {
  const bh = { kind: M.BLACK_HOLE_KIND, item: { key: "b", x: 0, y: 0, z: 0 } };
  const body = mark("marker-habitable", Object.assign(sys(5, "Planet", 0, 1), { bodyRef: "planet:1" }));
  const token = mark("marker-inhabited", Object.assign(sys(6, "Region", 0, 1), { regionKey: "9.0.0.0" }));
  const star = mark("marker-habitable", sys(7, "Star", 0, 1));
  const merged = M.mergeMarks([bh, body, token, star], atX, 22);
  assert.equal(merged.length, 4);
  assert.ok(!merged.some((e) => e.item.merged));
});

test("an item with no place on the screen is left as it is", () => {
  const chosen = [mark("marker-habitable", sys(1, "A", 0, 1)), mark("marker-habitable", sys(2, "B", 0, 1))];
  assert.equal(M.mergeMarks(chosen, () => null, 22).length, 2);
});

test("a kind pulled forward is served first from the budget and claims a system's symbol, its entries marked", () => {
  const eden = { id: 5, name: "Eden", x: 1, y: 2, z: 3, reasons: { habitable: { n: 2, score: 0.82, tier: "Green" }, "gas-giant": { n: 1 } } };
  const others = Array.from({ length: 6 }, (_, n) => ({ id: 20 + n, name: "G" + n, x: n, y: 0, z: 0, reasons: { "gas-giant": { n: 1 } } }));
  const node = { interests: [eden].concat(others) };
  const pick = (pulled) => M.selectMarkers({ kinds: M.INTEREST_KINDS, isOn: () => true, isPulled: (id) => pulled.includes(id),
    nodes: [{ key: "1.0.0.0", node }], budget: 3 });
  const plain = pick([]);
  assert.equal(plain.length, 3);
  assert.equal(plain.find((e) => e.item.name === "Eden").kind.interest, "habitable", "the registry's order decides without a pull");
  const forward = pick(["marker-gas-giant"]);
  assert.equal(forward.length, 3);
  assert.ok(forward.every((e) => e.kind.interest === "gas-giant" && e.pulled === true), "the pulled kind has the whole budget and the system's symbol");
  assert.ok(plain.every((e) => !e.pulled));
});

test("a crowd takes the strongest state of its marks: a pulled one gives it its symbol and its place in front", () => {
  const chosen = [
    mark("marker-habitable", sys(1, "A", 0, 1)),
    Object.assign(mark("marker-gas-giant", sys(2, "B", 1, undefined)), { pulled: true }),
  ];
  const [crowd] = M.mergeMarks(chosen, atX, 22);
  assert.equal(crowd.kind.id, "marker-gas-giant");
  assert.equal(crowd.pulled, true);
  const [plain] = M.mergeMarks(chosen.map((e) => ({ kind: e.kind, item: e.item })), atX, 22);
  assert.equal(plain.kind.id, "marker-habitable");
  assert.ok(!plain.pulled);
});

// MAP.339: the near-field rule for black holes.
const view = (extra) => Object.assign({ centre: [0, 0, 0], distance: 10, width: 1200, height: 600, densityPx: 120 }, extra);
const massive = (id, mass, x, extra) => ({ type: "black_hole", id, name: "BH " + id, mass_solar: mass, x, y: 0, z: 0, ...extra });

test("a black hole within 3 view distances is a candidate; one beyond needs its mass to pass the rising floor", () => {
  const rule = M.nearFieldRule(view());
  assert.equal(rule.accepts(massive(1, 6, 29.9)), true, "inside the reach");
  assert.equal(rule.accepts(massive(2, 5, 31)), false, "the lightest hole just past it: the floor has begun to rise");
  assert.equal(rule.accepts(massive(3, 5 * (60 / 30) ** 2 + 0.01, 60)), true, "heavy enough for twice the reach");
  assert.equal(rule.accepts(massive(4, 5 * (60 / 30) ** 2 - 0.5, 60)), false);
  assert.equal(rule.accepts(massive(5, 1e6, 3000)), true, "a supermassive hole stays");
  assert.equal(rule.accepts(massive(6, 1e4, 3000)), false, "an intermediate one needs the square of the distance");
  assert.equal(rule.accepts({ type: "black_hole", id: 7, x: 100, y: 0, z: 0 }), false, "a hole that gives no mass passes only inside the reach");
  assert.equal(rule.accepts({ type: "black_hole", id: 8, x: 5, y: 0, z: 0 }), true);
});

test("the reach and the floor scale with the view distance, so the rule is the same at every zoom", () => {
  const near = M.nearFieldRule(view({ distance: 10 }));
  const far = M.nearFieldRule(view({ distance: 1000 }));
  assert.equal(near.accepts(massive(1, 6, 50)), false);
  assert.equal(far.accepts(massive(1, 6, 50)), true);
  assert.equal(M.nearFieldRule(view({ distance: 10 })).accepts(massive(2, 20, 90)), M.nearFieldRule(view({ distance: 1000 })).accepts(massive(2, 20, 9000)));
});

test("black holes rank by mass over distance squared from the view centre", () => {
  const rule = M.nearFieldRule(view());
  assert.ok(rule.rank(massive(1, 6, 2)) > rule.rank(massive(2, 6, 4)), "nearer wins at the same mass");
  assert.ok(rule.rank(massive(3, 18, 4)) > rule.rank(massive(4, 6, 4)), "heavier wins at the same distance");
  assert.ok(Number.isFinite(rule.rank(massive(5, 6, 0))), "one at the very centre has a rank");
});

test("the cap is one mark per density square of the screen, at least one, and the default square is 120 px", () => {
  assert.equal(M.BLACK_HOLE_DENSITY_PX, 120);
  assert.equal(M.nearFieldRule(view()).cap, 50);
  assert.equal(M.nearFieldRule(view({ densityPx: 60 })).cap, 200);
  assert.equal(M.nearFieldRule(view({ width: 50, height: 50 })).cap, 1);
  assert.equal(M.nearFieldRule(view({ densityPx: 0 })).cap, 50, "no usable density reads as the default");
});

test("selectMarkers keeps the accepted black holes the rule ranks first, up to the cap, and leaves other kinds alone", () => {
  const registry = kinds();
  const points = [];
  for (let i = 0; i < 80; i++) points.push(massive(100 + i, 6 + (i % 5), 1 + i * 0.3));    // all within the reach
  points.push(massive(900, 6, 500), massive(901, 1e6, 3000));
  const nodes = [{ key: "9.1.0.0", node: { points, test: [{ id: "t1", x: 0, y: 0, z: 0, rank: 1 }] } }];
  const near = M.nearFieldRule(view());
  const chosen = M.selectMarkers({ kinds: registry.kinds(), isOn: () => true, nodes, near });
  const holes = chosen.filter((entry) => entry.kind.id === "marker-black-hole").map((entry) => entry.item);
  assert.equal(holes.length, near.cap);
  assert.ok(holes.some((item) => item.id === 901), "the supermassive hole is among them");
  assert.ok(!holes.some((item) => item.id === 900), "a light far hole is not");
  const ranks = holes.map((item) => near.rank(item));
  assert.deepEqual(ranks, ranks.slice().sort((p, q) => q - p), "in the rule's order");
  assert.equal(chosen.filter((entry) => entry.kind.id === "marker-test").length, 1);
  // Without a near-field rule the kind keeps its old ranking.
  const plain = M.selectMarkers({ kinds: registry.kinds(), isOn: () => true, nodes });
  assert.ok(plain.some((entry) => entry.item.id === 900));
});

// MAP.340: crowds of black holes.
const bhMark = (id, mass, x) => ({ kind: M.BLACK_HOLE_KIND, item: { key: String(id), id, name: "BH " + id, type: "black_hole", mass_solar: mass, x, y: 0, z: 0 } });

test("black holes within the merge distance are one mark with their count, the heaviest first, and spread out they split", () => {
  const chosen = [bhMark(1, 6, 0), bhMark(2, 12, 1), bhMark(3, 8, 2), bhMark(4, 7, 200)];
  const merged = M.mergeMarks(chosen, atX, 22);
  assert.equal(merged.length, 2);
  const crowd = merged.find((entry) => entry.item.merged);
  assert.equal(crowd.item.count, 3);
  assert.equal(crowd.item.holes, true);
  assert.deepEqual(crowd.item.members.map((item) => item.id), [2, 3, 1], "the heaviest first");
  assert.equal(crowd.item.key, "merge:2");
  assert.deepEqual([crowd.item.x, crowd.item.y, crowd.item.z], [1, 0, 0], "the middle of the holes");
  assert.equal(crowd.kind.id, "marker-black-hole");
  assert.equal(crowd.kind.describe(crowd.item), "3 black holes");
  assert.ok(merged.some((entry) => entry.item.id === 4 && !entry.item.merged), "the far one stays itself");
  const spread = M.mergeMarks(chosen, (item) => ({ x: item.x * 100, y: 0 }), 22);
  assert.equal(spread.length, 4);
  assert.ok(!spread.some((entry) => entry.item.merged));
});

test("a crowd of black holes and a crowd of systems on the same spot stay two marks", () => {
  const chosen = [bhMark(1, 6, 0), bhMark(2, 7, 0.3), mark("marker-habitable", sys(5, "A", 0, 1)), mark("marker-habitable", sys(6, "B", 0.1, 1))];
  const merged = M.mergeMarks(chosen, atX, 22);
  assert.equal(merged.length, 2);
  assert.deepEqual(merged.map((entry) => entry.item.holes === true).sort(), [false, true]);
  assert.deepEqual(merged.map((entry) => entry.item.count).sort(), [2, 2]);
});

test("the same crowd of black holes merges the same way whatever the order they come in", () => {
  const chosen = [bhMark(1, 6, 0), bhMark(2, 6, 1), bhMark(3, 6, 2)];
  const forward = M.mergeMarks(chosen, atX, 22);
  const backward = M.mergeMarks(chosen.slice().reverse(), atX, 22);
  assert.equal(forward[0].item.key, backward[0].item.key);
  assert.deepEqual(forward[0].item.members.map((item) => item.id), backward[0].item.members.map((item) => item.id));
});

test("black holes at very different depths on one line of sight are not a crowd, and ones at about the same depth are", () => {
  const chosen = [bhMark(1, 6, 0), bhMark(2, 7, 1), bhMark(3, 8, 2)];
  const depths = { 1: 100, 2: 400, 3: 120 };
  const merged = M.mergeMarks(chosen, (item) => ({ x: item.x * 10, y: 0, depth: depths[item.id] }), 22);
  assert.equal(merged.length, 2);
  const crowd = merged.find((entry) => entry.item.merged);
  assert.deepEqual(crowd.item.members.map((item) => item.id).sort(), [1, 3]);
  const same = M.mergeMarks(chosen, (item) => ({ x: item.x * 10, y: 0, depth: 100 + item.id }), 22);
  assert.equal(same.length, 1);
});
