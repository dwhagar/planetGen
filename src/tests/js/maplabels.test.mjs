// tests/js/maplabels.test.mjs -- static/maplabels.js (MAP.229): the label placement engine and its DOM layer.

import { test } from "node:test";
import assert from "node:assert/strict";

import { installDom, STATIC_URL } from "./fakedom.mjs";

installDom("http://localhost/galaxy");
const M = await import(new URL("maplabels.js", STATIC_URL).href);

const view = { width: 800, height: 600 };
// A flat projection: the anchor's x and y are its screen position, z its depth.
const project = (x, y, z) => (z < 0 ? null : { x, y, depth: z });
const measure = (text) => text.length * 7;

function label(id, x, y, extra) {
  return Object.assign({ id, kind: "sector", text: id, x, y, z: 1, priority: 1 }, extra || {});
}

function place(candidates, extra) {
  return M.placeLabels(candidates, Object.assign({ project, measure, viewport: view }, extra || {}));
}

function overlap(a, b) {
  return a.left < b.left + b.width && b.left < a.left + a.width && a.top < b.top + b.height && b.top < a.top + a.height;
}

test("labels outside the view or behind the camera are dropped", () => {
  const placed = place([label("in", 400, 300), label("left", -20, 300), label("below", 400, 700), label("behind", 400, 300, { z: -1 })]);
  assert.deepEqual(placed.map((p) => p.candidate.id), ["in"]);
});

test("no two placed labels overlap, and the higher priority wins the spot", () => {
  const placed = place([label("minor", 400, 300, { priority: 1 }), label("major", 405, 302, { priority: 5 }), label("far", 100, 100)]);
  const ids = placed.map((p) => p.candidate.id);
  assert.deepEqual(ids, ["major", "far"], "the lower priority one near it is the one left out");
  for (let i = 0; i < placed.length; i++) for (let j = i + 1; j < placed.length; j++) assert.equal(overlap(placed[i], placed[j]), false);
});

test("many labels never overlap and the budget caps how many are placed", () => {
  const many = [];
  for (let i = 0; i < 400; i++) many.push(label("l" + i, (i * 37) % 800, (i * 91) % 600, { priority: i % 7 }));
  const placed = place(many, { budget: 25 });
  assert.ok(placed.length > 0 && placed.length <= 25, String(placed.length));
  for (let i = 0; i < placed.length; i++) for (let j = i + 1; j < placed.length; j++) assert.equal(overlap(placed[i], placed[j]), false);
  assert.equal(place(many).length <= M.DEFAULT_BUDGET, true);
});

test("the size range hides a label that is too small or too large, and a shown one keeps showing a little past the edge", () => {
  const nebula = (pixels) => label("neb", 400, 300, { kind: "nebula", pixels, minPixels: 40 });
  assert.equal(place([nebula(39)]).length, 0, "just under the range is not shown");
  assert.equal(place([nebula(41)]).length, 1);
  assert.equal(place([nebula(36)], { previous: new Set(["neb"]) }).length, 1, "hysteresis for a label already shown");
  assert.equal(place([nebula(30)], { previous: new Set(["neb"]) }).length, 0, "but not far past it");
  assert.equal(place([label("big", 400, 300, { pixels: 900, maxPixels: 600 })]).length, 0);
});

test("the same input gives the same labels in the same places, and a shown label wins an equal tie", () => {
  const input = [label("a", 400, 300), label("b", 402, 300)];
  const first = place(input);
  assert.deepEqual(place(input).map((p) => [p.candidate.id, p.left, p.top]), first.map((p) => [p.candidate.id, p.left, p.top]));
  assert.equal(first[0].candidate.id, "a");
  assert.equal(place(input, { previous: new Set(["b"]) })[0].candidate.id, "b", "stable between frames");
});

test("a label near an edge is kept inside the view", () => {
  const [near] = place([label("edge", 790, 5, { text: "A rather long name" })]);
  assert.ok(near.left >= 0 && near.left + near.width <= view.width && near.top >= 0);
});

test("the layer writes names as text, fades labels in, and takes dropped ones away", () => {
  const { document } = installDom("http://localhost/galaxy");
  const container = document.createElement("div");
  document.body.appendChild(container);
  const layer = M.createLabelLayer({
    container, project, measure, viewport: () => view, reducedMotion: () => true,
  });
  layer.update([label("<b>x</b>", 400, 300, { text: "<b>x</b>" })]);
  assert.equal(container.children.length, 1);
  assert.equal(container.children[0].textContent, "<b>x</b>", "markup in a name is only text");
  assert.equal(container.children[0].className, "map-label");
  assert.deepEqual(layer.shown().map((s) => s.text), ["<b>x</b>"]);
  layer.update([]);
  assert.equal(container.children.length, 0, "gone at once with reduced motion");
  assert.deepEqual(layer.shown(), []);
});

test("placeBoxes keeps a box that clears, shifts one that clashes, slides one back in and hides one that cannot clear", () => {
  const box = (x, y) => ({ x, y, w: 40, h: 10 });
  const first = M.placeBoxes([
    { id: "a", box: box(10, 10) },
    { id: "b", box: box(12, 12), shifts: [0, 14] },
    { id: "c", box: box(12, 12) },
    { id: "edge", box: box(-5, 60) },
  ], { bounds: { x: 0, y: 0, w: 100, h: 100 } });
  assert.deepEqual(first.map((p) => [p.id, p.hidden]), [["a", false], ["b", false], ["c", true], ["edge", false]]);
  assert.equal(first[0].dy, 0);
  assert.equal(first[1].dy, 14, "b moved down a line to clear a");
  assert.equal(first[3].dx, 7, "slid back in by the margin (2) and the overhang (5)");
});

test("placeBoxes keeps clear of marker circles but not of the label's own", () => {
  const items = [{ id: "n", box: { x: 10, y: 10, w: 30, h: 10 }, own: ["mine"] }];
  const mine = { id: "mine", x: 20, y: 15, r: 6 };
  const other = { id: "other", x: 20, y: 15, r: 6 };
  assert.equal(M.placeBoxes(items, { obstacles: [mine] })[0].hidden, false);
  assert.equal(M.placeBoxes(items, { obstacles: [other] })[0].hidden, true);
});
