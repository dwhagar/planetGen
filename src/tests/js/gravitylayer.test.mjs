// tests/js/gravitylayer.test.mjs -- static/gravitylayer.js, the Sector Map's gravity layer (MAP.168): the colors (percentile range, log scale, ramp, transparency), the
// pixels of a slab in the scene's frame, the words of the legend, and the controller's loading, retrying, drawing, URL and saddle picking.

import { test } from "node:test";
import assert from "node:assert/strict";

import { h, installDom, STATIC_URL } from "./fakedom.mjs";

installDom("http://localhost/sector/1");
const THREE = await import(new URL("vendor/three.module.min.js", STATIC_URL).href);
const G = await import(new URL("gravitylayer.js", STATIC_URL).href);

const ZONES = 16;

// A grid whose zone (i, j, k) reads `f(i, j, k)`.
function gridOf(f, extra) {
  const values = [];
  for (let i = 0; i < ZONES; i++) for (let j = 0; j < ZONES; j++) for (let k = 0; k < ZONES; k++) values.push(f(i, j, k));
  return Object.assign({
    mode: "pull", field: "centered", unit: "m/s^2", zones: ZONES, zone_edge_pc: 0.25, edge_pc: 4, center_pc: [100, 200, 300], origin_pc: [98, 198, 298],
    values: values, minimum: Math.min.apply(null, values), maximum: Math.max.apply(null, values), smooth: { level: 1.6e-10, span: 1e-13 }, saddles: [], stamp: "s",
  }, extra || {});
}

test("the ramp runs from dark to bright and is clamped", () => {
  assert.deepEqual(G.rampColor(0), G.RAMP[0]);
  assert.deepEqual(G.rampColor(1), G.RAMP[G.RAMP.length - 1]);
  assert.deepEqual(G.rampColor(-3), G.RAMP[0]);
  assert.deepEqual(G.rampColor(9), G.RAMP[G.RAMP.length - 1]);
  const mid = G.rampColor(0.5);
  assert.deepEqual(mid, G.RAMP[2]);
});

test("a percentile interpolates between the sorted values", () => {
  assert.equal(G.percentile([5, 1, 3, 2, 4], 0), 1);
  assert.equal(G.percentile([5, 1, 3, 2, 4], 100), 5);
  assert.equal(G.percentile([5, 1, 3, 2, 4], 50), 3);
  assert.equal(G.percentile([1, 2], 50), 1.5);
  assert.equal(G.percentile([], 50), 0);
});

test("the colors' range is the sector's own 2nd to 98th percentile, and a sector with no variation is flat", () => {
  const values = Array.from({ length: 1000 }, (_unused, n) => 1 + n / 1000);
  values[0] = 1000;                              // one zone at a star's core does not stretch the range
  const range = G.colorRange(values);
  assert.ok(!range.flat && range.hi < 2.1 && range.lo > 0.99 && range.lo < 1.1);
  assert.ok(G.colorRange(new Array(50).fill(3)).flat);
  assert.ok(G.colorRange(new Array(50).fill(0)).flat);
});

test("a shade runs from 0 at the bottom of the range to 1 at its top on a log scale, clipped outside it", () => {
  const range = { lo: 1, hi: 100, flat: false };
  assert.equal(G.shade(1, range), 0);
  assert.equal(G.shade(100, range), 1);
  assert.ok(Math.abs(G.shade(10, range) - 0.5) < 1e-9);
  assert.equal(G.shade(0.1, range), 0);
  assert.equal(G.shade(1e6, range), 1);
  assert.equal(G.shade(50, { lo: 1, hi: 100, flat: true }), 0);
  // a depth range that starts at 0 spans at most DECADES decades
  assert.ok(Math.abs(G.shade(1, { lo: 0, hi: 1e4, flat: false }) - 0) < 1e-9);
  assert.ok(Math.abs(G.shade(1e2, { lo: 0, hi: 1e4, flat: false }) - 0.5) < 1e-9);
});

test("no zone is filled past half opacity and the bottom of the range is clear", () => {
  assert.equal(G.alphaOf(0), 0);
  assert.ok(Math.abs(G.alphaOf(1) - G.MAX_ALPHA) < 1e-12 && G.MAX_ALPHA <= 0.5);
  assert.ok(G.alphaOf(0.5) > 0 && G.alphaOf(0.5) < G.alphaOf(1));
  // sixteen volume slabs one behind another stay under half opacity
  assert.ok(1 - Math.pow(1 - G.VOLUME_ALPHA, ZONES) <= 0.5);
});

test("a slab's pixels put the zone's i in the column and turn j over when the scene runs y the other way", () => {
  const values = gridOf((i, j) => (i === 3 && j === 5 ? 100 : 1)).values;
  const range = { lo: 1, hi: 100, flat: false };
  const plain = G.slabPixels(values, ZONES, 0, range, { flipY: false });
  const turned = G.slabPixels(values, ZONES, 0, range, { flipY: true });
  assert.equal(plain.length, ZONES * ZONES * 4);
  assert.equal(plain[(5 * ZONES + 3) * 4 + 3], Math.round(G.MAX_ALPHA * 255));          // the bright zone: row j, column i
  assert.equal(turned[((ZONES - 1 - 5) * ZONES + 3) * 4 + 3], Math.round(G.MAX_ALPHA * 255));
  assert.equal(plain[3], 0);                                                           // an ordinary zone is clear
  assert.deepEqual(Array.from(plain.slice((5 * ZONES + 3) * 4, (5 * ZONES + 3) * 4 + 3)), G.RAMP[G.RAMP.length - 1]);
});

test("a step along the slabs stays between the first and the last", () => {
  assert.equal(G.stepSlab(0, -1, ZONES), 0);
  assert.equal(G.stepSlab(15, 1, ZONES), 15);
  assert.equal(G.stepSlab(7, 1, ZONES), 8);
});

test("a value is written to three digits with its unit", () => {
  assert.equal(G.valueText("pull", 1.5534e-10), "1.55e-10 m/s²");
  assert.equal(G.valueText("tidal", 3.1e-27), "3.10e-27 s⁻²");
  assert.equal(G.valueText("depth", 1234.5), "1230 m²/s²");
  assert.equal(G.valueText("pull", 0), "0 m/s²");
});

test("a point of the grid goes to the scene's frame, parsecs from the center", () => {
  assert.deepEqual(G.localPoint([101, 198, 301], [100, 200, 300], false), [1, -2, 1]);
  assert.deepEqual(G.localPoint([101, 198, 301], [100, 200, 300], true), [1, 2, 1]);
});

test("the legend says what the colors are and that the galaxy's part is left out", () => {
  const grid = gridOf(() => 1);
  const pull = G.legendText({ mode: "pull", slab: 7, volume: false, grid: grid, range: { lo: 1.55e-10, hi: 1.66e-10, flat: false }, building: false, failed: false });
  assert.match(pull, /^Pull \(slab 8 of 16\): 1\.55e-10 m\/s² to 1\.66e-10 m\/s², the sector's own 2% to 98%, log scale\./);
  assert.match(pull, /Total pull 1\.60e-10 m\/s² from the smooth galaxy at the center/);
  const depth = G.legendText({ mode: "depth", slab: 0, volume: true, grid: grid, range: { lo: 1, hi: 1e4, flat: false }, building: false, failed: false });
  assert.match(depth, /all slabs/);
  assert.match(depth, /the galaxy's own slope across the sector \(1\.00e-13 m²\/s²\) is left out/);
  assert.match(G.legendText({ mode: "pull", slab: 0, volume: false, grid: grid, range: { lo: 1, hi: 1, flat: true }, building: false, failed: false }), /no variation across this sector/);
  assert.match(G.legendText({ mode: "pull", building: true, grid: null }), /working out this sector's grid/);
  assert.match(G.legendText({ mode: "pull", failed: true }), /could not be loaded/);
});

// --- The controller -----------------------------------------------------------------------------------------------------------------------------

function rig(responses) {
  const calls = [];
  const timers = [];
  const group = new THREE.Group();
  const parent = new THREE.Group();
  parent.add(group);
  const layers = [];
  const changes = [];
  const controller = G.createGravityController(THREE, {
    flipY: true,
    fetchGrid: (id, mode, field) => {
      calls.push([mode, field]);
      const next = responses.shift();
      return typeof next === "function" ? Promise.resolve().then(() => next(mode, field)) : Promise.resolve(next);
    },
    onChange: () => changes.push(1),
    setTimeout: (fn, ms) => { timers.push({ fn, ms }); return timers.length; },
    clearTimeout: () => {},
    picker: { addLayer: (layer) => { layers.push(layer); return layer; }, removeLayer: (layer) => { const at = layers.indexOf(layer); if (at >= 0) layers.splice(at, 1); } },
  });
  return { controller, calls, timers, group, layers, changes, flush: () => new Promise((resolve) => setImmediate(resolve)) };
}

const centered = (extra) => gridOf((i, j, k) => 1 + (i === 4 ? 9 : 0) + k * 0.01, extra);
const total = () => gridOf(() => 2, { field: "total" });

test("turning the layer on loads the grid for colors, then the total for numbers, and draws the slab", async () => {
  const r = rig([centered(), total()]);
  r.controller.attach({ id: 7, group: r.group });
  assert.deepEqual(r.calls, []);                                       // off: nothing asked
  r.controller.setOn(true);
  await r.flush();
  await r.flush();
  assert.deepEqual(r.calls, [["pull", "centered"], ["pull", "total"]]);
  const drawn = r.group.children[0];
  assert.equal(drawn.name, "gravity-layer");
  const planes = drawn.children.filter((c) => c.userData.slab !== undefined);
  assert.equal(planes.length, 1);
  assert.equal(planes[0].userData.slab, G.DEFAULT_SLAB);
  assert.ok(Math.abs(planes[0].position.z - ((G.DEFAULT_SLAB + 0.5) * 0.25 - 2)) < 1e-9);
  assert.equal(planes[0].material.opacity, 1);
  assert.ok(r.changes.length > 0);
});

test("the volume draws every slab, fainter, and stepping moves the one slab", async () => {
  const r = rig([centered(), total()]);
  r.controller.attach({ id: 7, group: r.group });
  r.controller.setOn(true);
  await r.flush();
  await r.flush();
  const slabs = () => r.group.children[0].children.filter((c) => c.userData.slab !== undefined).map((c) => c.userData.slab);
  r.controller.step(1);
  assert.deepEqual(slabs(), [G.DEFAULT_SLAB + 1]);
  r.controller.setVolume(true);
  assert.deepEqual(slabs(), Array.from({ length: ZONES }, (_unused, n) => n));
  r.controller.setVolume(false);
  r.controller.setSlab(0);
  r.controller.step(-1);
  assert.deepEqual(slabs(), [0]);
});

test("a grid still being worked out is asked for again until it is there", async () => {
  const r = rig([{ building: true }, { building: true }, centered(), total()]);
  r.controller.attach({ id: 7, group: r.group });
  r.controller.setOn(true);
  await r.flush();
  assert.equal(r.timers.length, 1);
  assert.equal(r.timers[0].ms, G.RETRY_MS);
  assert.match(r.controller.legend(true).text, /working out/);
  assert.equal(r.group.children.length, 0);
  r.timers[0].fn();
  await r.flush();
  assert.equal(r.timers.length, 2);
  r.timers[1].fn();
  await r.flush();
  await r.flush();
  assert.equal(r.group.children.length, 1);
  assert.equal(r.calls.length, 4);
});

test("a grid that never arrives, or a failed ask, says so and draws nothing", async () => {
  const r = rig([() => { throw new Error("HTTP 500"); }]);
  r.controller.attach({ id: 7, group: r.group });
  r.controller.setOn(true);
  await r.flush();
  await r.flush();
  assert.match(r.controller.legend(true).text, /could not be loaded/);
  assert.equal(r.group.children.length, 0);
});

test("a sector with no variation draws nothing and says so", async () => {
  const r = rig([gridOf(() => 5), total()]);
  r.controller.attach({ id: 7, group: r.group });
  r.controller.setOn(true);
  await r.flush();
  await r.flush();
  assert.equal(r.group.children.length, 0);
  assert.match(r.controller.legend(true).text, /no variation across this sector/);
});

test("switching mode loads that mode's grid and leaving the sector removes the layer", async () => {
  const r = rig([centered(), total(), centered({ mode: "tidal" }), total()]);
  r.controller.attach({ id: 7, group: r.group });
  r.controller.setOn(true);
  await r.flush();
  await r.flush();
  r.controller.setMode("tidal");
  await r.flush();
  await r.flush();
  assert.deepEqual(r.calls.slice(2), [["tidal", "centered"], ["tidal", "total"]]);
  assert.equal(r.group.children.length, 1);
  r.controller.setMode("sound");                                        // not a mode: nothing changes
  assert.equal(r.controller.state.mode, "tidal");
  r.controller.detach();
  assert.equal(r.group.children.length, 0);
});

test("the saddle points are marked, put in the scene's frame and picked with their total value", async () => {
  const spot = { zone: [7, 8, 8], position_pc: [100.1, 200.2, 300.3] };
  const r = rig([centered({ saddles: [spot] }), gridOf(() => 3.5e-10, { field: "total" })]);
  r.controller.attach({ id: 7, group: r.group });
  r.controller.setOn(true);
  await r.flush();
  await r.flush();
  assert.equal(r.layers.length, 1);
  const marks = r.layers[0].meshes();
  assert.equal(marks.length, 1);
  assert.ok(Math.abs(marks[0].position.x - 0.1) < 1e-9 && Math.abs(marks[0].position.y + 0.2) < 1e-9 && Math.abs(marks[0].position.z - 0.3) < 1e-9);
  assert.equal(r.layers[0].entryOf({ object: marks[0] }), spot);
  assert.equal(r.layers[0].tooltip(spot), "Saddle point (zone 8, 9, 9): total 3.50e-10 m/s²");
  r.controller.setOn(false);
  assert.equal(r.layers.length, 0);
  assert.equal(r.group.children.length, 0);
});

test("the URL keeps the mode and the slab while the layer is on", () => {
  const r = rig([]);
  const off = new URLSearchParams("gravity=tidal&gslab=3&x=1");
  r.controller.writeUrl(off);
  assert.equal(off.get("gravity"), null);
  assert.equal(off.get("gslab"), null);
  assert.equal(off.get("x"), "1");
  r.controller.readUrl(new URLSearchParams("gravity=tidal&gslab=3"));
  assert.equal(r.controller.state.mode, "tidal");
  assert.equal(r.controller.state.slab, 2);
  r.controller.setOn(true);
  const on = new URLSearchParams();
  r.controller.writeUrl(on);
  assert.equal(on.get("gravity"), "tidal");
  assert.equal(on.get("gslab"), "3");
  r.controller.readUrl(new URLSearchParams("gravity=nonsense&gslab=all"));
  assert.equal(r.controller.state.mode, "tidal");
  assert.equal(r.controller.state.volume, true);
  const all = new URLSearchParams();
  r.controller.writeUrl(all);
  assert.equal(all.get("gslab"), "all");
});

test("the menu controls follow the state and drive the controller", async () => {
  const r = rig([centered(), total()]);
  r.controller.attach({ id: 7, group: r.group });
  const container = h("div", {});
  r.controller.settings(container);
  const modes = container.querySelectorAll("[data-gravity-mode]");
  assert.equal(modes.length, 3);
  assert.equal(modes.find((b) => b.dataset.gravityMode === "pull").getAttribute("aria-pressed"), "true");
  assert.equal(modes.find((b) => b.dataset.gravityMode === "tidal").getAttribute("aria-pressed"), "false");
  assert.match(container.querySelector(".gravity-slab").textContent, /slab 9 of 16/);
  container.querySelector("[data-gravity=next]").dispatchEvent({ type: "click" });
  assert.match(container.querySelector(".gravity-slab").textContent, /slab 10 of 16/);
  modes.find((b) => b.dataset.gravityMode === "depth").dispatchEvent({ type: "click" });
  assert.equal(r.controller.state.mode, "depth");
  assert.equal(modes.find((b) => b.dataset.gravityMode === "depth").getAttribute("aria-pressed"), "true");
  const box = container.querySelector("[data-gravity=volume]");
  box.checked = true;
  box.dispatchEvent({ type: "change" });
  assert.equal(r.controller.state.volume, true);
  assert.match(container.querySelector(".gravity-slab").textContent, /all slabs/);
  assert.equal(container.querySelector("[data-gravity=next]").disabled, true);
});
