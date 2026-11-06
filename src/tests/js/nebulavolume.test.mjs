// tests/js/nebulavolume.test.mjs -- static/nebulavolume.js (MAP.365): the packed points of a node's nebulae, which nodes ask for a
// volume, the capped sum of overlapping colors, and the store that fetches them.

import { test } from "node:test";
import assert from "node:assert/strict";

import { installDom, importPage } from "./fakedom.mjs";

installDom();
const V = await importPage("nebulavolume.js");

// Packs the way planetgen.galaxy.nebula_volume.pack does.
function pack(points, spacing, level) {
  const n = points.length;
  const bytes = new Uint8Array(14 + 14 * n);
  const view = new DataView(bytes.buffer);
  bytes.set([78, 86, 79, 76, 1, level]);
  view.setUint32(6, n, true);
  view.setFloat32(10, spacing, true);
  points.forEach((p, i) => {
    view.setInt32(14 + 4 * i, p[0], true);
    view.setInt32(14 + 4 * n + 4 * i, p[1], true);
    view.setInt32(14 + 8 * n + 4 * i, p[2], true);
    bytes[14 + 12 * n + i] = p[3];
    bytes[14 + 13 * n + i] = p[4];
  });
  return bytes;
}

test("a packed volume decodes to its columns and its positions", () => {
  const volume = V.decodeVolume(pack([[1, -2, 3, 200, 4], [-5, 6, 0, 17, 1]], 2.5, 3));
  assert.equal(volume.count, 2);
  assert.equal(volume.level, 3);
  assert.equal(volume.spacingPc, 2.5);
  assert.deepEqual(Array.from(volume.ix), [1, -5]);
  assert.deepEqual(Array.from(volume.iz), [3, 0]);
  assert.deepEqual(Array.from(volume.strength), [200, 17]);
  assert.deepEqual(Array.from(volume.family), [4, 1]);
  assert.deepEqual(Array.from(V.positionsOf(volume, [1, 1, 1])), [1.5, -6, 6.5, -13.5, 14, -1]);
  assert.deepEqual(V.centerOf(volume).map((v) => Math.round(v * 100) / 100), [-5, 5, 3.75]);
});

test("bytes that are not a volume are refused", () => {
  assert.throws(() => V.decodeVolume(new Uint8Array(20)), /not a nebula volume/);
  const good = pack([[0, 0, 0, 5, 0]], 1, 0);
  assert.throws(() => V.decodeVolume(good.slice(0, good.length - 1)), /wrong length/);
  const wrong = pack([], 1, 0);
  wrong[4] = 9;
  assert.throws(() => V.decodeVolume(wrong), /version 9/);
});

test("overlapping colors add and bend toward the cap, keeping the hue and never passing it", () => {
  const pink = [0.9, 0.3, 0.4];
  const blue = [0.2, 0.4, 0.9];
  assert.deepEqual(V.cappedSum([], V.CAP), [0, 0, 0]);
  const faint = V.cappedSum([pink.map((v) => v * 0.01)], V.CAP);
  assert.ok(Math.abs(faint[0] - 0.009) < 0.0002, "a faint sum is the sum");
  const two = V.cappedSum([pink.map((v) => v * 0.05), blue.map((v) => v * 0.05)], V.CAP);
  const plain = [0.055, 0.035, 0.065];
  two.forEach((v, i) => assert.ok(Math.abs(v - plain[i]) < 0.005, "two faint nebulae add"));
  for (const crowd of [5, 50, 5000]) {
    const terms = [];
    for (let i = 0; i < crowd; i++) terms.push(i % 2 ? pink : blue);
    const sum = V.cappedSum(terms, V.CAP);
    assert.ok(Math.max(...sum) <= V.CAP, "a crowd of " + crowd + " never passes the cap");
    assert.ok(Math.max(...sum) > 0.9 * V.CAP * (1 - Math.exp(-crowd * 0.01)), "and gets near it");
    const ratio = sum.map((v, i) => v / (blue[i] * Math.ceil(crowd / 2) + pink[i] * Math.floor(crowd / 2)));
    assert.ok(Math.abs(ratio[0] - ratio[1]) < 1e-9 && Math.abs(ratio[1] - ratio[2]) < 1e-9, "every channel scaled alike: the hue stays");
  }
  const white = V.cappedSum(Array.from({ length: 1000 }, () => [1, 1, 1]), V.CAP);
  assert.ok(Math.max(...white) <= V.CAP);
  assert.ok(V.CAP < 0.6, "the cap leaves the layer well short of white");
});

test("the soft cap is monotone and bounded", () => {
  let last = -1;
  for (let m = 0; m < 20; m += 0.25) {
    const bent = V.softCap(m, 0.5);
    assert.ok(bent >= last && bent <= 0.5);
    last = bent;
  }
  assert.ok(Math.abs(V.softCap(0.001, 0.5) - 0.001) < 1e-5);
});

test("a splat is as many pixels as its steps are, between its limits", () => {
  assert.equal(V.splatPx(1, 1e9, 800), V.MIN_SPLAT_PX);
  assert.equal(V.splatPx(100, 0.001, 800), V.MAX_SPLAT_PX);
  assert.ok(Math.abs(V.splatPx(2, 100, 1000) - V.SPLAT_STEPS * 20) < 1e-9);
});

test("a node is asked for only once a nebula in it is wide enough, the biggest first", () => {
  const near = V.sourceOf({ key: "27.3.0.5", block: {}, node: { clouds: [{ x: 0, y: 0, z: 0, radius_pc: 30 }] } });
  const far = V.sourceOf({ key: "27.4.0.1", block: {}, node: { nebula_regions: [{ center: [5000, 0, 0], radius_pc: 30 }] } });
  const empty = V.sourceOf({ key: "27.5.0.2", block: {}, node: { clouds: [] } });
  assert.equal(empty, null);
  assert.equal(V.sourceOf(null), null);
  const camera = [0, 0, 600];
  const found = V.wanted([far, near], camera, 800);
  assert.deepEqual(found.map((f) => f.source.key), ["27.3.0.5"]);
  assert.ok(found[0].radiusPx >= V.MIN_RADIUS_PX);
  assert.equal(V.wanted([far, near], camera, 800, () => false).length, 0);
  const many = Array.from({ length: 40 }, (_, i) => V.sourceOf({ key: "k" + i, block: {}, node: { clouds: [{ x: i, y: 0, z: 0, radius_pc: 30 }] } }));
  assert.equal(V.wanted(many, camera, 800).length, V.MAX_LOADED);
});

test("a volume shows from a pixel width and only for the share the cut gives its node", () => {
  assert.equal(V.volumeAlpha(V.MIN_RADIUS_PX, 1), 0);
  assert.equal(V.volumeAlpha(V.FULL_RADIUS_PX * 4, 1), 1);
  assert.equal(V.volumeAlpha(V.FULL_RADIUS_PX * 4, 0.25), 0.25);
  assert.equal(V.volumeAlpha(V.FULL_RADIUS_PX * 4, 0), 0);
  assert.ok(V.volumeAlpha(21, 1) > 0 && V.volumeAlpha(21, 1) < 1);
});

test("the store asks once for a key, three at a time, and tries a failed one again later", async () => {
  let time = 0;
  const asked = [];
  const waiting = [];
  const changes = [];
  const store = V.createVolumeStore((key) => {
    asked.push(key);
    return new Promise((resolve, reject) => waiting.push({ key, resolve, reject }));
  }, (key) => changes.push(key), () => time);
  ["a", "b", "c", "d", "a"].forEach((key) => store.request(key));
  await Promise.resolve();
  await Promise.resolve();
  assert.deepEqual(asked, ["a", "b", "c"]);
  waiting[0].resolve(pack([[0, 0, 0, 9, 0]], 1, 0));
  await new Promise((resolve) => setTimeout(resolve, 5));
  assert.equal(store.get("a").status, "ready");
  assert.equal(store.get("a").volume.count, 1);
  assert.deepEqual(asked, ["a", "b", "c", "d"]);
  waiting[1].reject(new Error("503"));
  await new Promise((resolve) => setTimeout(resolve, 5));
  assert.equal(store.get("b").status, "failed");
  store.request("b");
  assert.equal(store.get("b").status, "failed", "not at once");
  time += V.RETRY_MS;
  assert.equal(store.request("b").status !== "failed", true);
  store.release(new Set(["a"]));
  assert.equal(store.get("a").status, "ready");
  store.release(new Set());
  assert.equal(store.get("a"), null, "a key not kept is let go");
  assert.equal(store.get("c").status, "loading", "but not one on its way");
});

test("a star behind a nebula takes its hue, more the denser it is, and a little of its light is lost", () => {
  const white = [1, 1, 1];
  assert.deepEqual(V.tinted(white, [0, 0, 0]), white, "no nebula, no change");
  const faint = V.tinted(white, [0.02, 0.005, 0.01]);
  const dense = V.tinted(white, [0.8, 0.2, 0.4]);
  assert.ok(faint[0] > faint[1] && faint[0] > faint[2], "the star leans toward the nebula's hue");
  assert.ok(dense[0] - dense[1] > faint[0] - faint[1], "a denser core tints more");
  assert.ok(V.tintOf([5, 1, 2]).strength <= V.TINT_MAX && V.tintOf([5, 1, 2]).strength > 0.99 * V.TINT_MAX, "and the tint levels off, so it stays subtle");
  assert.ok(Math.max(...dense) < 1 && dense[0] > 0.7 * 1, "some light is lost but the star is still there");
  const blue = V.tinted([0.4, 0.5, 1.0], [0.1, 0.1, 0.5]);
  assert.ok(blue[2] > blue[0], "a blue star behind a blue nebula stays blue");
  assert.ok(V.NEBULA_TINT_GLSL.includes("nebulaTinted") && V.NEBULA_TINT_GLSL.includes(V.TINT_MAX.toFixed(3)), "the shader carries the same constants");
});

test("a cloud given by its sphere is its own key, the same for the same sphere, and the key gives the sphere back (MAP.392)", () => {
  const sphere = { x: 1.23456, y: -2, z: 3.5, radiusPc: 4.25 };
  const source = V.sphereSource("supernova-remnant", sphere, 0.5);
  assert.equal(source.key, "sphere:supernova-remnant:1.235,-2.000,3.500,4.250");
  assert.equal(source.single, true);
  assert.equal(source.block, null);
  assert.equal(source.weight, 0.5);
  assert.deepEqual(source.reaches, [sphere]);
  assert.equal(V.sphereSource("supernova-remnant", { ...sphere }, 1).key, source.key, "one request for one cloud, however often the scene is built");
  assert.deepEqual(V.sphereOf(source.key), { family: "supernova-remnant", x: 1.235, y: -2, z: 3.5, radiusPc: 4.25 });
  assert.equal(V.sphereOf("27.3.0.5"), null, "a node key is not one");
  assert.equal(V.sphereOf("nebula:5"), null, "nor is a stored nebula's");
  assert.equal(V.singleIdOf(source.key), null, "and a sphere key is not a stored nebula's id");
  assert.equal(V.sphereSource("diffuse", sphere, 9).weight, 1, "held to 0..1");
  assert.equal(V.wanted([source], [0, 0, 0], 800).length, 1, "it is asked for like any source once wide enough");
});

test("a single nebula is its own key, has no node block, and shows at the weight it is given (MAP.142)", () => {
  const source = V.singleSource("00001000000000-4000000-000", { x: 1, y: 2, z: 3, radiusPc: 4 }, 0.5);
  assert.deepEqual(source, { key: "nebula:00001000000000-4000000-000", single: true, block: null, weight: 0.5, reaches: [{ x: 1, y: 2, z: 3, radiusPc: 4 }] });
  assert.equal(V.singleIdOf(source.key), "00001000000000-4000000-000");
  assert.equal(V.singleIdOf("27.3.0.5"), null, "a node key is not one");
  assert.equal(V.singleSource(5, { x: 0, y: 0, z: 0, radiusPc: 1 }, 7).weight, 1, "held to 0..1");
  assert.equal(V.singleSource(5, { x: 0, y: 0, z: 0, radiusPc: 1 }, -2).weight, 0);
  assert.equal(V.singleSource(5, { x: 0, y: 0, z: 0, radiusPc: 1 }).weight, 1);
  const found = V.wanted([source], [0, 0, 0], 800);
  assert.equal(found.length, 1, "it is asked for like any source once wide enough");
});
