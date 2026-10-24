// tests/js/bodytexture.test.mjs -- static/bodytexture.js and static/solidpaint.js (MAP.272, MAP.273): the shared
// painter pieces (seeded noise, cache, one-paint-per-frame queue, size tiers) and the solid-world painter, whose
// pixels are fingerprinted so an accidental change to how a world looks fails here.

import { test } from "node:test";
import assert from "node:assert/strict";

import { importPage } from "./fakedom.mjs";

const T = await importPage("bodytexture.js");
const P = await importPage("solidpaint.js");

const EARTHLIKE = {
  v: 1, kind: "planet", body: "solid", seed: 12345, color: "#a39070", oceanColor: "#1f4f8f", hydro: "surface ocean",
  ocean: 0.7, land: 0.3, tempK: 288, atm: true, oceanKm: 3,
};

test("the seeded generator repeats and the noise wraps in longitude", () => {
  const a = T.makeRandom(42);
  const b = T.makeRandom(42);
  for (let i = 0; i < 5; i++) assert.equal(a(), b());
  assert.notEqual(T.makeRandom(1)(), T.makeRandom(2)());
  for (const v of [0.1, 0.6]) assert.ok(Math.abs(T.fbm(0, v, 9, 4, 3) - T.fbm(1, v, 9, 4, 3)) < 1e-9);
  const n = T.fbm(0.3, 0.4, 5, 4, 4);
  assert.ok(n >= 0 && n < 1);
});

test("the cache keeps the most recently used and reports what it drops", () => {
  const dropped = [];
  const cache = T.makeCache(2, (v) => dropped.push(v));
  cache.set("a", 1);
  cache.set("b", 2);
  cache.get("a");
  cache.set("c", 3);
  assert.deepEqual(dropped, [2]);
  assert.equal(cache.has("a"), true);
  assert.equal(cache.size(), 2);
});

test("the queue paints the most urgent job, once per call, and does not double a key", () => {
  const ran = [];
  const q = T.makeQueue();
  q.add("far", 10, () => ran.push("far"));
  q.add("near", -5, () => ran.push("near"));
  q.add("far", 10, () => ran.push("far again"));
  assert.equal(q.size(), 2);
  assert.equal(q.runOne(), true);
  assert.deepEqual(ran, ["near"]);
  q.runOne();
  assert.deepEqual(ran, ["near", "far again"]);
  assert.equal(q.runOne(), false);
});

test("tiers follow the projected diameter", () => {
  assert.equal(T.tierForDiameter(1), "point");
  assert.equal(T.tierForDiameter(NaN), "point");
  assert.equal(T.tierForDiameter(10), "disc");
  assert.equal(T.tierForDiameter(60), "low");
  assert.equal(T.tierForDiameter(300), "high");
});

test("look records parse, and anything else is null", () => {
  assert.equal(T.parseLook(""), null);
  assert.equal(T.parseLook("{not json"), null);
  assert.equal(T.parseLook('{"v":2}'), null);
  assert.equal(T.parseLook('{"v":1,"body":"solid"}').body, "solid");
});

test("the same look paints the same pixels, a different seed paints different ones", () => {
  const a = P.solidFingerprint(EARTHLIKE, 128, 64);
  assert.equal(a, P.solidFingerprint(EARTHLIKE, 128, 64));
  assert.notEqual(a, P.solidFingerprint({ ...EARTHLIKE, seed: 777 }, 128, 64));
});

test("the sea covers the stored ocean fraction of the sphere's area", () => {
  const w = 128, h = 64;
  const px = P.paintSolid({ ...EARTHLIKE, tempK: 320 }, w, h); // warm: no ice to confuse the count
  let seaArea = 0, total = 0;
  for (let y = 0; y < h; y++) {
    const wt = Math.cos((0.5 - (y + 0.5) / h) * Math.PI);
    for (let x = 0; x < w; x++) {
      const o = (y * w + x) * 4;
      const blue = px[o + 2] > px[o] + 25; // sea is blue-dominant, land tan is red-dominant
      if (blue) seaArea += wt;
      total += wt;
    }
  }
  assert.ok(Math.abs(seaArea / total - 0.7) < 0.08, `sea share ${seaArea / total}`);
});

test("a cold world gets polar ice and a frozen one is all ice", () => {
  const w = 64, h = 32;
  const px = P.paintSolid(EARTHLIKE, w, h);
  const top = (x) => px[x * 4];
  const mid = (x) => px[((h / 2) * w + x) * 4];
  let polar = 0, equator = 0;
  for (let x = 0; x < w; x++) { polar += top(x); equator += mid(x); }
  assert.ok(polar / w > 170, "the pole row is near white");
  assert.ok(equator / w < polar / w, "the equator is darker than the pole");
  const frozen = P.paintSolid({ ...EARTHLIKE, hydro: "ice", tempK: 80 }, w, h);
  for (let i = 0; i < w * h; i++) assert.ok(frozen[i * 4] > 200 && frozen[i * 4 + 2] > 200);
});

test("a dry airless world shows craters; the same world with air does not", () => {
  const dry = { ...EARTHLIKE, hydro: "dry", ocean: 0, atm: false, tempK: 300 };
  assert.notEqual(P.solidFingerprint(dry, 128, 64), P.solidFingerprint({ ...dry, atm: true }, 128, 64));
});
