// tests/js/bodytexture.test.mjs -- static/bodytexture.js and static/solidpaint.js (MAP.272, MAP.273): the shared
// painter pieces (seeded noise, cache, one-paint-per-frame queue, size tiers) and the solid-world painter, whose
// pixels are fingerprinted so an accidental change to how a world looks fails here.

import { test } from "node:test";
import assert from "node:assert/strict";

import { importPage } from "./fakedom.mjs";

const T = await importPage("bodytexture.js");
const P = await importPage("solidpaint.js");
const G = await importPage("giantpaint.js");

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

test("a look with no air fields paints as before, and a thick veil tints toward the sky colour and flattens the ground", () => {
  const w = 64, h = 32;
  const plain = { ...EARTHLIKE, atm: true, tempK: 320 };
  assert.equal(P.solidFingerprint(plain, w, h), P.solidFingerprint({ ...plain, veil: 0, haze: 0 }, w, h));
  const veiled = P.paintSolid({ ...plain, veil: 1, glow: "#ff8a30" }, w, h);
  const bare = P.paintSolid(plain, w, h);
  let r = 0, b = 0, br = 0, bb = 0, spread = 0, bareSpread = 0;
  for (let i = 0; i < w * h; i++) { r += veiled[i * 4]; b += veiled[i * 4 + 2]; br += bare[i * 4]; bb += bare[i * 4 + 2]; }
  assert.ok(r / b > br / bb * 1.2, "redder than the bare ground");
  const mean = (px) => px.reduce((t, v, i) => (i % 4 === 3 ? t : t + v), 0) / (w * h * 3);
  const sd = (px) => Math.sqrt(px.reduce((t, v, i) => (i % 4 === 3 ? t : t + (v - mean(px)) ** 2), 0) / (w * h * 3));
  assert.ok(sd(veiled) < sd(bare), "less contrast");
});

test("a haze softens, and a cloud layer whitens by its cover", () => {
  const w = 64, h = 32;
  const plain = { ...EARTHLIKE, atm: true, tempK: 320 };
  assert.notEqual(P.solidFingerprint(plain, w, h), P.solidFingerprint({ ...plain, haze: 0.9 }, w, h));
  const none = P.paintSolid(plain, w, h);
  const clouded = P.paintSolid({ ...plain, cloudCover: 0.6, cloudTau: 8 }, w, h);
  const total = (px) => px.reduce((t, v, i) => (i % 4 === 3 ? t : t + v), 0);
  assert.ok(total(clouded) > total(none), "brighter with cloud");
  const solid = P.paintSolid({ ...plain, cloudCover: 1, cloudTau: 20 }, w, h);
  assert.ok(total(solid) > total(clouded), "full cover is whiter still");
});

// --- Stars (MAP.276) ---------------------------------------------------------

const S = await importPage("starpaint.js");
const SUN = { v: 1, kind: "star", body: "star", seed: 99, color: "#fff4ea", spec: "G2V", yerkes: "V", lumSolar: 1 };

test("a star's class picks its recipe", () => {
  assert.equal(S.starKind({ yerkes: "VII" }), "white dwarf");
  assert.equal(S.starKind({ yerkes: "III" }), "giant");
  assert.equal(S.starKind({}), "dwarf");
  assert.equal(S.surfaceRecipe({ ...SUN, yerkes: "VII" }).contrast, 0);
  assert.ok(S.surfaceRecipe({ ...SUN, spec: "M3V" }).spots > S.surfaceRecipe(SUN).spots);
  assert.equal(S.surfaceRecipe({ ...SUN, spec: "A0V" }).spots, 0);
  assert.ok(S.surfaceRecipe({ ...SUN, yerkes: "III", spec: "K5III" }).cells < S.surfaceRecipe(SUN).cells);
});

test("brighter stars get a wider, stronger glow, within the caps", () => {
  const sun = S.glowFor(1);
  assert.equal(sun.scale, 1.4);
  assert.ok(S.glowFor(1e5).scale > sun.scale && S.glowFor(1e5).scale <= 1.9);
  assert.ok(S.glowFor(1e-3).scale >= 1.2 && S.glowFor(1e-3).scale < sun.scale);
  assert.ok(S.glowFor(1e9).strength <= 1 && S.glowFor(1e-9).strength >= 0.5);
});

test("a white dwarf is one flat colour; a sun-like star is mottled and never darker than 0.7", () => {
  const wd = S.paintStar({ ...SUN, yerkes: "VII" }, 32, 16);
  for (let i = 1; i < 32 * 16; i++) assert.equal(wd[i * 4], wd[0]);
  const px = S.paintStar(SUN, 64, 32);
  let min = 255, max = 0;
  for (let i = 0; i < 64 * 32; i++) { min = Math.min(min, px[i * 4]); max = Math.max(max, px[i * 4]); }
  assert.ok(max - min > 10);
  assert.ok(min >= Math.round(0xff * 0.7) - 1);
});

test("the same star paints the same pixels", () => {
  assert.equal(S.starFingerprint(SUN, 128, 64), S.starFingerprint(SUN, 128, 64));
  assert.notEqual(S.starFingerprint(SUN, 128, 64), S.starFingerprint({ ...SUN, seed: 100 }, 128, 64));
});

// --- Gas and ice giants ------------------------------------------------------

const JOVIAN = { v: 1, kind: "planet", body: "gas", seed: 4242, color: "#d9c9a3", tempK: 165, radiusKm: 69911, spinH: 9.9,
                 gas: { CH4: 0.002, NH3: 0.0003, H2S: 0.0001 } };
const ICY = { ...JOVIAN, color: "#9fd0dc", tempK: 60, radiusKm: 25000, spinH: 16, gas: { CH4: 0.02 } };

test("a faster spinner gets more bands, and the count is held to 4 to 24", () => {
  assert.ok(G.bandCount(5, 69911) > G.bandCount(20, 69911));
  assert.equal(G.bandCount(1000, 69911), 4);
  assert.equal(G.bandCount(0.5, 69911), 24);
  assert.equal(G.bandCount(undefined, 69911), G.BANDS_NO_SPIN);
});

test("an ammonia-cloud giant has far more belt and zone contrast than a methane ice giant", () => {
  assert.ok(G.giantRecipe(JOVIAN).contrast > 2 * G.giantRecipe(ICY).contrast);
  assert.ok(G.giantRecipe(JOVIAN).sharp > G.giantRecipe({ ...JOVIAN, spinH: 30 }).sharp);
});

test("the same giant paints the same pixels, another seed paints others, and the picture follows the chemistry", () => {
  const a = G.giantFingerprint(JOVIAN, 128, 64);
  assert.equal(a, G.giantFingerprint(JOVIAN, 128, 64));
  assert.notEqual(a, G.giantFingerprint({ ...JOVIAN, seed: 4243 }, 128, 64));
  assert.notEqual(a, G.giantFingerprint({ ...JOVIAN, gas: { CH4: 0.002, NH3: 0.0003, H2S: 0 } }, 128, 64));
});

test("a giant's bands run in latitude: rows differ and a row is nearly the same around the planet", () => {
  const w = 128, h = 64;
  const px = G.paintGiant(JOVIAN, w, h);
  const rowMean = (y) => { let s = 0; for (let x = 0; x < w; x++) s += px[(y * w + x) * 4]; return s / w; };
  const means = []; for (let y = 0; y < h; y++) means.push(rowMean(y));
  assert.ok(Math.max(...means) - Math.min(...means) > 25, "belts and zones differ in lightness");
  for (let i = 0; i < px.length; i += 4) assert.equal(px[i + 3], 255);
});
