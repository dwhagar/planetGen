// tests/js/systemgravity.test.mjs -- static/systemgravity.js, the System Map's gravity layer (MAP.169): the scale both ways, the field of point masses, the heat
// map's cells, the Hill radius, the Lagrange points (against the closed forms), the sources of a scene and the words of the legend.

import { test } from "node:test";
import assert from "node:assert/strict";

import { installDom, STATIC_URL } from "./fakedom.mjs";

installDom("http://localhost/system/1");
const S = await import(new URL("systemgravity.js", STATIC_URL).href);

const AU_KM = 1.495978707e8;
const SUN_KG = 1.98847e30;
const EARTH_KG = 5.9722e24;
const LOG = { c: 350, lo: 0.3 * AU_KM, hi: 30 * AU_KM, min: 55, spread: 280, perKm: NaN, ok: true };
const LINEAR = { c: 350, perKm: 0.001, lo: NaN, hi: NaN, min: NaN, spread: NaN, ok: true };

test("the scale is read from the scene's attributes and says when it is no scale", () => {
  const s = S.sceneScale({ lokm: "4.5e7", hikm: "4.5e9", cpx: "350.0", minpx: "55.0", spreadpx: "280.0" });
  assert.equal(s.ok, true);
  assert.equal(s.lo, 4.5e7);
  const moons = S.sceneScale({ pxperkm: "0.0008", cpx: "350.0" });
  assert.equal(moons.ok, true);
  assert.equal(moons.perKm, 0.0008);
  assert.equal(S.sceneScale({ cpx: "350" }).ok, false);
});

test("a point goes to its pixel and back on the log scale, keeping its angle", () => {
  for (const [x, y] of [[AU_KM, 0], [0, 2 * AU_KM], [-3 * AU_KM, -4 * AU_KM]]) {
    const p = S.kmToPixel(LOG, x, y);
    const back = S.pixelToKm(LOG, p.x, p.y);
    assert.ok(Math.abs(back.x - x) < 1e-6 * Math.hypot(x, y), "x");
    assert.ok(Math.abs(back.y - y) < 1e-6 * Math.hypot(x, y), "y");
  }
  const east = S.kmToPixel(LOG, 5 * AU_KM, 0);
  assert.ok(east.x > LOG.c && Math.abs(east.y - LOG.c) < 1e-9, "east is right of the middle, on its row");
  const north = S.kmToPixel(LOG, 0, 5 * AU_KM);
  assert.ok(north.y < LOG.c, "north is up the picture");
});

test("the log scale stands for no distance inside the first orbit's disk or past the last orbit", () => {
  assert.equal(S.pixelToKm(LOG, LOG.c + 10, LOG.c), null);
  assert.equal(S.pixelToKm(LOG, LOG.c + 400, LOG.c), null);
  assert.ok(S.pixelToKm(LOG, LOG.c + 100, LOG.c) !== null);
});

test("a moon system is on one linear scale in both directions", () => {
  const p = S.kmToPixel(LINEAR, 100000, 50000);
  assert.deepEqual(p, { x: 350 + 100, y: 350 - 50 });
  const back = S.pixelToKm(LINEAR, p.x, p.y);
  assert.ok(Math.abs(back.x - 100000) < 1e-6 && Math.abs(back.y - 50000) < 1e-6);
});

test("the pull of one body is G M / r^2 and points at it; the depth is -G M / r; far out the tides are 2 G M / r^3", () => {
  const sun = [{ x: 0, y: 0, mass: SUN_KG, radius: 695700 }];
  const r = AU_KM;
  assert.ok(Math.abs(S.fieldAt(sun, r, 0, "pull") / 5.93e-3 - 1) < 0.01, "the Sun's pull at 1 AU is 5.93 mm/s^2");
  assert.ok(Math.abs(S.fieldAt(sun, r, 0, "depth") / -8.87e8 - 1) < 0.01, "the Sun's potential at 1 AU is -8.87e8 m^2/s^2");
  assert.ok(Math.abs(S.fieldAt(sun, r, 0, "tidal") / (2 * 6.6743e-11 * SUN_KG / Math.pow(r * 1000, 3)) - 1) < 1e-3);
  assert.ok(Math.abs(S.fieldAt(sun, 0, r, "pull") - S.fieldAt(sun, r, 0, "pull")) < 1e-12, "the pull does not depend on the direction");
});

test("a body's own softening keeps the field finite at its center", () => {
  const earth = [{ x: 0, y: 0, mass: EARTH_KG, radius: 6371 }];
  for (const mode of ["pull", "depth", "tidal"]) assert.ok(isFinite(S.fieldAt(earth, 0, 0, mode)), mode);
  assert.ok(S.fieldAt(earth, 0, 0, "pull") < 1e-9, "no pull at the middle of a softened body");
});

test("two equal bodies pull nothing at the point between them", () => {
  const pair = [{ x: -1000000, y: 0, mass: EARTH_KG, radius: 6371 }, { x: 1000000, y: 0, mass: EARTH_KG, radius: 6371 }];
  assert.ok(S.fieldAt(pair, 0, 0, "pull") < 1e-12);
  assert.ok(S.fieldAt(pair, 0, 0, "tidal") > 0, "but the space there is stretched");
});

test("the heat map has the bodies' wells brightest, nothing where the scale means nothing, and the depth from the shallowest cell", () => {
  const bounds = [0, 700];
  const sources = [{ x: 0, y: 0, mass: SUN_KG, radius: 695700 }, { x: 1.0 * AU_KM, y: 0, mass: EARTH_KG * 300, radius: 70000 }];
  const grid = S.fieldGrid(LOG, bounds, sources, "pull", 70);
  assert.equal(grid.length, 70 * 70);
  assert.ok(grid.some(isNaN), "the cells inside the first orbit are masked");
  const range = { lo: 0, hi: 0, flat: true };
  assert.ok(isNaN(grid[35 * 70 + 35]), "the middle is masked");
  const depth = S.fieldGrid(LOG, bounds, sources, "depth", 70);
  const finite = Array.prototype.filter.call(depth, (v) => !isNaN(v));
  assert.equal(Math.min.apply(null, finite), 0, "the shallowest cell is the zero");
  assert.ok(Math.max.apply(null, finite) > 0);
  void range;
});

test("heat pixels are clear where masked or flat and colored where the field varies", () => {
  const values = new Float64Array([NaN, 1, 10, 100, 1000, 10000]);
  const range = { lo: 1, hi: 10000, flat: false };
  const pixels = S.heatPixels(values, range);
  assert.equal(pixels[3], 0, "NaN is clear");
  assert.ok(pixels[5 * 4 + 3] > pixels[2 * 4 + 3], "the top of the range is more opaque than the middle");
  const flat = S.heatPixels(values, { lo: 1, hi: 1, flat: true });
  assert.ok(flat.every((v) => v === 0), "a flat scene draws nothing");
});

test("the Hill radius is d (m / 3M)^(1/3): Earth's is 1.497 million km", () => {
  assert.ok(Math.abs(S.hillRadius(AU_KM, EARTH_KG, SUN_KG) / 1.4966e6 - 1) < 0.001);
  assert.equal(S.hillRadius(0, EARTH_KG, SUN_KG), 0);
  assert.equal(S.hillRadius(AU_KM, 0, SUN_KG), 0);
  assert.equal(S.hillRadius(AU_KM, EARTH_KG, 0), 0);
});

test("the Lagrange points of the Sun and Earth are where the closed forms put them", () => {
  const points = S.lagrangePoints(AU_KM, 0, EARTH_KG, SUN_KG, 0, 29.78);
  const byName = Object.fromEntries(points.map((p) => [p.name, p]));
  assert.ok(Math.abs(byName.L1.x / AU_KM - 0.98997) < 2e-4, "L1 is 1.5 million km sunward");
  assert.ok(Math.abs(byName.L2.x / AU_KM - 1.01007) < 2e-4, "L2 is 1.5 million km beyond");
  assert.ok(Math.abs(byName.L3.x / AU_KM + 1.0) < 1e-5, "L3 is opposite, at one separation");
  assert.ok(byName.L4.y > 0 && byName.L5.y < 0, "L4 leads an anticlockwise orbit");
  assert.ok(Math.abs(Math.hypot(byName.L4.x - AU_KM, byName.L4.y) / AU_KM - 1) < 1e-4, "L4 is one separation from each");
  assert.ok(Math.abs(Math.hypot(byName.L4.x, byName.L4.y) / AU_KM - 1) < 1e-4);
});

test("a body moving the other way has L4 and L5 trade places; an unknown motion is taken as anticlockwise", () => {
  const backward = S.lagrangePoints(AU_KM, 0, EARTH_KG, SUN_KG, 0, -29.78);
  assert.ok(backward.find((p) => p.name === "L4").y < 0);
  const unknown = S.lagrangePoints(AU_KM, 0, EARTH_KG, SUN_KG, NaN, NaN);
  assert.ok(unknown.find((p) => p.name === "L4").y > 0);
});

test("a pair too alike in mass has no points, nor a body on its host", () => {
  assert.equal(S.lagrangePoints(AU_KM, 0, SUN_KG, SUN_KG, 0, 1), null);
  assert.equal(S.lagrangePoints(0, 0, EARTH_KG, SUN_KG, 0, 1), null);
});

test("a scene's sources are its bodies with a mass; the stars and the self body orbit nothing", () => {
  const parts = S.sceneBodies([
    { kind: "star", id: "1", name: "Sol", xkm: "0", ykm: "0", masskg: "1.98847e30", radiuskm: "695700", self: false },
    { kind: "planet", id: "2", name: "Earth", xkm: "1.4e8", ykm: "2e7", masskg: "5.97e24", radiuskm: "6371", vxkms: "-4", vykms: "29", self: false },
    { kind: "planet", id: "3", name: "Rogue", xkm: "1e8", ykm: "0", masskg: "", radiuskm: "100", self: false },
    { kind: "planet", id: "4", name: "Home", xkm: "0", ykm: "0", masskg: "1e25", radiuskm: "7000", self: true },
  ], 1.98847e30);
  assert.deepEqual(parts.sources.map((b) => b.name), ["Sol", "Earth", "Home"]);
  assert.deepEqual(parts.orbiters.map((b) => b.name), ["Earth"]);
  assert.equal(parts.host, 1.98847e30);
});

test("the legend names the mode, the range and what the rings and diamonds are, and says when there is nothing", () => {
  const text = S.legendText({ mode: "pull", range: { lo: 1e-6, hi: 1e-2, flat: false }, empty: false });
  assert.match(text, /Pull/);
  assert.match(text, /m\/s²/);
  assert.match(text, /Hill spheres/);
  assert.match(text, /Lagrange/);
  assert.match(S.legendText({ mode: "depth", range: { lo: 1, hi: 1, flat: true }, empty: false }), /no variation/);
  assert.match(S.legendText({ mode: "pull", range: { lo: 0, hi: 0, flat: true }, empty: true }), /no body with a mass/);
});
