// tests/js/galaxystrip.test.mjs -- static/galaxystrip.js (MAP.375), the Galaxy Map's linked side strip for height:
// which stars the wedge holds, the frame and its scales (heights drawn 2.2 times true), picking, and the height readout.

import { test } from "node:test";
import assert from "node:assert/strict";

import { importPage } from "./fakedom.mjs";

const S = await importPage("galaxystrip.js");

const near = (a, b, eps = 1e-9) => Math.abs(a - b) <= eps;
const star = (x, y, z, lum = 1) => ({ x, y, z, luminosity_sol: lum });

test("azimuth runs 0 to 360 from the x axis and differences take the short way round", () => {
  assert.equal(S.azimuthDeg(1, 0), 0);
  assert.ok(near(S.azimuthDeg(0, 1), 90));
  assert.ok(near(S.azimuthDeg(-1, 0), 180));
  assert.ok(near(S.azimuthDeg(0, -1), 270));
  assert.ok(near(S.angleDifferenceDeg(350, 10), 20));
  assert.ok(near(S.angleDifferenceDeg(10, 350), -20));
  assert.ok(near(Math.abs(S.angleDifferenceDeg(0, 180)), 180));
});

test("a wedge of 60 degrees holds 30 either side of its centre, across the 0 and 360 seam, and everything with none", () => {
  const wedge = { centreDeg: 350, widthDeg: 60 };
  assert.equal(S.inWedge(1, 0, wedge), true, "10 degrees from 350");
  assert.equal(S.inWedge(Math.cos(Math.PI / 18 * 2.4), Math.sin(Math.PI / 18 * 2.4), wedge), false, "43 degrees");
  assert.equal(S.inWedge(Math.cos(-Math.PI / 6), Math.sin(-Math.PI / 6), { centreDeg: 0, widthDeg: 60 }), true, "the edge counts");
  assert.equal(S.inWedge(-1, 0, wedge), false);
  assert.equal(S.inWedge(-1, 0, null), true);
  assert.equal(S.inWedge(0, 0, wedge), true, "the centre is in every wedge");
});

test("a wedge about the focus points at the focus, and at the fallback when the focus is the centre", () => {
  assert.ok(near(S.wedgeAbout(0, 5).centreDeg, 90));
  assert.equal(S.wedgeAbout(0, 5).widthDeg, S.WEDGE_DEG);
  assert.equal(S.wedgeAbout(0, 0, 60, 123).centreDeg, 123);
  assert.equal(S.wedgeAbout(0, 0).centreDeg, 0);
});

test("the strip's entries are the wedge's stars, brightest first, cut to the limit, with R and height", () => {
  const stars = [star(1000, 0, 50, 1), star(0, 1000, 10, 5), star(-1000, 0, 0, 9), star(2000, 100, -30, 3)];
  const found = S.stripEntries(stars, { centreDeg: 0, widthDeg: 60 });
  assert.deepEqual(found.map((e) => e.star.luminosity_sol), [3, 1]);
  assert.ok(near(found[1].r, 1000) && found[1].h === 50);
  assert.ok(near(found[0].r, Math.hypot(2000, 100)) && found[0].h === -30);
  assert.equal(S.stripEntries(stars, null).length, 4);
  assert.equal(S.stripEntries(stars, null, 2).length, 2);
  assert.deepEqual(S.stripEntries(stars, null, 2).map((e) => e.star.luminosity_sol), [9, 5]);
});

test("the frame fits the distance across and draws heights 2.2 times taller than true", () => {
  const entries = S.stripEntries([star(1000, 0, 0), star(9000, 0, 0)], null);
  const f = S.stripFrame(entries, 856, 160);
  assert.equal(f.rMin, 1000);
  assert.equal(f.rMax, 9000);
  assert.ok(near(f.x(1000), S.STRIP_PAD_PX));
  assert.ok(near(f.x(9000), 856 - S.STRIP_PAD_PX));
  const horizontal = f.x(2000) - f.x(1000);
  const vertical = f.y(0) - f.y(1000);
  assert.ok(near(vertical / horizontal, S.HEIGHT_EXAGGERATION), "1000 pc up is 2.2 times 1000 pc across");
  assert.equal(f.y(0), 80);
  assert.ok(f.y(100) < f.y(0) && f.y(-100) > f.y(0), "up is above the plane");
  assert.ok(near(f.r(f.x(4321)), 4321) && near(f.h(f.y(-77)), -77));
});

test("a frame never covers less than a few hundred parsecs, takes the selection in, and is something with no stars", () => {
  const one = S.stripFrame(S.stripEntries([star(8000, 0, 5)], null), 600, 160);
  assert.ok(one.rMax - one.rMin >= S.MIN_ACROSS_PC - 1e-9);
  assert.ok(one.x(8000) > 0 && one.x(8000) < 600);
  const withSelection = S.stripFrame(S.stripEntries([star(8000, 0, 5)], null), 600, 160, [{ x: 12000, y: 0, z: 0 }]);
  assert.equal(withSelection.rMax, 12000);
  const empty = S.stripFrame([], 600, 160);
  assert.ok(Number.isFinite(empty.pxPerPc) && empty.pxPerPc > 0);
});

test("a pick takes the nearest star within reach, the brightest of those level, and nothing beyond it", () => {
  const entries = S.stripEntries([star(1000, 0, 0, 1), star(1010, 0, 0, 50), star(5000, 0, 100, 2)], null);
  const f = S.stripFrame(entries, 600, 160);
  const far = S.pickInStrip(entries, f, f.x(5000) + 4, f.y(100) - 3, 10);
  assert.equal(far.star.luminosity_sol, 2);
  assert.equal(S.pickInStrip(entries, f, f.x(5000) + 30, f.y(100), 10), null, "beyond the reach");
  assert.equal(S.pickInStrip(entries, f, f.x(5000) + 15, f.y(100), 22).star.luminosity_sol, 2, "a finger reaches farther");
  const level = S.pickInStrip(entries, f, f.x(1005), f.y(0), 22);
  assert.equal(level.star.luminosity_sol, 50, "two level stars: the brighter wins");
  const near1 = S.pickInStrip(entries, f, f.x(1000) - 1, f.y(0), 22);
  assert.equal(near1.star.luminosity_sol, 1, "a clearly nearer star wins over a brighter one");
});

test("the height readout is exact and says how many scale heights of each disc it is", () => {
  const shape = { disk_scale_height_pc: 300, thick_disk_scale_height_pc: 900 };
  const r = S.heightReadout(-450, shape);
  assert.equal(r.heightPc, -450);
  assert.ok(near(r.thin, 1.5) && near(r.thick, 0.5));
  assert.equal(S.heightReadout(10, null).thin, null);
  assert.equal(S.readoutText(8200, 134, shape), "Height +134 pc above the plane (0.4 thin-disc and 0.1 thick-disc scale heights), 8.2 kpc from the centre");
  assert.equal(S.readoutText(450, -450, shape), "Height -450 pc below the plane (1.5 thin-disc and 0.5 thick-disc scale heights), 450 pc from the centre");
  assert.equal(S.readoutText(900, 0, null), "Height 0 pc on the plane, 900 pc from the centre");
});

test("the distance ticks are round numbers inside the range", () => {
  const ticks = S.distanceTicks(3000, 13000);
  assert.ok(ticks.length >= 3 && ticks.length <= 12);
  for (const t of ticks) assert.ok(t.r >= 3000 && t.r <= 13000 && /kpc$/.test(t.label), JSON.stringify(t));
  assert.equal(S.niceStep(10000, 6), 2000);
  const small = S.distanceTicks(100, 500);
  assert.ok(small.every((t) => /pc$/.test(t.label) && !/kpc/.test(t.label)));
});
