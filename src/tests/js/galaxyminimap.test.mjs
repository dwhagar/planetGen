// tests/js/galaxyminimap.test.mjs -- static/galaxyminimap.js (MAP.75): the mini map's square, its locked camera and its drawing.

import { test } from "node:test";
import assert from "node:assert/strict";

import { importPage } from "./fakedom.mjs";

const M = await importPage("galaxyminimap.js");

test("the square sits in the bottom right corner, a share of the short side held to its bounds", () => {
  const rect = M.miniRect(1200, 800);
  assert.equal(rect.side, Math.round(800 * M.MINI_SHARE));
  assert.equal(rect.x, 1200 - rect.side - M.MINI_MARGIN);
  assert.equal(rect.y, M.MINI_MARGIN);
  assert.equal(M.miniRect(4000, 3000).side, M.MINI_MAX_SIDE);
  assert.equal(M.miniRect(M.MINI_MIN_CANVAS, 2000).side, M.MINI_MIN_SIDE);
});

test("a canvas too small for it has no mini map", () => {
  assert.equal(M.miniRect(M.MINI_MIN_CANVAS - 1, 900), null);
  assert.equal(M.miniRect(900, 200), null);
});

test("the camera stands where a sphere of the galaxy's radius fills its share of the square, and the locked policy keeps it there", () => {
  const radius = 15000;
  const distance = M.miniFitDistance(radius, 50);
  // The sphere's silhouette angle at that distance: asin(R / d) must be the fitted share of half the field of view.
  const half = (50 * Math.PI) / 360;
  assert.ok(Math.abs(Math.asin(radius / distance) - Math.asin(0.86 * Math.sin(half))) < 1e-9);
  assert.equal(M.MINI_POLICY.kind, "locked");
});

test("the eye looks along the main camera's own direction at the fitted distance", () => {
  const eye = M.miniEyeAt({ theta: 0.7, phi: 1.1 }, 1000);
  assert.ok(Math.abs(Math.hypot(...eye) - 1000) < 1e-9);
  assert.ok(Math.abs(eye[2] - 1000 * Math.cos(1.1)) < 1e-9);
  assert.ok(Math.abs(Math.atan2(eye[1], eye[0]) - 0.7) < 1e-9);
});

test("the outline holds three rings, the disk's two edge rings and two axes; the marker a cross and a line to the eye", () => {
  const outline = M.outlineSegments(100, 16);
  assert.equal(outline.length / 6, 5 * 16 + 2);
  const radii = [];
  for (let i = 0; i < 16; i++) radii.push(Math.hypot(outline[i * 6], outline[i * 6 + 1]));
  assert.ok(radii.every((r) => Math.abs(r - 100 / 3) < 1e-9));
  const marker = M.markerSegments([10, 20, 30], [1, 2, 3], 4);
  assert.equal(marker.length / 6, 4);
  assert.deepEqual(marker.slice(18), [10, 20, 30, 1, 2, 3]);
});
