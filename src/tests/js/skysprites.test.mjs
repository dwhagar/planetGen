// tests/js/skysprites.test.mjs -- static/skysprites.js (VIEW.36): which nebulae and remnants of the sky are drawn and how, the quads, the shader text.

import { test } from "node:test";
import assert from "node:assert/strict";

import { importPage } from "./fakedom.mjs";

const S = await importPage("skysprites.js");

const near = (a, b, eps = 1e-9) => Math.abs(a - b) <= eps * Math.max(1, Math.abs(a), Math.abs(b));

function record(over) {
  return {
    id: 1, family: "emission", direction: [1, 0, 0], angularRadiusRad: 0.05, inside: false, extinctionMag: 0, surfaceBrightness: 19,
    visible: true, seen: true, ...over,
  };
}

test("the glow is nothing at the limit and its full strength five magnitudes brighter", () => {
  assert.equal(S.glowStrength(S.SB_LIMIT), 0);
  assert.equal(S.glowStrength(S.SB_LIMIT + 2), 0);
  assert.ok(near(S.glowStrength(S.SB_LIMIT - S.SB_RAMP), S.GLOW_MAX));
  assert.ok(near(S.glowStrength(10), S.GLOW_MAX));
  assert.equal(S.glowStrength(NaN), 0);
  assert.equal(S.glowStrength(null), 0);
  assert.ok(S.glowStrength(21) > 0 && S.glowStrength(21) < S.GLOW_MAX);
});

test("the disc follows the chord of a ball and its mean over the disc is one", () => {
  assert.ok(near(S.chordMean(0), 2 / 3));
  assert.ok(near(S.profile(0, 0) , 1 / S.chordMean(0)));
  assert.equal(S.profile(1, 0), 0);
  let sum = 0;
  const n = 4000;
  for (let i = 0; i < n; i++) {
    const r = (i + 0.5) / n;
    sum += S.profile(r, 0) * 2 * r / n;
  }
  assert.ok(near(sum, 1, 1e-3));
});

test("a shell is a ring: nothing at the centre, brightest at the inner edge, and its mean over the disc is one", () => {
  const inner = 1 - S.SHELL_THICKNESS;
  assert.ok(S.profile(0, inner) < S.profile(inner * 0.999, inner));
  assert.ok(S.profile(0.5, inner) < 0.5 * S.profile(0.99, inner) || S.profile(0.5, inner) < S.profile(0.99, inner));
  let sum = 0;
  const n = 20000;
  for (let i = 0; i < n; i++) {
    const r = (i + 0.5) / n;
    sum += S.profile(r, inner) * 2 * r / n;
  }
  assert.ok(near(sum, 1, 2e-3));
});

test("a cloud's opacity is 1 - 10^(-0.4 A) through its middle and nothing at its edge", () => {
  assert.ok(near(S.darkOpacity(5, 0), 1 - Math.pow(10, -2)));
  assert.equal(S.darkOpacity(5, 1), 0);
  assert.ok(S.darkOpacity(2, 0.9) < S.darkOpacity(2, 0.1));
});

test("a nebula under the limit is a glow in its family's colour", () => {
  const look = S.spriteLook(record({}));
  assert.equal(look.kind, "glow");
  assert.equal(look.inner, 0);
  assert.ok(look.strength > 0);
  assert.equal(look.colour.length, 3);
  assert.deepEqual(look.colour, S.spriteLook(record({ family: "emission" })).colour);
  assert.notDeepEqual(look.colour, S.spriteLook(record({ family: "reflection" })).colour);
});

test("a remnant that has been seen is a ring; one not yet seen is not drawn", () => {
  const ring = S.spriteLook(record({ family: "supernova-remnant", surfaceBrightness: 20 }));
  assert.equal(ring.kind, "glow");
  assert.ok(near(ring.inner, 1 - S.SHELL_THICKNESS));
  assert.equal(S.spriteLook(record({ family: "supernova-remnant", seen: false })).kind, "none");
});

test("a faint nebula, one that holds the observer and one too wide are not drawn", () => {
  assert.equal(S.spriteLook(record({ visible: false, surfaceBrightness: 25 })).kind, "none");
  assert.equal(S.spriteLook(record({ surfaceBrightness: 25 })).kind, "none");
  assert.equal(S.spriteLook(record({ inside: true })).kind, "none");
  assert.equal(S.spriteLook(record({ angularRadiusRad: S.MAX_ANGLE + 0.1 })).kind, "none");
  assert.equal(S.spriteLook(record({ angularRadiusRad: 0 })).kind, "none");
  assert.equal(S.spriteLook(null).kind, "none");
});

test("a cloud with no light of its own dims what lies behind it", () => {
  const dark = S.spriteLook(record({ family: "dark", surfaceBrightness: null, visible: false, extinctionMag: 3 }));
  assert.equal(dark.kind, "dark");
  assert.equal(dark.extinctionMag, 3);
  assert.equal(S.spriteLook(record({ family: "dark", surfaceBrightness: null, extinctionMag: 0.001 })).kind, "none");
});

test("the tangent basis is three unit vectors at right angles", () => {
  for (const d of [[1, 0, 0], [0, 0, 1], [0.3, -0.5, 0.8], [0, 0, -2]]) {
    const { n, e1, e2 } = S.tangentBasis(d);
    const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    assert.ok(near(dot(n, n), 1));
    assert.ok(near(dot(e1, e1), 1));
    assert.ok(near(dot(e2, e2), 1));
    assert.ok(Math.abs(dot(n, e1)) < 1e-9 && Math.abs(dot(n, e2)) < 1e-9 && Math.abs(dot(e1, e2)) < 1e-9);
  }
});

test("a quad is four corners at tan(radius) from its centre on the tangent plane", () => {
  const r = record({ direction: [0, 0, 2], angularRadiusRad: 0.1 });
  const q = S.quads([{ record: r, look: S.spriteLook(r) }]);
  assert.equal(q.position.length, 12);
  assert.equal(q.index.length, 6);
  const k = Math.tan(0.1);
  for (let c = 0; c < 4; c++) {
    assert.ok(near(Math.hypot(q.offset[2 * c], q.offset[2 * c + 1]), k * Math.SQRT2, 1e-6));
    assert.ok(near(q.position[3 * c + 2], 1, 1e-6)); // plane at unit distance along +z
    assert.ok(near(q.theta[c], 0.1, 1e-6));
  }
});

test("the shaders hold the chord, the shell and the extinction laws", () => {
  assert.match(S.GLOW_FRAGMENT, /sqrt\(1\.0 - r \* r\)/);
  assert.match(S.GLOW_FRAGMENT, /inner/);
  assert.match(S.DARK_FRAGMENT, /pow\(10\.0, -0\.4/);
});

test("createSkySprites makes a glow mesh and a dark one from the records that are drawn, none for the rest", () => {
  const made = [];
  class Attr { constructor(a, n) { this.array = a; this.n = n; } }
  class Geo { setAttribute(k, v) { (this.attrs ||= {})[k] = v; } setIndex(v) { this.index = v; } dispose() {} }
  class Mat { constructor(o) { Object.assign(this, o); } dispose() {} }
  class Mesh { constructor(g, m) { this.geometry = g; this.material = m; this.userData = {}; this.position = { copy() {} }; made.push(this); } updateMatrixWorld() {} }
  const THREE = { BufferGeometry: Geo, BufferAttribute: Attr, ShaderMaterial: Mat, Mesh, AdditiveBlending: 2, NormalBlending: 1, DoubleSide: 2 };
  const records = [
    record({ id: 1 }), record({ id: 2, family: "supernova-remnant" }), record({ id: 3, visible: false, surfaceBrightness: 26 }),
    record({ id: 4, family: "dark", surfaceBrightness: null, visible: false, extinctionMag: 2 }),
  ];
  const sprites = S.createSkySprites(THREE, records, {});
  assert.deepEqual(sprites.count, { glow: 2, dark: 1 });
  assert.equal(sprites.objects.length, 2);
  assert.ok(sprites.glow.renderOrder < sprites.dark.renderOrder);
  assert.ok(sprites.dark.renderOrder < 0);
  assert.equal(S.createSkySprites(THREE, [], {}).objects.length, 0);
  sprites.update({ near: 0.1, far: 1000, position: {} });
  assert.ok(sprites.glow.material.uniforms.spriteRadius.value > 0);
});
