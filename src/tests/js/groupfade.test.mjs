// tests/js/groupfade.test.mjs -- static/groupfade.js (MAP.152): a scene group drawn at a share of its strength.

import { test } from "node:test";
import assert from "node:assert/strict";

import { importPage } from "./fakedom.mjs";

const G = await importPage("groupfade.js");

function fakeGroup(materials) {
  const objects = materials.map((material) => ({ material }));
  return { visible: true, traverse: (visit) => { visit({}); objects.forEach(visit); } };
}

const plain = () => ({ opacity: 0.8, transparent: false, userData: {} });
const shader = () => ({
  isShaderMaterial: true, userData: {},
  uniforms: { fade: { value: 1 }, maxAlpha: { value: 0.4 }, glowStrength: { value: 2 }, other: { value: 7 } },
});

test("every strength is scaled from its built value, and the shares do not compound", () => {
  const a = plain();
  const b = shader();
  const group = fakeGroup([a, b]);
  assert.equal(G.fadeGroup(group, 0.5), 0.5);
  assert.equal(a.opacity, 0.4);
  assert.equal(a.transparent, true);
  assert.equal(b.uniforms.fade.value, 0.5);
  assert.equal(b.uniforms.maxAlpha.value, 0.2);
  assert.equal(b.uniforms.glowStrength.value, 1);
  assert.equal(b.uniforms.other.value, 7, "an unrelated uniform is left alone");
  G.fadeGroup(group, 0.5);
  assert.equal(a.opacity, 0.4);
  assert.equal(b.uniforms.maxAlpha.value, 0.2);
});

test("0 hides the group, 1 brings back what was built", () => {
  const a = plain();
  const b = shader();
  const group = fakeGroup([a, b]);
  G.fadeGroup(group, 0);
  assert.equal(group.visible, false);
  G.fadeGroup(group, 1);
  assert.equal(group.visible, true);
  assert.equal(a.opacity, 0.8);
  assert.equal(a.transparent, false, "opaque again once whole");
  assert.equal(b.uniforms.maxAlpha.value, 0.4);
  assert.equal(b.uniforms.fade.value, 1);
});

test("the share is held to 0..1", () => {
  const group = fakeGroup([plain()]);
  assert.equal(G.fadeGroup(group, 7), 1);
  assert.equal(G.fadeGroup(group, -3), 0);
});
