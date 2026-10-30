// tests/js/spin.test.mjs -- static/spin.js, the slow turn of a rendered body on its axis (MAP.328).

import { test } from "node:test";
import assert from "node:assert/strict";

import { installDom, STATIC_URL } from "./fakedom.mjs";

installDom("http://localhost/");
const THREE = await import(new URL("vendor/three.module.min.js", STATIC_URL).href);
const S = await import(new URL("spin.js", STATIC_URL).href);

const near = (a, b, eps = 1e-9) => assert.ok(Math.abs(a - b) <= eps * Math.max(1, Math.abs(b)), `${a} is not ${b}`);

test("a day of a few hours shows about a turn in half a minute, never faster than the cap", () => {
  near(S.secondsPerTurn(S.FASTEST_PERIOD_HOURS), S.FASTEST_SECONDS_PER_TURN);
  near(S.secondsPerTurn(0.01), S.FASTEST_SECONDS_PER_TURN, 1e-9);
  const jupiter = S.secondsPerTurn(9.9);
  assert.ok(jupiter > 20 && jupiter < 60, `${jupiter}`);
});

test("the slowest period shows one turn in ten minutes and nothing is slower", () => {
  near(S.secondsPerTurn(S.SLOWEST_PERIOD_HOURS), S.SLOWEST_SECONDS_PER_TURN);
  near(S.secondsPerTurn(S.SLOWEST_PERIOD_HOURS * 3), S.SLOWEST_SECONDS_PER_TURN, 1e-9);
  assert.equal(S.FASTEST_SECONDS_PER_TURN, 20);
  assert.equal(S.SLOWEST_SECONDS_PER_TURN, 600);
});

test("one law for every body keeps the order of speeds: a longer period never turns faster", () => {
  const periods = [0.5, 2, 5, 10, 24, 100, 600, 1000, 5000, 20000];
  const rates = periods.map((p) => S.rateOf({ spinH: p }));
  for (let n = 1; n < rates.length; n += 1) assert.ok(rates[n] <= rates[n - 1], `${periods[n]} h is faster than ${periods[n - 1]} h`);
  assert.ok(rates[3] > rates[4] && rates[4] > rates[5], "inside the range the order is strict");
});

test("a body with no stored period is left still, never given a made-up one", () => {
  for (const look of [null, {}, { spinH: null }, { spinH: 0 }, { spinH: -5 }, { spinH: NaN }]) {
    assert.equal(S.rateOf(look), 0);
    assert.equal(S.secondsPerTurn(look && look.spinH), null);
    near(S.angleAt(look, 1000), S.phaseOf(look));
  }
});

test("the angle grows at the rate from the phase of the seed", () => {
  const look = { spinH: 24, seed: 4242 };
  near(S.angleAt(look, 10) - S.angleAt(look, 0), S.rateOf(look) * 10);
  near(S.angleAt(look, 0), S.phaseOf(look));
  assert.notEqual(S.phaseOf({ seed: 1 }), S.phaseOf({ seed: 2 }), "bodies do not start on the same face");
});

test("the axis is the stored one made a unit vector, else straight up", () => {
  assert.deepEqual(S.axisOf({ axis: [0, 3, 4] }).map((v) => Math.round(v * 1000) / 1000), [0, 0.6, 0.8]);
  assert.deepEqual(S.axisOf({}), [0, 1, 0]);
  assert.deepEqual(S.axisOf({ axis: [0, 0, 0] }), [0, 1, 0]);
});

test("the picture's axis leans with the axial tilt, and past 90 degrees points down", () => {
  assert.deepEqual(S.pictureAxisOf({}), [0, 1, 0]);
  const leaning = S.pictureAxisOf({ tiltDeg: 23.5 });
  near(leaning[0], Math.sin((23.5 * Math.PI) / 180));
  assert.ok(leaning[1] > 0);
  assert.ok(S.pictureAxisOf({ tiltDeg: 177 })[1] < 0, "retrograde: the axis points down");
});

// Where a point on the equator of a turned sphere has gone: +1 if it went counter-clockwise seen from above the axis.
function sweep(look, axis) {
  const mesh = new THREE.Mesh(new THREE.SphereGeometry(1, 8, 6));
  S.turn(THREE, mesh, look, 0, axis);
  const a = new THREE.Vector3(1, 0, 0).applyQuaternion(mesh.quaternion);
  S.turn(THREE, mesh, look, 1, axis);
  const b = new THREE.Vector3(1, 0, 0).applyQuaternion(mesh.quaternion);
  return new THREE.Vector3().crossVectors(a, b).y;
}

test("a retrograde body turns the other way about the same view axis", () => {
  const forward = sweep({ spinH: 24, seed: 0, tiltDeg: 20 }, S.pictureAxisOf({ tiltDeg: 20 }));
  const backward = sweep({ spinH: 24, seed: 0, tiltDeg: 160 }, S.pictureAxisOf({ tiltDeg: 160 }));
  assert.ok(forward > 0 && backward < 0, `${forward} ${backward}`);
});

test("turn lays the sphere's poles on the axis and turns it about it", () => {
  const mesh = new THREE.Mesh(new THREE.SphereGeometry(1, 8, 6));
  const look = { spinH: 24, seed: 0, axis: [1, 0, 0] };
  S.turn(THREE, mesh, look, 5);
  const pole = new THREE.Vector3(0, 1, 0).applyQuaternion(mesh.quaternion);
  near(pole.x, 1, 1e-9);
  const q0 = mesh.quaternion.clone();
  S.turn(THREE, mesh, look, 6);
  const relative = mesh.quaternion.clone().multiply(q0.invert());
  near(2 * Math.acos(Math.min(1, Math.abs(relative.w))), S.rateOf(look), 1e-6);
});

test("spin time runs only while bodies turn, and a long gap counts for a second at most", () => {
  S.setSwitch(true);
  const start = S.tick(1000);
  near(S.tick(1500) - start, 0.5);
  near(S.tick(1500 + 60000) - start, 1.5, 1e-9);
  S.setSwitch(false);
  assert.equal(S.active(), false);
  const stopped = S.tick(70000);
  assert.equal(S.tick(70500), stopped, "switched off, spin time stands still");
  S.setSwitch(true);
  assert.equal(S.switchedOn(), true);
  near(S.tick(71000) - stopped, 0.5);
});

test("reduced motion stops the bodies whatever the switch says", () => {
  S.setSwitch(true);
  globalThis.matchMedia = (query) => ({ matches: /reduce/.test(query) });
  assert.equal(S.reducedMotion(), true);
  assert.equal(S.active(), false);
  const held = S.tick(200000);
  assert.equal(S.tick(201000), held);
  delete globalThis.matchMedia;
  assert.equal(S.active(), true);
});

test("a point of a body's picture is found in the body's own frame, its axis leaning by its tilt", () => {
  // Upright: the middle of the picture is longitude 0, latitude 0; the right-hand edge is 90 degrees east; the top is the pole.
  const [lon, lat] = S.bodyPoint({}, 0, 0, 1);
  near(lon, 0); near(lat, 0);
  near(S.bodyPoint({}, 1, 0, 0)[0], Math.PI / 2);
  near(S.bodyPoint({}, 0, 1, 0)[1], Math.PI / 2);
  // Tilted 90 degrees the axis lies along the picture's right-hand side: that edge is the pole, the top is on the equator.
  near(S.bodyPoint({ tiltDeg: 90 }, 1, 0, 0)[1], Math.PI / 2);
  near(S.bodyPoint({ tiltDeg: 90 }, 0, 1, 0)[1], 0);
  // Past 90 degrees the picture shows the south pole at the top: a body that turns the other way.
  near(S.bodyPoint({ tiltDeg: 180 }, 0, 1, 0)[1], -Math.PI / 2);
});

test("a slow rotator still turns visibly and a year-long day is told apart from a 47-day one (UX.244)", () => {
  const forty7 = S.secondsPerTurn(46.8 * 24);
  const venus = S.secondsPerTurn(243 * 24);
  assert.ok(forty7 > 100 && forty7 < 350, `${forty7}`);
  assert.ok(venus > forty7 && venus < S.SLOWEST_SECONDS_PER_TURN, "not clamped to the slowest turn");
  assert.ok(S.rateOf({ spinH: 46.8 * 24 }) > 2 * Math.PI / 400, "a quarter of a degree a second at least");
});

test("a locked body turns to keep its longitude 0 toward what it orbits, wherever it is, about any axis", () => {
  const look = { spinH: 1123, seed: 7, locked: true };
  for (const axis of [[0, 1, 0], [0.3, 0.9, 0.2], [1, 0, 0], [0, -1, 0], [0.6, -0.6, 0.5]]) {
    for (const toward of [[1, 0, 0], [0, 0, -1], [-0.5, 0.2, 0.7], [0.2, 0.9, -0.1]]) {
      const mesh = new THREE.Mesh(new THREE.SphereGeometry(1, 8, 6));
      S.turn(THREE, mesh, look, 12345, axis, toward);
      const face = new THREE.Vector3(1, 0, 0).applyQuaternion(mesh.quaternion);
      const pole = new THREE.Vector3(0, 1, 0).applyQuaternion(mesh.quaternion);
      const ax = new THREE.Vector3(...axis).normalize();
      near(pole.dot(ax), 1, 1e-9);
      // The face is the direction to the parent as far as the axis allows: its part off the axis lies along the parent's.
      const t = new THREE.Vector3(...toward).normalize();
      const wanted = t.clone().sub(ax.clone().multiplyScalar(t.dot(ax)));
      if (wanted.length() < 1e-6) continue;
      near(face.clone().normalize().dot(wanted.normalize()), 1, 1e-9);
    }
  }
});

test("a locked body is not turned by the display law, and an unlocked one ignores the direction given", () => {
  const mesh = new THREE.Mesh(new THREE.SphereGeometry(1, 8, 6));
  S.turn(THREE, mesh, { spinH: 1123, seed: 7, locked: true }, 0, null, [0, 0, 1]);
  const first = new THREE.Vector3(1, 0, 0).applyQuaternion(mesh.quaternion);
  S.turn(THREE, mesh, { spinH: 1123, seed: 7, locked: true }, 500, null, [0, 0, 1]);
  const later = new THREE.Vector3(1, 0, 0).applyQuaternion(mesh.quaternion);
  near(first.distanceTo(later), 0, 1e-9);
  const free = { spinH: 24, seed: 7 };
  S.turn(THREE, mesh, free, 0, null, [0, 0, 1]);
  const a = new THREE.Vector3(1, 0, 0).applyQuaternion(mesh.quaternion);
  S.turn(THREE, mesh, free, 50, null, [0, 0, 1]);
  assert.ok(a.distanceTo(new THREE.Vector3(1, 0, 0).applyQuaternion(mesh.quaternion)) > 0.1);
  assert.equal(S.isLocked({ locked: true }), true);
  assert.equal(S.isLocked({}), false);
  assert.equal(S.isLocked(null), false);
});
