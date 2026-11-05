// tests/js/systemsky.test.mjs -- static/systemsky.js (VIEW.35): the preference, the URLs, the status lines, the turn from the galaxy's axes to the scene's, and the
// controller's flow with a fake view and fake routes.

import { test } from "node:test";
import assert from "node:assert/strict";

import { importPage } from "./fakedom.mjs";

const Y = await importPage("systemsky.js");

const near = (a, b, eps = 1e-12) => Math.abs(a - b) <= eps;

function store(initial) {
  const data = new Map(initial ? [[Y.STORAGE_KEY, initial]] : []);
  return { getItem: (k) => (data.has(k) ? data.get(k) : null), setItem: (k, v) => data.set(k, v), data };
}

test("the sky is on unless the browser remembers off, and a broken store means on", () => {
  assert.equal(Y.readPreference(null), true);
  assert.equal(Y.readPreference(store()), true);
  assert.equal(Y.readPreference(store("on")), true);
  assert.equal(Y.readPreference(store("off")), false);
  assert.equal(Y.readPreference({ getItem() { throw new Error("blocked"); } }), true);
  const s = store();
  Y.writePreference(s, false);
  assert.equal(s.data.get(Y.STORAGE_KEY), "off");
  Y.writePreference({ setItem() { throw new Error("blocked"); } }, true);   // no throw
});

test("a narrow window is a phone", () => {
  assert.equal(Y.isPhone({ innerWidth: 390 }), true);
  assert.equal(Y.isPhone({ innerWidth: 700 }), true);
  assert.equal(Y.isPhone({ innerWidth: 1280 }), false);
  assert.equal(Y.isPhone(null), false);
});

test("the sky and band URLs add their query to a path that may have one", () => {
  assert.equal(Y.skyUrlFor("/system/A/sky", 7.5), "/system/A/sky?limit=7.5");
  assert.equal(Y.skyUrlFor("/system/A/sky?db=x", 6.5), "/system/A/sky?db=x&limit=6.5");
  assert.equal(Y.bandUrlFor({ href: "/system/A/sky/band?limit=7.5", key: "abc123" }), "/system/A/sky/band?limit=7.5&key=abc123");
  assert.equal(Y.bandUrlFor({ href: "/b", key: "k" }), "/b?key=k");
});

test("the status line says what is drawn", () => {
  assert.match(Y.noteFor({ kind: "off" }), /off/);
  assert.match(Y.noteFor({ kind: "loading" }), /Loading/);
  assert.match(Y.noteFor({ kind: "ready", count: 24512 }), /24,512 stars/);
  assert.match(Y.noteFor({ kind: "filling", count: 900, report: { built_to_pc: 180 } }), /900 stars, built to 180 pc, filling in/);
  assert.match(Y.noteFor({ kind: "filling", count: 900, report: { built_to_pc: 2500 } }), /2\.5 kpc/);
  assert.match(Y.noteFor({ kind: "unavailable", status: 404 }), /no place in the galaxy/);
  assert.match(Y.noteFor({ kind: "unavailable", status: 502 }), /plain background/);
});

test("the galaxy's axes turn into the scene's: north is up, x stays x", () => {
  const [qx, qy, qz, qw] = Y.GALAXY_TO_SCENE;
  const rotate = ([x, y, z]) => {                   // v' = q v q^-1
    const ix = qw * x + qy * z - qz * y, iy = qw * y + qz * x - qx * z, iz = qw * z + qx * y - qy * x, iw = -qx * x - qy * y - qz * z;
    return [ix * qw + iw * -qx + iy * -qz - iz * -qy, iy * qw + iw * -qy + iz * -qx - ix * -qz, iz * qw + iw * -qz + ix * -qy - iy * -qx];
  };
  const close = (a, b) => a.every((v, i) => near(v, b[i], 1e-9));
  assert.ok(close(rotate([1, 0, 0]), [1, 0, 0]));
  assert.ok(close(rotate([0, 0, 1]), [0, 1, 0]));
  assert.ok(close(rotate([0, 1, 0]), [0, 0, -1]));
  assert.ok(near(qx * qx + qw * qw, 1));
});

// A fake page: routes, a view that records what is attached, and a THREE stand-in that the (fake) layer builder never touches.
function harness(routes, extra) {
  const log = { attached: [], detached: 0, notes: [], timers: [] };
  const view = { attachSky(layer) { log.attached.push(layer); }, detachSky() { log.detached += 1; } };
  const fetched = [];
  const fetch = async (url) => {
    fetched.push(url);
    const answer = routes(url);
    return { ok: answer.status === undefined, status: answer.status || 200, arrayBuffer: async () => answer.body };
  };
  return { log, view, fetch, fetched, extra };
}

test("a sky that cannot be had leaves the plain background and says why", async () => {
  const h = harness(() => ({ status: 404, body: new ArrayBuffer(0) }));
  const sky = Y.createSystemSky({}, { skyUrl: "/system/A/sky", view: h.view, fetch: h.fetch, storage: store(), onState: (s, n) => h.log.notes.push([s.kind, n]) });
  const state = await sky.reload();
  assert.equal(state.kind, "unavailable");
  assert.equal(state.status, 404);
  assert.equal(h.log.attached.length, 0);
  assert.match(h.log.notes.at(-1)[1], /no place in the galaxy/);
  assert.deepEqual(h.fetched, ["/system/A/sky?limit=7.5"]);
});

test("a phone asks for the lower limit, a bad file is unavailable too", async () => {
  const h = harness(() => ({ body: new Uint8Array([1, 2, 3, 4]).buffer }));
  const sky = Y.createSystemSky({}, { skyUrl: "/s", view: h.view, fetch: h.fetch, storage: store(), phone: true });
  const state = await sky.reload();
  assert.deepEqual(h.fetched, ["/s?limit=6.5"]);
  assert.equal(state.kind, "unavailable");
});

test("turning the sky off removes it and is remembered; turning it on asks again", async () => {
  const h = harness(() => ({ status: 502, body: new ArrayBuffer(0) }));
  const s = store();
  const states = [];
  const sky = Y.createSystemSky({}, { skyUrl: "/s", view: h.view, fetch: h.fetch, storage: s, onState: (st) => states.push(st.kind) });
  assert.equal(sky.enabled(), true);
  sky.disable();
  assert.equal(sky.enabled(), false);
  assert.equal(s.data.get(Y.STORAGE_KEY), "off");
  assert.equal(states.at(-1), "off");
  await sky.enable();
  assert.equal(s.data.get(Y.STORAGE_KEY), "on");
  assert.deepEqual(states, ["off", "loading", "unavailable"]);
  const again = Y.createSystemSky({}, { skyUrl: "/s", view: h.view, fetch: h.fetch, storage: store("off") });
  assert.equal(again.enabled(), false);
  assert.equal(again.state().kind, "off");
});
