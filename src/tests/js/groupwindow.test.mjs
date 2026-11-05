// tests/js/groupwindow.test.mjs -- static/groupwindow.js (MAP.222): the group window's bounds, members and framing.

import { test } from "node:test";
import assert from "node:assert/strict";

import { installDom, importPage } from "./fakedom.mjs";

installDom();
const G = await importPage("groupwindow.js");

const cloud = (id, x, y, z, r, extra) => Object.assign({ type: "nebula", id, name: "N" + id, x, y, z, radius_pc: r }, extra);

test("the bounds hold every member sphere", () => {
  const spheres = [cloud(1, 0, 0, 0, 10), cloud(2, 30, 0, 0, 5), cloud(3, 0, -20, 4, 8)];
  const b = G.boundsOf(spheres);
  spheres.forEach((s) => {
    const d = Math.hypot(s.x - b.centre[0], s.y - b.centre[1], s.z - b.centre[2]);
    assert.ok(d + s.radius_pc <= b.radiusPc + 1e-9);
  });
  assert.equal(G.boundsOf([]), null);
  const one = G.boundsOf([cloud(1, 5, 6, 7, 3)]);
  assert.deepEqual(one.centre, [5, 6, 7]);
  assert.equal(one.radiusPc, 3);
});

test("a nebula's window holds it and the nebulae that reach it, and says so in words", () => {
  const subject = cloud(1, 0, 0, 0, 10);
  const near = cloud(2, 14, 0, 0, 6);
  const far = cloud(3, 200, 0, 0, 6);
  const win = G.createWindow("nebula", subject, [subject, near, far]);
  assert.equal(win.kind, "nebula");
  assert.equal(win.key, "nebula:1");
  assert.equal(win.title, "N1");
  assert.equal(win.count, 2);
  assert.equal(win.summary, "2 nebulae");
  assert.ok(G.contains(win, 14, 0, 0) && !G.contains(win, 200, 0, 0));
  assert.equal(G.createWindow("nebula", subject, [subject]).summary, "1 nebula");
});

test("a window never grows past three times the subject's radius", () => {
  const subject = cloud(1, 0, 0, 0, 2);
  const huge = cloud(2, 3, 0, 0, 40);
  const win = G.createWindow("nebula", subject, [huge]);
  assert.equal(win.count, 1);
  assert.ok(win.radiusPc <= 3 * 2 + 1e-9);
});

test("a window crossing a sector border is bounded by its own sphere, not by sectors", () => {
  const subject = cloud(1, 3.9, 0, 0, 10);    // sectors are 4 pc here; the nebula reaches into several
  const win = G.createWindow("nebula", subject, []);
  assert.ok(G.contains(win, -5, 0, 0) && G.contains(win, 12, 0, 0) && !G.contains(win, 40, 0, 0));
});

test("framing is the bounds with the margin; no window contains everything; unknown kinds make none", () => {
  const win = G.createWindow("nebula", cloud(1, 1, 2, 3, 10), []);
  const f = G.framing(win);
  assert.deepEqual(f.center, [1, 2, 3]);
  assert.ok(Math.abs(f.radius - 10 * G.FRAME_MARGIN) < 1e-9);
  assert.equal(G.contains(null, 1e9, 0, 0), true);
  assert.equal(G.createWindow("constellation", cloud(1, 0, 0, 0, 1), []), null);
});

test("a second kind plugs in", () => {
  G.registerKind("association", {
    label: (c) => c.name, members: (c) => [c], noun: (n) => n + " associations",
  });
  assert.ok(G.kindNames().includes("association") && G.kindNames().includes("nebula"));
  const win = G.createWindow("association", cloud(9, 0, 0, 0, 5), []);
  assert.equal(win.kind, "association");
  assert.equal(win.summary, "1 associations");
});

test("a cluster's window is its own sphere and says what it is (MAP.223)", () => {
  const cluster = { ref: "globular_cluster:FE-1", name: "Globular cluster FE-1", x: 10, y: 20, z: 30, radius_pc: 25 };
  assert.ok(G.kindNames().includes("cluster"));
  const win = G.createWindow("cluster", cluster, [cloud(1, 12, 20, 30, 5)]);
  assert.equal(win.kind, "cluster");
  assert.equal(win.key, "globular_cluster:FE-1");
  assert.equal(win.title, "Globular cluster FE-1");
  assert.equal(win.summary, "a globular cluster");
  assert.equal(win.count, 1, "a nebula that reaches it does not join");
  assert.deepEqual(win.centre, [10, 20, 30]);
  assert.equal(win.radiusPc, 25);
  assert.ok(G.contains(win, 10, 20, 54) && !G.contains(win, 10, 20, 60));
  assert.equal(G.createWindow("cluster", { ...cluster, name: null }, []).title, "Globular cluster");
});
