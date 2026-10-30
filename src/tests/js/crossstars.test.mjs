// tests/js/crossstars.test.mjs -- static/crossstars.js (MAP.329): the systems of the sectors round the open one.

import { test } from "node:test";
import assert from "node:assert/strict";

import { importPage } from "./fakedom.mjs";

const C = await importPage("crossstars.js");

const EDGE = 4;
const rec = (id, x) => ({ id, shell: 1, band: 0, wedge: id, x });
const centreOf = (r) => [r.x, 0, 0];

test("the sectors read are the ones near the system, nearest first, never the open one", () => {
  const records = [rec(1, 0), rec(2, 4), rec(3, -4), rec(4, 20), rec(5, 5.5)];
  const got = C.neighbourRecords(records, centreOf, [1.8, 0, 0], EDGE, 1).map((r) => r.id);
  assert.deepEqual(got, [2, 5, 3], "the far sector (4) and the open one (1) are left out");
});

test("no more than CROSS_SECTORS are read, and a record with no id is skipped", () => {
  const records = Array.from({ length: 20 }, (_, i) => rec(i + 10, (i % 3) * 0.1));
  records.push({ shell: 1, band: 0, wedge: 0, x: 0 });
  assert.equal(C.neighbourRecords(records, centreOf, [0, 0, 0], EDGE, null).length, C.CROSS_SECTORS);
});

test("a scene's systems are placed as the scene is built, with a binary's companion folded into its primary", () => {
  const data = {
    centerPc: [10, 20, 30], halfEdgePc: 2, sceneHalfPx: 160,
    stars: [
      { x: 80, y: 40, z: -160, name: "A", endpoint: "system:7", href: "/system/7" },
      { x: 81, y: 41, z: -160, name: "A2", endpoint: "system:7", href: "/system/7", companionOf: 0 },
      { x: 0, y: 0, z: 0, name: "B", endpoint: "system:8", href: "/system/8" },
      { x: 5, y: 5, z: 5, name: "no link" },
    ],
  };
  const got = C.starsOfScene(data, { shell: 3, band: 1, wedge: 2 });
  assert.equal(got.length, 2);
  const unit = 2 / 160;
  assert.deepEqual([got[0].x, got[0].y, got[0].z], [10 + 80 * unit, 20 - 40 * unit, 30 - 160 * unit]);
  assert.equal(got[0].endpoint, "system:7");
  assert.deepEqual(got[0].sector, { shell: 3, band: 1, wedge: 2 });
  assert.equal(got[1].name, "B");
  assert.deepEqual(C.starsOfScene(null, {}), []);
  assert.deepEqual(C.starsOfScene({ centerPc: [0, 0, 0], halfEdgePc: 0, stars: [] }, {}), []);
});
