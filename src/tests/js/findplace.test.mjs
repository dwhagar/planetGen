// tests/js/findplace.test.mjs -- static/findplace.js (MAP.336): what a line typed in the find-a-place box can be.

import { test } from "node:test";
import assert from "node:assert/strict";

import { STATIC_URL } from "./fakedom.mjs";

const F = await import(new URL("findplace.js", STATIC_URL).href);
const S = await import(new URL("galaxystages.js", STATIC_URL).href);

const parseSector = S.parseSectorDesignation;
const SECTOR = S.sectorDesignation(30, 0, 1);
const SYSTEM = SECTOR + "-0000005-000";
const MOON = SECTOR + "-0000005-003";
const blobs = [
  { key: "a", name: "Fardoo Beez Region", x: 10, y: 20, z: 30, radiusPc: 12, block: { m: 9, shell: 5, band: 0, wedge: 1 }, parent: "Outer Reach" },
  { key: "b", name: "Beez Reach", x: 1, y: 2, z: 3, radiusPc: 15, block: { m: 9, shell: 6, band: 0, wedge: 2 }, parent: null },
  { key: "c", name: null, x: 0, y: 0, z: 0, radiusPc: 12, block: { m: 9, shell: 7, band: 0, wedge: 3 } },
];
const memory = () => {
  const values = {};
  return { getItem: (k) => (k in values ? values[k] : null), setItem: (k, v) => { values[k] = String(v); }, removeItem: (k) => { delete values[k]; } };
};

test("an object ID names its sector, system and body", () => {
  const system = F.objectIdOf(SYSTEM, parseSector);
  assert.deepEqual(system.sector, parseSector(SECTOR));
  assert.equal(system.system, SYSTEM);
  assert.equal(system.object, null);
  const moon = F.objectIdOf("moon:" + MOON.toLowerCase(), parseSector);
  assert.equal(moon.system, SYSTEM);
  assert.equal(moon.object, "moon:" + MOON);
  assert.equal(F.objectIdOf(SECTOR, parseSector), null);
  assert.equal(F.objectIdOf("Kanbee 3", parseSector), null);
  assert.equal(F.objectIdOf("12-34", parseSector), null);
});

test("regions are found by part of a name, the names that start with it first", () => {
  const found = F.regionResults(blobs, "beez");
  assert.deepEqual(found.map((r) => r.label), ["Beez Reach", "Fardoo Beez Region"]);
  assert.equal(found[1].detail, "region in Outer Reach");
  assert.deepEqual(found[1].target.region.center, [10, 20, 30]);
  assert.equal(found[1].target.region.radius, 12);
  assert.deepEqual(F.regionResults(blobs, ""), []);
  assert.equal(F.regionResults(blobs, "region c")[0].label, "Region c");
});

test("the visit log keeps the last places, newest first, once each, and survives storage that throws", () => {
  const log = F.createRecents(memory());
  for (let n = 0; n < 12; n++) log.add("Place " + n, { center: [n, 0, 0], radius: 1 });
  log.add("place 5", { center: [5, 5, 5], radius: 1 });
  const list = log.list();
  assert.equal(list.length, F.RECENT_LIMIT);
  assert.equal(list[0].label, "place 5");
  assert.equal(list.filter((r) => r.label.toLowerCase() === "place 5").length, 1);
  log.clear();
  assert.deepEqual(log.list(), []);
  const broken = F.createRecents({ getItem() { throw new Error("blocked"); }, setItem() { throw new Error("blocked"); }, removeItem() { throw new Error("blocked"); } });
  broken.add("Here", { center: [0, 0, 0], radius: 1 });
  assert.equal(broken.list()[0].label, "Here");
});

test("a typed line gives groups in order: this place, recent, regions, names, the galaxy", () => {
  const recents = [{ label: "Beez Station", target: { center: [1, 1, 1], radius: 1 } }];
  const names = [{ kind: "system", id: SYSTEM, name: "Beezlebub", shell_index: 5, band_index: 0, wedge_index: 1, parents: [] }];
  const groups = F.buildGroups({
    text: "beez", parsed: { name: "beez" }, id: null, regions: blobs, recents: recents, names: names,
    home: { label: "Whole galaxy", target: { home: true } },
  });
  assert.deepEqual(groups.map((g) => g.key), ["recent", "region", "name"]);
  const name = groups[2].items[0];
  assert.equal(name.target.open, true);
  assert.equal(name.target.system, SYSTEM);
  assert.deepEqual(name.target.sector, { shell: 5, band: 0, wedge: 1 });
});

test("coordinates, a designation and an object ID each give one result to go to", () => {
  const point = F.buildGroups({ text: "10 -20 30.5", parsed: { point: [10, -20, 30.5], sector: { shell: 1, band: 0, wedge: 0 } } });
  assert.equal(point[0].key, "here");
  assert.deepEqual(point[0].items[0].target.center, [10, -20, 30.5]);
  const designation = F.buildGroups({ text: "7/-2/14", parsed: { sector: { shell: 7, band: -2, wedge: 14 } } });
  assert.deepEqual(designation[0].items[0].target.sector, { shell: 7, band: -2, wedge: 14 });
  const id = F.objectIdOf("moon:" + MOON, parseSector);
  const object = F.buildGroups({ text: "moon:" + MOON, parsed: { name: "x" }, id: id });
  assert.equal(object[0].items[0].target.object, "moon:" + MOON);
  assert.equal(F.asksServer({ name: "x" }, "x"), false);
  assert.equal(F.asksServer({ name: "ka" }, "ka"), true);
  assert.equal(F.asksServer({ sector: {} }, "7/1/1"), false);
});

test("a blank line lists the visit log and the way home", () => {
  const groups = F.buildGroups({
    text: "", parsed: {}, regions: blobs, recents: [{ label: "A", target: { home: true } }], home: { label: "Whole galaxy", target: { home: true } },
  });
  assert.deepEqual(groups.map((g) => g.key), ["recent", "home"]);
});

test("a blank line also offers the charted regions with the most sectors, biggest first (MAP.357)", () => {
  const charted = [
    { key: "a", name: "Aye", block: { shell: 1 }, x: 1, y: 2, z: 3, radiusPc: 40, sectorCount: 2 },
    { key: "b", name: "Bee", block: { shell: 2 }, x: 4, y: 5, z: 6, radiusPc: 50, sectorCount: 9 },
    { key: "c", name: "Sea", block: { shell: 3 }, x: 7, y: 8, z: 9, radiusPc: 60, sectorCount: 0 },
  ];
  const groups = F.buildGroups({ text: "", parsed: {}, regions: charted, recents: [], home: { label: "Whole galaxy", target: { home: true } } });
  assert.deepEqual(groups.map((g) => g.key), ["charted", "home"]);
  assert.deepEqual(groups[0].items.map((i) => i.label), ["Bee", "Aye"]);
  assert.deepEqual(groups[0].items[0].target.region, { block: { shell: 2 }, center: [4, 5, 6], radius: 50 });
  assert.equal(F.GROUP_TITLES.charted, "Charted places");
});
