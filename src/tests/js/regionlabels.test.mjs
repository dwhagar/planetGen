// tests/js/regionlabels.test.mjs -- static/regionlabels.js (MAP.314): which charted-region names the nodes of a view offer, by zoom.

import { test } from "node:test";
import assert from "node:assert/strict";

import { STATIC_URL } from "./fakedom.mjs";

const R = await import(new URL("regionlabels.js", STATIC_URL).href);

const region = (key, name, parent, centre, sectors, percent = sectors / 7.8, radius = 12) =>
  ({ key, name, parent, centre, radius, sectors, percent });
const nine = (name, centre, sectors = 100) => ({ name, parent: "Vang Reach", extent: { centre, radius: 12 }, charted: { sectors, percent: sectors / 7.8 } });

test("the level of a node key is its first number", () => {
  assert.equal(R.levelOf("81.3.0.1"), 81);
  assert.equal(R.levelOf("9.33.0.0"), 9);
  assert.ok(Number.isNaN(R.levelOf("galaxy")));
});

test("a charted region's own node offers its name at the middle of its sphere, ranked by the share charted, first of all", () => {
  const [label] = R.regionLabels([{ key: "9.4.0.0", node: nine("Vang Zarh Region", [10, 20, 30], 78) }]);
  assert.equal(label.id, "region:9.4.0.0");
  assert.equal(label.text, "Vang Zarh Region");
  assert.deepEqual([label.x, label.y, label.z], [10, 20, 30]);
  assert.deepEqual(label.ref, { kind: "region", key: "9.4.0.0", extent: { centre: [10, 20, 30], radius: 12 } });
  assert.deepEqual(label.slice, { domain: "region", shareCharted: 0.1, first: true });
});

test("a 9-block that is not a charted region offers nothing, nor does a node not loaded", () => {
  assert.deepEqual(R.regionLabels([
    { key: "9.5.0.0", node: { aggregate: {}, stars: [] } },                         // no name: no generated sector, or no naming key
    { key: "9.6.0.0", node: { name: "Odd", parent: "Odd Reach" } },                 // a name with no extent is not drawn
    { key: "9.7.0.0", node: undefined },                                            // still on its way
    { key: "9.8.0.0" },
  ]), []);
});

test("a 27-block offers each region it lists, once, even when the region's own node is in the view too", () => {
  const listed = [region("9.4.0.0", "Vang Zarh Region", "Vang Reach", [10, 0, 0], 100), region("9.10.0.0", "Vang Zarl Region", "Vang Reach", [50, 0, 0], 30)];
  const labels = R.regionLabels([{ key: "27.1.0.0", node: { regions: listed } },
    { key: "9.4.0.0", node: nine("Vang Zarh Region", [10, 0, 0], 100) }]);
  assert.deepEqual(labels.map((label) => label.text).sort(), ["Vang Zarh Region", "Vang Zarl Region"]);
  assert.equal(labels.find((label) => label.id === "region:9.10.0.0").ref.extent.radius, 12);
  assert.ok(labels.every((label) => label.slice.first));
});

test("a region listed with no radius is framed by its block", () => {
  const [label] = R.regionLabels([{ key: "27.1.0.0", node: { regions: [{ key: "9.4.0.0", name: "Vang Zarh Region", parent: "Vang Reach", centre: [1, 2, 3], sectors: 5, percent: 1 }] } }]);
  assert.equal(label.ref.extent, null);
  assert.equal(label.ref.key, "9.4.0.0");
});

test("an 81-block and a 243-block each offer one family name, the family with the most sectors, at the centroid of its regions", () => {
  const regions = [region("9.4.0.0", "Vang Zarh Region", "Vang Reach", [0, 0, 0], 100), region("9.5.0.0", "Vang Zarl Region", "Vang Reach", [30, 0, 0], 300),
    region("9.90.0.0", "Orr Beez Region", "Orr Reach", [500, 0, 0], 150)];
  for (const key of ["81.1.0.0", "243.0.0.0"]) {
    const labels = R.regionLabels([{ key, node: { regions } }]);
    assert.equal(labels.length, 1, key);
    assert.equal(labels[0].text, "Vang Reach");
    assert.equal(labels[0].id, "family:" + key);
    assert.ok(Math.abs(labels[0].x - 22.5) < 1e-9 && labels[0].y === 0 && labels[0].z === 0, "weighted by sectors, not the middle of the block");
    assert.equal(labels[0].ref.extent, null);
    assert.equal(labels[0].ref.key, key);
  }
});

test("a tie between families goes to the name that sorts first, and a region with no family adds nothing", () => {
  const labels = R.regionLabels([{ key: "81.1.0.0", node: { regions: [region("9.4.0.0", "B Region", "Bee Reach", [0, 0, 0], 10), region("9.5.0.0", "A Region", "Ann Reach", [9, 9, 9], 10),
    { key: "9.6.0.0", name: "Lone", centre: [1, 1, 1], sectors: 99, percent: 9 }] } }]);
  assert.deepEqual(labels.map((label) => label.text), ["Ann Reach"]);
});

test("a node of any other size offers nothing, whatever it carries, and a block with no charted region offers nothing", () => {
  const regions = [region("9.4.0.0", "Vang Zarh Region", "Vang Reach", [0, 0, 0], 100)];
  assert.deepEqual(R.regionLabels([{ key: "3.12.0.0", node: { regions } }, { key: "1.40.2.0", node: { regions } },
    { key: "729.0.0.0", node: { regions } }, { key: "2187.0.0.0", node: { regions } }, { key: "81.0.0.0", node: { stars: [] } },
    { key: "27.0.0.0", node: { regions: [] } }, { key: "243.0.0.0", node: { regions: [] } }]), []);
  assert.deepEqual(R.regionLabels([]), []);
});
