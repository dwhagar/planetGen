// tests/js/systemgroups.test.mjs -- static/systemgroups.js (MAP.428): the stars of one system shown as one system.

import { test } from "node:test";
import assert from "node:assert/strict";

import { importPage } from "./fakedom.mjs";

const G = await importPage("systemgroups.js");

const star = (id, systemId, name, lum, x = 0) => ({ id, system_id: systemId, name, luminosity_sol: lum, x, y: 0, z: 0 });
const OFRIGGOR = [star(1, 7, "Ofriggor A", 0.591, 0), star(2, 7, "Ofriggor B", 0.01, 1), star(3, 7, "Ofriggor C", 0.012, 1), star(4, 7, "Ofriggor D", 0.0352, 8)];

test("only a system with two or more stars is a group, brightest member first", () => {
  const groups = G.groupBySystem([...OFRIGGOR, star(5, 9, "Loner", 1), star(6, null, "Scattered", 2), { id: 7, system_id: 3, phenomenon: true }]);
  assert.equal(groups.length, 1);
  assert.equal(groups[0].systemId, 7);
  assert.deepEqual(groups[0].members.map((s) => s.id), [1, 4, 3, 2]);
  assert.ok(Math.abs(groups[0].luminositySol - 0.6482) < 1e-9);
});

test("the system's name is the members' name without the component letter", () => {
  const [group] = G.groupBySystem(OFRIGGOR);
  assert.equal(G.systemNameOf(group), "Ofriggor");
  const [carried] = G.groupBySystem([star(1, 1, "Ofriggor", 1), star(2, 1, "Ofriggor B", 0.1)]);
  assert.equal(G.systemNameOf(carried), "Ofriggor", "the primary already carries the system's name");
  const [unnamed] = G.groupBySystem([star(1, 2, null, 1), star(2, 2, "", 0.1)]);
  assert.equal(G.systemNameOf(unnamed), "");
  const [plain] = G.groupBySystem([star(1, 3, "Alpha Centauri", 1), star(2, 3, "Proxima", 0.1)]);
  assert.equal(G.systemNameOf(plain), "Proxima", "no letter to take off: the shortest name stands");
});

test("a name consisting of one letter after a word is the only thing taken off", () => {
  const [group] = G.groupBySystem([star(1, 4, "Star X", 1), star(2, 4, "Star X B", 0.1)]);
  assert.equal(G.systemNameOf(group), "Star");
});

test("links run from the brightest member to the others while they are apart on screen", () => {
  const [group] = G.groupBySystem(OFRIGGOR);
  const project = (s) => ({ x: 100 + s.x * 30, y: 50, visible: true });
  const lines = G.linksOf(group, project);
  assert.equal(lines.length, 3);
  assert.ok(lines.every((line) => line.from.x === 100 && line.from.y === 50));
  assert.deepEqual(lines.map((line) => line.to.x).sort((a, b) => a - b), [130, 130, 340]);
});

test("no link while the members are one dot, or too far apart, or when fewer than two are on screen", () => {
  const [group] = G.groupBySystem(OFRIGGOR);
  assert.deepEqual(G.linksOf(group, (s) => ({ x: 100 + s.x * 0.1, y: 5, visible: true })), []);
  assert.deepEqual(G.linksOf(group, (s) => ({ x: 100 + s.x * 5000, y: 5, visible: true })), []);
  assert.deepEqual(G.linksOf(group, (s) => ({ x: 100 + s.x * 30, y: 5, visible: s.id === 1 })), []);
});

test("the shared name is anchored at the middle of the members", () => {
  const [group] = G.groupBySystem([star(1, 5, "A", 1, 0), star(2, 5, "B", 1, 4)]);
  assert.deepEqual(G.groupCenter(group), [2, 0, 0]);
});
