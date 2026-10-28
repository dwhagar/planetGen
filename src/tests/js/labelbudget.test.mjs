// tests/js/labelbudget.test.mjs -- static/labelbudget.js (MAP.317): the domain-sliced label budget.

import { test } from "node:test";
import assert from "node:assert/strict";

import { STATIC_URL } from "./fakedom.mjs";

const B = await import(new URL("labelbudget.js", STATIC_URL).href);

const star = (id, radiusSol, luminositySol, massSolar) => ({ id: "s" + id, domain: "star", radiusSol, luminositySol, massSolar });
const world = (id, category, habitability, population) => ({ id: "w" + id, domain: "inhabited", category, habitability, population });
const rock = (id, category, mass, volume) => ({ id: "u" + id, domain: "uninhabited", category, mass, volume });
const cloud = (id, volume) => ({ id: "n" + id, domain: "nebula", volume });
const region = (id, shareCharted, population) => ({ id: "r" + id, domain: "region", shareCharted, population });
const domainsOf = (chosen) => chosen.map((c) => c.domain);
const count = (chosen, domain) => chosen.filter((c) => c.domain === domain).length;

test("the budget is the screen area over the area one name takes", () => {
  assert.equal(B.labelBudget({ width: 1200, height: 600 }), Math.floor(1200 * 600 / (120 * 24)));
  assert.equal(B.labelBudget({ width: 1, height: 1 }), 1);
  assert.equal(B.labelBudget({ width: 240, height: 48 }, { width: 120, height: 24 }), 4);
});

test("a sparse view names everything it draws, and nothing undrawn", () => {
  const items = [star(1, 1, 1, 1), world(1, "planet", 7), rock(1, "gas giant", 0.001), cloud(1, 50), region(1, 0.14),
    Object.assign(star(2, 2, 3, 1.5), { drawn: false })];
  const chosen = B.chooseLabels(items, { budget: 100 });
  assert.equal(chosen.length, 5);
  assert.ok(!chosen.some((c) => c.id === "s2"), "a star that is not drawn is not named");
  assert.deepEqual(chosen.map((c) => c.rank), [0, 1, 2, 3, 4]);
});

test("a dense view gives each domain its share, not only the brightest stars", () => {
  const items = [];
  for (let i = 0; i < 400; i++) items.push(star(i, 1 + i / 100, 1 + i, 1 + i / 50));
  for (let i = 0; i < 200; i++) items.push(world(i, "planet", i % 10));
  for (let i = 0; i < 200; i++) items.push(rock(i, "rocky world", i / 1000));
  for (let i = 0; i < 100; i++) items.push(cloud(i, 10 + i));
  const chosen = B.chooseLabels(items, { budget: 100 });
  assert.equal(chosen.length, 100);
  // Shares 30:30:25:10 (the regions have no items, so their 5 percent goes round).
  const shares = ["inhabited", "star", "uninhabited", "nebula"].map((d) => count(chosen, d));
  assert.ok(shares[0] >= 30 && shares[1] >= 30 && shares[2] >= 25 && shares[3] >= 10, shares.join());
  assert.ok(shares[0] <= 35 && shares[1] <= 35 && shares[2] <= 30 && shares[3] <= 14, shares.join());
});

test("a domain with less than its share lends the rest to the others", () => {
  const items = [world(1, "planet", 5)];
  for (let i = 0; i < 200; i++) items.push(star(i, 1 + i / 100, 1 + i));
  const chosen = B.chooseLabels(items, { budget: 60 });
  assert.equal(chosen.length, 60);
  assert.equal(count(chosen, "inhabited"), 1);
  assert.equal(count(chosen, "star"), 59);
});

test("stars are taken a third each by size, brightness and mass, once each", () => {
  // Star 1 is the largest, 2 the brightest, 3 the heaviest, 4 second in all three.
  const items = [star(1, 100, 1, 1), star(2, 1, 1000, 1), star(3, 1, 1, 50), star(4, 50, 500, 25), star(5, 0.1, 0.1, 0.1)];
  const chosen = B.chooseLabels(items, { budget: 4 });
  assert.deepEqual(chosen.slice(0, 3).map((c) => c.id), ["s1", "s2", "s3"]);
  assert.equal(chosen[3].id, "s4");
  assert.equal(new Set(chosen.map((c) => c.id)).size, 4);
});

test("a star without a mass is ranked by the main-sequence estimate of its brightness", () => {
  assert.equal(B.starMass({ massSolar: 2 }), 2);
  assert.ok(Math.abs(B.starMass({ luminositySol: 1 }) - 1) < 1e-12);
  assert.ok(B.starMass({ luminositySol: 100 }) > B.starMass({ luminositySol: 10 }));
  assert.equal(B.starMass({}), 0);
});

test("inhabited worlds share per type, populated ones by population before habitable ones by score", () => {
  const items = [world(1, "planet", 9), world(2, "planet", 4, 5e6), world(3, "planet", 6, 2e9), world(4, "station", 1, 800), world(5, "moon", 8)];
  const chosen = B.chooseLabels(items, { budget: 5 });
  // Types in turn (moon, planet, station), the best of each first.
  assert.deepEqual(chosen.slice(0, 3).map((c) => c.id).sort(), ["w3", "w4", "w5"]);
  assert.deepEqual(chosen.slice(3).map((c) => c.id), ["w2", "w1"]);
});

test("uninhabited bodies share per category by mass, a belt by the volume it occupies", () => {
  const items = [rock(1, "gas giant", 1e-3), rock(2, "gas giant", 5e-4), rock(3, "rocky world", 3e-6), rock(4, "rocky world", 1e-6),
    rock(5, "asteroid belt", undefined, 4e6), rock(6, "asteroid belt", undefined, 9e6), rock(7, "star-like", 2.0)];
  const chosen = B.chooseLabels(items, { budget: 4 });
  assert.deepEqual(chosen.map((c) => c.category).sort(), ["asteroid belt", "gas giant", "rocky world", "star-like"]);
  assert.equal(chosen.find((c) => c.category === "asteroid belt").id, "u6");
  assert.equal(chosen.find((c) => c.category === "gas giant").id, "u1");
});

test("nebulae are taken by volume and regions by share charted, ties by population", () => {
  const nebulae = B.chooseLabels([cloud(1, 5), cloud(2, 900), cloud(3, 40)], { budget: 3 });
  assert.deepEqual(nebulae.map((c) => c.id), ["n2", "n3", "n1"]);
  const regions = B.chooseLabels([region(1, 0.14, 0), region(2, 0.5, 0), region(3, 0.14, 1e6)], { budget: 3 });
  assert.deepEqual(regions.map((c) => c.id), ["r2", "r3", "r1"]);
});

test("a first item (a region at the farthest zoom) goes before everything", () => {
  const items = [star(1, 1, 1, 1), Object.assign(region(1, 0.1), { first: true }), world(1, "planet", 7)];
  const chosen = B.chooseLabels(items, { budget: 3 });
  assert.equal(chosen[0].id, "r1");
  assert.deepEqual(domainsOf(B.chooseLabels(items, { budget: 1 })), ["region"]);
});

test("the choice is the same for the same items in any order", () => {
  const items = [];
  for (let i = 0; i < 60; i++) items.push(star(i, 1 + (i * 7) % 13, 1 + (i * 5) % 17, 1 + (i * 3) % 11), world(i, i % 2 ? "planet" : "moon", i % 10), rock(i, "rocky world", (i % 9) / 100));
  const a = B.chooseLabels(items, { budget: 40 }).map((c) => c.id);
  const b = B.chooseLabels(items.slice().reverse(), { budget: 40 }).map((c) => c.id);
  assert.deepEqual(a, b);
});

test("an inhabited world is named before a habitable one, and its type is served first (MAP.321)", () => {
  const home = Object.assign(world(1, "inhabited world", 0.2, 0), { inhabited: true });
  const good = world(2, "habitable world", 1.4, 0);
  const better = world(3, "habitable world", 1.6, 0);
  const chosen = B.chooseLabels([good, better, home], { budget: 100 });
  assert.deepEqual(chosen.map((c) => c.id), ["w1", "w3", "w2"]);
});

test("inside a type an inhabited world with no score yet still goes before a habitable one with a high score", () => {
  const first = Object.assign(world(1, "planet", 0, 0), { inhabited: true });
  const chosen = B.chooseLabels([world(2, "planet", 1.9, 0), first], { budget: 100 });
  assert.deepEqual(chosen.map((c) => c.id), ["w1", "w2"]);
});

const sector = (id, massKg) => ({ id: "sec" + id, domain: "region", first: true, massKg });

test("sectors are named heaviest first by their total mass, and a short budget drops the lightest (MAP.322)", () => {
  const items = [sector(1, 3e30), sector(2, 9e32), sector(3, 0), sector(4, 5e31)];
  assert.deepEqual(B.chooseLabels(items, { budget: 100 }).map((c) => c.id), ["sec2", "sec4", "sec1", "sec3"]);
  assert.deepEqual(B.chooseLabels(items, { budget: 2 }).map((c) => c.id), ["sec2", "sec4"]);
});

test("things that go first with no mass stay ahead of the sectors in their own order, and equal masses keep the list's order", () => {
  const farRegion = { id: "r1", domain: "region", first: true, shareCharted: 0.1 };
  const items = [sector(1, 2e30), farRegion, sector(2, 2e30), { id: "r2", domain: "region", first: true }];
  assert.deepEqual(B.chooseLabels(items, { budget: 100 }).map((c) => c.id), ["r1", "r2", "sec1", "sec2"]);
});

test("a sparse view names every sector whatever its mass", () => {
  const items = [sector(1, 1e30), sector(2, 2e30), sector(3, 3e30)];
  assert.equal(B.chooseLabels(items, { budget: 50 }).length, 3);
});
