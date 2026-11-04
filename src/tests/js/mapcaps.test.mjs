// tests/js/mapcaps.test.mjs -- MAP.367: the caps on what the Galaxy Map draws at once are about half what they were, and a sector opened in place
// keeps the old numbers (everything is shown).

import { test } from "node:test";
import assert from "node:assert/strict";

import { installDom, importPage } from "./fakedom.mjs";

installDom();
const MAG = await importPage("starmagnitude.js");
const MARKERS = await importPage("galaxymarkers.js");
const NEBULAE = await importPage("galaxynebulae.js");
const BUDGET = await importPage("labelbudget.js");
const CUT = await importPage("galaxynodecut.js");

// What each cap was before MAP.367; "about half" is 0.45 to 0.55 of it.
const WAS = { stars: 20000, starsPhone: 8000, markers: 240, regions: 600, named: 400, nodes: 1200 };
const about = (now, was) => now >= 0.45 * was && now <= 0.55 * was;

test("the stars held on the screen are about half what they were", () => {
  assert.ok(about(MAG.TARGET_STARS, WAS.stars), String(MAG.TARGET_STARS));
  assert.ok(about(MAG.TARGET_STARS_PHONE, WAS.starsPhone), String(MAG.TARGET_STARS_PHONE));
});

test("an opened sector keeps the old star numbers", () => {
  assert.equal(MAG.TARGET_STARS_OPENED, WAS.stars);
  assert.equal(MAG.TARGET_STARS_OPENED_PHONE, WAS.starsPhone);
});

test("the marks drawn at once are about half what they were, and an opened sector's budget is the larger", () => {
  assert.ok(about(MARKERS.MARKER_BUDGET, WAS.markers), String(MARKERS.MARKER_BUDGET));
  assert.ok(MARKERS.SECTOR_MARKER_BUDGET >= 10 * MARKERS.MARKER_BUDGET);
});

test("the nebula regions drawn at once are about half what they were", () => {
  assert.ok(about(NEBULAE.REGION_BUDGET, WAS.regions), String(NEBULAE.REGION_BUDGET));
});

test("the named stars and objects offered to the label engine are about half what they were", () => {
  assert.ok(about(BUDGET.NAMED_LIMIT, WAS.named), String(BUDGET.NAMED_LIMIT));
});

test("the nodes a cut holds stay at 1200: halving them (tried) made the blocks hand over in bursts, see starfade.test.mjs", () => {
  assert.equal(CUT.MAX_NODES, WAS.nodes);
});

test("the marker selection holds to the budget it is given, and the default is the smaller one", () => {
  const kinds = MARKERS.createMarkerKinds();
  kinds.register({ id: "marker-test", label: "T", glyph: "black-hole", color: "#fff", maxBlock: Infinity, defaultOn: true, source: (node) => node.items, rank: (item) => item.r });
  const nodes = [{ key: "1.0.0.0", block: { m: 1 }, node: { items: Array.from({ length: 1000 }, (_, i) => ({ key: String(i), x: i, y: 0, z: 0, r: i })) } }];
  const chosen = (budget) => MARKERS.selectMarkers({ kinds: kinds.kinds(), isOn: () => true, nodes, budget }).length;
  assert.ok(chosen(undefined) <= MARKERS.MARKER_BUDGET);
  assert.ok(chosen(MARKERS.SECTOR_MARKER_BUDGET) > MARKERS.MARKER_BUDGET);
});
