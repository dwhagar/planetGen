// tests/js/bodymarkers.test.mjs -- static/bodymarkers.js (MAP.321): the symbols of the bodies of the open system.

import test from "node:test";
import assert from "node:assert/strict";

const B = await import("../../html/static/bodymarkers.js");
const M = await import("../../html/static/galaxymarkers.js");

test("a body's reason follows the order of the stars' marks: inhabited, gas giant, habitable, barren, belt", () => {
  const planet = (extra) => Object.assign({ kind: "planet", body_type: "t", inhabited: false, tier: null }, extra);
  assert.deepEqual(B.reasonOfBody(planet({ inhabited: true, tier: "Green" })), { reason: "inhabited", tier: 1 });
  assert.deepEqual(B.reasonOfBody(planet({ body_type: "g" })), { reason: "gas-giant", tier: null });
  assert.deepEqual(B.reasonOfBody(planet({ tier: "Blue" })), { reason: "habitable", tier: 0 });
  assert.deepEqual(B.reasonOfBody(planet({ tier: "Yellow" })), { reason: "habitable", tier: 2 });
  assert.deepEqual(B.reasonOfBody(planet({ tier: "Red" })), { reason: "uninhabitable", tier: 3 });
  assert.equal(B.reasonOfBody(planet({})), null, "an unscored terrestrial body has no symbol");
  assert.deepEqual(B.reasonOfBody({ kind: "moon", body_type: "t", tier: "Green" }), { reason: "habitable", tier: 1 });
  assert.deepEqual(B.reasonOfBody({ kind: "belt" }), { reason: "asteroid-field", tier: null });
  assert.equal(B.reasonOfBody({ kind: "star" }), null);
  assert.equal(B.reasonOfBody({ kind: "comet" }), null);
  assert.equal(B.reasonOfBody(null), null);
});

test("a body has every reason it qualifies for, strongest first, from the stored life scores and class as they are (MAP.343)", () => {
  const planet = (extra) => Object.assign({ kind: "planet", body_type: "t", inhabited: false, tier: null }, extra);
  const reasons = (body) => B.reasonsOfBody(body).map((one) => one.reason);
  assert.deepEqual(reasons(planet({ life_microbial: 1 })), ["life"]);
  assert.deepEqual(reasons(planet({ life_microbial: 0, life_flora: 0, life_fauna: null, life_sentient: null })), [], "a zero score is none");
  assert.deepEqual(reasons(planet({ life_sentient: 2 })), ["sentient", "life"]);
  assert.deepEqual(reasons(planet({ life_stage: "technological_civilization", life_sentient: 4, inhabited: true, tier: "Green", planet_class: "M" })),
    ["inhabited", "civilization", "sentient", "life", "habitable", "class-M"]);
  assert.deepEqual(B.reasonsOfBody(planet({ inhabited: true, tier: "Green" }))[0], { reason: "inhabited", tier: 1 });
  assert.deepEqual(B.reasonsOfBody(planet({ life_microbial: 2, tier: "Green" })).map((one) => one.tier), [null, 1], "only the world reasons carry the habitability tier");
  assert.deepEqual(reasons(planet({ body_type: "g", planet_class: "J" })), ["gas-giant", "class-J"]);
  assert.deepEqual(reasons(planet({ planet_class: "R" })), [], "a letter that is not a class has no symbol");
  assert.deepEqual(reasons({ kind: "moon", body_type: "t", tier: "Red", planet_class: "D", life_microbial: 1 }), ["life", "uninhabitable", "class-D"]);
  assert.deepEqual(reasons({ kind: "star", life_microbial: 9 }), []);
  assert.deepEqual(reasons(null), []);
  assert.deepEqual(B.reasonOfBody(planet({ life_flora: 1 })), { reason: "life", tier: null });
  assert.deepEqual(reasons(planet({ life_microbial: 3, planet_class: "M", life_class_id: "4b" })), ["life", "class-M", "lifeclass-4b"], "the life class that lives on a world is its last reason (MAP.351)");
  assert.deepEqual(reasons(planet({ life_class_id: "zz" })), [], "an id that is not a life class has no symbol");
});

test("the symbol of a body is its strongest reason whose layer is on, and a layer that is off lets the next one name it (MAP.343)", () => {
  const body = { body_type: "t", tier: "Green", inhabited: true, life_microbial: 2, life_sentient: 1, life_stage: "technological_civilization", planet_class: "M" };
  const home = { ref: "planet:1", kind: "planet", name: "Home", x: 0, y: 0, z: 0, body };
  const on = (...ids) => B.bodyMarkers(M.INTEREST_KINDS, (id) => ids.indexOf(id) >= 0, [home]).map((c) => [c.kind.id, c.item.tier]);
  assert.deepEqual(on("marker-inhabited", "marker-life"), [["marker-inhabited", 1]]);
  assert.deepEqual(on("marker-life", "marker-habitable"), [["marker-life", undefined]]);
  assert.deepEqual(on("marker-habitable"), [["marker-habitable", 1]], "the inhabited layer being off, the habitable one names it");
  assert.deepEqual(on("marker-sentient", "marker-civilization"), [["marker-civilization", undefined]]);
  assert.deepEqual(on("marker-class-M"), [["marker-class-M", undefined]]);
  assert.deepEqual(on("marker-class-K"), []);
  assert.equal(B.bodyMarkers(M.INTEREST_KINDS, (id) => id === "marker-class-M", [home])[0].kind.describe({ name: "Home" }), "Home, class M world");
});

test("the markers for the bodies follow the bodies and only the layers that are on draw", () => {
  const earth = { ref: "planet:1", kind: "planet", name: "Earthlike", x: 1, y: 2, z: 3, body: { body_type: "t", tier: "Green", inhabited: true } };
  const rock = { ref: "planet:2", kind: "planet", name: "Rock", x: 4, y: 5, z: 6, body: { body_type: "t", tier: "Red" } };
  const belt = { ref: "belt:3", kind: "belt", name: "Belt", x: 7, y: 8, z: 9, body: {} };
  const star = { ref: "star:4", kind: "star", name: "Sun", x: 0, y: 0, z: 0, body: {} };
  const entries = [earth, rock, belt, star];
  const all = B.bodyMarkers(M.INTEREST_KINDS, () => true, entries);
  assert.deepEqual(all.map((c) => [c.item.name, c.kind.id]), [
    ["Earthlike", "marker-inhabited"], ["Rock", "marker-uninhabitable"], ["Belt", "marker-asteroid-field"]]);
  assert.equal(all[0].item.follow, earth, "the symbol follows the body's entry");
  assert.equal(all[0].item.bodyRef, "planet:1");
  assert.equal(all[0].item.tier, 1);
  assert.equal(all[1].item.tier, 3);
  assert.equal(all[2].item.tier, undefined);
  assert.equal(all[0].kind.describe(all[0].item), "Earthlike, inhabited world, Green");
  assert.equal(all[2].kind.describe(all[2].item), "Belt, asteroid belt");
  earth.x = 10;
  assert.equal(all[0].item.follow.x, 10);
  const onlyBelts = B.bodyMarkers(M.INTEREST_KINDS, (id) => id === "marker-asteroid-field", entries);
  assert.deepEqual(onlyBelts.map((c) => c.item.name), ["Belt"]);
  assert.deepEqual(B.bodyMarkers(M.INTEREST_KINDS, () => false, entries), []);
  assert.deepEqual(B.bodyMarkers(M.INTEREST_KINDS, () => true, null), []);
});

test("the hand-off has a gap, so a camera resting between the two distances does not flicker", () => {
  const handoff = B.createHandoff();
  const radius = 10;
  assert.equal(handoff.update(10 * B.HANDOFF_OFF, radius, true), false, "far out the stars carry the symbols");
  assert.equal(handoff.update(10 * (B.HANDOFF_ON + B.HANDOFF_OFF) / 2, radius, true), false, "in the gap it stays as it was");
  assert.equal(handoff.update(10 * B.HANDOFF_ON, radius, true), true, "close in the bodies carry them");
  assert.equal(handoff.update(10 * (B.HANDOFF_ON + B.HANDOFF_OFF) / 2, radius, true), true, "and in the gap it stays so");
  assert.equal(handoff.update(10 * B.HANDOFF_OFF, radius, true), false, "beyond the far limit it hands back");
  assert.equal(handoff.update(10 * B.HANDOFF_ON, radius, true), true);
  assert.equal(handoff.update(0, radius, false), false, "a system that is closed hands back at once");
  assert.equal(handoff.isOn(), false);
  assert.equal(handoff.update(1, 0, true), false, "a system with no size yet is not close");
  assert.ok(B.HANDOFF_ON < B.HANDOFF_OFF);
});
