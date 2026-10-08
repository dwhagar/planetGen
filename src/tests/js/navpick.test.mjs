// tests/js/navpick.test.mjs -- static/navpick.js, the NAV course pick the
// maps share (NAV.29, NAV.33): the Start Here and End Here buttons, the pick
// moving on in place after the first end, the URL query and the banner words.

import { test } from "node:test";
import assert from "node:assert/strict";

import { STATIC_URL } from "./fakedom.mjs";

const { createNavPick } = await import(new URL("navpick.js", STATIC_URL).href);

test("with no pick an object offers Start Here and End Here", () => {
  const pick = createNavPick({ navUrl: "/nav" });
  assert.equal(pick.active(), false);
  assert.equal(pick.query(), "");
  assert.deepEqual(pick.actionsFor("system:1", "Sol").map((a) => a.label), ["Start Here", "End Here"]);
  assert.deepEqual(pick.actionsFor(null, "Nothing"), []);
});

test("Start Here keeps the user on the map and asks for the destination", () => {
  let changes = 0;
  const pick = createNavPick({ navUrl: "/nav", onChange: () => { changes += 1; } });
  pick.actionsFor("system:1", "Sol")[0].onClick();
  assert.equal(changes, 1);
  assert.equal(pick.pick(), "to");
  assert.equal(pick.query(), "?pick=to&from=system:1");
  const next = pick.actionsFor("nebula:3", "Veil");
  assert.equal(next.length, 1);
  assert.equal(next[0].label, "End Here");
  assert.equal(next[0].href, "/nav?from=system:1&to=nebula:3");
  assert.deepEqual(pick.bannerParts("star")[0], "Start: Sol. Choosing a destination");
  assert.equal(pick.cancelUrl(), null, "a pick begun on the map is cancelled in place");
});

test("End Here asks for the start, and the course keeps the order from, to", () => {
  const pick = createNavPick({ navUrl: "/nav" });
  pick.actionsFor("system:2", "Vega")[1].onClick();
  assert.equal(pick.pick(), "from");
  assert.equal(pick.query(), "?pick=from&to=system:2");
  assert.equal(pick.actionsFor("system:1", "Sol")[0].href, "/nav?from=system:1&to=system:2");
});

test("a pick the NAV page began with one end set goes straight back to NAV", () => {
  const pick = createNavPick({ navUrl: "/nav", pick: "to", other: "system:12", cancelUrl: "/nav?from=system:12" });
  assert.equal(pick.query(), "?pick=to&from=system:12");
  assert.deepEqual(pick.actionsFor("system:5", "Alpha"), [
    { label: "End Here", href: "/nav?from=system:12&to=system:5", primary: true },
  ]);
  assert.equal(pick.cancelUrl(), "/nav?from=system:12");
  assert.equal(pick.bannerParts("system")[0], "Choosing a destination");
});

test("a pick the NAV page began with no end set moves on in place", () => {
  const pick = createNavPick({ navUrl: "/nav", pick: "from" });
  assert.equal(pick.query(), "?pick=from");
  const only = pick.actionsFor("system:7", "Rigel");
  assert.equal(only.length, 1);
  assert.equal(only[0].label, "Start Here");
  only[0].onClick();
  assert.equal(pick.query(), "?pick=to&from=system:7");
});

test("clearing ends the pick", () => {
  let changes = 0;
  const pick = createNavPick({ navUrl: "/nav", onChange: () => { changes += 1; } });
  pick.begin("from", "system:1", "Sol");
  pick.clear();
  assert.equal(changes, 2);
  assert.equal(pick.active(), false);
  assert.equal(pick.query(), "");
  assert.equal(pick.bannerParts(), null);
});
