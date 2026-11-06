// tests/js/viewmenu.test.mjs -- static/viewmenu.js (UX.278): the View menu's shared builder.

import { test } from "node:test";
import assert from "node:assert/strict";

import { h, installDom, FakeEvent, STATIC_URL } from "./fakedom.mjs";

installDom("http://localhost/system");
const L = await import(new URL("maplayers.js", STATIC_URL).href);
const V = await import(new URL("viewmenu.js", STATIC_URL).href);

function storage() {
  const kept = new Map();
  return {
    kept,
    getItem: (key) => (kept.has(key) ? kept.get(key) : null),
    setItem: (key, value) => { kept.set(key, String(value)); },
  };
}

function registry() {
  installDom("http://localhost/system");
  const reg = L.createLayerRegistry("galaxy");
  const base = { modes: ["galaxy"], urlKey: "hide", urlKeyOn: "show" };
  reg.register({ ...base, id: "star", label: "Stars", group: "Show on the map", heading: "Show", keywords: ["star"], common: true });
  reg.register({ ...base, id: "nebula", label: "Nebulae", group: "Show on the map", heading: "Show", keywords: ["nebula"], data: { kind: "nebula" } });
  reg.register({ ...base, id: "rogue", label: "Rogue planets", group: "Show on the map", heading: "Show", keywords: ["rogue"], defaultOn: false });
  reg.register({ ...base, id: "grid", label: "Cell grid", heading: "Color and fill", keywords: ["grid"], defaultOn: false });
  reg.register({ ...base, id: "spin", label: "Spin", heading: "Move and Settings", keywords: ["spin"], data: { setting: "spin" } });
  return reg;
}

function change(box, checked) {
  box.checked = checked;
  box.dispatchEvent(new FakeEvent("change"));
}

test("the plan has a heading in the registry's order for each heading that has a layer, and groups in first-seen order", () => {
  const reg = registry();
  const plan = V.menuPlan(reg.layers(), reg.headings);
  assert.deepEqual(plan.map((p) => p.heading), ["Show", "Color and fill", "Move and Settings"], "Chart and side views has nothing, so it is left out");
  assert.deepEqual(plan[0].groups.map((g) => [g.name, g.layers.map((l) => l.id)]), [["Show on the map", ["star", "nebula", "rogue"]]]);
  assert.deepEqual(plan[1].groups.map((g) => g.name), [null], "a layer with no group sits in a group with no name");
  const odd = V.menuPlan([{ id: "x", label: "X", heading: "Nowhere" }], reg.headings);
  assert.deepEqual(odd.map((p) => p.heading), ["Show"], "a layer that names no known heading goes under the first");
});

test("a group is all, some or none, and a count says how many are on and how many differ from the start", () => {
  const reg = registry();
  const layers = reg.layers().slice(0, 3);
  assert.equal(V.groupState(layers, reg.isOn), "some");
  reg.setOn("rogue", true);
  assert.equal(V.groupState(layers, reg.isOn), "all");
  reg.setGroup("Show on the map", "hidden");
  assert.equal(V.groupState(layers, reg.isOn), "none");
  assert.deepEqual(V.tally(layers, reg.isOn, reg.isPulled), { on: 0, total: 3, changed: 2 }, "stars and nebulae are off, rogue planets are as they start");
  assert.equal(V.countText({ on: 4, total: 12, changed: 0 }), "4 of 12 on");
  assert.equal(V.countText({ on: 4, total: 12, changed: 2 }), "4 of 12 on, 2 changed");
});

test("the menu keeps the layer selectors, shows a count per heading, and a checkbox switches the layer", () => {
  const reg = registry();
  const container = h("div", {});
  V.renderViewMenu(reg, container, { storage: storage() });
  assert.deepEqual(container.querySelectorAll("details").map((d) => d.dataset.heading), ["Show", "Color and fill", "Move and Settings"]);
  const star = container.querySelector('input[data-layer="star"]');
  assert.equal(star.type, "checkbox");
  assert.equal(container.querySelector('input[data-layer="nebula"]').dataset.kind, "nebula", "the layer's data-* stay on the control");
  assert.equal(container.querySelector('input[data-layer="spin"]').dataset.setting, "spin");
  assert.deepEqual(container.querySelectorAll("input[data-layer]").map((b) => b.checked), [true, true, false, false, true]);
  assert.equal(container.querySelectorAll("summary .view-menu-count")[0].textContent, "2 of 3 on");
  change(star, false);
  assert.equal(reg.isOn("star"), false);
  assert.equal(container.querySelectorAll("summary .view-menu-count")[0].textContent, "1 of 3 on, 1 changed");
  reg.setOn("grid", true);
  assert.equal(container.querySelector('input[data-layer="grid"]').checked, true, "a change made elsewhere shows in the box");
});

test("a layer that is not as it starts is marked in words", () => {
  const reg = registry();
  const container = h("div", {});
  V.renderViewMenu(reg, container, { storage: storage() });
  const mark = (id) => container.querySelector('[data-row="' + id + '"] .view-menu-changed');
  assert.equal(mark("star").hidden, true);
  reg.setOn("star", false);
  assert.equal(mark("star").hidden, false);
  assert.equal(mark("star").textContent, "changed");
  reg.setOn("star", true);
  assert.equal(mark("star").hidden, true);
});

test("a group of two or more has a box that is checked, mixed or clear, and pressing it sets them all", () => {
  const reg = registry();
  const container = h("div", {});
  V.renderViewMenu(reg, container, { storage: storage() });
  const box = container.querySelector('input[data-group-box="Show on the map"]');
  assert.equal(box.checked, false);
  assert.equal(box.indeterminate, true, "some are on");
  change(box, true);
  assert.deepEqual(["star", "nebula", "rogue"].map((id) => reg.isOn(id)), [true, true, true]);
  assert.equal(box.indeterminate, false);
  change(box, false);
  assert.deepEqual(["star", "nebula", "rogue"].map((id) => reg.isOn(id)), [false, false, false]);
  assert.equal(container.querySelectorAll("input[data-group-box]").length, 1, "a layer with no group has no box");
});

test("a small menu starts open; a long one folds all but the first heading, and the choice is kept", () => {
  const reg = registry();
  const small = h("div", {});
  const kept = storage();
  V.renderViewMenu(reg, small, { storage: kept });
  assert.deepEqual(small.querySelectorAll("details").map((d) => d.open), [true, true, true], "five layers: nothing to fold away");
  const long = h("div", {});
  V.renderViewMenu(reg, long, { storage: kept, collapseAbove: 2 });
  assert.deepEqual(long.querySelectorAll("details").map((d) => d.open), [true, false, false]);
  const second = long.querySelectorAll("details")[1];
  second.open = true;
  second.dispatchEvent(new FakeEvent("toggle", { bubbles: false }));
  assert.equal(kept.kept.get("planetgen.viewMenu.open.galaxy.Color and fill"), "1");
  const again = h("div", {});
  V.renderViewMenu(reg, again, { storage: kept, collapseAbove: 2 });
  assert.deepEqual(again.querySelectorAll("details").map((d) => d.open), [true, true, false], "the next page remembers it");
});

test("a menu with no storage still works", () => {
  const reg = registry();
  const container = h("div", {});
  V.renderViewMenu(reg, container, { storage: null });
  assert.equal(container.querySelectorAll("details").length, 3);
  const broken = { getItem() { throw new Error("blocked"); }, setItem() { throw new Error("blocked"); } };
  const other = h("div", {});
  V.renderViewMenu(reg, other, { storage: broken, collapseAbove: 1 });
  const second = other.querySelectorAll("details")[1];
  second.open = true;
  second.dispatchEvent(new FakeEvent("toggle", { bubbles: false }));
  assert.equal(other.querySelectorAll("details").length, 3);
});

test("a layer that is not available has no row until it is, and a headed layer's settings hook draws under its row", () => {
  installDom("http://localhost/system");
  const reg = L.createLayerRegistry("galaxy");
  let there = false;
  let drawn = 0;
  reg.register({ id: "star", label: "Stars", heading: "Show", keywords: ["star"] });
  reg.register({ id: "nebula", label: "Nebulae", heading: "Show", keywords: ["nebula"], available: () => there });
  reg.register({ id: "grid", label: "Cell grid", heading: "Color and fill", keywords: ["grid"],
    settings: (place) => { drawn += 1; place.textContent = "opacity slider"; } });
  const container = h("div", {});
  V.renderViewMenu(reg, container, { storage: storage() });
  assert.deepEqual(container.querySelectorAll("input[data-layer]").map((b) => b.dataset.layer), ["star", "grid"]);
  assert.equal(container.querySelector('[data-settings="grid"]').textContent, "opacity slider");
  assert.equal(drawn, 1);
  there = true;
  reg.refresh();
  assert.deepEqual(container.querySelectorAll("input[data-layer]").map((b) => b.dataset.layer), ["star", "nebula", "grid"]);
});

test("a forward layer has a Forward checkbox that pulls it forward and is greyed while the layer is off", () => {
  installDom("http://localhost/system");
  const reg = L.createLayerRegistry("galaxy");
  reg.register({ id: "marker-a", label: "Habitable worlds", heading: "Show", keywords: ["habitable"], forward: true });
  const container = h("div", {});
  V.renderViewMenu(reg, container, { storage: storage() });
  const front = container.querySelector('input[data-forward="marker-a"]');
  assert.equal(front.disabled, false);
  change(front, true);
  assert.equal(reg.isPulled("marker-a"), true);
  assert.equal(container.querySelector('[data-row="marker-a"] .view-menu-changed').hidden, false, "pulled forward is not as it starts");
  reg.setOn("marker-a", false);
  assert.equal(front.checked, false);
  assert.equal(front.disabled, true);
  assert.equal(front.getAttribute("aria-label"), "Pull habitable worlds forward");
});
