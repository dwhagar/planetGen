// tests/js/maplayers.test.mjs -- static/maplayers.js (MAP.217): the layer registry the three maps share.

import { test } from "node:test";
import assert from "node:assert/strict";

import { h, installDom, STATIC_URL } from "./fakedom.mjs";

installDom("http://localhost/galaxy");
const L = await import(new URL("maplayers.js", STATIC_URL).href);

function registry(mode) {
  installDom("http://localhost/galaxy");
  const reg = L.createLayerRegistry(mode || "galaxy");
  reg.register({ id: "star", label: "Stars", modes: ["galaxy", "sector"], group: "Show on the map",
    urlKey: "hide", urlKeyOn: "show", legend: "Stars are points of light" });
  reg.register({ id: "rogue", label: "Rogue planets", modes: ["galaxy", "sector"], group: "Show on the map",
    defaultOn: false, urlKey: "hide", urlKeyOn: "show" });
  reg.register({ id: "course", label: "Course", modes: ["galaxy", "sector", "system"], fadedWithSectors: false,
    urlKey: "hide", urlKeyOn: "show", legend: { text: "The plotted course", swatch: "#f90" } });
  reg.register({ id: "facility", label: "Facilities", modes: ["system"], urlKey: "hide", urlKeyOn: "show" });
  return reg;
}

test("a layer needs an id and a label, a known mode and urlKeyOn when it begins off", () => {
  const reg = L.createLayerRegistry("galaxy");
  assert.throws(() => reg.register({ id: "x" }), /id and a label/);
  assert.throws(() => reg.register({ id: "x", label: "X", modes: ["moon"] }), /unknown map mode/);
  assert.throws(() => reg.register({ id: "x", label: "X", defaultOn: false, urlKey: "hide" }), /urlKeyOn/);
  reg.register({ id: "x", label: "X" });
  assert.throws(() => reg.register({ id: "x", label: "X" }), /already registered/);
  assert.throws(() => L.createLayerRegistry("moon"), /unknown map mode/);
});

test("a registry holds only the layers of its mode", () => {
  assert.deepEqual(registry("galaxy").layers().map((l) => l.id), ["star", "rogue", "course"]);
  assert.deepEqual(registry("system").layers().map((l) => l.id), ["course", "facility"]);
});

test("a layer is on by default unless it says otherwise, and can be set", () => {
  const reg = registry();
  assert.equal(reg.isOn("star"), true);
  assert.equal(reg.isOn("rogue"), false);
  const seen = [];
  reg.onChange((id, on) => seen.push([id, on]));
  assert.equal(reg.setOn("star", false), false);
  assert.equal(reg.toggle("rogue"), true);
  reg.setOn("rogue", true);
  assert.deepEqual(seen, [["star", false], ["rogue", true]], "setting what it already is does not announce");
});

test("the URL lists the layers turned off under urlKey and those turned on under urlKeyOn", () => {
  const reg = registry();
  assert.equal(reg.writeUrl(new URLSearchParams()).toString(), "");
  reg.setOn("star", false);
  reg.setOn("course", false);
  reg.setOn("rogue", true);
  const params = reg.writeUrl(new URLSearchParams("sector=5&hide=old"));
  assert.equal(params.get("hide"), "star,course");
  assert.equal(params.get("show"), "rogue");
  assert.equal(params.get("sector"), "5", "other parameters stay");
  reg.setOn("star", true);
  reg.setOn("course", true);
  reg.setOn("rogue", false);
  assert.equal(reg.writeUrl(params).toString(), "sector=5");
});

test("the URL round trip, including a kind registered after the URL was read", () => {
  const reg = registry();
  reg.readUrl("hide=star,black-hole&show=rogue");
  assert.equal(reg.isOn("star"), false);
  assert.equal(reg.isOn("rogue"), true);
  const written = reg.writeUrl(new URLSearchParams()).toString().replace(/%2C/gi, ",");
  assert.equal(written, "hide=star%2Cblack-hole&show=rogue".replace(/%2C/gi, ","));
  reg.register({ id: "black-hole", label: "Black holes", modes: ["galaxy"], urlKey: "hide", urlKeyOn: "show" });
  assert.equal(reg.isOn("black-hole"), false, "the hiding the URL asked for applies when the layer comes");
  reg.unregister("black-hole");
  assert.equal(reg.writeUrl(new URLSearchParams()).get("hide"), "star,black-hole", "a sector without it keeps the setting");
});

test("a registry that forgets unknown ids drops the ones no layer is registered for", () => {
  const reg = L.createLayerRegistry("galaxy", { forgetUnknown: true });
  reg.register({ id: "star", label: "Stars", urlKey: "hide" });
  reg.readUrl("hide=star,nonsense");
  assert.equal(reg.isOn("star"), false);
  assert.equal(reg.writeUrl(new URLSearchParams()).get("hide"), "star");
});

test("reading a URL puts layers it does not name back to their defaults and announces the changes", () => {
  const reg = registry();
  reg.setOn("star", false);
  const seen = [];
  reg.onChange((id) => seen.push(id));
  reg.readUrl("");
  assert.equal(reg.isOn("star"), true);
  assert.deepEqual(seen, ["star"]);
});

test("the legend lists the layers that are on and have one, with their swatch", () => {
  const reg = registry();
  assert.deepEqual(reg.legend(), [
    { id: "star", text: "Stars are points of light", swatch: null, swatchClass: null },
    { id: "course", text: "The plotted course", swatch: "#f90", swatchClass: null },
  ]);
  reg.setOn("star", false);
  assert.deepEqual(reg.legend().map((line) => line.id), ["course"]);
  const withFunction = L.createLayerRegistry("system");
  withFunction.register({ id: "f", label: "F", legend: (on) => (on ? "on text" : "off text") });
  assert.equal(withFunction.legend()[0].text, "on text");
});

test("the picker gets the declarations of the layers that are on", () => {
  const reg = L.createLayerRegistry("sector");
  const picked = { points: () => [] };
  reg.register({ id: "a", label: "A", pickerLayer: picked });
  reg.register({ id: "b", label: "B" });
  assert.deepEqual(reg.pickerLayers(), [picked]);
  reg.setOn("a", false);
  assert.deepEqual(reg.pickerLayers(), []);
});

test("the blocker fade hides a layer with the sectors unless the layer says it does not", () => {
  const reg = registry();
  assert.equal(reg.visibleShare("star", 0.25), 0.25);
  assert.equal(reg.visibleShare("star", 7), 1);
  assert.equal(reg.visibleShare("course", 0), 1, "a course layer is not faded");
  reg.setOn("course", false);
  assert.equal(reg.visibleShare("course", 1), 0, "an off layer shows none");
  assert.equal(reg.visibleShare("nothing", 1), 0);
});

test("draw calls every layer of the mode with the scene and whether it is on", () => {
  const reg = L.createLayerRegistry("sector");
  const calls = [];
  reg.register({ id: "a", label: "A", draw: (scene, on) => calls.push(["a", scene, on]) });
  reg.register({ id: "b", label: "B", modes: ["galaxy"], draw: () => calls.push(["b"]) });
  reg.setOn("a", false);
  reg.draw("scene");
  assert.deepEqual(calls, [["a", "scene", false]]);
});

test("the menu has a heading per group and a pressed button per layer, and a click toggles the layer", () => {
  const reg = registry();
  const container = h("div", {});
  reg.renderMenu(container);
  const buttons = container.querySelectorAll("button");
  assert.deepEqual(buttons.map((b) => b.dataset.layer), ["star", "rogue", "course"]);
  assert.deepEqual(buttons.map((b) => b.getAttribute("aria-pressed")), ["true", "false", "true"]);
  assert.equal(container.querySelectorAll(".galaxy-kinds-heading").length, 1);
  buttons[0].click();
  assert.equal(reg.isOn("star"), false);
  assert.equal(buttons[0].getAttribute("aria-pressed"), "false");
  reg.setOn("rogue", true);
  assert.equal(buttons[1].getAttribute("aria-pressed"), "true", "a change made elsewhere shows on the button");
});

test("a menu drawn again in the same container replaces the first", () => {
  const reg = registry();
  const container = h("div", {});
  reg.renderMenu(container);
  reg.renderMenu(container);
  assert.equal(container.querySelectorAll("button").length, 3);
});

test("the legend block shows a line per entry, follows the layers and hides when empty", () => {
  const reg = registry("system");
  reg.register({ id: "swatchy", label: "S", legend: { text: "a swatch", swatch: "red" } });
  const container = h("div", {});
  reg.renderLegend(container);
  assert.equal(container.hidden, false);
  assert.equal(container.querySelectorAll(".map-layer-legend-line").length, 2);
  reg.setOn("course", false);
  assert.equal(container.querySelectorAll(".map-layer-legend-line").length, 1);
  reg.setOn("swatchy", false);
  assert.equal(container.hidden, true);
});

test("a layer can be left out of the menu while it has nothing to show, and its title can change", () => {
  const reg = L.createLayerRegistry("sector");
  let present = false;
  let count = 3;
  reg.register({ id: "nebula", label: "Nebulae", available: () => present, title: () => "Show the nebulae (" + count + ")",
    data: { kind: "nebula" }, urlKey: "hide" });
  reg.register({ id: "star", label: "Stars" });
  const container = h("div", {});
  reg.renderMenu(container);
  assert.deepEqual(container.querySelectorAll("button").map((b) => b.dataset.layer), ["star"]);
  reg.setOn("nebula", false);
  assert.equal(reg.writeUrl(new URLSearchParams()).get("hide"), "nebula", "an absent layer keeps its URL entry");
  present = true;
  count = 5;
  reg.refresh();
  const button = container.querySelector('button[data-layer="nebula"]');
  assert.equal(button.dataset.kind, "nebula");
  assert.equal(button.title, "Show the nebulae (5)");
  assert.equal(button.getAttribute("aria-pressed"), "false");
  assert.deepEqual(container.querySelectorAll("button").map((b) => b.dataset.layer), ["nebula", "star"]);
});

test("a legend swatch can be a CSS class", () => {
  const reg = L.createLayerRegistry("system");
  reg.register({ id: "facility", label: "Facilities", legend: { text: "marks a facility", swatchClass: "sysmap-legend-facility" } });
  const container = h("div", {});
  reg.renderLegend(container);
  const swatch = container.querySelector(".map-layer-legend-swatch");
  assert.ok(swatch.classList.contains("sysmap-legend-facility"));
});
