// tests/js/mappick.test.mjs -- static/mappick.js, the picking, hover and
// info-panel layer the Galaxy Map and the Sector Map share (MAP.65): which
// layer wins under the pointer, the tooltip kept inside the map, and the
// info panel's rows, NAV links, pick button and ☆.

import { test } from "node:test";
import assert from "node:assert/strict";

import { h, installDom, STATIC_URL } from "./fakedom.mjs";

installDom("http://localhost/sector/1");
const THREE = await import(new URL("vendor/three.module.min.js", STATIC_URL).href);
const MP = await import(new URL("mappick.js", STATIC_URL).href);

const SIZE = 400;

// A camera 100 units up the z axis looking at the origin, on a 400 px
// canvas.
function scene() {
  installDom("http://localhost/sector/1");
  const canvasEl = h("canvas", {});
  canvasEl.rect = { left: 0, top: 0, width: SIZE, height: SIZE };
  canvasEl.clientWidth = SIZE;
  canvasEl.clientHeight = SIZE;
  const camera = new THREE.PerspectiveCamera(45, 1, 0.1, 1000);
  camera.position.set(0, 0, 100);
  camera.lookAt(0, 0, 0);
  camera.updateMatrixWorld();
  return { canvasEl, camera, picker: MP.createPicker(camera, canvasEl) };
}

// A square plate facing the camera at height z.
function plate(z, size) {
  const mesh = new THREE.Mesh(new THREE.PlaneGeometry(size || 40, size || 40), new THREE.MeshBasicMaterial());
  mesh.position.set(0, 0, z);
  mesh.updateMatrixWorld();
  return mesh;
}

const CENTER = SIZE / 2;

test("a point of light beats the see-through layer it sits in", () => {
  const s = scene();
  const point = { x: 0, y: 0, z: 0, name: "star" };
  const block = plate(10);
  s.picker.addLayer({ name: "blocks", meshes: () => [block], entryOf: () => "block" });
  s.picker.addLayer({ name: "points", points: () => [point], reach: () => 8 });
  // Added after, but the blocks layer was first: it wins.
  assert.equal(s.picker.pick(CENTER, CENTER).layer.name, "blocks");
  const t = scene();
  t.picker.addLayer({ name: "points", points: () => [point], reach: () => 8 });
  t.picker.addLayer({ name: "blocks", meshes: () => [block], entryOf: () => "block" });
  const found = t.picker.pick(CENTER, CENTER);
  assert.equal(found.layer.name, "points");
  assert.equal(found.entry, point);
  assert.equal(t.picker.pick(CENTER + 30, CENTER).layer.name, "blocks", "off the point, the block");
});

test("priority puts a layer after the others whatever order they came in", () => {
  const s = scene();
  const block = plate(10);
  s.picker.addLayer({ name: "blocks", priority: 10, meshes: () => [block], entryOf: () => 0 });
  s.picker.addLayer({ name: "points", points: () => [{ x: 0, y: 0, z: 0 }], reach: () => 8 });
  assert.equal(s.picker.pick(CENTER, CENTER).layer.name, "points");
  assert.equal(s.picker.pick(CENTER + 30, CENTER).entry, 0, "an entry of 0 is still an entry");
});

test("a solid layer in front hides a point behind it, but not one in front", () => {
  const s = scene();
  const behind = { x: 0, y: 0, z: -20 };
  const body = plate(0, 10);
  s.picker.addLayer({ name: "points", points: () => [behind], reach: () => 8 });
  s.picker.addLayer({ name: "bodies", occludes: true, meshes: () => [body], entryOf: () => "body" });
  assert.equal(s.picker.pick(CENTER, CENTER).entry, "body");
  behind.z = 20;
  assert.equal(s.picker.pick(CENTER, CENTER).layer.name, "points");
});

test("a disabled layer and one left out by `only` are never picked", () => {
  const s = scene();
  let on = false;
  s.picker.addLayer({ name: "points", enabled: () => on, points: () => [{ x: 0, y: 0, z: 0 }], reach: () => 8 });
  const blocks = s.picker.addLayer({ name: "blocks", meshes: () => [plate(0)], entryOf: () => "block" });
  assert.equal(s.picker.pick(CENTER, CENTER).layer.name, "blocks");
  on = true;
  assert.equal(s.picker.pick(CENTER, CENTER).layer.name, "points");
  assert.equal(s.picker.pick(CENTER, CENTER, (layer) => layer === blocks).layer.name, "blocks");
  assert.equal(s.picker.pick(SIZE - 2, SIZE - 2), null, "nothing there");
});

test("a custom layer answers with its own entry and distance", () => {
  const s = scene();
  s.picker.addLayer({
    name: "clouds",
    pick: (ctx) => (ctx.ray.direction.z < 0 ? { entry: "cloud", distance: 50 } : null),
  });
  assert.equal(s.picker.pick(CENTER, CENTER).entry, "cloud");
});

test("the tooltip shows beside the pointer, stays inside its map, and hides", () => {
  installDom("http://localhost/");
  const tip = h("div", { hidden: true });
  const viewport = h("div", {}, [tip]);
  viewport.rect = { left: 0, top: 0, width: 300, height: 200 };
  tip.offsetWidth = 100;
  tip.offsetHeight = 20;
  document.body.append(viewport);
  const tooltip = MP.createTooltip(tip);
  tooltip.show("Sol, G2V", 10, 10);
  assert.equal(tip.hidden, false);
  assert.equal(tip.textContent, "Sol, G2V");
  assert.equal(tip.style.left, "22px");
  tooltip.show("Sol, G2V", 290, 190);
  assert.equal(tip.style.left, "196px", "pulled back inside on the right");
  assert.equal(tip.style.top, "176px", "and at the bottom");
  tooltip.show("", 0, 0);
  assert.equal(tip.hidden, true);
  MP.createTooltip(null).show("no element", 1, 1);
});

test("the info panel shows the rows it has a value for, then the actions", () => {
  installDom("http://localhost/");
  const panelEl = h("aside", {});
  document.body.append(panelEl, h("span", { "data-bookmark-db": "test" }));
  const panel = MP.infoPanelOf(panelEl);
  assert.equal(MP.infoPanelOf(panelEl), panel, "one panel per element");
  panel.show({
    title: "Sol",
    fields: [["Star type", "G2V"], ["Octant", null], ["Systems", 0]],
    nav: [{ label: "Start Here", href: "/nav?from=system:1" }, { label: "End Here", href: "/nav?to=system:1" }],
    bookmark: MP.endpointBookmark("system:1", "Sol", "/system/1"),
    links: [{ href: "/system/1", label: "View system →" }],
    hint: "Click another star.",
  });
  assert.equal(panelEl.querySelector("h3").textContent, "Sol");
  assert.deepEqual(panelEl.querySelectorAll("dt").map((dt) => dt.textContent), ["Star type", "Systems"]);
  assert.ok(panelEl.querySelector('a[href="/nav?from=system:1"]'));
  assert.ok(panelEl.querySelector('a[href="/nav?to=system:1"]'));
  assert.ok(panelEl.querySelector('a[href="/system/1"]'));
  assert.equal(panelEl.querySelector(".map-info-bookmark").textContent, "☆");
  assert.equal(panelEl.querySelector("p.hint").textContent, "Click another star.");
});

test("while a NAV end is picked the panel offers only the pick button", () => {
  installDom("http://localhost/");
  const panelEl = h("aside", {});
  document.body.append(panelEl);
  MP.infoPanelOf(panelEl).show({
    title: "Sol",
    nav: [{ label: "Start Here", href: "/nav?from=system:1&to=system:2", primary: true }],
  });
  const links = panelEl.querySelectorAll("a");
  assert.equal(links.length, 1);
  assert.equal(links[0].textContent, "Start Here");
  assert.equal(links[0].getAttribute("href"), "/nav?from=system:1&to=system:2");
});

test("NAV buttons without a link are buttons that act in place", () => {
  installDom("http://localhost/");
  const panelEl = h("aside", {});
  document.body.append(panelEl);
  const clicks = [];
  MP.infoPanelOf(panelEl).show({
    title: "Sol",
    nav: [
      { label: "Start Here", onClick: () => clicks.push("start") },
      { label: "End Here", onClick: () => clicks.push("end") },
    ],
  });
  const buttons = panelEl.querySelectorAll("button");
  assert.deepEqual([...buttons].map((b) => b.textContent), ["Start Here", "End Here"]);
  buttons[1].click();
  assert.deepEqual(clicks, ["end"]);
});

test("a NAV endpoint's bookmark entry takes its kind from the endpoint", () => {
  assert.deepEqual(MP.endpointBookmark("nebula:3", "Veil", "/phenomenon/nebula/3"),
    { kind: "nebula", value: "nebula:3", name: "Veil", url: "/phenomenon/nebula/3" });
  assert.equal(MP.endpointBookmark("", "x"), null);
  assert.equal(MP.endpointBookmark(undefined, "x"), null);
});
