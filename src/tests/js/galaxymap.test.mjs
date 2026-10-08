// tests/js/galaxymap.test.mjs -- the Galaxy Map's own logic (TEST.57):
// static/galaxystageview.js (the drill-down: stages, breadcrumb, the
// map's Back and Forward through browser history, the address bar, and
// the free view's zoom, pan and tilt clamps) and static/galaxymap3d.js
// (sector designations, addresses and the control buttons' handlers).
//
// The stage view runs for real against three.js, galaxyblocks.js and the
// galaxy's density shape from Python (test_js_unit.py's fixtures), with
// a fake host in place of the WebGL page: picking by pointer is real
// raycasting against the drawn blocks; only nothing is rendered.

import { test } from "node:test";
import assert from "node:assert/strict";

import { FakeEvent, h, installDom, settle, STATIC_URL } from "./fakedom.mjs";

const F = JSON.parse(process.env.PLANETGEN_JS_FIXTURES || "null");
if (!F) throw new Error("run through test_js_unit.py (PLANETGEN_JS_FIXTURES is not set)");

installDom("http://localhost/galaxy");
const THREE = await import(new URL("vendor/three.module.min.js", STATIC_URL).href);
const S = await import(new URL("galaxystages.js", STATIC_URL).href);
const GB = await import(new URL("galaxyblocks.js", STATIC_URL).href);
// Without a query, so they share one three.js with this file (the
// stage view gets it from the host anyway).
const SV = await import(new URL("galaxystageview.js", STATIC_URL).href);
const G3 = await import(new URL("galaxymap3d.js", STATIC_URL).href);

const WIDTH = 800;
const HEIGHT = 600;
const FOV = 50;

// --- A fake page around the stage view --------------------------------------

function makeBlockMesh(part, translucent) {
  const geometry = new THREE.BufferGeometry();
  geometry.setAttribute("position", new THREE.BufferAttribute(part.positions, 3));
  geometry.setIndex(new THREE.BufferAttribute(part.indices, 1));
  const material = new THREE.ShaderMaterial({ uniforms: { fade: { value: 1 }, gridEdges: { value: 1 } }, transparent: translucent });
  const mesh = new THREE.Mesh(geometry, material);
  mesh.frustumCulled = false;
  mesh.userData.part = part;
  return mesh;
}

function controlButtons() {
  return F.galaxyControls.admin.map((action) => h("button", { type: "button", "data-action": action }, action));
}

// The stage view with a fake host, on a fresh window at `href`. `server`
// answers the stage API: (query) -> {children, sectors}.
function open(href, options) {
  options = options || {};
  let win = options.reload;
  if (win) {
    // The same tab reloaded: its history stays, the page is new.
    document.body.replaceChildren();
  } else {
    win = installDom(href || "http://localhost/galaxy");
  }
  const scene = new THREE.Scene();
  const camera = new THREE.PerspectiveCamera(FOV, WIDTH / HEIGHT, 0.01, 1e7);
  camera.up.set(0, 0, 1);
  const canvasEl = h("canvas", { tabindex: "0" });
  canvasEl.rect = { left: 0, top: 0, width: WIDTH, height: HEIGHT };
  canvasEl.clientWidth = WIDTH;
  canvasEl.clientHeight = HEIGHT;
  const tooltip = h("div", { hidden: true });
  const viewport = h("div", { class: "starmap-viewport" }, [canvasEl, tooltip]);
  viewport.rect = canvasEl.rect;
  const els = {
    crumbs: h("nav"), slabs: h("div"), tooltip: tooltip, notice: h("p", { hidden: true }),
    address: h("form", {}, [h("input", { type: "search" })]), matches: h("div", { hidden: true }),
    controls: h("div", { id: "galaxymap3d-controls" }, controlButtons()),
  };
  document.body.append(h("div", { "data-bookmark-db": "planetgen" }), viewport, els.crumbs, els.slabs, els.notice,
    els.address, els.matches, els.controls);
  const calls = { cameras: [], fetches: [], blockInfo: [], placed: [], cells: [], hints: [], sizes: [], clips: [], locates: [], scenes: [], infos: [], hidden: [] };
  const server = options.server || (() => ({ children: [], sectors: [] }));
  const host = {
    THREE, scene, camera, canvasEl, edgePc: F.edgePc, shape: F.shape, galaxyRadius: F.galaxyRadius,
    reducedMotion: true, accentColor: "#4f5fe8", canGenerate: !!options.canGenerate, courseSectors: options.courseSectors || [],
    blockScene: GB.createBlockScene({
      edgePc: F.edgePc, galaxyRadius: F.galaxyRadius, shape: F.shape,
      palette: { dim: [0.2, 0.2, 0.3], accent: [0.3, 0.4, 0.9], hot: [1, 1, 1], placedLow: [0.3, 0.25, 0.2], placedHigh: [1, 0.96, 0.87] },
    }),
    makeBlockMesh,
    setCamera(v) {
      camera.quaternion.fromArray(v.quaternion);
      const back = new THREE.Vector3(0, 0, v.dist).applyQuaternion(camera.quaternion);
      camera.position.set(v.target[0] + back.x, v.target[1] + back.y, v.target[2] + back.z);
      camera.updateMatrixWorld();
      // The tilt from straight down and the screen's up direction on the
      // galaxy plane, for the checks below.
      const up = new THREE.Vector3(0, 1, 0).applyQuaternion(camera.quaternion);
      calls.cameras.push({ ...v, tilt: Math.acos(Math.max(-1, Math.min(1, back.z / v.dist))), heading: Math.atan2(up.y, up.x) });
    },
    fetchStage(query) {
      calls.fetches.push(query);
      try {
        return Promise.resolve(server(query));
      } catch (err) {
        return Promise.reject(err);
      }
    },
    showBlockInfo: (info) => calls.blockInfo.push(info),
    showPlacedInfo: (entry) => calls.placed.push(entry),
    showCellInfo: (cell) => calls.cells.push(cell),
    showHint: (text) => calls.hints.push(text),
    setWedgeClip: (clip) => calls.clips.push(clip),
    sectorUrl: (id) => "/sector/" + id,
    // The sector's scene JSON: the Sector Map fixture, placed.
    fetchSectorScene(id) {
      calls.scenes.push(id);
      const data = structuredClone(F.sectorMap.data);
      data.centerPc = data.centerPc || [1000, 2000, 30];
      data.halfEdgePc = data.halfEdgePc || 5;
      return Promise.resolve(options.sectorScene ? options.sectorScene(data) : data);
    },
    lightBackground: false,
    pixelRatio: () => 1,
    showInfo: (spec) => calls.infos.push(spec),
    hideCell: (bounds) => calls.hidden.push(bounds),
    locate(name) {
      calls.locates.push(name);
      return options.locate ? options.locate(name) : Promise.resolve([]);
    },
    els,
  };
  const view = SV.createStageView(host);
  return { win, view, host, els, calls, camera };
}

// Lets the stage's data load and any flight finish.
async function arrive(m) {
  for (let n = 0; n < 4; n++) {
    await settle();
    m.view.step(performance.now() + 1e9);
  }
}

async function start(href, options) {
  const m = open(href, options);
  m.view.openFromLocation();
  m.view.setActive(true);
  await arrive(m);
  return m;
}

function key(m, name) {
  const event = new FakeEvent("keydown", { key: name });
  m.view.onKey(event);
  return event;
}

function button(m, action) {
  return m.els.controls.querySelector(`[data-action="${action}"]`);
}

function wire(m, extra) {
  const ctx = Object.assign({
    stageView: m.view,
    setTerritories() {},
    territoriesWanted: () => false,
  }, extra || {});
  G3.wireMapControls(m.els.controls, G3.mapControlHandlers(ctx));
  return ctx;
}

function crumbLabels(m) {
  return m.els.crumbs.querySelectorAll("li").map((li) => li.textContent);
}

function expectedCrumbs(m) {
  const outline = S.galaxyOutline(F.edgePc, F.shape, F.galaxyRadius);
  return S.crumbs(m.view.stage(), outline, F.edgePc).map((c) => c.label);
}

// Picks the first choice by keyboard and takes it; true when the stage
// moved on.
async function drillOnce(m) {
  const before = JSON.stringify(m.view.stage());
  key(m, "ArrowRight");
  key(m, "Enter");
  await arrive(m);
  return JSON.stringify(m.view.stage()) !== before;
}

function lastCamera(m) {
  return m.calls.cameras[m.calls.cameras.length - 1];
}

// --- Pure helpers --------------------------------------------------------------

test("clampZoom keeps the distance between 1/8 and 2.5 times the fit", () => {
  assert.equal(SV.clampZoom(1, 100), 100 * SV.MIN_ZOOM);
  assert.equal(SV.clampZoom(100, 100), 100);
  assert.equal(SV.clampZoom(1e9, 100), 100 * SV.MAX_ZOOM);
  assert.equal(SV.MIN_ZOOM, 1 / 8);
  assert.equal(SV.MAX_ZOOM, 2.5);
});

test("clampPan pulls the view's middle back onto the reach sphere", () => {
  const fit = [10, 20, 30];
  assert.deepEqual(SV.clampPan([11, 20, 30], fit, 5), [11, 20, 30], "inside the reach: untouched");
  const pulled = SV.clampPan([10 + 30, 20 + 40, 30], fit, 5);
  assert.ok(Math.abs(pulled[0] - 13) < 1e-9 && Math.abs(pulled[1] - 24) < 1e-9 && pulled[2] === 30, String(pulled));
  assert.deepEqual(SV.clampPan([10, 20, 30], fit, 0), [10, 20, 30]);
});

test("sectorDesignation packs like Python's provisional_sector_designation", () => {
  assert.ok(F.designations.length > 50);
  for (const d of F.designations) {
    assert.equal(G3.sectorDesignation(d.ring, d.layer, d.slot), d.designation, JSON.stringify(d));
    assert.equal(S.sectorDesignation(d.ring, d.layer, d.slot), d.designation, JSON.stringify(d));
    assert.deepEqual(S.parseSectorDesignation(d.designation), { ring: d.ring, layer: d.layer, slot: d.slot });
  }
  // Past 32 bits of packing: BigInt, not a Number shift (which wraps).
  const wide = F.designations.find((d) => d.ring >= 4);
  assert.ok(wide.designation.length > 8);
});

test("formatAddress names ring, layer and slot", () => {
  assert.equal(G3.formatAddress(312, -3, 1042), "ring 312 layer -3 slot 1042");
});

test("every control button the panel draws has a handler", () => {
  const handlers = G3.mapControlHandlers({});
  for (const action of new Set(F.galaxyControls.public.concat(F.galaxyControls.admin))) {
    assert.equal(typeof handlers[action], "function", `the ${action} button would do nothing`);
  }
});

test("the panel has no Wedges or Whole galaxy button (MAP.55, MAP.85)", () => {
  const actions = F.galaxyControls.public.concat(F.galaxyControls.admin);
  assert.equal(actions.includes("wedges"), false);
  for (const action of ["back", "forward", "up", "reset", "reset-view", "charted-only"]) {
    assert.ok(actions.includes(action), action);
  }
});

test("Escape or a press outside closes the controls' Menu", () => {
  installDom("http://localhost/galaxy");
  const summary = h("summary", {}, "Menu");
  const inside = h("button", { "data-action": "reset-view" });
  const menu = h("details", {}, [summary, h("div", {}, [inside])]);
  document.body.append(menu);
  G3.wireMapMenu(menu);
  menu.open = true;
  const other = new FakeEvent("keydown", { key: "a", bubbles: true });
  inside.dispatchEvent(other);
  assert.equal(menu.open, true);
  const escape = new FakeEvent("keydown", { key: "Escape", bubbles: true });
  inside.dispatchEvent(escape);
  assert.equal(menu.open, false);
  assert.ok(escape.defaultPrevented);
  assert.equal(document.activeElement, summary);
  // A press inside leaves it open; one anywhere else closes it.
  const outside = h("p");
  document.body.append(outside);
  menu.open = true;
  inside.dispatchEvent(new FakeEvent("pointerdown", { bubbles: true }));
  assert.equal(menu.open, true);
  outside.dispatchEvent(new FakeEvent("pointerdown", { bubbles: true }));
  assert.equal(menu.open, false);
});

test("the toggle buttons flip aria-pressed and what they control", () => {
  installDom("http://localhost/galaxy");
  const calls = [];
  let wanted = false;
  const controls = h("div", {}, [
    h("button", { "data-action": "charted-only", "aria-pressed": "false" }),
    h("button", { "data-action": "territories", "aria-pressed": "false" }),
    h("button", { "data-action": "no-such-action" }),
  ]);
  G3.wireMapControls(controls, G3.mapControlHandlers({
    setChartedOnly: (on) => calls.push(["charted-only", on]),
    setTerritories: (on, btn) => { calls.push(["territories", on, btn.dataset.action]); wanted = on; },
    territoriesWanted: () => wanted,
  }));
  const [only, territories, unknown] = controls.children;
  only.click();
  assert.equal(only.getAttribute("aria-pressed"), "true");
  only.click();
  territories.click();
  assert.equal(territories.getAttribute("aria-pressed"), "true");
  unknown.click();
  assert.deepEqual(calls, [["charted-only", true], ["charted-only", false], ["territories", true, "territories"]]);
});

// --- The drill-down ------------------------------------------------------------

test("opens on the whole galaxy with nothing to go back to", async () => {
  const m = await start();
  assert.deepEqual(m.view.stage(), { at: null, picks: [] });
  assert.deepEqual(crumbLabels(m), ["Galaxy"]);
  assert.equal(button(m, "back").disabled, true);
  assert.equal(button(m, "forward").disabled, true);
  assert.equal(button(m, "up").disabled, true);
  assert.equal(button(m, "reset-view").disabled, false, "the whole galaxy turns (MAP.85)");
  assert.equal(m.calls.fetches[0], "", "the galaxy's own counts");
  const cam = lastCamera(m);
  assert.ok(cam.tilt < 1e-6, "the galaxy starts top-down (MAP.97)");
});

// Every drawn block mesh's gridEdges uniform.
function gridEdges(m) {
  const values = [];
  m.host.scene.traverse((object) => {
    if (object.isMesh && object.material.uniforms && object.material.uniforms.gridEdges) {
      values.push(object.material.uniforms.gridEdges.value);
    }
  });
  return Array.from(new Set(values));
}

test("the whole galaxy has no grid lines; an arc shows only its slabs' lines (MAP.85, MAP.77)", async () => {
  const m = await start();
  assert.deepEqual(gridEdges(m), [0]);
  await drillOnce(m);
  assert.equal(m.view.stage().picks[0].kind, "arc");
  // Picking a slab: no lines between the blocks inside the slabs, only
  // the boundaries between slabs.
  assert.deepEqual(gridEdges(m), [0]);
  const lines = m.host.canvasEl.galaxyLines();
  assert.equal(lines.kind, "layer");
  assert.ok(lines.slabLines > 0, "the slabs' boundaries are drawn");
});

test("the whole galaxy turns and zooms in to half its fit, never out past it", async () => {
  const m = await start();
  const fit = lastCamera(m);
  const wheel = (deltaY, deltaMode) => m.view.onWheel(new FakeEvent("wheel", { deltaY, deltaMode: deltaMode || 0 }));
  assert.equal(wheel(100), true, "the wheel zooms");
  assert.ok(Math.abs(lastCamera(m).dist - fit.dist) < 1e-9 * fit.dist, "no further out than the fit");
  for (let n = 0; n < 100; n++) wheel(-3, 1);
  assert.ok(Math.abs(lastCamera(m).dist - fit.dist * SV.GALAXY_MIN_ZOOM) < 1e-6 * fit.dist, "twice as close at most");
  m.view.onPointerDown(new FakeEvent("pointerdown", { clientX: 400, clientY: 300, pointerId: 1, button: 0, pointerType: "mouse" }));
  m.view.onPointerMove(new FakeEvent("pointermove", { clientX: 480, clientY: 260, pointerId: 1, pointerType: "mouse" }));
  m.view.onPointerUp(new FakeEvent("pointerup", { clientX: 480, clientY: 260, pointerId: 1, button: 0, pointerType: "mouse" }));
  const turned = lastCamera(m);
  assert.ok(Math.abs(turned.heading - fit.heading) > 0.05 && turned.tilt > 0.05, "dragging turns and tilts it");
  assert.deepEqual(m.view.stage(), { at: null, picks: [] }, "a drag picks nothing");
});

test("each step down by keyboard pushes a history entry and a breadcrumb", async () => {
  const m = await start();
  const seen = [];
  for (let n = 0; n < 12; n++) {
    if (!(await drillOnce(m))) break;
    seen.push(m.view.stage());
    assert.equal(m.win.location.search, S.stageQuery(m.view.stage()), "the URL names the stage");
    // A sector picked in a cube of them is selected too, after the stage.
    const labels = crumbLabels(m);
    if (labels[labels.length - 1].startsWith("Sector ")) labels.pop();
    assert.deepEqual(labels, expectedCrumbs(m));
    assert.equal(m.win.history.length, seen.length + 1);
    assert.equal(button(m, "back").disabled, false);
    assert.equal(button(m, "up").disabled, false);
  }
  assert.ok(seen.length >= 5, `only ${seen.length} steps down`);
  assert.ok(seen.some((s) => s.at && s.at.m === 3), "reached a level-3 block of sectors");
});

test("Back, Forward, Up and Reset buttons move through the stages", async () => {
  const m = await start();
  wire(m);
  await drillOnce(m);
  const first = m.view.stage();
  await drillOnce(m);
  const second = m.view.stage();
  assert.notDeepEqual(first, second);

  button(m, "back").click();
  await arrive(m);
  assert.deepEqual(m.view.stage(), first);
  assert.equal(m.win.location.search, S.stageQuery(first));
  assert.equal(button(m, "forward").disabled, false);

  button(m, "forward").click();
  await arrive(m);
  assert.deepEqual(m.view.stage(), second);
  assert.equal(button(m, "forward").disabled, true);

  button(m, "up").click();
  await arrive(m);
  assert.deepEqual(m.view.stage(), first, "up goes to the parent stage");
  assert.equal(button(m, "forward").disabled, true, "a new step clears Forward");

  button(m, "reset").click();
  await arrive(m);
  assert.deepEqual(m.view.stage(), { at: null, picks: [] });
  assert.deepEqual(crumbLabels(m), ["Galaxy"]);
});

test("Back past the map's first stage is never taken", async () => {
  const m = await start();
  m.view.travel(-1);
  await arrive(m);
  assert.equal(m.win.history.index, 0);
  m.view.travel(1);
  await arrive(m);
  assert.equal(m.win.history.index, 0);
});

test("a reload keeps the map's place in its Back and Forward history", async () => {
  const m = await start();
  await drillOnce(m);
  const first = m.view.stage();
  await drillOnce(m);
  m.view.travel(-1);
  await arrive(m);
  assert.deepEqual(m.view.stage(), first);
  m.view.setActive(false);

  const again = await start(undefined, { reload: m.win });
  assert.deepEqual(again.view.stage(), first);
  assert.equal(button(again, "back").disabled, false);
  assert.equal(button(again, "forward").disabled, false, "the stage after this one is still ahead");
  again.view.travel(1);
  await arrive(again);
  assert.equal(again.win.history.index, 2);
});

test("a breadcrumb button goes back to that stage", async () => {
  const m = await start();
  await drillOnce(m);
  await drillOnce(m);
  await drillOnce(m);
  const crumbButtons = m.els.crumbs.querySelectorAll("button.galaxy-crumb");
  assert.equal(crumbButtons.length, 3);
  crumbButtons[0].click();
  await arrive(m);
  assert.deepEqual(m.view.stage(), { at: null, picks: [] });
  assert.equal(m.els.crumbs.querySelectorAll('[aria-current="location"]').length, 1);
});

test("Escape goes up and Home goes back to the galaxy", async () => {
  const m = await start();
  await drillOnce(m);
  const first = m.view.stage();
  await drillOnce(m);
  assert.ok(key(m, "Escape").defaultPrevented);
  await arrive(m);
  assert.deepEqual(m.view.stage(), first);
  await drillOnce(m);
  key(m, "Home");
  await arrive(m);
  assert.deepEqual(m.view.stage(), { at: null, picks: [] });
});

test("a stage URL opens that stage; a bad one opens the galaxy and says why", async () => {
  const deep = await start();
  for (let n = 0; n < 4; n++) await drillOnce(deep);
  const query = S.stageQuery(deep.view.stage());

  const m = await start("http://localhost/galaxy" + query);
  assert.deepEqual(m.view.stage(), deep.view.stage());
  assert.deepEqual(crumbLabels(m), crumbLabels(deep));

  const bad = await start("http://localhost/galaxy?at=nonsense");
  assert.deepEqual(bad.view.stage(), { at: null, picks: [] });
  assert.equal(bad.els.notice.hidden, false);
  assert.match(bad.els.notice.textContent, /There is no block nonsense/);
});

test("?sector= opens the stage around that sector with it selected", async () => {
  const d = { ring: 400, layer: 0, slot: 1000 };
  const designation = S.sectorDesignation(d.ring, d.layer, d.slot);
  const m = await start("http://localhost/galaxy?sector=" + designation);
  const outline = S.galaxyOutline(F.edgePc, F.shape, F.galaxyRadius);
  assert.deepEqual(m.view.stage(), S.sectorStage(d.ring, d.layer, d.slot, outline, F.edgePc));
  const crumbs = crumbLabels(m);
  assert.equal(crumbs[crumbs.length - 1], "Sector " + S.blockLabel({ m: 1, ring: d.ring, wedge: d.slot, slab: d.layer }));
  assert.equal(m.calls.cells.length, 1, "the sector's panel is shown");
  assert.deepEqual(m.calls.cells[0].address, d);
});

test("the address bar flies to a sector, or says what is wrong", async () => {
  const m = await start();
  const input = m.els.address.querySelector("input");
  input.value = "400/0/1000";
  const submit = new FakeEvent("submit");
  m.els.address.dispatchEvent(submit);
  assert.ok(submit.defaultPrevented, "the form isn't sent");
  await arrive(m);
  assert.equal(m.win.location.search, "?sector=" + S.sectorDesignation(400, 0, 1000));
  assert.equal(m.els.notice.hidden, true);

  input.value = "0/0/99";
  m.els.address.dispatchEvent(new FakeEvent("submit"));
  assert.equal(m.els.notice.textContent, "Ring 0 has no slot 99.");

  input.value = "";
  m.els.address.dispatchEvent(new FakeEvent("submit"));
  assert.match(m.els.notice.textContent, /^Type a designation/);

  const far = 1e6;
  input.value = `${far}, 0, 0`;
  m.els.address.dispatchEvent(new FakeEvent("submit"));
  await arrive(m);
  assert.match(m.els.notice.textContent, /is outside the galaxy\.$/);
});

test("the address bar looks names up: one match flies there, several are listed", async () => {
  const sector = { kind: "sector", name: "Bead", sector_id: 7, ring: 400, layer: 0, slot: 1000 };
  const answers = {
    Bead: [sector],
    Nowhere: [],
    Many: [Object.assign({}, sector, { name: "Many A" }), { kind: "system", name: "Many B", sector_id: 8, sector_name: "Far", ring: 401, layer: 1, slot: 10 }],
  };
  const m = await start(undefined, {
    locate: (name) => (name === "Broken" ? Promise.reject(new Error("503")) : Promise.resolve(answers[name])),
  });
  const input = m.els.address.querySelector("input");
  const ask = async (text) => {
    input.value = text;
    m.els.address.dispatchEvent(new FakeEvent("submit"));
    await arrive(m);
  };

  await ask("Bead");
  assert.deepEqual(m.calls.locates, ["Bead"]);
  assert.equal(m.win.location.search, "?sector=" + S.sectorDesignation(400, 0, 1000));

  await ask("Nowhere");
  assert.equal(m.els.notice.textContent, "Nothing is named like Nowhere.");

  await ask("Broken");
  assert.equal(m.els.notice.textContent, "The lookup failed. Please try again shortly.");

  await ask("Many");
  assert.equal(m.els.notice.textContent, "2 names match; pick one.");
  assert.equal(m.els.matches.hidden, false);
  const choices = m.els.matches.querySelectorAll("button");
  assert.deepEqual(choices.map((b) => b.textContent), ["Many A (sector)", "Many B (system in Far)"]);
  choices[1].click();
  await arrive(m);
  assert.equal(m.els.matches.hidden, true);
  assert.equal(m.win.location.search, "?sector=" + S.sectorDesignation(401, 1, 10));
});

// --- The free view -------------------------------------------------------------

async function freeStage() {
  const m = await start();
  // An arc of the galaxy, seen at the isometric slant.
  await drillOnce(m);
  assert.equal(m.view.stage().picks[0].kind, "arc");
  assert.ok(Math.abs(lastCamera(m).tilt - SV.ISO_TILT) < 1e-9, "an arc is seen at a slant");
  assert.equal(button(m, "reset-view").disabled, false);
  return m;
}

test("the wheel zooms a free stage, within 1/8 to 2.5 times its fit", async () => {
  const m = await freeStage();
  const fit = lastCamera(m).dist;
  const wheel = (deltaY, deltaMode) => m.view.onWheel(new FakeEvent("wheel", { deltaY, deltaMode: deltaMode || 0 }));
  assert.equal(wheel(100), true);
  assert.ok(lastCamera(m).dist > fit, "scrolling down zooms out");
  for (let n = 0; n < 50; n++) wheel(200);
  assert.ok(Math.abs(lastCamera(m).dist - fit * SV.MAX_ZOOM) < 1e-6 * fit);
  for (let n = 0; n < 100; n++) wheel(-3, 1);
  assert.ok(Math.abs(lastCamera(m).dist - fit * SV.MIN_ZOOM) < 1e-6 * fit);
});

test("dragging turns a free stage any way, past edge-on and under the plane (MAP.96)", async () => {
  const m = await freeStage();
  const before = lastCamera(m);
  m.view.onPointerDown(new FakeEvent("pointerdown", { clientX: 400, clientY: 300, pointerId: 1, button: 0, pointerType: "mouse" }));
  m.view.onPointerMove(new FakeEvent("pointermove", { clientX: 450, clientY: 300, pointerId: 1, pointerType: "mouse" }));
  const turned = lastCamera(m);
  assert.ok(Math.abs(turned.quaternion[2] - before.quaternion[2]) + Math.abs(turned.quaternion[3] - before.quaternion[3]) > 0.01, "it turned");
  let most = 0;
  for (let n = 1; n <= 40; n++) {
    m.view.onPointerMove(new FakeEvent("pointermove", { clientX: 450, clientY: 300 - 20 * n, pointerId: 1, pointerType: "mouse" }));
    most = Math.max(most, lastCamera(m).tilt);
  }
  assert.ok(most > Math.PI / 2 + 0.1, `turns past edge-on to under the plane (${most})`);
  m.view.onPointerUp(new FakeEvent("pointerup", { clientX: 450, clientY: -500, pointerId: 1, button: 0, pointerType: "mouse" }));
  const stage = m.view.stage();
  await arrive(m);
  assert.deepEqual(m.view.stage(), stage, "a drag picks nothing");
});

test("right-dragging moves a free stage, no further than 1.5 fits", async () => {
  const m = await freeStage();
  const fit = lastCamera(m);
  m.view.onPointerDown(new FakeEvent("pointerdown", { clientX: 400, clientY: 300, pointerId: 1, button: 2, pointerType: "mouse" }));
  for (let n = 1; n <= 40; n++) {
    m.view.onPointerMove(new FakeEvent("pointermove", { clientX: 400 + 200 * n, clientY: 300, pointerId: 1, pointerType: "mouse" }));
  }
  const moved = lastCamera(m);
  const off = Math.hypot(moved.target[0] - fit.target[0], moved.target[1] - fit.target[1], moved.target[2] - fit.target[2]);
  const vertical = THREE.MathUtils.degToRad(FOV) / 2;
  const half = Math.min(vertical, Math.atan(Math.tan(vertical) * (WIDTH / HEIGHT)));
  const reach = SV.PAN_REACH * fit.dist * Math.tan(half);
  assert.ok(off > 0.9 * reach && off <= reach * (1 + 1e-9), `moved ${off} of ${reach}`);
  assert.deepEqual(moved.quaternion, fit.quaternion, "moving doesn't turn");
  m.view.onPointerUp(new FakeEvent("pointerup", { clientX: 8400, clientY: 300, pointerId: 1, button: 2, pointerType: "mouse" }));
});

test("Reset view flies back to the stage's own view", async () => {
  const m = await freeStage();
  wire(m);
  const fit = lastCamera(m);
  m.view.onWheel(new FakeEvent("wheel", { deltaY: 200 }));
  m.view.onPointerDown(new FakeEvent("pointerdown", { clientX: 400, clientY: 300, pointerId: 1, button: 0, pointerType: "mouse" }));
  m.view.onPointerMove(new FakeEvent("pointermove", { clientX: 500, clientY: 250, pointerId: 1, pointerType: "mouse" }));
  m.view.onPointerUp(new FakeEvent("pointerup", { clientX: 500, clientY: 250, pointerId: 1, button: 0, pointerType: "mouse" }));
  assert.notDeepEqual(lastCamera(m), fit);
  button(m, "reset-view").click();
  await arrive(m);
  const back = lastCamera(m);
  assert.ok(Math.abs(back.dist - fit.dist) < 1e-9 * fit.dist, "dist");
  for (let i = 0; i < 4; i++) assert.ok(Math.abs(Math.abs(back.quaternion[i]) - Math.abs(fit.quaternion[i])) < 1e-9, "quaternion " + i);
  assert.deepEqual(back.target, fit.target);
});

test("a pinch zooms a free stage within the same limits", async () => {
  const m = await freeStage();
  const fit = lastCamera(m).dist;
  const touch = (type, id, x, y) => m.view["onPointer" + type](new FakeEvent("pointer" + type.toLowerCase(), { clientX: x, clientY: y, pointerId: id, button: 0, pointerType: "touch" }));
  touch("Down", 1, 300, 300);
  touch("Down", 2, 500, 300);
  touch("Move", 2, 5000, 300);
  assert.ok(Math.abs(lastCamera(m).dist - fit * SV.MIN_ZOOM) < 1e-6 * fit, "spread apart: closest zoom");
  touch("Move", 2, 301, 300);
  assert.ok(Math.abs(lastCamera(m).dist - fit * SV.MAX_ZOOM) < 1e-6 * fit, "pinched together: furthest zoom");
  touch("Up", 2, 301, 300);
  touch("Up", 1, 300, 300);
});

// --- Slab buttons (MAP.110) -----------------------------------------------------

test("slabOrder sorts by slab number, either way, whatever order it is given", () => {
  const rows = [3, -1, 0, 7, 2].map((lo) => ({ pick: { lo: lo, hi: lo } }));
  const los = (list) => list.map((r) => r.pick.lo);
  assert.deepEqual(los(S.slabOrder(rows, true, (r) => r.pick)), [7, 3, 2, 0, -1]);
  assert.deepEqual(los(S.slabOrder(rows, false, (r) => r.pick)), [-1, 0, 2, 3, 7]);
  assert.deepEqual(los(rows), [3, -1, 0, 7, 2], "the rows given are left alone");
});

test("nearestOnPieces keeps a line's end at or below the one above it", () => {
  // An outline slanting up to the right: its point nearest a button at
  // the right is high up, above where the line above it ended.
  const pieces = [[[0, 100], [200, 20]]];
  const free = S.nearestOnPieces(pieces, [300, 30], null);
  assert.ok(free.y < 60, `nearest is high (${free.y})`);
  const kept = S.nearestOnPieces(pieces, [300, 30], null, 60);
  assert.ok(Math.abs(kept.y - 60) < 1e-9 && Math.abs(kept.x - 100) < 1e-9, JSON.stringify(kept));
  // Nothing of the outline reaches below: the nearest point of it all.
  assert.deepEqual(S.nearestOnPieces(pieces, [300, 30], null, 500), free);
  assert.equal(S.nearestOnPieces([], [0, 0], null), null);
});

// The slab buttons' numbers, top first, in each column.
function slabColumns(m) {
  const lists = [m.els.slabs, m.els.slabs.parentElement.querySelector("#galaxymap3d-slabs-side")].filter(Boolean)
    .map((box) => box.querySelector("ol.galaxy-slab-buttons")).filter(Boolean);
  return lists.map((list) => list.querySelectorAll(".galaxy-slab-button").map((b) => Number(b.dataset.slab)))
    .filter((column) => column.length);
}

function inNumberOrder(column) {
  const up = column.every((n, i) => i === 0 || n > column[i - 1]);
  const down = column.every((n, i) => i === 0 || n < column[i - 1]);
  return up ? "up" : down ? "down" : null;
}

test("the slab buttons stay in slab-number order however the view turns (MAP.110)", async () => {
  const m = await freeStage();
  // The map's row has a size, so the lines are drawn (and the buttons
  // ordered) as the view turns.
  m.els.slabs.parentElement.rect = { left: 0, top: 0, width: 1200, height: 700 };
  m.view.step(performance.now() + 1e9);
  let columns = slabColumns(m);
  assert.ok(columns.flat().length >= 2, `slab buttons drawn: ${JSON.stringify(columns)}`);
  columns.forEach((column) => assert.equal(inNumberOrder(column), "down", `top slab first from above: ${column}`));
  const seen = new Set();
  let most = 0;
  m.view.onPointerDown(new FakeEvent("pointerdown", { clientX: 400, clientY: 300, pointerId: 1, button: 0, pointerType: "mouse" }));
  for (let n = 1; n <= 40; n++) {
    m.view.onPointerMove(new FakeEvent("pointermove", { clientX: 400 + 7 * n, clientY: 300 - 20 * n, pointerId: 1, pointerType: "mouse" }));
    m.view.step(performance.now() + 1e9);
    most = Math.max(most, lastCamera(m).tilt);
    columns = slabColumns(m);
    columns.forEach((column) => {
      const way = inNumberOrder(column);
      assert.ok(way, `out of number order after turn ${n}: ${column}`);
      seen.add(way);
    });
  }
  m.view.onPointerUp(new FakeEvent("pointerup", { clientX: 680, clientY: -500, pointerId: 1, button: 0, pointerType: "mouse" }));
  assert.ok(most > Math.PI / 2 + 0.1, `the view turned under the plane (${most})`);
  assert.ok(seen.has("up"), "seen from below, the bottom slab (lowest number) shows on top and leads");
});

// --- Picking by pointer --------------------------------------------------------

// A screen point over a choice, found by moving the pointer over a grid
// and reading the tooltip: {x, y, text}.
function pointOver(m, wanted) {
  for (let y = 20; y < HEIGHT; y += 20) {
    for (let x = 20; x < WIDTH; x += 20) {
      m.view.onPointerMove(new FakeEvent("pointermove", { clientX: x, clientY: y, pointerId: 9, pointerType: "mouse" }));
      const tip = m.els.tooltip;
      if (!tip.hidden && tip.textContent && (!wanted || wanted.test(tip.textContent))) return { x, y, text: tip.textContent };
    }
  }
  return null;
}

function clickAt(m, x, y, pointerType) {
  const base = { clientX: x, clientY: y, pointerId: 3, button: 0, pointerType: pointerType || "mouse" };
  m.view.onPointerDown(new FakeEvent("pointerdown", base));
  m.view.onPointerUp(new FakeEvent("pointerup", base));
}

test("clicking an arc on the map drills into it", async () => {
  const m = await start();
  const spot = pointOver(m, /^Arc /);
  assert.ok(spot, "some arc shows a tooltip under the pointer");
  clickAt(m, spot.x, spot.y);
  await arrive(m);
  const stage = m.view.stage();
  assert.equal(stage.picks.length, 1);
  assert.equal(stage.picks[0].kind, "arc");
  assert.ok(spot.text.startsWith(crumbLabels(m)[1]), `${spot.text} vs ${crumbLabels(m)}`);
  assert.match(m.win.location.search, /^\?p=a[0-2]\.\d+$/, "the URL names the arc by band and bearing");
});

// The outline lines drawn over the map: [opacity, line count].
function outlines(m) {
  const found = [];
  m.host.scene.traverse((object) => {
    if (object.isLineSegments) found.push(object.material.opacity);
  });
  return found.sort((p, q) => q - p);
}

test("hovering an arc outlines it, its neighbors faintly, and nothing else (MAP.85)", async () => {
  const m = await start();
  assert.deepEqual(outlines(m), [], "nothing is outlined until the pointer is over an arc");
  const spot = pointOver(m, / \(middle\),/);
  assert.ok(spot, "a middle arc under the pointer");
  const lines = outlines(m);
  assert.equal(lines[0], 1, "the arc itself in full");
  // Either side in its band, and the arcs in and out of it.
  assert.equal(lines.length, 5, String(lines));
  assert.ok(lines.slice(1).every((opacity) => opacity > 0 && opacity < 1), "its neighbors faintly");
  m.view.onPointerLeave();
  assert.deepEqual(outlines(m), []);
});

test("on touch, the first tap highlights and the second takes it", async () => {
  const m = await start();
  const spot = pointOver(m, /^Arc /);
  m.view.onPointerLeave();
  clickAt(m, spot.x, spot.y, "touch");
  await arrive(m);
  assert.deepEqual(m.view.stage(), { at: null, picks: [] }, "one tap only shows what it would pick");
  assert.equal(m.els.tooltip.hidden, false);
  clickAt(m, spot.x, spot.y, "touch");
  await arrive(m);
  assert.equal(m.view.stage().picks.length, 1);
});

test("a generated sector opens in place; one not generated is selected", async () => {
  let known = null;
  const server = (query) => {
    if (!known) return { children: [], sectors: [] };
    return { children: [], sectors: [Object.assign({ id: 42, name: "Bead", system_count: 3 }, known)] };
  };
  const m = await start(undefined, { server });
  for (let n = 0; n < 12; n++) {
    if (S.parseStageQuery(m.win.location.search).stage.at?.m === 3 && !(await drillOnce(m))) break;
    if (!S.parseStageQuery(m.win.location.search).stage.at || S.parseStageQuery(m.win.location.search).stage.at.m !== 3) {
      await drillOnce(m);
    }
  }
  // In a level-3 block, down to one layer of sectors.
  for (let n = 0; n < 3; n++) await drillOnce(m);
  key(m, "ArrowRight");
  const tip = m.els.tooltip.textContent;
  const match = /^Sector ([\d,]+)·(-?\d+)·([\d,]+),.*not generated$/.exec(tip);
  assert.ok(match, "the tooltip names a sector: " + tip);
  const sector = { ring: Number(match[1].replace(/,/g, "")), layer: Number(match[2]), slot: Number(match[3].replace(/,/g, "")) };

  key(m, "Enter");
  await arrive(m);
  assert.equal(m.win.location.assigned.length, 0);
  const crumbs = crumbLabels(m);
  assert.equal(crumbs[crumbs.length - 1], "Sector " + S.blockLabel({ m: 1, ring: sector.ring, wedge: sector.slot, slab: sector.layer }));
  assert.deepEqual(m.calls.cells[m.calls.cells.length - 1].address, sector);

  known = sector;
  m.view.invalidate();
  await arrive(m);
  key(m, "Enter");
  await arrive(m);
  assert.deepEqual(m.win.location.assigned, [], "no page load: the sector opens in the map");
  assert.deepEqual(m.calls.scenes, [42]);
  assert.equal(S.parseStageQuery(m.win.location.search).open, true, "the URL says the sector is open: " + m.win.location.search);
  assert.deepEqual(S.parseStageQuery(m.win.location.search).sector, sector);
  assert.equal(m.calls.hidden.length, 1, "the galaxy's own stars leave the cell");
  assert.ok(m.calls.hidden[0].r1 > 0);
  assert.ok(m.view.stage(), "still on the stage that holds it");
  // Up closes it, keeping the sector selected and the cell's stars back.
  m.view.up();
  await arrive(m);
  assert.equal(m.calls.hidden[m.calls.hidden.length - 1], null);
  assert.equal(S.parseStageQuery(m.win.location.search).open, false);
  assert.deepEqual(S.parseStageQuery(m.win.location.search).sector, sector);
});
