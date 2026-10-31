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
const P = await import(new URL("galaxyprisms.js", STATIC_URL).href);
const GB = await import(new URL("galaxyblocks.js", STATIC_URL).href);
// Without a query, so they share one three.js with this file (the
// stage view gets it from the host anyway).
const SV = await import(new URL("galaxystageview.js", STATIC_URL).href);
const G3 = await import(new URL("galaxymap3d.js", STATIC_URL).href);
const NF = await import(new URL("nearfield.js", STATIC_URL).href);

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
  if (part.origin) mesh.position.set(part.origin[0], part.origin[1], part.origin[2]);   // MAP.240, as galaxymap3d.js
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
    crumbs: h("nav"), where: h("p", { hidden: true }), tooltip: tooltip, notice: h("p", { hidden: true }),
    address: h("form", {}, [h("input", { type: "search" })]), matches: h("div", { hidden: true }),
    controls: h("div", { id: "galaxymap3d-controls" }, controlButtons()),
  };
  document.body.append(h("div", { "data-bookmark-db": "planetgen" }), viewport, els.crumbs, els.where, els.notice,
    els.address, els.matches, els.controls);
  const calls = { cameras: [], fetches: [], blockInfo: [], placed: [], cells: [], hints: [], sizes: [], clips: [], locates: [], scenes: [], infos: [], hidden: [], pages: [], uncharted: [] };
  const server = options.server || (() => ({ children: [], sectors: [] }));
  const near = { D: 1, R: 0 };
  const host = {
    THREE, scene, camera, canvasEl, edgePc: F.edgePc, shape: F.shape, galaxyRadius: F.galaxyRadius,
    reducedMotion: true, accentColor: "#4f5fe8", canGenerate: !!options.canGenerate, courseSectors: options.courseSectors || [],
    blockScene: GB.createBlockScene({
      edgePc: F.edgePc, galaxyRadius: F.galaxyRadius, shape: F.shape,
      palette: { dim: [0.2, 0.2, 0.3], accent: [0.3, 0.4, 0.9], hot: [1, 1, 1], placedLow: [0.3, 0.25, 0.2], placedHigh: [1, 0.96, 0.87] },
    }),
    makeBlockMesh,
    // The near field of the camera last set (MAP.149), as galaxymap3d.js works it.
    nearFieldAt: (x, y, z) => NF.nearFieldAtWorld(camera.matrixWorldInverse.elements, x, y, z, near.D, near.R),
    setCamera(v) {
      near.D = v.dist;
      near.R = v.focusRadius || 0;
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
    // A window of blocks (MAP.151): the server answers it as it does a stage.
    fetchWindow(query) {
      calls.fetches.push("window" + query);
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
    sectorUrl: (id) => "/sector/" + id,
    pinned: options.pinned || null,
    openSectorPage: (id) => calls.pages.push(id),
    // The sector's scene JSON: the Sector Map fixture, placed.
    fetchSectorScene(id) {
      calls.scenes.push(id);
      const data = structuredClone(F.sectorMap.data);
      data.centerPc = data.centerPc || [1000, 2000, 30];
      data.halfEdgePc = data.halfEdgePc || 5;
      return Promise.resolve(options.sectorScene ? options.sectorScene(data) : data);
    },
    // MAP.162: a cell nothing was generated in, opened from what the scatters left there.
    fetchUnchartedScene: options.unchartedScene ? (at) => {
      calls.uncharted.push(at);
      const data = structuredClone(F.sectorMap.data);
      data.centerPc = [1000, 2000, 30];
      data.halfEdgePc = 5;
      Object.assign(data, { uncharted: true, designation: S.sectorDesignation(at.shell, at.band, at.wedge) });
      try {
        return Promise.resolve(options.unchartedScene(data, at));
      } catch (err) {
        return Promise.reject(err);
      }
    } : undefined,
    viewport: options.unchartedScene ? viewport : undefined,
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

// A query without its ?cam= (the camera rides the address once the view is still, MAP.228).
function withoutCam(search) {
  const params = new URLSearchParams(search);
  params.delete("cam");
  const text = params.toString().replace(/%2C/gi, ",");
  return text ? "?" + text : "";
}

function crumbLabels(m) {
  return m.els.crumbs.querySelectorAll("li").map((li) => li.textContent);
}

function expectedCrumbs(m) {
  const outline = S.galaxyOutline(F.edgePc, F.shape, F.galaxyRadius);
  return S.crumbs(m.view.stage(), outline, F.edgePc).map((c) => c.label);
}

// Opens the window one size down round the first block of the stage (what a wheel run into it carries the view to:
// the blocks themselves take no click since MAP.297); true when the stage moved on. A window of sectors has none.
async function drillOnce(m) {
  const before = JSON.stringify(m.view.stage());
  const outline = S.galaxyOutline(F.edgePc, F.shape, F.galaxyRadius);
  const stage = m.view.stage();
  const resolved = S.resolveStage(stage, outline, F.edgePc);
  const next = resolved.options && resolved.options.length ? S.windowInto(stage, resolved.options[0].blocks[0], null, F.edgePc) : null;
  if (!next) return false;
  m.view.flyToStage(next);
  await arrive(m);
  return JSON.stringify(m.view.stage()) !== before;
}

// A flight to a window moves the camera only (MAP.337): the map is reopened on the address it wrote, which opens the display on the window
// that holds the camera.
async function reopenOnAddress(m) {
  m.view.setActive(false);
  return start(undefined, { reload: m.win });
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

test("sectorDesignation packs like Python's provisional_sector_designation", () => {
  assert.ok(F.designations.length > 50);
  for (const d of F.designations) {
    assert.equal(G3.sectorDesignation(d.shell, d.band, d.wedge), d.designation, JSON.stringify(d));
    assert.equal(S.sectorDesignation(d.shell, d.band, d.wedge), d.designation, JSON.stringify(d));
    assert.deepEqual(S.parseSectorDesignation(d.designation), { shell: d.shell, band: d.band, wedge: d.wedge });
  }
  // Past 32 bits of packing: BigInt, not a Number shift (which wraps).
  const wide = F.designations.find((d) => d.shell >= 4000);
  assert.ok(wide.designation.length > 8);
});

test("formatAddress names shell, band and wedge", () => {
  assert.equal(G3.formatAddress(312, -3, 1042), "shell 312 band -3 wedge 1042");
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
    h("button", { "data-action": "cell-grid", "aria-pressed": "false" }),
    h("button", { "data-action": "territories", "aria-pressed": "false" }),
    h("button", { "data-action": "no-such-action" }),
  ]);
  G3.wireMapControls(controls, G3.mapControlHandlers({
    stageView: { setCellGrid: (on) => calls.push(["cell-grid", on]) },
    setChartedOnly: (on) => calls.push(["charted-only", on]),
    setTerritories: (on, btn) => { calls.push(["territories", on, btn.dataset.action]); wanted = on; },
    territoriesWanted: () => wanted,
  }));
  const [only, grid, territories, unknown] = controls.children;
  only.click();
  assert.equal(only.getAttribute("aria-pressed"), "true");
  only.click();
  grid.click();
  assert.equal(grid.getAttribute("aria-pressed"), "true");
  territories.click();
  assert.equal(territories.getAttribute("aria-pressed"), "true");
  unknown.click();
  assert.deepEqual(calls, [["charted-only", true], ["charted-only", false], ["cell-grid", true], ["territories", true, "territories"]]);
});

// --- The drill-down ------------------------------------------------------------

test("opens on the whole galaxy with nothing to go back to", async () => {
  const m = await start();
  assert.deepEqual(m.view.stage(), { win: null });
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

test("the whole galaxy has no grid lines; a window shows the lines between its blocks (MAP.85)", async () => {
  const first = await start();
  assert.deepEqual(gridEdges(first), [0]);
  await drillOnce(first);
  const m = await reopenOnAddress(first);
  assert.equal(m.view.stage().win.centre.m, 81);
  const rest = gridEdges(m);
  assert.ok(rest.length === 1 && rest[0] > 0 && rest[0] <= 0.15, "nearly invisible at rest (MAP.202): " + rest);
  assert.deepEqual(Object.keys(m.host.canvasEl.galaxyLines()).sort(), ["blockEdges", "blocksDrawn", "chartedLines", "hovered"]);
});

test("the whole galaxy turns and zooms in to half its fit, never out past it", async () => {
  const m = await start();
  const fit = lastCamera(m);
  const wheel = (deltaY, deltaMode) => m.view.onWheel(new FakeEvent("wheel", { deltaY, deltaMode: deltaMode || 0 }));
  assert.equal(wheel(100), true, "the wheel zooms");
  assert.ok(Math.abs(lastCamera(m).dist - fit.dist) < 1e-9 * fit.dist, "no further out than the fit");
  // Up to the stage's hand-over; one more tick there takes the block under the pointer (MAP.152), tested below. No wall stops the tick that
  // reaches it (MAP.334): it may carry the camera a tick's share past (a tick of this wheel is 0.78 of the distance).
  for (let n = 0; n < 100 && lastCamera(m).dist > fit.dist * SV.GALAXY_MIN_ZOOM * (1 + 1e-6); n++) wheel(-3, 1);
  const reached = lastCamera(m).dist / fit.dist;
  assert.ok(reached <= SV.GALAXY_MIN_ZOOM * (1 + 1e-6) && reached > SV.GALAXY_MIN_ZOOM * 0.75, "twice as close, by a tick at most more");
  m.view.onPointerDown(new FakeEvent("pointerdown", { clientX: 400, clientY: 300, pointerId: 1, button: 0, pointerType: "mouse" }));
  m.view.onPointerMove(new FakeEvent("pointermove", { clientX: 480, clientY: 260, pointerId: 1, pointerType: "mouse" }));
  m.view.onPointerUp(new FakeEvent("pointerup", { clientX: 480, clientY: 260, pointerId: 1, button: 0, pointerType: "mouse" }));
  const turned = lastCamera(m);
  assert.ok(Math.abs(turned.heading - fit.heading) > 0.05 && turned.tilt > 0.05, "dragging turns and tilts it");
  assert.deepEqual(m.view.stage(), { win: null }, "a drag picks nothing");
});

test("each window opened down the ladder pushes a history entry and a breadcrumb, a window one size down each time", async () => {
  const m = await start();
  const seen = [];
  for (let n = 0; n < 12; n++) {
    if (!(await drillOnce(m))) break;
    seen.push(m.view.stage());
    assert.match(m.win.location.search, /^\?cam=[-\d.e+,]+$/, "the URL holds the camera, not a stage");
    assert.ok(S.parseCamera(new URLSearchParams(m.win.location.search).get("cam")), "a camera that reads back");
    assert.deepEqual(crumbLabels(m), expectedCrumbs(m));
    assert.equal(m.win.history.length, seen.length + 1);
    assert.equal(button(m, "back").disabled, false);
    assert.equal(button(m, "up").disabled, false);
  }
  assert.deepEqual(seen.map((stage) => stage.win.centre.m), [81, 27, 9, 3, 1], "each step goes one level down; sectors are the last");
  assert.ok(seen.every((stage) => stage.win.g === 9));
});

test("Back, Forward, Up and Reset buttons move through the stages", async () => {
  const m = await start();
  wire(m);
  await drillOnce(m);
  const first = m.view.stage();
  const firstSearch = m.win.location.search;
  await drillOnce(m);
  const second = m.view.stage();
  assert.notDeepEqual(first, second);

  button(m, "back").click();
  await arrive(m);
  assert.deepEqual(m.view.stage(), first);
  const sameCamera = (a, b) => {
    const x = S.parseCamera(new URLSearchParams(a).get("cam"));
    const y = S.parseCamera(new URLSearchParams(b).get("cam"));
    return [...x.target, x.dist, ...x.quat].every((value, k) => Math.abs(value - [...y.target, y.dist, ...y.quat][k]) <= 1e-8 * Math.max(1, Math.abs(value)));
  };
  assert.ok(sameCamera(m.win.location.search, firstSearch), "back to the entry's own address");
  assert.equal(button(m, "forward").disabled, false);

  button(m, "forward").click();
  await arrive(m);
  assert.deepEqual(m.view.stage(), second);
  assert.equal(button(m, "forward").disabled, true);

  const outline = S.galaxyOutline(F.edgePc, F.shape, F.galaxyRadius);
  button(m, "up").click();
  await arrive(m);
  assert.deepEqual(m.view.stage(), S.parentStage(second, outline, F.edgePc), "up goes to the parent stage");
  assert.equal(m.view.stage().win.centre.m, 243 / 3, "the window one size up, centred on the block holding the middle");
  assert.equal(button(m, "forward").disabled, true, "a new step clears Forward");

  button(m, "reset").click();
  await arrive(m);
  assert.deepEqual(m.view.stage(), { win: null });
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
  const crumbButtons = m.els.crumbs.querySelectorAll("button.crumb");
  assert.equal(crumbButtons.length, 4, "the galaxy and a window at each level down to the sectors");
  crumbButtons[0].click();
  await arrive(m);
  assert.deepEqual(m.view.stage(), { win: null });
  assert.equal(m.els.crumbs.querySelectorAll('[aria-current="location"]').length, 1);
});

test("Escape goes up and Home goes back to the galaxy", async () => {
  const m = await start();
  const outline = S.galaxyOutline(F.edgePc, F.shape, F.galaxyRadius);
  await drillOnce(m);
  await drillOnce(m);
  const second = m.view.stage();
  assert.ok(key(m, "Escape").defaultPrevented);
  await arrive(m);
  assert.deepEqual(m.view.stage(), S.parentStage(second, outline, F.edgePc));
  await drillOnce(m);
  key(m, "Home");
  await arrive(m);
  assert.deepEqual(m.view.stage(), { win: null });
});

test("Up from a window is the window one size up, centred on the block holding its middle (MAP.150)", async () => {
  const m = await start();
  const outline = S.galaxyOutline(F.edgePc, F.shape, F.galaxyRadius);
  await drillOnce(m);
  const window27 = m.view.stage();
  await drillOnce(m);
  const window3 = m.view.stage();
  m.view.up();
  await arrive(m);
  const up = m.view.stage();
  assert.deepEqual(up.win.centre, S.stripBlock(P.drillParent(window3.win.centre)));
  assert.equal(up.win.g, 9);
  assert.notDeepEqual(up, window27, "recentred on where the click was, not back where the window was");
  m.view.up();
  await arrive(m);
  assert.equal(m.view.stage().win.centre.m, 243);
  m.view.up();
  await arrive(m);
  assert.deepEqual(m.view.stage(), { win: null });
  assert.equal(button(m, "up").disabled, true, "nothing above the galaxy");
  assert.equal(S.parentStage(m.view.stage(), outline, F.edgePc), null);
});

test("a stage URL opens that stage; a bad one opens the galaxy and says why", async () => {
  const deep = await start();
  for (let n = 0; n < 6; n++) await drillOnce(deep);
  const query = deep.win.location.search;

  const m = await start("http://localhost/galaxy" + query);
  assert.deepEqual(m.view.stage(), deep.view.stage());
  assert.deepEqual(crumbLabels(m), crumbLabels(deep));

  const bad = await start("http://localhost/galaxy?cam=nonsense");
  assert.deepEqual(bad.view.stage(), { win: null });
  assert.equal(bad.els.notice.hidden, false);
  assert.match(bad.els.notice.textContent, /There is no camera nonsense/);

  // The window tokens of before MAP.228 are not read.
  const old = await start("http://localhost/galaxy?win=81.1.0.1&g=9");
  assert.deepEqual(old.view.stage(), { win: null });
});

test("?sector= opens the stage around that sector with it selected", async () => {
  const d = { shell: 400, band: 0, wedge: 1000 };
  const designation = S.sectorDesignation(d.shell, d.band, d.wedge);
  const m = await start("http://localhost/galaxy?sector=" + designation);
  const outline = S.galaxyOutline(F.edgePc, F.shape, F.galaxyRadius);
  assert.deepEqual(m.view.stage(), S.sectorStage(d.shell, d.band, d.wedge));
  const crumbs = crumbLabels(m);
  assert.equal(crumbs[crumbs.length - 1], "Sector " + S.blockLabel(d) + " (1)");
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
  assert.match(m.win.location.search, new RegExp("^\\?sector=" + S.sectorDesignation(400, 0, 1000) + "(&cam=[-\\d.e+,]+)?$"));
  assert.equal(m.els.notice.hidden, true);

  input.value = "0/0/99";
  m.els.address.dispatchEvent(new FakeEvent("submit"));
  assert.equal(m.els.notice.textContent, "Band 0 of shell 0 has no wedge 99.");

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
  const sector = { kind: "sector", name: "Bead", sector_id: 7, shell_index: 400, band_index: 0, wedge_index: 1000 };
  const answers = {
    Bead: [sector],
    Nowhere: [],
    Many: [Object.assign({}, sector, { name: "Many A" }), { kind: "system", name: "Many B", sector_id: 8, sector_name: "Far", shell_index: 401, band_index: 1, wedge_index: 10 }],
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
  assert.equal(withoutCam(m.win.location.search), "?sector=" + S.sectorDesignation(400, 0, 1000));

  await ask("Nowhere");
  assert.equal(m.els.notice.textContent, "Nothing is named like Nowhere.");

  await ask("Broken");
  assert.equal(m.els.notice.textContent, "The lookup failed. Please try again shortly.");

  await ask("Many");
  assert.equal(m.els.notice.textContent, "2 names match; pick one.");
  assert.equal(m.els.matches.hidden, false);
  const choices = m.els.matches.querySelectorAll('[role="option"]');
  assert.deepEqual(choices.map((b) => b.textContent), ["Many Asector", "Many Bsystem"]);
  choices[1].click();
  await arrive(m);
  assert.equal(m.els.matches.hidden, true);
  // MAP.336: a system is opened in its sector (no printed ID in this answer, so no system of it).
  assert.equal(withoutCam(m.win.location.search), "?sector=" + S.sectorDesignation(401, 1, 10) + "&open=1");
});

// --- The free view -------------------------------------------------------------

async function freeStage() {
  const first = await start();
  // A window of blocks one size down, seen at the isometric slant.
  await drillOnce(first);
  const m = await reopenOnAddress(first);
  assert.equal(m.view.stage().win.centre.m, 81);
  assert.ok(Math.abs(lastCamera(m).tilt - SV.ISO_TILT) < 1e-9, "a window is seen at a slant");
  assert.equal(button(m, "reset-view").disabled, false);
  return m;
}

test("the wheel zooms a free stage, within 1/8 to 2.5 times its fit", async () => {
  const m = await freeStage();
  const fit = lastCamera(m).dist;
  const wheel = (deltaY, deltaMode) => m.view.onWheel(new FakeEvent("wheel", { deltaY, deltaMode: deltaMode || 0 }));
  assert.equal(wheel(100), true);
  assert.ok(lastCamera(m).dist > fit, "scrolling down zooms out");
  // Up to each limit.
  // MAP.334: no wall stops the tick that reaches a hand-over: it may carry the camera a tick's share past (a tick of 200 is 1.65 times the
  // distance, one of -3 lines 0.78 of it).
  for (let n = 0; n < 50 && lastCamera(m).dist < fit * SV.MAX_ZOOM * (1 - 1e-6); n++) wheel(200);
  assert.ok(lastCamera(m).dist >= fit * SV.MAX_ZOOM * (1 - 1e-6) && lastCamera(m).dist < fit * SV.MAX_ZOOM * 1.7);
  for (let n = 0; n < 100 && lastCamera(m).dist > fit * SV.MIN_ZOOM * (1 + 1e-6); n++) wheel(-3, 1);
  assert.ok(lastCamera(m).dist <= fit * SV.MIN_ZOOM * (1 + 1e-6) && lastCamera(m).dist > fit * SV.MIN_ZOOM * 0.75);
});

test("the wheel zooms toward the point under the pointer, which stays under it (MAP.225)", async () => {
  const m = await freeStage();
  const screenOf = (point) => {
    const ndc = new THREE.Vector3(point[0], point[1], point[2]).project(m.camera);
    return [ndc.x, ndc.y];
  };
  for (const [x, y] of [[560, 210], [260, 380]]) {
    const before = lastCamera(m);
    const anchor = m.view.anchorUnder(x, y);
    assert.ok(anchor, "something is under the pointer");
    const where = screenOf(anchor);
    m.view.onWheel(new FakeEvent("wheel", { deltaY: -120, clientX: x, clientY: y }));
    const after = lastCamera(m);
    assert.ok(after.dist < before.dist, "closer");
    assert.notDeepEqual(after.target, before.target, "the center moved toward the pointer");
    const now = screenOf(anchor);
    assert.ok(Math.abs(now[0] - where[0]) < 1e-6 && Math.abs(now[1] - where[1]) < 1e-6, "the point stays where it is on the screen");
    m.view.onWheel(new FakeEvent("wheel", { deltaY: 120, clientX: x, clientY: y }));
    const back = lastCamera(m);
    const off = Math.hypot(back.target[0] - before.target[0], back.target[1] - before.target[1], back.target[2] - before.target[2]);
    assert.ok(off < 1e-6 * before.dist, "zooming out about it brings the center back: " + off);
  }
});

test("a wheel without a position zooms about the view's middle, as before", async () => {
  const m = await freeStage();
  const before = lastCamera(m);
  m.view.onWheel(new FakeEvent("wheel", { deltaY: -120 }));
  assert.deepEqual(lastCamera(m).target, before.target);
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

test("right-dragging moves a free stage anywhere (MAP.263: no pan clamp)", async () => {
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
  assert.ok(off > 1.5 * fit.dist * Math.tan(half), `moved ${off}: past where the clamp held it`);
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
  // (The address the map was reopened on carries the camera to ten digits.)
  assert.ok(Math.abs(back.dist - fit.dist) < 1e-6 * fit.dist, "dist");
  for (let i = 0; i < 4; i++) assert.ok(Math.abs(Math.abs(back.quaternion[i]) - Math.abs(fit.quaternion[i])) < 1e-6, "quaternion " + i);
  for (let i = 0; i < 3; i++) assert.ok(Math.abs(back.target[i] - fit.target[i]) < 1e-6 * fit.dist, "target " + i);
});

test("a pinch zooms a free stage in and out past the stage's own range, the zoom being the camera's alone (MAP.337)", async () => {
  const m = await freeStage();
  const fit = lastCamera(m).dist;
  const touch = (type, id, x, y) => m.view["onPointer" + type](new FakeEvent("pointer" + type.toLowerCase(), { clientX: x, clientY: y, pointerId: id, button: 0, pointerType: "touch" }));
  touch("Down", 1, 300, 300);
  touch("Down", 2, 500, 300);
  touch("Move", 2, 5000, 300);
  assert.ok(lastCamera(m).dist <= fit * SV.MIN_ZOOM * (1 + 1e-6), "spread apart: in to the stage's closest zoom at least");
  touch("Move", 2, 301, 300);
  assert.ok(lastCamera(m).dist >= fit * SV.MAX_ZOOM * (1 - 1e-6), "pinched together: out to the stage's furthest zoom at least (MAP.337: the zoom is the camera's alone, so no end of range holds it)");
  touch("Up", 2, 301, 300);
  touch("Up", 1, 300, 300);
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

test("the Galaxy Map's blocks take no hover, click, touch or outline (MAP.297)", async () => {
  const m = await start();
  assert.equal(pointOver(m, /^Block /), null, "no block names itself under the pointer");
  assert.deepEqual(outlines(m), [], "nothing is outlined");
  for (const [x, y] of [[400, 300], [300, 250], [500, 350]]) {
    clickAt(m, x, y);
    clickAt(m, x, y, "touch");
    clickAt(m, x, y, "touch");
  }
  await arrive(m);
  assert.deepEqual(m.view.stage(), { win: null }, "a click or tap on the galaxy opens no window");
  assert.deepEqual(outlines(m), [], "and outlines nothing");
});

// The outline lines drawn over the map: [opacity, line count].
function outlines(m) {
  const found = [];
  m.host.scene.traverse((object) => {
    if (object.isLineSegments) found.push(object.material.opacity);
  });
  return found.sort((p, q) => q - p);
}

test("the galaxy's thousands of blocks are one mesh; a window's are one each (MAP.150)", async () => {
  const m = await start();
  const meshes = () => {
    const found = [];
    m.host.scene.traverse((object) => { if (object.isMesh) found.push(object); });
    return found.length;
  };
  const outline = S.galaxyOutline(F.edgePc, F.shape, F.galaxyRadius);
  assert.ok(S.resolveStage({ win: null }, outline, F.edgePc).options.length > 800, "more blocks than a mesh each");
  assert.ok(meshes() <= 2, `the galaxy is one solid and one glass mesh (${meshes()})`);
  await drillOnce(m);
  const window = await reopenOnAddress(m);
  const windowMeshes = [];
  window.host.scene.traverse((object) => { if (object.isMesh) windowMeshes.push(object); });
  assert.ok(windowMeshes.length > 20, "a window has a mesh or two for each block");
});

// The address of a sector whose map opens in place (the sector blocks are not picked on the Galaxy Map since
// MAP.297: a sector opens from its address, the wheel's arrival or a double click on one of its stars).
const SECTOR = { shell: 3, band: 0, wedge: 4 };
const openAddress = (sector) => "http://localhost/galaxy?sector=" + S.sectorDesignation(sector.shell, sector.band, sector.wedge) + "&open=1";

test("a generated sector opens in place from its address", async () => {
  const server = () => ({ children: [], sectors: [Object.assign({ id: 42, name: "Bead", system_count: 3 }, SECTOR)] });
  const m = await start(openAddress(SECTOR), { server });
  assert.deepEqual(m.win.location.assigned, [], "no page load: the sector opens in the map");
  assert.deepEqual(m.calls.scenes, [42]);
  assert.equal(S.parseStageQuery(m.win.location.search).open, true, "the URL says the sector is open: " + m.win.location.search);
  assert.deepEqual(S.parseStageQuery(m.win.location.search).sector, SECTOR);
  assert.match(m.win.location.search, /^\?sector=\w+&open=1$/, "the sector link has no window token");
  assert.equal(m.calls.hidden.length, 1, "the galaxy's own stars leave the cell");
  assert.ok(m.calls.hidden[0].r1 > 0);
  // Up closes it, keeping the sector selected and the cell's stars back.
  m.view.up();
  await arrive(m);
  assert.equal(m.calls.hidden[m.calls.hidden.length - 1], null);
  assert.equal(S.parseStageQuery(m.win.location.search).open, false);
  assert.deepEqual(S.parseStageQuery(m.win.location.search).sector, SECTOR);
});

test("a sector nothing was generated in opens in place with what the scatters left, marked uncharted (MAP.162)", async () => {
  const m = await start(openAddress(SECTOR), { unchartedScene: (data) => data });
  assert.deepEqual(m.calls.uncharted, [SECTOR]);
  assert.deepEqual(m.calls.scenes, [], "there is no sector to fetch");
  assert.deepEqual(m.win.location.assigned, [], "no page load: it opens in the map");
  const query = S.parseStageQuery(m.win.location.search);
  assert.equal(query.open, true, "the URL says it is open: " + m.win.location.search);
  assert.deepEqual(query.sector, SECTOR);
  assert.ok(m.calls.hidden[m.calls.hidden.length - 1].r1 > 0, "the galaxy's own stars leave the cell");
  // The title, the info panel and the frame all carry the mark.
  const crumbs = crumbLabels(m);
  assert.match(crumbs[crumbs.length - 1], /^Sector .* \(uncharted\)$/);
  const cell = m.calls.cells[m.calls.cells.length - 1];
  assert.deepEqual(cell.address, SECTOR);
  assert.equal(cell.uncharted.designation, S.sectorDesignation(SECTOR.shell, SECTOR.band, SECTOR.wedge));
  assert.ok(cell.uncharted.stars > 0);
  const badge = m.win.document.querySelector(".sector-uncharted-badge");
  assert.ok(badge && /Uncharted sector/.test(badge.textContent), "the frame carries a badge");
  // Up closes it: the mark goes with it and the cell is selected again.
  m.view.up();
  await arrive(m);
  assert.equal(m.win.document.querySelector(".sector-uncharted-badge"), null);
  assert.equal(m.calls.hidden[m.calls.hidden.length - 1], null);
  assert.equal(S.parseStageQuery(m.win.location.search).open, false);
  assert.ok(!/uncharted/.test(crumbLabels(m).join(" ")));
});

test("a cell the scatters left empty is only selected, not opened (MAP.162)", async () => {
  const m = await start(openAddress(SECTOR), { unchartedScene: (data) => Object.assign(data, { stars: [], clouds: [] }) });
  assert.deepEqual(m.calls.uncharted, [SECTOR]);
  assert.equal(m.calls.hidden.filter((cell) => cell).length, 0, "the galaxy's own stars stay: nothing opened");
  assert.equal(m.win.document.querySelector(".sector-uncharted-badge"), null);
  const cell = m.calls.cells[m.calls.cells.length - 1];
  assert.deepEqual(cell.address, SECTOR);
  assert.equal(cell.uncharted, null);
  assert.ok(!/uncharted/.test(crumbLabels(m).join(" ")));
});

test("a failed uncharted scene leaves the cell selected (MAP.162)", async () => {
  const m = await start(openAddress(SECTOR), { unchartedScene: () => { throw new Error("HTTP 502"); } });
  assert.equal(m.win.document.querySelector(".sector-uncharted-badge"), null);
  assert.deepEqual(m.calls.cells[m.calls.cells.length - 1].address, SECTOR);
  assert.equal(m.els.notice.hidden, true, "no error notice for what is only an extra");
});

test("?sector=...&open=1 opens an uncharted sector from the address (MAP.162)", async () => {
  const m = await start("http://localhost/galaxy?sector=" + S.sectorDesignation(3, 0, 4) + "&open=1", {
    unchartedScene: (data) => data,
  });
  assert.deepEqual(m.calls.uncharted, [{ shell: 3, band: 0, wedge: 4 }]);
  assert.ok(m.win.document.querySelector(".sector-uncharted-badge"));
});

test("a sector page's map opens its sector, keeps to it and sends another sector to its page (MAP.68)", async () => {
  const pinned = { shell: 3, band: 0, wedge: 4 };
  const server = () => ({ children: [], sectors: [
    Object.assign({ id: 42, name: "Pinned", system_count: 3 }, pinned),
    Object.assign({ id: 43, name: "Beside", system_count: 2 }, { shell: 3, band: 0, wedge: 5 }),
  ] });
  const m = await start("http://localhost/sector/42", { server, pinned });
  assert.deepEqual(m.calls.scenes, [42], "its sector is opened with nothing asked of the URL");
  assert.equal(m.win.location.search, "", "no history entry or query of its own");
  const stage = JSON.stringify(m.view.stage());
  // Up, Reset, Escape and Home leave nothing to go back to.
  m.view.up();
  m.view.home();
  key(m, "Escape");
  key(m, "Home");
  await arrive(m);
  assert.equal(JSON.stringify(m.view.stage()), stage);
  assert.equal(m.calls.hidden.filter((cell) => cell === null).length, 0, "the sector stays open");
  assert.equal(m.win.location.search, "");
  // Another generated sector is that sector's page; the open one is not.
  for (let n = 0; n < 40 && !m.calls.pages.length; n++) {
    key(m, "ArrowRight");
    key(m, "Enter");
  }
  await arrive(m);
  assert.deepEqual(m.calls.pages, [43]);
  assert.deepEqual(m.calls.scenes, [42]);
  assert.equal(m.win.location.search, "");
});

test("kinds of object are hidden per sector, kept in the URL and read back from it (MAP.79)", async () => {
  const pinned = { shell: 3, band: 0, wedge: 4 };
  const server = () => ({ children: [], sectors: [Object.assign({ id: 42, name: "Pinned", system_count: 3 }, pinned)] });
  const m = await start("http://localhost/sector/42?hide=nebula,nonsense", { server, pinned });
  assert.equal(m.view.kindHidden("nebula"), true, "the URL's choice is on the opened sector");
  assert.equal(m.view.kindHidden("nonsense"), false, "a kind that doesn't exist is ignored");
  const kinds = m.view.kinds();
  assert.ok(kinds.some((k) => k.kind === "star") && kinds.some((k) => k.kind === "nebula"), JSON.stringify(kinds));
  m.view.setKindHidden("star", true);
  assert.equal(m.win.location.search, "?hide=nebula,star");
  m.view.setKindHidden("nebula", false);
  assert.equal(m.win.location.search, "?hide=star");
  m.view.setKindHidden("star", false);
  assert.equal(m.win.location.search, "", "nothing hidden, nothing in the URL");
  assert.equal(m.view.kindHidden("star"), false);
});

test("on the Galaxy Map a hidden kind rides every stage URL and stays hidden after Up (MAP.79)", async () => {
  const sector = { shell: 3, band: 0, wedge: 4 };
  const server = () => ({ children: [], sectors: [Object.assign({ id: 42, name: "Bead", system_count: 3 }, sector)] });
  const m = await start("http://localhost/galaxy?sector=" + S.sectorDesignation(3, 0, 4) + "&open=1&hide=star", { server });
  assert.deepEqual(m.calls.scenes, [42]);
  assert.equal(m.view.kindHidden("star"), true);
  m.view.up();
  await arrive(m);
  assert.match(withoutCam(m.win.location.search), /^\?sector=\w+&hide=star$/, m.win.location.search);
  assert.equal(m.view.kindHidden("star"), true, "a kind stays hidden after Up");
});

test("a stage URL names the system opened in place (MAP.125)", () => {
  const CODE = S.sectorDesignation(3, 0, 4);
  const asked = S.parseStageQuery("?sector=" + CODE + "&open=1&system=702");
  assert.equal(asked.system, "702");
  assert.equal(asked.open, true);
  // MAP.196: the printed object ID (hex digits in dashed groups) is read too, not only digits.
  const printed = S.parseStageQuery("?sector=" + CODE + "&open=1&system=0008000000-4000000-00d&object=planet:0008000000-4000000-00e");
  assert.equal(printed.system, "0008000000-4000000-00D");
  assert.equal(printed.object, "planet:0008000000-4000000-00E");
  assert.equal(S.parseStageQuery("?sector=" + CODE + "&open=1&system=FE81000A2B-0000005-000").system, "FE81000A2B-0000005-000");
  assert.equal(S.parseStageQuery("?sector=" + CODE + "&open=1&system=0008000000-4000000-00D-1").system, null);
  assert.equal(S.parseStageQuery("?sector=" + CODE + "&open=1&system=12&object=sector:0008000000-4000000-00D").object, null);
  assert.equal(S.parseStageQuery("?sector=" + CODE + "&open=1&system=x").system, null);
  assert.equal(S.parseStageQuery("?sector=" + CODE + "&open=1").system, null);
  assert.equal(S.parseStageQuery("?sector=" + CODE + "&open=1&system=702&object=planet:12").object, "planet:12");
  assert.equal(S.parseStageQuery("?sector=" + CODE + "&open=1&system=702&object=nope").object, null);
});

test("the breadcrumb and the header name the cell the camera's middle is in (MAP.227)", async () => {
  const d = { shell: 400, band: 0, wedge: 1000 };
  const m = await start("http://localhost/galaxy?sector=" + S.sectorDesignation(d.shell, d.band, d.wedge));
  const outline = S.galaxyOutline(F.edgePc, F.shape, F.galaxyRadius);
  const stage = m.view.stage();
  const here = P.sectorPositionPc(d.shell, d.band, d.wedge, F.edgePc);
  assert.equal(S.sectorAtPoint(here, F.edgePc).wedge, d.wedge);
  assert.equal(S.stageAtPoint(stage, here, outline, F.edgePc), stage, "in its own cell the stage stands");
  assert.equal(m.els.where.textContent, "In Sector " + S.blockLabel(d));
  assert.equal(m.els.where.hidden, false);

  const next = { shell: d.shell, band: d.band, wedge: d.wedge + 1 };
  const there = P.sectorPositionPc(next.shell, next.band, next.wedge, F.edgePc);
  const moved = S.stageAtPoint(stage, there, outline, F.edgePc);
  assert.deepEqual(moved.win.centre, { m: 1, ...next }, "the cell the middle moved to");
  m.host.canvasEl.galaxyCentreOn(there[0], there[1], there[2]);
  await arrive(m);
  assert.equal(m.els.where.textContent, "In Sector " + S.blockLabel(next));
  const crumbs = crumbLabels(m);
  assert.equal(crumbs[crumbs.length - 1], "Sector " + S.blockLabel(next) + " (1)", "the trail ends at the cell it is in");
  assert.deepEqual(m.view.stage(), stage, "the data stage is the one it was opened at");
});

test("the address holds the camera: it is written in place once the view is still, and opening it puts the camera there (MAP.228)", async () => {
  const m = await freeStage();
  const entries = m.win.history.length;
  m.view.onWheel(new FakeEvent("wheel", { deltaY: -200, clientX: 560, clientY: 210 }));
  m.view.onPointerDown(new FakeEvent("pointerdown", { clientX: 400, clientY: 300, pointerId: 1, button: 0, pointerType: "mouse" }));
  m.view.onPointerMove(new FakeEvent("pointermove", { clientX: 440, clientY: 320, pointerId: 1, pointerType: "mouse" }));
  m.view.onPointerUp(new FakeEvent("pointerup", { clientX: 440, clientY: 320, pointerId: 1, pointerType: "mouse" }));
  const moved = lastCamera(m);
  await arrive(m);
  assert.equal(m.win.history.length, entries, "the camera replaces the entry, it adds none");
  const camera = S.parseCamera(new URLSearchParams(m.win.location.search).get("cam"));
  assert.ok(camera, m.win.location.search);
  const close = (a, b, scale) => Math.abs(a - b) <= 1e-8 * scale;
  assert.ok(close(camera.dist, moved.dist, moved.dist) && camera.target.every((v, k) => close(v, moved.target[k], moved.dist)));
  assert.ok(camera.quat.every((v, k) => close(v, moved.quaternion[k], 1)), "the turn is in the address too");

  const again = await start("http://localhost/galaxy" + m.win.location.search);
  assert.equal(again.view.stage().win.centre.m, m.view.stage().win.centre.m, "the stage is the one the camera looks at");
  const there = lastCamera(again);
  assert.ok(close(there.dist, moved.dist, moved.dist), "the same distance: " + there.dist + " vs " + moved.dist);
  assert.ok(there.target.every((v, k) => close(v, moved.target[k], moved.dist)), "the same middle");
  assert.ok(there.quaternion.every((v, k) => close(v, moved.quaternion[k], 1)), "the same turn");
});

test("the stage a camera opens is the finest window that holds its view, and the galaxy for a view of the whole (MAP.228)", async () => {
  const m = await freeStage();
  const outline = S.galaxyOutline(F.edgePc, F.shape, F.galaxyRadius);
  assert.ok(outline);
  const whole = await start();
  const galaxy = lastCamera(whole);
  const open = await start("http://localhost/galaxy?cam=" + [...galaxy.target, galaxy.dist, ...galaxy.quaternion].join(","));
  assert.deepEqual(open.view.stage(), { win: null }, "a view of the whole galaxy opens it");

  const close = lastCamera(m);
  const nearer = await start("http://localhost/galaxy?cam=" + [...close.target, close.dist / 3, ...close.quaternion].join(","));
  assert.ok(nearer.view.stage().win.centre.m <= m.view.stage().win.centre.m, "a nearer view opens the same window or a finer one");
});
