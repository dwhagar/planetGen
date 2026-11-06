// html/static/gravitylayer.js
//
// MAP.168: the Sector Map's gravity layer. The sector's 16 x 16 x 16 grid (GET /sector/<id>/gravity, which is API.24's /api/sectors/<id>/gravity) drawn as translucent
// colored zones over one slab at a time or over every slab (the volume), in one of three modes (pull, well depth, tidal strength), with the saddle points marked.
// Design: docs/design/orbits/gravity-map.md; the choice of colors: reports/gravity-range-2026-11-06/.
//
// - The numbers are the TOTAL field (the generated stars and the smooth galaxy, Boss 2026-11-05). The colors are chosen over each sector's own range after the smooth
//   galaxy is held fixed at the sector's center (`field=centered`), from the sector's 2nd to its 98th percentile, on a log scale: the smooth galaxy is 99.9 percent of the
//   pull and 700 to 900 times the stars' wells in the depth, so colors over the plain total are one flat color or one ramp. The legend says so.
// - The ramp is viridis, a color-blind-safe sequential ramp; transparency follows strength and no fill goes above `MAX_ALPHA` (50 percent, as MAP.131) so the stars and
//   routes stay readable; the volume's slabs are fainter (`VOLUME_ALPHA`) so sixteen of them stacked stay under it too.
// - The zone grid is a cube about the sector's center in the galaxy's axes (x, y, z) in parsecs; the scene's y runs the other way (`flipY`), so the picture of a slab is the
//   grid turned over to match the stars drawn with it. Slabs step along z.
//
// The first half is plain functions (tests/js/gravitylayer.test.mjs); `createGravityController` is the drawing half, a part of the Sector Map's stage (galaxysector.js).
// No DOM at import. Plain DOM calls only, never markup from strings.

export const MODES = [["pull", "Pull"], ["depth", "Well depth"], ["tidal", "Tidal strength"]];
export const MODE_NAMES = MODES.map(function (pair) { return pair[0]; });
export const UNIT_TEXT = { pull: "m/s²", depth: "m²/s²", tidal: "s⁻²" };
// Viridis, dark to bright.
export const RAMP = [[68, 1, 84], [59, 82, 139], [33, 145, 140], [94, 201, 98], [253, 231, 37]];
export const MAX_ALPHA = 0.5;
export const VOLUME_ALPHA = 0.04;
export const LOW_PERCENTILE = 2;
export const HIGH_PERCENTILE = 98;
// The colors span at most this many decades below the top (a depth starts at 0, which has no log).
export const DECADES = 4;
export const RETRY_MS = 3000;
export const MAX_TRIES = 60;
export const DEFAULT_SLAB = 8;

// The color of position `t` (0 to 1) on the ramp, [r, g, b].
export function rampColor(t) {
  const x = Math.min(1, Math.max(0, t)) * (RAMP.length - 1);
  const i = Math.min(RAMP.length - 2, Math.floor(x));
  const f = x - i;
  return RAMP[i].map(function (v, n) { return Math.round(v + (RAMP[i + 1][n] - v) * f); });
}

// The `p`th percentile (0 to 100) of `values`, by linear interpolation between the sorted values.
export function percentile(values, p) {
  const sorted = Array.prototype.slice.call(values).sort(function (a, b) { return a - b; });
  if (!sorted.length) return 0;
  const x = (p / 100) * (sorted.length - 1);
  const i = Math.floor(x);
  const f = x - i;
  return i + 1 < sorted.length ? sorted[i] + (sorted[i + 1] - sorted[i]) * f : sorted[i];
}

// The colors' range for a grid's numbers: {lo, hi, flat}. `flat` is true when the sector has no variation to show (every zone reads the same to within a
// part in a million), and then nothing is drawn.
export function colorRange(values) {
  const lo = percentile(values, LOW_PERCENTILE);
  const hi = percentile(values, HIGH_PERCENTILE);
  const flat = !(hi > lo * (1 + 1e-6) + (lo === 0 ? 1e-300 : 0));
  return { lo: lo, hi: hi, flat: flat };
}

// Position on the ramp (0 to 1) of `value` over `range`, on a log scale from `lo` (at least `DECADES` decades under `hi`) to `hi`; below the range 0, above 1.
export function shade(value, range) {
  if (range.flat || !(range.hi > 0)) return 0;
  const lo = Math.max(range.lo, range.hi * Math.pow(10, -DECADES));
  if (!(value > lo)) return 0;
  return Math.min(1, (Math.log(value) - Math.log(lo)) / (Math.log(range.hi) - Math.log(lo)));
}

// Opacity of a zone at ramp position `t`: nothing at the bottom of the range, `cap` at the top.
export function alphaOf(t, cap) {
  return (cap == null ? MAX_ALPHA : cap) * Math.pow(Math.min(1, Math.max(0, t)), 0.8);
}

// The zone at grid position (i, j, k): `values[(i * zones + j) * zones + k]` (i along x).
export function zoneIndex(zones, i, j, k) {
  return (i * zones + j) * zones + k;
}

// The pixels (RGBA bytes, `zones` x `zones`, row 0 at the bottom of the picture) of slab `k`. A pixel's column is the zone's i; its row is j, or `zones - 1 - j` when the
// scene runs y the other way (`flipY`).
export function slabPixels(values, zones, k, range, options) {
  const o = options || {};
  const cap = o.cap == null ? MAX_ALPHA : o.cap;
  const pixels = new Uint8ClampedArray(zones * zones * 4);
  for (let i = 0; i < zones; i++) {
    for (let j = 0; j < zones; j++) {
      const t = shade(values[zoneIndex(zones, i, j, k)], range);
      const color = rampColor(t);
      const row = o.flipY ? zones - 1 - j : j;
      const at = (row * zones + i) * 4;
      pixels[at] = color[0];
      pixels[at + 1] = color[1];
      pixels[at + 2] = color[2];
      pixels[at + 3] = Math.round(alphaOf(t, cap) * 255);
    }
  }
  return pixels;
}

// A slab number kept to 0 .. zones - 1 after a step of `delta`.
export function stepSlab(slab, delta, zones) {
  return Math.min(zones - 1, Math.max(0, slab + delta));
}

// Three significant digits and the mode's unit: "1.55e-10 m/s²".
export function valueText(mode, value) {
  if (!(value > 0)) return "0 " + UNIT_TEXT[mode];
  const text = value >= 0.01 && value < 10000 ? String(Number(value.toPrecision(3))) : value.toExponential(2).replace("e+", "e");
  return text + " " + UNIT_TEXT[mode];
}

// A saddle's place in the scene's frame, parsecs from the sector's center: the galaxy's x and z, and y turned over when the scene runs y the other way.
export function localPoint(position, center, flipY) {
  return [position[0] - center[0], (flipY ? -1 : 1) * (position[1] - center[1]), position[2] - center[2]];
}

// The words of the legend line for the layer's state: {mode, range, grid, building, failed, slab, volume}.
export function legendText(state) {
  if (state.failed) return "Gravity: the grid could not be loaded for this sector.";
  if (state.building || !state.grid) return "Gravity: working out this sector's grid. This can take a minute the first time.";
  const label = MODES.find(function (pair) { return pair[0] === state.mode; })[1];
  const smooth = state.grid.smooth || { level: 0, span: 0 };
  const where = state.volume ? "all slabs" : "slab " + (state.slab + 1) + " of " + state.grid.zones;
  if (state.range.flat) return label + " (" + where + "): no variation across this sector, so nothing is drawn.";
  const span = valueText(state.mode, state.range.lo) + " to " + valueText(state.mode, state.range.hi);
  const base = state.mode === "depth"
    ? "Colors show how deep the stars' wells are below the shallowest zone; the galaxy's own slope across the sector (" + valueText(state.mode, smooth.span) + ") is left out."
    : "Total " + label.toLowerCase() + " " + valueText(state.mode, smooth.level) + " from the smooth galaxy at the center; colors show the variation across this sector.";
  return label + " (" + where + "): " + span + ", the sector's own 2% to 98%, log scale. " + base;
}

// The drawing half. `deps`: `fetchGrid(sectorId, mode, field) -> Promise<grid | {building: true}>` (field "centered" or "total"), `flipY`, `onChange()` (the URL and the menu follow),
// `setTimeout`/`clearTimeout` (default the window's), `picker` (mappick.js: the saddle points are a layer of it), `invalidate()` (draw again).
export function createGravityController(THREE, deps) {
  const flipY = !!deps.flipY;
  const later = deps.setTimeout || function (fn, ms) { return setTimeout(fn, ms); };
  const cancel = deps.clearTimeout || function (id) { clearTimeout(id); };
  const state = { on: false, mode: "pull", slab: DEFAULT_SLAB, volume: false, grids: {}, building: false, failed: false };
  let attached = null; // {id, group}
  let drawn = null; // {group, planes: [mesh], saddles: [mesh], zones}
  let pickerLayer = null;
  let token = 0;
  let retry = null;
  let view = null; // the settings' refresh

  function grid() {
    const held = state.grids[state.mode];
    return held && held.centered ? held : null;
  }

  function currentState() {
    const held = grid();
    return {
      mode: state.mode, slab: state.slab, volume: state.volume, grid: held ? held.centered : null,
      range: held ? held.range : null, building: state.building, failed: state.failed,
    };
  }

  function changed() {
    if (view) view();
    if (deps.onChange) deps.onChange();
  }

  function clearDrawn() {
    if (!drawn) return;
    if (drawn.group.parent) drawn.group.parent.remove(drawn.group);
    drawn.planes.concat(drawn.saddles).forEach(function (mesh) {
      if (mesh.material.map) mesh.material.map.dispose();
      mesh.material.dispose();
      mesh.geometry.dispose();
    });
    drawn = null;
    if (pickerLayer && deps.picker) deps.picker.removeLayer(pickerLayer);
    pickerLayer = null;
  }

  function plane(held, k, cap) {
    const g = held.centered;
    const texture = new THREE.DataTexture(slabPixels(held.centered.values, g.zones, k, held.range, { cap: cap, flipY: flipY }), g.zones, g.zones, THREE.RGBAFormat);
    texture.magFilter = THREE.NearestFilter;
    texture.minFilter = THREE.NearestFilter;
    texture.needsUpdate = true;
    const mesh = new THREE.Mesh(new THREE.PlaneGeometry(g.edge_pc, g.edge_pc),
      new THREE.MeshBasicMaterial({ map: texture, transparent: true, depthWrite: false, side: THREE.DoubleSide }));
    mesh.position.set(0, 0, (k + 0.5) * g.zone_edge_pc - g.edge_pc / 2);
    mesh.renderOrder = 1;
    mesh.userData.slab = k;
    return mesh;
  }

  function saddleText(spot, mode) {
    const total = state.grids[mode] && state.grids[mode].total;
    const g = state.grids[mode].centered;
    const at = total ? total.values[zoneIndex(g.zones, spot.zone[0], spot.zone[1], spot.zone[2])] : null;
    return "Saddle point (zone " + spot.zone.map(function (n) { return n + 1; }).join(", ") + ")" + (at != null ? ": total " + valueText(mode, at) : "");
  }

  function draw() {
    clearDrawn();
    if (!attached || !state.on) return;
    const held = grid();
    if (!held || held.range.flat) return;
    const g = held.centered;
    const group = new THREE.Group();
    group.name = "gravity-layer";
    const planes = [];
    if (state.volume) {
      for (let k = 0; k < g.zones; k++) planes.push(plane(held, k, VOLUME_ALPHA));
    } else {
      planes.push(plane(held, Math.min(state.slab, g.zones - 1), MAX_ALPHA));
    }
    planes.forEach(function (mesh) { group.add(mesh); });
    const saddles = (g.saddles || []).map(function (spot) {
      const mesh = new THREE.Mesh(new THREE.OctahedronGeometry(g.zone_edge_pc * 0.35), new THREE.MeshBasicMaterial({ color: 0xffffff, transparent: true, opacity: 0.9, depthWrite: false }));
      const p = localPoint(spot.position_pc, g.center_pc, flipY);
      mesh.position.set(p[0], p[1], p[2]);
      mesh.renderOrder = 2;
      mesh.userData.saddle = spot;
      group.add(mesh);
      return mesh;
    });
    attached.group.add(group);
    drawn = { group: group, planes: planes, saddles: saddles, zones: g.zones };
    if (deps.picker && saddles.length) {
      pickerLayer = deps.picker.addLayer({
        meshes: function () { return saddles; },
        entryOf: function (hit) { return hit.object.userData.saddle; },
        tooltip: function (spot) { return saddleText(spot, state.mode); },
        priority: 5,
      });
    }
    if (deps.invalidate) deps.invalidate();
  }

  function receive(mode, field, body, mine) {
    if (mine !== token) return false;
    if (body && body.building) return null;
    const held = state.grids[mode] || (state.grids[mode] = {});
    held[field] = body;
    if (field === "centered") held.range = colorRange(body.values);
    return true;
  }

  // Loads the grid of the open sector for the current mode (the colors' field, then the numbers' field), asking again while it is being worked out.
  function load(tries) {
    if (!attached || !state.on || grid() && state.grids[state.mode].total) {
      state.building = false;
      draw();
      changed();
      return;
    }
    const mine = token;
    const mode = state.mode;
    state.building = true;
    state.failed = false;
    const wanted = grid() ? "total" : "centered";
    deps.fetchGrid(attached.id, mode, wanted).then(function (body) {
      const result = receive(mode, wanted, body, mine);
      if (result === false) return;
      if (result === null) {
        if ((tries || 0) >= MAX_TRIES) {
          state.building = false;
          state.failed = true;
          changed();
          return;
        }
        retry = later(function () { retry = null; if (mine === token) load((tries || 0) + 1); }, RETRY_MS);
        changed();
        return;
      }
      load(0);
    }, function () {
      if (mine !== token) return;
      state.building = false;
      state.failed = true;
      changed();
    });
  }

  function start() {
    token += 1;
    if (retry) { cancel(retry); retry = null; }
    load(0);
  }

  const controller = {
    state: state,
    // The sector stage opened a sector: its scene `group` takes the layer's.
    attach: function (info) {
      attached = info;
      state.grids = {};
      state.failed = false;
      if (state.on) start();
    },
    detach: function () {
      token += 1;
      if (retry) { cancel(retry); retry = null; }
      clearDrawn();
      attached = null;
      state.grids = {};
      state.building = false;
      changed();
    },
    setOn: function (on) {
      on = !!on;
      if (on === state.on) return;
      state.on = on;
      if (on) {
        if (attached) start();
      } else {
        token += 1;
        if (retry) { cancel(retry); retry = null; }
        state.building = false;
        clearDrawn();
      }
      changed();
    },
    setMode: function (mode) {
      if (MODE_NAMES.indexOf(mode) < 0 || mode === state.mode) return;
      state.mode = mode;
      if (state.on && attached) start();
      changed();
    },
    setVolume: function (volume) {
      state.volume = !!volume;
      draw();
      changed();
    },
    setSlab: function (slab) {
      const zones = grid() ? grid().centered.zones : 16;
      state.slab = stepSlab(Number(slab) || 0, 0, zones);
      draw();
      changed();
    },
    step: function (delta) {
      const zones = grid() ? grid().centered.zones : 16;
      state.slab = stepSlab(state.slab, delta, zones);
      draw();
      changed();
    },
    isOn: function () { return state.on; },
    // The URL: `gravity=<mode>` while the layer is on (the registry's `show` carries on and off), `gslab=<n>` or `all`.
    readUrl: function (params) {
      const mode = params.get("gravity");
      if (mode && MODE_NAMES.indexOf(mode) >= 0) state.mode = mode;
      const slab = params.get("gslab");
      if (slab === "all") state.volume = true;
      else if (slab && /^\d+$/.test(slab)) { state.slab = Math.max(0, Math.min(15, Number(slab) - 1)); state.volume = false; }
    },
    writeUrl: function (params) {
      if (state.on) {
        params.set("gravity", state.mode);
        params.set("gslab", state.volume ? "all" : String(state.slab + 1));
      } else {
        params.delete("gravity");
        params.delete("gslab");
      }
    },
    legend: function (on) {
      if (!on) return null;
      return { text: legendText(currentState()), swatchClass: "gravity-ramp" };
    },
    // The layer's own controls, under its row of the View menu (maplayers.js `settings`).
    settings: function (container) {
      const doc = container.ownerDocument;
      function button(label, handler, name) {
        const b = doc.createElement("button");
        b.type = "button";
        b.className = "gravity-button";
        b.textContent = label;
        b.addEventListener("click", handler);
        if (name) b.dataset.gravity = name;
        container.appendChild(b);
        return b;
      }
      container.textContent = "";
      const modes = MODES.map(function (pair) {
        const b = button(pair[1], function () { controller.setMode(pair[0]); }, "mode");
        b.dataset.gravityMode = pair[0];
        return b;
      });
      const previous = button("Slab −", function () { controller.step(-1); }, "previous");
      const label = doc.createElement("span");
      label.className = "gravity-slab";
      container.appendChild(label);
      const next = button("Slab +", function () { controller.step(1); }, "next");
      const volumeLabel = doc.createElement("label");
      volumeLabel.className = "gravity-volume";
      const box = doc.createElement("input");
      box.type = "checkbox";
      box.dataset.gravity = "volume";
      box.addEventListener("change", function () { controller.setVolume(box.checked); });
      volumeLabel.appendChild(box);
      volumeLabel.appendChild(doc.createTextNode(" Show every slab"));
      container.appendChild(volumeLabel);
      view = function () {
        modes.forEach(function (b) { b.setAttribute("aria-pressed", b.dataset.gravityMode === state.mode ? "true" : "false"); });
        const zones = grid() ? grid().centered.zones : 16;
        label.textContent = state.volume ? "all slabs" : "slab " + (state.slab + 1) + " of " + zones;
        previous.disabled = state.volume || state.slab <= 0;
        next.disabled = state.volume || state.slab >= zones - 1;
        box.checked = state.volume;
      };
      view();
    },
  };
  return controller;
}
