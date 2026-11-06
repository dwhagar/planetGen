// html/static/systemgravity.js
//
// The System Map's Gravity layer (MAP.169; docs/design/orbits/gravity-map.md): the orbital plane of the scene on show as a heat map of the pull, the well depth or
// the tidal strength of the bodies drawn in it, the Hill sphere of each body outlined, and the five Lagrange points of each body and what it orbits marked.
//
// The System Map is a true-position diagram on a shared log radial scale (planetgen/web/maps/systemmap.py `_radial_px`; the scene carries its bounds as
// data-lokm/data-hikm), or on one linear scale in a moon system (data-pxperkm). This layer works in kilometres and maps every point the way the markers were
// placed: a pixel of the picture is turned back into the kilometres it stands for (`pixelToKm`), the field is worked out there, and a Hill sphere or a point
// is sampled in kilometres and mapped forward (`kmToPixel`). So the heat map shows the field where the body's own markers say it is; on the log scale it is
// drawn between the innermost and the outermost orbit only, where the scale stands for a distance.
//
// The numbers: the bodies are point masses (softened by their own radius, so a zone holding one stays finite) at the places their markers stand for, G in SI.
// The colors are the Sector Map layer's (gravitylayer.js): the scene's own 2nd to 98th percentile, log scale, one ramp. The pure functions are exported for
// tests/js/systemgravity.test.mjs; the controller at the end draws.

import { MODES, MODE_NAMES, UNIT_TEXT, alphaOf, colorRange, rampColor, shade, valueText } from "./gravitylayer.js";

export const BIG_G = 6.6743e-11;
export const GRID = 140; // cells across the scene's frame
export const HILL_STEPS = 72;
export const MAX_ALPHA = 0.55;
export const POINT_NAMES = ["L1", "L2", "L3", "L4", "L5"];
const SVG_NS = "http://www.w3.org/2000/svg";

// The scene's scale from its body layer's attributes: {c, lo, hi, min, spread, perKm, ok}.
export function sceneScale(attrs) {
  const number = function (name) { return parseFloat(attrs[name]); };
  const scale = {
    lo: number("lokm"), hi: number("hikm"), c: number("cpx"), perKm: number("pxperkm"), min: number("minpx"), spread: number("spreadpx"),
  };
  scale.ok = isFinite(scale.c) && (scale.perKm > 0 || (scale.lo > 0 && scale.hi > 0 && isFinite(scale.min) && isFinite(scale.spread)));
  return scale;
}

// Kilometres (x east, y north, from the scene's middle) to the picture's pixel: the log radial scale keeps the real angle and maps the distance.
export function kmToPixel(scale, x, y) {
  const d = Math.hypot(x, y);
  if (d <= 0) return { x: scale.c, y: scale.c };
  if (scale.perKm > 0) return { x: scale.c + x * scale.perKm, y: scale.c - y * scale.perKm };
  let r;
  if (scale.hi <= scale.lo) {
    r = scale.min + scale.spread;
  } else {
    const clamped = Math.max(scale.lo, Math.min(scale.hi, d));
    r = scale.min + (Math.log10(clamped) - Math.log10(scale.lo)) / (Math.log10(scale.hi) - Math.log10(scale.lo)) * scale.spread;
  }
  return { x: scale.c + r * x / d, y: scale.c - r * y / d };
}

// The pixel back to kilometres, or null where the scale stands for no distance (inside the innermost orbit's disk, past the outermost).
export function pixelToKm(scale, px, py) {
  const dx = px - scale.c;
  const dy = scale.c - py;
  if (scale.perKm > 0) return { x: dx / scale.perKm, y: dy / scale.perKm };
  const r = Math.hypot(dx, dy);
  if (scale.hi <= scale.lo || r < scale.min || r > scale.min + scale.spread) return null;
  const d = Math.pow(10, Math.log10(scale.lo) + (r - scale.min) / scale.spread * (Math.log10(scale.hi) - Math.log10(scale.lo)));
  return r > 0 ? { x: dx / r * d, y: dy / r * d } : { x: 0, y: 0 };
}

// The three modes at kilometres (x, y) in the orbital plane from `sources` [{x, y, mass, radius}] (km, kg, km): the pull in m/s^2 (the size of the acceleration),
// the potential in m^2/s^2 (negative; the depth is measured from the shallowest cell), or the tidal strength in 1/s^2 (the largest eigenvalue of d a / d x, the
// third axis being the plane's normal, where each body only squeezes).
export function fieldAt(sources, x, y, mode) {
  let ax = 0, ay = 0, potential = 0, txx = 0, tyy = 0, txy = 0, tzz = 0;
  for (let n = 0; n < sources.length; n++) {
    const s = sources[n];
    const dx = (x - s.x) * 1000;
    const dy = (y - s.y) * 1000;
    const eps = Math.max(s.radius || 0, 1) * 1000;
    const rho2 = dx * dx + dy * dy + eps * eps;
    const rho = Math.sqrt(rho2);
    const gm = BIG_G * s.mass;
    if (mode === "tidal") {
      const inv3 = gm / (rho2 * rho);
      const inv5 = gm / (rho2 * rho2 * rho);
      txx += 3 * dx * dx * inv5 - inv3;
      tyy += 3 * dy * dy * inv5 - inv3;
      txy += 3 * dx * dy * inv5;
      tzz += -inv3;
    } else if (mode === "depth") {
      potential += -gm / rho;
    } else {
      ax += -gm * dx / (rho2 * rho);
      ay += -gm * dy / (rho2 * rho);
    }
  }
  if (mode === "tidal") {
    const mean = (txx + tyy) / 2;
    const half = Math.sqrt(Math.pow((txx - tyy) / 2, 2) + txy * txy);
    return Math.max(mean + half, tzz);
  }
  if (mode === "depth") return potential;
  return Math.hypot(ax, ay);
}

// The heat map's numbers: a GRID x GRID array over the scene's frame (`bounds` [lo, size], the viewBox), row 0 the top, NaN where the scale stands for no
// distance, and the 'depth' measured from the shallowest cell (a well's depth, not the potential).
export function fieldGrid(scale, bounds, sources, mode, grid) {
  const n = grid || GRID;
  const values = new Float64Array(n * n);
  const step = bounds[1] / n;
  let top = -Infinity;
  for (let row = 0; row < n; row++) {
    for (let col = 0; col < n; col++) {
      const at = pixelToKm(scale, bounds[0] + (col + 0.5) * step, bounds[0] + (row + 0.5) * step);
      const v = at ? fieldAt(sources, at.x, at.y, mode) : NaN;
      values[row * n + col] = v;
      if (v > top) top = v;
    }
  }
  if (mode === "depth") {
    for (let i = 0; i < values.length; i++) if (!isNaN(values[i])) values[i] = top - values[i];
  }
  return values;
}

// The picture's RGBA bytes for `values` and the range the colors use; the NaN cells and the flat scene stay clear.
export function heatPixels(values, range) {
  const pixels = new Uint8ClampedArray(values.length * 4);
  if (range.flat) return pixels;
  for (let i = 0; i < values.length; i++) {
    if (isNaN(values[i])) continue;
    const t = shade(values[i], range);
    const color = rampColor(t);
    pixels[i * 4] = color[0];
    pixels[i * 4 + 1] = color[1];
    pixels[i * 4 + 2] = color[2];
    pixels[i * 4 + 3] = Math.round(alphaOf(t, MAX_ALPHA) * 255);
  }
  return pixels;
}

// The Hill radius of a body of `mass` at `distance` from what it orbits (`host` mass): d (m / 3M)^(1/3), the same unit as the distance.
export function hillRadius(distance, mass, host) {
  if (!(distance > 0) || !(mass > 0) || !(host > 0)) return 0;
  return distance * Math.cbrt(mass / (3 * host));
}

// The x of L1, L2 and L3 in units of the separation, the barycenter at 0, the host at -mu and the body at 1 - mu (the roots of the force balance, by bisection).
export function collinearOffsets(mu) {
  const force = function (x) { return x - (1 - mu) * (x + mu) / Math.pow(Math.abs(x + mu), 3) - mu * (x - 1 + mu) / Math.pow(Math.abs(x - 1 + mu), 3); };
  const root = function (low, high) {
    let fLow = force(low);
    for (let i = 0; i < 200; i++) {
      const mid = 0.5 * (low + high);
      const fMid = force(mid);
      if ((fMid > 0) === (fLow > 0)) { low = mid; fLow = fMid; } else { high = mid; }
    }
    return 0.5 * (low + high);
  };
  const eps = 1e-12;
  return [root(-mu + eps, 1 - mu - eps), root(1 - mu + eps, 3), root(-3, -mu - eps)];
}

// The five Lagrange points of a body at (x, y) km with `mass` about a host at the middle with mass `host`, for a body moving with (vx, vy) (prograde, anticlockwise
// seen from above, when the motion is not known): {name, x, y} in km. Null for a pair too alike in mass for the restricted problem (mu >= 0.5).
export function lagrangePoints(x, y, mass, host, vx, vy) {
  const a = Math.hypot(x, y);
  const mu = mass / (mass + host);
  if (!(a > 0) || !(mu > 0) || !(mu < 0.5)) return null;
  const ux = x / a, uy = y / a;
  let sense = 1; // the orbit's normal along +z (anticlockwise) or -z
  if (isFinite(vx) && isFinite(vy) && (x * vy - y * vx) < 0) sense = -1;
  const wx = -sense * uy, wy = sense * ux;
  const offsets = collinearOffsets(mu);
  const out = [];
  [["L1", offsets[0]], ["L2", offsets[1]], ["L3", offsets[2]]].forEach(function (pair) {
    out.push({ name: pair[0], x: a * pair[1] * ux, y: a * pair[1] * uy });
  });
  [["L4", 1], ["L5", -1]].forEach(function (pair) {
    out.push({ name: pair[0], x: a * ((0.5 - mu) * ux + pair[1] * Math.sqrt(3) / 2 * wx), y: a * ((0.5 - mu) * uy + pair[1] * Math.sqrt(3) / 2 * wy) });
  });
  return out;
}

// The sources and the orbiting bodies of a scene from its markers' attributes `elements` [{kind, id, name, xkm, ykm, masskg, radiuskm, vxkms, vykms, self, el}] and the
// mass `host` the scene's bodies orbit (0 when unknown, then there are no Hill spheres and points). A body with no mass is no source.
export function sceneBodies(elements, host) {
  const sources = [];
  const orbiters = [];
  elements.forEach(function (e) {
    const mass = parseFloat(e.masskg);
    if (!(mass > 0)) return;
    const x = parseFloat(e.xkm) || 0;
    const y = parseFloat(e.ykm) || 0;
    const radius = parseFloat(e.radiuskm) || 0;
    const body = { kind: e.kind, x: x, y: y, mass: mass, radius: radius, name: e.name, vx: parseFloat(e.vxkms), vy: parseFloat(e.vykms), id: e.id, el: e.el };
    sources.push(body);
    if (e.kind !== "star" && !e.self) orbiters.push(body);
  });
  return { sources: sources, orbiters: orbiters, host: host };
}

export function legendText(state) {
  if (state.empty) return "Gravity: this view has no body with a mass to work the field from.";
  const label = MODES.find(function (pair) { return pair[0] === state.mode; })[1];
  if (state.range.flat) return label + ": no variation across this view, so nothing is drawn.";
  const span = valueText(state.mode, state.range.lo) + " to " + valueText(state.mode, state.range.hi);
  return label + ": " + span + ", this view's own 2% to 98%, log scale, from the bodies drawn here (the stars, planets and moons of this view, as points; the galaxy's own pull is too small to count). "
    + "Dashed rings are Hill spheres; diamonds are the Lagrange points L1 to L5 of each body and what it orbits.";
}

function svgElement(doc, name, attrs) {
  const e = doc.createElementNS(SVG_NS, name);
  Object.keys(attrs || {}).forEach(function (k) { e.setAttribute(k, attrs[k]); });
  return e;
}

// The layer: `deps` {root, onChange()}. It draws into the orbit layer of the scene on show (the decorative, aria-hidden svg under the markers).
export function createSystemGravity(deps) {
  const state = { on: false, mode: "pull", range: { lo: 0, hi: 0, flat: true }, empty: false };
  let view = null;

  function visibleLayers() {
    const doc = deps.root.ownerDocument;
    const body = doc.querySelector(".sysmap-svg:not(.sysmap-orbits-layer):not(.sysmap-hidden)");
    if (!body) return null;
    const orbits = Array.prototype.find.call(doc.querySelectorAll(".sysmap-orbits-layer"), function (layer) { return layer.dataset.scene === body.dataset.scene; });
    return orbits ? { body: body, orbits: orbits } : null;
  }

  function clear() {
    deps.root.ownerDocument.querySelectorAll(".sysmap-gravity-group").forEach(function (g) { g.remove(); });
  }

  function elementsOf(body) {
    return Array.prototype.map.call(body.querySelectorAll('[data-kind="star"],[data-kind="planet"],[data-kind="moon"]'), function (el) {
      const d = el.dataset;
      return { kind: d.kind, id: d.id, name: d.name, xkm: d.xkm, ykm: d.ykm, masskg: d.masskg, radiuskm: d.radiuskm, vxkms: d.vxkms, vykms: d.vykms, self: d.self === "true", el: el };
    });
  }

  function draw() {
    clear();
    state.empty = false;
    if (!state.on) { changed(); return; }
    const layers = visibleLayers();
    if (!layers) { changed(); return; }
    const doc = deps.root.ownerDocument;
    const scale = sceneScale(layers.body.dataset);
    const box = layers.body.viewBox.baseVal;
    const bounds = [box.x, box.width];
    const host = parseFloat(layers.body.dataset.hostkg) || 0;
    const parts = sceneBodies(elementsOf(layers.body), host);
    // In a star's scene the star is at the middle and a source like any body; in a moon scene the planet is the self body.
    if (!scale.ok || !parts.sources.length) { state.empty = true; changed(); return; }
    const group = svgElement(doc, "g", { class: "sysmap-gravity-group", "aria-hidden": "true" });
    const values = fieldGrid(scale, bounds, parts.sources, state.mode);
    state.range = colorRange(Array.prototype.filter.call(values, function (v) { return !isNaN(v); }));
    const canvas = doc.createElement("canvas");
    canvas.width = GRID;
    canvas.height = GRID;
    const context = canvas.getContext("2d");
    if (context) {
      context.putImageData(new ImageData(heatPixels(values, state.range), GRID, GRID), 0, 0);
      group.appendChild(svgElement(doc, "image", {
        class: "sysmap-gravity-image", x: String(bounds[0]), y: String(bounds[0]), width: String(bounds[1]), height: String(bounds[1]),
        href: canvas.toDataURL("image/png"), preserveAspectRatio: "none",
      }));
    }
    if (host > 0) drawBodies(doc, group, scale, parts);
    layers.orbits.insertBefore(group, layers.orbits.firstChild);
    changed();
  }

  // The Hill sphere and the Lagrange points of each body that orbits something, in the places the body's own marker stands (a marker nudged for room takes its ring and points along).
  function drawBodies(doc, group, scale, parts) {
    parts.orbiters.forEach(function (b) {
      const d = Math.hypot(b.x, b.y);
      const radius = hillRadius(d, b.mass, parts.host);
      const home = kmToPixel(scale, b.x, b.y);
      const marker = b.el ? deps.markerCenter(b.el) : null;
      const shift = marker ? { x: marker.x - home.x, y: marker.y - home.y } : { x: 0, y: 0 };
      if (radius > 0) {
        const points = [];
        for (let i = 0; i < HILL_STEPS; i++) {
          const angle = 2 * Math.PI * i / HILL_STEPS;
          const p = kmToPixel(scale, b.x + radius * Math.cos(angle), b.y + radius * Math.sin(angle));
          points.push((p.x + shift.x).toFixed(1) + "," + (p.y + shift.y).toFixed(1));
        }
        const ring = svgElement(doc, "polygon", { class: "sysmap-hill", points: points.join(" ") });
        const title = svgElement(doc, "title", {});
        title.textContent = "Hill sphere of " + (b.name || "this body") + ": " + Math.round(radius).toLocaleString("en-US") + " km";
        ring.appendChild(title);
        group.appendChild(ring);
      }
      const found = lagrangePoints(b.x, b.y, b.mass, parts.host, b.vx, b.vy);
      if (!found) return;
      found.forEach(function (point) {
        const p = kmToPixel(scale, point.x, point.y);
        const x = p.x + shift.x, y = p.y + shift.y;
        const mark = svgElement(doc, "path", { class: "sysmap-lagrange", d: "M" + x.toFixed(1) + "," + (y - 4).toFixed(1) + " l4,4 l-4,4 l-4,-4 z", "data-point": point.name });
        const title = svgElement(doc, "title", {});
        title.textContent = (b.name || "Body") + "-" + point.name;
        mark.appendChild(title);
        group.appendChild(mark);
        const label = svgElement(doc, "text", { class: "sysmap-lagrange-label", x: (x + 6).toFixed(1), y: (y + 3).toFixed(1) });
        label.textContent = point.name;
        group.appendChild(label);
      });
    });
  }

  function changed() {
    if (view) view();
    if (deps.onChange) deps.onChange();
  }

  const controller = {
    state: state,
    setOn: function (on) {
      on = !!on;
      if (on === state.on) return;
      state.on = on;
      draw();
    },
    setMode: function (mode) {
      if (MODE_NAMES.indexOf(mode) < 0 || mode === state.mode) return;
      state.mode = mode;
      if (state.on) draw(); else changed();
    },
    // A scene was shown (the map switched between the system and a moon system, or between a moon system's sizes).
    sceneShown: function () { if (state.on) draw(); },
    readUrl: function (params) {
      const mode = params.get("gravity");
      if (mode && MODE_NAMES.indexOf(mode) >= 0) state.mode = mode;
    },
    writeUrl: function (params) {
      if (state.on) params.set("gravity", state.mode); else params.delete("gravity");
    },
    legend: function (on) {
      if (!on) return null;
      return { text: legendText({ mode: state.mode, range: state.range, empty: state.empty }), swatchClass: "gravity-ramp" };
    },
    settings: function (container) {
      const doc = container.ownerDocument;
      container.textContent = "";
      const modes = MODES.map(function (pair) {
        const b = doc.createElement("button");
        b.type = "button";
        b.className = "gravity-button";
        b.textContent = pair[1];
        b.dataset.gravity = "mode";
        b.dataset.gravityMode = pair[0];
        b.addEventListener("click", function () { controller.setMode(pair[0]); });
        container.appendChild(b);
        return b;
      });
      view = function () {
        modes.forEach(function (b) { b.setAttribute("aria-pressed", b.dataset.gravityMode === state.mode ? "true" : "false"); });
      };
      view();
    },
    unitText: UNIT_TEXT,
  };
  return controller;
}
