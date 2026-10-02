// html/static/mapcore.js
//
// The helpers the Galaxy Map (galaxymap3d.js), the Sector Map
// (sectormap.js) and the System Map (systemmap.js) share (MAP.63), so
// each is written once: reading the page's scene JSON, theme colors, the
// info panel's fields, a sector's address, the highlight ring's texture,
// the scale bar's nice numbers and screen span, fitting the renderer to
// its canvas, and picking a point of light on screen. No visible change
// from the copies they replace.

// Imported with this module's own `?v=` query, as systemmap.js explains.
const VERSION_QUERY = new URL(import.meta.url).search;
const THREE = await import(`./vendor/three.module.min.js${VERSION_QUERY}`);

// The scene data a panel's `<script type="application/json">` block
// carries, or null when it is missing or not JSON.
export function readSceneData(dataEl) {
  if (!dataEl) {
    return null;
  }
  try {
    return JSON.parse(dataEl.textContent);
  } catch (err) {
    return null;
  }
}

// A CSS custom property on the page root (the theme's colors), or
// `fallback` when it isn't set.
export function cssVar(name, fallback) {
  var value = getComputedStyle(document.documentElement).getPropertyValue(name).trim();
  return value || fallback;
}

// Whether the maps' background (--bg-subtle, the light or dark theme) is
// light, from its computed color.
export function isLightBackground() {
  var probe = document.createElement("span");
  probe.style.color = cssVar("--bg-subtle", "#000");
  document.body.appendChild(probe);
  var rgb = (getComputedStyle(probe).color.match(/[\d.]+/g) || [0, 0, 0]).map(Number);
  probe.remove();
  return 0.2126 * rgb[0] + 0.7152 * rgb[1] + 0.0722 * rgb[2] > 140;
}

// One "label: value" row of an info panel's <dl>, left out when there's
// no value (0 is a value). Plain text, never markup: names are database
// content.
export function addField(dl, label, value) {
  if (!value && value !== 0) {
    return;
  }
  var dt = document.createElement("dt");
  dt.textContent = label;
  var dd = document.createElement("dd");
  dd.textContent = value;
  dl.appendChild(dt);
  dl.appendChild(dd);
}

// A sector's grid address as the info panels show it.
export function formatAddress(ringIndex, layerIndex, slotIndex) {
  return "ring " + ringIndex + " layer " + layerIndex + " slot " + slotIndex;
}

// A thin ring in `color` (the selection highlight, a marked rogue
// planet), as a sprite texture.
export function makeRingTexture(color) {
  var size = 64;
  var canvasEl = document.createElement("canvas");
  canvasEl.width = canvasEl.height = size;
  var ctx = canvasEl.getContext("2d");
  var r = size / 2;
  ctx.beginPath();
  ctx.arc(r, r, r - 4, 0, Math.PI * 2);
  ctx.lineWidth = 3;
  ctx.strokeStyle = color;
  ctx.stroke();
  var texture = new THREE.CanvasTexture(canvasEl);
  texture.colorSpace = THREE.SRGBColorSpace;
  return texture;
}

// --- Scale bar -------------------------------------------------------------------

// Snaps a positive value to the nearest "nice" 1, 2 or 5 x 10^n, the usual
// map-scale convention, so a label reads "5 ly" or "20 ly" rather than
// "6.283 ly"; 0 for anything not a positive number.
export function niceScaleValue(raw) {
  if (!isFinite(raw) || raw <= 0) {
    return 0;
  }
  var magnitude = Math.pow(10, Math.floor(Math.log10(raw)));
  var mantissa = raw / magnitude;
  var niceMantissa;
  if (mantissa < 1.5) niceMantissa = 1;
  else if (mantissa < 3.5) niceMantissa = 2;
  else if (mantissa < 7.5) niceMantissa = 5;
  else niceMantissa = 10;
  return niceMantissa * magnitude;
}

// How many world units one screen pixel spans at `distance` from a
// perspective camera, on a canvas `heightPx` tall.
export function worldUnitsPerPixel(camera, distance, heightPx) {
  return (2 * distance * Math.tan(THREE.MathUtils.degToRad(camera.fov) / 2)) / (heightPx || 1);
}

// --- Canvas size -----------------------------------------------------------------

// Sizes the renderer's drawing buffer to the canvas as laid out and the
// camera's aspect to match.
export function fitRendererToCanvas(renderer, camera, canvasEl) {
  var width = canvasEl.clientWidth || 1;
  var height = canvasEl.clientHeight || 1;
  renderer.setSize(width, height, false);
  camera.aspect = width / height;
  camera.updateProjectionMatrix();
}

// Calls `resize` now and whenever `target` (or the window) changes size.
export function watchResize(target, resize) {
  if (typeof ResizeObserver !== "undefined" && target) {
    new ResizeObserver(resize).observe(target);
  }
  resize();
  window.addEventListener("resize", resize);
}

// --- Picking on screen -----------------------------------------------------------

// The entry (anything with x, y, z) drawn nearest a screen point, within
// `options.reach(entry)` pixels of its center, as {entry, px}, or null.
// Points of light are a fixed number of pixels across at any distance, so
// they are picked on screen rather than by raycast. `options.accept(entry)`
// leaves an entry out when it returns false; behind the camera or past
// its far plane never counts. Between two at exactly the same distance
// the first wins, or the last with `options.lastWins`.
var projected = new THREE.Vector3();
export function nearestOnScreen(entries, camera, rect, clientX, clientY, options) {
  if (!rect.width || !rect.height) {
    return null;
  }
  var best = null;
  var bestPx = Infinity;
  entries.forEach(function (entry) {
    if (options.accept && !options.accept(entry)) {
      return;
    }
    projected.set(entry.x, entry.y, entry.z).project(camera);
    if (projected.z < -1 || projected.z > 1) {
      return;
    }
    var px = Math.hypot(
      rect.left + ((projected.x + 1) / 2) * rect.width - clientX,
      rect.top + ((1 - projected.y) / 2) * rect.height - clientY);
    if (px <= options.reach(entry) && (px < bestPx || (options.lastWins && px === bestPx))) {
      best = entry;
      bestPx = px;
    }
  });
  return best ? { entry: best, px: bestPx } : null;
}
