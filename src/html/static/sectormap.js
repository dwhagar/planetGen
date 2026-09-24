// html/static/sectormap.js
//
// Renders the 3D sector map built by `lib/starmap.py` as a real WebGL
// scene (three.js, vendored at `static/vendor/three.module.min.js` -- see
// that directory's `THIRD_PARTY_NOTICES.txt` for why it's vendored rather
// than loaded from a CDN) instead of the CSS `transform-style:
// preserve-3d` scene this file used to drive directly. `#starmap-data`
// (a `<script type="application/json">` block `starmap.py` writes) is the
// only thing read from the page -- every position/size/color/label for
// every star and phenomenon cloud, plus the compass arrow, is data
// `starmap.py` already computed; this file only ever turns that data
// into real 3D bodies/lines and wires up drag-to-rotate, scroll/button-
// to-zoom, and click/keyboard-for-info, the same interaction set the old
// CSS version had (a real perspective camera now does the projection/
// occlusion a browser's `preserve-3d` compositor used to).
//
// Every star and phenomenon cloud (nebula/asteroid field/black hole/
// neutron star) is a real textured, glowing `THREE.Mesh` sphere -- not a
// flat camera-facing `THREE.Sprite` the way this file used to draw them
// (a label and the highlight ring still are sprites, which suits them
// fine: flat 2D content with no "seen from any angle" concern a body's
// own true 3D shape has). `./bodyRendering.js` supplies the shared
// glow-shell shader and star granulation texture, the same building
// blocks `static/systemmap.js` uses for its own per-marker spheres.
//
// Built with plain DOM calls (never innerHTML/textContent-with-markup)
// when filling the info panel, same discipline the old version had --
// every star/system/phenomenon name is still database content (a system
// name can contain arbitrary characters via `--name`).

import * as THREE from "./vendor/three.module.min.js";
import { makeGlowMaterial, makeStarSurfaceTexture } from "./bodyRendering.js";

var canvas = document.getElementById("starmap-canvas");
var dataEl = document.getElementById("starmap-data");

function readSceneData() {
  if (!dataEl) {
    return null;
  }
  try {
    return JSON.parse(dataEl.textContent);
  } catch (err) {
    return null;
  }
}

var sceneData = readSceneData();

function cssVar(name, fallback) {
  var value = getComputedStyle(document.documentElement).getPropertyValue(name).trim();
  return value || fallback;
}

function addField(dl, label, value) {
  if (!value) {
    return;
  }
  var dt = document.createElement("dt");
  dt.textContent = label;
  var dd = document.createElement("dd");
  dd.textContent = value;
  dl.appendChild(dt);
  dl.appendChild(dd);
}

// A cloud entry carries `kind` (its phenomenon-texture recipe, see
// `CLOUD_KIND_RECIPES` below); a star entry never does; a neighboring-
// sector indicator carries `isNeighbor` -- that alone is enough to tell
// all three apart, unlike the old version's explicit
// `data-kind="phenomenon"` marker.
function showObjectInfo(entry) {
  var panel = document.getElementById("starmap-info");
  if (!panel || !entry) {
    return;
  }
  panel.textContent = "";

  if (entry.isNeighbor) {
    showNeighborInfo(panel, entry);
    return;
  }

  var heading = document.createElement("h3");
  heading.textContent = entry.name || "Unknown";
  panel.appendChild(heading);

  var dl = document.createElement("dl");
  if (entry.kind) {
    addField(dl, "Type", entry.typeLabel);
    addField(dl, "Radius", entry.radiusText);
    addField(dl, "Distance", entry.distanceText);
    panel.appendChild(dl);
    panel.appendChild(navLink(entry, "View phenomenon →"));
    return;
  }
  addField(dl, "Star type", entry.starType);
  addField(dl, "Temperature", entry.temp);
  addField(dl, "Octant", entry.quadrant);
  addField(dl, "Location", entry.location);
  panel.appendChild(dl);
  panel.appendChild(navLink(entry, "View system →"));
}

// A neighboring sector's own address -- shared display convention with
// `static/galaxymap3d.js`'s identically-named helpers for its own
// "planned" (not-yet-generated) sector addresses, since this is the same
// underlying concept one level in: a `(shell_index, shell_slot_index)`
// address, generated or not.
function formatAddress(shellIndex, slotIndex) {
  return "shell " + shellIndex + " slot " + slotIndex;
}

function cliSnippet(shellIndex, slotIndex) {
  return "generate.py galaxy --shell " + shellIndex + " --slot " + slotIndex;
}

function makeCopyButton(text) {
  var button = document.createElement("button");
  button.type = "button";
  button.className = "btn";
  button.textContent = "Copy CLI command";
  button.addEventListener("click", function () {
    var restore = button.textContent;
    var onDone = function () {
      button.textContent = "Copied!";
      setTimeout(function () {
        button.textContent = restore;
      }, 1500);
    };
    var onFail = function () {
      // Clipboard API unavailable (insecure context, permissions, older
      // browser) -- fall back to a selectable readonly field the visitor
      // can copy by hand, rather than silently doing nothing.
      var input = document.createElement("input");
      input.type = "text";
      input.readOnly = true;
      input.value = text;
      input.className = "galaxymap3d-address-field";
      button.insertAdjacentElement("afterend", input);
      input.focus();
      input.select();
    };
    if (navigator.clipboard && navigator.clipboard.writeText) {
      navigator.clipboard.writeText(text).then(onDone, onFail);
    } else {
      onFail();
    }
  });
  return button;
}

// A neighboring sector that already exists just links straight to it
// (same as a star/cloud entry); one that doesn't yet shows its address
// and a copyable `generate.py galaxy --shell K --slot N` command instead
// -- the same "not yet generated" info panel shape
// `static/galaxymap3d.js`'s own `showPlannedInfo` already uses for its
// "planned" tier, one galaxy-map zoom level up from this sector-level view.
function showNeighborInfo(panel, entry) {
  var heading = document.createElement("h3");
  heading.textContent = entry.exists ? entry.name || "Unnamed sector" : "Not yet generated";
  panel.appendChild(heading);

  var dl = document.createElement("dl");
  addField(dl, "Address", formatAddress(entry.shellIndex, entry.shellSlotIndex));
  addField(dl, "Designation", entry.designation);
  panel.appendChild(dl);

  if (entry.exists) {
    panel.appendChild(navLink(entry, "View sector →"));
    return;
  }

  var code = document.createElement("code");
  code.className = "galaxymap3d-cli-snippet";
  code.textContent = cliSnippet(entry.shellIndex, entry.shellSlotIndex);
  panel.appendChild(code);
  panel.appendChild(makeCopyButton(cliSnippet(entry.shellIndex, entry.shellSlotIndex)));
}

// A real, focusable `<a>` carrying `data-nav-target`/`data-nav-params`
// instead of an `href` query string -- `static/navform.js`'s document-
// level click handler (loaded on every page, see `lib/page.py`'s
// `render`) is what actually follows it, by posting a throwaway hidden
// form, the same convention every other in-app link now uses
// (`lib/fmt.py`'s `post_link` builds the non-JS-required `<form>` version
// of the same idea; a `<form>` can't be dynamically inserted into this
// panel's own DOM update flow as conveniently as a plain `<a>` can, so
// this stays in the `data-nav-target` camp like the map's own SVG-era
// markers already had to for the same "can't nest a form" reason).
function navLink(entry, label) {
  var link = document.createElement("a");
  link.href = "#";
  link.className = "btn";
  link.dataset.navTarget = entry.navTarget;
  link.dataset.navParams = JSON.stringify(entry.navParams || {});
  link.textContent = label;
  return link;
}

// --- Body textures -------------------------------------------------------
//
// A scene this size (a sector holds "only a handful of systems" -- see
// spaceSector.py) is nowhere near enough markers for a fresh canvas
// texture + fresh sphere geometry/material per instance to matter -- no
// shared texture atlas/instancing here because there's no need for one
// at this scale, the same reasoning this file already applied to its old
// sprite textures.

function makeRingTexture(color) {
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

// A neighboring-sector indicator's own small filled dot -- same recipe
// `static/galaxymap3d.js`'s own `makeDotTexture` uses for its "placed"/
// "planned" sector dots, one galaxy-map zoom level up from this
// sector-level view.
function makeDotTexture(fillColor, strokeColor) {
  var size = 64;
  var canvasEl = document.createElement("canvas");
  canvasEl.width = canvasEl.height = size;
  var ctx = canvasEl.getContext("2d");
  var r = size / 2;
  ctx.beginPath();
  ctx.arc(r, r, r - 3, 0, Math.PI * 2);
  ctx.fillStyle = fillColor;
  ctx.fill();
  if (strokeColor) {
    ctx.lineWidth = 3;
    ctx.strokeStyle = strokeColor;
    ctx.stroke();
  }
  var texture = new THREE.CanvasTexture(canvasEl);
  texture.colorSpace = THREE.SRGBColorSpace;
  return texture;
}

function makeSimpleRadialTexture(size, stops) {
  var canvasEl = document.createElement("canvas");
  canvasEl.width = canvasEl.height = size;
  var ctx = canvasEl.getContext("2d");
  var r = size / 2;
  var gradient = ctx.createRadialGradient(r, r, 0, r, r, r);
  stops.forEach(function (stop) {
    gradient.addColorStop(stop[0], stop[1]);
  });
  ctx.fillStyle = gradient;
  ctx.fillRect(0, 0, size, size);
  var texture = new THREE.CanvasTexture(canvasEl);
  texture.colorSpace = THREE.SRGBColorSpace;
  return texture;
}

function makeNebulaTexture(coreColor, edgeColor) {
  return makeSimpleRadialTexture(128, [
    [0, coreColor],
    [0.55, edgeColor],
    [0.78, "rgba(0,0,0,0)"],
  ]);
}

// Reproduces lib/starmap.py's old (now retired) `_ASTEROID_FIELD_BACKGROUND`
// -- a soft tan base disc under three dark "clump" splotches -- as three
// canvas radial-gradient fills instead of four stacked CSS ones. Not
// pixel-identical (CSS's own unsized `radial-gradient(circle at X% Y%, ...)`
// scales each clump to that *element's* own farthest-corner distance from
// its center, not a fixed fraction of the sprite's radius the way this
// does), just the same "mottled rocky scatter" read at a glance.
function makeAsteroidTexture() {
  var size = 128;
  var canvasEl = document.createElement("canvas");
  canvasEl.width = canvasEl.height = size;
  var ctx = canvasEl.getContext("2d");

  function radialDisc(cx, cy, radius, color) {
    var gradient = ctx.createRadialGradient(cx, cy, 0, cx, cy, radius);
    gradient.addColorStop(0, color);
    gradient.addColorStop(0.85, color);
    gradient.addColorStop(1, "rgba(0,0,0,0)");
    ctx.fillStyle = gradient;
    ctx.beginPath();
    ctx.arc(cx, cy, radius, 0, Math.PI * 2);
    ctx.fill();
  }

  var base = ctx.createRadialGradient(size / 2, size / 2, 0, size / 2, size / 2, size / 2);
  base.addColorStop(0, "#b89a6ea0");
  base.addColorStop(0.55, "#b89a6e50");
  base.addColorStop(0.78, "rgba(0,0,0,0)");
  ctx.fillStyle = base;
  ctx.fillRect(0, 0, size, size);

  // Drawn back-to-front relative to the CSS recipe's own stacking order
  // (its first-listed layer is topmost) -- painted last here instead.
  radialDisc(size * 0.68, size * 0.58, size * 0.14, "#00000060");
  radialDisc(size * 0.42, size * 0.78, size * 0.1, "#00000055");
  radialDisc(size * 0.3, size * 0.32, size * 0.12, "#00000070");

  var texture = new THREE.CanvasTexture(canvasEl);
  texture.colorSpace = THREE.SRGBColorSpace;
  return texture;
}

// Every cloud kind besides "nebula" is a fixed recipe (never a function of
// per-instance data beyond which kind/descriptor it is) -- lib/starmap.py
// only ever sends `kind`, these draw what it used to mean by
// `_BLACK_HOLE_ACCRETING_BACKGROUND`/`_BLACK_HOLE_QUIESCENT_BACKGROUND`/
// `_NEUTRON_STAR_BACKGROUND` (also now retired from there).
var CLOUD_KIND_RECIPES = {
  asteroidField: makeAsteroidTexture,
  blackHoleAccreting: function () {
    return makeSimpleRadialTexture(96, [
      [0, "#000000f5"], [0.34, "#000000f5"], [0.55, "#ff9d4dc0"], [0.72, "#ff9d4d30"], [0.86, "rgba(0,0,0,0)"],
    ]);
  },
  blackHoleQuiescent: function () {
    return makeSimpleRadialTexture(96, [
      [0, "#000000f5"], [0.55, "#000000f5"], [0.78, "#4b2f6660"], [0.9, "rgba(0,0,0,0)"],
    ]);
  },
  neutronStar: function () {
    return makeSimpleRadialTexture(96, [
      [0, "#ffffff"], [0.35, "#cfe8ffe0"], [0.6, "#8fc7ff80"], [0.82, "rgba(0,0,0,0)"],
    ]);
  },
};

function textureForCloud(cloud) {
  if (cloud.kind === "nebula") {
    return makeNebulaTexture(cloud.coreColor, cloud.edgeColor);
  }
  var recipe = CLOUD_KIND_RECIPES[cloud.kind];
  return recipe ? recipe() : makeNebulaTexture("#c9a8e090", "#c9a8e030");
}

// Extra fresnel glow shell layered over each body's own core texture
// above (see ./bodyRendering.js's own makeGlowMaterial) -- a true 3D
// "brightest at the silhouette, from any angle" glow the flat radial
// textures above can't give on their own once mapped onto a sphere
// (their own gradient varies by latitude there, not by view angle).
// power/strength/scale tuned per kind so each reads as "this kind of
// object glows this much" -- a neutron star's tight, intense beam vs. a
// quiescent black hole's dim, contained one.
var CLOUD_GLOW_RECIPES = {
  asteroidField: { color: "#b89a6e", power: 3.0, strength: 0.5, scale: 1.12 },
  blackHoleAccreting: { color: "#ff9d4d", power: 1.8, strength: 1.6, scale: 1.4 },
  blackHoleQuiescent: { color: "#4b2f66", power: 2.5, strength: 0.7, scale: 1.2 },
  neutronStar: { color: "#8fc7ff", power: 1.2, strength: 2.4, scale: 1.45 },
};

function glowRecipeForCloud(cloud) {
  if (cloud.kind === "nebula") {
    // cloud.coreColor is an 8-digit RGBA hex (alpha baked in, see
    // lib/starmap.py's own data["coreColor"]) -- stripped to a plain
    // 6-digit color for the glow shader's own uniform, which controls
    // opacity itself via glowStrength rather than a texture alpha
    // channel.
    var core = (cloud.coreColor || "#c9a8e090").slice(0, 7);
    return { color: core, power: 2.2, strength: 1.3, scale: 1.3 };
  }
  return CLOUD_GLOW_RECIPES[cloud.kind] || CLOUD_GLOW_RECIPES.asteroidField;
}

// A real 3D body: a textured core sphere plus a fresnel glow shell
// (./bodyRendering.js), replacing the flat always-facing-camera Sprite
// this file used to draw for every star/cloud. `coreTexture` is
// whatever this kind's own recipe already built above (a star's own
// granulation via makeStarSurfaceTexture, or one of the existing
// nebula/asteroid-field/compact-remnant textures reused unchanged);
// `coreColorHex` tints it (0xffffff/no-op when the texture already has
// its own baked-in color, as every one of these does). Returns
// `{core, glow}` -- the caller adds `core` to the raycastable
// interactiveGroup and `glow` directly to the scene (see this file's own
// module docstring on why the glow shell is never a raycast target of
// its own).
function makeBodySpheres(radius, coreTexture, coreColorHex, glow) {
  var coreMaterial = new THREE.MeshBasicMaterial({
    map: coreTexture || null, color: coreColorHex, transparent: !!coreTexture,
  });
  var core = new THREE.Mesh(new THREE.SphereGeometry(1, 32, 24), coreMaterial);
  core.scale.setScalar(radius);

  var glowMaterial = makeGlowMaterial(THREE, glow.color, glow.power, glow.strength);
  var glowMesh = new THREE.Mesh(new THREE.SphereGeometry(1, 32, 24), glowMaterial);
  glowMesh.scale.setScalar(radius * glow.scale);

  return { core: core, glow: glowMesh };
}

function makeTextSprite(text, color) {
  var measuring = document.createElement("canvas").getContext("2d");
  var font = "600 28px system-ui, -apple-system, Segoe UI, Roboto, sans-serif";
  measuring.font = font;
  var textWidth = measuring.measureText(text).width;

  var paddingX = 14;
  var canvasEl = document.createElement("canvas");
  canvasEl.width = Math.ceil(textWidth) + paddingX * 2;
  canvasEl.height = 40;
  var ctx = canvasEl.getContext("2d");
  ctx.font = font;
  ctx.fillStyle = color;
  ctx.textBaseline = "middle";
  ctx.textAlign = "left";
  ctx.fillText(text, paddingX, canvasEl.height / 2);

  var texture = new THREE.CanvasTexture(canvasEl);
  texture.colorSpace = THREE.SRGBColorSpace;
  var sprite = new THREE.Sprite(new THREE.SpriteMaterial({ map: texture, transparent: true, depthWrite: false }));
  var worldHeight = 22;
  sprite.scale.set(worldHeight * (canvasEl.width / canvasEl.height), worldHeight, 1);
  return sprite;
}

// --- Scene setup ---------------------------------------------------------

function initStarmap(canvasEl, data) {
  var viewport = canvasEl.closest(".starmap-viewport");

  var renderer = new THREE.WebGLRenderer({ canvas: canvasEl, antialias: true, alpha: true, logarithmicDepthBuffer: true });
  renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
  renderer.setClearColor(0x000000, 0);
  if (THREE.SRGBColorSpace) {
    renderer.outputColorSpace = THREE.SRGBColorSpace;
  }

  var scene = new THREE.Scene();

  var FOV_DEG = 45;
  var camera = new THREE.PerspectiveCamera(FOV_DEG, 1, 1, 5e6);

  // The world-unit distance at which a sphere of radius `sceneHalfPx`
  // (lib/starmap.py's own scale reference -- every position/radius in
  // `data` is expressed in these units) exactly fills the frame
  // vertically -- this is what "zoom = 1" means for a real camera, the
  // direct replacement for the old CSS version's `scale(1)`.
  var sceneHalfPx = data.sceneHalfPx || 160;
  var referenceDistance = sceneHalfPx / Math.tan(THREE.MathUtils.degToRad(FOV_DEG / 2));

  var MIN_ZOOM = 0.2;
  var MAX_ZOOM = 2.5;
  var ZOOM_STEP = 0.15;
  var WHEEL_ZOOM_STEP = 0.08;
  var KEY_ROTATE_STEP = THREE.MathUtils.degToRad(6);
  var ROTATE_SENSITIVITY = THREE.MathUtils.degToRad(0.4); // radians per pixel of drag
  var DRAG_CLICK_THRESHOLD_PX = 4;
  var MIN_POLAR = THREE.MathUtils.degToRad(2);
  var MAX_POLAR = THREE.MathUtils.degToRad(178);

  var defaultZoom = data.defaultZoom > 0 && data.defaultZoom <= 1 ? data.defaultZoom : 1;
  var defaultAzimuth = THREE.MathUtils.degToRad(-32);
  var defaultPolar = THREE.MathUtils.degToRad(90 - 18);

  var spherical = new THREE.Spherical(referenceDistance / defaultZoom, defaultPolar, defaultAzimuth);

  function applyCamera() {
    camera.position.setFromSpherical(spherical);
    camera.lookAt(0, 0, 0);
  }
  applyCamera();

  function currentZoom() {
    return referenceDistance / spherical.radius;
  }

  function setZoom(zoom) {
    spherical.radius = referenceDistance / Math.max(MIN_ZOOM, Math.min(MAX_ZOOM, zoom));
    applyCamera();
    updateScaleBar();
  }

  function resetView() {
    spherical.set(referenceDistance / defaultZoom, defaultPolar, defaultAzimuth);
    applyCamera();
    updateScaleBar();
  }

  var accentColor = cssVar("--accent", "#4f5fe8");

  if (data.compass) {
    var tip = data.compass.tip;
    var arrowGeometry = new THREE.BufferGeometry().setFromPoints([
      new THREE.Vector3(0, 0, 0),
      new THREE.Vector3(tip[0], tip[1], tip[2]),
    ]);
    scene.add(new THREE.Line(arrowGeometry, new THREE.LineBasicMaterial({ color: new THREE.Color(accentColor) })));

    // Plain "N" at the arrow's own tip (the standard compass-rose
    // convention) -- no extra arrow glyph appended to the text itself,
    // since the line already drawn above IS the arrow.
    var label = makeTextSprite(data.compass.label, accentColor);
    label.position.set(tip[0], tip[1], tip[2]);
    scene.add(label);
  }

  // Every clickable/focusable marker -- raycasting and the accessible
  // fallback button list both only ever need to search this, not the
  // compass (which carries no `data-*`-equivalent info of its own).
  var interactiveGroup = new THREE.Group();
  scene.add(interactiveGroup);
  var entryByObject = new Map();

  // A star's own glow: bright/wide, like static/systemmap.js's identical
  // STAR_GLOW_POWER/_STRENGTH/_SCALE tuning for its own star spheres --
  // shared visual language for "this is a star" across both maps.
  var STAR_GLOW = { power: 1.5, strength: 2.2, scale: 1.4 };

  (data.stars || []).forEach(function (star) {
    var bodies = makeBodySpheres(
      star.r, makeStarSurfaceTexture(THREE, star.fill), 0xffffff,
      { color: star.fill, power: STAR_GLOW.power, strength: STAR_GLOW.strength, scale: STAR_GLOW.scale },
    );
    bodies.core.position.set(star.x, star.y, star.z);
    bodies.glow.position.set(star.x, star.y, star.z);
    interactiveGroup.add(bodies.core);
    scene.add(bodies.glow);
    entryByObject.set(bodies.core, star);
  });

  (data.clouds || []).forEach(function (cloud) {
    var bodies = makeBodySpheres(cloud.r, textureForCloud(cloud), 0xffffff, glowRecipeForCloud(cloud));
    bodies.core.position.set(cloud.x, cloud.y, cloud.z);
    bodies.glow.position.set(cloud.x, cloud.y, cloud.z);
    interactiveGroup.add(bodies.core);
    scene.add(bodies.glow);
    entryByObject.set(bodies.core, cloud);
  });

  // Neighboring-sector indicators: a small flat dot at the scene's own
  // edge, in the real direction of that neighbor (see lib/starmap.py's
  // own `_neighbor_indicator_data`) -- a plain billboard `THREE.Sprite`
  // (not a real 3D body like a star/cloud above) suits these fine, the
  // same reasoning the highlight ring and compass label below are
  // sprites too: flat 2D markers with no "seen from any angle" concern.
  // Bright accent color for an already-generated neighbor (clickable,
  // navigable, like every other marker in this scene); a muted blue for
  // one that isn't yet -- the same "planned" color
  // `static/galaxymap3d.js`'s own dots use for the identical not-yet-
  // generated concept one galaxy-map zoom level up.
  var NEIGHBOR_EXISTS_COLOR = accentColor;
  var NEIGHBOR_MISSING_FILL = "#7fa8d9";
  var NEIGHBOR_MISSING_STROKE = "#3f5f80";
  var neighborExistsTexture = null;
  var neighborMissingTexture = null;

  (data.neighbors || []).forEach(function (neighbor) {
    if (!neighborExistsTexture) {
      neighborExistsTexture = makeDotTexture(NEIGHBOR_EXISTS_COLOR, null);
    }
    if (!neighborMissingTexture) {
      neighborMissingTexture = makeDotTexture(NEIGHBOR_MISSING_FILL, NEIGHBOR_MISSING_STROKE);
    }
    var marker = new THREE.Sprite(
      new THREE.SpriteMaterial({
        map: neighbor.exists ? neighborExistsTexture : neighborMissingTexture,
        transparent: true, depthWrite: false,
      })
    );
    marker.position.set(neighbor.x, neighbor.y, neighbor.z);
    marker.scale.set(neighbor.r * 2, neighbor.r * 2, 1);
    scene.add(marker);
    interactiveGroup.add(marker);
    entryByObject.set(marker, neighbor);

    var labelText = neighbor.exists ? neighbor.name : neighbor.designation;
    if (labelText) {
      var label = makeTextSprite(labelText, neighbor.exists ? NEIGHBOR_EXISTS_COLOR : NEIGHBOR_MISSING_FILL);
      label.position.set(neighbor.x, neighbor.y - neighbor.r * 2.4, neighbor.z);
      scene.add(label);
    }
  });

  var highlightSprite = new THREE.Sprite(
    new THREE.SpriteMaterial({ map: makeRingTexture(accentColor), transparent: true, depthWrite: false })
  );
  highlightSprite.visible = false;
  scene.add(highlightSprite);

  function highlightEntry(entry) {
    highlightSprite.position.set(entry.x, entry.y, entry.z);
    var r = entry.r || 8;
    highlightSprite.scale.set(r * 2.6, r * 2.6, 1);
    highlightSprite.visible = true;
  }

  function selectEntry(entry) {
    if (!entry) {
      return;
    }
    showObjectInfo(entry);
    highlightEntry(entry);
  }

  // --- Accessible fallback list -----------------------------------------
  //
  // A canvas has no focusable children of its own the way the old CSS
  // version's real per-star `<div role="button">`s were, so this is what
  // keeps every star/cloud/neighbor reachable by keyboard/screen reader
  // without needing 3D hit-testing or focus management inside the canvas
  // itself -- a visually hidden button per entry, in the same list order
  // the scene data arrived in.
  function entryLabel(entry) {
    if (entry.isNeighbor) {
      return entry.exists ? entry.name || "Unnamed sector" : "Not yet generated (" + entry.designation + ")";
    }
    return entry.name || "Unknown";
  }

  if (viewport) {
    var list = document.createElement("ul");
    list.className = "starmap-sr-list sr-only";
    (data.stars || []).concat(data.clouds || []).concat(data.neighbors || []).forEach(function (entry) {
      var item = document.createElement("li");
      var button = document.createElement("button");
      button.type = "button";
      button.textContent = entryLabel(entry);
      button.addEventListener("click", function () {
        selectEntry(entry);
      });
      item.appendChild(button);
      list.appendChild(item);
    });
    viewport.appendChild(list);
  }

  // --- Pointer/keyboard interaction --------------------------------------

  var dragging = false;
  var dragDistance = 0;
  var lastClientX = 0;
  var lastClientY = 0;
  var suppressNextClick = false;

  canvasEl.addEventListener("pointerdown", function (event) {
    dragging = true;
    dragDistance = 0;
    lastClientX = event.clientX;
    lastClientY = event.clientY;
    try {
      canvasEl.setPointerCapture(event.pointerId);
    } catch (err) {
      // Pointer capture isn't essential -- dragging still works via
      // ordinary pointermove bubbling if the browser refuses it.
    }
  });

  canvasEl.addEventListener("pointermove", function (event) {
    if (!dragging) {
      return;
    }
    var deltaX = event.clientX - lastClientX;
    var deltaY = event.clientY - lastClientY;
    dragDistance += Math.abs(deltaX) + Math.abs(deltaY);
    lastClientX = event.clientX;
    lastClientY = event.clientY;

    spherical.theta -= deltaX * ROTATE_SENSITIVITY;
    spherical.phi = Math.max(MIN_POLAR, Math.min(MAX_POLAR, spherical.phi - deltaY * ROTATE_SENSITIVITY));
    applyCamera();
    updateScaleBar();
  });

  function endDrag(event) {
    dragging = false;
    if (dragDistance > DRAG_CLICK_THRESHOLD_PX) {
      suppressNextClick = true;
      // Safety net: a pointerup isn't always followed by a click (e.g.
      // pointercancel) -- don't leave this suppressing some unrelated
      // later click if one never arrives to consume and clear it.
      setTimeout(function () {
        suppressNextClick = false;
      }, 0);
    }
    dragDistance = 0;
    try {
      canvasEl.releasePointerCapture(event.pointerId);
    } catch (err) {
      // Already released/invalid -- nothing to clean up.
    }
  }
  canvasEl.addEventListener("pointerup", endDrag);
  canvasEl.addEventListener("pointercancel", endDrag);

  canvasEl.addEventListener(
    "wheel",
    function (event) {
      event.preventDefault();
      setZoom(currentZoom() + (event.deltaY < 0 ? WHEEL_ZOOM_STEP : -WHEEL_ZOOM_STEP));
    },
    { passive: false }
  );

  canvasEl.addEventListener("keydown", function (event) {
    var key = event.key;
    if (key !== "ArrowLeft" && key !== "ArrowRight" && key !== "ArrowUp" && key !== "ArrowDown") {
      return;
    }
    event.preventDefault();
    if (key === "ArrowLeft") spherical.theta += KEY_ROTATE_STEP;
    if (key === "ArrowRight") spherical.theta -= KEY_ROTATE_STEP;
    if (key === "ArrowUp") spherical.phi = Math.max(MIN_POLAR, spherical.phi - KEY_ROTATE_STEP);
    if (key === "ArrowDown") spherical.phi = Math.min(MAX_POLAR, spherical.phi + KEY_ROTATE_STEP);
    applyCamera();
    updateScaleBar();
  });

  var raycaster = new THREE.Raycaster();

  function entryAtClientPoint(clientX, clientY) {
    var rect = canvasEl.getBoundingClientRect();
    if (rect.width === 0 || rect.height === 0) {
      return null;
    }
    var ndc = new THREE.Vector2(
      ((clientX - rect.left) / rect.width) * 2 - 1,
      -((clientY - rect.top) / rect.height) * 2 + 1
    );
    raycaster.setFromCamera(ndc, camera);
    var hits = raycaster.intersectObjects(interactiveGroup.children, false);
    return hits.length ? entryByObject.get(hits[0].object) || null : null;
  }

  canvasEl.addEventListener("click", function (event) {
    if (suppressNextClick) {
      suppressNextClick = false;
      return;
    }
    selectEntry(entryAtClientPoint(event.clientX, event.clientY));
  });

  var controlsEl = document.getElementById("starmap-controls");
  if (controlsEl) {
    controlsEl.querySelectorAll("[data-action]").forEach(function (button) {
      button.addEventListener("click", function () {
        var action = button.dataset.action;
        if (action === "zoom-in") setZoom(currentZoom() + ZOOM_STEP);
        else if (action === "zoom-out") setZoom(currentZoom() - ZOOM_STEP);
        else if (action === "reset") resetView();
      });
    });
  }

  // --- Scale bar -----------------------------------------------------------

  var scaleEl = document.getElementById("starmap-scale");
  var scaleBarEl = document.getElementById("starmap-scale-bar");
  var scaleLabelEl = document.getElementById("starmap-scale-label");
  var lyPerWorldUnit = data.lyPerPxAtZoom1 || 0;
  var SCALE_BAR_TARGET_PX = 70;

  // Snaps an arbitrary positive value to the nearest "nice" 1/2/5 * 10^n
  // -- the standard map-scale-bar convention, so the label reads "5 ly"
  // or "20 ly" rather than an ugly "6.283 ly".
  function niceScaleValue(raw) {
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

  function formatLy(value) {
    if (value >= 100) return Math.round(value) + " ly";
    if (value >= 1) return Math.round(value * 10) / 10 + " ly";
    return Math.round(value * 1000) / 1000 + " ly";
  }

  // Unlike the old CSS version (a fixed 320px scene scaled by a flat CSS
  // `zoom` factor, so ly-per-pixel was that one ratio divided by `zoom`),
  // a real perspective camera's screen-pixels-per-world-unit depends on
  // both camera distance *and* the canvas's own live rendered size (this
  // panel is a responsive `min(100%, 22rem)` box, not a fixed 320px one)
  // -- so this recomputes it from first principles every time instead.
  function worldUnitsPerScreenPixel() {
    var fovRad = THREE.MathUtils.degToRad(camera.fov);
    var heightPx = canvasEl.clientHeight || 1;
    return (2 * spherical.radius * Math.tan(fovRad / 2)) / heightPx;
  }

  function updateScaleBar() {
    if (!scaleEl || !scaleBarEl || !scaleLabelEl || !lyPerWorldUnit) {
      return;
    }
    var lyPerScreenPx = worldUnitsPerScreenPixel() * lyPerWorldUnit;
    var niceLy = niceScaleValue(SCALE_BAR_TARGET_PX * lyPerScreenPx);
    if (!niceLy) {
      return;
    }
    scaleBarEl.style.width = (niceLy / lyPerScreenPx).toFixed(1) + "px";
    scaleLabelEl.textContent = formatLy(niceLy);
  }

  // --- Resize/render loop ----------------------------------------------

  function resize() {
    var width = canvasEl.clientWidth || 1;
    var height = canvasEl.clientHeight || 1;
    renderer.setSize(width, height, false);
    camera.aspect = width / height;
    camera.updateProjectionMatrix();
    updateScaleBar();
  }

  if (typeof ResizeObserver !== "undefined" && viewport) {
    new ResizeObserver(resize).observe(viewport);
  }
  resize();
  window.addEventListener("resize", resize);

  (function animate() {
    requestAnimationFrame(animate);
    renderer.render(scene, camera);
  })();
}

if (canvas && sceneData) {
  initStarmap(canvas, sceneData);
}
