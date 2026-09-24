// html/static/galaxymap3d.js
//
// Renders the interactive 3D Galaxy Map (lib/galaxymap3d.py) as a real
// WebGL scene (three.js, vendored at ./vendor/three.module.min.js -- see
// sectormap.js's own docstring for why it's vendored rather than loaded
// from a CDN). Unlike sectormap.js's scene (every star/cloud baked into
// one JSON block at page load, camera always orbiting a fixed origin),
// this camera's own orbit TARGET moves freely through the galaxy -- most
// of what's drawn is fetched live from galaxy_tiles.py as the camera
// moves, one fixed cube of space ("tile") at a time from galaxy_tiles.py,
// never baked into the page beyond the very first frame
// (#galaxymap3d-data's own "initial" payload) -- see the "Cube tiles"
// section below.
//
// Coordinate convention: every position here is the design doc's own
// galaxy-frame parsecs (docs/design/galaxy-coordinate-system.md), passed
// straight through as three.js world units with NO axis flip/remap --
// x/y span the galactic plane, +z is galactic north. `camera.up` is set
// to (0, 0, 1) once, at scene setup, specifically so that convention
// reads as "up" on screen without needing any per-position sign flip the
// way this map's own former flat SVG projection (a legacy "+y is down on
// screen" convention) needed one.
//
// Three content tiers (each tile's placed/planned lists plus a separate
// density cloud -- see stellarObjects.galaxyViewport's module docstring):
//   - placed:  real, already-generated sectors -- bright sprites, sized
//              by system_count, click selects + navigates.
//   - planned: real, not-yet-generated qualifying addresses -- small dim
//              sprites, click selects (shows the copyable designation/
//              CLI snippet) but never navigates.
//   - density: an illustrative cloud of soft, translucent, additively-
//              blended spheres (THREE.InstancedMesh, one shared low-poly
//              SphereGeometry -- not individual Sprite/Mesh objects,
//              which wouldn't scale to a thousand-plus instances as
//              cheaply) -- not interactive, no identity to track. Sized
//              in real world-space parsecs (unlike placed/planned's own
//              constant-screen-pixel markers below), scaled per-instance
//              by that point's own relative_density AND by the camera's
//              current orbit radius (see DENSITY_RADIUS_FRACTION) so the
//              cloud keeps reading as roughly the same relative size
//              across zoom levels. Additive blending lets overlapping
//              spheres brighten where they overlap rather than simply
//              occluding each other -- the cheap way many soft
//              transparent blobs merge into continuous-looking shading
//              along a spiral arm instead of reading as a sparse
//              scatter-plot of discrete dots.
//
// Click-to-zoom is LOGARITHMIC, not a flat factor: clickZoomFactor()
// below interpolates between lib/galaxymap3d.py's own
// clickZoomFactorMin/Max by the camera's CURRENT distance in log space,
// recomputed fresh on every click -- big multiplicative jumps while
// zoomed out over the whole galaxy, small fine ones once close to a
// single sector, so a handful of clicks crosses the galaxy without ever
// needing so fine a step that reaching the core takes dozens of clicks,
// and without a step so coarse up close that it blows past the one
// sector being approached.

import * as THREE from "./vendor/three.module.min.js";

var canvas = document.getElementById("galaxymap3d-canvas");
var dataEl = document.getElementById("galaxymap3d-data");

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

// A real, focusable <a> carrying data-nav-target/data-nav-params instead
// of an href query string -- static/navform.js's document-level click
// handler (loaded on every page) is what actually follows it. Same
// convention sectormap.js's own navLink uses.
function navLink(navTarget, navParams, label) {
  var link = document.createElement("a");
  link.href = "#";
  link.className = "btn";
  link.dataset.navTarget = navTarget;
  link.dataset.navParams = JSON.stringify(navParams || {});
  link.textContent = label;
  return link;
}

function formatAddress(shellIndex, slotIndex) {
  return "shell " + shellIndex + " slot " + slotIndex;
}

function cliSnippet(shellIndex, slotIndex) {
  return "generate.py galaxy --shell " + shellIndex + " --slot " + slotIndex;
}

// --- Info panel ----------------------------------------------------------

function showPlacedInfo(entry) {
  var panel = document.getElementById("galaxymap3d-info");
  if (!panel) {
    return;
  }
  panel.textContent = "";

  var heading = document.createElement("h3");
  heading.textContent = entry.name || "Unnamed sector";
  panel.appendChild(heading);

  var relativeDensity = placedRelativeDensity(entry, sceneData.referenceDensityPerLy3);

  var dl = document.createElement("dl");
  addField(dl, "Systems", entry.system_count != null ? entry.system_count : 0);
  addField(dl, "Density", relativeDensity != null ? relativeDensity.toFixed(2) + "× local average" : null);
  addField(dl, "Distance from core", entry.galactic_radius_pc != null ? Math.round(entry.galactic_radius_pc) + " pc" : null);
  addField(dl, "Address", entry.shell_index != null ? formatAddress(entry.shell_index, entry.shell_slot_index) : null);
  addField(dl, "Designation", entry.designation);
  panel.appendChild(dl);
  panel.appendChild(navLink("sector.py", { db: sceneData.db, id: entry.id }, "View sector →"));
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

function showPlannedInfo(entry) {
  var panel = document.getElementById("galaxymap3d-info");
  if (!panel) {
    return;
  }
  panel.textContent = "";

  var heading = document.createElement("h3");
  heading.textContent = "Not yet generated";
  panel.appendChild(heading);

  var dl = document.createElement("dl");
  addField(dl, "Designation", entry.designation);
  addField(dl, "Address", formatAddress(entry.shell_index, entry.shell_slot_index));
  if (entry.predicted_star_count != null) {
    addField(dl, "Predicted systems", Math.max(1, Math.round(entry.predicted_star_count)));
  }
  panel.appendChild(dl);

  var code = document.createElement("code");
  code.className = "galaxymap3d-cli-snippet";
  code.textContent = cliSnippet(entry.shell_index, entry.shell_slot_index);
  panel.appendChild(code);
  panel.appendChild(makeCopyButton(cliSnippet(entry.shell_index, entry.shell_slot_index)));
}

// --- Sprite textures -------------------------------------------------------

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

// --- Per-tier visual formulas -----------------------------------------
//
// Unlike lib/starmap.py (size/color baked server-side, once), these run
// client-side for EVERY fetch (initial payload and every live re-fetch
// alike) -- see galaxymap3d.py's own module docstring for why: almost
// everything drawn here arrives through a live fetch that never passes
// through that Python module again, so a server-computed style would
// only ever apply to the first frame.

// Every marker below is sized in constant SCREEN pixels, not world
// parsecs -- recomputed every frame from each object's own live distance
// to the camera (updateMarkerScales, in the render loop below), the
// standard "billboard with constant screen size" technique. A world-unit
// sprite radius (three.js's own Sprite default) would perspective-shrink
// with distance like any other object: correct for something meant to
// represent real physical size, but wrong for a point-of-interest marker
// that needs to stay visible/clickable from a full-galaxy overview
// thousands of parsecs out, the same way a system marker on a game's own
// galaxy map stays a legible dot regardless of camera distance.
var PLACED_MIN_PX = 3.0;
var PLACED_MAX_PX = 7.0;
// Near-white/gray, not a fixed hue -- placedColor() below tints each
// marker's own sprite material to its own real-density color via plain
// multiplication (THREE.SpriteMaterial's own `color` * texture), which
// only stays clean (scales brightness/saturation) against a white/gray
// base; a colored base texture would shift hue unpredictably instead.
var PLACED_CORE_BASE_FILL = "#ffffff";
var PLACED_CORE_BASE_STROKE = "#c4c4c4";
var PLACED_HALO_BASE_FILL = "#ffffff";

function placedScreenRadiusPx(systemCount) {
  var count = systemCount || 0;
  return Math.max(PLACED_MIN_PX, Math.min(PLACED_MAX_PX, PLACED_MIN_PX + 0.6 * Math.sqrt(count)));
}

// A placed sector's own REAL stellar density (system_count / edge_ly^3),
// relative to physical_constants.LOCAL_STELLAR_DENSITY_LY3 (the real
// local-neighborhood average this whole generator already calibrates
// against -- see lib/galaxymap3d.py's own referenceDensityPerLy3
// comment) -- 1.0 means exactly average, >1 denser, <1 sparser. `null`
// when edge_ly isn't available (a sector placed before per-sector edge
// tracking existed) rather than a false 0, so placedDensityColor below
// can fall back to a neutral mid-tone instead of reading as "empty".
var PLACED_LOW_DENSITY_COLOR = "#4a3f2e";
var PLACED_HIGH_DENSITY_COLOR = "#fff6df";

function placedRelativeDensity(entry, referenceDensityPerLy3) {
  if (!entry.edge_ly || !referenceDensityPerLy3) {
    return null;
  }
  var densityPerLy3 = (entry.system_count || 0) / Math.pow(entry.edge_ly, 3);
  return densityPerLy3 / referenceDensityPerLy3;
}

// Same "log2(x+1)/3" shape densityIntensity (below, for the illustrative
// cloud) uses, for a consistent dim-to-bright response curve -- fed real
// per-sector density here instead of the server's own illustrative
// model, and mapped through a warm bronze-to-gold range (rather than the
// cloud's cooler dim-to-accent one) so a real, already-generated
// sector's own marker stays visually distinct from illustrative shading
// at a glance, exactly like it already was before density coloring.
function placedDensityColor(entry, referenceDensityPerLy3) {
  var relative = placedRelativeDensity(entry, referenceDensityPerLy3);
  var t = relative == null ? 0.5 : Math.max(0, Math.min(1, Math.log2(relative + 1) / 3));
  var low = new THREE.Color(PLACED_LOW_DENSITY_COLOR);
  var high = new THREE.Color(PLACED_HIGH_DENSITY_COLOR);
  return low.clone().lerp(high, t);
}

var PLANNED_PX = 3.0;
var PLANNED_FILL = "#7fa8d9";
var PLANNED_STROKE = "#3f5f80";

var placedCoreTexture = null;
var placedHaloTexture = null;
var plannedTexture = null;
var highlightTexture = null;

function ensureTextures(accentColor) {
  if (!placedCoreTexture) {
    placedCoreTexture = makeDotTexture(PLACED_CORE_BASE_FILL, PLACED_CORE_BASE_STROKE);
    placedHaloTexture = makeDotTexture(PLACED_HALO_BASE_FILL, null);
    plannedTexture = makeDotTexture(PLANNED_FILL, PLANNED_STROKE);
    highlightTexture = makeRingTexture(accentColor);
  }
}

// --- Scene setup -----------------------------------------------------------

function initGalaxyMap3d(canvasEl, data) {
  var viewport = canvasEl.closest(".starmap-viewport");
  var accentColor = cssVar("--accent", "#4f5fe8");
  ensureTextures(accentColor);

  var renderer = new THREE.WebGLRenderer({ canvas: canvasEl, antialias: true, alpha: true, logarithmicDepthBuffer: true });
  renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
  renderer.setClearColor(0x000000, 0);
  if (THREE.SRGBColorSpace) {
    renderer.outputColorSpace = THREE.SRGBColorSpace;
  }

  var scene = new THREE.Scene();

  var FOV_DEG = 50;
  var farPlane = Math.max(data.maxViewRadiusPc * 6, 1000);
  var camera = new THREE.PerspectiveCamera(FOV_DEG, 1, Math.max(data.minViewRadiusPc / 50, 0.001), farPlane);
  camera.up.set(0, 0, 1);

  var MIN_RADIUS = data.minViewRadiusPc;
  var MAX_RADIUS = data.maxViewRadiusPc;
  var CLICK_FACTOR_MIN = data.clickZoomFactorMin || 1.15;
  var CLICK_FACTOR_MAX = data.clickZoomFactorMax || 4.0;
  var WHEEL_ZOOM_RATIO = 1.1;
  var KEY_ROTATE_STEP = THREE.MathUtils.degToRad(6);
  var ROTATE_SENSITIVITY = THREE.MathUtils.degToRad(0.4);
  var DRAG_CLICK_THRESHOLD_PX = 4;
  var MIN_POLAR = THREE.MathUtils.degToRad(2);
  var MAX_POLAR = THREE.MathUtils.degToRad(178);

  var target = new THREE.Vector3(data.initialCenter[0], data.initialCenter[1], data.initialCenter[2]);
  var initialTarget = target.clone();
  var initialRadius = data.initialRadiusPc;

  // theta: azimuth from +x in the xy-plane; phi: polar angle from +z --
  // the same (r, theta, phi) convention docs/design/galaxy-coordinate-
  // system.md and stellarObjects.galaxyGeometry.sector_position_pc use,
  // deliberately not THREE.Spherical (which assumes a +y-up world).
  var orbit = { radius: initialRadius, theta: THREE.MathUtils.degToRad(-32), phi: THREE.MathUtils.degToRad(60) };
  var initialOrbit = { radius: orbit.radius, theta: orbit.theta, phi: orbit.phi };

  function offsetFromOrbit(o) {
    var sinPhi = Math.sin(o.phi);
    return new THREE.Vector3(
      o.radius * sinPhi * Math.cos(o.theta),
      o.radius * sinPhi * Math.sin(o.theta),
      o.radius * Math.cos(o.phi)
    );
  }

  function applyCamera() {
    var offset = offsetFromOrbit(orbit);
    camera.position.set(target.x + offset.x, target.y + offset.y, target.z + offset.z);
    camera.lookAt(target);
  }
  applyCamera();

  // --- Logarithmic click-zoom factor ------------------------------------
  //
  // See this file's own module docstring. Computed fresh from whatever
  // radius is passed in (always the CURRENT orbit.radius, before it
  // changes) -- never cached, since the whole point is that it shrinks
  // as the camera gets closer.
  function clickZoomFactor(currentRadius) {
    if (MAX_RADIUS <= MIN_RADIUS) {
      return CLICK_FACTOR_MIN;
    }
    var clamped = Math.max(MIN_RADIUS, Math.min(MAX_RADIUS, currentRadius));
    var t = (Math.log(clamped) - Math.log(MIN_RADIUS)) / (Math.log(MAX_RADIUS) - Math.log(MIN_RADIUS));
    return CLICK_FACTOR_MIN + t * (CLICK_FACTOR_MAX - CLICK_FACTOR_MIN);
  }

  function setRadius(radius) {
    orbit.radius = Math.max(MIN_RADIUS, Math.min(MAX_RADIUS, radius));
    applyCamera();
    updateScaleBar();
  }

  function resetView() {
    target.copy(initialTarget);
    orbit.radius = initialOrbit.radius;
    orbit.theta = initialOrbit.theta;
    orbit.phi = initialOrbit.phi;
    applyCamera();
    updateScaleBar();
    scheduleFetch(true);
  }

  // --- Content groups ------------------------------------------------------

  var interactiveGroup = new THREE.Group();
  scene.add(interactiveGroup);
  var entryByObject = new Map();
  var placedSpritesByKey = new Map();
  var plannedSpritesByKey = new Map();

  // See this file's own module docstring for why this is an InstancedMesh
  // of real, additively-blended spheres rather than the flat THREE.Points
  // scatter this used to be. DENSITY_INSTANCE_CAPACITY is a safety margin
  // above stellarObjects.galaxyViewport.DENSITY_TILE_SAMPLE_COUNT (1600) --
  // InstancedMesh.count (set per-update in applyDensity) can render fewer
  // instances than this capacity with no reallocation, but never more, so
  // this stays a headroom margin rather than an exact mirror of that
  // server-side constant.
  var DENSITY_INSTANCE_CAPACITY = 2000;
  // Each sphere's world-space radius is this fraction of the camera's
  // CURRENT orbit radius (recomputed in applyDensity, which reruns on
  // every live re-fetch as the camera moves) -- not a fixed parsec value,
  // since a fixed size would either vanish at the full-galaxy starting
  // view or dwarf the scene once zoomed in close. DENSITY_MIN/MAX_SCALE
  // then varies that base size per-instance by the point's own
  // relative_density (denser regions read as visibly bigger/brighter
  // blobs, not just differently colored ones).
  var DENSITY_RADIUS_FRACTION = 0.05;
  var DENSITY_MIN_SCALE = 0.6;
  var DENSITY_MAX_SCALE = 2.2;

  var densityGeometry = new THREE.SphereGeometry(1, 12, 10);
  var densityMaterial = new THREE.MeshBasicMaterial({
    vertexColors: true, transparent: true, opacity: 0.4,
    depthWrite: false, blending: THREE.AdditiveBlending,
  });
  var densityMesh = new THREE.InstancedMesh(densityGeometry, densityMaterial, DENSITY_INSTANCE_CAPACITY);
  densityMesh.count = 0;
  densityMesh.frustumCulled = false;
  scene.add(densityMesh);

  var highlightSprite = new THREE.Sprite(
    new THREE.SpriteMaterial({ map: highlightTexture, transparent: true, depthWrite: false })
  );
  highlightSprite.visible = false;
  highlightSprite.userData.screenRadiusPx = PLACED_MAX_PX;
  scene.add(highlightSprite);

  function highlightPosition(x, y, z, screenRadiusPx) {
    highlightSprite.position.set(x, y, z);
    highlightSprite.userData.screenRadiusPx = (screenRadiusPx || PLACED_MIN_PX) * 1.3;
    highlightSprite.visible = true;
  }

  function makePlacedSprite(entry) {
    var group = new THREE.Group();
    var color = placedDensityColor(entry, data.referenceDensityPerLy3);
    var halo = new THREE.Sprite(new THREE.SpriteMaterial({
      map: placedHaloTexture, color: color, transparent: true, depthWrite: false, opacity: 0.45,
    }));
    var core = new THREE.Sprite(new THREE.SpriteMaterial({
      map: placedCoreTexture, color: color, transparent: true, depthWrite: false,
    }));
    group.add(halo);
    group.add(core);
    group.position.set(entry.x, entry.y, entry.z);
    group.userData.screenRadiusPx = placedScreenRadiusPx(entry.system_count);
    return group;
  }

  function makePlannedSprite(entry) {
    var sprite = new THREE.Sprite(new THREE.SpriteMaterial({ map: plannedTexture, transparent: true, depthWrite: false, opacity: 0.75 }));
    sprite.position.set(entry.x, entry.y, entry.z);
    sprite.userData.screenRadiusPx = PLANNED_PX;
    return sprite;
  }

  // --- Constant-screen-size billboard scaling -----------------------------
  //
  // Recomputed every frame (cheap: a handful of thousand simple vector
  // ops) from each marker's own live distance to the camera, so a
  // marker's ON-SCREEN size stays whatever its own userData.screenRadiusPx
  // says regardless of how far the camera currently is -- see the
  // PLACED_MIN_PX/MAX_PX comment above for why. worldPerScreenPixel here
  // is the same "world units per screen pixel at distance 1" factor
  // updateScaleBar's own worldUnitsPerScreenPixel divides out at the
  // camera's CURRENT orbit radius; this multiplies it back in per-object
  // by that object's own real distance instead.
  function updateMarkerScales() {
    var heightPx = canvasEl.clientHeight || 1;
    var fovRad = THREE.MathUtils.degToRad(camera.fov);
    var worldPerScreenPixelPerUnitDistance = (2 * Math.tan(fovRad / 2)) / heightPx;

    // `anchor` is the world position to measure distance from -- a placed
    // sector's halo/core sprites sit at (0, 0, 0) INSIDE their group, so
    // their own `.position` is not where they are; measuring from it gave
    // the camera's distance to the galactic origin instead, which blew a
    // sector far from the core up into a sphere big enough to fill the
    // whole view once zoomed in on it.
    function screenSizedScale(object3d, screenRadiusPx, sizeMultiplier, anchor) {
      var distance = camera.position.distanceTo(anchor || object3d.position);
      var worldDiameter = screenRadiusPx * 2 * worldPerScreenPixelPerUnitDistance * distance;
      var scaled = worldDiameter * (sizeMultiplier || 1);
      object3d.scale.set(scaled, scaled, 1);
    }

    placedSpritesByKey.forEach(function (group) {
      var px = group.userData.screenRadiusPx;
      screenSizedScale(group.children[0], px, 1.2, group.position); // halo
      screenSizedScale(group.children[1], px, 1.0, group.position); // core
    });
    plannedSpritesByKey.forEach(function (sprite) {
      screenSizedScale(sprite, sprite.userData.screenRadiusPx, 1.0);
    });
    if (highlightSprite.visible) {
      screenSizedScale(highlightSprite, highlightSprite.userData.screenRadiusPx, 1.0);
    }
  }

  // Diffs the live tier against what's already on screen -- keyed so a
  // sector already visible doesn't get torn down and rebuilt (which
  // would flicker/lose its raycast-hit continuity) just because a
  // debounced re-fetch landed while the camera barely moved.
  function syncTier(entries, byKey, keyOf, factory, kind) {
    var seen = new Set();
    entries.forEach(function (entry) {
      var key = keyOf(entry);
      seen.add(key);
      if (byKey.has(key)) {
        return;
      }
      var object3d = factory(entry);
      interactiveGroup.add(object3d);
      byKey.set(key, object3d);
      entryByObject.set(object3d, Object.assign({ kind: kind }, entry));
      object3d.traverse(function (child) {
        if (child !== object3d) {
          entryByObject.set(child, Object.assign({ kind: kind }, entry));
        }
      });
    });
    byKey.forEach(function (object3d, key) {
      if (!seen.has(key)) {
        interactiveGroup.remove(object3d);
        entryByObject.delete(object3d);
        object3d.traverse(function (child) {
          entryByObject.delete(child);
        });
        byKey.delete(key);
      }
    });
  }

  // Shared by densityColor (below) and applyDensity's own per-instance
  // scale -- both want the same "how dense is this, on a 0..1 scale"
  // number, just mapped to a color and a size respectively.
  function densityIntensity(relativeDensity) {
    return Math.max(0, Math.min(1, Math.log2((relativeDensity || 0) + 1) / 3));
  }

  function densityColor(relativeDensity) {
    // Dim, cool color at low density warming toward the accent color at
    // high density -- purely illustrative (see lib/galaxyViewport.py's
    // own docstring), so this is a display choice, not derived physics.
    var base = new THREE.Color(accentColor);
    var dim = new THREE.Color(0x3a3f55);
    return dim.clone().lerp(base, densityIntensity(relativeDensity));
  }

  var densityMatrix = new THREE.Matrix4();
  var densityInstanceColor = new THREE.Color();

  function applyDensity(points) {
    var count = Math.min(points.length, DENSITY_INSTANCE_CAPACITY);
    var baseRadius = Math.max(orbit.radius * DENSITY_RADIUS_FRACTION, 1e-6);
    for (var i = 0; i < count; i++) {
      var point = points[i];
      var t = densityIntensity(point.relative_density);
      var scale = baseRadius * (DENSITY_MIN_SCALE + (DENSITY_MAX_SCALE - DENSITY_MIN_SCALE) * t);
      densityMatrix.makeScale(scale, scale, scale);
      densityMatrix.setPosition(point.x, point.y, point.z);
      densityMesh.setMatrixAt(i, densityMatrix);
      densityMesh.setColorAt(i, densityInstanceColor.copy(densityColor(point.relative_density)));
    }
    densityMesh.count = count;
    densityMesh.instanceMatrix.needsUpdate = true;
    if (densityMesh.instanceColor) {
      densityMesh.instanceColor.needsUpdate = true;
    }
  }

  var PLACED_KEY_OF = function (e) { return "p" + e.id; };
  var PLANNED_KEY_OF = function (e) { return e.shell_index + ":" + e.shell_slot_index; };

  // See pinnedEntry's own comment above selectEntry -- re-inserts the
  // pinned entry into a freshly-fetched tier list if the live query
  // itself didn't happen to include it, so syncTier never drops it.
  function withPinned(list, keyOf, kind) {
    if (!pinnedEntry || pinnedEntry.kind !== kind) {
      return list;
    }
    var pinnedKey = keyOf(pinnedEntry);
    for (var i = 0; i < list.length; i++) {
      if (keyOf(list[i]) === pinnedKey) {
        return list;
      }
    }
    return list.concat([pinnedEntry]);
  }


  // --- Cube tiles ----------------------------------------------------------
  //
  // The map asks for fixed cubes of space ("tiles") rather than "everything
  // within R of the camera target" -- see stellarObjects.galaxyViewport's
  // "Cube tiles" section. Space is an octree: level 0 is one cube
  // tileRootEdgePc on a side centered on the galactic origin, each level
  // halves the edge, and a tile's key is "level/ix/iy/iz". neededTiles()
  // picks the smallest level whose tiles are at least the view radius
  // across (so at most 27 tiles cover the view), plus the smallest tiles
  // (the only ones that list planned slots) around the target when zoomed
  // in. lib/galaxymap3d.py's initial_tile_request does the same in Python
  // for the first frame.
  //
  // A tile's contents depend only on its key and the database's content
  // stamp, so tiles are cached three ways: in memory here, in this
  // browser's localStorage (so a revisit or reload doesn't refetch them),
  // and on the server's disk (lib/tilecache.py) -- only tiles in none of
  // those reach the API and database. Every response carries the current
  // stamp; when it changes (new sectors generated), every cached tile is
  // dropped and refetched.

  var TILE_ROOT = data.tileRootEdgePc || 65536;
  var TILE_MAX_LEVEL = data.tileMaxLevel != null ? data.tileMaxLevel : 12;
  var PLANNED_TILE_EDGE = data.plannedTileMaxEdgePc || 16;
  var PLANNED_VIEW_RADIUS = data.plannedViewRadiusPc || 20;
  var PLANNED_MAX_VIEW_RADIUS = data.plannedMaxViewRadiusPc || 200;
  var FETCH_RADIUS_FACTOR = data.fetchRadiusFactor || 1.6;
  var MAX_TILES_PER_REQUEST = data.maxTilesPerRequest || 128;
  var hasShape = !!data.hasShape;

  function tileLevelForRadius(radius) {
    if (!(radius > 0)) {
      return TILE_MAX_LEVEL;
    }
    var level = Math.floor(Math.log2(TILE_ROOT / radius));
    return Math.max(0, Math.min(TILE_MAX_LEVEL, level));
  }

  function tileEdge(level) {
    return TILE_ROOT / Math.pow(2, level);
  }

  // Nearest first, matching galaxyViewport.tiles_intersecting_sphere.
  function tilesIntersectingSphere(level, center, radius) {
    var edge = tileEdge(level);
    var origin = -TILE_ROOT / 2;
    var span = Math.pow(2, level);
    var c = [center.x, center.y, center.z];
    var lo = [];
    var hi = [];
    for (var axis = 0; axis < 3; axis++) {
      lo.push(Math.max(0, Math.floor((c[axis] - radius - origin) / edge)));
      hi.push(Math.min(span - 1, Math.floor((c[axis] + radius - origin) / edge)));
    }
    var found = [];
    for (var ix = lo[0]; ix <= hi[0]; ix++) {
      for (var iy = lo[1]; iy <= hi[1]; iy++) {
        for (var iz = lo[2]; iz <= hi[2]; iz++) {
          var index = [ix, iy, iz];
          var distanceSq = 0;
          for (var a = 0; a < 3; a++) {
            var boxLo = origin + index[a] * edge;
            var nearest = Math.min(Math.max(c[a], boxLo), boxLo + edge);
            distanceSq += (c[a] - nearest) * (c[a] - nearest);
          }
          if (distanceSq <= radius * radius) {
            found.push({ d: distanceSq, key: level + "/" + ix + "/" + iy + "/" + iz });
          }
        }
      }
    }
    found.sort(function (p, q) { return p.d - q.d; });
    return found.map(function (f) { return f.key; });
  }

  function tileContaining(level, point) {
    var edge = tileEdge(level);
    var origin = -TILE_ROOT / 2;
    var span = Math.pow(2, level);
    var p = [point.x, point.y, point.z];
    var index = p.map(function (v) {
      return Math.max(0, Math.min(span - 1, Math.floor((v - origin) / edge)));
    });
    return level + "/" + index.join("/");
  }

  function neededTiles() {
    var viewRadius = Math.max(MIN_RADIUS, Math.min(MAX_RADIUS * FETCH_RADIUS_FACTOR, orbit.radius * FETCH_RADIUS_FACTOR));
    var level = tileLevelForRadius(viewRadius);
    var keys = tilesIntersectingSphere(level, target, viewRadius);
    var plannedKeys = [];
    var plannedRadius = 0;
    if (viewRadius <= PLANNED_MAX_VIEW_RADIUS) {
      plannedRadius = Math.min(viewRadius, PLANNED_VIEW_RADIUS);
      plannedKeys = tilesIntersectingSphere(tileLevelForRadius(PLANNED_TILE_EDGE), target, plannedRadius);
      plannedKeys.forEach(function (key) {
        if (keys.indexOf(key) < 0) {
          keys.push(key);
        }
      });
    }
    var densityKey = hasShape && viewRadius > PLANNED_MAX_VIEW_RADIUS ? tileContaining(level, target) : null;
    return { keys: keys, plannedKeys: plannedKeys, plannedRadius: plannedRadius, densityKey: densityKey };
  }

  // --- Tile caches ---------------------------------------------------------

  var TILE_MEMORY_MAX = 2000;
  var DENSITY_MEMORY_MAX = 24;
  var STORAGE_PREFIX = "planetgen:tile:" + data.db + ":";
  var currentStamp = (data.initial && data.initial.stamp) || "";
  var tileMemory = new Map();
  var densityMemory = new Map();

  // Map iteration order is insertion order, so re-inserting on every hit
  // and evicting from the front makes these LRU caches.
  function memoryGet(cache, key) {
    if (!cache.has(key)) {
      return undefined;
    }
    var value = cache.get(key);
    cache.delete(key);
    cache.set(key, value);
    return value;
  }

  function memorySet(cache, key, value, max) {
    cache.delete(key);
    cache.set(key, value);
    while (cache.size > max) {
      cache.delete(cache.keys().next().value);
    }
  }

  // localStorage can be missing, full, or throw on any access (private
  // browsing, blocked site data) -- every use is wrapped, and the map
  // works the same without it, just refetching from the server cache.
  function storage() {
    try {
      return window.localStorage || null;
    } catch (err) {
      return null;
    }
  }

  // Drops this database's stored tiles, except the current stamp's when
  // keepCurrent is set.
  function purgeStoredTiles(keepCurrent) {
    var store = storage();
    if (!store) {
      return;
    }
    try {
      var doomed = [];
      for (var i = 0; i < store.length; i++) {
        var name = store.key(i);
        if (name && name.indexOf(STORAGE_PREFIX) === 0) {
          if (!keepCurrent || name.indexOf(STORAGE_PREFIX + currentStamp + ":") !== 0) {
            doomed.push(name);
          }
        }
      }
      doomed.forEach(function (name) { store.removeItem(name); });
    } catch (err) {
      // Nothing more to do -- storage just stays as it was.
    }
  }

  function storageName(key) {
    return STORAGE_PREFIX + currentStamp + ":" + key;
  }

  function storedTile(key) {
    var store = storage();
    if (!store || !currentStamp) {
      return undefined;
    }
    try {
      var raw = store.getItem(storageName(key));
      return raw ? JSON.parse(raw) : undefined;
    } catch (err) {
      return undefined;
    }
  }

  function storeTile(key, tile) {
    var store = storage();
    if (!store || !currentStamp) {
      return;
    }
    var raw = JSON.stringify(tile);
    try {
      store.setItem(storageName(key), raw);
    } catch (err) {
      // Probably full: drop other stamps' leftovers and try once more,
      // then drop this database's tiles entirely and try a last time.
      purgeStoredTiles(true);
      try {
        store.setItem(storageName(key), raw);
      } catch (err2) {
        purgeStoredTiles(false);
        try {
          store.setItem(storageName(key), raw);
        } catch (err3) {
          // Leave it memory-only.
        }
      }
    }
  }

  function getTile(key) {
    var tile = memoryGet(tileMemory, key);
    if (tile === undefined) {
      tile = storedTile(key);
      if (tile !== undefined) {
        memorySet(tileMemory, key, tile, TILE_MEMORY_MAX);
      }
    }
    return tile;
  }

  function putTile(key, tile) {
    memorySet(tileMemory, key, tile, TILE_MEMORY_MAX);
    storeTile(key, tile);
  }

  // Density clouds (~100 KB each) stay in memory only; they would crowd
  // tiles out of localStorage's few-megabyte budget, and the server's
  // disk cache already serves them without touching the database.
  function getDensity(key) {
    return memoryGet(densityMemory, key);
  }

  function putDensity(key, points) {
    memorySet(densityMemory, key, points, DENSITY_MEMORY_MAX);
  }

  function adoptStamp(stamp) {
    if (!stamp || stamp === currentStamp) {
      return false;
    }
    currentStamp = stamp;
    tileMemory.clear();
    densityMemory.clear();
    purgeStoredTiles(true);
    return true;
  }

  function absorb(payload) {
    if (!payload) {
      return false;
    }
    var stampChanged = adoptStamp(payload.stamp);
    if (payload.has_shape != null) {
      hasShape = !!payload.has_shape;
    }
    Object.keys(payload.tiles || {}).forEach(function (key) {
      putTile(key, payload.tiles[key]);
    });
    if (payload.density && payload.density.key) {
      putDensity(payload.density.key, payload.density.points || []);
    }
    return stampChanged;
  }

  purgeStoredTiles(true);
  absorb(data.initial);

  // --- Drawing from tiles --------------------------------------------------

  function withinSphere(entry, radius) {
    var dx = entry.x - target.x;
    var dy = entry.y - target.y;
    var dz = entry.z - target.z;
    return dx * dx + dy * dy + dz * dz <= radius * radius;
  }

  // Draws whatever of the needed tiles is already cached; returns what's
  // still missing.
  function renderFromCache(need) {
    var placed = [];
    var planned = [];
    var missing = [];
    var plannedSet = new Set(need.plannedKeys);
    need.keys.forEach(function (key) {
      var tile = getTile(key);
      if (tile === undefined) {
        missing.push(key);
        return;
      }
      Array.prototype.push.apply(placed, tile.placed || []);
      if (plannedSet.has(key)) {
        (tile.planned || []).forEach(function (entry) {
          if (withinSphere(entry, need.plannedRadius)) {
            planned.push(entry);
          }
        });
      }
    });
    syncTier(withPinned(placed, PLACED_KEY_OF, "placed"), placedSpritesByKey, PLACED_KEY_OF, makePlacedSprite, "placed");
    syncTier(withPinned(planned, PLANNED_KEY_OF, "planned"), plannedSpritesByKey, PLANNED_KEY_OF, makePlannedSprite, "planned");

    var densityMissing = false;
    if (!need.densityKey) {
      applyDensity([]);
    } else {
      var points = getDensity(need.densityKey);
      if (points === undefined) {
        // Keep the previous cloud on screen until the new one arrives.
        densityMissing = true;
      } else {
        applyDensity(points);
      }
    }
    return { tiles: missing, density: densityMissing ? need.densityKey : null };
  }

  // --- Live tile fetching --------------------------------------------------

  var FETCH_DEBOUNCE_MS = 300;
  var fetchTimer = null;
  var activeAbort = null;

  // Caps how often a click/double-click/zoom-button interaction (every
  // caller that passes scheduleFetch(true) -- see each one below) can
  // actually trigger an IMMEDIATE live fetch: at most MAX_CLICKS_PER_SECOND
  // per second. `doFetch`'s own activeAbort.abort() only stops the
  // BROWSER from waiting on a superseded response -- it doesn't reliably
  // stop the server from finishing a request it already started, so
  // rapid clicking still costs the server work per click even when every
  // earlier response gets thrown away client-side. A click/double-click's
  // own visual effect (the camera recentering/zooming, via centerOn/
  // centerAndZoom) is never throttled here, only the network fetch that
  // follows it -- clicking faster than the cap still feels instant, it
  // just falls back to the standard debounced delay below instead of
  // firing right away, so a rapid burst still settles on exactly one
  // fetch shortly after it stops. (With tiles, an aborted request's work
  // isn't wasted either: the server caches whatever it computed.)
  var MAX_CLICKS_PER_SECOND = 4;
  var MIN_MS_BETWEEN_IMMEDIATE_FETCHES = 1000 / MAX_CLICKS_PER_SECOND;
  var lastImmediateFetchAt = 0;

  function scheduleFetch(immediate) {
    if (fetchTimer) {
      clearTimeout(fetchTimer);
    }
    var delay = FETCH_DEBOUNCE_MS;
    if (immediate) {
      var now = Date.now();
      if (now - lastImmediateFetchAt >= MIN_MS_BETWEEN_IMMEDIATE_FETCHES) {
        lastImmediateFetchAt = now;
        delay = 0;
      }
      // else: rate-limited -- falls through to the debounced delay above
      // instead of a bare no-op, so state still eventually syncs.
    }
    fetchTimer = setTimeout(doFetch, delay);
  }

  function doFetch() {
    fetchTimer = null;
    var need = neededTiles();
    var missing = renderFromCache(need);
    if (!missing.tiles.length && !missing.density) {
      return;
    }
    if (activeAbort) {
      activeAbort.abort();
    }
    var controller = typeof AbortController !== "undefined" ? new AbortController() : null;
    activeAbort = controller;
    var params = new URLSearchParams({
      db: data.db, tiles: missing.tiles.slice(0, MAX_TILES_PER_REQUEST).join(","),
    });
    if (missing.density) {
      params.set("density", missing.density);
    }
    fetch(data.fetchPath + "?" + params.toString(), controller ? { signal: controller.signal } : undefined)
      .then(function (response) {
        if (!response.ok) {
          throw new Error("galaxy_tiles.py returned " + response.status);
        }
        return response.json();
      })
      .then(function (payload) {
        if (activeAbort === controller) {
          activeAbort = null;
        }
        var stampChanged = absorb(payload);
        var stillMissing = renderFromCache(neededTiles());
        // A new stamp dropped every cached tile, and a view needing more
        // than one request's worth of tiles has more to fetch -- either
        // way, go again (tiles already fetched are cached by now).
        if (stampChanged || stillMissing.tiles.length || stillMissing.density) {
          scheduleFetch(false);
        }
      })
      .catch(function (err) {
        if (err && err.name === "AbortError") {
          return;
        }
        // A transient fetch failure just leaves the currently-drawn
        // content in place -- the next camera move retries automatically,
        // and there's no useful place to surface a network error inside
        // this canvas.
      });
  }

  // First frame: normally every tile is already cached (the page embeds
  // them), so this only fetches if the page's own tile fetch came up short.
  var initialMissing = renderFromCache(neededTiles());
  if (initialMissing.tiles.length || initialMissing.density) {
    scheduleFetch(false);
  }

  // --- Pointer/keyboard interaction --------------------------------------

  var dragging = false;
  var dragDistance = 0;
  var lastClientX = 0;
  var lastClientY = 0;
  var suppressNextClick = false;

  canvasEl.addEventListener("pointerdown", function (event) {
    if (event.button !== 0) {
      return;
    }
    dragging = true;
    dragDistance = 0;
    lastClientX = event.clientX;
    lastClientY = event.clientY;
    try {
      canvasEl.setPointerCapture(event.pointerId);
    } catch (err) {
      // Not essential -- dragging still works via ordinary bubbling.
    }
  });

  canvasEl.addEventListener("pointermove", function (event) {
    hideTooltip();
    if (!dragging) {
      scheduleTooltipCheck(event.clientX, event.clientY);
      return;
    }
    var deltaX = event.clientX - lastClientX;
    var deltaY = event.clientY - lastClientY;
    dragDistance += Math.abs(deltaX) + Math.abs(deltaY);
    lastClientX = event.clientX;
    lastClientY = event.clientY;

    orbit.theta -= deltaX * ROTATE_SENSITIVITY;
    orbit.phi = Math.max(MIN_POLAR, Math.min(MAX_POLAR, orbit.phi - deltaY * ROTATE_SENSITIVITY));
    applyCamera();
  });

  function endDrag(event) {
    if (!dragging) {
      return;
    }
    dragging = false;
    if (dragDistance > DRAG_CLICK_THRESHOLD_PX) {
      suppressNextClick = true;
      setTimeout(function () {
        suppressNextClick = false;
      }, 0);
      scheduleFetch(false);
    }
    dragDistance = 0;
    try {
      canvasEl.releasePointerCapture(event.pointerId);
    } catch (err) {
      // Already released/invalid.
    }
  }
  canvasEl.addEventListener("pointerup", endDrag);
  canvasEl.addEventListener("pointercancel", endDrag);

  canvasEl.addEventListener(
    "wheel",
    function (event) {
      event.preventDefault();
      setRadius(orbit.radius * (event.deltaY < 0 ? 1 / WHEEL_ZOOM_RATIO : WHEEL_ZOOM_RATIO));
      scheduleFetch(false);
    },
    { passive: false }
  );

  canvasEl.addEventListener("keydown", function (event) {
    var key = event.key;
    if (key !== "ArrowLeft" && key !== "ArrowRight" && key !== "ArrowUp" && key !== "ArrowDown") {
      return;
    }
    event.preventDefault();
    if (key === "ArrowLeft") orbit.theta += KEY_ROTATE_STEP;
    if (key === "ArrowRight") orbit.theta -= KEY_ROTATE_STEP;
    if (key === "ArrowUp") orbit.phi = Math.max(MIN_POLAR, orbit.phi - KEY_ROTATE_STEP);
    if (key === "ArrowDown") orbit.phi = Math.min(MAX_POLAR, orbit.phi + KEY_ROTATE_STEP);
    applyCamera();
    scheduleFetch(false);
  });

  var raycaster = new THREE.Raycaster();

  function ndcFromClientPoint(clientX, clientY) {
    var rect = canvasEl.getBoundingClientRect();
    if (rect.width === 0 || rect.height === 0) {
      return null;
    }
    return new THREE.Vector2(
      ((clientX - rect.left) / rect.width) * 2 - 1,
      -((clientY - rect.top) / rect.height) * 2 + 1
    );
  }

  function entryAtClientPoint(clientX, clientY) {
    var ndc = ndcFromClientPoint(clientX, clientY);
    if (!ndc) {
      return null;
    }
    raycaster.setFromCamera(ndc, camera);
    var hits = raycaster.intersectObjects(interactiveGroup.children, true);
    for (var i = 0; i < hits.length; i++) {
      var entry = entryByObject.get(hits[i].object);
      if (entry) {
        return { entry: entry, point: hits[i].point };
      }
    }
    return null;
  }

  // Empty-space click target: intersects an invisible sphere of the
  // camera's OWN current orbit radius, centered on the current target --
  // an approximation of "the depth the camera is presently looking at",
  // so clicking past empty space still lands somewhere reasonable in 3D
  // rather than needing a ground plane a free-flying galaxy scene has no
  // natural equivalent of.
  function depthPointAtClientPoint(clientX, clientY) {
    var ndc = ndcFromClientPoint(clientX, clientY);
    if (!ndc) {
      return null;
    }
    raycaster.setFromCamera(ndc, camera);
    var sphere = new THREE.Sphere(target, Math.max(orbit.radius, MIN_RADIUS));
    var hitPoint = new THREE.Vector3();
    var hit = raycaster.ray.intersectSphere(sphere, hitPoint);
    return hit ? hitPoint : null;
  }

  // selectEntry alone (no target/radius change) just updates the
  // highlight/info panel -- shared by centerOn/centerAndZoom below so a
  // click/double-click on empty space (entry === null) leaves whatever
  // was last selected showing, the same "recentering doesn't clear your
  // selection" behavior the old single zoomToward function had.
  //
  // Also PINS the selected entry (see pinnedEntry/renderFromCache above): the
  // live view's own tile set shrinks as orbit.radius shrinks
  // (neededTiles' own FETCH_RADIUS_FACTOR), so once you're centering/
  // zooming in toward one specific sector, a click/double-click that
  // landed even slightly off that sector's own exact stored position
  // (easy to do from far out, where its marker is only a handful of
  // screen pixels wide) could otherwise cause a later, smaller-radius
  // re-fetch to legitimately no longer include it -- and syncTier
  // removes anything not present in a fresh fetch, making the very
  // sector you just selected and are flying toward vanish outright.
  // Pinning keeps it in the scene regardless of what the live viewport
  // query returns, for as long as it stays selected.
  var pinnedEntry = null;

  function selectEntry(point, entry) {
    if (!entry) {
      return;
    }
    pinnedEntry = entry;
    highlightPosition(point.x, point.y, point.z, entry.kind === "placed" ? placedScreenRadiusPx(entry.system_count) : PLANNED_PX);
    if (entry.kind === "placed") {
      showPlacedInfo(entry);
    } else {
      showPlannedInfo(entry);
    }
  }

  // Single click: re-centers the view on the clicked point (or selects/
  // navigates-info for a clicked dot) WITHOUT zooming -- deliberately not
  // "click to zoom" any more (see this file's own module docstring's
  // former "Click-to-zoom is LOGARITHMIC" note, now double-click's own
  // job below). Centering alone, with no zoom commitment, is what makes
  // it possible to walk the camera across the galaxy toward a small/
  // distant dot over several clicks without a bad click also zooming
  // into empty space you didn't mean to approach.
  function centerOn(point, entry) {
    target.copy(point);
    applyCamera();
    updateScaleBar();
    selectEntry(point, entry);
    scheduleFetch(true);
  }

  // Double click: the old single-click behavior -- centers AND zooms in
  // by one clickZoomFactor step, in the same motion a click used to.
  function centerAndZoom(point, entry) {
    var factor = clickZoomFactor(orbit.radius);
    target.copy(point);
    orbit.radius = Math.max(MIN_RADIUS, orbit.radius / factor);
    applyCamera();
    updateScaleBar();
    selectEntry(point, entry);
    scheduleFetch(true);
  }

  // Resolves a click/double-click's target point the same way for both:
  // a hit dot's own position, or (empty space) the depth-sphere fallback
  // -- shared so centerOn/centerAndZoom above never have to duplicate the
  // raycast-then-fall-back logic.
  function resolveClickTarget(clientX, clientY) {
    var hit = entryAtClientPoint(clientX, clientY);
    if (hit) {
      return hit;
    }
    var depthPoint = depthPointAtClientPoint(clientX, clientY);
    return depthPoint ? { point: depthPoint, entry: null } : null;
  }

  canvasEl.addEventListener("click", function (event) {
    if (suppressNextClick) {
      suppressNextClick = false;
      return;
    }
    var resolved = resolveClickTarget(event.clientX, event.clientY);
    if (resolved) {
      centerOn(resolved.point, resolved.entry);
    }
  });

  canvasEl.addEventListener("dblclick", function (event) {
    if (suppressNextClick) {
      suppressNextClick = false;
      return;
    }
    var resolved = resolveClickTarget(event.clientX, event.clientY);
    if (resolved) {
      centerAndZoom(resolved.point, resolved.entry);
    }
  });

  // No right-click action any more (see this file's own module
  // docstring) -- the browser's own default context menu is left alone,
  // rather than preventDefault()-ing it for nothing.

  var controlsEl = document.getElementById("galaxymap3d-controls");
  if (controlsEl) {
    controlsEl.querySelectorAll("[data-action]").forEach(function (button) {
      button.addEventListener("click", function () {
        var action = button.dataset.action;
        if (action === "zoom-in") {
          setRadius(orbit.radius / clickZoomFactor(orbit.radius));
          scheduleFetch(true);
        } else if (action === "zoom-out") {
          setRadius(orbit.radius * clickZoomFactor(orbit.radius));
          scheduleFetch(true);
        } else if (action === "reset") {
          resetView();
        }
      });
    });
  }

  // --- Hover tooltip ---------------------------------------------------

  var tooltipEl = null;
  function ensureTooltipEl() {
    if (!tooltipEl && viewport) {
      tooltipEl = document.createElement("div");
      tooltipEl.className = "galaxymap3d-tooltip";
      tooltipEl.hidden = true;
      viewport.appendChild(tooltipEl);
    }
    return tooltipEl;
  }

  function hideTooltip() {
    if (tooltipEl) {
      tooltipEl.hidden = true;
    }
  }

  var tooltipPending = false;
  function scheduleTooltipCheck(clientX, clientY) {
    if (tooltipPending) {
      return;
    }
    tooltipPending = true;
    requestAnimationFrame(function () {
      tooltipPending = false;
      var hit = entryAtClientPoint(clientX, clientY);
      var el = ensureTooltipEl();
      if (!hit || !el || !viewport) {
        hideTooltip();
        return;
      }
      var entry = hit.entry;
      el.textContent =
        entry.kind === "placed"
          ? (entry.name || "Unnamed sector") + " — " + (entry.system_count != null ? entry.system_count : 0) + " system" + (entry.system_count === 1 ? "" : "s")
          : (entry.designation || "Not yet generated") + " — not yet generated";
      var rect = viewport.getBoundingClientRect();
      el.style.left = (clientX - rect.left + 14) + "px";
      el.style.top = (clientY - rect.top + 14) + "px";
      el.hidden = false;
    });
  }

  canvasEl.addEventListener("pointerleave", hideTooltip);

  // --- Scale bar -----------------------------------------------------------

  var scaleEl = document.getElementById("galaxymap3d-scale");

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

  function formatPc(value) {
    if (value >= 100) return Math.round(value).toLocaleString() + " pc";
    if (value >= 1) return Math.round(value * 10) / 10 + " pc";
    return Math.round(value * 1000) / 1000 + " pc";
  }

  function updateScaleBar() {
    if (!scaleEl) {
      return;
    }
    var fovRad = THREE.MathUtils.degToRad(camera.fov);
    var heightPx = canvasEl.clientHeight || 1;
    var pcPerScreenPx = (2 * orbit.radius * Math.tan(fovRad / 2)) / heightPx;
    var nicePc = niceScaleValue(70 * pcPerScreenPx);
    if (!nicePc) {
      return;
    }
    scaleEl.textContent = "≈ " + formatPc(nicePc) + " (reference)";
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
    updateMarkerScales();
    renderer.render(scene, camera);
  })();
}

if (canvas && sceneData) {
  initGalaxyMap3d(canvas, sceneData);
}
