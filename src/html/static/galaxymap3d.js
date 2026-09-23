// html/static/galaxymap3d.js
//
// Renders the interactive 3D Galaxy Map (lib/galaxymap3d.py) as a real
// WebGL scene (three.js, vendored at ./vendor/three.module.min.js -- see
// sectormap.js's own docstring for why it's vendored rather than loaded
// from a CDN). Unlike sectormap.js's scene (every star/cloud baked into
// one JSON block at page load, camera always orbiting a fixed origin),
// this camera's own orbit TARGET moves freely through the galaxy -- most
// of what's drawn is fetched live from galaxy_view.py as the camera
// moves, never baked into the page beyond the very first frame
// (#galaxymap3d-data's own "initial" payload).
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
// Three content tiers per live fetch (queryDb.galaxy_view's own
// placed/planned/density lists -- see stellarObjects.galaxyViewport's
// module docstring):
//   - placed:  real, already-generated sectors -- bright sprites, sized
//              by system_count, click selects + navigates.
//   - planned: real, not-yet-generated qualifying addresses -- small dim
//              sprites, click selects (shows the copyable designation/
//              CLI snippet) but never navigates.
//   - density: an illustrative point cloud (THREE.Points, not individual
//              sprites -- not interactive, no identity to track).
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

  var dl = document.createElement("dl");
  addField(dl, "Systems", entry.system_count != null ? entry.system_count : 0);
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
var PLACED_MIN_PX = 4.0;
var PLACED_MAX_PX = 16.0;
var PLACED_CORE_FILL = "#fff6df";
var PLACED_CORE_STROKE = "#caa54d";
var PLACED_HALO_FILL = "#ffd88a";

function placedScreenRadiusPx(systemCount) {
  var count = systemCount || 0;
  return Math.max(PLACED_MIN_PX, Math.min(PLACED_MAX_PX, PLACED_MIN_PX + 2.5 * Math.sqrt(count)));
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
    placedCoreTexture = makeDotTexture(PLACED_CORE_FILL, PLACED_CORE_STROKE);
    placedHaloTexture = makeDotTexture(PLACED_HALO_FILL, null);
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

  // sizeAttenuation: false -- a constant SCREEN-pixel point size
  // regardless of camera distance (three.js's PointsMaterial supports
  // this natively, unlike Sprite -- see updateMarkerScales below for how
  // placed/planned markers get the same effect). This illustrative cloud
  // needs to read as a recognizable galaxy shape from any zoom level,
  // including the full-galaxy starting view thousands of parsecs out,
  // where a world-space point size would shrink to sub-pixel and vanish.
  var densityGeometry = new THREE.BufferGeometry();
  var densityMaterial = new THREE.PointsMaterial({
    size: 2.2, sizeAttenuation: false, vertexColors: true,
    transparent: true, opacity: 0.55, depthWrite: false,
  });
  var densityPoints = new THREE.Points(densityGeometry, densityMaterial);
  scene.add(densityPoints);

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
    var halo = new THREE.Sprite(new THREE.SpriteMaterial({ map: placedHaloTexture, transparent: true, depthWrite: false, opacity: 0.45 }));
    var core = new THREE.Sprite(new THREE.SpriteMaterial({ map: placedCoreTexture, transparent: true, depthWrite: false }));
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

    function screenSizedScale(object3d, screenRadiusPx, sizeMultiplier) {
      var distance = camera.position.distanceTo(object3d.position);
      var worldDiameter = screenRadiusPx * 2 * worldPerScreenPixelPerUnitDistance * distance;
      var scaled = worldDiameter * (sizeMultiplier || 1);
      object3d.scale.set(scaled, scaled, 1);
    }

    placedSpritesByKey.forEach(function (group) {
      var px = group.userData.screenRadiusPx;
      screenSizedScale(group.children[0], px, 1.2); // halo
      screenSizedScale(group.children[1], px, 1.0); // core
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

  function densityColor(relativeDensity) {
    // Dim, cool color at low density warming toward the accent color at
    // high density -- purely illustrative (see lib/galaxyViewport.py's
    // own docstring), so this is a display choice, not derived physics.
    var t = Math.max(0, Math.min(1, Math.log2((relativeDensity || 0) + 1) / 3));
    var base = new THREE.Color(accentColor);
    var dim = new THREE.Color(0x3a3f55);
    return dim.clone().lerp(base, t);
  }

  function applyDensity(points) {
    var positions = new Float32Array(points.length * 3);
    var colors = new Float32Array(points.length * 3);
    for (var i = 0; i < points.length; i++) {
      positions[i * 3] = points[i].x;
      positions[i * 3 + 1] = points[i].y;
      positions[i * 3 + 2] = points[i].z;
      var color = densityColor(points[i].relative_density);
      colors[i * 3] = color.r;
      colors[i * 3 + 1] = color.g;
      colors[i * 3 + 2] = color.b;
    }
    densityGeometry.setAttribute("position", new THREE.BufferAttribute(positions, 3));
    densityGeometry.setAttribute("color", new THREE.BufferAttribute(colors, 3));
    densityGeometry.computeBoundingSphere();
  }

  function applyView(view) {
    syncTier(view.placed || [], placedSpritesByKey, function (e) { return "p" + e.id; }, makePlacedSprite, "placed");
    syncTier(
      view.planned || [], plannedSpritesByKey,
      function (e) { return e.shell_index + ":" + e.shell_slot_index; },
      makePlannedSprite, "planned",
    );
    applyDensity(view.density || []);
  }

  applyView(data.initial || { placed: [], planned: [], density: [] });

  // --- Live viewport fetching ---------------------------------------------

  var FETCH_DEBOUNCE_MS = 300;
  var FETCH_RADIUS_FACTOR = 1.6;
  var fetchTimer = null;
  var activeAbort = null;

  function scheduleFetch(immediate) {
    if (fetchTimer) {
      clearTimeout(fetchTimer);
    }
    fetchTimer = setTimeout(doFetch, immediate ? 0 : FETCH_DEBOUNCE_MS);
  }

  function doFetch() {
    if (activeAbort) {
      activeAbort.abort();
    }
    var controller = typeof AbortController !== "undefined" ? new AbortController() : null;
    activeAbort = controller;
    var radius = Math.max(MIN_RADIUS, Math.min(MAX_RADIUS, orbit.radius * FETCH_RADIUS_FACTOR));
    var params = new URLSearchParams({
      db: data.db, cx: target.x, cy: target.y, cz: target.z, radius_pc: radius,
    });
    fetch(data.fetchPath + "?" + params.toString(), controller ? { signal: controller.signal } : undefined)
      .then(function (response) {
        if (!response.ok) {
          throw new Error("galaxy_view.py returned " + response.status);
        }
        return response.json();
      })
      .then(applyView)
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

  function zoomToward(point, entry) {
    var factor = clickZoomFactor(orbit.radius);
    var newRadius = Math.max(MIN_RADIUS, orbit.radius / factor);
    target.copy(point);
    orbit.radius = newRadius;
    applyCamera();
    updateScaleBar();
    if (entry) {
      highlightPosition(point.x, point.y, point.z, entry.kind === "placed" ? placedScreenRadiusPx(entry.system_count) : PLANNED_PX);
      if (entry.kind === "placed") {
        showPlacedInfo(entry);
      } else {
        showPlannedInfo(entry);
      }
    }
    scheduleFetch(true);
  }

  canvasEl.addEventListener("click", function (event) {
    if (suppressNextClick) {
      suppressNextClick = false;
      return;
    }
    var hit = entryAtClientPoint(event.clientX, event.clientY);
    if (hit) {
      zoomToward(hit.point, hit.entry);
      return;
    }
    var depthPoint = depthPointAtClientPoint(event.clientX, event.clientY);
    if (depthPoint) {
      zoomToward(depthPoint, null);
    }
  });

  canvasEl.addEventListener("contextmenu", function (event) {
    event.preventDefault();
    var factor = clickZoomFactor(orbit.radius);
    setRadius(orbit.radius * factor);
    scheduleFetch(true);
  });

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
