// html/static/sectormap.js
//
// The sector page's 3D Sector Map: planetgen/web/maps/starmap.py's
// `#starmap-data` (a `<script type="application/json">` block) drawn as a
// real WebGL scene (three.js, vendored at `static/vendor/three.module.min.js`
// -- see that directory's `THIRD_PARTY_NOTICES.txt` for why it's vendored
// rather than loaded from a CDN). The sector itself (its stars, phenomena,
// neighbors, outline and compass) is built by ./sectorscene.js, which the
// Galaxy Map uses too for its sector stage (MAP.66); this page adds an
// orbit camera round it with drag-to-rotate, scroll/button-to-zoom, the
// scale bar, the controls, the screen-reader list and the "Show on map"
// buttons, and picks, hovers and fills the info panel through
// ./mappick.js (MAP.65).

// Sibling modules are imported with this module's own `?v=<version>`
// query (planetgen/web/lib/fmt.py's `static_url`), so they are cached and
// refreshed with the page's script. A plain static `import "./x.js"`
// would drop the query: an update could then leave a stale copy cached,
// and a page that also loaded the same file by its versioned URL would
// get a second, separate instance of it.
const VERSION_QUERY = new URL(import.meta.url).search;
const THREE = await import(`./vendor/three.module.min.js${VERSION_QUERY}`);
const { formatDistanceLy } = await import(`./distance.js${VERSION_QUERY}`);
const MC = await import(`./mapcontrol.js${VERSION_QUERY}`);
const {
  cssVar, fitRendererToCanvas, isLightBackground, niceScaleValue, readSceneData, watchResize, worldUnitsPerPixel,
} = await import(`./mapcore.js${VERSION_QUERY}`);
const { createPicker, createRing, createTooltip, infoPanelOf } = await import(`./mappick.js${VERSION_QUERY}`);
const {
  buildSectorScene, entryLabel, infoSpec, POINT_CLOSE_GROWTH, tooltipText,
} = await import(`./sectorscene.js${VERSION_QUERY}`);

var canvas = document.getElementById("starmap-canvas");
var dataEl = document.getElementById("starmap-data");

var sceneData = readSceneData(dataEl);

// --- Scene setup ---------------------------------------------------------

// `options.renderer` stands in for the WebGL renderer: the page never
// passes it; tests/js/sectormap.test.mjs does, where there is no WebGL.
export function initStarmap(canvasEl, data, options) {
  var viewport = canvasEl.closest(".starmap-viewport");

  var renderer = (options && options.renderer)
    || new THREE.WebGLRenderer({ canvas: canvasEl, antialias: true, alpha: true, logarithmicDepthBuffer: true });
  renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
  renderer.setClearColor(0x000000, 0);
  if (THREE.SRGBColorSpace) {
    renderer.outputColorSpace = THREE.SRGBColorSpace;
  }

  var scene = new THREE.Scene();

  var FOV_DEG = 45;
  var camera = new THREE.PerspectiveCamera(FOV_DEG, 1, 1, 5e6);

  // The world-unit distance at which a sphere of radius `sceneHalfPx`
  // (planetgen/web/maps/starmap.py's own scale reference -- every position/radius in
  // `data` is expressed in these units) exactly fills the frame
  // vertically -- this is what "zoom = 1" means for a real camera, the
  // direct replacement for the old CSS version's `scale(1)`.
  var sceneHalfPx = data.sceneHalfPx || 160;
  var referenceDistance = sceneHalfPx / Math.tan(THREE.MathUtils.degToRad(FOV_DEG / 2));

  var defaultZoom = data.defaultZoom > 0 && data.defaultZoom <= 1 ? data.defaultZoom : 1;
  // Always room to zoom out past the opening view: a crowded sector opens
  // at starmap.py's 0.2 floor, where a fixed 0.2 minimum left the - button
  // and Reset view with nothing to do (MAP.114).
  var MIN_ZOOM = Math.min(0.2, defaultZoom / 2);
  var MAX_ZOOM = 2.5;
  var ZOOM_STEP = 0.15;
  var WHEEL_ZOOM_STEP = 0.08;
  var KEY_ROTATE_STEP = THREE.MathUtils.degToRad(6);
  var ROTATE_SENSITIVITY = THREE.MathUtils.degToRad(0.4); // radians per pixel of drag
  var DRAG_CLICK_THRESHOLD_PX = 4;
  var MIN_POLAR = THREE.MathUtils.degToRad(2);
  var MAX_POLAR = THREE.MathUtils.degToRad(178);
  // The zoom policy (mapcontrol.js): a short range, zoom MIN_ZOOM to
  // MAX_ZOOM, as camera distances from the zoom-1 distance.
  var ZOOM_POLICY = MC.zoomPolicy(MC.ZOOM_RANGE, 1 / MAX_ZOOM, 1 / MIN_ZOOM);

  function clampPolar(phi) {
    return Math.max(MIN_POLAR, Math.min(MAX_POLAR, phi));
  }

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
    spherical.radius = MC.clampDistance(ZOOM_POLICY, referenceDistance / Math.max(zoom, 1e-6), referenceDistance);
    applyCamera();
    updateScaleBar();
  }

  function resetView() {
    spherical.set(referenceDistance / defaultZoom, defaultPolar, defaultAzimuth);
    applyCamera();
    updateScaleBar();
  }

  var accentColor = cssVar("--accent", "#4f5fe8");

  // How much the points have grown with the camera's zoom: their own size
  // up to default zoom, POINT_CLOSE_GROWTH times it at MAX_ZOOM.
  function pointSizeScale() {
    var share = THREE.MathUtils.clamp((currentZoom() - 1) / (MAX_ZOOM - 1), 0, 1);
    return 1 + (POINT_CLOSE_GROWTH - 1) * share;
  }

  var sector = buildSectorScene(data, {
    accentColor: accentColor, lightBackground: isLightBackground(), pixelRatio: renderer.getPixelRatio(),
    sizeScale: pointSizeScale,
  });
  scene.add(sector.group);

  // The selection ring and, fainter, the hover ring (mappick.js): a
  // point of light keeps its ring a size on screen, a body or cloud gets
  // one round it in the scene.
  var selectionRing = createRing(scene, camera, canvasEl, accentColor);
  var hoverRing = createRing(scene, camera, canvasEl, accentColor, { opacity: 0.45 });
  var highlightedPoint = null;
  var hoveredPoint = null;
  var infoPanel = infoPanelOf(document.getElementById("starmap-info"));
  var tooltip = createTooltip(document.getElementById("starmap-tooltip"));

  function ringAround(ring, entry) {
    ring.at(entry.x, entry.y, entry.z, sector.ringSize(entry));
    return entry.light ? entry : null;
  }

  function highlightEntry(entry) {
    highlightedPoint = ringAround(selectionRing, entry);
  }

  // A point's ring grows with it as the view closes in.
  function updatePointHighlight() {
    if (highlightedPoint) selectionRing.setPx(sector.ringSize(highlightedPoint).px);
    if (hoveredPoint) hoverRing.setPx(sector.ringSize(hoveredPoint).px);
  }

  // "Mark rogue planets": the rings grow with the bigger points.
  function setRoguesMarked(marked) {
    sector.setRoguesMarked(marked);
    updatePointHighlight();
  }

  function selectEntry(entry) {
    if (!entry) {
      return;
    }
    if (infoPanel) infoPanel.show(infoSpec(entry, data));
    highlightEntry(entry);
  }

  function hoverEntry(entry, event) {
    if (!entry) {
      hoverRing.hide();
      hoveredPoint = null;
      tooltip.hide();
      canvasEl.style.cursor = "";
      return;
    }
    hoveredPoint = ringAround(hoverRing, entry);
    tooltip.show(tooltipText(entry), event.clientX, event.clientY);
    canvasEl.style.cursor = "pointer";
  }

  // --- Accessible fallback list -----------------------------------------
  //
  // A canvas has no focusable children of its own the way the old CSS
  // version's real per-star `<div role="button">`s were, so this is what
  // keeps every star/cloud/neighbor reachable by keyboard/screen reader
  // without needing 3D hit-testing or focus management inside the canvas
  // itself -- a visually hidden button per entry, in the same list order
  // the scene data arrived in.
  if (viewport) {
    var list = document.createElement("ul");
    list.className = "starmap-sr-list sr-only";
    sector.entries.forEach(function (entry) {
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

  canvasEl.addEventListener(
    "wheel",
    function (event) {
      event.preventDefault();
      setZoom(currentZoom() + (event.deltaY < 0 ? WHEEL_ZOOM_STEP : -WHEEL_ZOOM_STEP));
    },
    { passive: false }
  );

  canvasEl.addEventListener("keydown", function (event) {
    if (!MC.orbitByKey(spherical, event.key, KEY_ROTATE_STEP, clampPolar)) {
      return;
    }
    event.preventDefault();
    applyCamera();
    updateScaleBar();
  });

  // Picking (mappick.js), with the sector's own layers (sectorscene.js).
  var picker = createPicker(camera, canvasEl);
  sector.layers.forEach(picker.addLayer);

  function entryAtClientPoint(clientX, clientY) {
    var found = picker.pick(clientX, clientY);
    return found ? found.entry : null;
  }

  // Any button's drag turns the view from its first move; a click that
  // travelled no more than DRAG_CLICK_THRESHOLD_PX picks what's under it.
  MC.createPointerControl(canvasEl, {
    attach: true,
    dragClickPx: DRAG_CLICK_THRESHOLD_PX,
    measure: "path",
    turnAtOnce: true,
    onDragStart: function () {
      hoverEntry(null);
    },
    onDrag: function (dx, dy) {
      MC.orbitByDrag(spherical, dx, dy, ROTATE_SENSITIVITY, clampPolar);
      applyCamera();
      updateScaleBar();
    },
    clickOn: "click",
    onClick: function (event) {
      selectEntry(entryAtClientPoint(event.clientX, event.clientY));
    },
    onHover: function (event) {
      if (event.pointerType === "touch") return;
      hoverEntry(entryAtClientPoint(event.clientX, event.clientY), event);
    },
    onLeave: function () {
      hoverEntry(null);
    },
  });

  var controlsEl = document.getElementById("starmap-controls");
  if (controlsEl) {
    controlsEl.querySelectorAll("[data-action]").forEach(function (button) {
      button.addEventListener("click", function () {
        var action = button.dataset.action;
        if (action === "zoom-in") setZoom(currentZoom() + ZOOM_STEP);
        else if (action === "zoom-out") setZoom(currentZoom() - ZOOM_STEP);
        else if (action === "reset") resetView();
        else if (action === "toggle-rogue-markers") {
          setRoguesMarked(!sector.roguesMarked());
          button.setAttribute("aria-pressed", sector.roguesMarked() ? "true" : "false");
        }
      });
      if (button.dataset.action === "toggle-rogue-markers") {
        setRoguesMarked(button.getAttribute("aria-pressed") === "true");
      }
    });
  }

  // The sector page's Contents "Show on map" buttons (MAP.46) name a
  // cloud by its `key`; selecting it shows its details and its ring, and
  // brings the map into view.
  document.querySelectorAll("[data-map-target]").forEach(function (button) {
    var entry = sector.entryByKey.get(button.dataset.mapTarget);
    if (!entry) {
      button.hidden = true;
      return;
    }
    button.hidden = false;
    button.addEventListener("click", function () {
      selectEntry(entry);
      if (viewport) {
        viewport.scrollIntoView({ block: "center" });
      }
      canvasEl.focus({ preventScroll: true });
    });
  });

  // --- Scale bar -----------------------------------------------------------

  var scaleEl = document.getElementById("starmap-scale");
  var scaleBarEl = document.getElementById("starmap-scale-bar");
  var scaleLabelEl = document.getElementById("starmap-scale-label");
  var lyPerWorldUnit = data.lyPerPxAtZoom1 || 0;
  var SCALE_BAR_TARGET_PX = 70;

  // Unlike the old CSS version (a fixed 320px scene scaled by a flat CSS
  // `zoom` factor, so ly-per-pixel was that one ratio divided by `zoom`),
  // a real perspective camera's screen-pixels-per-world-unit depends on
  // both camera distance *and* the canvas's own live rendered size (this
  // panel is a responsive `min(100%, 22rem)` box, not a fixed 320px one)
  // -- so this recomputes it from first principles every time instead.
  function worldUnitsPerScreenPixel() {
    return worldUnitsPerPixel(camera, spherical.radius, canvasEl.clientHeight);
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
    scaleLabelEl.textContent = formatDistanceLy(niceLy);
  }

  // --- Resize/render loop ----------------------------------------------

  watchResize(viewport, function () {
    fitRendererToCanvas(renderer, camera, canvasEl);
    updateScaleBar();
  });

  (function animate() {
    requestAnimationFrame(animate);
    sector.update();
    updatePointHighlight();
    selectionRing.update();
    hoverRing.update();
    renderer.render(scene, camera);
  })();
}

if (canvas && sceneData) {
  initStarmap(canvas, sceneData);
}
