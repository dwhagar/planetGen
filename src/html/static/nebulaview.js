// html/static/nebulaview.js
//
// A nebula page's 3D view (MAP.105, MAP.142): the nebula drawn from its own cloud of points (`GET /galaxy/nebula/<id>/volume`, the
// same points and the same layer, nebulavolumelayer.js, that draw it on the Galaxy Map and the Sector Map: soft splats that fade
// unevenly toward an irregular edge) with the brightest stars of the galaxy round it, dimmed, there only for reference
// (`GET /galaxy/nebula/<id>/surroundings`, which also says where the nebula is). The scene is in the galaxy's own frame, in parsecs,
// and the camera turns round the nebula's centre. Drag turns the view, the wheel (or pinch) zooms, and it turns slowly by itself
// until the visitor touches it, unless they prefer reduced motion. The panel's still drawing stays until this swaps in the canvas,
// and stays for good without WebGL (or a float render target, which the layer needs), for a nebula never placed in the galaxy, or when
// a fetch fails.
// (planetgen/web/maps/phenomenonrender.py's `render_nebula_view_panel`.)

const VERSION_QUERY = new URL(import.meta.url).search;
const THREE = await import(`./vendor/three.module.min.js${VERSION_QUERY}`);
const { createNebulaVolumeLayer } = await import(`./nebulavolumelayer.js${VERSION_QUERY}`);
const { singleSource } = await import(`./nebulavolume.js${VERSION_QUERY}`);
const { FAMILIES, familyColour } = await import(`./galaxynebulae.js${VERSION_QUERY}`);
const { NEBULA_LOOKS } = await import(`./nebulalook.js${VERSION_QUERY}`);

var viewport = document.getElementById("nebulaview");
var canvas = document.getElementById("nebulaview-canvas");

var reducedMotion = typeof window.matchMedia === "function"
  && window.matchMedia("(prefers-reduced-motion: reduce)").matches;

// How dim the reference stars are against the nebula.
var STAR_OPACITY = 0.5;
var STAR_PX = 2.2;
var TURN_PER_PX = 0.01;
var AUTO_TURN_PER_S = 0.12;
var TILT = 0.45;
var FOV_DEG = 40;
var ZOOM_RANGE = [0.35, 6];

function readView() {
  try {
    return JSON.parse(viewport.getAttribute("data-view"));
  } catch (err) {
    return null;
  }
}

function getJson(url) {
  return fetch(url, { headers: { Accept: "application/json" } }).then(function (response) {
    if (!response.ok) {
      throw new Error(url + " " + response.status);
    }
    return response.json();
  });
}

// A star's color from its temperature: orange below 4000 K to blue-white
// above 9000 K, kept pale because the points are dimmed anyway.
function starColor(temperatureK) {
  var t = Math.max(0, Math.min(1, ((temperatureK || 5800) - 3000) / 9000));
  var color = new THREE.Color();
  color.setHSL(0.08 + 0.5 * t, 0.55, 0.78);
  return color;
}

function starPoints(stars) {
  var positions = new Float32Array(stars.length * 3);
  var colors = new Float32Array(stars.length * 3);
  stars.forEach(function (star, i) {
    positions[3 * i] = star.x;
    positions[3 * i + 1] = star.y;
    positions[3 * i + 2] = star.z;
    var color = starColor(star.temperature_k);
    colors[3 * i] = color.r;
    colors[3 * i + 1] = color.g;
    colors[3 * i + 2] = color.b;
  });
  var geometry = new THREE.BufferGeometry();
  geometry.setAttribute("position", new THREE.BufferAttribute(positions, 3));
  geometry.setAttribute("color", new THREE.BufferAttribute(colors, 3));
  return new THREE.Points(geometry, new THREE.PointsMaterial({
    size: STAR_PX, sizeAttenuation: false, vertexColors: true, transparent: true, opacity: STAR_OPACITY,
    depthWrite: false,
  }));
}

function start(view, around) {
  var center = around.center_pc;
  if (!center) {
    return; // Never placed in the galaxy: no points to draw, the still drawing stays.
  }
  var renderer;
  try {
    renderer = new THREE.WebGLRenderer({ canvas: canvas, antialias: true });
  } catch (err) {
    return; // No WebGL: the still drawing stays.
  }
  var radius = around.radius_pc;
  var scene = new THREE.Scene();
  scene.background = new THREE.Color(0x000000);

  var camera = new THREE.PerspectiveCamera(FOV_DEG, 1, radius * 0.01, around.half_width_pc * 12);
  camera.up.set(0, 0, 1);   // the galaxy's own frame: z is up
  var layer = createNebulaVolumeLayer(scene, camera, canvas, renderer, {
    colours: FAMILIES.map(function (family) { return familyColour(family, NEBULA_LOOKS); }),
    fetchBytes: function () {
      return fetch(view.volumePath).then(function (response) {
        if (!response.ok) {
          throw new Error(view.volumePath + " " + response.status);
        }
        return response.arrayBuffer();
      });
    },
    onChange: function () { draw(); },
  });
  if (!layer) {
    renderer.dispose();
    return; // No float render target: the still drawing stays.
  }
  var source = singleSource(view.nebulaId, { x: center[0], y: center[1], z: center[2], radiusPc: radius }, 1);

  // The reference stars sit at their offsets from the nebula's centre, the group at the centre.
  var stars = starPoints(around.stars || []);
  stars.position.set(center[0], center[1], center[2]);
  scene.add(stars);

  var azimuth = 0;
  var zoom = 1;
  function placeCamera() {
    var distance = radius * 3.2 * zoom;
    camera.position.set(
      center[0] + Math.cos(TILT) * Math.cos(azimuth) * distance,
      center[1] + Math.cos(TILT) * Math.sin(azimuth) * distance,
      center[2] + Math.sin(TILT) * distance);
    camera.lookAt(center[0], center[1], center[2]);
    camera.updateMatrixWorld();
  }
  placeCamera();

  canvas.hidden = false;
  viewport.classList.add("phenomrender-live");

  function resize() {
    var w = viewport.clientWidth;
    var h = viewport.clientHeight;
    if (!w || !h) {
      return;
    }
    renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
    renderer.setSize(w, h, false);
    camera.aspect = w / h;
    camera.updateProjectionMatrix();
  }

  function draw() {
    var focalPx = (canvas.clientHeight || 1) / 2 / Math.tan(THREE.MathUtils.degToRad(FOV_DEG) / 2);
    layer.update([source], null, focalPx);
    layer.prepare();
    renderer.render(scene, camera);
    canvas.dataset.volumeDrawn = String(layer.drawn());   // the tests read this: how many volumes are up
  }

  var touched = reducedMotion;
  var dragging = null;
  canvas.addEventListener("pointerdown", function (event) {
    dragging = { x: event.clientX };
    touched = true;
    canvas.setPointerCapture(event.pointerId);
  });
  canvas.addEventListener("pointermove", function (event) {
    if (!dragging) {
      return;
    }
    azimuth -= (event.clientX - dragging.x) * TURN_PER_PX;
    placeCamera();
    dragging = { x: event.clientX };
    draw();
  });
  function endDrag() {
    dragging = null;
  }
  canvas.addEventListener("pointerup", endDrag);
  canvas.addEventListener("pointercancel", endDrag);
  canvas.addEventListener("wheel", function (event) {
    event.preventDefault();
    touched = true;
    zoom = Math.max(ZOOM_RANGE[0], Math.min(ZOOM_RANGE[1], zoom * Math.exp(event.deltaY * 0.001)));
    placeCamera();
    draw();
  }, { passive: false });

  // MAP.138: C (or the button) puts the view back about the nebula's center
  // at its opening angle and zoom.
  function recenter() {
    azimuth = 0;
    zoom = 1;
    touched = true;
    placeCamera();
    draw();
  }
  canvas.tabIndex = 0;
  canvas.addEventListener("keydown", function (event) {
    if ((event.key === "c" || event.key === "C") && !event.ctrlKey && !event.metaKey && !event.altKey) {
      event.preventDefault();
      recenter();
    }
  });
  var recenterButton = document.getElementById("nebulaview-recenter");
  if (recenterButton) {
    recenterButton.hidden = false;
    recenterButton.addEventListener("click", recenter);
  }

  if (typeof ResizeObserver === "function") {
    new ResizeObserver(function () {
      resize();
      draw();
    }).observe(viewport);
  }
  resize();
  draw();
  // The tests read this: the view is up, with this many reference stars.
  canvas.dataset.nebulaView = "ready";
  canvas.dataset.stars = String((around.stars || []).length);
  canvas.__nebulaVolumes = layer;     // the browser tests read what is drawn
  if (reducedMotion) {
    return;
  }

  var visible = true;
  if (typeof IntersectionObserver === "function") {
    new IntersectionObserver(function (entries) {
      visible = entries[0].isIntersecting;
    }).observe(viewport);
  }
  var last = performance.now();
  renderer.setAnimationLoop(function (now) {
    var dt = Math.min(0.1, (now - last) / 1000);
    last = now;
    if (visible && !document.hidden) {
      if (!touched) {
        azimuth -= AUTO_TURN_PER_S * dt;
        placeCamera();
      }
      draw();
    }
  });
}

var view = viewport && canvas ? readView() : null;
if (view) {
  getJson(view.surroundingsPath).then(function (around) {
    start(view, around);
  }, function () {});
}
