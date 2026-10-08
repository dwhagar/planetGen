// html/static/nebulaview.js
//
// A nebula page's 3D view (MAP.105): the nebula's own shape (the mesh
// `GET /galaxy/nebula/<id>/shape?lod=full` serves, drawn translucent)
// with the brightest stars of the galaxy round it, dimmed, there only for
// reference (`GET /galaxy/nebula/<id>/surroundings`). One parsec is one
// scene unit. Drag turns the view, the wheel (or pinch) zooms, and it
// turns slowly by itself until the visitor touches it, unless they prefer
// reduced motion. The panel's still drawing stays until this swaps in the
// canvas, and stays for good without WebGL or when a fetch fails.
// (planetgen/web/maps/phenomenonrender.py's `render_nebula_view_panel`.)

const VERSION_QUERY = new URL(import.meta.url).search;
const THREE = await import(`./vendor/three.module.min.js${VERSION_QUERY}`);
const { buildNebulaMesh } = await import(`./nebulamesh.js${VERSION_QUERY}`);

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

function start(view, shape, around) {
  var renderer;
  try {
    renderer = new THREE.WebGLRenderer({ canvas: canvas, antialias: true });
  } catch (err) {
    return; // No WebGL: the still drawing stays.
  }
  var radius = around.radius_pc;
  var scene = new THREE.Scene();
  scene.background = new THREE.Color(0x05070c);

  var pivot = new THREE.Group();
  scene.add(pivot);
  var mesh = buildNebulaMesh(THREE, shape, radius, [view.color, view.coreOpacity, view.edgeOpacity],
    { depthTest: true });
  // The surface alone is faint from outside: a brighter core helps read the shape.
  mesh.material.opacity = Math.max(mesh.material.opacity, 0.22);
  pivot.add(mesh);
  pivot.add(starPoints(around.stars || []));

  var camera = new THREE.PerspectiveCamera(40, 1, radius * 0.01, around.half_width_pc * 12);
  var zoom = 1;
  function placeCamera() {
    var distance = radius * 3.2 * zoom;
    camera.position.set(0, Math.sin(TILT) * distance, Math.cos(TILT) * distance);
    camera.lookAt(0, 0, 0);
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
    renderer.render(scene, camera);
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
    pivot.rotation.y += (event.clientX - dragging.x) * TURN_PER_PX;
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
        pivot.rotation.y += AUTO_TURN_PER_S * dt;
      }
      draw();
    }
  });
}

var view = viewport && canvas ? readView() : null;
if (view) {
  Promise.all([getJson(view.shapePath + "?lod=full"), getJson(view.surroundingsPath)]).then(function (results) {
    start(view, results[0], results[1]);
  }, function () {});
}
