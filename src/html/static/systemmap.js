// html/static/systemmap.js
//
// Click-for-info and drill-into-moons behavior for the "System Map" panel
// built by `lib/systemmap.py`. Unlike `sectormap.js` (which has to resolve
// clicks by geometry through a rotated 3D `preserve-3d` stack -- see that
// file's own comment for why), this map is flat, static, fixed-size SVG
// with no rotation/scroll/zoom, so a plain event-target lookup is all
// clicking needs.
//
// Every planet/moon/belt/star is one `<g data-kind="..." data-*="...">` --
// clicking (or Enter/Space on a focused one) either fills the info side
// panel from its `data-*` attributes, or -- for a planet with moons,
// marked with `data-scene="planet-<id>"` -- swaps which `<svg data-scene>`
// is visible so that planet takes the star's place with its own moons
// arranged around it. Built with plain DOM calls (never innerHTML with
// unescaped content), same as `sectormap.js`, since every data-* value is
// still database content.
//
// Every star/planet/moon marker in the currently visible scene also gets
// its own live-rendered 3D sphere on `#sysmap-spheres-canvas` (three.js,
// the same vendored build `sectormap.js` uses -- see
// `static/vendor/THIRD_PARTY_NOTICES.txt`), sized and positioned to
// exactly replace that marker's own flat SVG circle -- an appearance
// layer only (color/gas-giant banding+ring/atmosphere glow, all from that
// marker's own `data-*`), not a second position plot: the SVG scenes
// remain this map's actual true-position diagram, and a sphere that fails
// to render (no WebGL) just leaves that marker's flat circle showing.

import * as THREE from "./vendor/three.module.min.js";

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

function classField(el) {
  var cls = el.dataset.class;
  if (!cls) {
    return "";
  }
  return el.dataset.classdesc ? "Class " + cls + " -- " + el.dataset.classdesc : "Class " + cls;
}

// --- Per-marker body spheres (the one 3D layer on this page) --------------

var GLOW_VERTEX_SHADER = [
  "varying vec3 vNormal;",
  "varying vec3 vViewDir;",
  "void main() {",
  "  vNormal = normalize(normalMatrix * normal);",
  "  vec4 viewPosition = modelViewMatrix * vec4(position, 1.0);",
  "  vViewDir = normalize(-viewPosition.xyz);",
  "  gl_Position = projectionMatrix * viewPosition;",
  "}",
].join("\n");

// A standard fresnel rim-glow: brightest where the surface normal points
// away from the camera (the limb), near-zero head-on -- rendered on a
// slightly larger, back-face-only, additively-blended sphere around the
// body itself, the common cheap "planet atmosphere" trick (no real
// scattering simulation, just a glow that reads as one). glowPower/
// glowStrength are the same shell reused for a star's own corona too
// (see configureBody) -- a lower glowPower spreads the glow in across
// more of the disc instead of a thin limb-only rim, and a higher
// glowStrength brightens it, together reading as "this is a light
// source" rather than the same subtle atmosphere haze a planet gets.
var GLOW_FRAGMENT_SHADER = [
  "uniform vec3 glowColor;",
  "uniform float glowPower;",
  "uniform float glowStrength;",
  "varying vec3 vNormal;",
  "varying vec3 vViewDir;",
  "void main() {",
  "  float rim = 1.0 - max(dot(vNormal, vViewDir), 0.0);",
  "  float intensity = pow(rim, glowPower) * glowStrength;",
  "  gl_FragColor = vec4(glowColor, intensity);",
  "}",
].join("\n");

function hexToRgba(hex, alpha) {
  var c = new THREE.Color(hex);
  return "rgba(" + Math.round(c.r * 255) + "," + Math.round(c.g * 255) + "," + Math.round(c.b * 255) + "," + alpha + ")";
}

// A gas giant's banding is turbulent, not a mechanical stripe grid --
// layered sine waves at different frequencies/phases give each band a
// slightly irregular width/edge, closer to a real Jupiter/Saturn photo's
// look than evenly-spaced rings would. Height-only (width 4px, repeated)
// since SphereGeometry's default UVs vary banding by latitude (v), not
// longitude (u) -- a real horizontally-banded planet has no longitude
// dependence to texture in the first place.
function makeBandTexture(baseColorHex) {
  var width = 4;
  var height = 256;
  var canvasEl = document.createElement("canvas");
  canvasEl.width = width;
  canvasEl.height = height;
  var ctx = canvasEl.getContext("2d");
  var base = new THREE.Color(baseColorHex);
  for (var y = 0; y < height; y++) {
    var t = y / height;
    var shade = 0.85 + 0.15 * Math.sin(t * Math.PI * 10) + 0.08 * Math.sin(t * Math.PI * 23 + 1.7);
    var r = Math.min(255, Math.round(base.r * 255 * shade));
    var g = Math.min(255, Math.round(base.g * 255 * shade));
    var b = Math.min(255, Math.round(base.b * 255 * shade));
    ctx.fillStyle = "rgb(" + r + "," + g + "," + b + ")";
    ctx.fillRect(0, y, width, 1);
  }
  var texture = new THREE.CanvasTexture(canvasEl);
  texture.colorSpace = THREE.SRGBColorSpace;
  return texture;
}

// A ring's own inner-to-outer alpha profile, painted along
// RingGeometry's `v` axis (0 at the inner radius, 1 at the outer) -- a
// plain vertical gradient, not a 2D radial one, since that's the axis its
// polar UV mapping actually varies along. The alternating opaque/faint
// bands are a cheap "Cassini division"-style gap cue, not a real
// particle-density simulation.
function makeRingTexture(baseColorHex) {
  var width = 4;
  var height = 128;
  var canvasEl = document.createElement("canvas");
  canvasEl.width = width;
  canvasEl.height = height;
  var ctx = canvasEl.getContext("2d");
  var gradient = ctx.createLinearGradient(0, 0, 0, height);
  gradient.addColorStop(0.0, "rgba(0,0,0,0)");
  gradient.addColorStop(0.12, hexToRgba(baseColorHex, 0.85));
  gradient.addColorStop(0.45, hexToRgba(baseColorHex, 0.3));
  gradient.addColorStop(0.6, hexToRgba(baseColorHex, 0.75));
  gradient.addColorStop(0.85, hexToRgba(baseColorHex, 0.25));
  gradient.addColorStop(1.0, "rgba(0,0,0,0)");
  ctx.fillStyle = gradient;
  ctx.fillRect(0, 0, width, height);
  var texture = new THREE.CanvasTexture(canvasEl);
  texture.colorSpace = THREE.SRGBColorSpace;
  return texture;
}

// A star's photosphere reads as mottled granulation, not a clean flat
// disc -- unlike makeBandTexture's height-only banding (correct for a
// horizontally-banded gas giant), this varies in BOTH UV directions:
// several layered sine "octaves" at different frequencies/phases/axes,
// cross-modulated against each other, give an irregular blotchy pattern
// (bright granulation cells and darker starspot-like patches) rather
// than a regular grid a plain 2D sine product would produce. Built once
// per star marker (cached as marker.starTexture, same convention
// bandTexture/ringTexture already use), not per frame.
function makeStarTexture(baseColorHex) {
  var size = 128;
  var canvasEl = document.createElement("canvas");
  canvasEl.width = canvasEl.height = size;
  var ctx = canvasEl.getContext("2d");
  var base = new THREE.Color(baseColorHex);
  var imageData = ctx.createImageData(size, size);
  for (var y = 0; y < size; y++) {
    var v = y / size;
    for (var x = 0; x < size; x++) {
      var u = x / size;
      var n =
        0.5 +
        0.22 * Math.sin(u * Math.PI * 18 + Math.sin(v * 11) * 2.0) +
        0.18 * Math.sin(v * Math.PI * 14 + Math.cos(u * 9) * 2.4) +
        0.12 * Math.sin((u + v) * Math.PI * 23) +
        0.1 * Math.sin((u - v) * Math.PI * 27 + 1.3);
      n = Math.max(0, Math.min(1, n));
      // Kept bright overall (0.72-1.12x) -- a self-luminous surface, so
      // this is texture/granulation, never allowed to read as "shadow".
      var shade = 0.72 + 0.4 * n;
      var idx = (y * size + x) * 4;
      imageData.data[idx] = Math.min(255, Math.round(base.r * 255 * shade));
      imageData.data[idx + 1] = Math.min(255, Math.round(base.g * 255 * shade));
      imageData.data[idx + 2] = Math.min(255, Math.round(base.b * 255 * shade));
      imageData.data[idx + 3] = 255;
    }
  }
  ctx.putImageData(imageData, 0, 0);
  var texture = new THREE.CanvasTexture(canvasEl);
  texture.colorSpace = THREE.SRGBColorSpace;
  return texture;
}

function parseSurfaceTempK(text) {
  var match = /(-?[\d.]+)/.exec(text || "");
  return match ? parseFloat(match[1]) : NaN;
}

// Not a spectral simulation -- just enough of a temperature cue that a
// scorched Class N reads hazy/orange, a frigid world reads pale blue, and
// an Earth-like temperate one reads sky-blue, the same "gesture at the
// physics, not model it exactly" spirit as this module's own
// `_CLASS_COLORS` (see lib/systemmap.py).
function glowColorForTemp(tempK) {
  if (!isFinite(tempK)) return "#bcdfff";
  if (tempK >= 320) return "#ffb066";
  if (tempK <= 200) return "#bcd7ff";
  return "#bfe3ff";
}

// The whole diagram's fixed coordinate space -- must match
// `lib/systemmap.py`'s own `_VIEW_SIZE_PX`, since this is what turns a
// marker's `cx`/`cy`/`r` (in that same fixed viewBox) into real canvas
// pixels below.
var VIEW_SIZE_PX = 700;

// One shared WebGL context that renders every visible marker's own sphere
// via a scissored sub-viewport per marker -- not one `<canvas>`/context
// per body. A browser caps how many WebGL contexts can exist at once
// (commonly single digits to a couple dozen, silently dropping the
// oldest once exceeded), which a system with a dozen-plus planets/moons
// would blow through immediately; a single context drawing N small
// scissored regions in one render loop has no such ceiling and is also
// just cheaper (one GL context, one set of shared geometry/materials,
// reconfigured per marker before each of that marker's own draw calls).
// Returns `null` (leaving every marker's flat circle as its plain
// fallback) if this browser can't create a WebGL context at all.
function initSphereField(canvasEl) {
  var renderer;
  try {
    renderer = new THREE.WebGLRenderer({ canvas: canvasEl, antialias: true, alpha: true });
  } catch (err) {
    return null;
  }
  renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
  renderer.setScissorTest(true);
  if (THREE.SRGBColorSpace) {
    renderer.outputColorSpace = THREE.SRGBColorSpace;
  }

  var FOV_DEG = 32;
  var scene = new THREE.Scene();
  var camera = new THREE.PerspectiveCamera(FOV_DEG, 1, 0.1, 100);
  camera.aspect = 1; // every marker's own scissored viewport is square
  camera.updateProjectionMatrix();

  // How far back the camera sits is chosen per marker (see `renderFrame`)
  // rather than fixed: a gas giant's ring reaches much further from
  // center (out to `RING_OUTER_R`) than a bare sphere (`SPHERE_R` + a
  // little for the atmosphere glow shell) does, and a fixed framing tight
  // enough for the sphere alone clips the ring clean out of the frustum --
  // confirmed directly (an early version framed for the sphere only, and
  // the ring never appeared -- it was simply outside the visible frame,
  // not a visibility/material bug).
  var elevation = 0.12; // slight downward look, a hint of "looking at a globe" rather than dead-on
  function frameCamera(halfExtent) {
    var distance = halfExtent / Math.tan(THREE.MathUtils.degToRad(FOV_DEG / 2));
    camera.position.set(0, distance * elevation, distance);
    camera.lookAt(0, 0, 0);
  }

  scene.add(new THREE.AmbientLight(0xffffff, 0.35));
  var sun = new THREE.DirectionalLight(0xffffff, 1.15);
  sun.position.set(-3, 2, 4);
  scene.add(sun);

  var bodyGroup = new THREE.Group();
  scene.add(bodyGroup);

  var SPHERE_R = 1;
  var GLOW_R = 1.14;
  var RING_INNER_R = 1.25;
  var RING_OUTER_R = 1.7;

  // A planet/moon is externally lit (its own `sun` above); a star is
  // self-luminous, so it gets its own unlit material instead -- swapped
  // onto the one shared `sphere` mesh per marker (see `configureBody`)
  // rather than a second sphere mesh, since only one is ever drawn at a
  // time regardless.
  var planetMaterial = new THREE.MeshStandardMaterial({ color: 0xffffff, roughness: 0.9, metalness: 0.05 });
  var starMaterial = new THREE.MeshBasicMaterial({ color: 0xffffff });
  var sphere = new THREE.Mesh(new THREE.SphereGeometry(SPHERE_R, 48, 32), planetMaterial);
  bodyGroup.add(sphere);

  var ring = new THREE.Mesh(
    new THREE.RingGeometry(RING_INNER_R, RING_OUTER_R, 64),
    new THREE.MeshBasicMaterial({ transparent: true, side: THREE.DoubleSide, depthWrite: false })
  );
  ring.rotation.x = THREE.MathUtils.degToRad(70);
  bodyGroup.add(ring);

  var PLANET_GLOW_POWER = 2.5;
  var PLANET_GLOW_STRENGTH = 1.0;
  var STAR_GLOW_POWER = 1.5;
  var STAR_GLOW_STRENGTH = 2.2;
  var STAR_GLOW_SCALE = 1.45;

  var glowMaterial = new THREE.ShaderMaterial({
    uniforms: {
      glowColor: { value: new THREE.Color(0xbcdfff) },
      glowPower: { value: PLANET_GLOW_POWER },
      glowStrength: { value: PLANET_GLOW_STRENGTH },
    },
    vertexShader: GLOW_VERTEX_SHADER,
    fragmentShader: GLOW_FRAGMENT_SHADER,
    side: THREE.BackSide,
    blending: THREE.AdditiveBlending,
    transparent: true,
    depthWrite: false,
  });
  var glow = new THREE.Mesh(new THREE.SphereGeometry(GLOW_R, 48, 32), glowMaterial);
  bodyGroup.add(glow);

  // The markers of whichever scene is currently visible -- see
  // `gatherSphereMarkers` -- and, in lockstep, each one's own on-canvas
  // pixel rect (`recomputeRects`), kept separate so a resize alone (no
  // marker-set change) only has to redo the cheap geometry math, not
  // rebuild every marker's cached textures.
  var markers = [];
  var rects = [];

  function disposeMarker(marker) {
    if (marker.bandTexture) {
      marker.bandTexture.dispose();
    }
    if (marker.ringTexture) {
      marker.ringTexture.dispose();
    }
    if (marker.starTexture) {
      marker.starTexture.dispose();
    }
  }

  function recomputeRects() {
    var size = canvasEl.clientWidth || 0;
    var scale = size / VIEW_SIZE_PX;
    rects = markers.map(function (marker) {
      var cx = (parseFloat(marker.circle.getAttribute("cx")) || 0) * scale;
      var cy = (parseFloat(marker.circle.getAttribute("cy")) || 0) * scale;
      var r = (parseFloat(marker.circle.getAttribute("r")) || 0) * scale;
      // The local half-extent (`frameCamera`'s own units) that a gas
      // giant's ring or a plain glow shell needs to stay in-frame -- see
      // that function's own comment -- mapped so the body's own
      // `SPHERE_R` (1 local unit) lands at exactly this marker's own
      // on-screen radius, the same relationship the flat circle it
      // replaces already had to its neighbors. A star's own glow shell
      // is scaled up by STAR_GLOW_SCALE (see configureBody) -- framed
      // for that larger reach here too, or its corona would simply clip
      // outside the frustum the same way an under-framed gas giant ring
      // once did (see this function's own earlier fix for that).
      var halfExtent = marker.isGasGiant
        ? RING_OUTER_R * 1.15
        : marker.isStar
          ? GLOW_R * STAR_GLOW_SCALE * 1.15
          : GLOW_R * 1.15;
      return { cx: cx, cy: cy, side: 2 * halfExtent * r, halfExtent: halfExtent };
    });
  }

  function setMarkers(list) {
    markers.forEach(disposeMarker);
    markers = list;
    recomputeRects();
  }

  function resize() {
    var size = canvasEl.clientWidth || 160;
    renderer.setSize(size, size, false);
    recomputeRects();
  }
  if (typeof ResizeObserver !== "undefined") {
    new ResizeObserver(resize).observe(canvasEl);
  }
  resize();
  window.addEventListener("resize", resize);

  function configureBody(marker) {
    if (marker.isStar) {
      sphere.material = starMaterial;
      // marker.starTexture (see makeStarTexture) already bakes the
      // star's own color into its granulation pattern -- tinting
      // starMaterial.color on top of that too would double-multiply it
      // (darker/oversaturated), so the material's own color resets to
      // white whenever a texture is doing the coloring instead. Falls
      // back to the old flat-tinted-sphere look only if the texture
      // somehow isn't there.
      if (marker.starTexture) {
        starMaterial.color.set(0xffffff);
        starMaterial.map = marker.starTexture;
      } else {
        starMaterial.color.set(marker.color);
        starMaterial.map = null;
      }
      starMaterial.needsUpdate = true;
      ring.visible = false;
      // A star gets its own glow shell too, tinted to its own spectral
      // color rather than `glowColorForTemp` -- a cheap "it's a light
      // source" cue, not a real corona simulation. Bigger and brighter
      // than a planet's own subtle atmosphere rim (STAR_GLOW_SCALE/
      // _POWER/_STRENGTH vs. PLANET_*) -- it needs to read unmistakably
      // as "this is the light source", not just a faint haze.
      glow.visible = true;
      glow.scale.setScalar(STAR_GLOW_SCALE);
      glowMaterial.uniforms.glowColor.value.set(marker.color);
      glowMaterial.uniforms.glowPower.value = STAR_GLOW_POWER;
      glowMaterial.uniforms.glowStrength.value = STAR_GLOW_STRENGTH;
      return;
    }
    sphere.material = planetMaterial;
    planetMaterial.map = marker.isGasGiant ? marker.bandTexture : null;
    planetMaterial.color.set(marker.isGasGiant ? "#ffffff" : marker.color);
    planetMaterial.needsUpdate = true;

    ring.visible = marker.isGasGiant;
    if (marker.isGasGiant) {
      ring.material.map = marker.ringTexture;
      ring.material.needsUpdate = true;
    }

    glow.visible = marker.hasAtmosphere;
    if (marker.hasAtmosphere) {
      glow.scale.setScalar(1);
      glowMaterial.uniforms.glowColor.value.set(marker.glowColor);
      glowMaterial.uniforms.glowPower.value = PLANET_GLOW_POWER;
      glowMaterial.uniforms.glowStrength.value = PLANET_GLOW_STRENGTH;
    }
  }

  // Renders every current marker's own sphere into its own scissored
  // sub-viewport of the shared canvas, once per frame. The whole canvas
  // is explicitly cleared first (scissor test off) rather than relying on
  // each marker's own per-rect autoClear -- the canvas persists across
  // scene switches, so without this, a marker from a now-hidden scene
  // would leave its last-drawn sphere ghosted on screen forever, since
  // nothing else would ever touch those particular pixels again.
  function renderFrame() {
    var w = canvasEl.clientWidth || 0;
    var h = canvasEl.clientHeight || 0;
    if (!w || !h) {
      return;
    }
    renderer.setScissorTest(false);
    renderer.setViewport(0, 0, w, h);
    renderer.clear();
    renderer.setScissorTest(true);

    for (var i = 0; i < markers.length; i++) {
      var marker = markers[i];
      var rect = rects[i];
      if (!rect || rect.side <= 0) {
        continue;
      }
      marker.rotation += 0.006;
      bodyGroup.rotation.y = marker.rotation;
      configureBody(marker);

      // `setViewport`/`setScissor` both take the rect's bottom-left
      // corner (WebGL's own coordinate convention), so `cy`/`side` --
      // computed above in ordinary top-left DOM pixel space, same as
      // `cx`/`cy` on the marker's own SVG circle -- need flipping here.
      var x = rect.cx - rect.side / 2;
      var glY = h - (rect.cy - rect.side / 2) - rect.side;
      renderer.setViewport(x, glY, rect.side, rect.side);
      renderer.setScissor(x, glY, rect.side, rect.side);
      frameCamera(rect.halfExtent);
      renderer.render(scene, camera);
    }
  }

  (function animate() {
    requestAnimationFrame(animate);
    renderFrame();
  })();

  return { setMarkers: setMarkers };
}

var spheresCanvas = document.getElementById("sysmap-spheres-canvas");
var sphereField = spheresCanvas ? initSphereField(spheresCanvas) : null;

// Builds `sphereField`'s marker list from whichever scene `<svg>` is now
// visible -- every `.sysmap-body` (star/planet/moon; `.sysmap-belt` isn't
// a sphere and is left alone) with its own `<circle class="sysmap-body-
// fill">`, `data-color`/`data-bodytype`/`data-hasatmosphere`/`data-
// surfacetemp` read the same way the old single-body preview read them
// off the clicked marker. A gas giant's band/ring textures are built once
// here (not per frame -- `renderFrame` only re-reads the cached texture),
// and `sysmap-sphere-active` is added so `style.css` can drop that
// marker's own flat circle fill/ring silhouette in favor of the sphere
// now standing in for them -- only once WebGL is confirmed working
// (`sphereField` non-null), so an unsupported browser keeps the plain
// flat marker instead of an empty transparent hole.
function gatherSphereMarkers(sceneEl) {
  if (!sceneEl) {
    return [];
  }
  var markers = [];
  var els = sceneEl.querySelectorAll(".sysmap-body");
  for (var i = 0; i < els.length; i++) {
    var el = els[i];
    var circle = el.querySelector("circle.sysmap-body-fill");
    if (!circle) {
      continue;
    }
    var isStar = el.dataset.kind === "star";
    var isGasGiant = !isStar && el.dataset.bodytype === "Gas Giant";
    var marker = {
      circle: circle,
      isStar: isStar,
      isGasGiant: isGasGiant,
      color: el.dataset.color || "#8a8f9c",
      hasAtmosphere: !isStar && el.dataset.hasatmosphere === "true",
      glowColor: glowColorForTemp(parseSurfaceTempK(el.dataset.surfacetemp)),
      rotation: Math.random() * Math.PI * 2, // dephased so a scene's spheres don't all spin in lockstep
    };
    if (isGasGiant) {
      marker.bandTexture = makeBandTexture(marker.color);
      marker.ringTexture = makeRingTexture(marker.color);
    }
    if (isStar) {
      marker.starTexture = makeStarTexture(marker.color);
    }
    el.classList.add("sysmap-sphere-active");
    markers.push(marker);
  }
  return markers;
}

// --- Info panel / scene switching ------------------------------------------

function showInfo(el) {
  var panel = document.getElementById("sysmap-info");
  if (!panel) {
    return;
  }
  panel.textContent = "";

  var heading = document.createElement("h3");
  heading.textContent = el.dataset.name || "Unknown";
  panel.appendChild(heading);

  var dl = document.createElement("dl");
  var kind = el.dataset.kind;
  if (kind === "star") {
    addField(dl, "Role", el.dataset.role);
    addField(dl, "Star type", el.dataset.type);
    addField(dl, "Temperature", el.dataset.temp);
    addField(dl, "Mass", el.dataset.mass);
    addField(dl, "Radius", el.dataset.radius);
    addField(dl, "Luminosity", el.dataset.lum);
  } else if (kind === "belt") {
    addField(dl, "Density", el.dataset.density);
    addField(dl, "Distance", el.dataset.distance);
    addField(dl, "Composition", el.dataset.composition);
  } else {
    addField(dl, "Class", classField(el));
    addField(dl, "Type", el.dataset.bodytype);
    addField(dl, "Zone", el.dataset.zone);
    addField(dl, "Distance", el.dataset.distance);
    addField(dl, "Period", el.dataset.period);
    addField(dl, "Gravity", el.dataset.gravity);
    addField(dl, "Atmosphere", el.dataset.atmosphere);
    addField(dl, "Surface composition", el.dataset.composition);
    addField(dl, "Surface temperature", el.dataset.surfacetemp);
    addField(dl, "Life Chemistry", el.dataset.life);
    if (kind === "moon") {
      addField(dl, "Orbits", el.dataset.parent);
    } else if (el.dataset.moons) {
      addField(dl, "Moons", el.dataset.moons);
    }
  }
  panel.appendChild(dl);

  if (el.dataset.scene) {
    var hint = document.createElement("p");
    hint.className = "hint";
    hint.textContent = "Click again to view its moon system.";
    panel.appendChild(hint);
  }
}

function resetInfo(panel) {
  panel.textContent = "";
  var hint = document.createElement("p");
  hint.className = "hint";
  hint.textContent = "Click a star, planet, moon, or asteroid belt for details.";
  panel.appendChild(hint);
}

function initSystemMap(root) {
  var info = document.getElementById("sysmap-info");
  var crumb = document.getElementById("sysmap-crumb");
  if (!info) {
    return;
  }
  // `.sysmap-orbits-layer`s are excluded here (a separate list below) --
  // `lib/systemmap.py`'s `_scene_svg_pair` now emits TWO sibling `<svg
  // class="sysmap-svg" data-scene="...">`s per scene (an orbits-only
  // layer plus this, the body-marker layer -- see that function's own
  // docstring for why), and `active`/`self`/`sceneFocusLabel` below all
  // need the body-marker one specifically (the orbits layer has no
  // `[data-self]`/`[data-kind]` markers of its own to find).
  var scenes = Array.prototype.slice.call(root.querySelectorAll(".sysmap-svg:not(.sysmap-orbits-layer)"));
  var orbitLayers = Array.prototype.slice.call(root.querySelectorAll(".sysmap-orbits-layer"));

  // A drilled-into scene's own crumb label depends on what was drilled
  // into: a planet's own moons (`kind === "planet"`, from a scene
  // reached off a planet-with-moons marker) vs. a wide binary's
  // companion star and its own planets (`kind === "star"`, from
  // `lib/systemmap.py`'s `_render_wide_binary_scenes`) -- both scenes
  // share the same "self" marker convention, just with a different noun.
  function sceneFocusLabel(sceneEl) {
    var self = sceneEl.querySelector('[data-self="true"]');
    if (!self) {
      return "";
    }
    var noun = self.dataset.kind === "star" ? "Planets of " : "Moons of ";
    return noun + self.dataset.name;
  }

  function showScene(sceneId) {
    var active = null;
    scenes.forEach(function (scene) {
      var isActive = scene.dataset.scene === sceneId;
      scene.classList.toggle("sysmap-hidden", !isActive);
      if (isActive) {
        active = scene;
      }
    });
    // Kept in lockstep with the body-marker layer above by the same
    // data-scene value, not folded into that same loop -- an orbits
    // layer is never a candidate for `active` (see `scenes`'s own
    // comment above).
    orbitLayers.forEach(function (layer) {
      layer.classList.toggle("sysmap-hidden", layer.dataset.scene !== sceneId);
    });
    if (!active) {
      return;
    }

    if (sphereField) {
      sphereField.setMarkers(gatherSphereMarkers(active));
    }

    crumb.textContent = "";
    if (sceneId !== "system") {
      var backBtn = document.createElement("button");
      backBtn.type = "button";
      backBtn.className = "starmap-btn sysmap-back-btn";
      backBtn.textContent = "← Back to system";
      backBtn.addEventListener("click", function () {
        showScene("system");
      });
      crumb.appendChild(backBtn);
      var label = sceneFocusLabel(active);
      if (label) {
        var current = document.createElement("span");
        current.className = "sysmap-crumb-current hint";
        current.textContent = label;
        crumb.appendChild(current);
      }
    }

    var self = active.querySelector('[data-self="true"]');
    if (self) {
      showInfo(self);
    } else {
      resetInfo(info);
    }
  }

  root.addEventListener("click", function (event) {
    var el = event.target.closest("[data-kind]");
    if (!el) {
      return;
    }
    if (el.dataset.scene) {
      showScene(el.dataset.scene);
    } else {
      showInfo(el);
    }
  });

  root.addEventListener("keydown", function (event) {
    if (event.key !== "Enter" && event.key !== " ") {
      return;
    }
    var el = event.target.closest("[data-kind]");
    if (!el) {
      return;
    }
    event.preventDefault();
    if (el.dataset.scene) {
      showScene(el.dataset.scene);
    } else {
      showInfo(el);
    }
  });

  showScene("system");
}

var rootEl = document.getElementById("sysmap-root");
if (rootEl) {
  initSystemMap(rootEl);
}
