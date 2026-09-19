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
// A planet/moon selection also redraws `#sysmap-preview` -- a small
// rotating shaded sphere (three.js, the same vendored build `sectormap.js`
// uses -- see `static/vendor/THIRD_PARTY_NOTICES.txt`), the one genuinely
// 3D element on this otherwise flat-SVG page. It's an appearance preview
// only (color/gas-giant banding+ring/atmosphere glow from the clicked
// marker's own `data-*`), not a second position plot -- the SVG scenes
// above remain this map's actual true-position diagram.

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

// --- Body preview (the one 3D element on this page) -----------------------

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
// scattering simulation, just a glow that reads as one).
var GLOW_FRAGMENT_SHADER = [
  "uniform vec3 glowColor;",
  "varying vec3 vNormal;",
  "varying vec3 vViewDir;",
  "void main() {",
  "  float rim = 1.0 - max(dot(vNormal, vViewDir), 0.0);",
  "  float intensity = pow(rim, 2.5);",
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

function initPreview(canvasEl) {
  var renderer = new THREE.WebGLRenderer({ canvas: canvasEl, antialias: true, alpha: true });
  renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
  if (THREE.SRGBColorSpace) {
    renderer.outputColorSpace = THREE.SRGBColorSpace;
  }

  var FOV_DEG = 32;
  var scene = new THREE.Scene();
  var camera = new THREE.PerspectiveCamera(FOV_DEG, 1, 0.1, 100);

  // How far back the camera sits is chosen per body (see `frameCamera`,
  // called from `updatePreview`) rather than fixed: a gas giant's ring
  // reaches much further from center (out to `RING_OUTER_R`) than a bare
  // sphere (`SPHERE_R` + a little for the atmosphere glow shell) does, and
  // a fixed framing tight enough for the sphere alone clips the ring
  // clean out of the frustum -- confirmed directly (an early version framed
  // for the sphere only, and the ring never appeared -- it was simply
  // outside the visible frame, not a visibility/material bug).
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

  var sphere = new THREE.Mesh(
    new THREE.SphereGeometry(SPHERE_R, 48, 32),
    new THREE.MeshStandardMaterial({ color: 0xffffff, roughness: 0.9, metalness: 0.05 })
  );
  bodyGroup.add(sphere);

  var ring = new THREE.Mesh(
    new THREE.RingGeometry(RING_INNER_R, RING_OUTER_R, 64),
    new THREE.MeshBasicMaterial({ transparent: true, side: THREE.DoubleSide, depthWrite: false })
  );
  ring.rotation.x = THREE.MathUtils.degToRad(70);
  ring.visible = false;
  bodyGroup.add(ring);

  var glowMaterial = new THREE.ShaderMaterial({
    uniforms: { glowColor: { value: new THREE.Color(0xbcdfff) } },
    vertexShader: GLOW_VERTEX_SHADER,
    fragmentShader: GLOW_FRAGMENT_SHADER,
    side: THREE.BackSide,
    blending: THREE.AdditiveBlending,
    transparent: true,
    depthWrite: false,
  });
  var glow = new THREE.Mesh(new THREE.SphereGeometry(GLOW_R, 48, 32), glowMaterial);
  glow.visible = false;
  bodyGroup.add(glow);

  frameCamera(GLOW_R * 1.15); // a sensible default before anything's been clicked yet

  function resize() {
    var size = canvasEl.clientWidth || 160;
    renderer.setSize(size, size, false);
    camera.aspect = 1;
    camera.updateProjectionMatrix();
  }
  if (typeof ResizeObserver !== "undefined") {
    new ResizeObserver(resize).observe(canvasEl);
  }
  resize();
  window.addEventListener("resize", resize);

  (function animate() {
    requestAnimationFrame(animate);
    bodyGroup.rotation.y += 0.006;
    renderer.render(scene, camera);
  })();

  return {
    sphere: sphere, ring: ring, glow: glow, glowMaterial: glowMaterial,
    frameCamera: frameCamera, ringOuterR: RING_OUTER_R, glowR: GLOW_R,
  };
}

var previewCanvas = document.getElementById("sysmap-preview-canvas");
var previewContainer = document.getElementById("sysmap-preview");
var preview = previewCanvas ? initPreview(previewCanvas) : null;

function updatePreview(dataset) {
  if (!preview || !previewContainer) {
    return;
  }
  var color = dataset.color || "#8a8f9c";
  var isGasGiant = dataset.bodytype === "Gas Giant";

  // A gas giant's ring reaches out to `ringOuterR` in its own local X --
  // untouched by the ring's own tilt (only its Y-extent foreshortens, see
  // this file's own `frameCamera` comment) -- so it, not the sphere/glow,
  // is what the camera needs to fit back far enough for.
  preview.frameCamera((isGasGiant ? preview.ringOuterR : preview.glowR) * 1.15);

  if (preview.sphere.material.map) {
    preview.sphere.material.map.dispose();
  }
  preview.sphere.material.map = isGasGiant ? makeBandTexture(color) : null;
  preview.sphere.material.color.set(isGasGiant ? "#ffffff" : color);
  preview.sphere.material.needsUpdate = true;

  preview.ring.visible = isGasGiant;
  if (isGasGiant) {
    if (preview.ring.material.map) {
      preview.ring.material.map.dispose();
    }
    preview.ring.material.map = makeRingTexture(color);
    preview.ring.material.needsUpdate = true;
  }

  var hasAtmosphere = dataset.hasatmosphere === "true";
  preview.glow.visible = hasAtmosphere;
  if (hasAtmosphere) {
    preview.glowMaterial.uniforms.glowColor.value.set(glowColorForTemp(parseSurfaceTempK(dataset.surfacetemp)));
  }

  previewContainer.hidden = false;
}

function hidePreview() {
  if (previewContainer) {
    previewContainer.hidden = true;
  }
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

  if (kind === "planet" || kind === "moon") {
    updatePreview(el.dataset);
  } else {
    hidePreview();
  }
}

function resetInfo(panel) {
  panel.textContent = "";
  var hint = document.createElement("p");
  hint.className = "hint";
  hint.textContent = "Click a star, planet, moon, or asteroid belt for details.";
  panel.appendChild(hint);
  hidePreview();
}

function initSystemMap(root) {
  var info = document.getElementById("sysmap-info");
  var crumb = document.getElementById("sysmap-crumb");
  if (!info) {
    return;
  }
  var scenes = Array.prototype.slice.call(root.querySelectorAll(".sysmap-svg"));

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
    if (!active) {
      return;
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
