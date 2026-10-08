// html/static/systemmap.js
//
// Click-for-info and drill-into-moons behavior for the "System Map" panel
// built by `planetgen/web/maps/systemmap.py`. Unlike `sectorscene.js` (which has to resolve
// clicks by geometry through a rotated 3D `preserve-3d` stack -- see that
// file's own comment for why), this map is flat, static, fixed-size SVG
// with no rotation/scroll/zoom, so a plain event-target lookup is all
// clicking needs.
//
// Every planet/moon/belt/star/facility is one `<g data-kind="..." data-*="...">` --
// clicking (or Enter/Space on a focused one) either fills the info side
// panel from its `data-*` attributes, or -- for a planet with moons,
// marked with `data-scene="planet-<id>"` -- swaps which `<svg data-scene>`
// is visible so that planet takes the star's place with its own moons
// arranged around it. Built with plain DOM calls (never innerHTML with
// unescaped content), same as `sectorscene.js`, since every data-* value is
// still database content.
//
// Every star/planet/moon marker in the currently visible scene also gets
// its own live-rendered 3D sphere on `#sysmap-spheres-canvas` (three.js,
// the same vendored build `sectorscene.js` uses -- see
// `static/vendor/THIRD_PARTY_NOTICES.txt`), sized and positioned to
// exactly replace that marker's own flat SVG circle -- an appearance
// layer only (color/gas-giant banding+ring/atmosphere glow, all from that
// marker's own `data-*`), not a second position plot: the SVG scenes
// remain this map's actual true-position diagram, and a sphere that fails
// to render (no WebGL) just leaves that marker's flat circle showing.

// Sibling modules are imported with this module's own `?v=<version>`
// query (planetgen/web/lib/fmt.py's `static_url`), so they are cached and
// refreshed with the page's script. A plain static `import "./x.js"`
// would drop the query: an update could then leave a stale copy cached,
// and a page that also loaded the same file by its versioned URL would
// get a second, separate instance of it.
const VERSION_QUERY = new URL(import.meta.url).search;
const THREE = await import(`./vendor/three.module.min.js${VERSION_QUERY}`);
const { glowInnerRatio, makeGlowMaterial, makeStarSurfaceTexture } = await import(`./bodyRendering.js${VERSION_QUERY}`);
const { formatDistanceKm: formatLadderKm } = await import(`./distance.js${VERSION_QUERY}`);
const { addField } = await import(`./mapcore.js${VERSION_QUERY}`);

function classField(el) {
  var cls = el.dataset.class;
  if (!cls) {
    return "";
  }
  return el.dataset.classdesc ? "Class " + cls + " -- " + el.dataset.classdesc : "Class " + cls;
}

// --- Per-marker body spheres (the one 3D layer on this page) --------------
//
// The glow shader (GLOW_VERTEX_SHADER/GLOW_FRAGMENT_SHADER) and the star
// granulation texture now live in ./bodyRendering.js, shared with
// sectorscene.js's own star/nebula/asteroid-field/black-hole/neutron-star
// spheres -- imported above (makeGlowMaterial/makeStarSurfaceTexture)
// rather than duplicated here.

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


// The kelvin at the front of "1,230 K (957 °C, 1,754 °F)"
// (format.format_temperature_k), comma grouping and all.
function parseSurfaceTempK(text) {
  var match = /(-?[\d.,]+)/.exec(text || "");
  return match ? parseFloat(match[1].replace(/,/g, "")) : NaN;
}

// Not a spectral simulation -- just enough of a temperature cue that a
// scorched Class N reads hazy/orange, a frigid world reads pale blue, and
// an Earth-like temperate one reads sky-blue, the same "gesture at the
// physics, not model it exactly" spirit as this module's own
// `_CLASS_COLORS` (see planetgen/web/maps/systemmap.py).
function glowColorForTemp(tempK) {
  if (!isFinite(tempK)) return "#bcdfff";
  if (tempK >= 320) return "#ffb066";
  if (tempK <= 200) return "#bcd7ff";
  return "#bfe3ff";
}

// The diagram's coordinate space when a scene has no viewBox to read --
// `planetgen/web/maps/systemmap.py`'s own `_VIEW_SIZE_PX`. A scene's viewBox is usually
// that 700 px frame, but MAP.88 widens it around the center when anything
// drawn would run past the edge, so marker positions are turned into
// canvas pixels through each marker's own scene's viewBox.
var VIEW_SIZE_PX = 700;

function viewBoxOf(el) {
  var svg = el.ownerSVGElement;
  var box = svg && svg.viewBox && svg.viewBox.baseVal;
  return box && box.width ? box : { x: 0, y: 0, width: VIEW_SIZE_PX, height: VIEW_SIZE_PX };
}

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

  var PLANET_GLOW_POWER = 1.4;
  var PLANET_GLOW_STRENGTH = 0.7;
  var STAR_GLOW_POWER = 1.8;
  var STAR_GLOW_STRENGTH = 1.0;
  var STAR_GLOW_SCALE = 1.45;
  // The glow shell's radius relative to the body sphere's own, planet
  // (unscaled) and star (scaled up by STAR_GLOW_SCALE) -- the shader
  // starts its fade at the body's limb (see bodyRendering.js).
  var PLANET_GLOW_INNER = glowInnerRatio(GLOW_R / SPHERE_R);
  var STAR_GLOW_INNER = glowInnerRatio((GLOW_R * STAR_GLOW_SCALE) / SPHERE_R);

  var glowMaterial = makeGlowMaterial(THREE, 0xbcdfff, PLANET_GLOW_POWER, PLANET_GLOW_STRENGTH, GLOW_R / SPHERE_R);
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
    rects = markers.map(function (marker) {
      var view = viewBoxOf(marker.circle);
      var scale = size / view.width;
      var cx = ((parseFloat(marker.circle.getAttribute("cx")) || 0) - view.x) * scale;
      var cy = ((parseFloat(marker.circle.getAttribute("cy")) || 0) - view.y) * scale;
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
      glowMaterial.uniforms.innerRatio.value = STAR_GLOW_INNER;
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
      glowMaterial.uniforms.innerRatio.value = PLANET_GLOW_INNER;
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
      marker.starTexture = makeStarSurfaceTexture(THREE, marker.color);
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
  } else if (kind === "facility") {
    addField(dl, "Kind", el.dataset.facilitykind);
    addField(dl, "Host", el.dataset.host);
    addField(dl, "Placement", el.dataset.placement);
    addField(dl, "Orbital distance", el.dataset.distance);
    addField(dl, "Orbital period", el.dataset.period);
    addField(dl, "Orbital speed", el.dataset.speed);
  } else if (kind === "belt") {
    addField(dl, "Density", el.dataset.density);
    addField(dl, "Distance", el.dataset.distance);
    addField(dl, "Composition", el.dataset.composition);
  } else {
    addField(dl, "Class", classField(el));
    addField(dl, "Type", el.dataset.bodytype);
    addField(dl, "Radius", el.dataset.radius);
    addField(dl, "Mass", el.dataset.mass);
    addField(dl, "Zone", el.dataset.zone);
    addField(dl, "Distance", el.dataset.distance);
    addField(dl, "Period", el.dataset.period);
    addField(dl, "Gravity", el.dataset.gravity);
    addField(dl, "Atmosphere", el.dataset.atmosphere);
    addField(dl, "Surface composition", el.dataset.composition);
    addField(dl, "Surface temperature", el.dataset.surfacetemp);
    addField(dl, "Surface pressure", el.dataset.surfacepressure);
    addField(dl, "Life Chemistry", el.dataset.life);
    if (kind === "moon") {
      addField(dl, "Orbits", el.dataset.parent);
    } else if (el.dataset.moons) {
      addField(dl, "Moons", el.dataset.moons);
    }
    addField(dl, "Note", el.dataset.note);
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
  hint.textContent = "Click a star, planet, moon, asteroid belt or facility for details.";
  panel.appendChild(hint);
}

// --- Labels that never overlap ------------------------------------------
//
// planetgen/web/maps/systemmap.py places each name with an estimated text width (the
// no-script layout). Once a scene is shown, the real text is measured
// here: a label that overlaps another label or any marker but its own is
// nudged up or down a line, and hidden if neither clears (its name is
// still the marker's aria-label and shows on hover or focus). Star names
// are placed first. A label with a leader line only stays or hides, since
// moving it would pull it off its line.

function boxesOverlap(a, b) {
  return a.x < b.x + b.w && b.x < a.x + a.w && a.y < b.y + b.h && b.y < a.y + a.h;
}

function circleHitsBox(c, box) {
  var nx = Math.max(box.x, Math.min(c.x, box.x + box.w));
  var ny = Math.max(box.y, Math.min(c.y, box.y + box.h));
  return Math.hypot(c.x - nx, c.y - ny) < c.r;
}

function layoutLabels(sceneEl) {
  var labels = Array.prototype.slice.call(sceneEl.querySelectorAll("text.sysmap-label"));
  if (!labels.length) {
    return;
  }
  labels.sort(function (a, b) {
    return (b.classList.contains("sysmap-star-label") ? 1 : 0) - (a.classList.contains("sysmap-star-label") ? 1 : 0);
  });
  var fills = Array.prototype.slice.call(sceneEl.querySelectorAll("circle.sysmap-body-fill"));
  var circles = [];
  fills.forEach(function (circle) {
    try {
      var box = circle.getBBox();
      circles.push({ x: box.x + box.width / 2, y: box.y + box.height / 2, r: box.width / 2 + 1, el: circle });
    } catch (e) {
      // Not rendered: nothing to avoid.
    }
  });
  // Labels stay inside the scene's viewBox (MAP.50): a label crossing a
  // side edge slides back in sideways, one crossing the top or bottom
  // slides back in vertically, and a vertical nudge that would leave the
  // view isn't tried.
  var view = sceneEl.viewBox && sceneEl.viewBox.baseVal;
  var margin = 2;
  function slideIn(lo, hi, min, max) {
    if (lo < min) { return min - lo; }
    if (hi > max) { return max - hi; }
    return 0;
  }
  var kept = [];
  labels.forEach(function (label) {
    label.removeAttribute("transform");
    label.classList.remove("sysmap-label-hidden");
    var box;
    try {
      box = label.getBBox();
    } catch (e) {
      return;
    }
    if (!box.width) {
      return;
    }
    var marker = label.closest(".sysmap-body");
    var hasLeader = marker && marker.querySelector(".sysmap-label-leader");
    var dx = 0, dy = 0;
    if (view && view.width) {
      dx = slideIn(box.x, box.x + box.width, view.x + margin, view.x + view.width - margin);
      dy = slideIn(box.y, box.y + box.height, view.y + margin, view.y + view.height - margin);
    }
    var shifts = hasLeader ? [0] : [0, box.height * 0.9, -box.height * 0.9, box.height * 1.8, -box.height * 1.8];
    for (var i = 0; i < shifts.length; i++) {
      var trialY = box.y + dy + shifts[i];
      if (i && view && view.height &&
          (trialY < view.y + margin || trialY + box.height > view.y + view.height - margin)) {
        continue;
      }
      var trial = { x: box.x + dx - 1, y: trialY - 1, w: box.width + 2, h: box.height + 2 };
      var clash = kept.some(function (other) { return boxesOverlap(trial, other); }) ||
        circles.some(function (c) { return !(marker && marker.contains(c.el)) && circleHitsBox(c, trial); });
      if (!clash) {
        if (dx || dy + shifts[i]) {
          label.setAttribute("transform", "translate(" + dx.toFixed(1) + " " + (dy + shifts[i]).toFixed(1) + ")");
        }
        kept.push(trial);
        return;
      }
    }
    label.classList.add("sysmap-label-hidden");
  });
}

// --- Measure distance -------------------------------------------------
//
// Real straight-line distance between any two bodies (star, planet, or
// moon) in the currently visible scene, computed from their own raw
// `data-xkm`/`data-ykm` (set by `planetgen/web/maps/systemmap.py`'s `_planet_attrs`/
// `_wide_binary_star_attrs`/etc -- see that module's own comment on why
// this has to come from real km, not this scene's drawn pixel positions:
// the shared log radial scale (`_radial_px`) preserves real angle but not
// real distance). When the straight line between the two would pass
// through the scene's own center body (its star, or -- one level in, a
// moon scene -- the planet drilled into), also computes the shortest
// path that goes around it instead: two tangent line segments from each
// point to the obstacle circle, plus the arc between the two tangent
// points -- the standard "shortest path around a circular obstacle"
// construction (the same geometry as a belt wrapped partway around a
// pulley), not a literal spline curve fit, but the real minimum distance
// a route that has to clear the body would need to cover.

// Distances go through the shared ladder (static/distance.js).
function formatDistanceKm(km) {
  if (km == null || !isFinite(km)) {
    return "unknown";
  }
  return formatLadderKm(km);
}

// The distance from point (px, py) to the nearest point on segment AB.
function pointToSegmentDistance(px, py, ax, ay, bx, by) {
  var abx = bx - ax, aby = by - ay;
  var lengthSq = abx * abx + aby * aby;
  var t = lengthSq > 0 ? ((px - ax) * abx + (py - ay) * aby) / lengthSq : 0;
  t = Math.max(0, Math.min(1, t));
  var cx = ax + t * abx, cy = ay + t * aby;
  return Math.hypot(px - cx, py - cy);
}

// How far a route keeps from each body, as a multiple of its radius: a
// star gets a real safety margin, not just its surface; a planet or moon
// only needs clearing.
var STAR_CLEARANCE_RADII = 10;
var BODY_CLEARANCE_RADII = 1.2;
// Each keep-out circle becomes a regular polygon drawn just outside it, so
// the shortest route is a shortest path through a visibility graph of
// polygon corners (Dijkstra). 72 corners put the length within 0.07% of
// the true tangent-and-arc route.
var ROUTE_POLYGON_SIDES = 72;

function markerKm(el) {
  var x = parseFloat(el.dataset.xkm), y = parseFloat(el.dataset.ykm);
  return { x: isFinite(x) ? x : 0, y: isFinite(y) ? y : 0 };
}

// The keep-out circles for a route from elA to elB in `sceneEl`: every
// star, planet and moon in the scene except the two ends. A close pair's
// two stars become one circle around both, so no route threads between
// them, unless the route starts or ends at one of them.
function routeObstacles(sceneEl, elA, elB) {
  var obstacles = [];
  var stars = [];
  var bodies = sceneEl ? sceneEl.querySelectorAll('[data-kind="star"], [data-kind="planet"], [data-kind="moon"]') : [];
  for (var i = 0; i < bodies.length; i++) {
    var el = bodies[i];
    var radius = parseFloat(el.dataset.radiuskm);
    if (!isFinite(radius) || radius <= 0 || el.dataset.xkm == null) {
      continue;
    }
    var p = markerKm(el);
    var isStar = el.dataset.kind === "star";
    var obstacle = {
      x: p.x, y: p.y, bodyR: radius,
      r: radius * (isStar ? STAR_CLEARANCE_RADII : BODY_CLEARANCE_RADII),
      el: el, name: isStar ? "the star" : (el.dataset.name || "a body"),
    };
    // A drillable companion star (a wide pair's) is far away and stands
    // on its own; the stars of a close pair (no data-scene) group up.
    if (isStar && !el.dataset.scene) {
      stars.push(obstacle);
    } else if (el !== elA && el !== elB) {
      obstacles.push(obstacle);
    }
  }
  var endIsStar = stars.some(function (o) { return o.el === elA || o.el === elB; });
  if (stars.length >= 2 && !endIsStar) {
    var cx = 0, cy = 0;
    stars.forEach(function (o) { cx += o.x / stars.length; cy += o.y / stars.length; });
    var r = 0;
    stars.forEach(function (o) { r = Math.max(r, Math.hypot(o.x - cx, o.y - cy) + o.r); });
    obstacles.push({ x: cx, y: cy, r: r, bodyR: r, el: null, name: "the stars" });
  } else {
    stars.forEach(function (o) {
      if (o.el !== elA && o.el !== elB) {
        obstacles.push(o);
      }
    });
  }
  // A keep-out circle never swallows an end: shrink it to just short of
  // the end, but never inside the body itself.
  var ends = [markerKm(elA), markerKm(elB)];
  obstacles.forEach(function (o) {
    ends.forEach(function (e) {
      var d = Math.hypot(e.x - o.x, e.y - o.y);
      if (d < o.r) {
        o.r = Math.max(o.bodyR, d * 0.98);
      }
    });
  });
  return obstacles;
}

function segmentClear(a, b, obstacles, skip) {
  for (var i = 0; i < obstacles.length; i++) {
    var o = obstacles[i];
    if (o === skip) {
      continue;
    }
    if (pointToSegmentDistance(o.x, o.y, a.x, a.y, b.x, b.y) < o.r * (1 - 1e-9)) {
      return false;
    }
  }
  return true;
}

// The shortest route from `a` to `b` (km) that stays out of every
// obstacle circle: `{points: [{x, y, obstacle?}], km}`, or `null` when no
// route exists (an end sealed inside other bodies' keep-out zones).
function shortestRoute(a, b, obstacles) {
  if (segmentClear(a, b, obstacles)) {
    return { points: [a, b], km: Math.hypot(b.x - a.x, b.y - a.y) };
  }
  var nodes = [a, b];
  var outward = 1 / Math.cos(Math.PI / ROUTE_POLYGON_SIDES) * (1 + 1e-6);
  obstacles.forEach(function (o) {
    for (var k = 0; k < ROUTE_POLYGON_SIDES; k++) {
      var angle = (2 * Math.PI * k) / ROUTE_POLYGON_SIDES;
      var node = { x: o.x + o.r * outward * Math.cos(angle), y: o.y + o.r * outward * Math.sin(angle), obstacle: o };
      var inside = obstacles.some(function (other) {
        return other !== o && Math.hypot(node.x - other.x, node.y - other.y) < other.r;
      });
      if (!inside) {
        nodes.push(node);
      }
    }
  });
  var count = nodes.length;
  var dist = new Float64Array(count).fill(Infinity);
  var prev = new Int32Array(count).fill(-1);
  var done = new Uint8Array(count);
  dist[0] = 0;
  for (;;) {
    var u = -1;
    for (var i = 0; i < count; i++) {
      if (!done[i] && dist[i] < Infinity && (u < 0 || dist[i] < dist[u])) {
        u = i;
      }
    }
    if (u < 0 || u === 1) {
      break;
    }
    done[u] = 1;
    for (var v = 0; v < count; v++) {
      if (done[v] || v === u) {
        continue;
      }
      var step = Math.hypot(nodes[v].x - nodes[u].x, nodes[v].y - nodes[u].y);
      if (dist[u] + step >= dist[v]) {
        continue;
      }
      if (segmentClear(nodes[u], nodes[v], obstacles)) {
        dist[v] = dist[u] + step;
        prev[v] = u;
      }
    }
  }
  if (dist[1] === Infinity) {
    return null;
  }
  var points = [];
  for (var at = 1; at >= 0; at = prev[at]) {
    points.unshift(nodes[at]);
  }
  return { points: points, km: dist[1] };
}

// What the route bends around, for the result panel.
function routeAvoids(route) {
  var names = [];
  route.points.forEach(function (p) {
    if (p.obstacle && names.indexOf(p.obstacle.name) < 0) {
      names.push(p.obstacle.name);
    }
  });
  return names;
}

function measurableLabel(el) {
  return el.dataset.name || "Unknown";
}

// Straight-line distance plus, when that line passes through or too near
// any body in the scene (see `routeObstacles`), the shortest route that
// keeps clear of all of them. `sceneEl` is the active body-marker
// `<svg class="sysmap-svg">`.
function computeMeasurement(elA, elB, sceneEl) {
  var a = markerKm(elA), b = markerKm(elB);
  if (!isFinite(parseFloat(elA.dataset.xkm)) || !isFinite(parseFloat(elB.dataset.xkm))) {
    return null;
  }
  var result = { straightKm: Math.hypot(b.x - a.x, b.y - a.y), routeKm: null, routeLabel: null, route: null };
  var route = shortestRoute(a, b, routeObstacles(sceneEl, elA, elB));
  result.route = route || { points: [a, b], km: result.straightKm };
  if (route && route.points.length > 2) {
    result.routeKm = route.km;
    result.routeLabel = "Route clear of " + routeAvoids(route).join(", ");
  } else if (!route) {
    result.routeLabel = "No clear route";
  }
  return result;
}

// --- Drawing the measured route ----------------------------------------
//
// The markers sit on a log radial scale around the scene's center
// (planetgen/web/maps/systemmap.py `_radial_px`, whose bounds the scene carries as
// data-lokm/data-hikm), so a straight line in km is a curve on the map:
// each leg is sampled and every sample mapped the same way the markers
// were. The two ends snap to their markers' drawn centers (markers can be
// nudged apart for legibility), and a bend around a body keeps outside
// that body's drawn marker.

var SVG_NS = "http://www.w3.org/2000/svg";

function sceneScale(sceneEl) {
  var lo = parseFloat(sceneEl.dataset.lokm), hi = parseFloat(sceneEl.dataset.hikm);
  return {
    lo: lo, hi: hi, c: parseFloat(sceneEl.dataset.cpx),
    min: parseFloat(sceneEl.dataset.minpx), spread: parseFloat(sceneEl.dataset.spreadpx),
    ok: isFinite(lo) && isFinite(hi) && lo > 0,
  };
}

function kmToPx(scale, x, y) {
  var d = Math.hypot(x, y);
  if (d <= 0) {
    return { x: scale.c, y: scale.c };
  }
  var r;
  if (scale.hi <= scale.lo) {
    r = scale.min + scale.spread;
  } else {
    var clamped = Math.max(scale.lo, Math.min(scale.hi, d));
    r = scale.min + (Math.log10(clamped) - Math.log10(scale.lo)) / (Math.log10(scale.hi) - Math.log10(scale.lo)) * scale.spread;
  }
  return { x: scale.c + r * x / d, y: scale.c - r * y / d };
}

function markerCenterPx(el) {
  var shape = el.querySelector("circle") || el;
  try {
    var box = shape.getBBox();
    return { x: box.x + box.width / 2, y: box.y + box.height / 2, r: Math.max(box.width, box.height) / 2 };
  } catch (e) {
    return null;
  }
}

function routePathPx(sceneEl, elA, elB, route) {
  var scale = sceneScale(sceneEl);
  var startPx = markerCenterPx(elA), endPx = markerCenterPx(elB);
  if (!scale.ok || !startPx || !endPx) {
    return null;
  }
  var obstaclePx = new Map();
  function bendPx(p) {
    var o = p.obstacle;
    var mapped = kmToPx(scale, p.x, p.y);
    if (!o) {
      return mapped;
    }
    if (!obstaclePx.has(o)) {
      obstaclePx.set(o, o.el ? markerCenterPx(o.el) : kmToPx(scale, o.x, o.y));
    }
    var center = obstaclePx.get(o);
    if (!center) {
      return mapped;
    }
    // Keep the bend outside the body's drawn marker, in the real direction.
    var angle = Math.atan2(-(p.y - o.y), p.x - o.x);
    var reach = Math.max(Math.hypot(mapped.x - center.x, mapped.y - center.y), (center.r || 0) + 6);
    return { x: center.x + reach * Math.cos(angle), y: center.y + reach * Math.sin(angle) };
  }
  // Every drawn marker the route must not be seen crossing: the log scale
  // can put a line that clears a body in km right over its (much larger)
  // drawn marker, so samples inside one are pushed out to its edge.
  var glyphs = [];
  var bodies = sceneEl.querySelectorAll('[data-kind="star"], [data-kind="planet"], [data-kind="moon"]');
  for (var g = 0; g < bodies.length; g++) {
    if (bodies[g] !== elA && bodies[g] !== elB) {
      var glyph = markerCenterPx(bodies[g]);
      if (glyph) {
        glyphs.push(glyph);
      }
    }
  }
  function clearOfGlyphs(p) {
    for (var k = 0; k < glyphs.length; k++) {
      var c = glyphs[k], need = c.r + 5;
      var dx = p.x - c.x, dy = p.y - c.y, d = Math.hypot(dx, dy);
      if (d < need) {
        if (d < 1e-6) {
          dx = 1; dy = 0; d = 1;
        }
        p = { x: c.x + dx / d * need, y: c.y + dy / d * need };
      }
    }
    return p;
  }
  var points = route.points;
  var out = [startPx];
  for (var i = 0; i < points.length - 1; i++) {
    var p = points[i], q = points[i + 1];
    var aroundOne = p.obstacle && p.obstacle === q.obstacle;
    var steps = aroundOne ? 1 : 24;
    for (var s = 1; s <= steps; s++) {
      if (i === points.length - 2 && s === steps) {
        out.push(endPx);
      } else if (aroundOne || s === steps) {
        out.push(bendPx(s === steps ? q : p));
      } else {
        var t = s / steps;
        out.push(kmToPx(scale, p.x + (q.x - p.x) * t, p.y + (q.y - p.y) * t));
      }
    }
  }
  for (var j = 1; j < out.length - 1; j++) {
    out[j] = clearOfGlyphs(out[j]);
  }
  return out;
}

function clearMeasurePath(root) {
  var old = root.querySelectorAll(".sysmap-measure-path");
  for (var i = 0; i < old.length; i++) {
    old[i].remove();
  }
}

function drawMeasurePath(root, sceneEl, elA, elB, measurement) {
  clearMeasurePath(root);
  if (!measurement || !measurement.route) {
    return;
  }
  var pts = routePathPx(sceneEl, elA, elB, measurement.route);
  if (!pts) {
    return;
  }
  var line = document.createElementNS(SVG_NS, "polyline");
  line.setAttribute("class", "sysmap-measure-path");
  line.setAttribute("points", pts.map(function (p) { return p.x.toFixed(1) + "," + p.y.toFixed(1); }).join(" "));
  // In the scene's orbits layer, under the sphere canvas, so the ends
  // tuck under the two bodies' spheres (the layer is aria-hidden).
  var layer = root.querySelector('.sysmap-orbits-layer[data-scene="' + sceneEl.dataset.scene + '"]') || sceneEl;
  layer.appendChild(line);
}

function showMeasurementPrompt(panel, el) {
  panel.textContent = "";
  var heading = document.createElement("h3");
  heading.textContent = measurableLabel(el) + " selected";
  panel.appendChild(heading);
  var hint = document.createElement("p");
  hint.className = "hint";
  hint.textContent = "Click a second star, planet, or moon in this same view to measure the distance.";
  panel.appendChild(hint);
}

function showMeasurementResult(panel, elA, elB, measurement) {
  panel.textContent = "";
  var heading = document.createElement("h3");
  heading.textContent = "Distance";
  panel.appendChild(heading);

  var dl = document.createElement("dl");
  addField(dl, "Between", measurableLabel(elA) + " and " + measurableLabel(elB));
  addField(dl, "Straight-line", measurement ? formatDistanceKm(measurement.straightKm) : "unavailable");
  if (measurement && measurement.routeKm != null) {
    var detourKm = measurement.routeKm - measurement.straightKm;
    addField(dl, measurement.routeLabel, formatDistanceKm(measurement.routeKm) +
      (detourKm >= 1 ? " (" + formatDistanceKm(detourKm) + " longer)" : ""));
  } else if (measurement && measurement.routeLabel) {
    addField(dl, measurement.routeLabel, "every way passes too near a body");
  }
  panel.appendChild(dl);

  var hint = document.createElement("p");
  hint.className = "hint";
  hint.textContent = "Click a body to start a new measurement, or turn off Measure distance.";
  panel.appendChild(hint);
}

function initSystemMap(root) {
  var info = document.getElementById("sysmap-info");
  var crumb = document.getElementById("sysmap-crumb");
  if (!info) {
    return;
  }
  // `.sysmap-orbits-layer`s are excluded here (a separate list below) --
  // `planetgen/web/maps/systemmap.py`'s `_scene_svg_pair` now emits TWO sibling `<svg
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
  // `planetgen/web/maps/systemmap.py`'s `_render_wide_binary_scenes`) -- both scenes
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

    clearMeasureSelection();
    layoutLabels(active);

    var self = active.querySelector('[data-self="true"]');
    if (self) {
      showInfo(self);
    } else {
      resetInfo(info);
    }
  }

  // --- Measure distance ---------------------------------------------
  //
  // Off by default -- normal clicks keep their existing behavior (info
  // panel / drill into a planet's moons) until this is toggled on, at
  // which point clicking a star/planet/moon selects it as one of two
  // measurement endpoints instead (see `computeMeasurement` above).
  var measureBtn = document.getElementById("sysmap-measure-btn");
  var measureMode = false;
  var measureSelection = [];

  function clearMeasureSelection() {
    measureSelection.forEach(function (el) {
      el.classList.remove("sysmap-measure-selected");
    });
    measureSelection = [];
    clearMeasurePath(root);
  }

  function setMeasureMode(on) {
    measureMode = on;
    clearMeasureSelection();
    if (measureBtn) {
      measureBtn.setAttribute("aria-pressed", on ? "true" : "false");
      measureBtn.classList.toggle("starmap-btn-active", on);
    }
    if (on) {
      info.textContent = "";
      var hint = document.createElement("p");
      hint.className = "hint";
      hint.textContent = "Click two stars, planets, or moons in the same view to measure the distance between them.";
      info.appendChild(hint);
    } else {
      resetInfo(info);
    }
  }

  function isMeasurable(el) {
    var kind = el.dataset.kind;
    return (kind === "star" || kind === "planet" || kind === "moon") && el.dataset.xkm != null;
  }

  function handleMeasureClick(el) {
    if (measureSelection.length >= 2) {
      clearMeasureSelection();
    }
    if (measureSelection.indexOf(el) !== -1) {
      return;
    }
    measureSelection.push(el);
    el.classList.add("sysmap-measure-selected");
    if (measureSelection.length < 2) {
      showMeasurementPrompt(info, el);
      return;
    }
    var sceneEl = el.closest(".sysmap-svg");
    var measurement = computeMeasurement(measureSelection[0], measureSelection[1], sceneEl);
    showMeasurementResult(info, measureSelection[0], measureSelection[1], measurement);
    drawMeasurePath(root, sceneEl, measureSelection[0], measureSelection[1], measurement);
  }

  if (measureBtn) {
    measureBtn.addEventListener("click", function () {
      setMeasureMode(!measureMode);
    });
  }

  root.addEventListener("click", function (event) {
    var el = event.target.closest("[data-kind]");
    if (!el) {
      return;
    }
    if (measureMode) {
      if (isMeasurable(el)) {
        handleMeasureClick(el);
      }
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
    if (measureMode) {
      if (isMeasurable(el)) {
        handleMeasureClick(el);
      }
      return;
    }
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
