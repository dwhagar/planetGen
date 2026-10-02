// html/static/sectormap.js
//
// Renders the 3D sector map built by `lib/starmap.py` as a real WebGL
// scene (three.js, vendored at `static/vendor/three.module.min.js` -- see
// that directory's `THIRD_PARTY_NOTICES.txt` for why it's vendored rather
// than loaded from a CDN) instead of the CSS `transform-style:
// preserve-3d` scene this file used to drive directly. `#starmap-data`
// (a `<script type="application/json">` block `starmap.py` writes) is the
// only thing read from the page -- every position/size/color/label for
// every star and phenomenon cloud, plus the cell outline and compass
// arrow, is data `starmap.py` already computed; this file only ever turns that data
// into real 3D bodies/lines and wires up drag-to-rotate, scroll/button-
// to-zoom, and click/keyboard-for-info, the same interaction set the old
// CSS version had (a real perspective camera now does the projection/
// occlusion a browser's `preserve-3d` compositor used to).
//
// Every star, and every phenomenon that gives off light (quasar, neutron
// star, accreting black hole), is a point of light (MAP.15): a tiny bright
// core in a soft halo a fixed number of screen pixels across, the way
// `static/galaxymap3d.js` draws its bright stars, all in one `THREE.Points`
// (see "Points of light" below). A rogue planet is a faint point in the
// same set, bigger and ringed once marked (MAP.82 to MAP.84). A phenomenon
// with no light of its own (quiescent black hole, interstellar comet) or an asteroid
// field is a textured `THREE.Mesh` sphere with a fresnel glow shell from
// `./bodyRendering.js`, and a nebula or supernova remnant a see-through
// volume. Labels and the highlight ring are sprites.
//
// Built with plain DOM calls (never innerHTML/textContent-with-markup)
// when filling the info panel, same discipline the old version had --
// every star/system/phenomenon name is still database content (a system
// name can contain arbitrary characters via `--name`).

// Sibling modules are imported with this module's own `?v=<version>`
// query (html/lib/fmt.py's `static_url`), so they are cached and
// refreshed with the page's script. A plain static `import "./x.js"`
// would drop the query: an update could then leave a stale copy cached,
// and a page that also loaded the same file by its versioned URL would
// get a second, separate instance of it.
const VERSION_QUERY = new URL(import.meta.url).search;
const THREE = await import(`./vendor/three.module.min.js${VERSION_QUERY}`);
const { makeGlowMaterial } = await import(`./bodyRendering.js${VERSION_QUERY}`);
const { formatDistanceLy } = await import(`./distance.js${VERSION_QUERY}`);
const { generateButtons } = await import(`./generatebuttons.js${VERSION_QUERY}`);
const MC = await import(`./mapcontrol.js${VERSION_QUERY}`);
const {
  addField, cssVar, fitRendererToCanvas, formatAddress, isLightBackground, makeRingTexture, nearestOnScreen,
  niceScaleValue, readSceneData, watchResize, worldUnitsPerPixel,
} = await import(`./mapcore.js${VERSION_QUERY}`);

var canvas = document.getElementById("starmap-canvas");
var dataEl = document.getElementById("starmap-data");

var sceneData = readSceneData(dataEl);

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
    appendNavActions(panel, entry);
    if (!picking(entry)) {
      panel.appendChild(navLink(entry, "View phenomenon →"));
    }
    return;
  }
  addField(dl, "Star type", entry.starType);
  addField(dl, "Temperature", entry.temp);
  addField(dl, "Octant", entry.quadrant);
  addField(dl, "Location", entry.location);
  panel.appendChild(dl);
  appendNavActions(panel, entry);
  if (!picking(entry)) {
    panel.appendChild(navLink(entry, "View system →"));
  }
}

// Whether the page is choosing a NAV start or destination (NAV's "Pick
// on map", `entry.nav.pick`): the panel then offers only the pick button,
// no link that would leave the course being built (NAV.30).
function picking(entry) {
  return !!(entry.nav && entry.nav.pick);
}

// NAV links for a system or phenomenon (`entry.nav`, built by the sector
// page): in pick mode only a "Use as destination" (or start) button that
// lands on the plotted course, otherwise "Nav from here" and "Nav to
// here".
function appendNavActions(panel, entry) {
  var nav = entry.nav;
  if (!nav) {
    return;
  }
  if (nav.pick) {
    var pick = document.createElement("a");
    pick.href = nav.pick;
    pick.className = "btn starmap-pick";
    pick.textContent = nav.pickLabel;
    panel.appendChild(pick);
    return;
  }
  var links = document.createElement("p");
  links.className = "page-actions";
  [["from", "Nav from here"], ["to", "Nav to here"]].forEach(function (pair) {
    var link = document.createElement("a");
    link.href = nav[pair[0]];
    link.className = "btn btn-small";
    link.textContent = pair[1];
    links.appendChild(link);
  });
  panel.appendChild(links);
}

function showNeighborInfo(panel, entry) {
  var heading = document.createElement("h3");
  heading.textContent = entry.exists ? entry.name || "Unnamed sector" : "Not yet generated";
  panel.appendChild(heading);

  var dl = document.createElement("dl");
  addField(dl, "Address", formatAddress(entry.ringIndex, entry.layerIndex, entry.ringSlotIndex));
  addField(dl, "Designation", entry.designation);
  if (entry.brightStarCount) {
    // Bright stars the plan pre-placed here, waiting for the fill.
    var more = entry.brightStarCount - entry.brightStars.length;
    addField(dl, "Bright stars waiting",
      entry.brightStars.join("; ") + (more > 0 ? "; and " + more + " more" : ""));
  }
  panel.appendChild(dl);

  if (entry.exists) {
    panel.appendChild(navLink(entry, "View sector →"));
    return;
  }

  if (sceneData && sceneData.generate) {
    panel.appendChild(generateButtons(sceneData.generate, entry.ringIndex, entry.layerIndex, entry.ringSlotIndex,
      sceneData.edgeLy));
  }
}

// A plain `<a href>` to the entry's own page (`href`, built server-side by
// lib/starmap.py), so Back, open-in-new-tab and copy-link all work.
function navLink(entry, label) {
  var link = document.createElement("a");
  link.href = entry.href || "#";
  link.className = "btn";
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
  // A thin, bright expanding shell around a faint interior -- how a
  // supernova remnant (Cassiopeia A, the Veil) actually reads.
  supernovaRemnant: function () {
    return makeSimpleRadialTexture(128, [
      [0, "#ffb07a20"], [0.45, "#ffb07a30"], [0.62, "#ff8a5cb0"], [0.7, "#8fd6ffa0"], [0.8, "rgba(0,0,0,0)"],
    ]);
  },
  // A starless world lit only by its own internal heat: a cool violet,
  // bright enough to find on the dark scene (MAP.46) but no star color.
  roguePlanet: function () {
    return makeSimpleRadialTexture(96, [
      [0, "#e4dcff"], [0.45, "#a993f0f0"], [0.62, "#7a62c8c0"], [0.8, "rgba(0,0,0,0)"],
    ]);
  },
  // An icy nucleus inside a pale cyan coma.
  interstellarComet: function () {
    return makeSimpleRadialTexture(96, [
      [0, "#ffffff"], [0.2, "#e6fbffe0"], [0.5, "#8ff0e080"], [0.8, "rgba(0,0,0,0)"],
    ]);
  },
  // A blinding white-blue point ringed by a hot violet accretion glow --
  // brighter than anything else on the map, as a galaxy's nucleus is.
  quasar: function () {
    return makeSimpleRadialTexture(128, [
      [0, "#ffffff"], [0.18, "#f2f6ffff"], [0.4, "#a9c4ffd0"], [0.6, "#b98cff70"], [0.85, "rgba(0,0,0,0)"],
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
  supernovaRemnant: { color: "#ff8a5c", power: 2.0, strength: 1.2, scale: 1.2 },
  roguePlanet: { color: "#b59cff", power: 2.0, strength: 1.6, scale: 1.3 },
  interstellarComet: { color: "#8ff0e0", power: 1.4, strength: 1.8, scale: 1.4 },
  quasar: { color: "#c9d8ff", power: 1.0, strength: 3.0, scale: 1.7 },
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

// How much fainter a neighboring sector's cloud reaching into this one
// is drawn (lib/starmap.py marks it `neighbor`, MAP.45), so it reads as
// coming from outside rather than belonging to this sector.
var NEIGHBOR_CLOUD_DIM = 0.4;

// MAP.46: a marked rogue planet gets a ring that keeps one size on screen
// (a fraction of the canvas height) however far out the camera is, so a
// small dark world stays easy to spot zoomed out. A violet that reads on
// the light and the dark theme's map background.
var ROGUE_MARKER_COLOR = "#8f6df2";
var ROGUE_MARKER_SCREEN_SIZE = 0.03;
// MAP.84: how near (pixels) a click must land to pick a rogue planet:
// an unmarked one is a faint speck and only a click right on it takes
// it; marked, it is a big, easy target.
var ROGUE_PICK_PX = { unmarked: 3, marked: 14 };

// A real 3D body: a textured core sphere plus a fresnel glow shell
// (./bodyRendering.js), for every phenomenon that is neither a cloud
// volume nor a point of light. `coreTexture` is whatever this kind's own
// recipe already built above;
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

  var glowMaterial = makeGlowMaterial(THREE, glow.color, glow.power, glow.strength, glow.scale);
  var glowMesh = new THREE.Mesh(new THREE.SphereGeometry(1, 32, 24), glowMaterial);
  glowMesh.scale.setScalar(radius * glow.scale);

  return { core: core, glow: glowMesh };
}

// Nebulae and supernova remnants are gas, not bodies: each is drawn as a
// see-through volume whose opacity at a pixel grows with how much of the
// cloud the line of sight crosses (1 - exp(-density * path / diameter)),
// so it is densest through the middle and fades to nothing at the edge.
// The path is worked out per pixel by intersecting the view ray with the
// sphere, and only the far side of the mesh is drawn, so a cloud far
// larger than the sector still reads correctly with the camera inside it.
// A remnant is a shell: the path through its inner hollow (innerRatio of
// the radius) is left out, so it shows the bright rim the Veil does.
var CLOUD_VOLUME_VERTEX_SHADER = [
  "varying vec3 vWorld;",
  "void main() {",
  "  vec4 world = modelMatrix * vec4(position, 1.0);",
  "  vWorld = world.xyz;",
  "  gl_Position = projectionMatrix * viewMatrix * world;",
  "}",
].join("\n");

var CLOUD_VOLUME_FRAGMENT_SHADER = [
  "uniform vec3 cloudColor;",
  "uniform vec3 cloudCenter;",
  "uniform float cloudRadius;",
  "uniform float innerRatio;",
  "uniform float density;",
  "uniform float maxAlpha;",
  "varying vec3 vWorld;",
  "float pathThrough(vec3 origin, vec3 dir, float radius) {",
  "  vec3 oc = origin - cloudCenter;",
  "  float b = dot(dir, oc);",
  "  float h = b * b - (dot(oc, oc) - radius * radius);",
  "  if (h <= 0.0) return 0.0;",
  "  h = sqrt(h);",
  "  return max(-b + h, 0.0) - max(-b - h, 0.0);",
  "}",
  "void main() {",
  "  vec3 dir = normalize(vWorld - cameraPosition);",
  "  float path = pathThrough(cameraPosition, dir, cloudRadius);",
  "  if (innerRatio > 0.0) path -= pathThrough(cameraPosition, dir, cloudRadius * innerRatio);",
  "  float alpha = (1.0 - exp(-density * path / (2.0 * cloudRadius))) * maxAlpha;",
  "  if (alpha < 0.004) discard;",
  "  gl_FragColor = vec4(cloudColor, alpha);",
  "}",
].join("\n");

// Per kind: innerRatio (0 for a filled cloud), density, and the color and
// peak opacity (a nebula's come from lib/starmap.py's per-type coreColor).
var CLOUD_VOLUMES = {
  nebula: { innerRatio: 0, density: 1.6 },
  supernovaRemnant: { innerRatio: 0.82, density: 3.0, color: "#ff8a5c", alpha: 0.5 },
};

function makeCloudVolume(cloud, volume) {
  var color = volume.color || (cloud.coreColor || "#c9a8e090").slice(0, 7);
  var alpha = volume.alpha;
  if (alpha === undefined) {
    // coreColor's last two hex digits are its alpha (lib/starmap.py).
    var baked = parseInt((cloud.coreColor || "#c9a8e090").slice(7, 9), 16);
    alpha = Math.min(0.45, (isNaN(baked) ? 0x90 : baked) / 255);
  }
  if (cloud.neighbor) {
    alpha *= NEIGHBOR_CLOUD_DIM;
  }
  var material = new THREE.ShaderMaterial({
    uniforms: {
      cloudColor: { value: new THREE.Color(color) },
      cloudCenter: { value: new THREE.Vector3(cloud.x, cloud.y, cloud.z) },
      cloudRadius: { value: cloud.r },
      innerRatio: { value: volume.innerRatio },
      density: { value: volume.density },
      maxAlpha: { value: alpha },
    },
    vertexShader: CLOUD_VOLUME_VERTEX_SHADER,
    fragmentShader: CLOUD_VOLUME_FRAGMENT_SHADER,
    side: THREE.BackSide,
    transparent: true,
    depthWrite: false,
  });
  var mesh = new THREE.Mesh(new THREE.SphereGeometry(1, 48, 32), material);
  mesh.scale.setScalar(cloud.r);
  mesh.position.set(cloud.x, cloud.y, cloud.z);
  // Only its far wall is drawn, so a star inside the cloud (nearer than
  // that wall) is never tinted over and stays easy to pick out.
  return mesh;
}

// --- Points of light (MAP.15) ----------------------------------------------
//
// Every star, and every phenomenon lib/starmap.py sends a `light` for, is
// drawn the way static/galaxymap3d.js draws its bright stars: a point
// sprite a fixed number of pixels across (never sized by distance), a core
// `corePx` across in a soft halo `sizePx` across whose strength is `glow`,
// the core `bright` opaque and whitened toward its middle (by `whiten`, 1
// unless the entry says otherwise: an accreting black hole's core keeps
// its disc's orange). lib/starmap.py works every number out (from the
// star's radius, luminosity and temperature, or a fixed recipe for a
// phenomenon). The halo falls off a little more slowly than the Galaxy
// Map's, so even a red dwarf shows a bright aura on this smaller map (Boss:
// "realistic sizes with bright auras"). The only change with zoom: closing in grows a point up to
// POINT_CLOSE_GROWTH times its size. On the light theme's pale background
// a white core would vanish, so there the core is a darker shade of its
// own color instead, and the halo blends as usual.
var POINT_CLOSE_GROWTH = 1.5;
// A click within this many pixels of a point's center (or inside its
// core, for a big one) picks it.
var POINT_PICK_PX = 8;
// The highlight ring around a selected point, in pixels across at least.
var POINT_RING_MIN_PX = 20;

var POINT_VERTEX_SHADER = [
  "#include <common>",
  "#include <logdepthbuf_pars_vertex>",
  "attribute vec3 pointColor;",
  "attribute float pointSize;",
  "attribute float pointCore;",
  "attribute float pointGlow;",
  "attribute float pointBright;",
  "attribute float pointWhiten;",
  "uniform float pixelRatio;",
  "uniform float sizeScale;",
  "varying vec3 vColor;",
  "varying float vCore;",
  "varying float vGlow;",
  "varying float vBright;",
  "varying float vWhiten;",
  "void main() {",
  "  vColor = pointColor;",
  "  vCore = pointCore / pointSize;",
  "  vGlow = pointGlow;",
  "  vBright = pointBright;",
  "  vWhiten = pointWhiten;",
  "  gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);",
  "  gl_PointSize = pointSize * sizeScale * pixelRatio;",
  "  #include <logdepthbuf_vertex>",
  "}",
].join("\n");

var POINT_FRAGMENT_SHADER = [
  "#include <common>",
  "#include <logdepthbuf_pars_fragment>",
  "uniform float lightBackground;",
  "varying vec3 vColor;",
  "varying float vCore;",
  "varying float vGlow;",
  "varying float vBright;",
  "varying float vWhiten;",
  "void main() {",
  "  #include <logdepthbuf_fragment>",
  "  float r = length(gl_PointCoord * 2.0 - 1.0);",
  "  if (r > 1.0) discard;",
  "  float core = 1.0 - smoothstep(vCore * 0.5, vCore, r);",
  "  float halo = vGlow * exp(-r * r * 3.0) * (1.0 - r * r);",
  "  vec3 lit = mix(vColor, vec3(1.0), 0.6 * vBright * vWhiten);",
  "  vec3 coreColor = mix(lit, vColor * 0.55, lightBackground);",
  "  gl_FragColor = vec4(mix(vColor, coreColor, core), clamp(core * vBright + halo, 0.0, 1.0));",
  "}",
].join("\n");

// One THREE.Points holding every entry with a `light`; a neighboring
// sector's phenomenon (`neighbor`) is dimmed as its cloud would be.
function makePointsOfLight(entries, lightBackground) {
  var n = entries.length;
  var positions = new Float32Array(3 * n);
  var colors = new Float32Array(3 * n);
  var sizes = new Float32Array(n);
  var cores = new Float32Array(n);
  var glows = new Float32Array(n);
  var brights = new Float32Array(n);
  var whitens = new Float32Array(n);
  var color = new THREE.Color();
  entries.forEach(function (entry, i) {
    var light = entry.light;
    var dim = entry.neighbor ? NEIGHBOR_CLOUD_DIM : 1;
    positions.set([entry.x, entry.y, entry.z], 3 * i);
    color.set(light.color || "#ffffff");
    colors.set([color.r, color.g, color.b], 3 * i);
    sizes[i] = light.sizePx;
    cores[i] = light.corePx;
    glows[i] = light.glow * dim;
    brights[i] = light.bright * dim;
    whitens[i] = light.whiten != null ? light.whiten : 1;
  });
  var geometry = new THREE.BufferGeometry();
  geometry.setAttribute("position", new THREE.BufferAttribute(positions, 3));
  geometry.setAttribute("pointColor", new THREE.BufferAttribute(colors, 3));
  geometry.setAttribute("pointSize", new THREE.BufferAttribute(sizes, 1));
  geometry.setAttribute("pointCore", new THREE.BufferAttribute(cores, 1));
  geometry.setAttribute("pointGlow", new THREE.BufferAttribute(glows, 1));
  geometry.setAttribute("pointBright", new THREE.BufferAttribute(brights, 1));
  geometry.setAttribute("pointWhiten", new THREE.BufferAttribute(whitens, 1));
  var material = new THREE.ShaderMaterial({
    uniforms: {
      pixelRatio: { value: 1 },
      sizeScale: { value: 1 },
      lightBackground: { value: lightBackground ? 1 : 0 },
    },
    vertexShader: POINT_VERTEX_SHADER,
    fragmentShader: POINT_FRAGMENT_SHADER,
    transparent: true,
    depthWrite: false,
  });
  var points = new THREE.Points(geometry, material);
  // Drawn after the clouds, so a cloud's far wall never paints over one.
  points.renderOrder = 5;
  points.frustumCulled = false;
  return points;
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
  // The zoom policy (mapcontrol.js): a short range, zoom MIN_ZOOM to
  // MAX_ZOOM, as camera distances from the zoom-1 distance.
  var ZOOM_POLICY = MC.zoomPolicy(MC.ZOOM_RANGE, 1 / MAX_ZOOM, 1 / MIN_ZOOM);

  function clampPolar(phi) {
    return Math.max(MIN_POLAR, Math.min(MAX_POLAR, phi));
  }

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

  // The sector's own cell as a faint wireframe (starmap.py's
  // _outline_data): 12 edges, each a polyline, whose inner and outer ring
  // faces follow the ring's arc. Thin and translucent so the stars stay the
  // focus, in the muted text color so it reads in both themes.
  if (data.outline) {
    var outlineMaterial = new THREE.LineBasicMaterial({
      color: new THREE.Color(cssVar("--text-muted", "#5b6072")),
      transparent: true,
      opacity: 0.35,
      depthWrite: false,
    });
    data.outline.edges.forEach(function (edge) {
      var outlineGeometry = new THREE.BufferGeometry().setFromPoints(edge.map(function (point) {
        return new THREE.Vector3(point[0], point[1], point[2]);
      }));
      scene.add(new THREE.Line(outlineGeometry, outlineMaterial));
    });
  }

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

  // Cloud volumes (see makeCloudVolume): a click that reaches one only
  // through its far wall still picks a point of light in front of it.
  var volumeMeshes = new Set();

  // Every star and light-giving phenomenon, drawn as points of light
  // (MAP.15) and picked on screen (pointAtClientPoint), not by raycast.
  var pointEntries = (data.stars || []).concat((data.clouds || []).filter(function (cloud) {
    return cloud.light;
  }));
  var pointsOfLight = null;
  if (pointEntries.length) {
    pointsOfLight = makePointsOfLight(pointEntries, isLightBackground());
    pointsOfLight.material.uniforms.pixelRatio.value = renderer.getPixelRatio();
    scene.add(pointsOfLight);
  }

  // How much the points have grown with the camera's zoom: their own size
  // up to default zoom, POINT_CLOSE_GROWTH times it at MAX_ZOOM.
  function pointSizeScale() {
    var share = THREE.MathUtils.clamp((currentZoom() - 1) / (MAX_ZOOM - 1), 0, 1);
    return 1 + (POINT_CLOSE_GROWTH - 1) * share;
  }

  // Not in interactiveGroup: a click on a marker's empty middle should
  // still reach whatever is behind it. Rogue planets are points of light
  // (lib/starmap.py's _ROGUE_LIGHT); "Mark rogue planets" (off by
  // default) swaps in each one's markedLight and shows these rings.
  var rogueMarkers = new THREE.Group();
  rogueMarkers.visible = false;
  scene.add(rogueMarkers);
  var rogueMarkerMaterial = null;
  var roguesMarked = false;

  function setRoguesMarked(marked) {
    roguesMarked = marked;
    rogueMarkers.visible = marked;
    if (!pointsOfLight) {
      return;
    }
    var geometry = pointsOfLight.geometry;
    pointEntries.forEach(function (entry, i) {
      if (entry.kind !== "roguePlanet" || !entry.markedLight) {
        return;
      }
      var light = marked ? entry.markedLight : entry.unmarkedLight;
      var dim = entry.neighbor ? NEIGHBOR_CLOUD_DIM : 1;
      entry.light = light;
      var color = new THREE.Color(light.color || "#ffffff");
      geometry.getAttribute("pointColor").array.set([color.r, color.g, color.b], 3 * i);
      geometry.getAttribute("pointSize").array[i] = light.sizePx;
      geometry.getAttribute("pointCore").array[i] = light.corePx;
      geometry.getAttribute("pointGlow").array[i] = light.glow * dim;
      geometry.getAttribute("pointBright").array[i] = light.bright * dim;
      geometry.getAttribute("pointWhiten").array[i] = light.whiten != null ? light.whiten : 1;
    });
    ["pointColor", "pointSize", "pointCore", "pointGlow", "pointBright", "pointWhiten"].forEach(function (name) {
      geometry.getAttribute(name).needsUpdate = true;
    });
    if (highlightedPoint) {
      updatePointHighlight();
    }
  }

  (data.clouds || []).forEach(function (cloud) {
    if (cloud.kind === "roguePlanet") {
      cloud.unmarkedLight = cloud.light;
      if (!rogueMarkerMaterial) {
        rogueMarkerMaterial = new THREE.SpriteMaterial({
          map: makeRingTexture(ROGUE_MARKER_COLOR), transparent: true, opacity: 0.75, depthWrite: false,
          sizeAttenuation: false,
        });
      }
      var ring = new THREE.Sprite(rogueMarkerMaterial);
      ring.position.set(cloud.x, cloud.y, cloud.z);
      ring.scale.set(ROGUE_MARKER_SCREEN_SIZE, ROGUE_MARKER_SCREEN_SIZE, 1);
      rogueMarkers.add(ring);
    }
    if (cloud.light) {
      return;
    }
    var volume = CLOUD_VOLUMES[cloud.kind];
    if (volume) {
      var mesh = makeCloudVolume(cloud, volume);
      interactiveGroup.add(mesh);
      entryByObject.set(mesh, cloud);
      volumeMeshes.add(mesh);
      return;
    }
    var glow = glowRecipeForCloud(cloud);
    if (cloud.neighbor) {
      glow = Object.assign({}, glow, { strength: glow.strength * NEIGHBOR_CLOUD_DIM });
    }
    var bodies = makeBodySpheres(cloud.r, textureForCloud(cloud), 0xffffff, glow);
    if (cloud.neighbor) {
      bodies.core.material.opacity = NEIGHBOR_CLOUD_DIM;
    }
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

  var ringTexture = makeRingTexture(accentColor);
  var highlightSprite = new THREE.Sprite(
    new THREE.SpriteMaterial({ map: ringTexture, transparent: true, depthWrite: false })
  );
  highlightSprite.visible = false;
  scene.add(highlightSprite);
  // A point of light keeps its size on screen, so its ring does too
  // (sized each frame, updatePointHighlight).
  var pointHighlightSprite = new THREE.Sprite(
    new THREE.SpriteMaterial({ map: ringTexture, transparent: true, depthWrite: false, sizeAttenuation: false })
  );
  pointHighlightSprite.visible = false;
  pointHighlightSprite.renderOrder = 6;
  scene.add(pointHighlightSprite);
  var highlightedPoint = null;

  function highlightEntry(entry) {
    if (entry.light) {
      highlightSprite.visible = false;
      highlightedPoint = entry;
      pointHighlightSprite.position.set(entry.x, entry.y, entry.z);
      pointHighlightSprite.visible = true;
      updatePointHighlight();
      return;
    }
    highlightedPoint = null;
    pointHighlightSprite.visible = false;
    highlightSprite.position.set(entry.x, entry.y, entry.z);
    var r = Math.max(entry.r || 8, 4);
    highlightSprite.scale.set(r * 2.6, r * 2.6, 1);
    highlightSprite.visible = true;
  }

  // A sprite without size attenuation is scaled in units of the view's
  // height at distance 1, so a ring `px` pixels across needs
  // px / heightPx * 2 * tan(fov / 2).
  function updatePointHighlight() {
    if (!highlightedPoint) {
      return;
    }
    var light = highlightedPoint.light;
    var px = Math.max(POINT_RING_MIN_PX, light.corePx * 2 + 14) * pointSizeScale();
    var heightPx = canvasEl.clientHeight || 1;
    var scale = (px / heightPx) * 2 * Math.tan(THREE.MathUtils.degToRad(camera.fov) / 2);
    pointHighlightSprite.scale.set(scale, scale, 1);
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
    var point = pointAtClientPoint(clientX, clientY, rect);
    var hit = hits.length ? hits[0] : null;
    // A point of light wins unless a solid body sits in front of it (a
    // cloud volume is only ever hit at its far wall, see makeCloudVolume).
    if (point && (!hit || volumeMeshes.has(hit.object) || point.distance <= hit.distance)) {
      return point.entry;
    }
    return hit ? entryByObject.get(hit.object) || null : null;
  }

  // The point of light whose center is nearest a screen point, within
  // POINT_PICK_PX (or its own core, if bigger), as `{entry, distance}`
  // (distance from the camera, in world units), or null -- picked on
  // screen like static/galaxymap3d.js's stars, since a point is a fixed
  // number of pixels across whatever its distance.
  function pointAtClientPoint(clientX, clientY, rect) {
    var growth = pointSizeScale();
    var found = nearestOnScreen(pointEntries, camera, rect, clientX, clientY, {
      reach: function (entry) {
        return entry.kind === "roguePlanet"
          ? ROGUE_PICK_PX[roguesMarked ? "marked" : "unmarked"]
          : Math.max(POINT_PICK_PX, (entry.light.corePx / 2) * growth);
      },
    });
    if (!found) {
      return null;
    }
    var best = found.entry;
    return { entry: best, distance: camera.position.distanceTo(new THREE.Vector3(best.x, best.y, best.z)) };
  }

  // Any button's drag turns the view from its first move; a click that
  // travelled no more than DRAG_CLICK_THRESHOLD_PX picks what's under it.
  MC.createPointerControl(canvasEl, {
    attach: true,
    dragClickPx: DRAG_CLICK_THRESHOLD_PX,
    measure: "path",
    turnAtOnce: true,
    onDrag: function (dx, dy) {
      MC.orbitByDrag(spherical, dx, dy, ROTATE_SENSITIVITY, clampPolar);
      applyCamera();
      updateScaleBar();
    },
    clickOn: "click",
    onClick: function (event) {
      selectEntry(entryAtClientPoint(event.clientX, event.clientY));
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
          setRoguesMarked(!roguesMarked);
          button.setAttribute("aria-pressed", roguesMarked ? "true" : "false");
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
  var entryByKey = new Map();
  (data.clouds || []).forEach(function (cloud) {
    if (cloud.key) {
      entryByKey.set(cloud.key, cloud);
    }
  });
  document.querySelectorAll("[data-map-target]").forEach(function (button) {
    var entry = entryByKey.get(button.dataset.mapTarget);
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
    if (pointsOfLight) {
      pointsOfLight.material.uniforms.sizeScale.value = pointSizeScale();
    }
    updatePointHighlight();
    renderer.render(scene, camera);
  })();
}

if (canvas && sceneData) {
  initStarmap(canvas, sceneData);
}
