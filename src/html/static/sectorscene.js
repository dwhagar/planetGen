// The Sector Map's scene, as a part of any map (MAP.66, MAP.68): one
// sector's stars, phenomena, neighboring-sector markers, cell outline and
// compass, built from planetgen/web/maps/starmap.py's JSON into a
// THREE.Group, with the picker layers (mappick.js), info-panel contents
// and tooltip text for what is in it. The Sector Map draws it on its own
// page; the Galaxy Map opens it as the drill-down's last stage, in its
// place in the galaxy.
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
// volume. Labels are sprites.
//
// starmap.py's numbers are scene units around the sector's center with +y
// down the galaxy's y axis (the Sector Map's screen up). The builder
// places them: `options.origin` (where the sector's center is in the
// map's world), `options.unit` (world units per scene unit) and
// `options.flipY` (turn y back to the galaxy's own); with none of them the
// scene is the Sector Map's own. Each entry's x, y, z and r become world
// values (what the picker and the rings use), and `local` its place in the
// group.

const VERSION_QUERY = new URL(import.meta.url).search;
const THREE = await import(`./vendor/three.module.min.js${VERSION_QUERY}`);
const { makeGlowMaterial } = await import(`./bodyRendering.js${VERSION_QUERY}`);
const { generateButtons } = await import(`./generatebuttons.js${VERSION_QUERY}`);
const { cssVar, formatAddress, makeRingTexture } = await import(`./mapcore.js${VERSION_QUERY}`);
const { endpointBookmark } = await import(`./mappick.js${VERSION_QUERY}`);
const { buildNebulaMesh, createShapeLoader, disposeNebulaMesh } = await import(`./nebulamesh.js${VERSION_QUERY}`);

// A cloud entry carries `kind` (its phenomenon-texture recipe, see
// `CLOUD_KIND_RECIPES` below); a star entry never does; a neighboring-
// sector indicator carries `isNeighbor` -- that alone is enough to tell
// all three apart. The panel itself is mappick.js's, shared with the
// Galaxy Map (MAP.65); this says what goes in it.
// `navPick` (navpick.js) gives the NAV buttons and says whether a course is
// being picked; without one the panel has neither.
export function infoSpec(entry, data, navPick) {
  if (entry.isNeighbor) {
    return neighborSpec(entry, data);
  }
  var picking = !!(navPick && navPick.active());
  var spec = { title: entry.name || "Unknown", nav: [], links: [] };
  if (entry.kind) {
    spec.fields = [["Type", entry.typeLabel], ["Radius", entry.radiusText], ["Distance", entry.distanceText]];
    spec.bookmark = endpointBookmark(entry.key, entry.name, entry.href);
    if (navPick) spec.nav = navPick.actionsFor(entry.key, entry.name);
    if (!picking) {
      spec.links.push({ href: entry.href, label: "View phenomenon →" });
    }
    return spec;
  }
  spec.fields = [["Star type", entry.starType], ["Temperature", entry.temp], ["Octant", entry.quadrant],
    ["Location", entry.location]];
  spec.bookmark = endpointBookmark(entry.endpoint, entry.name, entry.href);
  if (navPick) spec.nav = navPick.actionsFor(entry.endpoint, entry.name);
  if (!picking) {
    spec.links.push({ href: entry.href, label: "View system →" });
  }
  return spec;
}

function neighborSpec(entry, data) {
  var fields = [["Address", formatAddress(entry.ringIndex, entry.layerIndex, entry.ringSlotIndex)],
    ["Designation", entry.designation]];
  if (entry.brightStarCount) {
    // Bright stars the plan pre-placed here, waiting for the fill.
    var more = entry.brightStarCount - entry.brightStars.length;
    fields.push(["Bright stars waiting", entry.brightStars.join("; ") + (more > 0 ? "; and " + more + " more" : "")]);
  }
  var spec = { title: entry.exists ? entry.name || "Unnamed sector" : "Uncharted", fields: fields };
  if (entry.exists) {
    spec.bookmark = {
      kind: "sector", value: entry.designation, name: entry.name || entry.designation, url: entry.href || null,
      sectorId: entry.sectorId != null ? entry.sectorId : null,
    };
    spec.links = [{ href: entry.href, label: "View sector →" }];
  } else if (data && data.generate) {
    spec.generate = generateButtons(data.generate, entry.ringIndex, entry.layerIndex, entry.ringSlotIndex, data.edgeLy);
  }
  return spec;
}

// The hover tooltip's text: what it is and its name.
export function tooltipText(entry) {
  if (entry.isNeighbor) {
    return (entry.exists ? entry.name || "Unnamed sector" : "Sector " + entry.designation + ", uncharted")
      + " (neighboring sector)";
  }
  if (entry.kind) {
    return (entry.name || "Unknown") + (entry.typeLabel ? ", " + entry.typeLabel : "");
  }
  return (entry.name || "Unknown") + (entry.starType ? ", " + entry.starType : "");
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

// Reproduces planetgen/web/maps/starmap.py's old (now retired) `_ASTEROID_FIELD_BACKGROUND`
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
// per-instance data beyond which kind/descriptor it is) -- planetgen/web/maps/starmap.py
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
    // planetgen/web/maps/starmap.py's own data["coreColor"]) -- stripped to a plain
    // 6-digit color for the glow shader's own uniform, which controls
    // opacity itself via glowStrength rather than a texture alpha
    // channel.
    var core = (cloud.coreColor || "#c9a8e090").slice(0, 7);
    return { color: core, power: 2.2, strength: 1.3, scale: 1.3 };
  }
  return CLOUD_GLOW_RECIPES[cloud.kind] || CLOUD_GLOW_RECIPES.asteroidField;
}

// How much fainter a neighboring sector's cloud reaching into this one
// is drawn (planetgen/web/maps/starmap.py marks it `neighbor`, MAP.45), so it reads as
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
// peak opacity (a nebula's come from planetgen/web/maps/starmap.py's per-type coreColor).
var CLOUD_VOLUMES = {
  nebula: { innerRatio: 0, density: 1.6 },
  supernovaRemnant: { innerRatio: 0.82, density: 3.0, color: "#ff8a5c", alpha: 0.5 },
};

function makeCloudVolume(cloud, volume) {
  var color = volume.color || (cloud.coreColor || "#c9a8e090").slice(0, 7);
  var alpha = volume.alpha;
  if (alpha === undefined) {
    // coreColor's last two hex digits are its alpha (planetgen/web/maps/starmap.py).
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
  mesh.position.fromArray(cloud.local);
  // Only its far wall is drawn, so a star inside the cloud (nearer than
  // that wall) is never tinted over and stays easy to pick out.
  return mesh;
}

// --- Points of light (MAP.15) ----------------------------------------------
//
// Every star, and every phenomenon planetgen/web/maps/starmap.py sends a `light` for, is
// drawn the way static/galaxymap3d.js draws its bright stars: a point
// sprite a fixed number of pixels across (never sized by distance), a core
// `corePx` across in a soft halo `sizePx` across whose strength is `glow`,
// the core `bright` opaque and whitened toward its middle (by `whiten`, 1
// unless the entry says otherwise: an accreting black hole's core keeps
// its disc's orange). planetgen/web/maps/starmap.py works every number out (from the
// star's radius, luminosity and temperature, or a fixed recipe for a
// phenomenon). The halo falls off a little more slowly than the Galaxy
// Map's, so even a red dwarf shows a bright aura on this smaller map (Boss:
// "realistic sizes with bright auras"). The only change with zoom: closing in grows a point up to
// POINT_CLOSE_GROWTH times its size. On the light theme's pale background
// a white core would vanish, so there the core is a darker shade of its
// own color instead, and the halo blends as usual.
export const POINT_CLOSE_GROWTH = 1.5;
// A highlighted kind of phenomenon (MAP.123) is drawn this many times larger.
var HIGHLIGHT_GROWTH = 2.2;
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
    positions.set(entry.local, 3 * i);
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

function makeTextSprite(text, color, worldHeight) {
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
  sprite.scale.set(worldHeight * (canvasEl.width / canvasEl.height), worldHeight, 1);
  return sprite;
}

// --- The kinds of object, shown or hidden (MAP.79) -----------------------------

// The kinds a map can show or hide, in the order its buttons come, with
// their button names. A star has no `kind`; an accreting and a quiet black
// hole are one kind here.
export var KINDS = [
  ["star", "Stars"], ["nebula", "Nebulae"], ["supernovaRemnant", "Supernova remnants"],
  ["asteroidField", "Asteroid fields"], ["blackHole", "Black holes"], ["neutronStar", "Neutron stars"],
  ["quasar", "Quasars"], ["roguePlanet", "Rogue planets"], ["interstellarComet", "Interstellar comets"],
  ["neighbor", "Neighboring sectors"],
];

// The kind (one of KINDS' keys) of an entry.
export function kindOf(entry) {
  if (entry.isNeighbor) return "neighbor";
  if (!entry.kind) return "star";
  if (entry.kind === "blackHoleAccreting" || entry.kind === "blackHoleQuiescent") return "blackHole";
  return entry.kind;
}

// The star classes the filter offers (MAP.123), "other" for the rest
// (white dwarfs and anything the type string does not start with a class).
export var STAR_CLASSES = ["O", "B", "A", "F", "G", "K", "M", "other"];

// The class (one of STAR_CLASSES) of a star entry.
export function starClassOf(entry) {
  var letter = ((entry && entry.starType) || "").charAt(0).toUpperCase();
  return STAR_CLASSES.indexOf(letter) >= 0 && letter !== "" ? letter : "other";
}

// --- Building a sector ------------------------------------------------------

// How big a neighboring sector's or compass label is, in scene units.
var LABEL_HEIGHT = 22;

// The scene for `data` (starmap.py's JSON, changed in place: see the top
// of this file). `options`:
//   origin, unit, flipY  where it goes (see the top of this file);
//   accentColor          the accent (default the page's --accent);
//   lightBackground      whether the map's background is light;
//   pixelRatio           the renderer's device pixels per CSS pixel;
//   sizeScale()          how much the points have grown with the zoom
//                        (default 1).
// Returns {group, pointEntries, entries (every entry, list order), layers
// (the picker layers, nearest first: points, bodies, volumes), ringSize(entry),
// setRoguesMarked(on), roguesMarked(), entryByKey, kinds() (the kinds in
// it, [{kind, label, count}]), setKindHidden(kind, hidden), kindHidden(kind),
// update() (every frame), nebulaMeshCount() (nebulae drawn from their shape
// so far), dispose()}. `options.onShape(count)` is called as each arrives.
export function buildSectorScene(data, options) {
  var o = options || {};
  var origin = o.origin || [0, 0, 0];
  var unit = o.unit || 1;
  var flip = o.flipY ? -1 : 1;
  var sizeScale = o.sizeScale || function () { return 1; };
  var accentColor = o.accentColor || cssVar("--accent", "#4f5fe8");

  function toLocal(x, y, z) {
    return [x * unit, flip * y * unit, z * unit];
  }

  function place(entry) {
    var local = toLocal(entry.x, entry.y, entry.z);
    entry.local = local;
    entry.x = origin[0] + local[0];
    entry.y = origin[1] + local[1];
    entry.z = origin[2] + local[2];
    if (entry.r != null) entry.r *= unit;
  }

  var group = new THREE.Group();
  group.position.set(origin[0], origin[1], origin[2]);

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
        return new THREE.Vector3().fromArray(toLocal(point[0], point[1], point[2]));
      }));
      group.add(new THREE.Line(outlineGeometry, outlineMaterial));
    });
  }

  if (data.compass) {
    var tip = toLocal(data.compass.tip[0], data.compass.tip[1], data.compass.tip[2]);
    var arrowGeometry = new THREE.BufferGeometry().setFromPoints([
      new THREE.Vector3(0, 0, 0),
      new THREE.Vector3().fromArray(tip),
    ]);
    group.add(new THREE.Line(arrowGeometry, new THREE.LineBasicMaterial({ color: new THREE.Color(accentColor) })));

    // Plain "N" at the arrow's own tip (the standard compass-rose
    // convention) -- no extra arrow glyph appended to the text itself,
    // since the line already drawn above IS the arrow.
    var compassLabel = makeTextSprite(data.compass.label, accentColor, LABEL_HEIGHT * unit);
    compassLabel.position.fromArray(tip);
    group.add(compassLabel);
  }

  var stars = data.stars || [];
  // A binary is one system to pick (MAP.136): a companion's pick is its primary's.
  stars.forEach(function (star) {
    if (star.companionOf != null && stars[star.companionOf]) star.pickAs = stars[star.companionOf];
  });
  var clouds = data.clouds || [];
  var neighbors = data.neighbors || [];
  stars.forEach(place);
  clouds.forEach(place);
  // A neighbor's label sits under its marker (scene units, before placing).
  var labelAt = new Map();
  neighbors.forEach(function (neighbor) {
    labelAt.set(neighbor, toLocal(neighbor.x, neighbor.y - neighbor.r * 2.4, neighbor.z));
    place(neighbor);
  });

  // Every clickable marker but the points of light, for the raycast.
  var bodyObjects = [];
  var volumeObjects = [];
  // A nebula is first the analytic sphere volume (above); once its shape
  // arrives (MAP.103) the shape's mesh takes its place for drawing and picking.
  var loadNebulaShape = data.nebulaShapePath ? createShapeLoader(data.nebulaShapePath) : null;
  var nebulaMeshes = 0;
  var disposed = false;

  function drawFromShape(cloud, sphere) {
    loadNebulaShape(cloud.nebulaId, "low").then(function (shape) {
      if (disposed) return;
      var core = (cloud.coreColor || "#c9a8e090").slice(0, 7);
      var edge = parseInt((cloud.edgeColor || "#c9a8e030").slice(7, 9), 16);
      var look = [core, 0, (isNaN(edge) ? 0x30 : edge) / 255 * (cloud.neighbor ? NEIGHBOR_CLOUD_DIM : 1)];
      var mesh = buildNebulaMesh(THREE, shape, cloud.r, look, { yScale: -flip, depthTest: true });
      mesh.position.fromArray(cloud.local);
      mesh.visible = sphere.visible;
      group.remove(sphere);
      group.add(mesh);
      volumeObjects[volumeObjects.indexOf(sphere)] = mesh;
      entryByObject.delete(sphere);
      entryByObject.set(mesh, cloud);
      var drawn = objectsOfEntry.get(cloud);
      drawn[drawn.indexOf(sphere)] = mesh;
      sphere.geometry.dispose();
      sphere.material.dispose();
      nebulaMeshes++;
      if (o.onShape) o.onShape(nebulaMeshes);
    }, function () {});
  }
  var entryByObject = new Map();
  // Everything drawn for an entry that isn't a point of light, so a hidden
  // kind (MAP.79) hides it all, its label and ring too.
  var objectsOfEntry = new Map();
  function drawnFor(entry, object) {
    if (!objectsOfEntry.has(entry)) objectsOfEntry.set(entry, []);
    objectsOfEntry.get(entry).push(object);
  }
  var hiddenKinds = new Set();
  // MAP.123: kinds of phenomenon drawn larger and brighter to find them.
  var markedKinds = new Set();
  // MAP.123: star classes left off, and the dimmest star shown (L☉; 0 shows all).
  var hiddenClasses = new Set();
  var minLuminosity = 0;

  // Whether `entry` is left off the map: its kind is hidden, or it is a star
  // of a class that is, or one dimmer than the floor.
  function hides(entry) {
    var kind = kindOf(entry);
    if (hiddenKinds.has(kind)) return true;
    if (kind !== "star") return false;
    if (hiddenClasses.has(starClassOf(entry))) return true;
    return minLuminosity > 0 && !((entry.luminositySol || 0) >= minLuminosity);
  }

  function anyHidden() {
    return hiddenKinds.size > 0 || hiddenClasses.size > 0 || minLuminosity > 0;
  }

  // Every star and light-giving phenomenon, drawn as points of light
  // (MAP.15) and picked on screen, not by raycast.
  var pointEntries = stars.concat(clouds.filter(function (cloud) {
    return cloud.light;
  }));
  var pointsOfLight = null;
  if (pointEntries.length) {
    pointsOfLight = makePointsOfLight(pointEntries, !!o.lightBackground);
    pointsOfLight.material.uniforms.pixelRatio.value = o.pixelRatio || 1;
    group.add(pointsOfLight);
  }

  // Not raycast: a click on a marker's empty middle should still reach
  // whatever is behind it. Rogue planets are points of light
  // (planetgen/web/maps/starmap.py's _ROGUE_LIGHT); "Mark rogue planets" (off by
  // default) swaps in each one's markedLight and shows these rings.
  var rogueMarkers = new THREE.Group();
  rogueMarkers.visible = false;
  group.add(rogueMarkers);
  var rogueMarkerMaterial = null;
  var roguesMarked = false;

  var POINT_ATTRIBUTES = ["pointColor", "pointSize", "pointCore", "pointGlow", "pointBright", "pointWhiten"];

  // Writes point `i`'s look: its entry's current light, or nothing at all
  // while its kind is hidden.
  function refreshPoint(i) {
    var geometry = pointsOfLight.geometry;
    var entry = pointEntries[i];
    var light = entry.light;
    var shown = !hides(entry);
    var dim = entry.neighbor ? NEIGHBOR_CLOUD_DIM : 1;
    var color = new THREE.Color(light.color || "#ffffff");
    geometry.getAttribute("pointColor").array.set([color.r, color.g, color.b], 3 * i);
    var mark = markedKinds.has(kindOf(entry)) ? HIGHLIGHT_GROWTH : 1;
    geometry.getAttribute("pointSize").array[i] = shown ? light.sizePx * mark : 0;
    geometry.getAttribute("pointCore").array[i] = shown ? light.corePx * mark : 0;
    geometry.getAttribute("pointGlow").array[i] = shown ? light.glow * dim * mark : 0;
    geometry.getAttribute("pointBright").array[i] = shown ? Math.min(1, light.bright * dim * mark) : 0;
    geometry.getAttribute("pointWhiten").array[i] = light.whiten != null ? light.whiten : 1;
  }

  function pointsChanged() {
    POINT_ATTRIBUTES.forEach(function (name) {
      pointsOfLight.geometry.getAttribute(name).needsUpdate = true;
    });
  }

  function setRoguesMarked(marked) {
    roguesMarked = marked;
    rogueMarkers.visible = marked;
    if (!pointsOfLight) {
      return;
    }
    pointEntries.forEach(function (entry, i) {
      if (entry.kind !== "roguePlanet" || !entry.markedLight) {
        return;
      }
      entry.light = marked ? entry.markedLight : entry.unmarkedLight;
      refreshPoint(i);
    });
    pointsChanged();
  }

  clouds.forEach(function (cloud) {
    if (cloud.kind === "roguePlanet") {
      cloud.unmarkedLight = cloud.light;
      if (!rogueMarkerMaterial) {
        rogueMarkerMaterial = new THREE.SpriteMaterial({
          map: makeRingTexture(ROGUE_MARKER_COLOR), transparent: true, opacity: 0.75, depthWrite: false,
          sizeAttenuation: false,
        });
      }
      var ring = new THREE.Sprite(rogueMarkerMaterial);
      ring.position.fromArray(cloud.local);
      ring.scale.set(ROGUE_MARKER_SCREEN_SIZE, ROGUE_MARKER_SCREEN_SIZE, 1);
      rogueMarkers.add(ring);
      drawnFor(cloud, ring);
    }
    if (cloud.light) {
      return;
    }
    var volume = CLOUD_VOLUMES[cloud.kind];
    if (volume) {
      var mesh = makeCloudVolume(cloud, volume);
      group.add(mesh);
      entryByObject.set(mesh, cloud);
      volumeObjects.push(mesh);
      drawnFor(cloud, mesh);
      if (cloud.kind === "nebula" && cloud.nebulaId != null && loadNebulaShape) {
        drawFromShape(cloud, mesh);
      }
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
    bodies.core.position.fromArray(cloud.local);
    bodies.glow.position.fromArray(cloud.local);
    group.add(bodies.core);
    group.add(bodies.glow);
    bodyObjects.push(bodies.core);
    entryByObject.set(bodies.core, cloud);
    drawnFor(cloud, bodies.core);
    drawnFor(cloud, bodies.glow);
  });

  // Neighboring-sector indicators: a small flat dot at the scene's own
  // edge, in the real direction of that neighbor (see planetgen/web/maps/starmap.py's
  // own `_neighbor_indicator_data`) -- a plain billboard `THREE.Sprite`
  // (not a real 3D body like a star/cloud above) suits these fine: flat
  // 2D markers with no "seen from any angle" concern.
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

  neighbors.forEach(function (neighbor) {
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
    marker.position.fromArray(neighbor.local);
    marker.scale.set(neighbor.r * 2, neighbor.r * 2, 1);
    group.add(marker);
    bodyObjects.push(marker);
    entryByObject.set(marker, neighbor);
    drawnFor(neighbor, marker);

    var labelText = neighbor.exists ? neighbor.name : neighbor.designation;
    if (labelText) {
      var label = makeTextSprite(labelText, neighbor.exists ? NEIGHBOR_EXISTS_COLOR : NEIGHBOR_MISSING_FILL,
        LABEL_HEIGHT * unit);
      label.position.fromArray(labelAt.get(neighbor));
      group.add(label);
      drawnFor(neighbor, label);
    }
  });

  // Picking (mappick.js): a point of light wins unless a solid body or a
  // neighbor's marker sits in front of it; a cloud volume, only ever hit
  // at its far wall (see makeCloudVolume), is picked when nothing else
  // is. A point is a fixed number of pixels across whatever its
  // distance, so it is picked on screen, within POINT_PICK_PX of its
  // center (or its own core, if bigger).
  var entryOf = function (hit) { return entryByObject.get(hit.object) || null; };
  // What of `objects` is not hidden with its kind (a hidden mesh isn't picked).
  function shown(objects) {
    return anyHidden() ? objects.filter(function (object) { return object.visible; }) : objects;
  }
  var layers = [
    {
      name: "sector-points",
      points: function () {
        return anyHidden() ? pointEntries.filter(function (entry) { return !hides(entry); }) : pointEntries;
      },
      reach: function (entry) {
        return entry.kind === "roguePlanet"
          ? ROGUE_PICK_PX[roguesMarked ? "marked" : "unmarked"]
          : Math.max(POINT_PICK_PX, (entry.light.corePx / 2) * sizeScale());
      },
    },
    { name: "sector-bodies", occludes: true, meshes: function () { return shown(bodyObjects); }, entryOf: entryOf },
    { name: "sector-volumes", meshes: function () { return shown(volumeObjects); }, entryOf: entryOf },
  ];

  // The selection or hover ring's size round an entry (mappick.js's
  // createRing): a point of light's a size on screen that grows with it,
  // a body's or a cloud's round it in the scene.
  function ringSize(entry) {
    if (entry.light) {
      return { px: Math.max(POINT_RING_MIN_PX, entry.light.corePx * 2 + 14) * sizeScale() };
    }
    return { radius: Math.max(entry.r || 8 * unit, 4 * unit) * 1.3 };
  }

  // The sector page's Contents "Show on map" buttons (MAP.46) name a
  // cloud by its `key`.
  var entryByKey = new Map();
  clouds.forEach(function (cloud) {
    if (cloud.key) {
      entryByKey.set(cloud.key, cloud);
    }
  });

  var allEntries = stars.concat(clouds).concat(neighbors);

  // The kinds of object in this sector, in KINDS' order, with how many.
  function kinds() {
    var counts = new Map();
    allEntries.forEach(function (entry) {
      var kind = kindOf(entry);
      counts.set(kind, (counts.get(kind) || 0) + 1);
    });
    return KINDS.filter(function (pair) { return counts.has(pair[0]); }).map(function (pair) {
      return { kind: pair[0], label: pair[1], count: counts.get(pair[0]) };
    });
  }

  // Puts every entry's objects and points as the hidden kinds, classes and
  // floor say (a hidden object cannot be hovered or picked).
  function applyHidden() {
    allEntries.forEach(function (entry) {
      var hidden = hides(entry);
      (objectsOfEntry.get(entry) || []).forEach(function (object) { object.visible = !hidden; });
    });
    if (pointsOfLight) {
      pointEntries.forEach(function (entry, i) { refreshPoint(i); });
      pointsChanged();
    }
  }

  // Hides or shows every object of a kind (MAP.79): its points, bodies,
  // clouds, rings and labels, none of which can then be hovered or picked.
  function setKindHidden(kind, hidden) {
    if (hidden) hiddenKinds.add(kind);
    else hiddenKinds.delete(kind);
    applyHidden();
  }

  // Draws every point of a phenomenon kind larger and brighter, or as it was (MAP.123).
  function setKindMarked(kind, marked) {
    if (marked) markedKinds.add(kind);
    else markedKinds.delete(kind);
    if (!pointsOfLight) return;
    pointEntries.forEach(function (entry, i) {
      if (kindOf(entry) === kind) refreshPoint(i);
    });
    pointsChanged();
  }

  // Hides or shows the stars of one class (MAP.123).
  function setStarClassHidden(starClass, hidden) {
    if (hidden) hiddenClasses.add(starClass);
    else hiddenClasses.delete(starClass);
    applyHidden();
  }

  // Shows only stars at least this luminous (L☉); 0 shows them all.
  function setMinLuminosity(value) {
    minLuminosity = value > 0 ? value : 0;
    applyHidden();
  }

  // The star classes in this sector, in STAR_CLASSES' order, with how many.
  function starClasses() {
    var counts = new Map();
    stars.forEach(function (entry) {
      var c = starClassOf(entry);
      counts.set(c, (counts.get(c) || 0) + 1);
    });
    return STAR_CLASSES.filter(function (c) { return counts.has(c); }).map(function (c) {
      return { starClass: c, count: counts.get(c) };
    });
  }

  // The dimmest and brightest luminosity (L\u2609) of the stars in this
  // sector, or null with none: the luminosity slider's ends.
  function luminosityRange() {
    var lo = Infinity, hi = 0;
    stars.forEach(function (entry) {
      var lum = entry.luminositySol;
      if (lum > 0) { lo = Math.min(lo, lum); hi = Math.max(hi, lum); }
    });
    return hi > 0 ? [lo, hi] : null;
  }

  function update() {
    if (pointsOfLight) {
      pointsOfLight.material.uniforms.sizeScale.value = sizeScale();
    }
  }

  function dispose() {
    disposed = true;
    group.traverse(function (object) {
      if (object.geometry) object.geometry.dispose();
      if (object.material) {
        if (object.material.map) object.material.map.dispose();
        object.material.dispose();
      }
    });
    if (group.parent) group.parent.remove(group);
  }

  return {
    group: group,
    pointEntries: pointEntries,
    entries: allEntries,
    layers: layers,
    ringSize: ringSize,
    setRoguesMarked: setRoguesMarked,
    roguesMarked: function () { return roguesMarked; },
    entryByKey: entryByKey,
    kinds: kinds,
    setKindHidden: setKindHidden,
    kindHidden: function (kind) { return hiddenKinds.has(kind); },
    setKindMarked: setKindMarked,
    kindMarked: function (kind) { return markedKinds.has(kind); },
    setStarClassHidden: setStarClassHidden,
    starClassHidden: function (starClass) { return hiddenClasses.has(starClass); },
    starClasses: starClasses,
    luminosityRange: luminosityRange,
    setMinLuminosity: setMinLuminosity,
    minLuminosity: function () { return minLuminosity; },
    // Whether an entry is left off the map by a hidden kind, class or the floor.
    hides: hides,
    update: update,
    nebulaMeshCount: function () { return nebulaMeshes; },
    dispose: dispose,
  };
}

// The screen-reader list's name for an entry.
export function entryLabel(entry) {
  if (entry.isNeighbor) {
    return entry.exists ? entry.name || "Unnamed sector" : "Uncharted (" + entry.designation + ")";
  }
  return entry.name || "Unknown";
}
