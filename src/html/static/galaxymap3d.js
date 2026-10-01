// html/static/galaxymap3d.js
// (Why three.js and not another renderer: docs/html-interface.md.)
//
// Renders the interactive 3D Galaxy Map (lib/galaxymap3d.py) as a real
// WebGL scene (three.js, vendored at ./vendor/three.module.min.js -- see
// sectormap.js's own docstring for why it's vendored rather than loaded
// from a CDN). Unlike sectormap.js's scene (every star/cloud baked into
// one JSON block at page load, camera always orbiting a fixed origin),
// this camera's own orbit TARGET moves freely through the galaxy -- most
// of what's drawn is fetched live from /galaxy/tiles (the JSON's
// fetchPath) as the camera moves, one fixed cube of space ("tile") at a time,
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
// One content tier: the galaxy's sector grid drawn as a solid of lit
// cylindrical segment blocks -- single sectors up close, blocks of m x m x
// ~m whole sectors (m a power of three) further out, sized so a block
// stays a few pixels across (see ./galaxyprisms.js). Density is computed
// right here from the galaxy's own analytic shape, not fetched, and colors
// each block. Unfilled space is translucent; each tile's `filled` summary
// says how many generated sectors every block holds, and a block grows
// more solid with that share, fully solid once every sector in it is
// generated. Clicking a block (or empty space, which is the sector there)
// shows its address or ranges, counts, coordinates and 8 corners; at one
// sector per block, a generated sector's panel links to its page.
//
// Nebulae and supernova remnants are drawn over the blocks as soft
// translucent spheres their real size (each tile lists the ones reaching
// into it); clicking one shows it and links to its page -- see "Clouds".
//
// Pre-placed bright stars (every star of 500 L☉ or more) are points of
// light a few pixels across with a big soft glow, the same size at every
// zoom; clicking one shows it -- see "Bright stars".
//
// Zooming stays smooth: the blocks for a view are built in a Web Worker
// (./galaxyblocks.js), built views are kept and the next zoom step's are
// prepared ahead, zoom steps glide instead of jumping, and a change of
// block size crossfades -- see "Block sets" and the zoom glide below.
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

// Sibling modules are imported with this module's own `?v=<version>`
// query (html/lib/fmt.py's `static_url`), so they are cached and
// refreshed with the page's script. A plain static `import "./x.js"`
// would drop the query: an update could then leave a stale copy cached,
// and a page that also loaded the same file by its versioned URL would
// get a second, separate instance of it.
const VERSION_QUERY = new URL(import.meta.url).search;
const THREE = await import(`./vendor/three.module.min.js${VERSION_QUERY}`);
const {
  blockAt, blockSectorCount, blockSectorRanges, blockSizeForScale, blockSlotRange, cellCoordinates, cellVertices,
  sectorAddressAt, sectorCellBounds, wedgeLines,
} = await import(`./galaxyprisms.js${VERSION_QUERY}`);
const { CELL_STRIDE, POINT_STRIDE, createBlockScene } = await import(`./galaxyblocks.js${VERSION_QUERY}`);
const { formatDistancePc } = await import(`./distance.js${VERSION_QUERY}`);
const { generateButtons } = await import(`./generatebuttons.js${VERSION_QUERY}`);

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

// A plain <a href> link (GET, bookmarkable, opens in a new tab).
function pageLink(href, label) {
  var link = document.createElement("a");
  link.href = href;
  link.className = "btn";
  link.textContent = label;
  return link;
}

// The sector page's URL: the server's template (sceneData.sectorUrl,
// built with page_url so it follows the sector page wherever it lives)
// with the id filled in.
function sectorUrl(id) {
  return String(sceneData.sectorUrl || "").replace("{id}", encodeURIComponent(id));
}

function formatAddress(ringIndex, layerIndex, slotIndex) {
  return "ring " + ringIndex + " layer " + layerIndex + " slot " + slotIndex;
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
  addField(dl, "Address", entry.ring_index != null ? formatAddress(entry.ring_index, entry.layer_index, entry.ring_slot_index) : null);
  addField(dl, "Designation", entry.designation);
  panel.appendChild(dl);
  if (sceneData.sectorUrl && entry.id != null) {
    panel.appendChild(pageLink(sectorUrl(entry.id), "View sector →"));
  }
}

// A share as a percentage, never rounding a nonzero one down to 0%.
function formatShare(share) {
  if (share > 0 && share < 0.001) {
    return "< 0.1%";
  }
  return (100 * share).toFixed(share < 0.1 ? 1 : 0) + "%";
}

function formatDeg(rad) {
  return ((rad * 180) / Math.PI).toFixed(3) + "°";
}

// provisional_sector_designation's hex code: ring << 33 | (layer + 4096)
// << 20 | slot (BigInt, since it runs past 32 bits).
function sectorDesignation(ring, layer, slot) {
  var packed = (BigInt(ring) << 33n) | (BigInt(layer + 4096) << 20n) | BigInt(slot);
  return packed.toString(16).toUpperCase();
}

// One block's info: a single sector (m == 1) or a block of sectors.
// `cell` is {m, bounds, address} for a sector or {m, bounds, ring, seg,
// slab, density, edgePc, shape} for a block, plus `filled` (generated
// sectors in it) when it came from the drawn blocks. A generated sector's
// own panel is showPlacedInfo.
function showCellInfo(cell) {
  var panel = document.getElementById("galaxymap3d-info");
  if (!panel) {
    return;
  }
  panel.textContent = "";
  var b = cell.bounds;
  var single = cell.m === 1;

  var heading = document.createElement("h3");
  heading.textContent = single ? "Sector cell" : "Sector block (" + cell.m + " sectors a side)";
  panel.appendChild(heading);

  var dl = document.createElement("dl");
  if (single) {
    var a = cell.address;
    addField(dl, "Address", formatAddress(a.ring, a.layer, a.slot));
    addField(dl, "Designation", sectorDesignation(a.ring, a.layer, a.slot));
    if (cell.filled != null) {
      addField(dl, "Generated", cell.filled > 0 ? "Yes" : "Not yet");
    }
  } else {
    var ranges = blockSectorRanges(cell.ring, cell.slab, cell.m);
    addField(dl, "Rings", ranges.ringFirst + "–" + ranges.ringLast);
    addField(dl, "Layers", ranges.layerFirst + "–" + ranges.layerLast);
    // Slot numbers restart in every ring, so give the innermost and
    // outermost member rings' ranges.
    [ranges.ringFirst, ranges.ringLast].forEach(function (ring) {
      var slots = blockSlotRange(cell.ring, cell.seg, cell.m, ring);
      addField(dl, "Slots in ring " + ring, slots.first + "–" + slots.last);
    });
    var total = blockSectorCount(cell.ring, cell.seg, cell.slab, cell.m, cell.edgePc, cell.shape);
    addField(dl, "Sectors", total.toLocaleString());
    if (cell.filled != null) {
      addField(dl, "Generated", cell.filled.toLocaleString() + (cell.filled > 0 && total > 0 ? " (" + formatShare(cell.filled / total) + ")" : ""));
    }
  }
  if (cell.density != null) {
    addField(dl, "Predicted density", cell.density.toFixed(2) + "× local average");
  }
  var coords = cellCoordinates(b);
  var c = coords.cartesian;
  addField(dl, "Center x, y, z", c.map(function (v) { return v.toFixed(1); }).join(", ") + " pc");
  addField(dl, "Cylindrical R, θ, z", coords.cylindrical[0].toFixed(1) + " pc, " + formatDeg(coords.cylindrical[1]) + ", " + coords.cylindrical[2].toFixed(1) + " pc");
  addField(dl, "Spherical r, θ, φ", coords.spherical[0].toFixed(1) + " pc, " + formatDeg(coords.spherical[1]) + ", " + formatDeg(coords.spherical[2]));
  addField(dl, "Radial width", formatDistancePc(b.r1 - b.r0));
  addField(dl, "Height", formatDistancePc(b.z1 - b.z0));
  addField(dl, "Mean arc length", formatDistancePc(((b.r0 + b.r1) / 2) * (b.t1 - b.t0)));
  panel.appendChild(dl);

  var corners = document.createElement("details");
  var summary = document.createElement("summary");
  summary.textContent = "8 corners (x, y, z pc)";
  corners.appendChild(summary);
  var list = document.createElement("ol");
  list.start = 0;
  cellVertices(b).forEach(function (v) {
    var item = document.createElement("li");
    item.textContent = v.map(function (n) { return n.toFixed(2); }).join(", ");
    list.appendChild(item);
  });
  corners.appendChild(list);
  panel.appendChild(corners);

  // A sector that isn't generated yet: for a logged-in admin
  // (sceneData.generate, set by lib/galaxymap3d.py only then), the same
  // Generate buttons the Sector Map gives a neighbor.
  if (single && !(cell.filled > 0) && sceneData.generate) {
    panel.appendChild(generateButtons(sceneData.generate, cell.address.ring, cell.address.layer, cell.address.slot));
  }
}

// --- Clouds: nebulae and supernova remnants --------------------------------

// A nebula's color by its type, the same hues the Sector Map uses
// (lib/starmap.py's _NEBULA_TYPE_COLORS and _NEBULA_TYPE_ALPHA, as core
// and edge opacity). A dark nebula is a near-black silhouette.
var NEBULA_LOOKS = {
  diffuse: ["#e3a6c8", 0.47, 0.14],
  emission: ["#ff6f91", 0.69, 0.25],
  reflection: ["#6fa8ff", 0.63, 0.22],
  planetary: ["#5be8c9", 0.66, 0.24],
  dark: ["#1c1c24", 0.91, 0.56],
};
var DEFAULT_NEBULA_LOOK = ["#c9a8e0", 0.56, 0.22];

function capitalize(text) {
  text = String(text || "");
  return text.charAt(0).toUpperCase() + text.slice(1);
}

// "Emission Nebula", "Supernova Remnant (Shell)".
function cloudTypeLabel(cloud) {
  if (cloud.type === "nebula") {
    return cloud.descriptor ? capitalize(cloud.descriptor) + " Nebula" : "Nebula";
  }
  return "Supernova Remnant" + (cloud.descriptor ? " (" + capitalize(cloud.descriptor) + ")" : "");
}

// The phenomenon page's URL: the server's template (sceneData.phenomenonUrl)
// with the type and id filled in.
function phenomenonUrl(cloud) {
  return String(sceneData.phenomenonUrl || "")
    .replace("{type}", encodeURIComponent(cloud.type))
    .replace("{id}", encodeURIComponent(cloud.id));
}

function showCloudInfo(cloud) {
  var panel = document.getElementById("galaxymap3d-info");
  if (!panel) {
    return;
  }
  panel.textContent = "";
  var heading = document.createElement("h3");
  heading.textContent = cloud.name || cloudTypeLabel(cloud);
  panel.appendChild(heading);
  var dl = document.createElement("dl");
  addField(dl, "Type", cloudTypeLabel(cloud));
  addField(dl, "Class", cloud.class);
  addField(dl, "Radius", formatDistancePc(cloud.radius_pc));
  addField(dl, "Center x, y, z", [cloud.x, cloud.y, cloud.z].map(function (v) { return v.toFixed(1); }).join(", ") + " pc");
  addField(dl, "Distance from core", formatDistancePc(Math.hypot(cloud.x, cloud.y, cloud.z)));
  panel.appendChild(dl);
  if (sceneData.phenomenonUrl) {
    panel.appendChild(pageLink(phenomenonUrl(cloud), "View phenomenon →"));
  }
}

// A star's color from its surface temperature: Tanner Helland's fit to
// blackbody colors, as [r, g, b] in 0..1 (red giants orange, O and B
// stars blue-white).
function starColor(temperatureK) {
  var t = Math.min(40000, Math.max(1000, temperatureK || 5800)) / 100;
  var r = t <= 66 ? 255 : 329.698727446 * Math.pow(t - 60, -0.1332047592);
  var g = t <= 66 ? 99.4708025861 * Math.log(t) - 161.1195681661 : 288.1221695283 * Math.pow(t - 60, -0.0755148492);
  var b = t >= 66 ? 255 : t <= 19 ? 0 : 138.5177312231 * Math.log(t - 10) - 305.0447927307;
  return [r, g, b].map(function (v) { return Math.min(255, Math.max(0, v)) / 255; });
}

// "12,300 L☉", "1.2 million L☉".
function formatLuminosity(sol) {
  if (sol >= 1e6) {
    return (sol / 1e6).toFixed(sol >= 1e7 ? 0 : 1) + " million L☉";
  }
  return Math.round(sol).toLocaleString("en-US") + " L☉";
}

function systemUrl(id) {
  return String(sceneData.systemUrl || "").replace("{id}", encodeURIComponent(id));
}

// A pre-placed bright star (queryDb.galaxy_bright_stars_in_box): what it
// is, and its system's page once its sector is filled.
function showStarInfo(star) {
  var panel = document.getElementById("galaxymap3d-info");
  if (!panel) {
    return;
  }
  panel.textContent = "";
  var heading = document.createElement("h3");
  heading.textContent = "Bright star";
  panel.appendChild(heading);
  var dl = document.createElement("dl");
  addField(dl, "Type", [star.star_type, star.yerkes_class].filter(Boolean).join(" "));
  addField(dl, "Luminosity", formatLuminosity(star.luminosity_sol));
  addField(dl, "Temperature", Math.round(star.temperature_k).toLocaleString("en-US") + " K");
  addField(dl, "Distance from core", formatDistancePc(Math.hypot(star.x, star.y, star.z)));
  addField(dl, "Sector", sectorDesignation(star.ring_index, star.layer_index, star.ring_slot_index));
  addField(dl, "Address", formatAddress(star.ring_index, star.layer_index, star.ring_slot_index));
  addField(dl, "System", star.system_id != null ? null : "Not generated yet (its sector isn't filled)");
  panel.appendChild(dl);
  if (sceneData.systemUrl && star.system_id != null) {
    panel.appendChild(pageLink(systemUrl(star.system_id), "View system →"));
  }
}

function rgba(hex, alpha) {
  var n = parseInt(hex.slice(1), 16);
  return "rgba(" + (n >> 16) + "," + ((n >> 8) & 255) + "," + (n & 255) + "," + alpha + ")";
}

// A cloud's sprite texture, drawn so the sphere's edge is the texture's
// edge (the sprite is scaled to the cloud's diameter): a nebula glows
// from a denser core out to nothing; a remnant is a thin bright shell
// around a faint interior, as Cassiopeia A or the Veil look.
var cloudTextures = new Map();

function cloudTexture(cloud) {
  var key = cloud.type === "nebula" ? "nebula:" + (NEBULA_LOOKS[cloud.descriptor] ? cloud.descriptor : "") : "remnant";
  var texture = cloudTextures.get(key);
  if (texture) {
    return texture;
  }
  var size = 128;
  var canvasEl = document.createElement("canvas");
  canvasEl.width = size;
  canvasEl.height = size;
  var ctx = canvasEl.getContext("2d");
  var gradient = ctx.createRadialGradient(size / 2, size / 2, 0, size / 2, size / 2, size / 2);
  if (cloud.type === "nebula") {
    var look = NEBULA_LOOKS[cloud.descriptor] || DEFAULT_NEBULA_LOOK;
    gradient.addColorStop(0, rgba(look[0], look[1]));
    gradient.addColorStop(0.35, rgba(look[0], look[1] * 0.8));
    gradient.addColorStop(0.75, rgba(look[0], look[2]));
    gradient.addColorStop(1, rgba(look[0], 0));
  } else {
    gradient.addColorStop(0, rgba("#ffb07a", 0.12));
    gradient.addColorStop(0.6, rgba("#ffb07a", 0.2));
    gradient.addColorStop(0.82, rgba("#ff8a5c", 0.7));
    gradient.addColorStop(0.9, rgba("#8fd6ff", 0.6));
    gradient.addColorStop(1, rgba("#8fd6ff", 0));
  }
  ctx.fillStyle = gradient;
  ctx.fillRect(0, 0, size, size);
  texture = new THREE.CanvasTexture(canvasEl);
  texture.colorSpace = THREE.SRGBColorSpace;
  cloudTextures.set(key, texture);
  return texture;
}

// --- Sprite textures -------------------------------------------------------

// A text label for a sprite, with its width / height in userData.aspect.
function makeLabelTexture(text, color) {
  var height = 64;
  var canvasEl = document.createElement("canvas");
  var ctx = canvasEl.getContext("2d");
  var font = "600 " + Math.round(height * 0.7) + "px system-ui, sans-serif";
  ctx.font = font;
  canvasEl.width = Math.ceil(ctx.measureText(text).width + height * 0.3);
  canvasEl.height = height;
  ctx.font = font;
  ctx.fillStyle = color;
  ctx.textAlign = "center";
  ctx.textBaseline = "middle";
  ctx.fillText(text, canvasEl.width / 2, height / 2);
  var texture = new THREE.CanvasTexture(canvasEl);
  texture.colorSpace = THREE.SRGBColorSpace;
  texture.userData.aspect = canvasEl.width / height;
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

// --- Filled-sector colors ---------------------------------------------
//
// Unlike lib/starmap.py (size/color baked server-side, once), these run
// client-side for EVERY fetch (initial payload and every live re-fetch
// alike) -- see galaxymap3d.py's own module docstring for why: almost
// everything drawn here arrives through a live fetch that never passes
// through that Python module again, so a server-computed style would
// only ever apply to the first frame.

// A placed sector's own REAL stellar density (system_count / edge_ly^3),
// relative to physical_constants.LOCAL_STELLAR_DENSITY_LY3 (the real
// local-neighborhood average this whole generator already calibrates
// against -- see lib/galaxymap3d.py's own referenceDensityPerLy3
// comment) -- 1.0 means exactly average, >1 denser, <1 sparser. `null`
// when edge_ly isn't available (a sector placed before per-sector edge
// tracking existed) rather than a false 0, so a generated sector's block
// (galaxyblocks.js) can fall back to a neutral mid-tone instead of
// reading as "empty". Generated sectors are colored on a warm
// bronze-to-gold range from PLACED_LOW_DENSITY_COLOR (sparse) to
// PLACED_HIGH_DENSITY_COLOR (dense), apart from the density shading's
// cooler dim-to-accent one, so they stand out at a glance.
var PLACED_LOW_DENSITY_COLOR = "#4a3f2e";
var PLACED_HIGH_DENSITY_COLOR = "#fff6df";

function placedRelativeDensity(entry, referenceDensityPerLy3) {
  if (!entry.edge_ly || !referenceDensityPerLy3) {
    return null;
  }
  var densityPerLy3 = (entry.system_count || 0) / Math.pow(entry.edge_ly, 3);
  return densityPerLy3 / referenceDensityPerLy3;
}

// The selection ring, drawn this many screen pixels across.
var HIGHLIGHT_PX = 9;

var highlightTexture = null;

function ensureTextures(accentColor) {
  if (!highlightTexture) {
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

  var FOV_DEG = data.fovDeg || 50;
  var farPlane = Math.max(data.maxViewRadiusPc * 6, 1000);
  var camera = new THREE.PerspectiveCamera(FOV_DEG, 1, Math.max(data.minViewRadiusPc / 50, 0.001), farPlane);
  camera.up.set(0, 0, 1);

  var MIN_RADIUS = data.minViewRadiusPc;
  var MAX_RADIUS = data.maxViewRadiusPc;
  var CLICK_FACTOR_MIN = data.clickZoomFactorMin || 1.15;
  var CLICK_FACTOR_MAX = data.clickZoomFactorMax || 4.0;
  // Wheel zoom scales with how far the wheel actually moved (normalized
  // to pixels): one ~100 px mouse-wheel notch is a ~1.28x step, while a
  // trackpad's stream of tiny deltas zooms smoothly instead of taking a
  // full step on every event. Pinch-zoom arrives as ctrl+wheel with even
  // smaller deltas, hence the boost.
  var WHEEL_ZOOM_PER_PX = 0.0025;
  var WHEEL_MAX_PX = 200;
  var PINCH_BOOST = 4;
  var GALAXY_RADIUS = data.galaxyRadiusPc || data.maxViewRadiusPc;
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
    // Raycasts read matrixWorld, which otherwise only updates at the next
    // render -- a click right after a move would pick against the old view.
    camera.updateMatrixWorld();
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
    zoomAnimation = null;
    target.copy(initialTarget);
    orbit.radius = initialOrbit.radius;
    orbit.theta = initialOrbit.theta;
    orbit.phi = initialOrbit.phi;
    applyCamera();
    updateScaleBar();
    scheduleFetch(true);
  }

  // --- Wedge lines ---------------------------------------------------------
  //
  // The master lines of the sector grid (galaxyprisms.wedgeLines): lines in
  // the galactic plane out past the edge, 3 from the core, then 6, 12, ...
  // each starting where its zone does. The coarsest are labelled with their
  // bearing (degrees counterclockwise from +X, the zero meridian), so a
  // view can be placed around the galaxy at a glance. Each zone's lines
  // show only while they are at least WEDGE_MIN_GAP_PX apart on screen
  // (updateWedgeLevels). Drawn over everything (no depth test) and never
  // picked. The Wedges button hides them.
  var wedgeColor = new THREE.Color(cssVar("--text", "#e6e8f0"));
  var wedgeGroup = new THREE.Group();
  wedgeGroup.renderOrder = 2;
  scene.add(wedgeGroup);
  var WEDGE_LABEL_PX = 16;
  // Zones with more master lines than this (15 degrees apart) go unlabelled.
  var WEDGE_LABEL_MAX_MASTERS = 24;
  var WEDGE_MIN_GAP_PX = 24;
  var wedgeLabels = [];
  // [{masters, r0, objects: [LineSegments, label sprites...]}], coarsest first.
  var wedgeLevels = [];
  (function buildWedgeLines() {
    var reach = GALAXY_RADIUS * 1.02;
    var byMasters = new Map();
    wedgeLines(data.edgePc || 1, GALAXY_RADIUS).forEach(function (line) {
      if (!byMasters.has(line.masters)) {
        byMasters.set(line.masters, []);
      }
      byMasters.get(line.masters).push(line);
    });
    var material = new THREE.LineBasicMaterial({
      color: wedgeColor, transparent: true, opacity: 0.85, depthTest: false, depthWrite: false,
    });
    byMasters.forEach(function (lines, masters) {
      var level = { masters: masters, r0: lines[0].r0, lines: lines, objects: [], segments: null };
      var points = new Float32Array(lines.length * 6);
      lines.forEach(function (line, n) {
        var cos = Math.cos(line.angleRad);
        var sin = Math.sin(line.angleRad);
        points.set([line.r0 * cos, line.r0 * sin, 0, reach * cos, reach * sin, 0], 6 * n);
        if (masters > WEDGE_LABEL_MAX_MASTERS) {
          return;
        }
        var label = new THREE.Sprite(new THREE.SpriteMaterial({
          map: makeLabelTexture(String(line.bearingDeg).padStart(3, "0"), "#" + wedgeColor.getHexString()),
          transparent: true, depthTest: false, depthWrite: false, sizeAttenuation: false,
        }));
        label.position.set(reach * 1.06 * cos, reach * 1.06 * sin, 0);
        label.renderOrder = 2;
        wedgeLabels.push(label);
        level.objects.push(label);
      });
      var geometry = new THREE.BufferGeometry();
      geometry.setAttribute("position", new THREE.BufferAttribute(points, 3));
      var segments = new THREE.LineSegments(geometry, material);
      segments.renderOrder = 2;
      segments.frustumCulled = false;
      level.objects.push(segments);
      level.segments = segments;
      level.objects.forEach(function (object) { wedgeGroup.add(object); });
      wedgeLevels.push(level);
    });
  })();

  // Shows each zone's lines while neighbours in that zone are at least
  // WEDGE_MIN_GAP_PX apart on screen.
  // - Labelled zones run out to the galaxy's edge. The gap is measured
  //   where the view reaches farthest out (the focus's radius plus half
  //   the screen's height), or at the zone's start.
  // - Finer zones are measured at the focus (or the zone's start) and are
  //   clipped to the view ball (fadeRadius around the target, in the
  //   plane): drawn out to the edge, their far ends would converge into a
  //   hatch toward the horizon.
  function updateWedgeLevels(pcPerPixel) {
    var focusR = Math.hypot(target.x, target.y);
    var reachR = Math.min(GALAXY_RADIUS, focusR + 0.5 * (canvasEl.clientHeight || 1) * pcPerPixel);
    wedgeLevels.forEach(function (level) {
      var labelled = level.masters <= WEDGE_LABEL_MAX_MASTERS;
      var gapPx = (2 * Math.PI * Math.max(level.r0, labelled ? reachR : focusR)) / level.masters / pcPerPixel;
      var show = level.masters === wedgeLevels[0].masters || gapPx >= WEDGE_MIN_GAP_PX;
      level.objects.forEach(function (object) { object.visible = show; });
      if (show && !labelled) {
        clipWedgeLevel(level);
      }
    });
  }

  // Each of a zone's lines as the part inside the circle where the view
  // ball meets the plane, from where the line starts.
  function clipWedgeLevel(level) {
    var radius = Math.sqrt(Math.max(0, fadeRadius * fadeRadius - target.z * target.z));
    var cx = target.x;
    var cy = target.y;
    var attribute = level.segments.geometry.getAttribute("position");
    var points = attribute.array;
    level.lines.forEach(function (line, n) {
      var cos = Math.cos(line.angleRad);
      var sin = Math.sin(line.angleRad);
      var along = cos * cx + sin * cy;
      var disc = along * along - (cx * cx + cy * cy) + radius * radius;
      var near = 0;
      var far = 0;
      if (disc > 0) {
        near = Math.max(line.r0, along - Math.sqrt(disc));
        far = Math.max(near, along + Math.sqrt(disc));
      }
      points.set([near * cos, near * sin, 0, far * cos, far * sin, 0], 6 * n);
    });
    attribute.needsUpdate = true;
  }

  // Keeps the labels WEDGE_LABEL_PX tall on screen: a sprite without size
  // attenuation spans scale / tan(fov / 2) half-heights of the view.
  function updateWedgeLabels() {
    var heightPx = canvasEl.clientHeight || 1;
    var h = (WEDGE_LABEL_PX * 2 * Math.tan(THREE.MathUtils.degToRad(camera.fov) / 2)) / heightPx;
    wedgeLabels.forEach(function (label) {
      label.scale.set(h * label.material.map.userData.aspect, h, 1);
    });
  }

  // --- The block solid -------------------------------------------------------
  //
  // Everything is drawn as one solid of blocks (galaxyprisms.js). Blocks
  // fill their whole cell, so the galaxy reads as one solid; a thin
  // brighter outline along each face's own edges (found from the face's
  // 0..1 uv, about a screen pixel wide) keeps neighboring blocks apart.
  //
  // One built view (a "block set", see "Block sets" below) is two meshes
  // sharing one shader:
  // - solid: blocks whose every sector is generated, opaque;
  // - glass: everything else, translucent (no depth writes), its blocks
  //   sorted back to front from the camera when built. Unfilled space is
  //   half to 80% opaque by density, and a block grows more solid with
  //   its filled share (galaxyblocks.blockOpacity).
  // Blocks holding filled sectors also get a warm tint on their faces and
  // warm, stronger face edges, in step with their filled share
  // (galaxyblocks.blockFillStep), so they can be picked out even where the
  // opacity step alone is small.
  // Two things depend on where the camera is, so the shader does them and
  // a set built once stays right while the camera zooms and turns:
  // - whole blocks close to the camera are dropped (nearCut, tested on
  //   each block's own center), so flying through the disk shows what's
  //   ahead instead of a wall of the nearest blocks;
  // - blocks darken toward the rim of the view ball (fadeCenter,
  //   fadeRadius: EDGE_FADE_START to EDGE_FADE_END of the way out), so the
  //   solid fades out instead of stopping at a hard spherical rim the eye
  //   reads as a ball.
  // `fade` scales a whole set's opacity, for crossfades. The logdepthbuf
  // chunks match the renderer's logarithmic depth buffer.
  var FILLED_TINT = new THREE.Color(PLACED_HIGH_DENSITY_COLOR);
  var EDGE_FADE_START = 0.6;
  var EDGE_FADE_END = 1.0;

  function makeBlockMaterial(translucent) {
    return new THREE.ShaderMaterial({
      uniforms: {
        nearCut: { value: 0 }, filledTint: { value: FILLED_TINT }, fade: { value: 1 },
        fadeCenter: { value: new THREE.Vector3() }, fadeRadius: { value: 0 },
      },
      transparent: translucent,
      depthWrite: !translucent,
      vertexShader: [
        "#include <common>",
        "#include <logdepthbuf_pars_vertex>",
        "uniform float nearCut;",
        "uniform vec3 fadeCenter;",
        "uniform float fadeRadius;",
        "attribute vec3 prismColor;",
        "attribute vec3 prismCenter;",
        "attribute float prismAlpha;",
        "attribute float prismFill;",
        "attribute vec2 faceUv;",
        "varying vec3 vColor;",
        "varying float vAlpha;",
        "varying float vFill;",
        "varying vec2 vUv;",
        "varying float vKeep;",
        "void main() {",
        "  float rim = 0.0;",
        "  if (fadeRadius > 0.0) {",
        "    rim = clamp((distance(prismCenter, fadeCenter) / fadeRadius - " + EDGE_FADE_START.toFixed(3) + ") / "
          + (EDGE_FADE_END - EDGE_FADE_START).toFixed(3) + ", 0.0, 1.0);",
        "  }",
        "  vColor = prismColor * (0.25 + 0.75 * (1.0 - rim * rim * (3.0 - 2.0 * rim)));",
        "  vAlpha = prismAlpha;",
        "  vFill = prismFill;",
        "  vUv = faceUv;",
        "  vKeep = step(nearCut, distance(cameraPosition, prismCenter));",
        "  vec4 mvPosition = modelViewMatrix * vec4(position, 1.0);",
        "  gl_Position = projectionMatrix * mvPosition;",
        "  #include <logdepthbuf_vertex>",
        "}",
      ].join("\n"),
      fragmentShader: [
        "#include <common>",
        "#include <logdepthbuf_pars_fragment>",
        "varying vec3 vColor;",
        "uniform vec3 filledTint;",
        "uniform float fade;",
        "varying float vAlpha;",
        "varying float vFill;",
        "varying vec2 vUv;",
        "varying float vKeep;",
        "void main() {",
        "  if (vKeep < 0.5) discard;",
        "  #include <logdepthbuf_fragment>",
        "  vec2 toEdge = min(vUv, 1.0 - vUv) / max(fwidth(vUv), vec2(1e-6));",
        "  float edge = 1.0 - smoothstep(0.5, 1.5, min(toEdge.x, toEdge.y));",
        "  vec3 face = mix(vColor, filledTint, 0.35 * vFill);",
        "  vec3 edgeColor = mix(vec3(1.0), filledTint, step(0.001, vFill));",
        "  gl_FragColor = vec4(mix(face, edgeColor, (0.18 + 0.7 * vFill) * edge), vAlpha * fade);",
        "  #include <colorspace_fragment>",
        "}",
      ].join("\n"),
    });
  }

  var highlightSprite = new THREE.Sprite(
    new THREE.SpriteMaterial({ map: highlightTexture, transparent: true, depthWrite: false, depthTest: false })
  );
  highlightSprite.visible = false;
  highlightSprite.renderOrder = 5;
  scene.add(highlightSprite);

  function highlightPosition(x, y, z) {
    highlightSprite.position.set(x, y, z);
    highlightSprite.visible = true;
  }

  // Keeps the selection ring HIGHLIGHT_PX across on screen whatever its
  // distance: a sprite spans its scale in world units, so the scale is the
  // world size of that many pixels at the sprite's own distance.
  function updateHighlightScale() {
    if (!highlightSprite.visible) {
      return;
    }
    var heightPx = canvasEl.clientHeight || 1;
    var worldPerPixel = (2 * Math.tan(THREE.MathUtils.degToRad(camera.fov) / 2) * camera.position.distanceTo(highlightSprite.position)) / heightPx;
    var size = HIGHLIGHT_PX * 2 * worldPerPixel;
    highlightSprite.scale.set(size, size, 1);
  }

  // The radius of the view ball around the target (the fetched tiles'
  // reach, which the solid fades out toward): set by renderFromCache.
  var fadeRadius = 0;

  // --- Filled sectors --------------------------------------------------------
  //
  // Every tile carries a `filled` summary (queryDb.galaxy_filled_in_box):
  // at g = 1 each generated sector's address, id, name and system count;
  // coarser, counts per cell g sectors a side. renderFromCache turns them
  // into points (a sector's center, or a cell's middle) with counts, and
  // galaxyblocks.js sums those into the blocks it draws. A cell never
  // spans more than one block ring or layer (g divides m), so only a
  // cell's wedge can straddle two blocks, and its count goes to the one
  // holding its middle.
  var filledPoints = [];
  var filledVersion = 0;

  function filledPointsOf(filled) {
    var out = [];
    if (!filled || !filled.cells) {
      return out;
    }
    var g = filled.g || 1;
    filled.cells.forEach(function (cell) {
      if (g === 1) {
        var bounds = sectorCellBounds(cell[0], cell[1], cell[2], edgePc);
        var c = cellCoordinates(bounds).cartesian;
        out.push({
          x: c[0], y: c[1], z: c[2], count: 1,
          sector: {
            id: cell[3], system_count: cell[4], name: cell[5],
            ring_index: cell[0], layer_index: cell[1], ring_slot_index: cell[2],
            designation: sectorDesignation(cell[0], cell[1], cell[2]),
            galactic_radius_pc: Math.hypot(c[0], c[1], c[2]), edge_ly: data.edgeLy,
          },
        });
        return;
      }
      var wedges = Math.max(3, Math.round(2 * Math.PI * (cell[0] + 0.5)));
      var r = (cell[0] + 0.5) * g * edgePc;
      var t = ((cell[2] + 0.5) * 2 * Math.PI) / wedges;
      out.push({ x: r * Math.cos(t), y: r * Math.sin(t), z: cell[1] * g * edgePc, count: cell[3] });
    });
    return out;
  }

  // The filled points as galaxyblocks.js takes them (POINT_STRIDE numbers
  // each).
  function packFilledPoints(points) {
    var packed = new Float64Array(points.length * POINT_STRIDE);
    points.forEach(function (point, n) {
      var o = n * POINT_STRIDE;
      packed[o] = point.x;
      packed[o + 1] = point.y;
      packed[o + 2] = point.z;
      packed[o + 3] = point.count;
      packed[o + 4] = point.sector ? 1 : 0;
      var relative = point.sector ? placedRelativeDensity(point.sector, data.referenceDensityPerLy3) : null;
      packed[o + 5] = relative == null ? NaN : relative;
    });
    return packed;
  }

  // With the slice on (the default), blocks above the focus's own layer
  // are left out, so the view looks down on the solid's cut face -- the
  // arms and whatever layer the focus is in -- instead of its outside.
  var sliceAtFocus = true;
  // Density shading's dim-to-accent-to-white ramp (galaxyblocks.js).
  var PRISM_DIM = new THREE.Color(0x1d2340);
  var PRISM_ACCENT = new THREE.Color(accentColor);
  var PRISM_HOT = new THREE.Color(0xeef0ff);
  // Blocks centered nearer the camera than this x the orbit radius (the
  // camera's distance to the target) aren't drawn.
  var NEAR_CUT = 0.4;
  var galaxyShape = data.densityShape || null;
  var edgePc = data.edgePc || 1;
  // How many sectors a side the drawn blocks are.
  var drawnSectorsPerPrism = 1;

  // Parsecs per screen pixel at the view's focus (the orbit target), for
  // the camera at orbit radius `radius` (default: where it is). Per CSS
  // pixel, on purpose, not per device pixel: the block size (a block at
  // least blockMinPx across) and the scale readout's "1 px" are about what
  // a person can see and click, which is the same on a 2x phone screen as
  // on a 1x monitor. The renderer draws at up to 2 device pixels per CSS
  // pixel (setPixelRatio) for sharpness only, and the bright stars' sizes
  // are CSS pixels too (scaled by pixelRatio in their shader).
  function pcPerPixelAtTarget(radius) {
    var heightPx = canvasEl.clientHeight || 1;
    return ((radius || orbit.radius) * 2 * Math.tan(THREE.MathUtils.degToRad(camera.fov) / 2)) / heightPx;
  }

  // --- Block sets ------------------------------------------------------------
  //
  // galaxyblocks.js builds a view's blocks; a Web Worker runs it, so a
  // zoom step never stalls the page. The first frame is built right here
  // (it must not wait on the worker), and so is everything when a worker
  // can't start.
  //
  // A built view is kept as a block set, keyed so that one set serves
  // every nearby camera position (blockViewFor):
  // - the block size m and the view ball's radius, rounded up to one of
  //   VIEW_BUCKETS_PER_DOUBLING steps per doubling;
  // - the ball's center, snapped to a grid VIEW_SNAP_FRACTION of that
  //   radius (the ball grows by the snap's reach, so it still covers the
  //   view);
  // - the slice's block layer and the block layer the camera is in (which
  //   decides the top and bottom faces it can see).
  // The most recent sets stay ready (up to BLOCK_SET_MAX_VERTICES in
  // all), so zooming back is instant. While the page is idle, the worker
  // prebuilds the views one zoom step in and out, and the nearest views
  // with the next finer and coarser block size. When the needed set isn't
  // built yet, the last one stays on screen until it is. A new set with a
  // different block size crossfades in over CROSSFADE_MS; one with the
  // same size (the same blocks, a little more or less of them) swaps.
  // When the generated sectors change (new tiles), a set already built
  // still shows at once and is rebuilt in the background.
  var VIEW_BUCKETS_PER_DOUBLING = 8;
  var VIEW_SNAP_FRACTION = 1 / 32;
  var BLOCK_SET_MAX_VERTICES = 4000000;
  var CROSSFADE_MS = 220;
  // A padded ball holds more blocks than the view's own, so its budget
  // (galaxyprisms.blocksForView's) grows by the same share of surface,
  // and the block size stays what the view's own ball would get.
  var VIEW_BUCKET_STEP = Math.pow(2, 1 / VIEW_BUCKETS_PER_DOUBLING);
  var VIEW_PAD = VIEW_BUCKET_STEP * (1 + VIEW_SNAP_FRACTION * Math.sqrt(3) / 2);
  var VIEW_BUDGET = Math.round((data.blockBudget || 60000) * VIEW_PAD * VIEW_PAD);

  var reducedMotion = typeof window.matchMedia === "function"
    && window.matchMedia("(prefers-reduced-motion: reduce)").matches;

  var blockSceneConfig = {
    edgePc: edgePc, galaxyRadius: GALAXY_RADIUS, shape: galaxyShape, minPx: data.blockMinPx, budget: data.blockBudget,
    palette: {
      dim: PRISM_DIM.toArray(), accent: PRISM_ACCENT.toArray(), hot: PRISM_HOT.toArray(),
      placedLow: new THREE.Color(PLACED_LOW_DENSITY_COLOR).toArray(),
      placedHigh: new THREE.Color(PLACED_HIGH_DENSITY_COLOR).toArray(),
    },
  };
  var localBlockScene = createBlockScene(blockSceneConfig);
  var localFilledVersion = -1;

  // The view ball for the camera at orbit radius `radius` around `center`
  // (the fetched tiles' reach).
  function viewRadiusFor(radius) {
    return Math.max(MIN_RADIUS, Math.min(MAX_RADIUS * FETCH_RADIUS_FACTOR, radius * FETCH_RADIUS_FACTOR));
  }

  // The block set the camera at orbit radius `radius` around `center`
  // (current angles) needs: {key, m, view}, view as galaxyblocks' build
  // takes it.
  function blockViewFor(radius, center) {
    var pcPerPixel = pcPerPixelAtTarget(radius);
    var m = blockSizeForScale(pcPerPixel, edgePc, data.blockMinPx, GALAXY_RADIUS);
    var size = m * edgePc;
    var bucket = Math.ceil(Math.log2(viewRadiusFor(radius)) * VIEW_BUCKETS_PER_DOUBLING - 1e-9);
    var ball = Math.pow(2, bucket / VIEW_BUCKETS_PER_DOUBLING);
    var snap = ball * VIEW_SNAP_FRACTION;
    var c = [Math.round(center.x / snap) * snap, Math.round(center.y / snap) * snap, Math.round(center.z / snap) * snap];
    var slice = sliceAtFocus ? Math.round(center.z / size) : null;
    var eye = offsetFromOrbit({ radius: radius, theta: orbit.theta, phi: orbit.phi }).add(center);
    var viewerLayer = Math.floor(eye.z / size + 0.5);
    return {
      key: [m, bucket, c.join(","), slice, viewerLayer].join(":"),
      m: m,
      view: {
        center: c, viewRadius: ball + snap * Math.sqrt(3) / 2, pcPerPixel: pcPerPixel,
        sliceZ: slice === null ? null : slice * size, viewerZ: viewerLayer * size, budget: VIEW_BUDGET,
        eye: [eye.x, eye.y, eye.z],
      },
    };
  }

  var blockSets = new Map();
  var blockSetVertices = 0;
  var shownSet = null;
  var fadingSet = null;
  var fadeStartedAt = 0;

  // A block set from galaxyblocks' build result: [solid mesh, glass mesh],
  // each keeping its packed part (for picking) and the filled points the
  // build counted (a sector's panel comes from them).
  function makeBlockSet(key, built, version, points) {
    var set = { key: key, m: built.m, version: version, points: points, meshes: [], vertexCount: 0 };
    [built.solid, built.glass].forEach(function (part, n) {
      var geometry = new THREE.BufferGeometry();
      geometry.setAttribute("position", new THREE.BufferAttribute(part.positions, 3));
      geometry.setAttribute("prismCenter", new THREE.BufferAttribute(part.centers, 3));
      geometry.setAttribute("prismColor", new THREE.BufferAttribute(part.colors, 3, true));
      geometry.setAttribute("faceUv", new THREE.BufferAttribute(part.uvs, 2, true));
      geometry.setAttribute("prismAlpha", new THREE.BufferAttribute(part.alphas, 1, true));
      geometry.setAttribute("prismFill", new THREE.BufferAttribute(part.fills, 1, true));
      geometry.setIndex(new THREE.BufferAttribute(part.indices, 1));
      var mesh = new THREE.Mesh(geometry, makeBlockMaterial(n === 1));
      mesh.frustumCulled = false;
      mesh.visible = part.vertexCount > 0;
      mesh.userData.part = part;
      mesh.userData.set = set;
      set.meshes.push(mesh);
      set.vertexCount += part.vertexCount;
    });
    return set;
  }

  function disposeSet(set) {
    set.meshes.forEach(function (mesh) {
      mesh.geometry.dispose();
      mesh.material.dispose();
    });
  }

  // Takes a set off screen, and frees it unless the cache still holds it.
  function retireSet(set) {
    set.meshes.forEach(function (mesh) { scene.remove(mesh); });
    if (blockSets.get(set.key) !== set) {
      disposeSet(set);
    }
  }

  function cacheSet(set) {
    var old = blockSets.get(set.key);
    if (old) {
      blockSets.delete(set.key);
      blockSetVertices -= old.vertexCount;
      if (old !== shownSet && old !== fadingSet) {
        disposeSet(old);
      }
    }
    blockSets.set(set.key, set);
    blockSetVertices += set.vertexCount;
    // Oldest first, never what's on screen.
    var keys = Array.from(blockSets.keys());
    for (var i = 0; i < keys.length && blockSetVertices > BLOCK_SET_MAX_VERTICES; i++) {
      var doomed = blockSets.get(keys[i]);
      if (doomed === set || doomed === shownSet || doomed === fadingSet) {
        continue;
      }
      blockSets.delete(keys[i]);
      blockSetVertices -= doomed.vertexCount;
      disposeSet(doomed);
    }
  }

  function cachedSet(key) {
    return memoryGet(blockSets, key);
  }

  // Puts a set on screen: crossfading from the last one when the block
  // size changed (and motion is welcome), else swapping.
  function showSet(set) {
    if (set === shownSet) {
      return;
    }
    if (fadingSet) {
      retireSet(fadingSet);
      fadingSet = null;
    }
    var old = shownSet;
    shownSet = set;
    set.meshes.forEach(function (mesh) { scene.add(mesh); });
    if (old) {
      if (old.m !== set.m && !reducedMotion) {
        fadingSet = old;
        fadeStartedAt = performance.now();
      } else {
        retireSet(old);
      }
    }
    stepCrossfade(performance.now());
    if (drawnSectorsPerPrism !== set.m) {
      drawnSectorsPerPrism = set.m;
      updateScaleBar();
    }
  }

  // Fades the shown set in and the last one out. While a set fades, even
  // its solid blends, and the outgoing one draws last and writes no depth.
  function stepCrossfade(now) {
    if (!shownSet) {
      return;
    }
    var t = fadingSet ? Math.min(1, (now - fadeStartedAt) / CROSSFADE_MS) : 1;
    setFade(shownSet, t, false);
    if (fadingSet) {
      if (t >= 1) {
        retireSet(fadingSet);
        fadingSet = null;
      } else {
        setFade(fadingSet, 1 - t, true);
      }
    }
  }

  function setFade(set, fade, outgoing) {
    set.meshes.forEach(function (mesh, n) {
      var translucent = n === 1;
      mesh.material.uniforms.fade.value = fade;
      mesh.material.transparent = translucent || fade < 1;
      mesh.material.depthWrite = !translucent && !outgoing;
      // Opaque first, then the translucent shell over it; an outgoing
      // set over both.
      mesh.renderOrder = (outgoing ? 2 : 0) + n;
    });
  }

  // The camera-dependent uniforms, on every set on screen.
  function updateBlockUniforms() {
    [shownSet, fadingSet].forEach(function (set) {
      if (!set) {
        return;
      }
      set.meshes.forEach(function (mesh) {
        var uniforms = mesh.material.uniforms;
        // Inside the solid, the slice is the main answer; this near cut
        // stays for views where the camera is below the cut (or the slice
        // is off).
        uniforms.nearCut.value = orbit.radius * NEAR_CUT;
        uniforms.fadeCenter.value.copy(target);
        uniforms.fadeRadius.value = fadeRadius;
      });
    });
  }

  var blockWorker = startBlockWorker();
  var workerFilledVersion = -1;
  var buildInFlight = null;
  var nextBuildId = 1;

  function startBlockWorker() {
    if (typeof Worker === "undefined") {
      return null;
    }
    var worker;
    try {
      worker = new Worker(new URL("./galaxyblocks.js" + VERSION_QUERY, import.meta.url), { type: "module" });
    } catch (err) {
      return null;
    }
    worker.onmessage = function (event) {
      var message = event.data;
      if (!buildInFlight || message.id !== buildInFlight.id) {
        return;
      }
      var job = buildInFlight;
      buildInFlight = null;
      if (message.type === "built") {
        var set = makeBlockSet(job.key, message.result, job.version, job.points);
        cacheSet(set);
        if (job.ahead && !zoomAnimation) {
          warmUp(set);
        }
      } else {
        abandonWorker();
      }
    };
    worker.onerror = function () {
      abandonWorker();
    };
    worker.postMessage({ type: "init", config: blockSceneConfig });
    return worker;
  }

  // From here on every set is built on the page.
  function abandonWorker() {
    if (blockWorker) {
      blockWorker.terminate();
    }
    blockWorker = null;
    buildInFlight = null;
  }

  function buildHere(wanted) {
    if (localFilledVersion !== filledVersion) {
      localFilledVersion = filledVersion;
      localBlockScene.setFilled(packFilledPoints(filledPoints));
    }
    var set = makeBlockSet(wanted.key, localBlockScene.build(wanted.view), filledVersion, filledPoints);
    cacheSet(set);
    return set;
  }

  // A set's geometry goes to the GPU the first time it's drawn, which for
  // a big set takes a frame or more: a set built ahead is drawn once into
  // a 1-pixel target right away, while the page is idle, so showing it
  // later costs nothing.
  var warmTarget = new THREE.WebGLRenderTarget(1, 1);
  var warmScene = new THREE.Scene();

  function warmUp(set) {
    set.meshes.forEach(function (mesh) { warmScene.add(mesh); });
    var previous = renderer.getRenderTarget();
    renderer.setRenderTarget(warmTarget);
    renderer.render(warmScene, camera);
    renderer.setRenderTarget(previous);
    set.meshes.forEach(function (mesh) { warmScene.remove(mesh); });
  }

  // Asks the worker for a set (`ahead`: one the camera may need next),
  // one at a time: by the time one comes back the camera may have moved
  // on, and the next frame asks for whatever is needed then.
  function buildInWorker(wanted, ahead) {
    if (buildInFlight) {
      return;
    }
    if (workerFilledVersion !== filledVersion) {
      workerFilledVersion = filledVersion;
      var packed = packFilledPoints(filledPoints);
      blockWorker.postMessage({ type: "filled", points: packed }, [packed.buffer]);
    }
    buildInFlight = { id: nextBuildId++, key: wanted.key, version: filledVersion, points: filledPoints, ahead: !!ahead };
    blockWorker.postMessage({ type: "build", id: buildInFlight.id, view: wanted.view });
  }

  // The orbit radii worth having ready: one view step and one click-zoom
  // step each way, and the nearest radii with the next finer and coarser
  // block size.
  function neighborRadii() {
    var r = orbit.radius;
    var factor = clickZoomFactor(r);
    var radii = [r / VIEW_BUCKET_STEP, r * VIEW_BUCKET_STEP, r / factor, r * factor];
    var m = blockViewFor(r, target).m;
    [1 / VIEW_BUCKET_STEP, VIEW_BUCKET_STEP].forEach(function (step) {
      for (var k = 1, at = r * step; k <= 4 * VIEW_BUCKETS_PER_DOUBLING; k++, at *= step) {
        if (at < MIN_RADIUS || at > MAX_RADIUS) {
          break;
        }
        if (blockViewFor(at, target).m !== m) {
          radii.push(at);
          break;
        }
      }
    });
    return radii.filter(function (radius) { return radius >= MIN_RADIUS && radius <= MAX_RADIUS; });
  }

  // Not while tiles are on their way: their sectors would make the set
  // the view needs wait behind a build that's already out of date.
  function prebuildNeighbors() {
    if (!blockWorker || buildInFlight || fetchTimer || activeAbort
        || blockSetVertices + (shownSet ? shownSet.vertexCount : 0) > BLOCK_SET_MAX_VERTICES) {
      return;
    }
    var radii = neighborRadii();
    for (var i = 0; i < radii.length; i++) {
      var wanted = blockViewFor(radii[i], target);
      if (!blockSets.has(wanted.key)) {
        buildInWorker(wanted, true);
        return;
      }
    }
  }

  // Every frame: shows the set the camera needs if it's built, and asks
  // for it (or, once it's showing and current, its neighbours) if not.
  function syncBlocks() {
    var wanted = blockViewFor(orbit.radius, target);
    var set = cachedSet(wanted.key);
    if (!set && (!shownSet || !blockWorker)) {
      if (blockWorker || !zoomAnimation) {
        set = buildHere(wanted);
      }
    }
    if (set) {
      showSet(set);
    }
    if (blockWorker) {
      if (!set || set.version !== filledVersion) {
        buildInWorker(wanted);
      } else {
        prebuildNeighbors();
      }
    } else if (set && set.version !== filledVersion && !zoomAnimation) {
      showSet(buildHere(wanted));
    }
    stepCrossfade(performance.now());
    updateBlockUniforms();
  }

  // --- Cube tiles ----------------------------------------------------------
  //
  // The map asks for fixed cubes of space ("tiles") rather than "everything
  // within R of the camera target" -- see stellarObjects.galaxyViewport's
  // "Cube tiles" section. Space is an octree: level 0 is one cube
  // tileRootEdgePc on a side centered on the galactic origin, each level
  // halves the edge, and a tile's key is "level/ix/iy/iz". neededTiles()
  // picks the smallest level whose tiles are at least the view radius
  // across (so at most 27 tiles cover the view). lib/galaxymap3d.py's
  // initial_tile_request does the same in Python for the first frame.
  // Once the view's own tiles are in, the tiles one click-zoom step in and
  // out are fetched too (prefetchTiles), so zooming never waits on them.
  //
  // A tile's contents depend only on its key and the database's content
  // stamp, so tiles are cached three ways: in memory here, in this
  // browser's localStorage (so a revisit or reload doesn't refetch them),
  // and on the server's disk (lib/tilecache.py) -- only tiles in none of
  // those reach the API and database. Every response carries the current
  // stamp; when it differs from ours, the response also lists which tiles
  // changed since (history), and only those are dropped and refetched. A
  // change the history can't account for (a deleted sector, a new
  // release, a stamp older than the history) drops every cached tile.

  var TILE_ROOT = data.tileRootEdgePc || 65536;
  var TILE_MAX_LEVEL = data.tileMaxLevel != null ? data.tileMaxLevel : 12;
  var FETCH_RADIUS_FACTOR = data.fetchRadiusFactor || 1.6;
  var MAX_TILES_PER_REQUEST = data.maxTilesPerRequest || 128;

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

  // The tiles the camera at orbit radius `radius` (default: where it is)
  // needs.
  function neededTiles(radius) {
    var viewRadius = viewRadiusFor(radius || orbit.radius);
    var level = tileLevelForRadius(viewRadius);
    var keys = tilesIntersectingSphere(level, target, viewRadius);
    return { keys: keys, viewRadius: viewRadius };
  }

  // --- Tile caches ---------------------------------------------------------

  var TILE_MEMORY_MAX = 2000;
  // Keyed by database name, the same as before the page moved to
  // /galaxy, so a visitor's stored tiles carry over.
  var STORAGE_PREFIX = "planetgen:tile:" + data.storageKey + ":";
  // Where the stamp our stored tiles are at is kept, and the generation
  // (the label our stored tiles are filed under since they last all went
  // stale).
  var STAMP_RECORD = STORAGE_PREFIX + "@stamp";
  var remembered = readStampRecord();
  var currentStamp = remembered.stamp || "";
  var currentGeneration = remembered.generation || "";
  var tileMemory = new Map();

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

  function readStampRecord() {
    var store = storage();
    try {
      var record = store ? JSON.parse(store.getItem(STAMP_RECORD) || "{}") : {};
      return record && typeof record === "object" ? record : {};
    } catch (err) {
      return {};
    }
  }

  function rememberStamp() {
    var store = storage();
    if (!store) {
      return;
    }
    try {
      store.setItem(STAMP_RECORD, JSON.stringify({ stamp: currentStamp, generation: currentGeneration }));
    } catch (err) {
      // Next visit just starts afresh.
    }
  }

  // Drops this database's stored tiles, except the current generation's
  // when keepCurrent is set.
  function purgeStoredTiles(keepCurrent) {
    var store = storage();
    if (!store) {
      return;
    }
    try {
      var doomed = [];
      for (var i = 0; i < store.length; i++) {
        var name = store.key(i);
        if (name && name.indexOf(STORAGE_PREFIX) === 0 && name !== STAMP_RECORD) {
          if (!keepCurrent || name.indexOf(STORAGE_PREFIX + currentGeneration + ":") !== 0) {
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
    return STORAGE_PREFIX + currentGeneration + ":" + key;
  }

  function removeStoredTile(key) {
    var store = storage();
    if (!store || !currentGeneration) {
      return;
    }
    try {
      store.removeItem(storageName(key));
    } catch (err) {
      // Nothing more to do.
    }
  }

  function storedTile(key) {
    var store = storage();
    if (!store || !currentGeneration) {
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
    if (!store || !currentGeneration) {
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

  // The tile keys that changed between fromStamp and toStamp, from a
  // response's history (checks oldest first, each {from, to, tiles}), or
  // null when the history doesn't reach back to fromStamp unbroken.
  function staleKeys(fromStamp, toStamp, history) {
    if (!fromStamp || !Array.isArray(history)) {
      return null;
    }
    var start = -1;
    for (var i = history.length - 1; i >= 0; i--) {
      if (history[i] && history[i].from === fromStamp) {
        start = i;
        break;
      }
    }
    if (start < 0) {
      return null;
    }
    var keys = [];
    var at = fromStamp;
    for (var j = start; j < history.length; j++) {
      var entry = history[j];
      if (!entry || entry.from !== at || !Array.isArray(entry.tiles)) {
        return null;
      }
      keys = keys.concat(entry.tiles);
      at = entry.to;
    }
    return at === toStamp ? keys : null;
  }

  function adoptStamp(payload) {
    var stamp = payload.stamp;
    if (!stamp || stamp === currentStamp) {
      return false;
    }
    var stale = currentGeneration ? staleKeys(currentStamp, stamp, payload.history) : null;
    currentStamp = stamp;
    if (stale) {
      stale.forEach(function (key) {
        tileMemory.delete(key);
        removeStoredTile(key);
      });
    } else {
      currentGeneration = payload.generation || stamp;
      tileMemory.clear();
      purgeStoredTiles(true);
    }
    rememberStamp();
    return true;
  }

  function absorb(payload) {
    if (!payload) {
      return false;
    }
    var stampChanged = adoptStamp(payload);
    Object.keys(payload.tiles || {}).forEach(function (key) {
      putTile(key, payload.tiles[key]);
    });
    return stampChanged;
  }

  absorb(data.initial);
  purgeStoredTiles(true);

  // --- Clouds ----------------------------------------------------------------
  //
  // Each tile lists the nebulae and supernova remnants that reach into it
  // (queryDb.galaxy_clouds_in_box); a cloud is a sprite the size of its
  // sphere, drawn over the blocks without a depth test (a flat sprite
  // through a cloud's middle would otherwise be cut in half by the blocks
  // around it, though the cloud fills that space). One that is only a pixel or two across
  // is hidden, and one the camera comes close to fades out so it never
  // fills the screen.
  var CLOUD_MIN_PX = 2;
  var CLOUD_OPACITY = 1;
  var CLOUD_NEAR_FADE = [1.2, 2.0];
  // Clicking picks a cloud only while it is small enough on screen to aim
  // at; up close, clicks go to the sectors inside it.
  var CLOUD_PICK_MAX_PX = 150;
  var cloudGroup = new THREE.Group();
  cloudGroup.renderOrder = 4;
  scene.add(cloudGroup);
  var cloudSprites = new Map();
  var cloudSignature = "";

  function cloudKey(cloud) {
    return cloud.type + ":" + cloud.id;
  }

  // Shows exactly `clouds` (one entry per type and id).
  function setClouds(clouds) {
    var wanted = new Map();
    clouds.forEach(function (cloud) {
      wanted.set(cloudKey(cloud), cloud);
    });
    cloudSprites.forEach(function (sprite, key) {
      if (!wanted.has(key)) {
        cloudGroup.remove(sprite);
        sprite.material.dispose();
        cloudSprites.delete(key);
      }
    });
    wanted.forEach(function (cloud, key) {
      if (cloudSprites.has(key)) {
        return;
      }
      var sprite = new THREE.Sprite(new THREE.SpriteMaterial({
        map: cloudTexture(cloud), transparent: true, depthWrite: false, depthTest: false, opacity: CLOUD_OPACITY,
      }));
      sprite.position.set(cloud.x, cloud.y, cloud.z);
      sprite.scale.set(2 * cloud.radius_pc, 2 * cloud.radius_pc, 1);
      sprite.renderOrder = 4;
      sprite.userData.cloud = cloud;
      cloudGroup.add(sprite);
      cloudSprites.set(key, sprite);
    });
  }

  // A cloud's radius in pixels from where the camera is.
  function cloudPixels(cloud, distance) {
    var heightPx = canvasEl.clientHeight || 1;
    return (cloud.radius_pc * heightPx) / (2 * Math.tan(THREE.MathUtils.degToRad(camera.fov) / 2) * Math.max(distance, 1e-6));
  }

  function updateClouds() {
    cloudSprites.forEach(function (sprite) {
      var cloud = sprite.userData.cloud;
      var distance = camera.position.distanceTo(sprite.position);
      var near = THREE.MathUtils.smoothstep(distance / cloud.radius_pc, CLOUD_NEAR_FADE[0], CLOUD_NEAR_FADE[1]);
      sprite.visible = near > 0 && cloudPixels(cloud, distance) >= CLOUD_MIN_PX;
      sprite.material.opacity = CLOUD_OPACITY * near;
    });
  }

  // The smallest visible cloud under a screen point that is still small
  // enough to aim at, as {cloud, core} (core: the point is within the
  // inner CLOUD_CORE of its radius), or null.
  var CLOUD_CORE = 0.5;

  function cloudAtClientPoint(clientX, clientY) {
    var ndc = ndcFromClientPoint(clientX, clientY);
    if (!ndc) {
      return null;
    }
    raycaster.setFromCamera(ndc, camera);
    var ray = raycaster.ray;
    var best = null;
    var bestOffset = 0;
    cloudSprites.forEach(function (sprite) {
      var cloud = sprite.userData.cloud;
      if (!sprite.visible || (best && best.radius_pc <= cloud.radius_pc)) {
        return;
      }
      var distance = camera.position.distanceTo(sprite.position);
      if (cloudPixels(cloud, distance) > CLOUD_PICK_MAX_PX) {
        return;
      }
      var offset = Math.sqrt(ray.distanceSqToPoint(sprite.position));
      if (offset <= cloud.radius_pc && ray.direction.dot(new THREE.Vector3().subVectors(sprite.position, ray.origin)) > 0) {
        best = cloud;
        bestOffset = offset;
      }
    });
    return best ? { cloud: best, core: bestOffset <= CLOUD_CORE * best.radius_pc } : null;
  }

  // --- Bright stars ----------------------------------------------------------
  //
  // Each tile lists its most luminous pre-placed stars (bright_stars, v43:
  // every star of 500 L☉ or more, placed before any sector is filled), so
  // the arms show before anything is generated. Boss: "no matter how far
  // the user zooms in they should be very small with a big glow" -- each
  // is a point of light a fixed number of pixels across at every zoom
  // (never sized by distance): a core of two or three pixels in a soft
  // halo up to STAR_MAX_PX wide, both a little bigger and brighter for a
  // brighter star. Halos blend normally rather than adding up, so a
  // crowded arm zoomed out glows in its stars' colors instead of burning
  // to white; stars are depth-tested against the blocks without hiding
  // them.
  var STAR_MIN_PX = 12;
  var STAR_MAX_PX = 30;
  var STAR_CORE_PX = [2.2, 3.2];
  var STAR_GLOW = [0.3, 0.6];
  var STAR_LOG_LUMINOSITY = [Math.log10(500), Math.log10(1e6)];
  // A click within this many pixels of a star's center picks it.
  var STAR_PICK_PX = 7;

  var starMaterial = new THREE.ShaderMaterial({
    uniforms: { pixelRatio: { value: renderer.getPixelRatio() } },
    vertexShader: [
      "#include <common>",
      "#include <logdepthbuf_pars_vertex>",
      "attribute float starSize;",
      "attribute float starCore;",
      "attribute float starGlow;",
      "attribute vec3 starColor;",
      "uniform float pixelRatio;",
      "varying vec3 vColor;",
      "varying float vCore;",
      "varying float vGlow;",
      "void main() {",
      "  vColor = starColor;",
      "  vCore = starCore / starSize;",
      "  vGlow = starGlow;",
      "  gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);",
      "  gl_PointSize = starSize * pixelRatio;",
      "  #include <logdepthbuf_vertex>",
      "}",
    ].join("\n"),
    fragmentShader: [
      "#include <common>",
      "#include <logdepthbuf_pars_fragment>",
      "varying vec3 vColor;",
      "varying float vCore;",
      "varying float vGlow;",
      "void main() {",
      "  #include <logdepthbuf_fragment>",
      "  float r = length(gl_PointCoord * 2.0 - 1.0);",
      "  if (r > 1.0) discard;",
      "  float core = 1.0 - smoothstep(vCore * 0.5, vCore, r);",
      "  float halo = vGlow * exp(-r * r * 4.0) * (1.0 - r);",
      "  gl_FragColor = vec4(mix(vColor, vec3(1.0), core * 0.75), clamp(core + halo, 0.0, 1.0));",
      "}",
    ].join("\n"),
    transparent: true,
    depthWrite: false,
  });
  var starPoints = new THREE.Points(new THREE.BufferGeometry(), starMaterial);
  starPoints.renderOrder = 5;
  starPoints.frustumCulled = false;
  scene.add(starPoints);
  var starList = [];

  // Draws exactly `stars` (one entry per id).
  function setStars(stars) {
    starList = stars;
    var n = stars.length;
    var positions = new Float32Array(3 * n);
    var colors = new Float32Array(3 * n);
    var sizes = new Float32Array(n);
    var cores = new Float32Array(n);
    var glows = new Float32Array(n);
    stars.forEach(function (star, i) {
      var t = THREE.MathUtils.clamp(
        (Math.log10(Math.max(star.luminosity_sol, 1)) - STAR_LOG_LUMINOSITY[0]) / (STAR_LOG_LUMINOSITY[1] - STAR_LOG_LUMINOSITY[0]), 0, 1);
      positions.set([star.x, star.y, star.z], 3 * i);
      colors.set(starColor(star.temperature_k), 3 * i);
      sizes[i] = THREE.MathUtils.lerp(STAR_MIN_PX, STAR_MAX_PX, t);
      cores[i] = THREE.MathUtils.lerp(STAR_CORE_PX[0], STAR_CORE_PX[1], t);
      glows[i] = THREE.MathUtils.lerp(STAR_GLOW[0], STAR_GLOW[1], t);
    });
    var geometry = new THREE.BufferGeometry();
    geometry.setAttribute("position", new THREE.BufferAttribute(positions, 3));
    geometry.setAttribute("starColor", new THREE.BufferAttribute(colors, 3));
    geometry.setAttribute("starSize", new THREE.BufferAttribute(sizes, 1));
    geometry.setAttribute("starCore", new THREE.BufferAttribute(cores, 1));
    geometry.setAttribute("starGlow", new THREE.BufferAttribute(glows, 1));
    starPoints.geometry.dispose();
    starPoints.geometry = geometry;
  }

  // The star whose center is nearest a screen point, within
  // STAR_PICK_PX, or null.
  function starAtClientPoint(clientX, clientY) {
    var rect = canvasEl.getBoundingClientRect();
    if (!starList.length || !rect.width || !rect.height) {
      return null;
    }
    var projected = new THREE.Vector3();
    var best = null;
    var bestPx = STAR_PICK_PX;
    starList.forEach(function (star) {
      projected.set(star.x, star.y, star.z).project(camera);
      if (projected.z < -1 || projected.z > 1) {
        return;
      }
      var px = Math.hypot(
        rect.left + ((projected.x + 1) / 2) * rect.width - clientX,
        rect.top + ((1 - projected.y) / 2) * rect.height - clientY);
      if (px <= bestPx) {
        best = star;
        bestPx = px;
      }
    });
    return best;
  }

  // --- Drawing from tiles --------------------------------------------------

  // Filled points per tile object (a refetched tile is a new object).
  var filledPointsByTile = new WeakMap();
  var filledSignature = "";
  var starSignature = "";

  // Takes the filled sectors from whatever of the needed tiles is already
  // cached (the blocks are rebuilt from them when they change); returns
  // what's still missing.
  function renderFromCache(need) {
    fadeRadius = need.viewRadius;
    var points = [];
    var present = [];
    var missing = [];
    var clouds = new Map();
    var stars = new Map();
    need.keys.forEach(function (key) {
      var tile = getTile(key);
      if (tile === undefined) {
        missing.push(key);
        return;
      }
      var tilePoints = filledPointsByTile.get(tile);
      if (!tilePoints) {
        tilePoints = filledPointsOf(tile.filled);
        filledPointsByTile.set(tile, tilePoints);
      }
      present.push(key + "=" + tilePoints.length);
      Array.prototype.push.apply(points, tilePoints);
      (tile.clouds || []).forEach(function (cloud) {
        clouds.set(cloudKey(cloud), cloud);
      });
      (tile.stars || []).forEach(function (star) {
        stars.set(star.id, star);
      });
    });
    var starKeys = currentStamp + "|" + Array.from(stars.keys()).sort().join(",");
    if (starKeys !== starSignature) {
      starSignature = starKeys;
      setStars(Array.from(stars.values()));
    }
    var cloudKeys = Array.from(clouds.keys()).sort().join(",");
    if (cloudKeys !== cloudSignature) {
      cloudSignature = cloudKeys;
      setClouds(Array.from(clouds.values()));
    }
    var signature = currentStamp + "|" + present.join(",");
    if (signature !== filledSignature) {
      filledSignature = signature;
      filledPoints = points;
      filledVersion++;
    }
    return missing;
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
  // zoomInOnTarget) is never throttled here, only the network fetch that
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
    // Mid-glide, the tiles where the zoom is headed too.
    if (zoomAnimation) {
      neededTiles(zoomAnimation.to).keys.forEach(function (key) {
        if (missing.indexOf(key) < 0 && getTile(key) === undefined) {
          missing.push(key);
        }
      });
    }
    if (!missing.length) {
      prefetchTiles();
      return;
    }
    if (activeAbort) {
      activeAbort.abort();
    }
    var controller = typeof AbortController !== "undefined" ? new AbortController() : null;
    activeAbort = controller;
    var params = new URLSearchParams({
      tiles: missing.slice(0, MAX_TILES_PER_REQUEST).join(","),
    });
    if (currentStamp) {
      params.set("stamp", currentStamp);
    }
    fetch(data.fetchPath + "?" + params.toString(), controller ? { signal: controller.signal } : undefined)
      .then(function (response) {
        if (!response.ok) {
          throw new Error(data.fetchPath + " returned " + response.status);
        }
        return response.json();
      })
      .then(function (payload) {
        if (activeAbort === controller) {
          activeAbort = null;
        }
        var stampChanged = absorb(payload);
        var stillMissing = renderFromCache(neededTiles());
        // A new stamp dropped changed (or all) cached tiles, and a view
        // needing more than one request's worth of tiles has more to
        // fetch -- either way, go again (tiles already fetched are
        // cached by now).
        if (stampChanged || stillMissing.length) {
          scheduleFetch(false);
        } else {
          prefetchTiles();
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

  // The tiles one click-zoom step in and out, fetched at low priority once
  // the view's own are in. One prefetch runs at a time and a real fetch
  // never waits on or cancels it (or the other way round); its tiles just
  // land in the cache.
  var prefetching = false;

  function prefetchTiles() {
    if (prefetching || activeAbort) {
      return;
    }
    var factor = clickZoomFactor(orbit.radius);
    var missing = [];
    [orbit.radius / factor, orbit.radius * factor].forEach(function (radius) {
      neededTiles(Math.max(MIN_RADIUS, Math.min(MAX_RADIUS, radius))).keys.forEach(function (key) {
        if (missing.indexOf(key) < 0 && getTile(key) === undefined) {
          missing.push(key);
        }
      });
    });
    if (!missing.length) {
      return;
    }
    var params = new URLSearchParams({ tiles: missing.slice(0, MAX_TILES_PER_REQUEST).join(",") });
    if (currentStamp) {
      params.set("stamp", currentStamp);
    }
    prefetching = true;
    fetch(data.fetchPath + "?" + params.toString(), { priority: "low" })
      .then(function (response) {
        if (!response.ok) {
          throw new Error(data.fetchPath + " returned " + response.status);
        }
        return response.json();
      })
      .then(function (payload) {
        prefetching = false;
        // A new stamp means the view's own tiles may be stale too.
        if (absorb(payload)) {
          scheduleFetch(false);
        } else if (missing.length > MAX_TILES_PER_REQUEST) {
          prefetchTiles();
        }
      })
      .catch(function () {
        // Only a head start: the real fetch gets these when they're needed.
        prefetching = false;
      });
  }

  // Normally every tile is already cached (the page embeds them), so this
  // only fetches if the page's own tile fetch came up short.
  // First frame's tiles, then (normally nothing is missing) the prefetch.
  renderFromCache(neededTiles());
  scheduleFetch(false);

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
    if (!dragging) {
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

  // Wheel, button and double-click zooms glide to the new radius over
  // ZOOM_MS (evenly in log space, easing out) instead of jumping; another
  // step mid-glide starts from where the camera is, toward the new goal.
  // With prefers-reduced-motion they jump.
  var ZOOM_MS = 160;
  var zoomAnimation = null;

  // Where the camera is headed: the end of the glide, or where it is.
  function zoomGoal() {
    return zoomAnimation ? zoomAnimation.to : orbit.radius;
  }

  function zoomTo(radius) {
    radius = Math.max(MIN_RADIUS, Math.min(MAX_RADIUS, radius));
    if (reducedMotion) {
      zoomAnimation = null;
      setRadius(radius);
      return;
    }
    zoomAnimation = { from: orbit.radius, to: radius, startedAt: performance.now() };
  }

  function stepZoom(now) {
    if (!zoomAnimation) {
      return;
    }
    var t = Math.min(1, (now - zoomAnimation.startedAt) / ZOOM_MS);
    var eased = 1 - Math.pow(1 - t, 3);
    setRadius(zoomAnimation.from * Math.pow(zoomAnimation.to / zoomAnimation.from, eased));
    if (t >= 1) {
      zoomAnimation = null;
    }
  }

  canvasEl.addEventListener(
    "wheel",
    function (event) {
      event.preventDefault();
      var deltaPx = event.deltaY;
      if (event.deltaMode === 1) {
        deltaPx *= 33;
      } else if (event.deltaMode === 2) {
        deltaPx *= canvasEl.clientHeight || 400;
      }
      if (event.ctrlKey) {
        deltaPx *= PINCH_BOOST;
      }
      deltaPx = Math.max(-WHEEL_MAX_PX, Math.min(WHEEL_MAX_PX, deltaPx));
      if (!deltaPx) {
        return;
      }
      zoomTo(zoomGoal() * Math.exp(deltaPx * WHEEL_ZOOM_PER_PX));
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

  // The block under a screen point, as showCellInfo's cell, plus the point
  // to recenter on (the mean of its filled sectors, or its middle) -- or
  // null. Along the ray, the nearest block holding filled sectors wins over
  // the translucent space in front of it, so clicking finds what's
  // generated. Blocks hidden by the near cut (NEAR_CUT) are skipped.
  function cellAtClientPoint(clientX, clientY) {
    var meshes = shownSet ? shownSet.meshes.filter(function (mesh) { return mesh.visible; }) : [];
    if (!meshes.length) {
      return null;
    }
    var ndc = ndcFromClientPoint(clientX, clientY);
    if (!ndc) {
      return null;
    }
    raycaster.setFromCamera(ndc, camera);
    var hits = raycaster.intersectObjects(meshes, false);
    var nearCut = orbit.radius * NEAR_CUT;
    var first = null;
    for (var i = 0; i < hits.length; i++) {
      var face = hits[i].face;
      if (!face) {
        continue;
      }
      var mesh = hits[i].object;
      var cellData = blockRecord(mesh.userData.set, mesh.userData.part, mesh.userData.part.owners[face.a]);
      if (!cellData) {
        continue;
      }
      var coords = cellCoordinates(cellData).cartesian;
      if (new THREE.Vector3(coords[0], coords[1], coords[2]).distanceTo(camera.position) < nearCut) {
        continue;
      }
      if (cellData.filled > 0) {
        return blockHit(cellData);
      }
      if (!first) {
        first = cellData;
      }
    }
    return first ? blockHit(first) : null;
  }

  // Block `index` of one mesh of a block set, the way the listing gives
  // blocks (bounds, density, meanDensity) plus its m, filled, total,
  // centroid and generated sector -- or null.
  function blockRecord(set, part, index) {
    var rec = part.cells;
    var o = index * CELL_STRIDE;
    if (!(o + CELL_STRIDE <= rec.length)) {
      return null;
    }
    var block = blockAt(rec[o], rec[o + 1], rec[o + 2], set.m, edgePc, null, null);
    block.m = set.m;
    block.density = isNaN(rec[o + 3]) ? null : rec[o + 3];
    block.meanDensity = isNaN(rec[o + 4]) ? null : rec[o + 4];
    block.filled = rec[o + 5];
    block.total = rec[o + 6];
    if (!isNaN(rec[o + 7])) {
      block.centroid = [rec[o + 7], rec[o + 8], rec[o + 9]];
    }
    if (!isNaN(rec[o + 10])) {
      var point = set.points[rec[o + 10]];
      block.sector = point && point.sector ? point.sector : null;
    }
    return block;
  }

  function blockHit(cellData) {
    var m = cellData.m;
    var at = cellData.centroid || cellCoordinates(cellData).cartesian;
    var cell = {
      m: m, bounds: cellData, density: cellData.density, ring: cellData.ring, seg: cellData.seg, slab: cellData.slab,
      edgePc: edgePc, shape: galaxyShape, filled: cellData.filled || 0, sector: cellData.sector || null,
    };
    if (m === 1) {
      cell.address = { ring: cellData.ring, layer: cellData.slab, slot: cellData.seg };
    }
    return { point: new THREE.Vector3(at[0], at[1], at[2]), cell: cell };
  }

  // The sector cell holding a galaxy-frame point -- every point in space
  // has one, whether or not anything was ever generated there.
  function sectorCellAt(point) {
    var a = sectorAddressAt(point.x, point.y, point.z, edgePc);
    return { m: 1, address: a, bounds: sectorCellBounds(a.ring, a.layer, a.slot, edgePc), density: null };
  }

  // Empty-space click target: where the click's ray meets the plane
  // through the current target parallel to the galactic disk -- so
  // clicking empty space over a spiral arm lands ON the arm, the point
  // you see under the cursor, rather than above or below it. When the
  // disk is seen nearly edge-on that ray can run off almost parallel to
  // the plane, so a hit farther than DISK_PICK_MAX_RADII orbit radii from
  // the target falls back to the plane through the target facing the
  // camera (the depth the camera is looking at). Either way the point is
  // kept inside the galaxy -- within GALAXY_RADIUS of the core across the
  // disk and a tenth of that above or below it -- so a click can never
  // recenter the view out in the void.
  //
  // (This used to intersect a sphere of the orbit radius around the
  // target -- but the camera sits ON that sphere, so the ray's nearest
  // hit was the camera itself or the sphere's far side, and a click
  // moved the view thousands of parsecs away from where it landed.)
  var DISK_PICK_MAX_RADII = 2;
  var GALAXY_HALF_THICKNESS = GALAXY_RADIUS / 10;

  function depthPointAtClientPoint(clientX, clientY) {
    var ndc = ndcFromClientPoint(clientX, clientY);
    if (!ndc) {
      return null;
    }
    raycaster.setFromCamera(ndc, camera);
    var ray = raycaster.ray;
    var hitPoint = new THREE.Vector3();
    var diskPlane = new THREE.Plane(new THREE.Vector3(0, 0, 1), -target.z);
    var hit = ray.intersectPlane(diskPlane, hitPoint);
    if (!hit || hitPoint.distanceTo(target) > orbit.radius * DISK_PICK_MAX_RADII) {
      var facing = new THREE.Vector3().subVectors(camera.position, target).normalize();
      var focalPlane = new THREE.Plane().setFromNormalAndCoplanarPoint(facing, target);
      hit = ray.intersectPlane(focalPlane, hitPoint);
    }
    if (!hit) {
      return null;
    }
    var across = Math.hypot(hitPoint.x, hitPoint.y);
    if (across > GALAXY_RADIUS) {
      hitPoint.x *= GALAXY_RADIUS / across;
      hitPoint.y *= GALAXY_RADIUS / across;
    }
    var maxZ = Math.max(GALAXY_HALF_THICKNESS, Math.abs(target.z));
    hitPoint.z = Math.max(-maxZ, Math.min(maxZ, hitPoint.z));
    return hitPoint;
  }

  // Shows a block's info and rings it: a generated sector's own panel (with
  // its link) at one sector per block, otherwise the block's.
  function selectCell(point, cell, cloud, star) {
    if (star) {
      highlightPosition(point.x, point.y, point.z);
      showStarInfo(star);
      return;
    }
    if (cloud) {
      highlightPosition(point.x, point.y, point.z);
      showCloudInfo(cloud);
      return;
    }
    if (!cell) {
      return;
    }
    highlightPosition(point.x, point.y, point.z);
    if (cell.sector) {
      showPlacedInfo(cell.sector);
    } else {
      showCellInfo(cell);
    }
  }

  // Single click: re-centers the view on the clicked block (and shows its
  // info) WITHOUT zooming -- deliberately not
  // "click to zoom" any more (see this file's own module docstring's
  // former "Click-to-zoom is LOGARITHMIC" note, now double-click's own
  // job below). Centering alone, with no zoom commitment, is what makes
  // it possible to walk the camera across the galaxy toward a small/
  // distant block over several clicks without a bad click also zooming
  // into empty space you didn't mean to approach.
  function centerOn(point, cell, cloud, star) {
    target.copy(point);
    applyCamera();
    updateScaleBar();
    selectCell(point, cell, cloud, star);
    scheduleFetch(true);
  }

  // Double click: centers AND zooms in by one clickZoomFactor step. A
  // double-click's first click has already centered on the point (see the
  // click handler), so this zooms in on that same point rather than
  // re-picking under the cursor, which by now is over something else.
  function zoomInOnTarget() {
    zoomTo(zoomGoal() / clickZoomFactor(zoomGoal()));
    scheduleFetch(true);
  }

  // Resolves a click/double-click's target point the same way for both:
  // a bright star clicked on (it is only a few pixels wide, so a click
  // that close means it), else a cloud small enough to aim at, clicked
  // near its middle (its center),
  // else a hit block holding generated sectors (their mean position, so a
  // double-click zooms in toward them), else such a cloud clicked
  // anywhere, else any hit block, or (empty space)
  // depthPointAtClientPoint -- shared so the click handlers below never
  // have to duplicate the raycast-then-fall-back logic.
  function resolveClickTarget(clientX, clientY) {
    var star = starAtClientPoint(clientX, clientY);
    if (star) {
      return { point: new THREE.Vector3(star.x, star.y, star.z), cell: null, cloud: null, star: star };
    }
    var hit = cellAtClientPoint(clientX, clientY);
    var found = cloudAtClientPoint(clientX, clientY);
    if (found && (found.core || !(hit && hit.cell.filled > 0))) {
      var cloud = found.cloud;
      return { point: new THREE.Vector3(cloud.x, cloud.y, cloud.z), cell: null, cloud: cloud };
    }
    if (hit) {
      return hit;
    }
    var depthPoint = depthPointAtClientPoint(clientX, clientY);
    return depthPoint ? { point: depthPoint, cell: galaxyShape ? sectorCellAt(depthPoint) : null } : null;
  }

  // Only a double-click's FIRST click recenters (event.detail counts the
  // clicks in a burst): the second one lands at the same screen spot,
  // which after the first recenter shows a different point, and
  // recentering again there walked the view off in that direction.
  var lastCenterClickAt = 0;

  canvasEl.addEventListener("click", function (event) {
    if (suppressNextClick) {
      suppressNextClick = false;
      return;
    }
    if (event.detail > 1) {
      return;
    }
    var resolved = resolveClickTarget(event.clientX, event.clientY);
    if (resolved) {
      centerOn(resolved.point, resolved.cell, resolved.cloud, resolved.star);
      lastCenterClickAt = Date.now();
    }
  });

  canvasEl.addEventListener("dblclick", function (event) {
    if (suppressNextClick) {
      suppressNextClick = false;
      return;
    }
    // Normally the first click just centered the view (see above); if it
    // didn't (a browser that skipped it), center here first.
    if (Date.now() - lastCenterClickAt > 1000) {
      var resolved = resolveClickTarget(event.clientX, event.clientY);
      if (resolved) {
        centerOn(resolved.point, resolved.cell, resolved.cloud, resolved.star);
      }
    }
    zoomInOnTarget();
  });

  // No right-click action any more (see this file's own module
  // docstring) -- the browser's own default context menu is left alone,
  // rather than preventDefault()-ing it for nothing.

  var controlsEl = document.getElementById("galaxymap3d-controls");
  if (controlsEl) {
    // Five buttons: let them wrap rather than run off the side panel.
    controlsEl.style.flexWrap = "wrap";
    controlsEl.querySelectorAll("[data-action]").forEach(function (button) {
      button.addEventListener("click", function () {
        var action = button.dataset.action;
        if (action === "zoom-in") {
          zoomTo(zoomGoal() / clickZoomFactor(zoomGoal()));
          scheduleFetch(true);
        } else if (action === "zoom-out") {
          zoomTo(zoomGoal() * clickZoomFactor(zoomGoal()));
          scheduleFetch(true);
        } else if (action === "reset") {
          resetView();
        } else if (action === "wedges") {
          wedgeGroup.visible = !wedgeGroup.visible;
          button.setAttribute("aria-pressed", String(wedgeGroup.visible));
        } else if (action === "slice") {
          // The next frame picks the block set for the new slice.
          sliceAtFocus = !sliceAtFocus;
          button.setAttribute("aria-pressed", String(sliceAtFocus));
        }
      });
    });
  }

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

  // Up to 3 significant figures, grouped: 0.0512, 3.4, 1,280.
  function formatCount(value) {
    if (!(value > 0)) return "0";
    if (value >= 100) return Math.round(value).toLocaleString();
    return String(Number(value.toPrecision(value >= 1 ? 3 : 2)));
  }

  function plural(count, word) {
    return count === 1 ? word : word + "s";
  }

  // One line of the readout: "<lead> <sectors> · <distance>", the distance
  // on the shared ladder (static/distance.js), which adds ly in
  // parentheses to parsec values.
  function scaleLine(lead, pc, suffix) {
    var sectors = pc / edgePc;
    var line = document.createElement("span");
    line.className = "starmap-scale-label";
    // Wraps on a phone instead of running off the canvas.
    line.style.whiteSpace = "normal";
    line.textContent = lead + " " + formatCount(sectors) + " " + plural(sectors, "sector") + (suffix || "")
      + " · " + formatDistancePc(pc);
    return line;
  }

  // Three lines: what one screen pixel spans at the focus, how big one
  // prism (block) is, and a bar of about SCALE_BAR_PX rounded to a nice
  // number of parsecs. All per CSS pixel, like the block size itself:
  // they're about what a person can see and click, not device pixels.
  var SCALE_BAR_PX = 70;
  if (scaleEl) {
    // A column instead of the other maps' one-line row (set here, not in
    // style.css: only this map's readout stacks).
    scaleEl.style.flexDirection = "column";
    scaleEl.style.alignItems = "flex-start";
    scaleEl.style.gap = "0.1rem";
    scaleEl.style.maxWidth = "calc(100% - 1.5rem)";
  }

  function updateScaleBar() {
    if (!scaleEl) {
      return;
    }
    var pcPerScreenPx = pcPerPixelAtTarget();
    var nicePc = niceScaleValue(SCALE_BAR_PX * pcPerScreenPx);
    if (!nicePc) {
      return;
    }
    var m = drawnSectorsPerPrism || 1;
    var block = scaleLine("1 block =", m * edgePc, m === 1 ? "" : " across (≈ " + formatCount(m * m * m) + ")");
    var barRow = document.createElement("span");
    barRow.style.display = "flex";
    barRow.style.alignItems = "center";
    barRow.style.gap = "0.4rem";
    var bar = document.createElement("span");
    bar.className = "starmap-scale-bar";
    bar.style.width = (nicePc / pcPerScreenPx).toFixed(1) + "px";
    barRow.appendChild(bar);
    barRow.appendChild(scaleLine("≈", nicePc));
    scaleEl.replaceChildren(scaleLine("1 px ≈", pcPerScreenPx), block, barRow);
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

  // Every frame: the zoom glide, then (when the view moved) the filled
  // sectors from cached tiles, a fetch for missing ones and the wedge
  // lines, then the block set for the view.
  var lastViewSignature = "";

  (function animate(now) {
    requestAnimationFrame(animate);
    stepZoom(now || performance.now());
    var viewSignature = [target.x, target.y, target.z, orbit.radius, canvasEl.clientHeight].join(",");
    if (viewSignature !== lastViewSignature) {
      lastViewSignature = viewSignature;
      var missing = renderFromCache(neededTiles());
      if (missing.length && !fetchTimer && !activeAbort) {
        scheduleFetch(false);
      }
      updateWedgeLevels(pcPerPixelAtTarget());
    }
    syncBlocks();
    updateClouds();
    updateHighlightScale();
    updateWedgeLabels();
    renderer.render(scene, camera);
  })();
}

if (canvas && sceneData) {
  initGalaxyMap3d(canvas, sceneData);
}
