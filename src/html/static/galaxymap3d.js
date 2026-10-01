// html/static/galaxymap3d.js
// TODO(galaxy-map #18): the renderer stays three.js (vendored r186).
// Babylon.js and deck.gl need a bundler or ship several MB; regl or raw
// WebGPU would mean rewriting picking, labels and lighting. The slow part
// is the JavaScript block listing, not drawing. When the map outgrows
// WebGL, try three's own WebGPURenderer first.
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
  blockAddressAt, blockAt, blockSectorCount, blockSectorRanges, blockSizeForScale, blockSlotRange, blocksForView,
  buildPrismGeometry, cellCoordinates, cellVertices, sectorAddressAt, sectorCellBounds, wedgeLines,
} = await import(`./galaxyprisms.js${VERSION_QUERY}`);
const { formatDistancePc } = await import(`./distance.js${VERSION_QUERY}`);

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

// TODO(sector-map #24): the Galaxy Map's unfilled-sector panel gets the
// same admin generate buttons as the sector map in place of this snippet.
function cliSnippet(ringIndex, layerIndex, slotIndex) {
  return "generate.py galaxy --ring " + ringIndex + " --layer " + layerIndex + " --slot " + slotIndex;
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

  if (single) {
    var snippet = cliSnippet(cell.address.ring, cell.address.layer, cell.address.slot);
    var code = document.createElement("code");
    code.className = "galaxymap3d-cli-snippet";
    code.textContent = snippet;
    panel.appendChild(code);
    panel.appendChild(makeCopyButton(snippet));
  }
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
// tracking existed) rather than a false 0, so placedDensityColor below
// can fall back to a neutral mid-tone instead of reading as "empty".
var PLACED_LOW_DENSITY_COLOR = "#4a3f2e";
var PLACED_HIGH_DENSITY_COLOR = "#fff6df";

function placedRelativeDensity(entry, referenceDensityPerLy3) {
  if (!entry.edge_ly || !referenceDensityPerLy3) {
    return null;
  }
  var densityPerLy3 = (entry.system_count || 0) / Math.pow(entry.edge_ly, 3);
  return densityPerLy3 / referenceDensityPerLy3;
}

// Colors a filled sector's own block at one sector per block: its real
// density on a "log2(x+1)/3" dim-to-bright curve, through a warm
// bronze-to-gold range (rather than the density shading's cooler
// dim-to-accent one) so a generated sector stands out from predicted
// space at a glance. Coarser blocks take the space's density color
// (prismShade), and their filled share only sets their opacity.
function placedDensityColor(entry, referenceDensityPerLy3) {
  var relative = placedRelativeDensity(entry, referenceDensityPerLy3);
  var t = relative == null ? 0.5 : Math.max(0, Math.min(1, Math.log2(relative + 1) / 3));
  var low = new THREE.Color(PLACED_LOW_DENSITY_COLOR);
  var high = new THREE.Color(PLACED_HIGH_DENSITY_COLOR);
  return low.clone().lerp(high, t);
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
  // Everything is drawn as one solid of blocks (galaxyprisms.js), rebuilt
  // whenever the view's block set changes (updatePrisms, below). Blocks
  // fill their whole cell, so the galaxy reads as one solid; a thin
  // brighter outline along each face's own edges (found from the face's
  // 0..1 uv, about a screen pixel wide) keeps neighboring blocks apart.
  //
  // Two meshes share one shader:
  // - solidMesh: blocks whose every sector is generated, opaque;
  // - glassMesh: everything else, translucent (no depth writes), its
  //   blocks sorted back to front from the camera on every rebuild. Unfilled
  //   space is BLOCK_OPACITY_SPARSE (sparsest) to BLOCK_OPACITY_DENSE
  //   (densest) opaque, and a block grows more solid with its filled share
  //   (blockOpacity).
  // Blocks holding filled sectors also get a warm tint on their faces and
  // warm, stronger face edges, in step with their filled share
  // (blockFillStep), so they can be picked out even where the opacity
  // step alone is small.
  // Whole blocks close to the camera are dropped (nearCut, tested on each
  // block's own center), so flying through the disk shows what's ahead
  // instead of a wall of the nearest blocks. The logdepthbuf chunks match
  // the renderer's logarithmic depth buffer.
  var FILLED_TINT = new THREE.Color(PLACED_HIGH_DENSITY_COLOR);

  function makeBlockMaterial(translucent) {
    return new THREE.ShaderMaterial({
      uniforms: { nearCut: { value: 0 }, filledTint: { value: FILLED_TINT } },
      transparent: translucent,
      depthWrite: !translucent,
      vertexShader: [
        "#include <common>",
        "#include <logdepthbuf_pars_vertex>",
        "uniform float nearCut;",
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
        "  vColor = prismColor;",
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
        "  gl_FragColor = vec4(mix(face, edgeColor, (0.18 + 0.7 * vFill) * edge), vAlpha);",
        "  #include <colorspace_fragment>",
        "}",
      ].join("\n"),
    });
  }
  var solidMaterial = makeBlockMaterial(false);
  var glassMaterial = makeBlockMaterial(true);
  var solidMesh = new THREE.Mesh(new THREE.BufferGeometry(), solidMaterial);
  var glassMesh = new THREE.Mesh(new THREE.BufferGeometry(), glassMaterial);
  [solidMesh, glassMesh].forEach(function (mesh) {
    mesh.frustumCulled = false;
    mesh.visible = false;
    mesh.userData.cells = [];
    mesh.userData.owners = null;
  });
  // Opaque first, then the translucent shell over it.
  solidMesh.renderOrder = 0;
  glassMesh.renderOrder = 1;
  scene.add(solidMesh);
  scene.add(glassMesh);

  var highlightSprite = new THREE.Sprite(
    new THREE.SpriteMaterial({ map: highlightTexture, transparent: true, depthWrite: false, depthTest: false })
  );
  highlightSprite.visible = false;
  highlightSprite.renderOrder = 3;
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

  // Blocks darken toward the edge of the view ball around the target, so
  // the solid fades out instead of stopping at a hard spherical rim the eye
  // reads as a ball. fadeRadius is the current view radius (set by
  // renderFromCache).
  var EDGE_FADE_START = 0.6;
  var EDGE_FADE_END = 1.0;
  var fadeRadius = 0;

  function edgeFade(x, y, z) {
    if (!(fadeRadius > 0)) {
      return 1;
    }
    var dx = x - target.x;
    var dy = y - target.y;
    var dz = z - target.z;
    var t = Math.sqrt(dx * dx + dy * dy + dz * dz) / fadeRadius;
    var u = Math.max(0, Math.min(1, (t - EDGE_FADE_START) / (EDGE_FADE_END - EDGE_FADE_START)));
    return 1 - u * u * (3 - 2 * u);
  }

  // --- Filled sectors --------------------------------------------------------
  //
  // Every tile carries a `filled` summary (queryDb.galaxy_filled_in_box):
  // at g = 1 each generated sector's address, id, name and system count;
  // coarser, counts per cell g sectors a side. renderFromCache turns them
  // into points (a sector's center, or a cell's middle) with counts, and
  // updatePrisms sums those into the blocks it draws. A cell never spans
  // more than one block ring or layer (g divides m), so only a cell's
  // wedge can straddle two blocks, and its count goes to the one holding
  // its middle.
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

  // Unfilled space: this opaque at the sparsest drawn density, rising to
  // BLOCK_OPACITY_DENSE at the densest (50% to 20% see-through).
  var BLOCK_OPACITY_SPARSE = 0.5;
  var BLOCK_OPACITY_DENSE = 0.8;
  // Any filled sector lifts a block at least this far from its unfilled
  // opacity toward solid, however tiny its share: at large m one filled
  // sector among 531,441 must still show. The rest of the way follows the
  // log of the filled count against the log of the block's total, so a
  // block turns solid only when every sector in it is filled.
  var FILLED_MIN_STEP = 0.3;

  // 0 for a block with nothing filled, FILLED_MIN_STEP for one filled
  // sector among many, 1 once every sector in it is filled.
  function blockFillStep(cell) {
    if (!(cell.filled > 0)) {
      return 0;
    }
    if (cell.filled >= cell.total) {
      return 1;
    }
    var share = Math.log(1 + cell.filled) / Math.log(1 + cell.total);
    return FILLED_MIN_STEP + (1 - FILLED_MIN_STEP) * share;
  }

  function blockOpacity(cell) {
    var base = BLOCK_OPACITY_SPARSE + (BLOCK_OPACITY_DENSE - BLOCK_OPACITY_SPARSE) * prismIntensity(cell.density || 0);
    return base + (1 - base) * blockFillStep(cell);
  }

  // The light the blocks' faces are shaded by: from galactic north, a
  // little off to one side, so tops, walls and sides all read apart.
  var PRISM_LIGHT = new THREE.Vector3(0.35, -0.3, 0.9).normalize();
  var PRISM_AMBIENT = 0.35;
  // With the slice on (the default), blocks above the focus's own layer
  // are left out, so the view looks down on the solid's cut face -- the
  // arms and whatever layer the focus is in -- instead of its outside.
  var sliceAtFocus = true;
  var lastPrismViewRadius = null;
  // Block shading runs over a wide density range: the drawing floor
  // (galaxyprisms.js's PRISM_MIN_DENSITY) up to the core, on a log scale,
  // through a dim-to-accent-to-white ramp.
  var PRISM_DENSITY_LOW = 0.02;
  var PRISM_DENSITY_HIGH = 100;
  var PRISM_DIM = new THREE.Color(0x1d2340);
  var PRISM_ACCENT = new THREE.Color(accentColor);
  var PRISM_HOT = new THREE.Color(0xeef0ff);

  // Density alone puts the arms (1 +/- arm_amplitude around the ring's
  // mean) on a sliver of that 0.02-100 ramp, so the spiral barely shows.
  // Each block's shade mixes its density with its arm factor (density /
  // its azimuthal mean, from galaxyprisms.js), stretched over the arm
  // model's own range: inter-arm troughs dim, arm crests bright.
  var PRISM_DENSITY_SHARE = 0.45;
  var ARM_AMPLITUDE = galaxyArmAmplitude(data.densityShape);

  function galaxyArmAmplitude(shape) {
    if (!shape || !(shape.arm_count > 0)) {
      return 0;
    }
    return Math.min(1, Math.abs(shape.arm_amplitude || 0));
  }

  // 0..1 shade for a block: density alone when the galaxy has no arms.
  function prismShade(cell) {
    var t = prismIntensity(cell.density);
    if (!(ARM_AMPLITUDE > 0) || !(cell.meanDensity > 0)) {
      return t;
    }
    var arm = (cell.density / cell.meanDensity - (1 - ARM_AMPLITUDE)) / (2 * ARM_AMPLITUDE);
    return PRISM_DENSITY_SHARE * t + (1 - PRISM_DENSITY_SHARE) * Math.max(0, Math.min(1, arm));
  }

  function prismIntensity(relativeDensity) {
    var t = Math.log(Math.max(relativeDensity, 1e-9) / PRISM_DENSITY_LOW) / Math.log(PRISM_DENSITY_HIGH / PRISM_DENSITY_LOW);
    return Math.max(0, Math.min(1, t));
  }

  function prismColor(t) {
    return t < 0.6 ? PRISM_DIM.clone().lerp(PRISM_ACCENT, t / 0.6) : PRISM_ACCENT.clone().lerp(PRISM_HOT, (t - 0.6) / 0.4);
  }
  // Blocks centered nearer the camera than this x the orbit radius (the
  // camera's distance to the target) aren't drawn.
  var NEAR_CUT = 0.4;
  var galaxyShape = data.densityShape || null;
  var edgePc = data.edgePc || 1;
  // Densities already worked out, per block size.
  var prismDensityCaches = new Map();
  var PRISM_CACHE_SIZES = 4;
  var prismSetKey = "";
  // How many sectors a side the drawn blocks are.
  var drawnSectorsPerPrism = 1;

  function prismDensityCache(m) {
    var cache = memoryGet(prismDensityCaches, m);
    if (!cache) {
      cache = new Map();
      memorySet(prismDensityCaches, m, cache, PRISM_CACHE_SIZES);
    }
    return cache;
  }

  // TODO(galaxy-map #19): this is per CSS pixel. The block size (#13)
  // should stay per CSS pixel (it is about what a person can see and
  // click), but say so in the readout's code.
  // Parsecs per screen pixel at the view's focus (the orbit target).
  function pcPerPixelAtTarget() {
    var heightPx = canvasEl.clientHeight || 1;
    return (orbit.radius * 2 * Math.tan(THREE.MathUtils.degToRad(camera.fov) / 2)) / heightPx;
  }

  // blockSectorCount per filled block (m/ring/seg/slab): at large m it
  // walks every member ring and layer, so it's kept between rebuilds.
  var blockTotals = new Map();
  var BLOCK_TOTALS_MAX = 20000;

  // The blocks to draw: the solid's surface (galaxyprisms.blocksForView),
  // plus every block in view holding filled sectors, interior ones
  // included (they show through the translucent shell). Each gets `filled`,
  // and filled blocks also `total` (blockSectorCount), `centroid` (the
  // mean position of their filled sectors, what a click centers on) and,
  // at one sector per block, `sector`. Without a galaxy shape only filled
  // blocks are drawn.
  function viewBlocks(center, viewRadius, pcPerPixel) {
    var options = {
      pcPerPixel: pcPerPixel, sliceZ: sliceAtFocus ? target.z : null, viewerZ: camera.position.z,
      minPx: data.blockMinPx, budget: data.blockBudget,
    };
    var view = galaxyShape
      ? blocksForView(center, viewRadius, edgePc, GALAXY_RADIUS, galaxyShape, prismDensityCache, options)
      : { m: blockSizeForScale(pcPerPixel, edgePc, data.blockMinPx, GALAXY_RADIUS), slice: null, blocks: [] };
    var m = view.m;
    var byKey = new Map();
    view.blocks.forEach(function (block) {
      byKey.set(block.ring + "/" + block.seg + "/" + block.slab, block);
    });
    var sums = new Map();
    filledPoints.forEach(function (point) {
      var a = blockAddressAt(point.x, point.y, point.z, m, edgePc);
      if (view.slice !== null && a.slab > view.slice) {
        return;
      }
      var key = a.ring + "/" + a.seg + "/" + a.slab;
      var sum = sums.get(key);
      if (!sum) {
        sum = { address: a, count: 0, x: 0, y: 0, z: 0, sector: null };
        sums.set(key, sum);
      }
      sum.count += point.count;
      sum.x += point.x * point.count;
      sum.y += point.y * point.count;
      sum.z += point.z * point.count;
      if (point.sector) {
        sum.sector = point.sector;
      }
    });
    var reach = viewRadius + m * edgePc;
    sums.forEach(function (sum, key) {
      var block = byKey.get(key);
      if (!block) {
        var a = sum.address;
        block = blockAt(a.ring, a.seg, a.slab, m, edgePc, galaxyShape, galaxyShape ? prismDensityCache(m) : null);
        var mid = cellCoordinates(block).cartesian;
        if (Math.hypot(mid[0] - center[0], mid[1] - center[1], mid[2] - center[2]) > reach) {
          return;
        }
        view.blocks.push(block);
      }
      block.filled = sum.count;
      var totalKey = m + "/" + key;
      var total = blockTotals.get(totalKey);
      if (total === undefined) {
        if (blockTotals.size > BLOCK_TOTALS_MAX) {
          blockTotals.clear();
        }
        total = blockSectorCount(block.ring, block.seg, block.slab, m, edgePc, galaxyShape);
        blockTotals.set(totalKey, total);
      }
      block.total = Math.max(sum.count, total);
      block.centroid = [sum.x / sum.count, sum.y / sum.count, sum.z / sum.count];
      if (m === 1) {
        block.sector = sum.sector;
      }
    });
    return view;
  }

  // One mesh's geometry from its blocks, colors and opacities (and each
  // block's fill step, for its tint).
  function setBlockGeometry(mesh, cells, colorsOf, alphasOf) {
    var fillsOf = cells.map(blockFillStep);
    var built = buildPrismGeometry(cells);
    var count = built.owners.length;
    var colors = new Float32Array(count * 3);
    var centers = new Float32Array(count * 3);
    var alphas = new Float32Array(count);
    var fills = new Float32Array(count);
    var normals = built.normals;
    var middles = cells.map(function (cell) { return cellCoordinates(cell).cartesian; });
    for (var v = 0; v < count; v++) {
      var owner = built.owners[v];
      var c = colorsOf[owner];
      var lit = normals[3 * v] * PRISM_LIGHT.x + normals[3 * v + 1] * PRISM_LIGHT.y + normals[3 * v + 2] * PRISM_LIGHT.z;
      var shade = PRISM_AMBIENT + (1 - PRISM_AMBIENT) * Math.max(0, lit);
      colors[3 * v] = c.r * shade;
      colors[3 * v + 1] = c.g * shade;
      colors[3 * v + 2] = c.b * shade;
      var mid = middles[owner];
      centers[3 * v] = mid[0];
      centers[3 * v + 1] = mid[1];
      centers[3 * v + 2] = mid[2];
      alphas[v] = alphasOf[owner];
      fills[v] = fillsOf[owner];
    }
    var geometry = new THREE.BufferGeometry();
    geometry.setAttribute("position", new THREE.BufferAttribute(built.positions, 3));
    geometry.setAttribute("prismColor", new THREE.BufferAttribute(colors, 3));
    geometry.setAttribute("prismCenter", new THREE.BufferAttribute(centers, 3));
    geometry.setAttribute("prismAlpha", new THREE.BufferAttribute(alphas, 1));
    geometry.setAttribute("prismFill", new THREE.BufferAttribute(fills, 1));
    geometry.setAttribute("faceUv", new THREE.BufferAttribute(built.uvs, 2));
    geometry.setIndex(new THREE.BufferAttribute(built.indices, 1));
    mesh.geometry.dispose();
    mesh.geometry = geometry;
    mesh.userData.cells = cells;
    mesh.userData.owners = built.owners;
    mesh.visible = cells.length > 0;
  }

  // TODO(galaxy-map #17): smooth zooming.
  //   - Hand the block listing and geometry to a worker, and keep a small
  //     cache of built meshes keyed by (m, slice, focus cell) for the
  //     current m and one step finer and coarser.
  //   - Once the current view is drawn, prerender the neighbours during
  //     idle time (requestIdleCallback, with a setTimeout fallback).
  //   - On a zoom that changes m, crossfade the old and new meshes over a
  //     few frames instead of swapping.
  // Edge cases:
  // - Cancel stale worker jobs when the camera moves again.
  // - Cap memory (drop the farthest cached mesh).
  // - The first frame must not wait on the worker: draw today's
  //   synchronous result once.
  function updatePrisms(viewRadius) {
    var center = [target.x, target.y, target.z];
    var pcPerPixel = pcPerPixelAtTarget();
    lastPrismViewRadius = viewRadius;
    // The camera's height only decides which top and bottom faces face it,
    // so it is keyed coarsely: a rebuild every 16 pixels' worth.
    var viewerBand = Math.round(camera.position.z / (16 * pcPerPixel));
    var key = [center.join(","), viewRadius, pcPerPixel.toPrecision(3), sliceAtFocus, viewerBand, filledVersion].join(":");
    if (key === prismSetKey) {
      return;
    }
    prismSetKey = key;
    updateWedgeLevels(pcPerPixel);
    var view = viewBlocks(center, viewRadius, pcPerPixel);
    drawnSectorsPerPrism = view.m;

    var solid = [];
    var glass = [];
    view.blocks.forEach(function (cell) {
      cell.opacity = blockOpacity(cell);
      (cell.opacity >= 1 ? solid : glass).push(cell);
    });
    // Back to front, so each translucent block blends over what's behind.
    var eye = camera.position;
    glass.forEach(function (cell) {
      var mid = cellCoordinates(cell).cartesian;
      cell.eyeDistance = Math.hypot(mid[0] - eye.x, mid[1] - eye.y, mid[2] - eye.z);
    });
    glass.sort(function (p, q) { return q.eyeDistance - p.eyeDistance; });

    function colorOf(cell) {
      var color = cell.sector
        ? placedDensityColor(cell.sector, data.referenceDensityPerLy3)
        : prismColor(galaxyShape ? prismShade(cell) : 0.5);
      var mid = cellCoordinates(cell).cartesian;
      return color.multiplyScalar(0.25 + 0.75 * edgeFade(mid[0], mid[1], mid[2]));
    }
    setBlockGeometry(solidMesh, solid, solid.map(colorOf), solid.map(function () { return 1; }));
    setBlockGeometry(glassMesh, glass, glass.map(colorOf), glass.map(function (cell) { return cell.opacity; }));
    updateScaleBar();
  }

  // Inside the solid, the slice is the main answer; this near cut stays
  // for views where the camera is below the cut (or the slice is off).
  function updateNearCut() {
    solidMaterial.uniforms.nearCut.value = orbit.radius * NEAR_CUT;
    glassMaterial.uniforms.nearCut.value = orbit.radius * NEAR_CUT;
  }


  // --- Cube tiles ----------------------------------------------------------
  //
  // The map asks for fixed cubes of space ("tiles") rather than "everything
  // within R of the camera target" -- see stellarObjects.galaxyViewport's
  // "Cube tiles" section. Space is an octree: level 0 is one cube
  // tileRootEdgePc on a side centered on the galactic origin, each level
  // halves the edge, and a tile's key is "level/ix/iy/iz". neededTiles()
  // picks the smallest level whose tiles are at least the view radius
  // across (so at most 27 tiles cover the view), plus the smallest tiles
  // (the only ones that list planned slots) around the target when zoomed
  // in. lib/galaxymap3d.py's initial_tile_request does the same in Python
  // for the first frame.
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

  // TODO(galaxy-map #17): also compute the tiles for the next zoom step in
  // and out (tile level +/- 1 around the same target) and fetch them at low
  // priority after the current view's tiles, so zooming never waits on the
  // network. Never let the prefetch abort a real fetch.
  function neededTiles() {
    var viewRadius = Math.max(MIN_RADIUS, Math.min(MAX_RADIUS * FETCH_RADIUS_FACTOR, orbit.radius * FETCH_RADIUS_FACTOR));
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

  // --- Drawing from tiles --------------------------------------------------

  // Filled points per tile object (a refetched tile is a new object).
  var filledPointsByTile = new WeakMap();
  var filledSignature = "";

  // Draws whatever of the needed tiles is already cached; returns what's
  // still missing.
  function renderFromCache(need) {
    fadeRadius = need.viewRadius;
    var points = [];
    var present = [];
    var missing = [];
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
    });
    var signature = currentStamp + "|" + present.join(",");
    if (signature !== filledSignature) {
      filledSignature = signature;
      filledPoints = points;
      filledVersion++;
    }
    updatePrisms(need.viewRadius);
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
    if (!missing.length) {
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

  // First frame: normally every tile is already cached (the page embeds
  // them), so this only fetches if the page's own tile fetch came up short.
  var initialMissing = renderFromCache(neededTiles());
  if (initialMissing.length) {
    scheduleFetch(false);
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

  // TODO(galaxy-map #17): animate zoom steps (wheel, buttons,
  // double-click) over ~150 ms toward the new radius instead of jumping,
  // so the prerendered neighbour level can fade in. Respect
  // prefers-reduced-motion by jumping as today.
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
      setRadius(orbit.radius * Math.exp(deltaPx * WHEEL_ZOOM_PER_PX));
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
  // generated. Blocks hidden by the near cut (updateNearCut) are skipped.
  function cellAtClientPoint(clientX, clientY) {
    var meshes = [solidMesh, glassMesh].filter(function (mesh) { return mesh.visible && mesh.userData.owners; });
    if (!meshes.length) {
      return null;
    }
    var ndc = ndcFromClientPoint(clientX, clientY);
    if (!ndc) {
      return null;
    }
    raycaster.setFromCamera(ndc, camera);
    var hits = raycaster.intersectObjects(meshes, false);
    var nearCut = solidMaterial.uniforms.nearCut.value;
    var first = null;
    for (var i = 0; i < hits.length; i++) {
      var face = hits[i].face;
      if (!face) {
        continue;
      }
      var mesh = hits[i].object;
      var cellData = mesh.userData.cells[mesh.userData.owners[face.a]];
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

  function blockHit(cellData) {
    var m = drawnSectorsPerPrism;
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
  function selectCell(point, cell) {
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
  function centerOn(point, cell) {
    target.copy(point);
    applyCamera();
    updateScaleBar();
    selectCell(point, cell);
    scheduleFetch(true);
  }

  // Double click: centers AND zooms in by one clickZoomFactor step. A
  // double-click's first click has already centered on the point (see the
  // click handler), so this zooms in on that same point rather than
  // re-picking under the cursor, which by now is over something else.
  function zoomInOnTarget() {
    setRadius(orbit.radius / clickZoomFactor(orbit.radius));
    scheduleFetch(true);
  }

  // Resolves a click/double-click's target point the same way for both:
  // a hit block (its filled sectors' mean position when it has any, so a
  // double-click zooms in toward them), or (empty space)
  // depthPointAtClientPoint -- shared so the click handlers below never
  // have to duplicate the raycast-then-fall-back logic.
  function resolveClickTarget(clientX, clientY) {
    var hit = cellAtClientPoint(clientX, clientY);
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
      centerOn(resolved.point, resolved.cell);
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
        centerOn(resolved.point, resolved.cell);
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
          setRadius(orbit.radius / clickZoomFactor(orbit.radius));
          scheduleFetch(true);
        } else if (action === "zoom-out") {
          setRadius(orbit.radius * clickZoomFactor(orbit.radius));
          scheduleFetch(true);
        } else if (action === "reset") {
          resetView();
        } else if (action === "wedges") {
          wedgeGroup.visible = !wedgeGroup.visible;
          button.setAttribute("aria-pressed", String(wedgeGroup.visible));
        } else if (action === "slice") {
          sliceAtFocus = !sliceAtFocus;
          button.setAttribute("aria-pressed", String(sliceAtFocus));
          if (lastPrismViewRadius !== null) {
            updatePrisms(lastPrismViewRadius);
          }
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

  (function animate() {
    requestAnimationFrame(animate);
    updateHighlightScale();
    updateWedgeLabels();
    updateNearCut();
    renderer.render(scene, camera);
  })();
}

if (canvas && sceneData) {
  initGalaxyMap3d(canvas, sceneData);
}
